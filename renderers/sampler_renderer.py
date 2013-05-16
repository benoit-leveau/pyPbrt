"""SamplerRenderer Class."""

import math

from core.pbrt import round_up_pow_2
from core.renderer import Renderer
from core.parallel import Task, num_system_cores
from core.rng import RNG
from core.spectrum import Spectrum
from core.sampler import Sample
from core.logger import logger
from core.intersection import Intersection


class SamplerRenderer(Renderer):

    """Renderer driven by a stream of samples."""

    def __init__(self, sampler, camera, surface_integrator, volume_integrator):
        """Default constructor for SamplerRenderer."""
        super(SamplerRenderer, self).__init__()
        self.sampler = sampler
        self.camera = camera
        self.surface_integrator = surface_integrator
        self.volume_integrator = volume_integrator
        
    def render(self, scene):
        """Render a scene."""

        # allow integrators to do preprocessing
        self.surface_integrator.preprocess(scene, self.camera, self)
        self.volume_integrator.preprocess(scene, self.camera, self)

        # allocate and initialize smaple
        self.sample = Sample.from_sampler(
            self.sampler,
            self.surface_integrator,
            self.volume_integrator,
            scene)

        # compute number of SamplerRendererTasks to create for rendering
        n_pixels = self.camera.film.x_resolution * self.camera.film.y_resolution
        n_tasks = max(32 * num_system_cores(), n_pixels / (16*16))
        n_tasks = round_up_pow_2(n_tasks)

        # create and launch SamplerRendererTasks for rendering image
        render_tasks = []
        for i in range(n_tasks):
            render_tasks.append(SamplerRendererTask(scene,
                                                    self,
                                                    self.camera,
                                                    self.sampler,
                                                    self.sample,
                                                    n_tasks-1-i,
                                                    n_tasks))

        # self.enqueue_tasks(render_tasks)
        # self.wait_for_all_tasks()
        for task in render_tasks:
            task.run()
        
        # clean up after rendering and store final image
        del(self.sample)
        self.camera.film.write_image()

    def Li(self, scene, ray, sample, rng, intersection=None, T=None):
        """Compute the incident radiance along a given ray."""
        # allocate local variables for isect and T if needed
        if not T:
            T = Spectrum(0.0)
        if not intersection:
            intersection = Intersection()
        Li = Spectrum(0.0)
        hit = scene.intersect(ray, intersection)
        if hit:
            Li = self.surface_integrator.Li(scene,
                                            self,
                                            ray,
                                            intersection,
                                            sample,
                                            rng)
        else:
            # handle ray that doesn't intersect any geometry
            for light in scene.lights:
                Li += light.Le(ray)
        
        # Lvi = self.volume_integrator.Li(scene, self, ray, sample, rng, T)
        # return T * Li + Lvi, intersection

        return Li, intersection, T
                                        

    def transmittance(self, scene, ray, sample, rng):
        """Compute the attenuation by volumetric scattering along a ray."""
        return self.volume_integrator.transmittance(scene,
                                                    self,
                                                    ray,
                                                    sample,
                                                    rng)


class SamplerRendererTask(Task):
    
    """Single Render Task for the SamplerRenderer."""

    def __init__(self, scene, renderer, camera,
                 sampler, sample, task_num, task_count):
        """Default constructor for SamplerRendererTask."""
        self.scene = scene
        self.renderer = renderer
        self.camera = camera
        self.main_sampler = sampler
        self.orig_sample = sample
        self.task_num = task_num
        self.task_count = task_count

    def run(self):
        """Execute the task."""

        print "executing task %d/%d" % (self.task_num, self.task_count)
        
        # get sub-sampler for SamplerRendererTask
        sampler = self.main_sampler.get_sub_sampler(self.task_num,
                                                    self.task_count)
        if not sampler:
            return
        
        # Declare local variables used for rendering loop
        rng = RNG(self.task_num)

        # allocate space for samples and intersections
        max_samples = sampler.maximum_sample_count()
        samples = self.orig_sample.duplicate(max_samples)
        rays = [] # RayDifferential[max_samples]
        Ls = [] # Spectrum[max_samples]
        Ts = [] # Spectrum[max_samples]
        isects = [] # Intersection[max_samples]
        for i in range(max_samples):
            isects.append(Intersection())

        # get samples from Sampler and update image
        while True:
            sample_count = sampler.get_more_samples(samples, rng)

            # if no more samples to compute, exit
            if sample_count <= 0:
                break
            
            # generate camera rays and compute radiance along rays
            for i in range(sample_count):
                # find camera ray for samples[i]
                ray_weight, ray_diff = self.camera.generate_ray_differential(samples[i])
                
                rays.append(ray_diff)
                coeff = 1.0 / math.sqrt(sampler.samples_per_pixel)
                ray_diff.scale_differentials(coeff)
                
                # evaluate radiance along camera ray
                if ray_weight > 0.0:
                    radiance, intersection, Ts_i = \
                              self.renderer.Li(self.scene,
                                               ray_diff,
                                               samples[i],
                                               rng,
                                               isects[i])
                    Ls_i = ray_weight * radiance
                else:
                    Ls_i = 0.0
                    Ts_i = 1.0

                Ls.append(Ls_i)
                
                # check for unexpected radiance values
                if Ls_i.has_nan():
                    logger.error(
                        "Not-a-number radiance value returned for image sample.  Setting to black.")
                    Ls_i = Spectrum(0.0)
                elif Ls_i.y() < -1e-5:
                    logger.error(
                        "Negative luminance value, %f, returned for image sample.  Setting to black." % Ls_i.y())
                    Ls_i = Spectrum(0.0)
                elif Ls_i.y() == float('inf'):
                    logger.error(
                        "Infinite luminance value returned for image sample.  Setting to black.")
                    Ls_i = Spectrum(0.0)

            # report sample results to Sampler, add contributions to image
            if sampler.report_results(samples, rays, Ls, isects, sample_count):
                for i in range(sample_count):
                    self.camera.film.add_sample(samples[i], Ls[i])

        # clean up after SamplerRendererTask is done with its image region
        pass
