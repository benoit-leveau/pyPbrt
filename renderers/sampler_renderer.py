"""SamplerRenderer Class."""


from core.pbrt import round_up_pow_2
from core.renderer import Renderer
from core.parallel import Task, num_system_cores
from core.rng import RNG
from core.spectrum import Spectrum

from core.logger import logger


class SamplerRenderer(Renderer):

    """Renderer driven by a stream of samples."""

    def __init__(self, sampler, camera, si, vi):
        """Default constructor for SamplerRenderer."""
        super(SamplererRenderer, self).__init__()
        self.sampler = sampler
        self.camera = camera
        self.surface_integrator = si
        self.volume_integrator = vi
        
    def render(self, scene):
        """Render a scene."""

        # allow integrators to do preprocessing
        self.surface_integrator.preprocess(scene, camera, self)
        self.volume_integrator.preprocess(scene, camera, self)

        # allocate and initialize smaple
        self.sample = Sample(self.sampler,
                             self.surface_integrator,
                             self.volume_integrator,
                             self.scene)

        # compute number of SamplerRendererTasks to create for rendering
        n_pixels = self.camera.film.x_resolution * self.camera.film.y_resolution
        n_tasks = max(32 * num_system_cores(), n_pixels / (16*16))
        n_tasks = round_up_pow_2(n_tasks)

        # create and launch SamplerRendererTasks for rendering image
        render_tasks = []
        for i in range(n_tasks):
            render_tasks.append(SamplerRendererTask(self.scene,
                                                    self,
                                                    self.camera,
                                                    self.sampler,
                                                    self.sample,
                                                    n_tasks-1-i,
                                                    n_tasks))

        self.enqueue_tasks(render_tasks)
        self.wait_for_all_tasks()
        
        # clean up after rendering and store final image
        del(self.sample)
        self.camera.film.write_image()

    def Li(self, scene, ray, sample, rng, isect=None, T=None):
        """Compute the incident radiance along a given ray."""
        # allocate local variables for isect and T if needed
        if not T:
            T = Spectrum()
        if not isect:
            isect = Intersection()
        Li = Spectrum(0.0)
        hit, isect = self.scene.intersect(ray)
        if hit:
            Li = self.surface_integrator.Li(self.scene,
                                            self,
                                            ray,
                                            isect,
                                            sample,
                                            rng)
        else:
            # handle ray that doesn't intersect any geometry
            for light in self.scene.lights:
                Li += light.Le(ray)
        
        Lvi = self.volume_integrator.Li(self.scene,
                                        self,
                                        ray,
                                        sample,
                                        rng,
                                        T,
                                        arena)

        return T * Li + Lvi
                                        

    def transmittance(self, scene, ray, sample, rng):
        """Compute the attenuation by volumetric scattering along a ray."""
        return self.volume_integrator.transmittance(self.scene,
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

        # get sub-sampler for SamplerRendererTask
        sampler = self.main_sampler.get_sub_sampler(self.task_num,
                                                    self.task_count)
        if not sampler:
            return
        
        # Declare local variables used for rendering loop
        rng = RNG(self.task_num)

        # allocate space for samples and intersections
        max_samples = self.sampler.maximum_sample_count()
        samples = self.orig_sample.duplicate(max_samples)
        rays = [] # RayDifferential[max_samples]
        Ls = [] # Spectrum[max_samples]
        Ts = [] # Spectrum[max_samples]
        isects = [] # Intersection[max_samples]

        # get samples from Sampler and update image
        while True:
            sample_count = sampler.get_more_samples(samples, rng)

            # if no more samples to compute, exit
            if sample_count <= 0:
                break
            
            # generate camera rays and compute radiance along rays
            for i in range(sample_count):
                # find camera ray for samples[i]
                ray_weight = self.camera.generate_ray_differential(samples[i],
                                                                   rays[i])
                coeff = 1.0 / math.sqrt(sampler.samples_per_pixel)
                rays[i].scale_differentials(coeff)
                
                # evaluate radiance along camera ray
                if ray_weight > 0.0:
                    radiance = self.renderer.Li(self.scene,
                                                rays[i],
                                                samples[i],
                                                rng,
                                                isects[i],
                                                Ts[i]
                                                )
                    Ls[i] = ray_weight * radiance
                else:
                    Ls[i] = 0.0
                    Ts[i] = 1.0

                # check for unexpected radiance values
                if Ls[i].has_nan():
                    logger.error(
                        "Not-a-number radiance value returned for image sample.  Setting to black.")
                    Ls[i] = Spectrum(0.0)
                elif Ls[i].y() < -1e-5:
                    logger.error(
                        "Negative luminance value, %f, returned for image sample.  Setting to black." % Ls[i].y())
                    Ls[i] = Spectrum(0.0)
                elif is_inf(Ls[i].y()):
                    logger.error(
                        "Infinite luminance value returned for image sample.  Setting to black.")
                    Ls[i] = Spectrum(0.0)
                
            # report sample results to Sampler, add contributions to image
            if sampler.report_results(samples, rays, Ls, isects, sample_count):
                for i in range(sample_count):
                    self.camera.film.add_sample(samples[i], Ls[i])

        # clean up after SamplerRendererTask is done with its image region
        pass
