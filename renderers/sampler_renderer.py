"""SamplerRenderer Class."""


from core.pbrt import round_up_pow_2
from core.renderer import Renderer
from core.parallel import Task, num_system_cores
from core.rng import Rng

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
        rnd = Rng(self.task_num)

        # allocate space for samples and intersections

        # get samples from Sampler and update image

        # clean up after SamplerRendererTask is done with its image region
        
        
