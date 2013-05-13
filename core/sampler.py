"""Interface classes for samplers."""

from abc import ABCMeta, abstractmethod

from core.pbrt import lerp


class Sampler(object):

    """Class describing a sampler."""

    __metaclass__ = ABCMeta
    
    def __init__(self, x_start, x_end, y_start, y_end, spp, s_open, s_close):
        """Default constructor for a sampler."""
        self.x_pixel_start = x_start
        self.x_pixel_end = x_end
        self.y_pixel_start = y_start
        self.y_pixel_end = y_end
        self.samples_per_pixel = spp
        self.shutter_open = s_open
        self.shutter_close = s_close

    @abstractmethod
    def get_more_samples(self, sample, rng):
        """Get more samples from the sampler."""
        pass

    @abstractmethod
    def maximum_sample_count(self):
        """Return the number of samples the sampler can return in one call."""
        pass

    def report_results(self, samples, rays, ls, intersections, count):
        """Read the samples back from the renderer.

        Return True to tell the renderer to keep the samples.
        Return False to tell the renderer to discard the samples.

        """
        return True

    @abstractmethod
    def get_sub_sampler(self, num, count):
        """Get a sub sampler."""
        pass

    @abstractmethod
    def round_size(self, size):
        """."""
        pass

    def compute_sub_window(self, num, count):
        """Compute the sub window position."""
        dx = self.x_pixel_end - self.x_pixel_start
        dy = self.y_pixel_end - self.y_pixel_start
        nx = count
        ny = 1
        while(((nx & 0x1) == 0) and ((2*dx*ny) < (dy*nx))):
            nx >>= 1
            ny <<= 1

        # compute x and y pixel sample range for sub-window
        xo = num % nx
        yo = num % ny
        tx0 = float(xo) / float(nx)
        tx1 = float(xo+1) / float(nx)
        ty0 = float(yo) / float(ny)
        ty1 = float(yo+1) / float(ny)

        x_start = int(lerp(tx0, self.x_pixel_start, self.x_pixel_end))
        x_end   = int(lerp(tx1, self.x_pixel_start, self.x_pixel_end))
        y_start = int(lerp(ty0, self.y_pixel_start, self.y_pixel_end))
        y_end   = int(lerp(ty1, self.y_pixel_start, self.y_pixel_end))
        
        return x_start, x_end, y_start, y_end


class CameraSample(object):

    """Class describing a camera sample."""

    def __init__(self,
                 image_x=0.0, image_y=0.0,
                 lens_u=0.0, lens_v=0.0,
                 time=0.0):
        """Default constructor for CameraSample."""
        self.image_x = image_x
        self.image_y = image_y
        self.lens_u = lens_u
        self.lens_v = lens_v
        self.time = time

    def __str__(self):
        """Return a string describing the sample."""
        return "CameraSample (img=(%f,%f) lns=(%f,%f) t=%f)" % \
               (self.image_x, self.image_y,
                self.lens_u, self.lens_v,
                self.time)


class Sample(CameraSample):

    """Class describing a generic sample."""

    def __init__(self):
        """Default constructor for Sample."""
        super(Sample, self).__init__()
        self.n1D = []
        self.n2D = []
        self.oneD = None
        self.twoD = None        

    @classmethod
    def from_sampler(cls,
                     sampler,
                     surface_integrator, volume_integrator,
                     scene):
        """Construct a Sample by querying the Sampler."""
        sample = cls()
        if surface_integrator:
           surface_integrator.request_samples(sampler, sample, scene)
        if volume_integrator:
            volume_integrator.request_samples(sampler, sample, scene)
        sample.allocate_sample_memory()
        return sample

    def add_1D(self, num):
        """Add a sample."""
        self.n1D.append(num)
        return len(n1D)-1

    def add_2D(self, num):
        """Add a sample."""
        self.n2D.append(num)
        return len(n2D)-1

    def duplicate(self, count):
        """Return a list of new duplicate samples."""
        ret = []
        for i in xrange(count):
            sample = Sample()
            sample.n1D = list(self.n1D)
            sample.n2D = list(self.n2D)
            ret.append(sample)
        return ret

    def allocate_sample_memory(self):
        """Reserve size for sample."""
        pass
