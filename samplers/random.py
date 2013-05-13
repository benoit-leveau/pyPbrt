"""RandomSampler Class."""

from core.sampler import Sampler
from core.rng import RNG


class RandomSampler(Sampler):

    """Class describing a random sampler."""

    def __init__(self,
                 x_start, x_end,
                 y_start, y_end,
                 n_samples,
                 s_open, s_close):
        """Default constructor for the random sampler."""
        super(RandomSampler, self).__init__(x_start, x_end,
                                            y_start, y_end,
                                            n_samples, s_open, s_close)
        self.x_pos = self.x_pixel_start
        self.y_pos = self.y_pixel_start
        self.n_samples = n_samples

        # get storage for a pixel's worth of stratified samples
        self.image_samples = [0.0] * (5*self.n_samples)

        rng = RNG(x_start + y_start * (x_end-x_start))
        for i in xrange(5*self.n_samples):
            self.image_samples[i] = rng.random_float()

        # shift image samples to pixel coordinates
        for o in xrange(0, 2*self.n_samples, 2):
            self.image_samples[o] += self.x_pos
            self.image_samples[o+1] += self.y_pos

        self.sample_pos = 0
        
    def get_more_samples(self, samples, rng):
        """Get more samples from the sampler."""
        if self.sample_pos == self.n_samples:
            # generate new samples
            if (self.x_pixel_start == self.x_pixel_end) or \
               (self.y_pixel_start == self.y_pixel_end):
                return 0
            self.x_pos += 1
            if self.x_pos == self.x_pixel_end:
                self.x_pos = self.x_pixel_start
                self.y_pos += 1
            if self.y_pos == self.y_pixel_end:
                return 0

            for i in range(5*self.n_samples):
                self.image_samples[i] = rng.random_float()

            # shift image samples to pixel coordinates
            for o in range(0, 2*self.n_samples, 2):
                self.image_samples[o] += self.x_pos
                self.image_samples[o+1] += self.y_pos

            self.sample_pos = 0

        # return next sample point
        sample = samples[0]
        sample.image_x = self.image_samples[2*self.sample_pos]
        sample.image_y = self.image_samples[2*self.sample_pos + 1]
        sample.lens_u  = self.image_samples[2*self.n_samples + \
                                            2*self.sample_pos]
        sample.lens_v  = self.image_samples[2*self.n_samples + \
                                            2*self.sample_pos + 1]
        sample.time    = self.image_samples[4*self.n_samples + \
                                            self.sample_pos]

        # generate stratified samples for integrators
        for i in range(len(sample.n1D)):
            for j in range(sample.n1D[i]):
                sample.oneD[i][j] = rng.random_float()
        for i in range(len(sample.n2D)):
            for j in range(2*sample.n2D[i]):
                sample.twoD[i][j] = rng.random_float()

        self.sample_pos += 1

        # print sample
        
        return 1

    def maximum_sample_count(self):
        """Return the number of samples the sampler can return in one call."""
        return 1

    def get_sub_sampler(self, num, count):
        """Get a sub sampler."""
        x0, x1, y0, y1 = self.compute_sub_window(num, count)
        if (x0 == x1) or (y0 == y1):
            return None
        return RandomSampler(x0, x1, y0, y1,
                             self.n_samples,
                             self.shutter_open, self.shutter_close)

    def round_size(self, size):
        """Return the size."""
        return size
