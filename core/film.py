"""Film Class."""

from abc import ABCMeta, abstractmethod


class Film(object):

    """Class describing a Film."""

    __metaclass__ = ABCMeta
    
    def __init__(self, x_res, y_res):
        """Default constructor for Film."""
        self.x_resolution = x_res
        self.y_resolution = y_res

    @abstractmethod
    def add_sample(self, sample, L):
        pass

    @abstractmethod
    def splat(self, sample, L):
        pass

    @abstractmethod
    def get_sample_extent(self):
        pass

    @abstractmethod
    def get_pixel_extent(self):
        pass
    
    def update_display(self, x0, y0, x1, y1, splat_scale=1.0):
        pass

    @abstractmethod
    def write_image(self, splat_scale=1.0):
        """Write the image to disk."""
        pass
