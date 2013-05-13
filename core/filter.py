"""Filter Class."""

from abc import ABCMeta, abstractmethod


class Filter(object):

    """Class describing a Filter."""

    __metaclass__ = ABCMeta
    
    def __init__(self, x_width, y_width):
        """Default constructor for Filter."""
        self.x_width = x_width
        self.y_width = y_width
        self.inv_x_width = 1.0 / self.x_width
        self.inv_y_width = 1.0 / self.y_width

    @abstractmethod
    def evaluate(self, x, y):
        """Evaluate filter at given position."""
        pass
