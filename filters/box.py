"""BoxFilter Class."""

from core.filter import Filter


class BoxFilter(Filter):

    """Class describing a Box Filter."""

    def __init__(self, x_width, y_width):
        """Default constructor for BoxFilter."""
        super(BoxFilter, self).__init__(x_width, y_width)

    def evaluate(self, x, y):
        """Evaluate filter at given position."""
        return 1.0

