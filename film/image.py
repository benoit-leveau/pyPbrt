"""ImageFilm Class."""

from core.film import Film


class ImageFilm(Film):

    """Class describing an ImageFilm."""

    def __init__(self, x_res, y_res, filter, crop, filename, open_window):
        """Default constructor for ImageFilm."""
        super(ImageFilm, self).__init__(x_res, y_res)

    def add_sample(self, sample, L):
        pass

    def splat(self, sample, L):
        pass

    def get_sample_extent(self):
        pass

    def get_pixel_extent(self):
        pass
    
    def update_display(self, x0, y0, x1, y1, splat_scale=1.0):
        pass

    def write_image(self, splat_scale=1.0):
        """Write the image to disk."""
        pass
