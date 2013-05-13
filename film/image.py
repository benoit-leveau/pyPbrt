"""ImageFilm Class."""

import math
from PIL.Image import new

from core.film import Film
from core.logger import logger
from core.spectrum import xyz_to_rgb


class ImageFilm(Film):

    """Class describing an ImageFilm."""

    def __init__(self, x_res, y_res, filter, crop_window, filename, open_window):
        """Default constructor for ImageFilm."""
        super(ImageFilm, self).__init__(x_res, y_res)
        self.filter = filter
        self.crop_window = list(crop_window)
        self.filename = filename

        # compute film image extent
        self.x_pixel_start = int(math.ceil(
            self.x_resolution * self.crop_window[0]))
        self.x_pixel_count = max(1, int(math.ceil(
            self.x_resolution * crop_window[1]) - self.x_pixel_start))
        self.y_pixel_start = int(math.ceil(
            self.y_resolution * self.crop_window[2]))
        self.y_pixel_count = max(1, int(math.ceil(
            self.y_resolution * crop_window[3]) - self.y_pixel_start))

        if self.x_pixel_count>self.x_resolution:
            logger.warning("ImageFilm.x_pixel_count is incorrect: %d" %
                           self.x_pixel_count)

        if self.y_pixel_count>self.y_resolution:
            logger.warning("ImageFilm.y_pixel_count is incorrect: %d" %
                           self.y_pixel_count)

        # allocate film image storage
        self.pixels = [0.0] * (self.x_resolution * self.y_resolution * 4)
        # for y in xrange(self.y_resolution):
        #    row = []
        #    for x in xrange(self.x_resolution):
        #        row.append(Pixel())
        #    self.pixels.append(row)

        # precompute filter weight table
        # pass

        # possibly open window for image display
        if open_window: # or pbrt_options.open_window:
            logger.warning("Support for opening image display window not available in this build.")
    
    def add_sample(self, sample, L):
        """Add a sample to the film."""
        d_image_x = sample.image_x - 0.5
        d_image_y = sample.image_y - 0.5
        x0 = int(math.ceil(d_image_x - self.filter.x_width))
        x1 = int(d_image_x + self.filter.x_width)
        y0 = int(math.ceil(d_image_y - self.filter.y_width))
        y1 = int(d_image_y + self.filter.y_width)
        x0 = max(x0, self.x_pixel_start)
        x1 = min(x1, self.x_pixel_start + self.x_resolution - 1)
        y0 = max(y0, self.y_pixel_start)
        y1 = min(y1, self.y_pixel_start + self.y_resolution - 1)
        if ((x1-x0)<0) or ((y1-y0)<0):
            logger.warning("Sample outside the image extent. %d %d %d %d" % (x0, x1, y0, y1))
            return
        # loop over filter support and add sample to pixel arrays
        x, y, z = L.to_xyz()

        for j in xrange(y0, y1+1):
            for i in xrange(x0, x1+1):
                filter_weight = self.filter.evaluate(i, j)
                self.pixels[j*self.x_resolution*4+i*4] = x
                self.pixels[j*self.x_resolution*4+i*4+1] = y
                self.pixels[j*self.x_resolution*4+i*4+2] = z
                self.pixels[j*self.x_resolution*4+i*4+3] = filter_weight

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
        n_pixels = self.x_resolution * self.y_resolution
        image = new("RGB", (self.x_resolution, self.y_resolution))
        offset = 0
        rgb = []
        for j in xrange(self.y_resolution):
            for i in xrange(self.x_resolution):
                # convert pixel XYZ color to RGB
                x = self.pixels[j*self.x_resolution*4+i*4]
                y = self.pixels[j*self.x_resolution*4+i*4+1]
                z = self.pixels[j*self.x_resolution*4+i*4+2]
                r, g, b = xyz_to_rgb(x, y, z)

                # normalize pixel with weight sum
                weight_sum = self.pixels[j*self.x_resolution*4+i*4+3]
                if weight_sum != 0.0:
                    inv_weight = 1.0 / weight_sum
                    r = max(0.0, r * inv_weight)
                    g = max(0.0, g * inv_weight)
                    b = max(0.0, b * inv_weight)

                if False: #splat_scale != 0.0:
                    # add splat value at pixel
                    x, y, z = self.pixels[j*self.x_resolution*4+i*4].splat_xyz
                    splat_r, splat_g, splat_b = xyz_to_rgb(x, y, z)
                    r += splat_scale * splat_r
                    g += splat_scale * splat_g
                    b += splat_scale * splat_b

                r = int(255.0 * r)
                g = int(255.0 * g)
                b = int(255.0 * b)

                image.im.putpixel((i, j), (r, g, b))

        image.save(self.filename)

class Pixel(object):

    """Class describing a Pixel."""

    def __init__(self):
        """Default constructor for Pixel."""
        self.Lxyz = [0.0, 0.0, 0.0]
        self.weight_sum = 0.0
        self.splat_xyz = [0.0, 0.0, 0.0]
        self.pad = 0.0

