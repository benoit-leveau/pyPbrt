"""ImageFilm Class."""

import math
from PIL.Image import new

from core.film import Film
from core.logger import logger
from core.spectrum import xyz_to_rgb


FILTER_TABLE_SIZE = 16


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
        self.pixels = [0.0] * (self.x_pixel_count * self.y_pixel_count * 8)

        # precompute filter weight table
        self.filter_table = [0.0] * (FILTER_TABLE_SIZE * FILTER_TABLE_SIZE)
        for y in xrange(FILTER_TABLE_SIZE):
            fy = (float(y) + 0.5) * self.filter.y_width / float(FILTER_TABLE_SIZE)
            for x in xrange(FILTER_TABLE_SIZE):
                fx = (float(x) + 0.5) * self.filter.x_width / float(FILTER_TABLE_SIZE)
                self.filter_table[y * FILTER_TABLE_SIZE + x] = \
                                    self.filter.evaluate(fx, fy)

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
        x1 = min(x1, self.x_pixel_start + self.x_pixel_count - 1)
        y0 = max(y0, self.y_pixel_start)
        y1 = min(y1, self.y_pixel_start + self.y_pixel_count - 1)
        if ((x1-x0)<0) or ((y1-y0)<0):
            logger.warning("Sample outside the image extent. %d %d %d %d" % (x0, x1, y0, y1))
            return

        # loop over filter support and add sample to pixel arrays
        L_x, L_y, L_z = L.to_xyz()

        # precompute x and y filter table offsets
        ifx = []
        for x in xrange(x0, x1+1):
            fx = abs((x - d_image_x) * self.filter.inv_x_width * FILTER_TABLE_SIZE)
            ifx.append(min(int(fx), FILTER_TABLE_SIZE-1))
        ify = []
        for y in xrange(y0, y1+1):
            fy = abs((y - d_image_y) * self.filter.inv_y_width * FILTER_TABLE_SIZE)
            ify.append(min(int(fy), FILTER_TABLE_SIZE-1))
        sync_needed = (self.filter.x_width > 0.5) or (self.filter.y_width > 0.5)
        for y in xrange(y0, y1+1):
            for x in xrange(x0, x1+1):
                # evaluate filter value at (x,y) pixel
                offset = ify[y-y0]*FILTER_TABLE_SIZE + ifx[x-x0]
                filter_weight = self.filter_table[offset]
                # update pixel values with filtered sample contribution

                if not sync_needed:
                    self.pixels[y*self.x_pixel_count*8+x*8] += filter_weight * L_x
                    self.pixels[y*self.x_pixel_count*8+x*8+1] += filter_weight * L_y
                    self.pixels[y*self.x_pixel_count*8+x*8+2] += filter_weight * L_z
                    self.pixels[y*self.x_pixel_count*8+x*8+3] += filter_weight
                else:
                    self.pixels[y*self.x_pixel_count*8+x*8] += filter_weight * L_x
                    self.pixels[y*self.x_pixel_count*8+x*8+1] += filter_weight * L_y
                    self.pixels[y*self.x_pixel_count*8+x*8+2] += filter_weight * L_z
                    self.pixels[y*self.x_pixel_count*8+x*8+3] += filter_weight

    def splat(self, sample, L):
        if L.has_nan():
            logger.warning("ImageFilm ignoring splatted spectrum with NaN values")
            return
        L_x, L_y, L_z = L.to_xyz()
        x = int(sample.image_x)
        y = int(sample.image_y)
        if (x < self.x_pixel_start) or (x - self.x_pixel_start >= self.x_pixel_count) or (y < self.y_pixel_start) or (y - self.y_pixel_start >= self.y_pixel_count):
            return
        y = y - self.y_pixel_start
        x = x - self.x_pixel_start
        self.pixels[y*self.x_pixel_count*8+8*x+4] = L_x
        self.pixels[y*self.x_pixel_count*8+8*x+5] = L_y
        self.pixels[y*self.x_pixel_count*8+8*x+6] = L_z

    def get_sample_extent(self):
        x_start = int(self.x_pixel_start + 0.5 - self.filter.x_width)
        x_end = int(math.ceil(self.x_pixel_start + 0.5 + self.x_pixel_count + self.filter.x_width))
        y_start = int(self.y_pixel_start + 0.5 - self.filter.y_width)
        y_end = int(math.ceil(self.y_pixel_start + 0.5 + self.y_pixel_count + self.filter.y_width))
        return x_start, x_end, y_start, y_end

    def get_pixel_extent(self):
        x_start = self.x_pixel_start
        x_end = self.x_pixel_start + self.x_pixel_count
        y_start = self.y_pixel_start
        y_end = self.y_pixel_start + self.y_pixel_count
        return x_start, x_end, y_start, y_end
    
    def update_display(self, x0, y0, x1, y1, splat_scale=1.0):
        pass

    def write_image(self, splat_scale=1.0):
        """Write the image to disk."""
        # convert image to RGB and compute final pixel values
        n_pixels = self.x_pixel_count * self.y_pixel_count
        image = new("RGB", (self.x_pixel_count, self.y_pixel_count))
        offset = 0
        for j in xrange(self.y_pixel_count):
            for i in xrange(self.x_pixel_count):
                # convert pixel XYZ color to RGB
                x = self.pixels[j*self.x_pixel_count*8+i*8]
                y = self.pixels[j*self.x_pixel_count*8+i*8+1]
                z = self.pixels[j*self.x_pixel_count*8+i*8+2]
                r, g, b = xyz_to_rgb(x, y, z)

                # normalize pixel with weight sum
                weight_sum = self.pixels[j*self.x_pixel_count*8+i*8+3]

                if weight_sum != 0.0:
                    inv_weight = 1.0 / weight_sum
                    r = max(0.0, r * inv_weight)
                    g = max(0.0, g * inv_weight)
                    b = max(0.0, b * inv_weight)

                # add splat value at pixel
                x = self.pixels[j*self.x_pixel_count*8+i*8+4]
                y = self.pixels[j*self.x_pixel_count*8+i*8+5]
                z = self.pixels[j*self.x_pixel_count*8+i*8+6]
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

