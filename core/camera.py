"""Camera Classes."""

from abc import ABCMeta, abstractmethod

from logger import logger
from core.transform import translate, scale, inverse
from core.geometry import Point, RayDifferential
from core.sampler import CameraSample


class Camera(object):

    """Interface Class for all cameras."""

    __metaclass__ = ABCMeta
    
    def __init__(self, camera_to_world, s_open, s_close, film):
        """Default constructor for Camera."""
        self.camera_to_world = camera_to_world
        self.shutter_open = s_open
        self.shutter_close = s_close
        self.film = film
        if camera_to_world.has_scale():
            logger.warning(
                "Scaling detected in world-to-camera transformation!\n" \
                "The system has numerous assumptions, implicit and explicit,\n" \
                "that this transform will have no scale factors in it.\n" \
                "Proceed at your own risk; your image may have errors or\n" \
                "the system may crash as a result of this.")
        
    @abstractmethod
    def generate_ray(self, sample):
        """Generate a Ray from the camera."""
        pass

    def generate_ray_differential(self, sample):
        """Generate a RayDifferential from the camera."""
        # generate the ray
        weight, ray = self.generate_ray(sample)
        ray_diff = RayDifferential.from_ray(ray)
        
        # find ray after shifting one pixel in the x direction
        sshift = CameraSample.from_sample(sample)
        sshift.image_x += 1
        weight_x, ray_x = self.generate_ray(sshift)
        ray_diff.rx_origin = ray_x.o
        ray_diff.rx_direction = ray_x.d
        
        # find ray after shifting one pixel in the y direction
        sshift.image_x -= 1
        sshift.image_y += 1
        weight_y, ray_y = self.generate_ray(sshift)
        ray_diff.ry_origin = ray_y.o
        ray_diff.ry_direction = ray_y.d

        if (weight_x == 0.0) or (weight_y == 0.0):
            return 0.0, ray_diff

        ray_diff.has_differentials = True
        return weight, ray_diff


class ProjectiveCamera(Camera):

    """Class descriing a ProjectiveCamera."""

    def __init__(self, camera_to_world, projection, screen_window,
                 s_open, s_close, lens_radius, focal_distance,
                 film):
        """Default constructor for ProjectiveCamera."""
        super(ProjectiveCamera, self).__init__(camera_to_world,
                                               s_open, s_close,
                                               film)
        # initialize depth of field parameters
        self.lens_radius = lens_radius
        self.focal_distance = focal_distance

        # compute projective camera transformations
        self.camera_to_screen = projection

        # compute projective camera screen transformations
        sc1 = scale(float(film.x_resolution), float(film.y_resolution), 1.0)
        sc2 = scale(1.0 / (screen_window[1] - screen_window[0]),
                    1.0 / (screen_window[2] - screen_window[3]),
                    1.0)
        tr = translate(Point(-screen_window[0], -screen_window[3], 0.0))
        self.screen_to_raster = sc1 * sc2 * tr
        self.raster_to_screen = inverse(self.screen_to_raster)
        self.raster_to_camera = inverse(self.camera_to_screen) * \
                                self.raster_to_screen
