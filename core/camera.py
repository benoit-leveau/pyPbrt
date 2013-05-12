"""Camera Classes."""

from abc import ABCMeta, abstractmethod

from logger import logger
from core.transform import translate, scale, inverse
from core.geometry import Vector


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
        pass

    def generate_ray_differentials(self, sample):
        raise NotImplementedError("TODO")


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
        tr = translate(Vector(-screen_window[0], -screen_window[3], 0.0))
        self.screen_to_raster = sc1 * sc2 * tr
        self.raster_to_screen = inverse(self.screen_to_raster)
        self.raster_to_camera = inverse(self.camera_to_screen) * \
                                self.raster_to_screen
