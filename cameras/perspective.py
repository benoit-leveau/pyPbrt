"""PerspectiveCamera Class."""

from core.camera import ProjectiveCamera
from core.geometry import Point
from core.transform import perspective


class PerspectiveCamera(ProjectiveCamera):

    """Class describing a Perspective Camera."""

    def __init__(self, camera_to_world, screen_window, s_open, s_close,
                 lens_radius, focal_distance, fov, film):
        
        super(PerspectiveCamera, self).__init__(
            camera_to_world,
            perspective(fov, 1e-2, 1000.0),
            screen_window,
            s_open, s_close,
            lens_radius, focal_distance,
            film)
        self.dx_camera = self.raster_to_camera(Point(1, 0, 0)) - \
                         self.raster_to_camera(Point(0, 0, 0))
        self.dy_camera = self.raster_to_camera(Point(0, 1, 0)) - \
                         self.raster_to_camera(Point(0, 0, 0))

    def generate_ray(self, sample):
        """Generate a Ray."""
        pass

    def generate_ray_differential(self, sample):
        """Generate a RayDifferential."""
        pass
