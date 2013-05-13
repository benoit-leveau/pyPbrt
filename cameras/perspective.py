"""PerspectiveCamera Class."""

from core.camera import ProjectiveCamera, Camera
from core.geometry import Point, Vector, Ray, RayDifferential, normalize
from core.transform import perspective
from core.monte_carlo import concentric_sample_disk

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
        """Generate a Ray from the camera."""
        # Generate raster and camera samples
        p_ras = Point(sample.image_x, sample.image_y, 0)
        p_camera = Point()
        p_camera = self.raster_to_camera(p_ras)
        ray = Ray(Point(0, 0, 0),
                  normalize(Vector.from_point(p_camera)),
                  0.0,
                  float('inf'))

        #  Modify ray for depth of field
        if self.lens_radius > 0.0:
            # Sample point on lens
            lens_u, lens_v = concentric_sample_disk(sample.lens_u,
                                                    sample.lens_v)
            lens_u *= self.lens_radius
            lens_v *= self.lens_radius

            # Compute point on plane of focus
            ft = self.focal_distance / ray.d.z
            p_focus = ray(ft)

            # Update ray for effect of lens
            ray.o = Point(lens_u, lens_v, 0.0)
            ray.d = normalize(p_focus - ray.o)

        ray.time = sample.time
        ray = self.camera_to_world(ray)

        return 1.0

    def generate_ray_differential(self, sample):
        """Generate a RayDifferential from the camera."""
        # Generate raster and camera samples
        p_ras = Point(sample.image_x, sample.image_y, 0)
        p_camera = Point()
        p_camera = self.raster_to_camera(p_ras)
        ray = RayDifferential(Point(0, 0, 0),
                              normalize(Vector.from_point(p_camera)),
                              0.0,
                              float('inf'))

        #  Modify ray for depth of field
        if self.lens_radius > 0.0:
            # Sample point on lens
            lens_u, lens_v = concentric_sample_disk(sample.lens_u,
                                                    sample.lens_v)
            lens_u *= self.lens_radius
            lens_v *= self.lens_radius

            # Compute point on plane of focus
            ft = self.focal_distance / ray.d.z
            p_focus = ray(ft)

            # Update ray for effect of lens
            ray.o = Point(lens_u, lens_v, 0.0)
            ray.d = normalize(p_focus - ray.o)

        ray.rx_origin = ray.o
        ray.ry_origin = ray.o
        ray.rx_direction = normalize(Vector.from_point(p_camera) + \
                                     self.dx_camera)
        ray.ry_direction = normalize(Vector.from_point(p_camera) + \
                                     self.dy_camera)
        ray.time = sample.time
        ray = self.camera_to_world(ray)
        ray.has_differentials = True
        return 1.0, ray
