"""Shape Interface class."""

from abc import ABCMeta, abstractmethod

from core.transform import Transform
from core.geometry import abs_dot, distance_squared


class Shape(object):
    
    """Class describing a Shape."""

    __metaclass__ = ABCMeta

    next_shape_id = 1
    
    def __init__(self, object_to_world, world_to_object, reverse_orientation):
        """Constructs a Shape."""

        self.object_to_world = Transform.from_transform(object_to_world)
        self.world_to_object = Transform.from_transform(world_to_object)

        self.reverse_orientation = reverse_orientation

        self.transform_swap_handedness = object_to_world.swap_handedness()

        self.shape_id = Shape._get_shape_id()

    @classmethod
    def _get_shape_id(cls):
        shape_id = cls.next_shape_id
        cls.next_shape_id += 1
        return shape_id

    @abstractmethod
    def object_bound(self):
        """Return bounding box in object space."""
        pass

    def world_bound(self):
        """Return bounding box in world space."""
        return object_to_world(self.object_bound())

    def can_intersect(self):
        """Return True if the shape can be intersected."""
        return True

    def refine(self):
        """Refine the current shape.

        Returns an array of shape.
        """
        raise NotImplementedError()

    def intersect(self, ray):
        """Intersect the ray with the shape."""
        raise NotImplementedError()

    def intersect_p(self, ray):
        """Return True if ray intersects the shape."""
        raise NotImplementedError()

    def get_shading_geometry(self, object_to_world, dg):
        """Get the shading geometry at the intersection point."""
        return dg

    @abstractmethod
    def area(self):
        """Return the area of the shape."""
        pass

    @abstractmethod
    def sample(self, u1, u2, ns):
        """Sample the shape."""
        pass

    def pdf(self, p):
        """Compute the probability distribution function."""
        return 1.0 / self.area()

    def sample_p(self, p, u1, u2, ns):
        """Sample at point p."""
        return self.sample(u1, u2, ns)

    def pdf_wi(self, p, wi):
        """Intersect sample ray with area light geometry."""
        dg_light = DifferentialGeoemtry()
        ray = Ray(p, wi, 1e-3)
        ray.depth = -1 # temp hack to ignore alpha mask
        intersect, t_hit, ray_epsilon, dg_light - self.intersect(ray)
        if not intersect:
            return 0.0
        # convert light sample weight to solid angle measure
        pdf = distance_squared(p, ray(t_hit)) / \
              (abs_dot(dg_light.nnm -wi) * self.area())
        if pdf == float('inf'):
            return 0.0
        return pdf

    @abstractmethod
    def __str__(self):
        """Return a string describing the shape."""
        pass
