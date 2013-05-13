"""Scene Class"""

from core.geometry import union

class Scene(object):

    """Class describing a pyPbrt scene."""

    def __init__(self, accel, lights=None, volume_region=None):
        """Default constructor for Scene."""

        # Primitive
        self.aggregate = accel

        # array of Light
        if lights:
            self.lights = lights
        else:
            self.lights = []

        # VolumeRegion
        self.volume_region = volume_region

        # initialize bounding box
        self.bound = self.aggregate.world_bound()
        if self.volume_region:
            self.bound = union(self.bound, self.volume_region.world_bound())

    def intersect(self, ray, intersection):
        """Intersect a ray with the scene."""
        hit = self.aggregate.intersect(ray, intersection)
        return hit

    def intersect_p(self, ray):
        """More efficient intersection function.

        Search for intersection is stopped as soon as one is found.

        """
        hit = aggregate.intersect_p(ray)
        return hit

    def world_bound(self):
        """Return the world bounding box of the scene."""
        return bound
