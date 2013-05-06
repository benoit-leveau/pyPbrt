"""Scene Class"""

from core.geometry import union

class Scene(object):

    """Class describing a pyPbrt scene."""

    def __init__(self, accel, lts=None, vr=None):
        """Default constructor for Scene."""

        # Primitive
        self.aggregate = accel

        # array of Light
        if lts:
            self.lights = lts
        else:
            self.lights = []

        # VolumeRegion
        self.volume_region = vr

        # initialize bounding box
        self.bound = self.aggregate.world_bound()
        if volumeRegion:
            self.bound = union(self.bound, volume_region.world_bound())

    def intersect(self, ray):
        """Intersect a ray with the scene."""
        hit, isect = aggregate.intersect(ray)
        return hit, isect

    def intersect_p(self, ray):
        """More efficient intersection function.

        Search for intersection is stopped as soon as one is found.

        """
        hit = aggregate.intersect_p(ray)
        return hit

    def world_bound(self):
        """Return the world bounding box of the scene."""
        return bound
