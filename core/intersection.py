"""Intersection Class."""

from core.transform import Transform
from core.diffgeom import DifferentialGeometry
from core.spectrum import Spectrum


class Intersection(object):

    """Class describing an intersection of a Primitive."""

    def __init__(self):
        """Default constructor for Intersection."""
        self.dg = DifferentialGeometry()

        self.primitive = None

        self.world_to_object = Transform()
        self.object_to_world = Transform()
        
        self.shape_id = 0
        self.primitive_id = 0

        self.ray_epsilon = 0.0

    def get_bsdf(self, ray):
        """Compute the BSDF."""
        self.dg.compute_differentials(ray)
        bsdf = self.primitive.get_bsdf(dg, self.object_to_world)
        return bsdf

    def get_bssrdf(self, ray):
        """Compute the BSSRDF."""
        self.dg.compute_differentials(ray)
        bssrdf = self.primitive.get_bssrdf(dg, self.object_to_world)
        return bssrdf

    def Le(self, w):
        """Return the light emitted by the object."""
        area = self.primitive.get_area_light()
        if area:
            return area.L(self.dg.p, self.dg.nn, w)
        else:
            return Spectrum(0.0)

    def __str__(self):
        """Return a string describing the intersection."""
        return "Intersection (prim='%s', dg='%s')" % (self.primitive,
                                                      self.dg)

