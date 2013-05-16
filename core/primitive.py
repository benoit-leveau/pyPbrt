"""Primitive Classes."""

from abc import ABCMeta, abstractmethod
from logger import logger

from core.transform import inverse, normalize


class Primitive(object):

    """Primitive Class."""

    __metaclass__ = ABCMeta

    next_primitive_id = 1
    
    def __init__(self):
        """Default constructor for Primitive."""
        self.primitive_id = Primitive._get_shape_id()

    @classmethod
    def _get_shape_id(cls):
        """Return a new unique primitive id."""
        primitive_id = cls.next_primitive_id
        cls.next_primitive_id += 1
        return primitive_id

    @abstractmethod
    def world_bound(self):
        """Return the bounding box of the primitive, in world space."""

    def can_intersect(self):
        """Return True if the primitive allow intersection calls."""
        return True

    @abstractmethod
    def intersect(self, ray, intersection):
        """Compute an intersection."""
        pass

    @abstractmethod
    def intersect_p(self, ray):
        """Return True if the submitted ray intersects the primitive."""
        pass

    def refine(self, refined):
        """Fill the list with refined shapes."""
        raise NotImplementedError()

    def fully_refine(self, refined):
        """Fill the list with all refined shapes, recursively."""
        to_process = []
        to_process.append(self)
        while(len(to_process)):
            # refine last primitive in to_process list
            primitive = to_process.pop()
            if primitive.can_intersect():
                refined.append(primitive)
            else:
                primitive.refine(to_process)

    @abstractmethod
    def get_area_light(self):
        """Return the AreaLight corresponding to the primitive."""
        pass
    
    @abstractmethod
    def get_bsdf(self, dg, object_to_world):
        """Compute the BSDF."""
        pass

    @abstractmethod
    def get_bssrdf(self, dg, object_to_world):
        """Compute the BSSRDF."""
        pass

    @abstractmethod
    def __str__(self):
        """Return a string describing the primitive."""
        pass


class GeometricPrimitive(Primitive):

    """GeometricPrimitive Class."""

    def __init__(self, shape, material, area_light=None):
        super(GeometricPrimitive, self).__init__()
        self.shape = shape
        self.material = material
        self.area_light = area_light

    def world_bound(self):
        """Return the bounding box of the primitive, in world space."""
        return self.shape.world_bound()

    def can_intersect(self):
        """Return True if the primitive allow intersection calls."""
        return self.shape.can_intersect()

    def refine(self, refined):
        """Fill the list with refined shapes."""
        shapes = []
        self.shape.refine(shapes)
        for shape in shapes:
            refined.append(GeometricPrimitive(shape,
                                              self.material,
                                              self.area_light))

    def intersect(self, ray, intersection):
        """Compute an intersection."""
        intersected, t_hit, ray_epsilon, dg = self.shape.intersect(ray)
        if not intersected:
            return False
        intersection.primitive = self
        intersection.world_to_object = self.shape.world_to_object
        intersection.object_to_world = self.shape.object_to_world
        intersection.shape_id = self.shape.shape_id
        intersection.primitive_id = self.primitive_id
        intersection.ray_epsilon = ray_epsilon
        ray.maxt = t_hit
        return True

    def intersect_p(self, ray):
        """Return True if the submitted ray intersects the primitive."""
        return self.shape.intersect_p(ray)

    def get_area_light(self):
        """Return the AreaLight corresponding to the primitive."""
        return self.area_light
    
    def get_bsdf(self, dg, object_to_world):
        """Compute the BSDF."""
        dgs = self.shape.get_shading_geometry(object_to_world, dg)
        return self.material.get_bsdf(dg, dgs)

    def get_bssrdf(self, dg, object_to_world):
        """Compute the BSSRDF."""
        dgs = self.shape.get_shading_geometry(object_to_world, dg)
        return self.material.get_bssrdf(dg, dgs)

    def __str__(self):
        """Return a string describing the geometric primitive."""
        return "GeometricPrimitive(shape='%s', material='%s', light='%s')" % \
               (self.shape, self.material, self.area_light)


class TransformedPrimitive(Primitive):

    """Class describing a TransformedPrimitive."""

    def __init__(self, primitive, world_to_primitive):
        """Default constructor for TransformedPrimitive."""
        super(TransformedPrimitive, self).__init__()
        self.primitive = primitive
        self.world_to_primitive = Transform.from_transform(world_to_primitive)

    def world_bound(self):
        """Return the bounding box of the primitive, in world space."""
        return self.world_to_primitive.motion_bounds(
            self.primitive.world_bound(), True)

    def intersect(self, ray, intersection):
        """Compute an intersection."""
        w2p = self.world_to_primitive.interpolate(ray.time)
        ray_primitive = w2p(ray)
        intersection = Intersection()
        found_intersect, ray_hit, ray_epsilon, dg = self.primitive.intersect(
            ray_primitive, intersection)
        if not found_intersect:
            return False
        ray.maxt = ray_primitive.maxt
        intersection.primitive_id = self.primitive_id
        if not w2p.is_identity():
            # Compute world-to-object transformation for instance
            intersection.world_to_object = intersection.world_to_object * w2p
            intersection.object_to_world = inverse(intersection.world_to_object)

            # Transform instance's differential geometry to world space
            p2w = inverse(w2p)
            intersection.dg.p = p2w(intersection.dg.p)
            intersection.dg.nn = normalize(p2w(intersection.dg.nn))
            intersection.dg.dpdu = p2w(intersection.dg.dpdu)
            intersection.dg.dpdv = p2w(intersection.dg.dpdv)
            intersection.dg.dndu = p2w(intersection.dg.dndu)
            intersection.dg.dndv = p2w(intersection.dg.dndv)
        return True

    def intersect_p(self, ray):
        """Return True if the submitted ray intersects the primitive."""
        return self.primitive.intersect_p(self.world_to_primitive(ray))

    def get_area_light(self):
        """Return the AreaLight corresponding to the primitive."""
        return None
    
    def get_bsdf(self, dg, object_to_world):
        """Compute the BSDF."""
        return None

    def get_bssrdf(self, dg, object_to_world):
        """Compute the BSSRDF."""
        return None

    def __str__(self):
        """Return a string describing the primitive."""
        return "TransformedPrimitive (prim='%s')" % self.primitive
    

class Aggregate(Primitive):

    """Class describing an Aggregate."""

    def get_area_light(self):
        """Return the AreaLight corresponding to the primitive."""
        logger.error("Aggregate.get_area_light() method called;"\
                     "should have gone to GeometricPrimitive.")
        return None
    
    def get_bsdf(self, dg, object_to_world):
        """Compute the BSDF."""
        logger.error("Aggregate.get_bsdf() method called;"\
                     "should have gone to GeometricPrimitive.")
        return None

    def get_bssrdf(self, dg, object_to_world):
        """Compute the BSSRDF."""
        logger.error("Aggregate.get_bssrdf() method called;"\
                     "should have gone to GeometricPrimitive.")
        return None

    def __str__(self):
        """Return a string describing the Aggregate."""
        return "Aggregate ()"
