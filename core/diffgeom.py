"""DifferentialGeometry class."""

from core.geometry import Point, Vector, Normal
from core.transform import normalize, cross

class DifferentialGeometry(object):
    
    """Class describing a DifferentialGeoemtry for intersections)."""

    def __init__(self):
        """Default constructor for DifferentialGeometry."""
        self.p = Point()
        self.nn = Normal()
        
        self.u = 0.0
        self.v = 0.0

        self.shape = None

        self.dp_du = Vector()
        self.dp_dv = Vector()

        self.dn_du = Normal()
        self.dn_dv = Normal()

        self.dp_dx = Vector()
        self.dp_dy = Vector()
        
        self.du_dx = 0.0
        self.dv_dx = 0.0

        self.du_dy = 0.0
        self.dv_dy = 0.0

    @classmethod
    def from_intersection(cls, p, dp_du, dp_dv, dn_du, dn_dv, uu, vv, shape=None):
        """Construct a DifferentialGeometry instance.

        Arguments:
        p : Point
        dp_du : Vector
        dp_dv : Vector
        dn_du : Normal
        dn_dv : Normal
        uu : float
        vv : float
        shape : optional Shape

        """
        # call default constructor
        diff_geom = cls()

        diff_geom.p = p
        diff_geom.dp_du = dp_du
        diff_geom.dp_dv = dp_dv
        diff_geom.dn_du = dn_du
        diff_geom.dn_dv = dn_dv

        diff_geom.nn = Normal.from_vector(normalize(cross(dp_du, dp_dv)))
        diff_geom.u = uu
        diff_geom.v = vv
        diff_geom.shape = shape

        if shape and \
           (shape.reverse_orientation ^ shape.transform_swap_handedness):
            diff_geom.nn *= -1.0

        return diff_geom
    
    def compute_differentials(self, raydiff):
        """Computes the differentials from a RayDifferential."""
        return cls()

    def __str__(self):
        """Return a string describing the diff geometry."""
        return "DifferentialGeometry (%s, %s)" % (self.p, self.nn)
