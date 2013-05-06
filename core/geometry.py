"""Geometry classes and utility functions."""

import math

from core.pbrt import lerp, eq


class Vector(object):
    
    """Class describing a 3D vector."""
    
    def __init__(self, x=0.0, y=0.0, z=0.0):
        """Construct a Vector instance."""
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    @classmethod
    def from_vector(cls, v):
        """Copy constructor."""
        return cls(v.x, v.y, v.z)
    
    @classmethod
    def from_normal(cls, n):
        """Construct a Vector from a Normal."""
        return cls(n.x, n.y, n.z)
    
    def length_squared(self):
        """Return the square value of the length of the vector."""
        return self.x*self.x + self.y*self.y + self.z*self.z
    
    def length(self):
        """Return the length of the vector."""
        return math.sqrt(self.length_squared())
        
    def normalize(self):
        """Returns a normalized (unit length) vector."""
        return self/self.length()
    
    def __neg__(self):
        """Overload the negation operator."""
        return Vector(-self.x,
                      -self.y,
                      -self.z)
    
    def __mul__(self, f):
        """Overload the multiplication operator."""
        return Vector(self.x * f,
                      self.y * f,
                      self.z * f)
    
    def __rmul__(self, f):
        """Overload the right multiplication operator."""
        return Vector(self.x * f,
                      self.y * f,
                      self.z * f)

    def __div__(self, f):
        """Overload the division operator."""
        inv = 1.0 / f
        return Vector(self.x * inv,
                      self.y * inv,
                      self.z * inv)

    def __eq__(self, v):
        """Overload the comparison operator."""
        return eq(self.x, v.x) and eq(self.y, v.y) and eq(self.z, v.z)

    def __ne__(self, v):
        """Overload the comparison operator."""
        return not eq(self.x, v.x) or not eq(self.y, v.y) or not (self.z, v.z)

    def __getitem__(self, index):
        """Overload the bracket operator.
            
        Example:
            v1[0] will return v1.x
            similarly for v1[1] and v1.y, v1[2] and v1.z
        
        """
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        raise IndexError("list index out of range")

    def __str__(self):
        """Return a string describing the vector."""
        return "Vector (%f, %f, %f)" % (self.x, self.y, self.z)


def dot(v1, v2):
    """Return the dot (scalar) product of two vectors.
    
    Formula:
        dot(v1, v2) = (v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z)
    
    """
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z


def abs_dot(v1, v2):
    """Return the absolute value of the dot (scalar) product of two vectors.
        
    Formula:
        abs_dot(v1, v2) = abs(dot(v1, v2))
    
    """
    return abs(v1.x * v2.x + v1.y * v2.y + v1.z * v2.z)


def cross(v1, v2):
    """Return the cross product between two vectors, as a new Vector.
    
    Formula:
        cross(v1, v2) = (v1.y * v2.z) - (v1.z * v2.y)
                        (v1.z * v2.x) - (v1.x * v2.z)
                        (v1.x * v2.y) - (v1.y * v2.x)

        
    """
    return Vector((v1.y * v2.z) - (v1.z * v2.y),
                  (v1.z * v2.x) - (v1.x * v2.z),
                  (v1.x * v2.y) - (v1.y * v2.x))


def normalize(v):
    """Return a the normalized vector of v."""
    return v/v.length()


def coordinate_system(v):
    """Return a coordinate system (defined by 3 vectors) based on the supplied vector."""
    v1 = Vector(v.x, v.y, v.z)
    if abs(v1.x) > abs(v1.y):
        invLen = 1.0 / math.sqrt(v1.x*v1.x + v1.z*v1.z)
        v2 = Vector(-v1.z * invLen, 0.0, v1.x * invLen)
    else:
        invLen = 1.0 / math.sqrt(v1.y*v1.y + v1.z*v1.z)
        v2 = Vector(-v1.z * invLen, 0.0, v1.x * invLen)
    v3 = cross(v1, v2)
    return v1, v2, v3


class Point(object):

    """3D Point"""

    def __init__(self, x=0.0, y=0.0, z=0.0):
        """Construct a Point instance."""
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

    @classmethod
    def from_point(cls, p):
        """Copy constructor."""
        return cls(p.x, p.y, p.z)

    def __add__(self, v_or_p):
        """Return a point translated by the supplied vector, or a sum of two points."""
        return Point(self.x+v_or_p.x,
                     self.y+v_or_p.y,
                     self.z+v_or_p.z)

    def __sub__(self, v_or_p):
        """Return a point translated by the supplied vector, or a vector if supplied with a point."""
        if isinstance(v_or_p, Point):
            return Vector(self.x-v_or_p.x,
                          self.y-v_or_p.y,
                          self.z-v_or_p.z)
        elif isinstance(v_or_p, Vector):
            return Point(self.x-v_or_p.x,
                         self.y-v_or_p.y,
                         self.z-v_or_p.z)

    def __mul__(self, f):
        """Overload the multiplication operator."""
        return Point(self.x * f,
                     self.y * f,
                     self.z * f)
    
    def __rmul__(self, f):
        """Overload the right multiplication operator."""
        return Point(self.x * f,
                     self.y * f,
                     self.z * f)

    def __eq__(self, p):
        """Overload the comparison operator."""
        return eq(self.x, p.x) and eq(self.y, p.y) and eq(self.z, p.z)
    
    def __ne__(self, p):
        """Overload the comparison operator."""
        return not eq(self.x, p.x) or not eq(self.y, p.y) or not (self.z, p.z)

    def __getitem__(self, index):
        """Overload the bracket operator.
            
        Example:
            v1[0] will return v1.x
            similarly for v1[1] and v1.y, v1[2] and v1.z
        
        """
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        raise IndexError("list index out of range")
    
    def __str__(self):
        """Return a string describing the point."""
        return "Point (%f, %f, %f)" % (self.x, self.y, self.z)
 

def distance(p1, p2):
    """Return the distance between two points."""
    return (p1-p2).length()


def distance_squared(p1, p2):
    """Return the squared distance between two points."""
    return (p1-p2).length_squared()


class Normal(object):
    
    """Class describing a 3D Normal."""
    
    def __init__(self, x=0.0, y=0.0, z=0.0):
        """Construct a Normal instance."""
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
    
    @classmethod
    def from_normal(cls, n):
        """Construct a Normal from another Normal."""
        return cls(n.x, n.y, n.z)

    @classmethod
    def from_vector(cls, v):
        """Construct a Normal from a Vector."""
        return cls(v.x, v.y, v.z)

    def length_squared(self):
        """Return the square value of the length of the normal."""
        return self.x*self.x + self.y*self.y + self.z*self.z
    
    def length(self):
        """Return the length of the normal."""
        return math.sqrt(self.length_squared())
    
    def normalize(self):
        """Returns a normalized (unit length) normal."""
        return self/self.length()
    
    def __neg__(self):
        """Overload the negation operator."""
        return Normal(-self.x,
                      -self.y,
                      -self.z)
    
    def __mul__(self, f):
        """Overload the multiplication operator."""
        return Normal(self.x * f,
                      self.y * f,
                      self.z * f)
    
    def __rmul__(self, f):
        """Overload the right multiplication operator."""
        return Normal(self.x * f,
                      self.y * f,
                      self.z * f)
    
    def __div__(self, f):
        """Overload the division operator."""
        inv = 1.0 / f
        return Normal(self.x * inv,
                      self.y * inv,
                      self.z * inv)
    
    def __eq__(self, n):
        """Overload the comparison operator."""
        return eq(self.x, n.x) and eq(self.y, n.y) and eq(self.z, n.z)
    
    def __ne__(self, n):
        """Overload the comparison operator."""
        return not eq(self.x, n.x) or not eq(self.y, n.y) or not (self.z, n.z)
    
    def __getitem__(self, index):
        """Overload the bracket operator.
            
            Example:
            n1[0] will return n1.x
            similarly for n1[1] and n1.y, n1[2] and n1.z
            
            """
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        raise IndexError("list index out of range")
    
    def __str__(self):
        """Return a string describing the normal."""
        return "Normal (%f, %f, %f)" % (self.x, self.y, self.z)

def face_forward(n, v):
    """Return the supplied normal flipped to the same hemisphere as v.""" 
    if dot(n, v) < 0.0:
        return -n
    else:
        return n


class Ray(object):

    """Class describing a 3D Ray."""

    def __init__(self, origin=None, direction=None, start=0.0, end=float('inf'), time=0.0, depth=0):
        """Constructor for a 3D Ray."""
        if origin:
            self.o = Point.from_point(origin)
        else:
            self.o = Point()
        if direction:
            self.d = Vector.from_vector(direction)
        else:
            self.d = Vector()

        # the following are used to constrain the ray to a segment
        self.mint = start
        self.maxt = end
        
        # time at which this ray is evaluated
        self.time = time
        
        # number of bounces this ray went through
        self.depth = depth

    @classmethod
    def from_ray(cls, ray):
        """Copy constructor for Ray."""
        return cls(origin=ray.o,
                   direction=ray.d,
                   start=ray.mint,
                   end=ray.maxt,
                   time=ray.time,
                   depth=ray.depth
                   )
        
    @classmethod
    def from_ray_parent(cls, origin, direction, parent, start, end=float('inf')):
        """Construct a (spawned) Ray from a parent Ray."""
        return cls(origin=origin,
                   direction=direction,
                   start=start,
                   end=end,
                   time=parent.time,
                   depth=parent.depth+1 # increment depth
                   )

    def __call__(self, t):
        """Override the operator().

        >>> r = Ray(Point(0, 0, 0), Vector(1, 2, 3))
        >>> p = r(1.7)

        """
        return self.o + self.d * t

    def __str__(self):
        """Return a string describing the ray."""
        return "Ray ( %s, %s, [%f-%f], %f, %d )" % (
            self.o, self.d,
            self.mint, self.maxt,
            self.time, self.depth
            )


class RayDifferential(Ray):

    """Class describing a 3D Ray Differential (adding dx & dy rays to Ray)."""

    def __init__(self, origin=None, direction=None, start=0.0, end=float('inf'), time=0.0, depth=0):
        """Constructor for a 3D RayDifferential."""
        super(RayDifferential, self).__init__(origin, direction, start, end, time, depth)
        self.has_differentials = False
        self.rx_origin = Point()
        self.ry_origin = Point()
        self.rx_direction = Vector()
        self.ry_direction = Vector()

    @classmethod
    def from_ray_differential(cls, ray):
        """Copy constructor."""
        ret = cls(origin=ray.o,
                  direction=ray.d,
                  start=ray.mint,
                  end=ray.maxt,
                  time=ray.time,
                  depth=ray.depth)
        ret.has_differentials = ray.has_differentials
        ret.rx_origin = ray.rx_origin
        ret.ry_origin = ray.ry_origin
        ret.rx_direction = ray.rx_direction
        ret.ry_direction = ray.ry_direction
        return ret

    @classmethod
    def from_ray_parent(cls, origin, direction, parent, start, end=float('inf')):
        """Construct a (spawned) Ray from a parent Ray."""
        return cls(origin=origin,
                   direction=direction,
                   start=start,
                   end=end,
                   time=parent.time,
                   depth=parent.depth+1)

    @classmethod
    def from_ray(cls, ray):
        """Construct a RayDifferential from a Ray."""
        return cls(origin=ray.o,
                   direction=ray.d,
                   start=ray.mint,
                   end=ray.maxt,
                   time=ray.time,
                   depth=ray.depth)
    
    def scale_differentials(self, s):
        """Scale the differentials rays to accomodate for different ray spacings."""
        self.rx_origin = self.o + (self.rx_origin - self.o) * s
        self.ry_origin = self.o + (self.ry_origin - self.o) * s
        self.rx_direction = self.d + (self.rx_direction - self.d) * s
        self.ry_direction = self.d + (self.ry_direction - self.d) * s
        
    def __str__(self):
        """Return a string describing the ray differential."""
        return "RayDifferential ( %s, %s, [%f-%f], %f, %d, %s)" % (
            self.o, self.d,
            self.mint, self.maxt,
            self.time, self.depth,
            str(self.has_differentials)
            )


class BBox(object):

    """Class for Axis-Aligned Bounding Box (AABB)."""

    def __init__(self, p1=None, p2=None):
        """Construct a BBox with optional points.

        >>> b = BBox(p)
        >>> b.pMin == b.pMax

        >>> b = BBox(p1, p2)
        >>> for i in range(3):
        >>>     b.pMin[i] == min(p1[i], p2[i])
        >>>     b.pMax[i] == max(p1[i], p2[i])

        """
        if p1 and p2:
            self.pMin = Point(min(p1.x, p2.x),
                              min(p1.y, p2.y),
                              min(p1.z, p2.z))
            self.pMax = Point(max(p1.x, p2.x),
                              max(p1.y, p2.y),
                              max(p1.z, p2.z))
        elif p1:
            self.pMin = Point.from_point(p1)
            self.pMax = Point.from_point(p1)
        elif p2:
            self.pMin = Point.from_point(p2)
            self.pMax = Point.from_point(p2)
        else:
            self.pMin = Point(float('inf'), float('inf'), float('inf'))
            self.pMax = Point(-float('inf'), -float('inf'), -float('inf'))

    @classmethod
    def from_bbox(cls, b):
        """Copy constructor."""
        return cls(b.pMin, b.pMax)

    def overlaps(self, b):
        """Return True if overlaps with specified bounding box."""
        x = (self.pMax.x >= b.pMin.x) and (self.pMin.x <= b.pMax.x)
        y = (self.pMax.y >= b.pMin.y) and (self.pMin.y <= b.pMax.y)
        z = (self.pMax.z >= b.pMin.z) and (self.pMin.z <= b.pMax.z)
        return x and y and z

    def inside(self, pt):
        """Return True if point is inside the box."""
        x = (pt.x >= self.pMin.x) and (pt.x <= self.pMax.x)
        y = (pt.y >= self.pMin.y) and (pt.y <= self.pMax.y)
        z = (pt.z >= self.pMin.z) and (pt.z <= self.pMax.z)
        return x and y and z

    def expand(self, delta):
        """Pad the bounding box by a constant factor."""
        self.pMin -= Vector(delta, delta, delta)
        self.pMax += Vector(delta, delta, delta)

    def surface_area(self):
        """Compute the surface area of the six faces of the box."""
        d = self.pMax - self.pMin
        return 2.0 * (d.x*d.y + d.x*d.z + d.y*d.z)
    
    def volume(self):
        """Compute the volume of the box."""
        d = self.pMax - self.pMin
        return d.x*d.y*d.z

    def maximum_extent(self):
        """Compute the volume of the box."""
        diag = self.pMax - self.pMin
        if diag.x>diag.y and diag.y>diag.z:
            return 0
        elif diag.y>diag.z:
            return 1
        else:
            return 2
    
    def lerp(self, tx, ty, tz):
        """Linear Interpolation between pMin and pMax."""
        return Point(lerp(tx, self.pMin.x, self.pMax.x),
                     lerp(ty, self.pMin.y, self.pMax.y))
                
    def offset(self, p):
        """Return position of a point relative to the corners."""
        return Vector((p.x - self.pMin.x) / (self.pMax.x - self.pMin.x),
                      (p.y - self.pMin.y) / (self.pMax.y - self.pMin.y),
                      (p.z - self.pMin.z) / (self.pMax.z - self.pMin.z))

    def bounding_sphere(self):
        """Return a sphere containing the entire bounding box."""
        center = 0.5 * (self.pMin+self.pMax)
        if self.inside(center):
            radius = distance(center, self.pMax)
        else:
            radius = 0.0
        return center, radius

    def __getitem__(self, index):
        """Overload the bracket operator.
        
        Example:
        bbox[0] will return bbox.pMin
        bbox[1] will return bbox.pMax
        
        """
        if index == 0:
            return self.pMin
        elif index == 1:
            return self.pMax
        raise IndexError("list index out of range")

    def __str__(self):
        """Return a string describing the bbox."""
        return "BBox (min='%s', max='%s')" % (str(self.pMin), str(self.pMax))


def union(b, b_or_p):
    """Return the union of a BBox and a Point/BBox."""
    ret = BBox()
    if isinstance(b_or_p, Point):
        ret.pMin.x = min(b.pMin.x, b_or_p.x)
        ret.pMin.y = min(b.pMin.y, b_or_p.y)
        ret.pMin.z = min(b.pMin.z, b_or_p.z)
        ret.pMax.x = max(b.pMax.x, b_or_p.x)
        ret.pMax.y = max(b.pMax.y, b_or_p.y)
        ret.pMax.z = max(b.pMax.z, b_or_p.z)
    else:
        ret.pMin.x = min(b.pMin.x, b_or_p.pMin.x)
        ret.pMin.y = min(b.pMin.y, b_or_p.pMin.y)
        ret.pMin.z = min(b.pMin.z, b_or_p.pMin.z)
        ret.pMax.x = max(b.pMax.x, b_or_p.pMax.x)
        ret.pMax.y = max(b.pMax.y, b_or_p.pMax.y)
        ret.pMax.z = max(b.pMax.z, b_or_p.pMax.z)
    return ret
