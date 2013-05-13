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
    def from_point(cls, p):
        """Construct a Vector from a Point."""
        return cls(p.x, p.y, p.z)
    
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

    def __add__(self, v):
        """Overload the addition operator."""
        return Vector(self.x + v.x,
                      self.y + v.y,
                      self.z + v.z)

    def __sub__(self, v):
        """Overload the subtraction operator."""
        return Vector(self.x - v.x,
                      self.y - v.y,
                      self.z - v.z)
    
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

    def __setitem__(self, index, value):
        """Overload the bracket operator for assignments.
            
        Example:
            v1[0] = 0.3 will set v1.x
            similarly for v1[1] and v1.y, v1[2] and v1.z
        
        """
        if index == 0:
            self.x = float(value)
        elif index == 1:
            self.y = float(value)
        elif index == 2:
            self.z = float(value)
        else:
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

    def __div__(self, f):
        """Overload the division operator."""
        inv = 1.0 / f
        return Point(self.x * inv,
                     self.y * inv,
                     self.z * inv)

    def __eq__(self, p):
        """Overload the comparison operator."""
        return eq(self.x, p.x) and eq(self.y, p.y) and eq(self.z, p.z)
    
    def __ne__(self, p):
        """Overload the comparison operator."""
        return not eq(self.x, p.x) or not eq(self.y, p.y) or not (self.z, p.z)

    def __getitem__(self, index):
        """Overload the bracket operator.
            
        Example:
            p1[0] will return p1.x
            similarly for p1[1] and p1.y, p1[2] and p1.z
        
        """
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        raise IndexError("list index out of range")

    def __setitem__(self, index, value):
        """Overload the bracket operator for assignments.
            
        Example:
            p1[0] = 0.3 will set p1.x
            similarly for p1[1] and p1.y, p1[2] and p1.z
        
        """
        if index == 0:
            self.x = float(value)
        elif index == 1:
            self.y = float(value)
        elif index == 2:
            self.z = float(value)
        else:
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

    def __setitem__(self, index, value):
        """Overload the bracket operator for assignments.
            
        Example:
            n1[0] = 0.3 will set n1.x
            similarly for n1[1] and n1.y, n1[2] and n1.z
        
        """
        if index == 0:
            self.x = float(value)
        elif index == 1:
            self.y = float(value)
        elif index == 2:
            self.z = float(value)
        else:
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
        >>> b.p_min == b.p_max

        >>> b = BBox(p1, p2)
        >>> for i in range(3):
        >>>     b.p_min[i] == min(p1[i], p2[i])
        >>>     b.p_max[i] == max(p1[i], p2[i])

        """
        if p1 and p2:
            self.p_min = Point(min(p1.x, p2.x),
                              min(p1.y, p2.y),
                              min(p1.z, p2.z))
            self.p_max = Point(max(p1.x, p2.x),
                              max(p1.y, p2.y),
                              max(p1.z, p2.z))
        elif p1:
            self.p_min = Point.from_point(p1)
            self.p_max = Point.from_point(p1)
        elif p2:
            self.p_min = Point.from_point(p2)
            self.p_max = Point.from_point(p2)
        else:
            self.p_min = Point(float('inf'), float('inf'), float('inf'))
            self.p_max = Point(-float('inf'), -float('inf'), -float('inf'))

    @classmethod
    def from_bbox(cls, b):
        """Copy constructor."""
        return cls(b.p_min, b.p_max)

    def overlaps(self, b):
        """Return True if overlaps with specified bounding box."""
        x = (self.p_max.x >= b.p_min.x) and (self.p_min.x <= b.p_max.x)
        y = (self.p_max.y >= b.p_min.y) and (self.p_min.y <= b.p_max.y)
        z = (self.p_max.z >= b.p_min.z) and (self.p_min.z <= b.p_max.z)
        return x and y and z

    def inside(self, pt):
        """Return True if point is inside the box."""
        x = (pt.x >= self.p_min.x) and (pt.x <= self.p_max.x)
        y = (pt.y >= self.p_min.y) and (pt.y <= self.p_max.y)
        z = (pt.z >= self.p_min.z) and (pt.z <= self.p_max.z)
        return x and y and z

    def expand(self, delta):
        """Pad the bounding box by a constant factor."""
        self.p_min -= Vector(delta, delta, delta)
        self.p_max += Vector(delta, delta, delta)

    def surface_area(self):
        """Compute the surface area of the six faces of the box."""
        d = self.p_max - self.p_min
        return 2.0 * (d.x*d.y + d.x*d.z + d.y*d.z)
    
    def volume(self):
        """Compute the volume of the box."""
        d = self.p_max - self.p_min
        return d.x*d.y*d.z

    def maximum_extent(self):
        """Compute the volume of the box."""
        diag = self.p_max - self.p_min
        if diag.x>diag.y and diag.y>diag.z:
            return 0
        elif diag.y>diag.z:
            return 1
        else:
            return 2
    
    def lerp(self, tx, ty, tz):
        """Linear Interpolation between p_min and p_max."""
        return Point(lerp(tx, self.p_min.x, self.p_max.x),
                     lerp(ty, self.p_min.y, self.p_max.y))
                
    def offset(self, p):
        """Return position of a point relative to the corners."""
        return Vector((p.x - self.p_min.x) / (self.p_max.x - self.p_min.x),
                      (p.y - self.p_min.y) / (self.p_max.y - self.p_min.y),
                      (p.z - self.p_min.z) / (self.p_max.z - self.p_min.z))

    def bounding_sphere(self):
        """Return a sphere containing the entire bounding box."""
        center = 0.5 * (self.p_min+self.p_max)
        if self.inside(center):
            radius = distance(center, self.p_max)
        else:
            radius = 0.0
        return center, radius

    def intersect_p(self, ray):
        """Compute intersections with the bounding box."""
        t0 = ray.mint
        t1 = ray.maxt
        for i in range(3):
            # Update interval for _i_th bounding box slab
            if ray.d[i] == 0.0:
                inv_ray_dir = float('inf')
            else:
                inv_ray_dir = 1.0 / ray.d[i]
            t_near = (self.p_min[i] - ray.o[i]) * inv_ray_dir
            t_far  = (self.p_max[i] - ray.o[i]) * inv_ray_dir

            # Update parametric interval from slab intersection $t$s
            if (t_near > t_far):
                t_near, t_far = t_far, t_near
            if t_near > t0:
                t0 = t_near
            if t_far < t1:
                t1 = t_far
            if t0 > t1:
                return False, float('inf'), float('inf')
        return True, t0, t1

    def __getitem__(self, index):
        """Overload the bracket operator.
        
        Example:
        bbox[0] will return bbox.p_min
        bbox[1] will return bbox.p_max
        
        """
        if index == 0:
            return self.p_min
        elif index == 1:
            return self.p_max
        raise IndexError("list index out of range")

    def __str__(self):
        """Return a string describing the bbox."""
        return "BBox (min='%s', max='%s')" % (str(self.p_min), str(self.p_max))


def union(b, elt):
    """Return the union of a BBox and a Point/BBox."""
    if isinstance(elt, Point):
        ret = BBox()
        ret.p_min.x = min(b.p_min.x, elt.x)
        ret.p_min.y = min(b.p_min.y, elt.y)
        ret.p_min.z = min(b.p_min.z, elt.z)
        ret.p_max.x = max(b.p_max.x, elt.x)
        ret.p_max.y = max(b.p_max.y, elt.y)
        ret.p_max.z = max(b.p_max.z, elt.z)
        return ret
    elif isinstance(elt, BBox):
        ret = BBox()
        ret.p_min.x = min(b.p_min.x, elt.p_min.x)
        ret.p_min.y = min(b.p_min.y, elt.p_min.y)
        ret.p_min.z = min(b.p_min.z, elt.p_min.z)
        ret.p_max.x = max(b.p_max.x, elt.p_max.x)
        ret.p_max.y = max(b.p_max.y, elt.p_max.y)
        ret.p_max.z = max(b.p_max.z, elt.p_max.z)
        return ret
    raise TypeError("second argument must be a Point or BBox")

