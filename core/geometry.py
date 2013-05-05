"""Geometry classes and utility functions."""

import math


class Vector(object):
    
    """Class holding a 3D vector."""
    
    def __init__(self, x=0.0, y=0.0, z=0.0):
        """Construct a Vector instance."""
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

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
        return self.x == v.x and self.y == v.y and self.z == v.z

    def __ne__(self, v):
        """Overload the comparison operator."""
        return self.x != v.x or self.y != v.y or self.z != v.z

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
        else:
            raise IndexError("list index out of range")

    def __str__(self):
        """Return a string describing the vector."""
        return "Vector " + '(' + str(self.x) + ',' + str(self.y) + ',' + str(self.z) + ')'


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
        return self.x == p.x and self.y == p.y and self.z == p.z
    
    def __ne__(self, p):
        """Overload the comparison operator."""
        return self.x != p.x or self.y != p.y or self.z != p.z
    
    def __str__(self):
        """Return a string describing the point."""
        return "Point " + '(' + str(self.x) + ',' + str(self.y) + ',' + str(self.z) + ')'
 

def distance(p1, p2):
    """Return the distance between two points."""
    return (p1-p2).length()


def distance_squared(p1, p2):
    """Return the squared distance between two points."""
    return (p1-p2).length_squared()


class Normal(object):
    
    """Class holding a 3D Normal."""
    
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
        return self.x == n.x and self.y == n.y and self.z == n.z
    
    def __ne__(self, n):
        """Overload the comparison operator."""
        return self.x != n.x or self.y != n.y or self.z != n.z
    
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
        else:
            raise IndexError("list index out of range")
    
    def __str__(self):
        """Return a string describing the vector."""
        return "Normal " + '(' + str(self.x) + ',' + str(self.y) + ',' + str(self.z) + ')'


def face_forward(n, v):
    """Return the supplied normal flipped to the same hemisphere as v.""" 
    if dot(n, v) < 0.0:
        return -n
    else:
        return n


class Ray(object):

    """Class hodling a 3D Ray."""

    def __init__(self, origin=None, direction=None, start=0.0, end=float('inf'), time=0.0, depth=0):
        """Constructor for a 3D Ray."""
        if origin:
            self.o = origin
        else:
            self.o = Point(0,0,0)
        if direction:
            self.d = direction
        else:
            self.d = Vector(0,0,0)

        # the following are used to constrain the ray to a segment
        self.mint = start
        self.maxt = end
        
        # time at which this ray is evaluated
        self.time = time
        
        # number of bounces this ray went through
        self.depth = depth
