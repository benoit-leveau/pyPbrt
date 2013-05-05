"""Geometry classes and utility functions."""

import math


class Vector(object):
    
    """Class holding a 3D vector."""
    
    def __init__(self, x=0.0, y=0.0, z=0.0):
        """Construct a Vector instance."""
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)

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
        return '(' + str(self.x) + ',' + str(self.y) + ',' + str(self.z) + ')'


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

    def __str__(self):
        """Return a string describing the point."""
        return '(' + str(self.x) + ',' + str(self.y) + ',' + str(self.z) + ')'
 

def distance(p1, p2):
    """Return the distance between two points."""
    return (p1-p2).length()


def distance_squared(p1, p2):
    """Return the squared distance between two points."""
    return (p1-p2).length_squared()

