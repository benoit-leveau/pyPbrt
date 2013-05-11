"""Quaternion Class."""

import math

from core.pbrt import clamp, eq
from core.geometry import Vector, dot
from core.transform import Matrix4x4, Transform, transpose

class Quaternion(object):

    """Class describing a Quaternion."""

    def __init__(self, v=None, w=1.0):
        """Default constructor for Quaternion."""
        if v is None:
            self.v = Vector(0.0, 0.0, 0.0)
        else:
            self.v = Vector.from_vector(v)
        self.w = w
        
    @classmethod
    def from_transform(cls, transform):
        """Constructor from a Transform."""
        m = transform.m
        trace = m.m[0][0] + m.m[1][1] + m.m[2][2]
        v = Vector()
        w = 1.0
        if trace > 0.0:
            # Compute w from matrix trace, then xyz
            # 4w^2 = m[0][0] + m[1][1] + m[2][2] + m[3][3] (but m[3][3] == 1)
            s = math.sqrt(trace + 1.0)
            w = s / 2.0
            s = 0.5 / s
            v.x = (m.m[2][1] - m.m[1][2]) * s
            v.y = (m.m[0][2] - m.m[2][0]) * s
            v.z = (m.m[1][0] - m.m[0][1]) * s
        else:
            # Compute largest of $x$, $y$, or $z$, then remaining components
            q = [0.0, 0.0, 0.0]
            i = 0
            if m.m[1][1] > m.m[0][0]:
                i = 1
            if m.m[2][2] > m.m[i][i]:
                i = 2
            j = (i+1) % 3
            k = (j+1) % 3
            s = math.sqrt((m.m[i][i] - (m.m[j][j] + m.m[k][k])) + 1.0)
            q[i] = s * 0.5
            if s != 0.0:
                s = 0.5 / s
            w = (m.m[k][j] - m.m[j][k]) * s
            q[j] = (m.m[j][i] + m.m[i][j]) * s
            q[k] = (m.m[k][i] + m.m[i][k]) * s
            v.x = q[0]
            v.y = q[1]
            v.z = q[2]
        return cls(v, w)
    
    def to_transform(self):
        """Return the equivalent Transform."""
        xx = self.v.x * self.v.x
        yy = self.v.y * self.v.y
        zz = self.v.z * self.v.z
        xy = self.v.x * self.v.y
        xz = self.v.x * self.v.z
        yz = self.v.y * self.v.z
        wx = self.v.x * self.w
        wy = self.v.y * self.w
        wz = self.v.z * self.w

        m = Matrix4x4()
        m.m[0][0] = 1.0 - 2.0 * (yy + zz)
        m.m[0][1] =       2.0 * (xy + wz)
        m.m[0][2] =       2.0 * (xz - wy)
        m.m[1][0] =       2.0 * (xy - wz)
        m.m[1][1] = 1.0 - 2.0 * (xx + zz)
        m.m[1][2] =       2.0 * (yz + wx)
        m.m[2][0] =       2.0 * (xz + wy)
        m.m[2][1] =       2.0 * (yz - wx)
        m.m[2][2] = 1.0 - 2.0 * (xx + yy)
            
        return Transform(transpose(m), m)
        
    def __add__(self, q):
        """Overload the addition operator."""
        return Quaternion(self.v + q.v,
                          self.w + q.w)

    def __sub__(self, q):
        """Overload the subtraction operator."""
        return Quaternion(self.v - q.v,
                          self.w - q.w)
                          
    def __mul__(self, f):
        """Overload the multiplication operator."""
        return Quaternion(self.v * f,
                          self.w * f)

    def __rmul__(self, f):
        """Overload the right multiplication operator."""
        return Quaternion(self.v * f,
                          self.w * f)

    def __div__(self, f):
        """Overload the division operator."""
        return Quaternion(self.v / f,
                          self.w / f)

    def __eq__(self, q):
        """Overload the comparison operator."""
        return self.v==q.v and eq(self.w, q.w)
    
    def __ne__(self, q):
        """Overload the coparison operator."""
        return (self.v!=q.v) or (not eq(self.w, q.w))

    def __str__(self):
        """Return a string describing the quaternion."""
        return "Quaternion (%s, %f)" % (self.v, self.w)
    
def dot_quaternions(q1, q2):
    """Dot product of two quaternions."""
    return dot(q1.v, q2.v) + q1.w * q2.w

def normalize(q):
    """Return a nornalized version of the quaternion."""
    return q / math.sqrt(dot_quaternions(q, q))

def slerp(t, q1, q2):
    """Spherical Linear Interpolation between two quaternions."""
    cos_theta = dot_quaternions(q1, q2)
    if cos_theta > 0.9995:
        return normalize((1.0 - t)*q1 + t*q2)
    else:
        theta = math.acos(clamp(cos_theta, -1.0, 1.0))
        theta_p = theta * t
        q_perp = normalize(q2 - q1*cos_theta)
        return q1*math.cos(theta_p) + q_perp*math.sin(theta_p)
