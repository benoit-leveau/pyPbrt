"""Transformation Classes."""

import math
import copy

from core.pbrt import eq, lerp
from core.geometry import Point, Vector, Normal, Ray, RayDifferential, BBox
from core.geometry import normalize, cross, union
from core.quaternion import Quaternion, slerp


class Matrix4x4(object):

    """Class describing a 4x4 transformation matrix."""

    def __init__(self,
                 t00=1.0, t01=0.0, t02=0.0, t03=0.0,
                 t10=0.0, t11=1.0, t12=0.0, t13=0.0,
                 t20=0.0, t21=0.0, t22=1.0, t23=0.0,
                 t30=0.0, t31=0.0, t32=0.0, t33=1.0,):
        """Default constructor for Matrix4x4."""
        self.m = [[t00, t01, t02, t03],
                  [t10, t11, t12, t13],
                  [t20, t21, t22, t23],
                  [t30, t31, t32, t33]]

    @classmethod
    def from_matrix4x4(cls, matrix):
        """Copy constructor."""
        return cls(
            matrix.m[0][0], matrix.m[0][1], matrix.m[0][2], matrix.m[0][3],
            matrix.m[1][0], matrix.m[1][1], matrix.m[1][2], matrix.m[1][3],
            matrix.m[2][0], matrix.m[2][1], matrix.m[2][2], matrix.m[2][3],
            matrix.m[3][0], matrix.m[3][1], matrix.m[3][2], matrix.m[3][3]
            )

    @classmethod
    def from_array(cls, array):
        """Constructor form 4x4 array."""
        return cls(
            array[0][0], array[0][1], array[0][2], array[0][3],
            array[1][0], array[1][1], array[1][2], array[1][3],
            array[2][0], array[2][1], array[2][2], array[2][3],
            array[3][0], array[3][1], array[3][2], array[3][3]
            )

    def __mul__(self, m2):
        """Overload the multiplication operator."""
        r = Matrix4x4()
        for i in range(4):
            for j in range(4):
                r.m[i][j] = self.m[i][0] * m2.m[0][j] + \
                            self.m[i][1] * m2.m[1][j] + \
                            self.m[i][2] * m2.m[2][j] + \
                            self.m[i][3] * m2.m[3][j]
        return r

    def __eq__(self, m2):
        """Overload the comparison operator."""
        for i in range(4):
            for j in range(4):
                if not eq(self.m[i][j], m2.m[i][j]):
                    return False
        return True

    def __ne__(self, m2):
        """Overload the comparison operator."""
        for i in range(4):
            for j in range(4):
                if not eq(self.m[i][j], m2.m[i][j]):
                    return True
        return False

    def __str__(self):
        """Return a string describing the matrix."""
        return "Matrix4x4 (%s)" % str(self.m)


def transpose(m):
    """Return the transpose of m."""
    if isinstance(m, Matrix4x4):
        return Matrix4x4(m.m[0][0], m.m[1][0], m.m[2][0], m.m[3][0],
                         m.m[0][1], m.m[1][1], m.m[2][1], m.m[3][1],
                         m.m[0][2], m.m[1][2], m.m[2][2], m.m[3][2],
                         m.m[0][3], m.m[1][3], m.m[2][3], m.m[3][3])
    elif isinstance(m, Transform):
        return Transform(transpose(m.m),
                         transpose(m.m_inv))

class Transform(object):

    """Class describing a 3D transformation."""

    def __init__(self, matrix=None, matrix_inverse=None):
        """Default constructor for Transform."""
        if matrix is None:
            self.m = Matrix4x4()
            self.m_inv = Matrix4x4()
        else:
            self.m = Matrix4x4.from_matrix4x4(matrix)
            if matrix_inverse is None:
                self.m_inv = inverse(self.m)
            else:
                self.m_inv = Matrix4x4.from_matrix4x4(matrix_inverse)

    @classmethod
    def from_transform(cls, transform):
        """Copy constructor."""
        return cls(transform.m, transform.m_inv)
    
    def inverse(self):
        """Return the inverse of the transform."""
        return Transform(self.m_inv, self.m)
        
    def is_identity(self):
        """Return True if it is the identity transform."""
        return self.m == Matrix4x4()

    def has_scale(self):
        """Return True if the transform has a scaling term."""
        la2 = self(Vector(1, 0, 0)).length_squared()
        lb2 = self(Vector(0, 1, 0)).length_squared()
        lc2 = self(Vector(0, 0, 1)).length_squared()
        not_one = lambda x: x<0.999 or x>1.001
        return not_one(la2) or not_one(lb2) or not_one(lc2)

    def swap_handedness(self):
        """Return True if matrix has changed handedness."""
        det = ((self.m.m[0][0] * \
                (self.m.m[1][1] * self.m.m[2][2] - \
                 self.m.m[1][2] * self.m.m[2][1])) - \
               (self.m.m[0][1] * \
                (self.m.m[1][0] * self.m.m[2][2] - \
                 self.m.m[1][2] * self.m.m[2][0])) + \
               (self.m.m[0][2] * \
                (self.m.m[1][0] * self.m.m[2][1] - \
                 self.m.m[1][1] * self.m.m[2][0])))
        return det < 0.0

    def __eq__(self, t):
        """Overload the comparison operator."""
        return self.m == t.m and self.m_inv == t.m_inv
    
    def __ne__(self, t):
        """Overload the comparison operator."""
        return self.m != t.m or self.m_inv != t.m_inv
    
    def __call__(self, elt):
        """Overload the operator().

        Supported operations:
        * Transform(Point)
        * Transform(Vector)
        * Transform(Normal)
        * Transform(Ray)
        * Transform(RayDifferential)
        * Transform(BBox)
        
        """
        if isinstance(elt, Point):
            x = elt.x
            y = elt.y
            z = elt.z
            xp = self.m.m[0][0]*x + self.m.m[0][1]*y + self.m.m[0][2]*z + self.m.m[0][3]
            yp = self.m.m[1][0]*x + self.m.m[1][1]*y + self.m.m[1][2]*z + self.m.m[1][3]
            zp = self.m.m[2][0]*x + self.m.m[2][1]*y + self.m.m[2][2]*z + self.m.m[2][3]
            wp = self.m.m[3][0]*x + self.m.m[3][1]*y + self.m.m[3][2]*z + self.m.m[3][3]
            if wp == 1.0:
                return Point(xp, yp, zp)
            else:
                return Point(xp, yp, zp)/wp
        elif isinstance(elt, Vector):
            x = elt.x
            y = elt.y
            z = elt.z
            xp = self.m.m[0][0]*x + self.m.m[0][1]*y + self.m.m[0][2]*z
            yp = self.m.m[1][0]*x + self.m.m[1][1]*y + self.m.m[1][2]*z
            zp = self.m.m[2][0]*x + self.m.m[2][1]*y + self.m.m[2][2]*z
            return Vector(xp, yp, zp)
        elif isinstance(elt, Normal):
            x = elt.x
            y = elt.y
            z = elt.z
            return Normal(self.m_inv.m[0][0]*x + self.m_inv.m[1][0]*y + self.m_inv.m[2][0]*z,
                          self.m_inv.m[0][1]*x + self.m_inv.m[1][1]*y + self.m_inv.m[2][1]*z,
                          self.m_inv.m[0][2]*x + self.m_inv.m[1][2]*y + self.m_inv.m[2][2]*z)
        elif isinstance(elt, RayDifferential):
            ray = RayDifferential.from_ray_differential(elt)
            ray.o = self(ray.o)
            ray.d = self(ray.d)
            ray.rx_origin = self(ray.rx_origin)
            ray.ry_origin = self(ray.ry_origin)
            ray.rx_direction = self(ray.rx_direction)
            ray.ry_direction = self(ray.ry_direction)
            return ray
        elif isinstance(elt, Ray):
            ray = Ray.from_ray(elt)
            ray.o = self(ray.o)
            ray.d = self(ray.d)
            return ray
        elif isinstance(elt, BBox):
            ret = BBox(self(Point(elt.p_min.x, elt.p_min.y, elt.p_min.z)))
            ret = union(ret, self(Point(elt.p_max.x, elt.p_min.y, elt.p_min.z)))
            ret = union(ret, self(Point(elt.p_min.x, elt.p_max.y, elt.p_min.z)))
            ret = union(ret, self(Point(elt.p_min.x, elt.p_min.y, elt.p_max.z)))
            ret = union(ret, self(Point(elt.p_min.x, elt.p_max.y, elt.p_max.z)))
            ret = union(ret, self(Point(elt.p_max.x, elt.p_max.y, elt.p_min.z)))
            ret = union(ret, self(Point(elt.p_max.x, elt.p_min.y, elt.p_max.z)))
            ret = union(ret, self(Point(elt.p_max.x, elt.p_max.y, elt.p_max.z)))
            return ret

    def __mul__(self, t):
        """Overload the multiplication operator for TxT."""
        m = self.m * t.m
        m_inv = t.m_inv * self.m_inv
        return Transform(m, m_inv)

    def __str__(self):
        """Return a string describing the transform."""
        return "Transform (m='%s', m_inv='%s')" % (str(self.m), str(self.m_inv))


def translate(delta):
    """Construct a Transform representing the translation by delta."""
    m = Matrix4x4(1.0, 0.0, 0.0, delta.x,
                  0.0, 1.0, 0.0, delta.y,
                  0.0, 0.0, 1.0, delta.z,
                  0.0, 0.0, 0.0, 1.0)
    m_inv = Matrix4x4(1.0, 0.0, 0.0, -delta.x,
                      0.0, 1.0, 0.0, -delta.y,
                      0.0, 0.0, 1.0, -delta.z,
                      0.0, 0.0, 0.0, 1.0)
    return Transform(m, m_inv)


def scale(x, y, z):
    """Construct a Transform representing a scale by (x, y, z)."""
    m = Matrix4x4(  x, 0.0, 0.0, 0.0,
                  0.0,   y, 0.0, 0.0,
                  0.0, 0.0,   z, 0.0,
                  0.0, 0.0, 0.0, 1.0)
    m_inv = Matrix4x4( 1.0/x, 0.0,   0.0,   0.0,
                       0.0,   1.0/y, 0.0,   0.0,
                       0.0,   0.0,   1.0/z, 0.0,
                       0.0,   0.0,   0.0,   1.0)
    return Transform(m, m_inv)


def rotate_x(angle):
    """Construct a Transform representing a rotation around x axis."""
    sin_t = math.sin(math.radians(angle))
    cos_t = math.cos(math.radians(angle))
    m = Matrix4x4(1.0, 0.0,      0.0, 0.0,
                  0.0, cos_t, -sin_t, 0.0,
                  0.0, sin_t,  cos_t, 0.0,
                  0.0, 0.0,      0.0, 1.0)
    return Transform(m, transpose(m))


def rotate_y(angle):
    """Construct a Transform representing a rotation around y axis."""
    sin_t = math.sin(math.radians(angle))
    cos_t = math.cos(math.radians(angle))
    m = Matrix4x4(cos_t,  0.0, sin_t, 0.0,
                  0.0,    1.0,   0.0, 0.0,
                  -sin_t, 0.0, cos_t, 0.0,
                  0.0,    0.0,   0.0, 1.0)
    return Transform(m, transpose(m))


def rotate_z(angle):
    """Construct a Transform representing a rotation around z axis."""
    sin_t = math.sin(math.radians(angle))
    cos_t = math.cos(math.radians(angle))
    m = Matrix4x4(cos_t, -sin_t, 0.0, 0.0,
                  sin_t,  cos_t, 0.0, 0.0,
                    0.0,    0.0, 1.0, 0.0,
                    0.0,    0.0, 0.0, 1.0)
    return Transform(m, transpose(m))


def rotate(angle, axis):
    """Construct a Transform representing a rotation around the specified axis."""
    a = normalize(axis)
    sin_t = math.sin(math.radians(angle))
    cos_t = math.cos(math.radians(angle))
    mat = Matrix4x4(a.x * a.x + (1.0 - a.x * a.x) * cos_t,
                    a.x * a.y * (1.0 - cos_t) - a.z * sin_t,
                    a.x * a.z * (1.0 - cos_t) + a.y * sin_t,
                    0.0,
                    a.x * a.y * (1.0 - cos_t) + a.z * sin_t,
                    a.y * a.y + (1.0 - a.y * a.y) * cos_t,
                    a.y * a.z * (1.0 - cos_t) - a.x * sin_t,
                    0.0,
                    a.x * a.z * (1.0 - cos_t) - a.y * sin_t,
                    a.y * a.z * (1.0 - cos_t) + a.x * sin_t,
                    a.z * a.z + (1.0 - a.z * a.z) * cos_t,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    1.0)
    return Transform(mat, transpose(mat))


def look_at(pos, look, up):
    """Construct a Transform corresponding to a viewpoint in world space."""
    cam_to_world = Matrix4x4()
    
    # initialize fourth column of the viewing matrix
    cam_to_world.m[0][3] = pos.x
    cam_to_world.m[1][3] = pos.y
    cam_to_world.m[2][3] = pos.z
    cam_to_world.m[3][3] = 1.0

    # construct the base
    dir = normalize(look-pos)
    left = normalize(cross(normalize(up), dir))
    new_up = cross(dir, left)
    
    # fill the other columns
    cam_to_world.m[0][0] = left.x
    cam_to_world.m[1][0] = left.y
    cam_to_world.m[2][0] = left.z
    cam_to_world.m[3][0] = 0.0
    cam_to_world.m[0][1] = new_up.x
    cam_to_world.m[1][1] = new_up.y
    cam_to_world.m[2][1] = new_up.z
    cam_to_world.m[3][1] = 0.0
    cam_to_world.m[0][2] = dir.x
    cam_to_world.m[1][2] = dir.y
    cam_to_world.m[2][2] = dir.z
    cam_to_world.m[3][2] = 0.0

    return Transform(inverse(cam_to_world), cam_to_world)


def inverse(m):
    """Return the inverse of matrix m."""
    if isinstance(m, Transform):
        return Transform(m.m_inv, m.m)
    indxc = [0, 0, 0, 0]
    indxr = [0, 0, 0, 0]
    ipiv  = [0, 0, 0, 0]
    minv = copy.deepcopy(m.m)
    for i in range(4):
        irow = -1
        icol = -1
        big = 0.0
        # Choose pivot
        for j in range(4):
            if (ipiv[j] != 1):
                for k in range(4):
                    if (ipiv[k] == 0):
                        if (abs(minv[j][k]) >= big):
                            big = float(abs(minv[j][k]))
                            irow = j
                            icol = k
                    elif (ipiv[k] > 1):
                        raise Exception("Singular matrix in MatrixInvert")
        ipiv[icol] += 1
        # Swap rows _irow_ and _icol_ for pivot
        if (irow != icol):
            for k in range(4):
                # swap
                minv[irow][k], minv[icol][k] = minv[icol][k], minv[irow][k]
                
        indxr[i] = irow
        indxc[i] = icol
        if (minv[icol][icol] == 0.0):
            raise Exception("Singular matrix in MatrixInvert")

        # Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
        pivinv = 1.0 / minv[icol][icol]
        minv[icol][icol] = 1.0
        for j in range(4):
            minv[icol][j] *= pivinv

        # Subtract this row from others to zero out their columns
        for j in range(4):
            if (j != icol):
                save = minv[j][icol]
                minv[j][icol] = 0
                for k in range(4):
                    minv[j][k] -= minv[icol][k]*save

    # Swap columns to reflect permutation
    for j in range(3,-1,-1):
        if (indxr[j] != indxc[j]):
            for k in range(4):
                # swap
                minv[k][indxr[j]], minv[k][indxc[j]] = \
                                   minv[k][indxc[j]], minv[k][indxr[j]]
    return Matrix4x4.from_array(minv)


class AnimatedTransform(object):

    """Class describing an Animated Transform."""

    def __init__(self, transform1, time1, transform2, time2):
        """Default constructor for AnimatedTransform."""
        self.start_time = float(time1)
        self.end_time = float(time2)
        self.start_transform = Transform.from_transform(transform1)
        self.end_transform = Transform.from_transform(transform2)
        self.t_0, self.r_o, self.s_o = decompose(self.start_transform.m)
        self.t_1, self.r_1, self.s_1 = decompose(self.end_transform.m)
        self.actually_animated = transform1 != transform2

    @classmethod
    def from_animatedtransform(cls, transform):
        """Copy constructor."""
        return cls(transform.start_transform, transform.start_time,
                   transform.end_transform, transform.end_time)
    
    def interpolate(self, time):
        """."""
        # Handle boundary conditions for matrix interpolation
        if (not self.actually_animated or time <= self.start_time):
            return Transform.from_transform(self.start_transform)
        if (time >= self.endTime):
            return self.end_transform

        dt = (time - self.start_time) / (self.end_ime - self.start_time)
        # Interpolate translation at _dt_
        trans = (1.0 - dt) * self.t_0 + dt * self.t_1

        # Interpolate rotation at _dt_
        rotate = slerp(dt, self.r_0, self.r_1)

        # Interpolate scale at _dt_
        scale = Matrix4x4()
        for i in range(3):
            for j in range(3):
                scale.m[i][j] = lerp(dt, self.s_o.m[i][j], self.s_1.m[i][j])

        # Compute interpolated matrix as product of interpolated components
        return translate(trans) * \
               rotate.to_transform() * \
               Transform(scale)

    def motion_bounds(self, bbox, use_inverse):
        """."""
        if (not self.actually_animated):
            return inverse(self.start_transform)(bbox)
        ret = BBox()
        n_steps = 128
        for i in range(n_steps):
            t = Transform()
            time = lerp(float(i)/float(n_steps-1),
                        self.start_time,
                        self.end_time)
            t = self.Interpolate(time)
            if (use_inverse):
                t = inverse(t)
            raise Exception("check_code_next_line")
            # ret = union(ret, t(b))
        return ret

    def has_scale(self):
        """Return True if the transform has a scale component."""
        return self.start_transform.has_scale() or \
               self.end_transform.has_scale()
    
    def __call__(self, *args):
        """Overload the operator().

        Supported operations:
        * AnimatedTransform(Ray)
        * AnimatedTransform(RayDifferential)
        * AnimatedTransform(float, Point)
        * AnimatedTransform(float, Vector)
        
        """
        
        if len(args) == 1:
            if isinstance(args[0], RayDifferential):
                return self.__call_raydiff(args[0])
            elif isinstance(args[0], Ray):
                return self.__call_ray(args[0])
        elif len(args) == 2:
            if isinstance(args[0], float) and isinstance(args[1], Point):
                return self.__call_float_point(args[0], args[1])
            elif isinstance(args[0], float) and isinstance(args[1], Vector):
                return self.__call_float_vector(args[0], args[1])
        raise TypeError("Invalid call.")

    def __call_ray(self, ray):
        """Implementation of AnimatedTransform(Ray)."""
        if (not self.actuall_animated) or ray.time <= self.start_time:
            ray_interpolated = self.start_transform(ray)
        elif ray.time >= self.end_time:
            ray_interpolated = self.end_transform(ray)
        else:
            transform = self.interpolate(ray.time)
            ray_interpolated = transform(ray)
        ray_interpolated.time = ray.time
        return ray_interpolated

    def __call_rayd(self, raydiff):
        """Implementation of AnimatedTransform(RayDifferential)."""
        if (not self.actuall_animated) or raydiff.time <= self.start_time:
            ray_interpolated = self.start_transform(raydiff)
        elif raydiff.time >= self.end_time:
            ray_interpolated = self.end_transform(raydiff)
        else:
            transform = self.interpolate(raydiff.time)
            ray_interpolated = transform(raydiff)
        ray_interpolated.time = raydiff.time
        return ray_interpolated

    def __call_float_point(self, time, point):
        """Implementation of AnimatedTransform(float, Point)."""
        if (not self.actually_animated) or time <= self.start_time:
            return self.start_transform(point)
        elif time >= self.end_time:
            return self.end_transform(point)
        t = self.interpolate(time)
        return t(point)

    def __call_float_vector(self, time, vector):
        """Implementation of AnimatedTransform(float, Vector)."""
        if (not self.actually_animated) or time <= self.start_time:
            return self.start_transform(vector)
        elif time >= self.end_time:
            return self.end_transform(vector)
        t = self.interpolate(time)
        return t(vector)

    def __str__(self):
        """Return a string describing the animated transform."""
        return "AnimatedTransform (s=%d, e=%f, m_s='%s', m_e='%s')" % (self.start_time, self.end_time, str(self.start_transform), str(self.end_transform))
        

def perspective(fov, near, far):
    """Build a Transform corresponding to a perspective field of view."""
    # Perform projective divide
    persp = Matrix4x4(1, 0,            0,               0,
                      0, 1,            0,               0,
                      0, 0,  far/ (far - near),  -far*near / (far - near),
                      0, 0,            1,               0);

    # Scale to canonical viewing volume
    inv_tan_ang = 1.0 / math.tan(math.radians(fov) / 2.0)
    return scale(inv_tan_ang, inv_tan_ang, 1.0) * Transform(persp)


def decompose(matrix):
    """Decompose a Matrix4x4 into T, R and S components.

    Returns:
    T , Vector
    R , Quaternion
    S , Matrix4x4

    """
    R = Quaternion()
    S = Matrix4x4()
    
    # Extract translation _T_ from transformation matrix
    T = Vector(matrix.m[0][3],
               matrix.m[1][3],
               matrix.m[2][3])

    # Compute new transformation matrix _M_ without translation
    M = Matrix4x4.from_matrix4x4(matrix)
    for i in range(3):
        M.m[i][3] = 0.0
        M.m[3][i] = 0.0
    M.m[3][3] = 1.0

    # Extract rotation _R_ from transformation matrix
    norm = 1.0
    count = 0
    R = Matrix4x4.from_matrix4x4(M)
    while (count<100 and norm>0.0001):
        # Compute next matrix _Rnext_ in series
        Rnext = Matrix4x4()
        Rit = inverse(transpose(R))
        for i in range(4):
            for j in range(4):
                Rnext.m[i][j] = 0.5 * (R.m[i][j] + Rit.m[i][j])

        # Compute norm of difference between _R_ and _Rnext_
        norm = 0.0
        for i in range(3):
            n = abs(R.m[i][0] - Rnext.m[i][0]) + \
                abs(R.m[i][1] - Rnext.m[i][1]) + \
                abs(R.m[i][2] - Rnext.m[i][2])
            norm = max(norm, n)

        R = Rnext
        count += 1

    # XXX TODO FIXME deal with flip...
    Rquat = Quaternion.from_transform(Transform(R))

    # Compute scale _S_ using rotation and original matrix
    S = inverse(R) * M

    return T, Rquat, S
