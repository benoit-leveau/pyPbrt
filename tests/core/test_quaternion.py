import unittest

from core.geometry import Vector, dot
from core.transform import Transform, rotate_z
from core.quaternion import Quaternion, dot_quaternions, normalize, slerp


class TestQuaternion(unittest.TestCase):
    
    def setUp(self):
        self.q = Quaternion()
        self.q1 = Quaternion(Vector(0.0, 0.0, 1.0), 0.5)
        self.q2 = Quaternion(Vector(1.0, 2.0, 3.0), 0.5)

    def test_constructor(self):
        self.assertEqual(self.q.v, Vector(0.0, 0.0, 0.0))
        self.assertEqual(self.q.w, 1.0)

    def test_transform(self):
        q3 = Quaternion.from_transform(self.q2.to_transform())
        self.assertEqual(q3, self.q2)
    
    def test_operators(self):
        q3 = self.q1 * 4.5
        self.assertEqual(q3.v, self.q1.v*4.5)
        self.assertEqual(q3.w, self.q1.w*4.5)

        q4 = self.q1 / 7.5
        self.assertEqual(q4.v, self.q1.v/7.5)
        self.assertEqual(q4.w, self.q1.w/7.5)

    def test_normalize(self):
        q = Quaternion.from_transform(rotate_z(30.0))
        self.assertEqual(normalize(q), q)

    def test_dot_w_zero(self):
        v1 = Vector(1.0, 3.0, 5.0)
        v2 = Vector(2.0, 3.0, 4.0)
        q1 = Quaternion(v1, 0.0)
        q2 = Quaternion(v2, 0.0)
        scalar = dot_quaternions(q1, q2)
        self.assertEqual(scalar, dot(v1, v2))

    def test_dot_perpendicular(self):
        v1 = Vector(1.0, 2.0, 0.0)
        v2 = Vector(0.0, 0.0, 5.0)
        q1 = Quaternion(v1, 3.0)
        q2 = Quaternion(v2, 4.0)
        scalar = dot_quaternions(q1, q2)
        self.assertEqual(scalar, 3.0*4.0)

    def test_slerp(self):
        transform1 = rotate_z(30.0)
        q1 = Quaternion.from_transform(transform1)
        transform2 = rotate_z(90.0)
        q2 = Quaternion.from_transform(transform2)
        q_interp = slerp(0.5, q1, q2)
        transform3 = rotate_z(60)
        q3 = Quaternion.from_transform(transform3)
        self.assertEqual(q_interp, q3)
        self.assertEqual(q_interp.to_transform(), transform3)

if __name__ == '__main__':
    unittest.main()
