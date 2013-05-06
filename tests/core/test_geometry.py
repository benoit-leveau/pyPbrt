import unittest
import math

from core.geometry import Vector, Point, Normal, Ray, RayDifferential
from core.geometry import face_forward

class TestGeometry(unittest.TestCase):
    
    def setUp(self):
        self.v1 = Vector(1, 1, 1)
        self.v2 = Vector(2, 2, 2)
        self.p1 = Point(1, 1, 1)
        self.p2 = Point(2, 2, 2)
    
    def test_point(self):
        # operator[]
        p = Point(1.0, 2.0, 3.0)
        self.assertEqual(p[0], 1.0)
        self.assertEqual(p[1], 2.0)
        self.assertEqual(p[2], 3.0)
        
        self.assertEqual(self.p1 * 4, 2 * self.p2)
        self.assertEqual(self.p2 - self.p1, self.v1)

    def test_vector(self):
        # operator[]
        v = Vector(1.0, 2.0, 3.0)
        self.assertEqual(v[0], 1.0)
        self.assertEqual(v[1], 2.0)
        self.assertEqual(v[2], 3.0)

        # misc operators
        self.assertEqual(self.v1 * 4, 2 * self.v2)
        self.assertEqual(self.p2 - self.p1, self.v1)

        # length functions
        self.assertEqual(self.v1.length_squared(), 3)
        self.assertEqual(self.v1.length(), math.sqrt(3))

    def test_normal(self):
        # operator[]
        n = Normal(1.0, 2.0, 3.0)
        self.assertEqual(n[0], 1.0)
        self.assertEqual(n[1], 2.0)
        self.assertEqual(n[2], 3.0)

        # face_forward
        n2 = Normal(1, 0, 0)
        v = Vector(-0.5, -0.1, 0.2)
        self.assertEqual(face_forward(n, v), -n)
    
    def test_ray(self):
        r = Ray(Point(0, 0, 0), Vector(1, 2, 3))
        
        # test copy constructor
        r2 = Ray.from_ray(r)
        self.assertEqual(r2.d, r.d)

        # test constructor from parent ray
        r3 = Ray.from_ray_parent(r.o, r.d, r, r.mint)
        self.assertEqual(r3.depth, r.depth+1)

        # test operator()
        p = r(1.7)        
        self.assertEqual(p, Point(1.7, 3.4, 5.1))        

    def test_ray_differential(self):
        r = Ray(Point(0, 0, 0), Vector(1, 2, 3))
        rd = RayDifferential(Point(0, 0, 0), Vector(1, 2, 3))
        
        # test copy constructor from Ray
        rd2 = RayDifferential.from_ray(r)
        self.assertEqual(rd2.d, r.d)
        self.assertEqual(rd2.has_differentials, False)

        # test constructor from parent ray
        rd3 = RayDifferential.from_ray_parent(r.o, r.d, r, r.mint)
        self.assertEqual(rd3.depth, r.depth+1)

        # test operator()
        p = rd(1.7)        
        self.assertEqual(p, Point(1.7, 3.4, 5.1))

if __name__ == '__main__':
    unittest.main()
