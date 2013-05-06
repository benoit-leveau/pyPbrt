import unittest
import math
import random

from core.geometry import Vector, Point, Normal, Ray, RayDifferential, BBox
from core.geometry import face_forward, union

class TestGeometry(unittest.TestCase):
    
    def setUp(self):
        self.v1 = Vector(1, 1, 1)
        self.v2 = Vector(2, 2, 2)
        self.p1 = Point(1, 1, 1)
        self.p2 = Point(2, 2, 2)

    def get_random_point(self):
        return Point(float(random.randint(0,100)) / 20.0,
                     float(random.randint(0,100)) / 20.0,
                     float(random.randint(0,100)) / 20.0)
    
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
        self.assertTrue(isinstance(r2, Ray))
        self.assertEqual(r2.d, r.d)

        # test constructor from parent ray
        r3 = Ray.from_ray_parent(r.o, r.d, r, r.mint)
        self.assertTrue(isinstance(r3, Ray))
        self.assertEqual(r3.depth, r.depth+1)

        # test operator()
        p = r(1.7)        
        self.assertEqual(p, Point(1.7, 3.4, 5.1))        

    def test_ray_differential(self):
        r = Ray(Point(0, 0, 0), Vector(1, 2, 3))
        rd = RayDifferential(Point(0, 0, 0), Vector(1, 2, 3))
        
        # test copy constructor from Ray
        rd.has_differentials = True
        rd1 = RayDifferential.from_ray_differential(rd)
        self.assertTrue(isinstance(rd1, RayDifferential))        
        self.assertEqual(rd1.o, rd.o)
        self.assertEqual(rd1.d, rd.d)
        self.assertEqual(rd1.rx_origin, rd.rx_origin)
        self.assertEqual(rd1.ry_origin, rd.ry_origin)
        self.assertEqual(rd1.rx_direction, rd.rx_direction)
        self.assertEqual(rd1.ry_direction, rd.ry_direction)
        self.assertEqual(rd1.has_differentials, rd.has_differentials)

        # test copy constructor from Ray
        rd2 = RayDifferential.from_ray(r)
        self.assertTrue(isinstance(rd2, RayDifferential))
        self.assertEqual(rd2.d, r.d)
        self.assertEqual(rd2.has_differentials, False)

        # test constructor from parent ray
        rd3 = RayDifferential.from_ray_parent(r.o, r.d, r, r.mint)
        self.assertTrue(isinstance(rd3, RayDifferential))
        self.assertEqual(rd3.depth, r.depth+1)

        # test operator()
        p = rd(1.7)        
        self.assertEqual(p, Point(1.7, 3.4, 5.1))

    def test_bounding_box_1(self):
        # test default constructor
        b = BBox()
        self.assertEqual(b.pMin,
                         Point(float('inf'), float('inf'), float('inf')))
        self.assertEqual(b.pMax,
                         Point(-float('inf'), -float('inf'), -float('inf')))

    def test_bounding_box_2(self):
        # test constructor from one point
        p = self.get_random_point()
        b1 = BBox(p)
        self.assertEqual(b1.pMin, p)
        self.assertEqual(b1.pMin, b1.pMax)

    def test_bounding_box_3(self):
        # test constructor from two points
        p1 = self.get_random_point()
        p2 = self.get_random_point()
        b2 = BBox(p1, p2)
        for i in range(3):
            self.assertEqual(b2.pMin[i], min(p1[i], p2[i]))
            self.assertEqual(b2.pMax[i], max(p1[i], p2[i]))

    def test_bounding_box_4(self):
        # test copy constructor
        bbox = BBox(Point(5, 5, 5),
                    Point(7, 7, 7))
        bbox2 = BBox.from_bbox(bbox)
        self.assertEqual(bbox.pMin, bbox2.pMin)
        self.assertEqual(bbox.pMax, bbox2.pMax)

        p1 = Point(6, 5.5, 7)
        p2 = Point(6, 7.5, 7)
        p3 = Point(6, 6.5, 4.5)
        
        # test methods
        self.assertEqual(bbox.inside(p1), True)
        self.assertEqual(bbox.inside(p2), False)
        self.assertEqual(bbox.inside(p3), False)
        
        bbox.expand(1)
        self.assertEqual(bbox.inside(p1), True)
        self.assertEqual(bbox.inside(p2), True)
        self.assertEqual(bbox.inside(p3), True)

    def test_bounding_box_5(self):
        bbox1 = BBox(Point(0, -2, 0),
                     Point(1, -1, 1))
        p1 = Point(-3,3,0.5)
        bbox2 = union(bbox1, p1)
        self.assertEqual(bbox2.pMin, Point(-3, -2, 0))
        self.assertEqual(bbox2.pMax, Point(1, 3, 1))

    def test_bounding_box_6(self):
        bbox1 = BBox(Point(0, -2, 0),
                     Point(1, -1, 1))
        bbox2 = BBox(Point(-2, 0, -2),
                     Point(-1, 1, -1))
        bbox3 = union(bbox1, bbox2)
        self.assertEqual(bbox3.pMin, Point(-2, -2, -2))
        self.assertEqual(bbox3.pMax, Point(1, 1, 1))

    def test_bounding_box_6(self):
        bbox1 = BBox(Point(0, 0, 0),
                     Point(2, 2, 2))

        bbox2 = BBox(Point(2.5, 2.5, 2.5),
                     Point(3, 3, 3))
        self.assertEqual(bbox1.overlaps(bbox2), False)

        bbox3 = BBox(Point(-1, -1, -1),
                     Point(0.5, 0.5, 0.5))
        self.assertEqual(bbox1.overlaps(bbox3), True)
        
if __name__ == '__main__':
    unittest.main()
