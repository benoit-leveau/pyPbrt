import unittest

from core.geometry import Point, Vector, Normal, Ray, RayDifferential, BBox
from core.transform import translate, scale, rotate_x, rotate_y, rotate_z, rotate


class TestGeometry(unittest.TestCase):
    
    def setUp(self):
        pass

    def test_rotate(self):
        self.assertEqual(rotate(40, Vector(1.0, 0.0, 0.0)),
                         rotate_x(40))
        self.assertEqual(rotate(20, Vector(0.0, 1.0, 0.0)),
                         rotate_y(20))
        self.assertEqual(rotate(70, Vector(0.0, 0.0, 1.0)),
                         rotate_z(70))

    def test_transform(self):
        p = Point(1, 2, 3)
        p2 = translate(Point(10, 20, 30))(p)
        self.assertTrue(isinstance(p2, Point))
        self.assertEqual(p2, Point(11, 22, 33))
        
        v = Vector(1, 2, 3)
        v2 = translate(Point(10, 20, 30))(v)
        self.assertTrue(isinstance(v2, Vector))
        self.assertEqual(v2, Vector(1, 2, 3))

        self.assertEqual(scale(2, 3, 4)(Point(1, 2, 3)),
                         Point(2, 6, 12))
        self.assertEqual(scale(2, 3, 4)(Vector(1, 2, 3)),
                         Vector(2, 6, 12))
        self.assertEqual(rotate(90, Vector(0, 1, 0))(Normal(1, 0, 0)),
                         Normal(0, 0, -1))
    
    def test_transform_ray(self):
        ray = Ray(origin=Point(1, 2, 3),
                  direction=Vector(10, 20, 30))
        ray_transformed = translate(Point(10, 20, 30))(ray)
        
        self.assertTrue(isinstance(ray_transformed, Ray))
        self.assertEqual(ray_transformed.o, Point(11, 22, 33))
        self.assertEqual(ray_transformed.d, Vector(10, 20, 30))
    
    def test_transform_ray(self):
        ray_differential = RayDifferential(origin=Point(1,2,3),
                                           direction=Vector(10,20,30))
        ray_differential.rx_origin = Point(4,5,6)
        ray_differential.ry_origin = Point(5,6,7)
        ray_differential.rx_direction = Vector(2,3,4)
        ray_differential.ry_direction = Vector(3,4,5)
        ray_transformed = translate(Point(10,20,30))(ray_differential)

        self.assertTrue(isinstance(ray_transformed, RayDifferential))
        self.assertEqual(ray_transformed.o, Point(11,22,33))
        self.assertEqual(ray_transformed.d, Vector(10,20,30))
        self.assertEqual(ray_transformed.rx_origin, Point(14,25,36))
        self.assertEqual(ray_transformed.ry_origin, Point(15,26,37))
        self.assertEqual(ray_transformed.rx_direction, Vector(2,3,4))
        self.assertEqual(ray_transformed.ry_direction, Vector(3,4,5))

    def test_transform_bbox(self):
        box = BBox(Point(-1, -2, 0), Point(0, 3, -4))
        box_transformed = translate(Point(10, 20, 30))(box)
        self.assertTrue(isinstance(box_transformed, BBox))
        self.assertEqual(box_transformed.pMin, Point(9, 18, 26))
        self.assertEqual(box_transformed.pMax, Point(10, 23, 30))
        
        box_transformed2 = scale(2, 3, 4)(box)
        self.assertTrue(isinstance(box_transformed2, BBox))
        self.assertEqual(box_transformed2.pMin, Point(-2, -6, -16))
        self.assertEqual(box_transformed2.pMax, Point(0, 9, 0))

    def test_transform_transform(self):
        m1 = scale(2, 3, 4) * translate(Point(10, 20, 30))
        m2 = translate(Point(20, 60, 120)) * scale(2, 3, 4)
        self.assertEqual(m1, m2)

    def test_transform_handedness(self):
        m1 = translate(Point(-17, 2, 31)) * scale(0.5, 6 , 1.4) * rotate(35, Vector(-15, 20, 0.2))
        self.assertFalse(m1.swap_handedness())

        m2 = translate(Point(5, 6, 7)) * scale(2, -3 , 4) * rotate(17, Vector(-1, 4, -2))
        self.assertTrue(m2.swap_handedness())


if __name__ == '__main__':
    unittest.main()
