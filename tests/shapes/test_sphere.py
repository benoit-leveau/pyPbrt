import unittest

from core.transform import Transform, translate, scale
from shapes.sphere import Sphere
from core.geometry import Ray, Vector, Point


class TestSphere(unittest.TestCase):
    
    def setUp(self):
        # create a transform
        o2w = translate(Vector(10,0,0)) * scale(1.3, 1.8, 2.0)
        w2o = o2w.inverse()

        # create the sphere
        self.sphere = Sphere(o2w, w2o, False, 1.0, -1.0, 1.0, 360)

    def test_intersect(self):
        # test an intersection
        ray = Ray(Point(20, 10, 10), Vector(-1, -1, -1))
        intersect, t_hit, ray_epsilon, dg = self.sphere.intersect(ray)
        self.assertTrue(intersect)
        intersect = self.sphere.intersect_p(ray)
        self.assertTrue(intersect)

        # test an intersection
        ray = Ray(Point(20, 10, 10), Vector(-1, 1, -1))
        intersect, t_hit, ray_epsilon, dg = self.sphere.intersect(ray)
        self.assertFalse(intersect)
        intersect = self.sphere.intersect_p(ray)
        self.assertFalse(intersect)

        # test an intersection
        ray = Ray(Point(10, 0, 0), Vector(3, 1, -2))
        intersect, t_hit, ray_epsilon, dg = self.sphere.intersect(ray)
        self.assertTrue(intersect)
        intersect = self.sphere.intersect_p(ray)
        self.assertTrue(intersect)


if __name__ == '__main__':
    unittest.main()
