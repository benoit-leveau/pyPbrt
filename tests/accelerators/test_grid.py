import unittest

from core.transform import Transform, translate, scale
from shapes.sphere import Sphere
from core.primitive import GeometricPrimitive
from core.geometry import Ray, Vector, Point
from accelerators.grid import GridAccel
from core.intersection import Intersection


class TestGrid(unittest.TestCase):
    
    def setUp(self):
        # create a transform
        o2w1 = translate(Vector(10,0,0)) * scale(1.3, 1.8, 2.0)
        w2o1 = o2w1.inverse()

        # create a sphere
        self.sphere1 = Sphere(o2w1, w2o1, False, 1.0, -1.0, 1.0, 360)
        self.primitive_sphere1 = GeometricPrimitive(self.sphere1,
                                                    None, None)
        
        # create a transform
        o2w2 = translate(Vector(-10,0,0)) * scale(1.3, 1.8, 2.0)
        w2o2 = o2w2.inverse()

        # create a sphere
        self.sphere2 = Sphere(o2w2, w2o2, False, 1.0, -1.0, 1.0, 360)
        self.primitive_sphere2 = GeometricPrimitive(self.sphere2,
                                                    None, None)

        primitives = [self.primitive_sphere1, self.primitive_sphere2]

        self.grid_accel = GridAccel(primitives, True)

    def test_intersect(self):
        # test an intersection
        ray = Ray(Point(0, 0, 10), Vector(0, 0, -1))
        intersection = Intersection()
        intersected = self.grid_accel.intersect(ray, intersection)
        self.assertFalse(intersected)
        ray.maxt = float('inf')
        intersected = self.grid_accel.intersect_p(ray)
        self.assertFalse(intersected)

        # test an intersection
        ray2 = Ray(Point(10, 0, 10), Vector(0, 0, -1))
        intersection = Intersection()
        intersected = self.grid_accel.intersect(ray2, intersection)
        self.assertTrue(intersected)
        self.assertEqual(intersection.primitive_id,
                        self.primitive_sphere1.primitive_id)
        self.assertEqual(intersection.shape_id,
                        self.sphere1.shape_id)
        ray2.maxt = float('inf')
        intersected = self.grid_accel.intersect_p(ray2)
        self.assertTrue(intersected)

        # test an intersection
        ray3 = Ray(Point(-10, 0, 10), Vector(0, 0, -1))
        intersection = Intersection()
        intersected = self.grid_accel.intersect(ray3, intersection)
        self.assertTrue(intersected)
        self.assertEqual(intersection.primitive_id,
                        self.primitive_sphere2.primitive_id)
        self.assertEqual(intersection.shape_id,
                        self.sphere2.shape_id)
        ray3.maxt = float('inf')
        intersected = self.grid_accel.intersect_p(ray3)
        self.assertTrue(intersected)


if __name__ == '__main__':
    unittest.main()
