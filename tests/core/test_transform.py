import unittest

from core.geometry import Point, Vector
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
        self.assertEqual(translate(Point(10,10,10))(Point(1,2,3)),
                         Point(11,12,13))
        self.assertEqual(translate(Point(10,10,10))(Vector(1,2,3)),
                         Vector(1,2,3))
        self.assertEqual(scale(2,3,4)(Point(1,2,3)),
                         Point(2,6,12))
        self.assertEqual(scale(2,3,4)(Vector(1,2,3)),
                         Vector(2,6,12))

if __name__ == '__main__':
    unittest.main()
