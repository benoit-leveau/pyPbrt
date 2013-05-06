import unittest

from core.geometry import Vector
from core.transform import rotate_x, rotate_y, rotate_z, rotate


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

if __name__ == '__main__':
    unittest.main()
