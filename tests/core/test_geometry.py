import core.geometry
import unittest

class TestGeometry(unittest.TestCase):
    
    def setUp(self):
        self.v1 = core.geometry.Vector(1,1,1)
        self.v2 = core.geometry.Vector(2,2,2)
        self.p1 = core.geometry.Point(1,1,1)
        self.p2 = core.geometry.Point(2,2,2)
    
    def test_misc(self):
        self.assertEqual(self.v1 * 4, 2 * self.v2)
        self.assertEqual(self.p2 - self.p1, self.v1)
        self.assertEqual(self.v1.length_squared(), 3)
        
if __name__ == '__main__':
    unittest.main()
