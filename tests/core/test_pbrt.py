import unittest

from core.pbrt import round_up_pow_2, lerp, eq, clamp, quadratic


class TestPBRT(unittest.TestCase):
    
    def test_round_up_pow_2(self):
        self.assertEqual(round_up_pow_2(11), 16)
        self.assertEqual(round_up_pow_2(0), 0)
        self.assertEqual(round_up_pow_2(1), 1)
        self.assertEqual(round_up_pow_2(2), 2)
        self.assertEqual(round_up_pow_2(3), 4)
        self.assertEqual(round_up_pow_2(112), 128)
        self.assertEqual(round_up_pow_2(255), 256)
        self.assertEqual(round_up_pow_2(256), 256)
        self.assertEqual(round_up_pow_2(pow(2,40)), pow(2,40))
        self.assertEqual(round_up_pow_2(pow(2,45)-103), pow(2,45))

    def test_lerp(self):
        v1 = 10.0
        v2 = 20.0
        self.assertEqual(lerp(0.0, v1, v2), v1)
        self.assertEqual(lerp(1.0, v1, v2), v2)
        self.assertEqual(lerp(0.5, v1, v2), 0.5*(v1+v2))

    def test_eq(self):
        self.assertTrue(eq(1.0, 1.0+1e-15))
        self.assertFalse(eq(1.0, 1.0+1e-5))
        self.assertFalse(eq(1.0, float('inf')))
        self.assertFalse(eq(1.0, -float('inf')))
        self.assertFalse(eq(1.0, 3.0))

    def test_clamp(self):
        low = -10.0
        high = 10.0
        self.assertEqual(clamp(3.4, low, high), 3.4)
        self.assertEqual(clamp(-100.0, low, high), low)
        self.assertEqual(clamp(100.0, low, high), high)

    def test_quadratic(self):
        # Ax2 + Bx + C = 0
        # use several known values to test it
        found, t0, t1 = quadratic(1.0, -3.0, -4.0)
        self.assertTrue(found)
        self.assertEqual(t0, -1.0)
        self.assertEqual(t1, 4.0)
        found, t0, t1 = quadratic(1.0, -3.0, 4.0)
        self.assertFalse(found)
        found, t0, t1 = quadratic(1.0, 0.0, -4.0)
        self.assertTrue(found)
        self.assertEqual(t0, -2.0)
        self.assertEqual(t1, 2.0)
        found, t0, t1 = quadratic(6.0, 11.0, -35.0)
        self.assertTrue(found)
        self.assertEqual(t0, -7.0/2.0)
        self.assertEqual(t1, 5.0/3.0)

if __name__ == '__main__':
    unittest.main()
