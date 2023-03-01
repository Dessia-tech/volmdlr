"""
Unittest for volmdlr.Vector2D
"""
import unittest
from volmdlr import Vector2D


class TestVector2D(unittest.TestCase):
    vector1 = Vector2D(1.1234567890, 2.873667521231412)
    vector2 = Vector2D(1.123456789, 2.873668521)

    def test_constructor(self):
        v1 = Vector2D(1.0, 2.0)
        self.assertEqual(v1.x, 1.0)
        self.assertEqual(v1.y, 2.0)

    def test_add(self):
        v1 = Vector2D(1.0, 2.0)
        v2 = Vector2D(2.0, 3.0)
        v3 = v1 + v2
        self.assertEqual(v3.x, 3.0)
        self.assertEqual(v3.y, 5.0)

    def test_sub(self):
        v1 = Vector2D(1.0, 2.0)
        v2 = Vector2D(2.0, 3.0)
        v3 = v1 - v2
        self.assertEqual(v3.x, -1.0)
        self.assertEqual(v3.y, -1.0)

    def test_mul(self):
        v1 = Vector2D(1.0, 2.0)
        self.assertEqual(v1 * 9, Vector2D(9, 18))

    def test_norm(self):
        v1 = Vector2D(3.0, 4.0)
        self.assertEqual(v1.norm(), 5.0)

    def test_normalize(self):
        v1 = Vector2D(3.0, 4.0)
        v1.normalize()
        self.assertAlmostEqual(v1.x, 0.6)
        self.assertAlmostEqual(v1.y, 0.8)

    def test_point_distance(self):
        v1 = Vector2D(1.0, 1.0)
        v2 = Vector2D(4.0, 5.0)
        self.assertAlmostEqual(v1.point_distance(v2), 5.0)


if __name__ == '__main__':
    unittest.main()
