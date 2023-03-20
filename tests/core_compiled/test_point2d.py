"""
Unittest for volmdlr.Point2D
"""
import unittest

import matplotlib

import volmdlr

import unittest
from volmdlr import Point2D, Vector2D, Vector3D


class TestPoint2D(unittest.TestCase):
    def setUp(self):
        self.point1 = Point2D(1.0, 2.0, "point1")
        self.point2 = Point2D(3.0, 4.0, "point2")
        self.point3 = Point2D(1.0000001, 2.0000001, "point3")
        self.vector1 = Vector2D(2.0, 3.0, "vector1")
        self.vector2 = Vector2D(2, 2, "vector2")
        self.vector3 = Vector2D(4, 4, "vector3")

    def test_add(self):
        result = self.point1 + self.vector1
        self.assertAlmostEqual(result.x, 3.0)
        self.assertAlmostEqual(result.y, 5.0)

    def test_neg(self):
        result = -self.point1
        self.assertAlmostEqual(result.x, -1.0)
        self.assertAlmostEqual(result.y, -2.0)

    def test_sub(self):
        result = self.point2 - self.point1
        self.assertAlmostEqual(result.x, 2.0)
        self.assertAlmostEqual(result.y, 2.0)

    def test_mul(self):
        result = self.point1 * 2.0
        self.assertAlmostEqual(result.x, 2.0)
        self.assertAlmostEqual(result.y, 4.0)

    def test_truediv(self):
        result = self.point1 / 2.0
        self.assertAlmostEqual(result.x, 0.5)
        self.assertAlmostEqual(result.y, 1.0)

    def test_truediv_zero(self):
        with self.assertRaises(ZeroDivisionError):
            self.point1 / 0

    def test_hash(self):
        self.assertIsInstance(hash(self.point1), int)

    def test_to_dict(self):
        result = self.point1.to_dict()
        self.assertIsInstance(result, dict)
        self.assertIn("object_class", result)
        self.assertEqual(result["object_class"], "volmdlr.Point2D")
        self.assertIn("x", result)
        self.assertEqual(result["x"], self.point1.x)
        self.assertIn("y", result)
        self.assertEqual(result["y"], self.point1.y)
        self.assertIn("name", result)
        self.assertEqual(result["name"], self.point1.name)

    def test_to_3d(self):
        plane_origin = Vector3D(1.0, 2.0, 3.0, "origin")
        vx = Vector3D(2.0, 0.0, 0.0, "vx")
        vy = Vector3D(0.0, 2.0, 0.0, "vy")
        result = self.point1.to_3d(plane_origin, vx, vy)
        self.assertAlmostEqual(result.x, 3.0)
        self.assertAlmostEqual(result.y, 6.0)
        self.assertAlmostEqual(result.z, 3.0)

    def test_to_vector(self):
        result = self.point1.to_vector()
        self.assertIsInstance(result, Vector2D)
        self.assertAlmostEqual(result.x, self.point1.x)
        self.assertAlmostEqual(result.y, self.point1.y)

    def test_plot(self):
        ax = self.point1.plot()
        self.assertIsInstance(ax, matplotlib.axes.Axes)
        lines = ax.lines
        self.assertEqual(len(lines), 1)
        line = lines[0]
        self.assertAlmostEqual(line.get_xdata()[0], self.point1.x)
        self.assertAlmostEqual(line.get_ydata()[0], self.point1.y)

    def test_point_distance(self):
        result = self.point1.point_distance(self.point2)
        self.assertAlmostEqual(result, 2.8284271247461903)

    def test_is_close(self):
        self.assertTrue(self.point1.is_close(self.point1))
        self.assertTrue(self.point1.is_close(self.point3))
        self.assertFalse(self.point1.is_close(self.point2))
        self.assertFalse(self.point1.is_close(self.vector1))
        self.assertFalse(self.point1.is_close(self.vector2))
        self.assertTrue(self.point1.is_close(self.point2, tol=3.6))

    def test_nearest_point(self):
        points = [Point2D(1, 2), Point2D(5, 6), Point2D(7, 8)]
        nearest = self.point2.nearest_point(points)
        self.assertEqual(nearest, Point2D(1, 2))


if __name__ == '__main__':
    unittest.main()
