import math
import unittest
from volmdlr import Point3D, Vector3D


class TestPoint3D(unittest.TestCase):
    def setUp(self):
        self.p1 = Point3D(1.0, 2.0, 3.0, "p1")
        self.p2 = Point3D(2.0, 3.0, 4.0, "p2")
        self.p3 = Point3D(-4, 5, -6)
        self.v = Vector3D(1.0, 1.0, 1.0, "v")

    def test_addition(self):
        p = self.p1 + self.v
        self.assertAlmostEqual(p.x, 2.0)
        self.assertAlmostEqual(p.y, 3.0)
        self.assertAlmostEqual(p.z, 4.0)
        self.assertEqual(p.name, "")

    def test_subtraction(self):
        v = self.p2 - self.p1
        self.assertAlmostEqual(v.x, 1.0)
        self.assertAlmostEqual(v.y, 1.0)
        self.assertAlmostEqual(v.z, 1.0)
        self.assertEqual(v.name, "")

    def test_negation(self):
        p = -self.p1
        self.assertAlmostEqual(p.x, -1.0)
        self.assertAlmostEqual(p.y, -2.0)
        self.assertAlmostEqual(p.z, -3.0)
        self.assertEqual(p.name, "")

    def test_multiplication(self):
        p = self.p1 * 2
        self.assertAlmostEqual(p.x, 2.0)
        self.assertAlmostEqual(p.y, 4.0)
        self.assertAlmostEqual(p.z, 6.0)
        self.assertEqual(p.name, "")

    def test_division(self):
        p = self.p1 / 2
        self.assertAlmostEqual(p.x, 0.5)
        self.assertAlmostEqual(p.y, 1.0)
        self.assertAlmostEqual(p.z, 1.5)
        self.assertEqual(p.name, "")

    def test_division_by_zero(self):
        with self.assertRaises(ZeroDivisionError):
            self.p1 / 0

    def test_hash(self):
        p1 = Point3D(1.0, 2.0, 3.0, "p1")
        p2 = Point3D(1.0, 2.0, 3.0, "p2")
        p3 = Point3D(1.0, 2.0, 4.0, "p3")
        self.assertEqual(hash(p1), hash(p2))
        self.assertNotEqual(hash(p1), hash(p3))

    def test_to_dict(self):
        point = Point3D(1, 2, 3, "test_point")
        point_dict = point.to_dict()
        self.assertEqual(point_dict["object_class"], "volmdlr.Point3D")
        self.assertEqual(point_dict["x"], 1)
        self.assertEqual(point_dict["y"], 2)
        self.assertEqual(point_dict["z"], 3)
        self.assertEqual(point_dict["name"], "test_point")

    def test_to_vector(self):
        expected_output = Vector3D(1, 2, 3)
        self.assertEqual(self.p1.to_vector(), expected_output)

    def test_point_distance(self):
        expected_output = math.sqrt((self.p3.x - self.p1.x) ** 2 +
                                    (self.p3.y - self.p1.y) ** 2 +
                                    (self.p3.z - self.p1.z) ** 2)
        self.assertAlmostEqual(self.p1.point_distance(self.p3), expected_output)

    def test_middle_point(self):
        expected_output = Point3D(-1.5, 3.5, -1.5)
        self.assertEqual(Point3D.middle_point(self.p1, self.p3), expected_output)

    def test_nearest_point(self):
        point_list = [Point3D(-3, 4, -5), Point3D(2, 0, 5), Point3D(-1, 2, 0)]
        expected_output = Point3D(2, 0, 5)
        self.assertEqual(self.p1.nearest_point(point_list), expected_output)

    def test_is_close(self):
        close_point = Point3D(1.0000001, 2.0000002, 3.0000001)
        self.assertTrue(self.p1.is_close(close_point, tol=1e-6))

        far_point = Point3D(1.001, 2.001, 3.001)
        self.assertFalse(self.p1.is_close(far_point, tol=1e-6))


if __name__ == "__main__":
    unittest.main()
