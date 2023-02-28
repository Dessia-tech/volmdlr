"""
Unittest for volmdlr.Vector3D
"""
import unittest

import volmdlr
from volmdlr import Vector3D, Point3D


class TestVector3D(unittest.TestCase):

    def setUp(self):
        self.v1 = Vector3D(1, 2, 3, "v1")
        self.v2 = Vector3D(3, 4, 5, "v2")
        self.v3 = Vector3D(1.00000001, 2.00000001, 3.00000001)
        self.v4 = Vector3D(2, 3, 4)
        self.v = Vector3D(1, 2, 3, "v1 bis")
        self.point1 = Point3D(2, 3, 4)
        self.point2 = Point3D(1, 2, 2)

    def tearDown(self) -> None:
        self.v1 = None
        self.v2 = None
        self.v3 = None
        self.v4 = None
        self.v = None
        self.point1 = None
        self.point2 = None

    def test___repr__(self):
        self.assertEqual(repr(self.v1), "Vector3D: [1.0, 2.0, 3.0]")
        self.assertEqual(repr(self.v2), "Vector3D: [3.0, 4.0, 5.0]")

    def test___setitem__(self):
        self.v1[0] = 4
        self.v1[1] = 5
        self.v1[2] = 6
        self.assertEqual(self.v1.x, 4)
        self.assertEqual(self.v1.y, 5)
        self.assertEqual(self.v1.z, 6)
        with self.assertRaises(IndexError):
            self.v1[3] = 7

    def test___getitem__(self):
        self.assertEqual(self.v1[0], 1)
        self.assertEqual(self.v1[1], 2)
        self.assertEqual(self.v1[2], 3)
        with self.assertRaises(IndexError):
            self.v1[3]

    def test___add__(self):
        result = self.v1 + self.v2
        self.assertEqual(result.x, 4)
        self.assertEqual(result.y, 6)
        self.assertEqual(result.z, 8)

    def test___neg__(self):
        result = -self.v1
        self.assertEqual(result.x, -1)
        self.assertEqual(result.y, -2)
        self.assertEqual(result.z, -3)

    def test___sub__(self):
        result = self.v1 - self.v2
        self.assertEqual(result.x, -2)
        self.assertEqual(result.y, -2)
        self.assertEqual(result.z, -2)

    def test___mul__(self):
        result = self.v1 * 2
        self.assertEqual(result.x, 2)
        self.assertEqual(result.y, 4)
        self.assertEqual(result.z, 6)

    def test___truediv__(self):
        result = self.v1 / 2
        self.assertEqual(result.x, 0.5)
        self.assertEqual(result.y, 1)
        self.assertEqual(result.z, 1.5)
        with self.assertRaises(ZeroDivisionError):
            self.v1 / 0

    def test___round__(self):
        result = round(self.v3, 1)
        self.assertEqual(result.x, 1.0)
        self.assertEqual(result.y, 2.0)
        self.assertEqual(result.z, 3.0)

    def test___hash__(self):
        self.assertEqual(hash(self.v1), hash(("vector", 1, 2, 3)))
        self.assertEqual(hash(self.v2), hash(("vector", 3, 4, 5)))

    def test___eq__(self):
        point1 = self.v1.to_point()
        self.assertTrue(self.v1 == self.v1)
        self.assertFalse(self.v1 == self.v2)
        self.assertFalse(self.v1 == point1)

    def test_is_close(self):
        self.assertTrue(self.v3.is_close(Vector3D(1.000000011, 2.000000009, 3.000000012)))

    def test_dot(self):
        self.assertEqual(self.v1.dot(self.v4), 20)

    def test_cross(self):
        expected_vector = Vector3D(-1, 2, -1)
        self.assertEqual(self.v1.cross(self.v4), expected_vector)

    def test_norm(self):
        self.assertAlmostEqual(self.v1.norm(), 3.7416573867739413)

    def test_normalize(self):
        self.v1.normalize()
        self.assertAlmostEqual(self.v1.norm(), 1.0)

    def test_point_distance(self):
        self.assertAlmostEqual(self.v1.point_distance(self.point1), 1.7320508075688772)

    def test_rotation(self):
        axis = Vector3D(1, 0, 0)
        angle = 3.141592653589793
        center = Point3D(0, 0, 0)
        rotated_vector = self.v1.rotation(center, axis, angle)
        expected_vector = Vector3D(1, -2, -3)
        self.assertAlmostEqual(rotated_vector.x, expected_vector.x)
        self.assertAlmostEqual(rotated_vector.y, expected_vector.y)
        self.assertAlmostEqual(rotated_vector.z, expected_vector.z)

    def test_vector_projection(self):
        v1 = volmdlr.Z3D
        v2 = Vector3D(1, 1, 1)
        v3 = Vector3D(1, 1, -1)
        v4 = volmdlr.X3D
        p1 = v2.vector_projection(v1)
        p2 = v3.vector_projection(v1)
        p3 = v4.vector_projection(v1)
        self.assertEqual(p1, Vector3D(0, 0, 1))
        self.assertEqual(p2, Vector3D(0, 0, -1))
        self.assertEqual(p3, Vector3D(0, 0, 0))


if __name__ == '__main__':
    unittest.main()
