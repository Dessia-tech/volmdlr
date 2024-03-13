"""
Unit tests for functions defined in volmdlr.core.py
"""
import unittest
import volmdlr
from volmdlr.core import delete_double_point


class TestDeleteDoublePoint(unittest.TestCase):
    def setUp(self):
        self.points = [
            volmdlr.Point3D(0.0, 0.0, 0.0),
            volmdlr.Point3D(1.0, 1.0, 1.0),
            volmdlr.Point3D(0.0, 0.0, 0.0),
            volmdlr.Point3D(2.0, 2.0, 2.0),
        ]

    def test_delete_double_point(self):
        result = delete_double_point(self.points)

        self.assertEqual(len(result), 3)
        self.assertTrue(volmdlr.Point3D(0.0, 0.0, 0.0) in result)
        self.assertTrue(volmdlr.Point3D(1.0, 1.0, 1.0) in result)
        self.assertTrue(volmdlr.Point3D(2.0, 2.0, 2.0) in result)


if __name__ == "__main__":
    unittest.main()
