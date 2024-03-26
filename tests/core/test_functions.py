"""
Unit tests for functions defined in volmdlr.core.py
"""
import unittest
import volmdlr
from volmdlr.core import delete_double_point, get_point_index_in_list
from volmdlr.models.contours import contour3d_all_edges, arc_ellipse


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

    def test_get_point_index_in_list(self):
        result = get_point_index_in_list(volmdlr.Point3D(0.0, 0.0, 0.0), self.points)
        self.assertEqual(result, 0)
        result = get_point_index_in_list(volmdlr.Point3D(1.0, 0.0, 0.0), self.points)
        self.assertIsNone(result)

    def test_get_point_index_in_list(self):
        result = get_point_index_in_list(arc_ellipse, contour3d_all_edges.primitives)
        self.assertEqual(result, 2)





if __name__ == "__main__":
    unittest.main()
