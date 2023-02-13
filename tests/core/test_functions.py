"""
Unit tests for functions defined in volmdlr core.py
"""
import unittest
import volmdlr
import numpy
from volmdlr.core import determinant, delete_double_point, step_ids_to_str


# class TestDeterminant(unittest.TestCase):
#     def test_determinant_1(self):
#         vec1 = numpy.ndarray([1, 2, 3])
#         vec2 = volmdlr.Vector3D(4, 5, 6)
#         vec3 = volmdlr.Vector3D(7, 8, 9)
#         result = determinant(vec1, vec2, vec3)
#         self.assertAlmostEqual(result, 0.0)
#
#     def test_determinant_2(self):
#         vec1 = volmdlr.Vector3D(3, 2, 1)
#         vec2 = volmdlr.Vector3D(1, 2, 3)
#         vec3 = volmdlr.Vector3D(4, 5, 6)
#         result = determinant(vec1, vec2, vec3)
#         self.assertAlmostEqual(result, -9.0)


class TestDeleteDoublePoint(unittest.TestCase):
    points = [
        volmdlr.Point3D(0, 0, 0),
        volmdlr.Point3D(1, 1, 1),
        volmdlr.Point3D(0, 0, 0),
        volmdlr.Point3D(2, 2, 2),
    ]

    def test_delete_double_point(self):
        result = delete_double_point(self.points)

        self.assertEqual(len(result), 3)
        self.assertTrue(volmdlr.Point3D(0, 0, 0) in result)
        self.assertTrue(volmdlr.Point3D(1, 1, 1) in result)
        self.assertTrue(volmdlr.Point3D(2, 2, 2) in result)


class TestStepIdsToStr(unittest.TestCase):
    ids_0 = [0]
    ids_1 = [1, 2, 3, 4, 5]
    ids_2 = [11, 22, 33, 44, 55]

    def test_step_ids_to_str(self):
        self.assertEqual(step_ids_to_str(self.ids_0), "#0")
        self.assertEqual(step_ids_to_str(self.ids_1), "#1,#2,#3,#4,#5")
        self.assertEqual(step_ids_to_str(self.ids_2), "#11,#22,#33,#44,#55")


if __name__ == "__main__":
    unittest.main()
