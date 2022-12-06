import unittest

import volmdlr
from volmdlr import edges


class TestLine2D(unittest.TestCase):

    def test_point_distance(self):
        line_2d = edges.Line2D(volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1))
        self.assertEqual(line_2d.point_distance(volmdlr.Point2D(2, 2)), 0.0)
        self.assertEqual(line_2d.point_distance(volmdlr.Point2D(1, 2)), 0.7071067811865475)

    def test_point_belongs(self):
        line = edges.Line2D(volmdlr.O2D, volmdlr.Point2D(1, 1))
        point1 = volmdlr.Point2D(2, 2)
        point2 = volmdlr.Point2D(1, 2)
        self.assertTrue(line.point_belongs(point1))
        self.assertFalse(line.point_belongs(point2))

    def test_sort_points_along_line(self):
        line2d = edges.Line2D(volmdlr.O2D, volmdlr.Point2D(1, 2))
        list_points2d = [volmdlr.Point2D(2, 4), volmdlr.Point2D(1.5, 3),
                         volmdlr.Point2D(4, 8), volmdlr.Point2D(2.5, 5)]
        sorted_points_along_line2d = line2d.sort_points_along_line(list_points2d)
        expected_sorted_points2d = [volmdlr.Point2D(1.5, 3), volmdlr.Point2D(2, 4),
                                    volmdlr.Point2D(2.5, 5), volmdlr.Point2D(4, 8)]
        for point, expected_point in zip(sorted_points_along_line2d, expected_sorted_points2d):
            self.assertEqual(point, expected_point)


if __name__ == '__main__':
    unittest.main()
