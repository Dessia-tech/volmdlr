import unittest

import volmdlr
from volmdlr import edges, curves


class TestLine3D(unittest.TestCase):
    line1 = curves.Line3D(volmdlr.O3D, volmdlr.Point3D(0, 1, 0), name='line1')
    line2 = curves.Line3D(volmdlr.Point3D(0, 0, 1), volmdlr.Point3D(1, 0, 1), name='line2')
    line3 = curves.Line3D(volmdlr.Point3D(0, 0, 0), volmdlr.Point3D(1, 0, 1), name='line3')
    line4 = curves.Line3D(volmdlr.Point3D(0, 2, 1), volmdlr.Point3D(1, 2, 1), name='line4')
    line5 = curves.Line3D(volmdlr.Point3D(0, 0, 0), volmdlr.Point3D(1, 1, 0), name='line5')
    line6 = curves.Line3D(volmdlr.Point3D(0, 0, 0), volmdlr.Point3D(0, 1, 1), name='line6')

    def test_to_2d(self):
        line_3d = curves.Line3D(volmdlr.O3D, volmdlr.Point3D(1, 1, 1))
        line_2d = line_3d.to_2d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D)
        self.assertEqual(line_2d.point1, volmdlr.Point2D(0, 0))
        self.assertEqual(line_2d.point2, volmdlr.Point2D(1, 1))

    def test_line_distance(self):
        self.assertEqual(self.line1.line_distance(self.line2), 1.0)
        self.assertEqual(self.line1.line_distance(self.line3), 0.0)
        self.assertEqual(self.line2.line_distance(self.line4), 2.0)

    def test_skew_to(self):
        self.assertTrue(self.line1.skew_to(self.line2))
        self.assertFalse(self.line1.skew_to(self.line3))
        self.assertFalse(self.line2.skew_to(self.line4))

    def test_intersection(self):
        self.assertIsNone(self.line1.intersection(self.line1))
        self.assertIsNone(self.line2.intersection(self.line4))
        self.assertIsNone(self.line3.intersection(self.line4))
        self.assertEqual(self.line1.intersection(self.line3), volmdlr.O3D)
        self.assertEqual(self.line1.intersection(self.line5), volmdlr.O3D)
        self.assertEqual(self.line1.intersection(self.line6), volmdlr.O3D)
        self.assertEqual(self.line2.intersection(self.line3), volmdlr.Point3D(1, 0, 1))
        self.assertEqual(self.line3.intersection(self.line5), volmdlr.O3D)
        self.assertEqual(self.line5.intersection(self.line6), volmdlr.O3D)

    def test_sort_points_along_line(self):
        line3d = curves.Line3D(volmdlr.O3D, volmdlr.Point3D(1, 2, 3))
        list_points_3d = [volmdlr.Point3D(0, 0, 0), volmdlr.Point3D(5, 10, 15),
                          volmdlr.Point3D(-1, -2, -3), volmdlr.Point3D(2, 4, 6)]
        sorted_points_along_line3d = line3d.sort_points_along_curve(list_points_3d)
        expected_sorted_points3d = [volmdlr.Point3D(-1, -2, -3), volmdlr.Point3D(0, 0, 0),
                                    volmdlr.Point3D(2, 4, 6), volmdlr.Point3D(5, 10, 15)]
        for point, expected_point in zip(sorted_points_along_line3d, expected_sorted_points3d):
            self.assertEqual(point, expected_point)

    def test_point_at_abscissa(self):
        line3d = curves.Line3D(volmdlr.Point3D(-0.16532959009, 0.669230747399, 0.6255868826019999),
                               volmdlr.Point3D(-0.1653592907723126, 0.6687168836734351, 0.6247295250680509))

        abscissa1 = -0.00019339887669772427
        self.assertTrue(line3d.point_at_abscissa(abscissa1).is_close(
            volmdlr.Point3D(-0.1653238460114036, 0.6693301280663, 0.6257526945859939)))

        abscissa2 = 0.00019337655078049366
        self.assertTrue(line3d.point_at_abscissa(abscissa2).is_close(
            volmdlr.Point3D(-0.16533533350550145, 0.6691313782041791, 0.6254210897592992)))


if __name__ == '__main__':
    unittest.main()
