"""
Unitests for wires.Ellipse2D
"""
import math
import unittest

import volmdlr
from volmdlr import edges, wires


class TestEllipse2D(unittest.TestCase):
    ellipse2d = wires.Ellipse2D(4, 2, volmdlr.O2D, volmdlr.Vector2D(1, 1))
    discretized_points = ellipse2d.discretization_points(number_points=10)

    def test_length(self):
        self.assertAlmostEqual(self.ellipse2d.length(), 19.376896441095354)

    def test_area(self):
        self.assertEqual(self.ellipse2d.area(), 25.132741228718345)

    def test_line_intersections(self):
        line = edges.Line2D(volmdlr.O2D, volmdlr.Point2D(2, 3))
        line_intersections = self.ellipse2d.line_intersections(line)
        self.assertEqual(len(line_intersections), 2)
        self.assertEqual(line_intersections[0], volmdlr.Point2D(-2.1009029257555607, -3.151354388633341))
        self.assertEqual(line_intersections[1], volmdlr.Point2D(2.1009029257555607, 3.151354388633341))

    def test_linesegment_intersections(self):
        line_segment = edges.LineSegment2D(volmdlr.O2D, volmdlr.Point2D(4, 4))
        linesegment_intersections = self.ellipse2d.linesegment_intersections(line_segment)
        self.assertEqual(len(linesegment_intersections), 1)
        self.assertEqual(linesegment_intersections[0], volmdlr.Point2D(2.82842712474619, 2.82842712474619))

    def test_discretization_points(self):
        self.assertEqual(self.discretized_points[0], volmdlr.Point2D(2.8284271247461903, 2.82842712474619))
        self.assertEqual(self.discretized_points[5], volmdlr.Point2D(-2.8284271247461903, -2.82842712474619))
        self.assertEqual(self.discretized_points[-2], volmdlr.Point2D(3.119499486825644, 1.4569917357158295))

    def test_point_over_ellipse(self):
        self.assertTrue(self.ellipse2d.point_over_ellipse(self.discretized_points[3]))
        self.assertFalse(self.ellipse2d.point_over_ellipse(volmdlr.Point2D(2, 2)))

    def test_abscissa(self):
        self.assertEqual(self.ellipse2d.abscissa(self.discretized_points[5]), 9.688448220547677)

    def test_point_angle_with_major_dir(self):
        point_angle_with_major_axis = self.ellipse2d.point_angle_with_major_dir(self.discretized_points[5])
        self.assertEqual(point_angle_with_major_axis, math.pi)

    def test_rotation(self):
        rotationed_ellipse = self.ellipse2d.rotation(volmdlr.O2D, math.pi / 4)
        rotationed_major_axis_point = (
            rotationed_ellipse.center + rotationed_ellipse.major_axis * rotationed_ellipse.major_dir
        )
        self.assertEqual(rotationed_major_axis_point, volmdlr.Point2D(0.0, 4.0))

    def test_translation(self):
        translated_ellipse = self.ellipse2d.translation(volmdlr.Vector2D(1, 0))
        translated_ellipse_major_axis_point = (
            translated_ellipse.center + translated_ellipse.major_axis * translated_ellipse.major_dir
        )
        self.assertEqual(translated_ellipse_major_axis_point, volmdlr.Point2D(3.8284271247461903, 2.8284271247461903))

    def test_frame_mapping(self):
        frame_mapped_ellipse = self.ellipse2d.frame_mapping(
            volmdlr.Frame2D(volmdlr.Point2D(1, 1), self.ellipse2d.major_dir, self.ellipse2d.minor_dir), "new"
        )
        frame_mapped_ellipse_major_axis_point = (
            frame_mapped_ellipse.center + frame_mapped_ellipse.major_axis * frame_mapped_ellipse.major_dir
        )
        self.assertEqual(frame_mapped_ellipse_major_axis_point, volmdlr.Point2D(2.585786437626905, 0.0))


if __name__ == "__main__":
    unittest.main(verbosity=0)
