"""
Unitests for wires.Ellipse2D
"""
import math
import unittest

import volmdlr
from volmdlr import edges, curves


class TestEllipse2D(unittest.TestCase):
    ellipse2d = curves.Ellipse2D(4, 2, volmdlr.Frame2D(volmdlr.O2D,
                                                       volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475),
                                                       volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)))
    discretized_points = ellipse2d.discretization_points(number_points=11)

    def test_length(self):
        self.assertAlmostEqual(self.ellipse2d.length(), 19.376896441095354)

    def test_area(self):
        self.assertEqual(self.ellipse2d.area(), 25.132741228718345)

    def test_line_intersections(self):
        line = curves.Line2D(volmdlr.O2D, volmdlr.Point2D(2, 3))
        line_intersections = self.ellipse2d.line_intersections(line)
        self.assertEqual(len(line_intersections), 2)
        self.assertTrue(line_intersections[1].is_close(volmdlr.Point2D(-2.1009029257555607,
                                                                       -3.151354388633341)))
        self.assertTrue(line_intersections[0].is_close(volmdlr.Point2D(2.1009029257555607,
                                                                       3.151354388633341)))

        ellipse2d = curves.Ellipse2D(1, 0.5,
                                     volmdlr.Frame2D(
                                         origin=volmdlr.Point2D(0.8660254037844388, -2.220446049250313e-16),
                                         u=volmdlr.Vector2D(1.0, -2.220446049250313e-16),
                                         v=volmdlr.Vector2D(-4.440892098500627e-16, 1.0)))
        line2d = curves.Line2D(volmdlr.Point2D(1.0, 0.0), volmdlr.Point2D(0.7660444431189781, 0.6427876096865393))
        intersections = ellipse2d.line_intersections(line2d)
        self.assertTrue(len(intersections), 2)
        self.assertTrue(intersections[0].is_close(volmdlr.Point2D(1.173187451788891, -0.47582861312286373)))
        self.assertTrue(intersections[1].is_close(volmdlr.Point2D(0.8182229267692979, 0.4994284040759036)))
        line2d = curves.Line2D(volmdlr.Point2D(.25, -1), volmdlr.Point2D(0.25, 1))
        intersections = ellipse2d.line_intersections(line2d)
        self.assertTrue(len(intersections), 2)
        self.assertTrue(intersections[0].is_close(volmdlr.Point2D(0.24999999999999933, 0.39386314307517345)))
        self.assertTrue(intersections[1].is_close(volmdlr.Point2D(0.24999999999999978, -0.39386314307517367)))

    def test_linesegment_intersections(self):
        line_segment = edges.LineSegment2D(volmdlr.O2D, volmdlr.Point2D(4, 4))
        linesegment_intersections = self.ellipse2d.linesegment_intersections(line_segment)
        self.assertEqual(len(linesegment_intersections), 1)
        self.assertTrue(linesegment_intersections[0].is_close(volmdlr.Point2D(2.82842712474619,
                                                                              2.82842712474619)))

        ellipse2d = curves.Ellipse2D(
            1.5029657596067132, 1.4999999999999993,
            volmdlr.Frame2D(origin=volmdlr.Point2D(-0.9980228269288582, -8.858169481883975e-16),
                            u=volmdlr.Vector2D(1.0, -8.840689907867823e-16), v=volmdlr.Vector2D(0.0, -1.0)))
        lineseg2d = edges.LineSegment2D(volmdlr.Point2D(0.17364817766693041, 0.984807753012208),
                                        volmdlr.Point2D(-0.4999999999999998, 0.8660254037844388))
        lineseg_intersections = ellipse2d.linesegment_intersections(lineseg2d)
        self.assertEqual(len(lineseg_intersections), 1)
        self.assertTrue(lineseg_intersections[0].is_close(volmdlr.Point2D(0.1406944277231128, 0.9789971177815928)))

    def test_discretization_points(self):
        self.assertTrue(self.discretized_points[0].is_close(volmdlr.Point2D(2.8284271247461903, 2.82842712474619)))
        self.assertTrue(self.discretized_points[5].is_close(volmdlr.Point2D(-2.8284271247461903, -2.82842712474619)))
        self.assertTrue(self.discretized_points[-2].is_close(volmdlr.Point2D(3.119499486825644, 1.4569917357158295)))

    def test_point_belongs(self):
        self.assertTrue(self.ellipse2d.point_belongs(self.discretized_points[3]))
        self.assertFalse(self.ellipse2d.point_belongs(volmdlr.Point2D(2, 2)))

    def test_point_inside(self):
        ellipse2d = curves.Ellipse2D(2, 1, volmdlr.Frame2D(volmdlr.O2D, volmdlr.X2D, -volmdlr.Y2D))

        point2d = volmdlr.Point2D(1.6, 0.7)
        self.assertFalse(ellipse2d.point_inside(point2d))
        point2d = volmdlr.Point2D(1.6, 0.5)
        self.assertTrue(ellipse2d.point_inside(point2d))

    def test_abscissa(self):
        self.assertEqual(self.ellipse2d.abscissa(self.discretized_points[5]), 9.688448220547677)

    def test_point_angle_with_major_dir(self):
        point_angle_with_major_axis = self.ellipse2d.point_angle_with_major_dir(
            self.discretized_points[5])
        self.assertEqual(point_angle_with_major_axis, math.pi)

    def test_rotation(self):
        rotationed_ellipse = self.ellipse2d.rotation(volmdlr.O2D, math.pi / 4)
        rotationed_major_axis_point = rotationed_ellipse.center +\
            rotationed_ellipse.major_axis * rotationed_ellipse.major_dir
        self.assertTrue(rotationed_major_axis_point.is_close(volmdlr.Point2D(0.0, 4.0)))

    def test_translation(self):
        translated_ellipse = self.ellipse2d.translation(volmdlr.Vector2D(1, 0))
        translated_ellipse_major_axis_point = translated_ellipse.center +\
            translated_ellipse.major_axis * translated_ellipse.major_dir
        self.assertTrue(translated_ellipse_major_axis_point.is_close(volmdlr.Point2D(3.8284271247461903,
                                                                                     2.8284271247461903)))

    def test_frame_mapping(self):
        frame_mapped_ellipse = self.ellipse2d.frame_mapping(
            volmdlr.Frame2D(volmdlr.Point2D(1, 1), self.ellipse2d.major_dir, self.ellipse2d.minor_dir), 'new')
        frame_mapped_ellipse_major_axis_point = frame_mapped_ellipse.center +\
            frame_mapped_ellipse.major_axis * frame_mapped_ellipse.major_dir
        self.assertTrue(frame_mapped_ellipse_major_axis_point.is_close(volmdlr.Point2D(2.585786437626905, 0.0)))

    def test_point_distance(self):
        ellipse2d = curves.Ellipse2D(2, 1, volmdlr.Frame2D(volmdlr.O2D, volmdlr.X2D, -volmdlr.Y2D))

        point2d = volmdlr.Point2D(1.6, 0.7)

        point_distance = ellipse2d.point_distance(point2d)
        self.assertAlmostEqual(point_distance, 0.08415399818595351)


if __name__ == '__main__':
    unittest.main(verbosity=0)
