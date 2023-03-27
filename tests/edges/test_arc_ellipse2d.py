import unittest

import volmdlr
from volmdlr import edges


class TestArcEllipse2D(unittest.TestCase):
    arc_ellipse2d = edges.ArcEllipse2D(start=volmdlr.Point2D(-0.125, -0.08416500663326211),
                                       interior=volmdlr.Point2D(-0.03543560762586048, -0.011930639375832372),
                                       end=volmdlr.Point2D(0.0, 0.125), center=volmdlr.Point2D(-0.15, 0.125),
                                       major_dir=volmdlr.Vector2D(0, 1))
    discretized_points = arc_ellipse2d.discretization_points(number_points=3)

    def test_length(self):
        self.assertAlmostEqual(self.arc_ellipse2d.length(), 0.2612813847745195)

    def test_discretization_points(self):
        expected_discretized_points = [volmdlr.Point2D(-0.1250000000001241, -0.0841650066332621),
                                       volmdlr.Point2D(-0.03543560762614463, -0.011930639376171975),
                                       volmdlr.Point2D(0.0, 0.12499999999999996)]
        for expected_point, point in zip(expected_discretized_points, self.discretized_points):
            self.assertEqual(expected_point, point)

    def test_point_belongs(self):
        self.assertTrue(self.arc_ellipse2d.point_belongs(self.discretized_points[1]))

    def test_abscissa(self):
        self.assertAlmostEqual(self.arc_ellipse2d.abscissa(self.discretized_points[1]), 0.08651454398796926)

    def reverse(self):
        reversed_arcellipse2d = self.arc_ellipse2d.reverse()
        self.assertEqual(reversed_arcellipse2d.start, self.arc_ellipse2d.end)
        self.assertEqual(reversed_arcellipse2d.end, self.arc_ellipse2d.start)
        reversed_arcellipse2d_length = reversed_arcellipse2d.length()
        self.assertEqual(reversed_arcellipse2d_length, self.arc_ellipse2d.length())
        self.assertEqual(reversed_arcellipse2d_length - reversed_arcellipse2d.abscissa(self.discretized_points[1]),
                         self.arc_ellipse2d.abscissa(self.discretized_points[1]))

    def test_bounding_rectangle(self):
        self.assertEqual(self.arc_ellipse2d.bounding_rectangle.bounds(),
                         (-0.1250000000001241, 0.0, -0.0841650066332621, 0.12499999999999996))

    def test_straight_line_area(self):
        self.assertEqual(self.arc_ellipse2d.straight_line_area(), 0.00663975840258857)

    def test_linesegment_intersections(self):
        line = edges.LineSegment2D(volmdlr.Point2D(0.0, -0.08), self.arc_ellipse2d.center)
        inters = self.arc_ellipse2d.linesegment_intersections(line)
        self.assertEqual(len(inters), 1)
        self.assertEqual(inters[0], volmdlr.Point2D(-0.04213625371737764, -0.022413786586250567))

    def test_translation(self):
        translated_ellipse = self.arc_ellipse2d.translation(volmdlr.Vector2D(1, 1))
        self.assertTrue(translated_ellipse.start.is_close(volmdlr.Point2D(0.875, 0.9158349933667379)))
        self.assertTrue(translated_ellipse.interior.is_close(volmdlr.Point2D(0.9645643923741395, 0.9880693606241676)))
        self.assertTrue(translated_ellipse.end.is_close(volmdlr.Point2D(1.0, 1.125)))

    def test_point_distance(self):
        pass

    def test_split(self):
        pass


if __name__ == '__main__':
    unittest.main()
