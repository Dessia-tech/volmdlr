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
        bounding_rectangle = self.arc_ellipse2d.bounding_rectangle()
        self.assertEqual(bounding_rectangle.bounds(), (-0.3, 0.0, -0.08713203435567374, 0.33713203435567374))

    def test_straight_line_area(self):
        self.assertEqual(self.arc_ellipse2d.straight_line_area(), 0.00663975840258857)


if __name__ == '__main__':
    unittest.main()
