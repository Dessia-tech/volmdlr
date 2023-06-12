import unittest
from itertools import product

import volmdlr
from volmdlr import edges


class TestArcEllipse2D(unittest.TestCase):
    arc_ellipse2d = edges.ArcEllipse2D(start=volmdlr.Point2D(-0.125, -0.08416500663326211),
                                       interior=volmdlr.Point2D(-0.03543560762586048, -0.011930639375832372),
                                       end=volmdlr.Point2D(0.0, 0.125), center=volmdlr.Point2D(-0.15, 0.125),
                                       major_dir=volmdlr.Vector2D(0, 1))
    discretized_points = arc_ellipse2d.discretization_points(number_points=3)

    def test_init(self):
        list_points = [volmdlr.Point2D(2, 0), volmdlr.Point2D(1.4142135623730951, 0.7071067811865475),
                       volmdlr.Point2D(0, 1), volmdlr.Point2D(-1.4142135623730951, 0.7071067811865475),
                       volmdlr.Point2D(-2, 0), volmdlr.Point2D(-1.4142135623730951, -0.7071067811865475),
                       volmdlr.Point2D(0, -1), volmdlr.Point2D(1.4142135623730951, -0.7071067811865475)]
        list_lengths = []
        for point1, point2, point3 in product(list_points, repeat=3):
            if point1 != point2 and point1 != point3 and point2 != point3:
                if not (abs(point3.x) == abs(point1.x) == abs(point2.x)):
                    arc_ellipse2d = edges.ArcEllipse2D(start=point1, interior=point2,
                                                       end=point3, center=volmdlr.Point2D(0, 0),
                                                       major_dir=volmdlr.Vector2D(1, 0))
                    list_lengths.append(round(arc_ellipse2d.length(), 7))
        list_expected_legths = [2.4221121, 3.8785604, 4.8442241, 5.8098879, 7.2663362, 8.7227845, 8.7227845, 3.8785604,
                                4.8442241, 5.8098879, 7.2663362, 8.7227845, 8.7227845, 7.2663362, 4.8442241, 5.8098879,
                                7.2663362, 8.7227845, 8.7227845, 7.2663362, 5.8098879, 5.8098879, 7.2663362, 8.7227845,
                                8.7227845, 7.2663362, 5.8098879, 4.8442241, 7.2663362, 8.7227845, 8.7227845, 7.2663362,
                                5.8098879, 4.8442241, 3.8785604, 8.7227845, 8.7227845, 7.2663362, 5.8098879, 4.8442241,
                                3.8785604, 2.4221121, 8.2319999, 6.7755516, 5.8098879, 4.8442241, 3.3877758, 1.9313275,
                                8.7227845, 2.9128966, 3.8785604, 4.8442241, 6.3006724, 7.7571207, 8.7227845, 8.2319999,
                                3.8785604, 6.3006724, 8.7227845, 8.2319999, 6.7755516, 4.8442241, 6.3006724, 7.7571207,
                                8.7227845, 8.2319999, 5.8098879, 6.3006724, 8.7227845, 8.2319999, 6.7755516, 5.8098879,
                                4.8442241, 7.7571207, 8.7227845, 8.2319999, 5.8098879, 3.3877758, 8.2319999, 8.2319999,
                                7.2663362, 6.3006724, 4.8442241, 3.3877758, 2.4221121, 8.2319999, 7.2663362, 6.3006724,
                                4.8442241, 3.3877758, 7.2663362, 8.2319999, 2.4221121, 3.3877758, 4.8442241, 6.3006724,
                                7.2663362, 8.2319999, 8.2319999, 3.3877758, 4.8442241, 6.3006724, 7.2663362, 8.2319999,
                                8.2319999, 7.2663362, 4.8442241, 6.3006724, 7.2663362, 8.2319999, 8.2319999, 7.2663362,
                                6.3006724, 6.3006724, 7.2663362, 8.2319999, 8.2319999, 7.2663362, 6.3006724, 4.8442241,
                                6.7755516, 8.2319999, 8.7227845, 7.7571207, 6.3006724, 4.8442241, 3.8785604, 8.2319999,
                                8.7227845, 6.3006724, 3.8785604, 2.9128966, 8.7227845, 7.7571207, 6.3006724, 4.8442241,
                                5.8098879, 6.7755516, 8.2319999, 1.9313275, 3.3877758, 4.8442241, 5.8098879, 8.2319999,
                                8.7227845, 3.3877758, 5.8098879, 6.7755516, 8.2319999, 8.7227845, 7.7571207, 4.8442241,
                                5.8098879, 8.2319999, 8.7227845, 6.3006724, 5.8098879, 7.2663362, 8.7227845, 8.7227845,
                                7.2663362, 5.8098879, 4.8442241, 7.2663362, 8.7227845, 8.7227845, 7.2663362, 5.8098879,
                                4.8442241, 3.8785604, 8.7227845, 8.7227845, 7.2663362, 5.8098879, 4.8442241, 3.8785604,
                                2.4221121, 8.7227845, 7.2663362, 5.8098879, 4.8442241, 5.8098879, 7.2663362, 8.7227845,
                                2.4221121, 3.8785604, 4.8442241, 5.8098879, 7.2663362, 8.7227845, 8.7227845, 3.8785604,
                                4.8442241, 5.8098879, 7.2663362, 8.7227845, 8.7227845, 7.2663362, 4.8442241, 6.3006724,
                                7.7571207, 8.7227845, 8.2319999, 6.7755516, 5.8098879, 6.3006724, 8.7227845, 8.2319999,
                                5.8098879, 4.8442241, 7.7571207, 8.7227845, 8.2319999, 6.7755516, 5.8098879, 3.3877758,
                                8.7227845, 8.2319999, 5.8098879, 4.8442241, 3.3877758, 1.9313275, 8.2319999, 6.7755516,
                                3.8785604, 4.8442241, 6.3006724, 7.7571207, 8.7227845, 2.9128966, 3.8785604, 6.3006724,
                                8.7227845, 8.2319999, 3.3877758, 4.8442241, 6.3006724, 7.2663362, 8.2319999, 8.2319999,
                                7.2663362, 4.8442241, 6.3006724, 7.2663362, 8.2319999, 8.2319999, 7.2663362, 6.3006724,
                                6.3006724, 7.2663362, 8.2319999, 8.2319999, 7.2663362, 6.3006724, 4.8442241, 7.2663362,
                                8.2319999, 8.2319999, 7.2663362, 6.3006724, 4.8442241, 3.3877758, 8.2319999, 8.2319999,
                                7.2663362, 6.3006724, 4.8442241, 3.3877758, 2.4221121, 8.2319999, 2.4221121, 3.3877758,
                                4.8442241, 6.3006724, 7.2663362, 8.2319999, 1.9313275, 3.3877758, 4.8442241, 5.8098879,
                                6.7755516, 8.2319999, 8.7227845, 3.3877758, 5.8098879, 8.2319999, 8.7227845, 7.7571207,
                                4.8442241, 5.8098879, 6.7755516, 8.2319999, 8.7227845, 6.3006724, 5.8098879, 8.2319999,
                                8.7227845, 7.7571207, 6.3006724, 4.8442241, 6.7755516, 8.2319999, 8.7227845, 6.3006724,
                                3.8785604, 8.2319999, 8.7227845, 7.7571207, 6.3006724, 4.8442241, 3.8785604, 2.9128966]
        for length, expected_length in zip(list_lengths, list_expected_legths):
            self.assertAlmostEqual(length, expected_length)
        with self.assertRaises(ValueError):
            edges.ArcEllipse2D(
                start=volmdlr.Point2D(1.4142135623730951, 0.7071067811865475),
                interior=volmdlr.Point2D(-1.4142135623730951, 0.7071067811865475),
                end=volmdlr.Point2D(-1.4142135623730951, -0.7071067811865475), center=volmdlr.Point2D(0, 0),
                major_dir=volmdlr.Vector2D(1, 0))

    def test_length(self):
        self.assertAlmostEqual(self.arc_ellipse2d.length(), 0.2612813847745195)

    def test_discretization_points(self):
        expected_discretized_points = [volmdlr.Point2D(-0.1250000000001241, -0.0841650066332621),
                                       volmdlr.Point2D(-0.03543560762614463, -0.011930639376171975),
                                       volmdlr.Point2D(0.0, 0.12499999999999996)]
        for expected_point, point in zip(expected_discretized_points, self.discretized_points):
            self.assertTrue(expected_point.is_close(point))

    def test_point_belongs(self):
        self.assertTrue(self.arc_ellipse2d.point_belongs(self.discretized_points[1]))

    def test_abscissa(self):
        self.assertAlmostEqual(self.arc_ellipse2d.abscissa(self.discretized_points[1]), 0.11816032940841212)

    def reverse(self):
        reversed_arcellipse2d = self.arc_ellipse2d.reverse()
        self.assertEqual(reversed_arcellipse2d.start, self.arc_ellipse2d.end)
        self.assertEqual(reversed_arcellipse2d.end, self.arc_ellipse2d.start)
        reversed_arcellipse2d_length = reversed_arcellipse2d.length()
        self.assertEqual(reversed_arcellipse2d_length, self.arc_ellipse2d.length())
        self.assertEqual(reversed_arcellipse2d_length - reversed_arcellipse2d.abscissa(self.discretized_points[1]),
                         self.arc_ellipse2d.abscissa(self.discretized_points[1]))

    def test_bounding_rectangle(self):
        expected_bounds = [-0.1250000000001241, 0.0, -0.0841650066332621, 0.12499999999999996]
        for i, bound in enumerate(self.arc_ellipse2d.bounding_rectangle.bounds()):
            self.assertAlmostEqual(bound, expected_bounds[i])

    def test_straight_line_area(self):
        self.assertAlmostEqual(self.arc_ellipse2d.straight_line_area(), 0.00663975840258857)

    def test_linesegment_intersections(self):
        line = edges.LineSegment2D(volmdlr.Point2D(0.0, -0.08), self.arc_ellipse2d.center)
        inters = self.arc_ellipse2d.linesegment_intersections(line)
        self.assertEqual(len(inters), 1)
        self.assertTrue(inters[0].is_close(volmdlr.Point2D(-0.04213625371737764, -0.022413786586250567)))

    def test_translation(self):
        translated_ellipse = self.arc_ellipse2d.translation(volmdlr.Vector2D(1, 1))
        self.assertTrue(translated_ellipse.start.is_close(volmdlr.Point2D(0.875, 0.9158349933667379)))
        self.assertTrue(translated_ellipse.interior.is_close(volmdlr.Point2D(0.9645643923741395, 0.9880693606241676)))
        self.assertTrue(translated_ellipse.end.is_close(volmdlr.Point2D(1.0, 1.125)))

    def test_point_distance(self):
        point = volmdlr.Point2D(0, 0)
        point_distance = self.arc_ellipse2d.point_distance(point)
        self.assertAlmostEqual(point_distance, 0.02594552758568249)

    def test_split(self):
        split1 = self.arc_ellipse2d.split(self.arc_ellipse2d.interior)
        split2 = self.arc_ellipse2d.split(self.arc_ellipse2d.start)
        split3 = self.arc_ellipse2d.split(self.arc_ellipse2d.end)
        self.assertTrue(split1[0].start.is_close(self.arc_ellipse2d.start))
        self.assertTrue(split1[0].end.is_close(self.arc_ellipse2d.interior))
        self.assertTrue(split1[1].start.is_close(self.arc_ellipse2d.interior))
        self.assertTrue(split1[1].end.is_close(self.arc_ellipse2d.end))
        self.assertIsNone(split2[0])
        self.assertEqual(split2[1], self.arc_ellipse2d)
        self.assertIsNone(split3[1])
        self.assertEqual(split3[0], self.arc_ellipse2d)

    def test_point_at_abscissa(self):
        list_abscissas = [0, self.arc_ellipse2d.length() / 4, self.arc_ellipse2d.length() / 3,
                          self.arc_ellipse2d.length() / 2, self.arc_ellipse2d.length() * 0.75,
                          self.arc_ellipse2d.length()]
        list_points = []
        for abscissa in list_abscissas:
            point_at_abscissa = self.arc_ellipse2d.point_at_abscissa(abscissa)
            list_points.append(point_at_abscissa)
        expected_points = [volmdlr.Point2D(-0.125, -0.08416500663326211),
                           volmdlr.Point2D(-0.06840405601123901, -0.053000572608877194),
                           volmdlr.Point2D(-0.053379764573816504, -0.037262319138834876),
                           volmdlr.Point2D(-0.029351056810099085, -0.0010462812393967902),
                           volmdlr.Point2D(-0.007169432517673158, 0.06020140444930423),
                           volmdlr.Point2D(0.0, 0.125)]
        for point, expected_point in zip(list_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_complementary(self):
        complementary = self.arc_ellipse2d.complementary()
        self.assertTrue(complementary.interior.is_close(volmdlr.Point2D(-0.26456439237413953, 0.2619306393758324)))
        self.assertAlmostEqual(complementary.length(), 0.884777951933013)


if __name__ == '__main__':
    unittest.main()
