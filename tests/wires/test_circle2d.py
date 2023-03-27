import unittest

from geomdl import utilities

import volmdlr.edges
import volmdlr.wires as vmw

circle = vmw.Circle2D(volmdlr.O2D, 0.50)
line = volmdlr.edges.Line2D(volmdlr.O2D, volmdlr.Point2D(0, 1))

arc1_validate = volmdlr.edges.Arc2D(volmdlr.Point2D(0, -0.5), volmdlr.Point2D(0.5, 0), volmdlr.Point2D(0, 0.5))
arc2_validate = volmdlr.edges.Arc2D(volmdlr.Point2D(0, -0.5), volmdlr.Point2D(-0.5, 0), volmdlr.Point2D(0, 0.5))


class TestCircle2D(unittest.TestCase):
    def test_split_by_line(self):
        list_arcs = circle.split_by_line(line=line)

        self.assertEqual(list_arcs[0].length(), arc1_validate.length())
        self.assertEqual(list_arcs[0].start, arc1_validate.start)
        self.assertEqual(list_arcs[0].end, arc1_validate.end)
        self.assertEqual(list_arcs[1].length(), arc2_validate.length())
        self.assertEqual(list_arcs[1].start, arc2_validate.start)
        self.assertEqual(list_arcs[1].end, arc2_validate.end)

    def test_point_distance(self):
        circle_ = vmw.Circle2D(volmdlr.Point2D(1.5, 0), 1)
        point1 = volmdlr.Point2D(0.5410304786421145, 0.2835091834608327)
        self.assertEqual(circle_.point_distance(point1), 0.0)
        point2 = volmdlr.Point2D(2, 1.5)
        self.assertEqual(circle_.point_distance(point2), 0.5811388300841898)

    def test_bspline_intersections(self):
        degree = 3
        points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1), volmdlr.Point2D(2, -1), volmdlr.Point2D(3, 0)]
        knotvector = utilities.generate_knot_vector(degree, len(points))
        knot_multiplicity = [1] * len(knotvector)
        bspline = volmdlr.edges.BSplineCurve2D(degree, points, knot_multiplicity, knotvector, None, False)
        circle_ = vmw.Circle2D(volmdlr.Point2D(1.5, 0), 1)
        circle_intersections = circle_.bsplinecurve_intersections(bspline)
        expected_intersections = [volmdlr.Point2D(0.5410304786421145, 0.2835091834608327),
                                  volmdlr.Point2D(2.4589695213572873, -0.2835091834628551)]
        for intersection, expected_intersection in zip(circle_intersections, expected_intersections):
            self.assertTrue(intersection.is_close(expected_intersection))

    def test_circle_intersections(self):
        circle1 = vmw.Circle2D(volmdlr.Point2D(0, 0), 1)
        circle2 = vmw.Circle2D(volmdlr.Point2D(1, 1), 1)
        circle_intersections = circle1.circle_intersections(circle2)
        self.assertEqual(len(circle_intersections), 2)
        self.assertTrue(circle_intersections[0].is_close(volmdlr.Point2D(1.0, 0.0)))
        self.assertTrue(circle_intersections[1].is_close(0.0, 1.0))


if __name__ == '__main__':
    unittest.main()
