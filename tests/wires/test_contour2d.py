import math
import unittest

from dessia_common.core import DessiaObject

import volmdlr
from volmdlr import edges, wires
from volmdlr.models.contours import contour2d_1, contour2d_2, contour1_cut_by_wire, contour2_cut_by_wire


class TestContour2D(unittest.TestCase):
    contour1 = wires.Contour2D([edges.FullArc2D(center=volmdlr.O2D, start_end=volmdlr.Point2D(0.029999999, 0))])
    not_ordered_contour = DessiaObject.load_from_file('wires/contour_not_ordered.json')
    ordered_contour = DessiaObject.load_from_file('wires/contour_ordered.json')
    contour_to_extract_from = contour = wires.Contour2D.from_points(
        [volmdlr.Point2D(-.15, .15), volmdlr.Point2D(-.15, -.15), volmdlr.Point2D(.15, -.15),
         volmdlr.Point2D(.15, .15), volmdlr.Point2D(-.15, .15)])
    point1_ = volmdlr.Point2D(0.12500000000000003, 0.15)
    point2_ = volmdlr.Point2D(0.12500000000000003, -0.15)
    point_to_extract_with = [(point1_, point2_), (volmdlr.Point2D(0.15, -0.05), volmdlr.Point2D(0.15, 0.05)),
                       (volmdlr.Point2D(-0.15, 0.15), point2_), (volmdlr.Point2D(-0.15, 0.15), point1_)]
    line_segment1 = edges.LineSegment2D(volmdlr.Point2D(1, -1), volmdlr.Point2D(1.5, 1))
    arc = edges.Arc2D(volmdlr.Point2D(1.5, 1), volmdlr.Point2D(1.3, 1.5), volmdlr.Point2D(0.5, 1.5))
    line_segment2 = edges.LineSegment2D(volmdlr.Point2D(0.5, 1.5), volmdlr.Point2D(-2, 1))
    line_segment3 = edges.LineSegment2D(volmdlr.Point2D(-2, 1), volmdlr.Point2D(-2, 0.7))
    lie_segment4 = edges.LineSegment2D(volmdlr.Point2D(-2, 0.7), volmdlr.Point2D(-1, 1))
    points2d = [volmdlr.Point2D(-1, 1), volmdlr.Point2D(2, 2), volmdlr.Point2D(-2, -2), volmdlr.Point2D(1, -1)]
    bspline = edges.BSplineCurve2D(3, points2d, knot_multiplicities=[4, 4], knots=[0.0, 1.0])
    contour2 = wires.Contour2D([bspline, line_segment1, arc, line_segment2, line_segment3, lie_segment4])
    contour3 = contour2.rotation(volmdlr.Point2D(0.5, 0.5), math.pi / 1.5)
    contour3 = contour3.translation(volmdlr.Vector2D(-0.3, 0))

    def test_point_belongs(self):
        point1 = volmdlr.Point2D(0.0144822, 0.00595264)
        point2 = volmdlr.Point2D(0.02, 0.02)
        self.assertTrue(self.contour1.point_belongs(point1))
        self.assertTrue(self.contour1.point_belongs(point2))

        point3 = volmdlr.Point2D(0, 0.013)
        self.assertTrue(contour2d_1.point_belongs(point3))
        self.assertFalse(contour2d_1.point_belongs(point1))

        point4 = volmdlr.Point2D(0.745, 0.0685)
        self.assertTrue(contour2d_2.point_belongs(point4))

    def test_is_ordered(self):
        self.assertTrue(self.ordered_contour.is_ordered())
        self.assertFalse(self.not_ordered_contour.is_ordered())

    def test_order_contour(self):
        ordered_contour = self.not_ordered_contour.order_contour()
        for previous_primitive, primitive in zip(ordered_contour.primitives, ordered_contour.primitives[1:] +
                                                                             [ordered_contour.primitives[0]]):
            self.assertEqual(previous_primitive.end, primitive.start)

    def test_cut_by_wire(self):
        results = contour1_cut_by_wire.cut_by_wire(contour2_cut_by_wire)
        results1 = contour2_cut_by_wire.cut_by_wire(contour1_cut_by_wire)
        list_expected_contour_lengths = [0.735061566418825, 0.7350615482415725, 0.7350613283792926, 0.7350615900834909,
                                         3.4786699386591753, 1.2716033189138256, 0.8996336040021796,
                                         1.2716033047094752, 0.8996337337370333, 3.4786699386591753]
        self.assertEqual(len(results)+len(results1), 10)
        for i, contour in enumerate(results + results1):
            self.assertAlmostEqual(contour.length(), list_expected_contour_lengths[i])

    def test_offset(self):
        contour_to_offset = DessiaObject.load_from_file('wires/contour_to_offset.json')
        stringer_contour_offset = contour_to_offset.offset(4)
        self.assertEqual(len(stringer_contour_offset.primitives), 10)
        self.assertAlmostEqual(stringer_contour_offset.area(), 546.1486677646163)

    def test_edge_crossings(self):
        points = [volmdlr.Point2D(-0.3, -0.2), volmdlr.Point2D(0.3, -0.2),
                  volmdlr.Point2D(0.2, 0.2), volmdlr.Point2D(0, 0.35), volmdlr.Point2D(-0.2, 0.2)]
        contour = wires.ClosedPolygon2D(points)
        lineseg = edges.LineSegment2D(volmdlr.Point2D(-.3, -.3), volmdlr.Point2D(.3, .3))
        lineseg1 = edges.LineSegment2D(volmdlr.Point2D(-.3, -.3), volmdlr.Point2D(.2, .2))
        edge_crossings = contour.edge_crossings(lineseg)
        edge_crossings1 = contour.edge_crossings(lineseg1)
        wire = wires.Wire2D([edges.LineSegment2D(p1, p2) for p1, p2 in (points[:2], points[1:3])])
        wire_edge_crossings = wire.edge_crossings(lineseg)
        self.assertEqual(len(edge_crossings), 2)
        self.assertTrue(edge_crossings[0].is_close(volmdlr.Point2D(-0.2, -0.2)))
        self.assertTrue(edge_crossings[1].is_close(volmdlr.Point2D(0.2, 0.2)))
        self.assertEqual(len(edge_crossings1), 1)
        self.assertTrue(edge_crossings1[0].is_close(volmdlr.Point2D(-0.2, -0.2)))
        self.assertEqual(len(wire_edge_crossings), 1)
        self.assertTrue(wire_edge_crossings[0].is_close(volmdlr.Point2D(-0.2, -0.2)))

    def test_crossings(self):
        contour_crossings = self.contour2.wire_crossings(self.contour3)
        expected_crossings = [volmdlr.Point2D(0.003455042474764269, 0.010365521626940934),
                              volmdlr.Point2D(-0.08592141662920678, -0.26441019657119236),
                              volmdlr.Point2D(1.3565385114925572, 0.4261540459702289),
                              volmdlr.Point2D(1.062649150008717, 1.6296941038180754),
                              volmdlr.Point2D(-0.8413589025398691, 1.0502025995770785),
                              volmdlr.Point2D(-0.16551404839994044, -0.5695209737217725),
                              volmdlr.Point2D(0.3497299962331255, -1.1558059537665162),
                              volmdlr.Point2D(0.5885392751147405, -1.116467059799847)]
        self.assertEqual(len(contour_crossings), len(expected_crossings))
        for crossing, expected_crossing in zip(contour_crossings, expected_crossings):
            self.assertTrue(crossing.is_close(expected_crossing))
        points = [volmdlr.Point2D(-0.3, -0.2), volmdlr.Point2D(0.3, -0.2),
                  volmdlr.Point2D(0.2, 0.2), volmdlr.Point2D(0, 0.3), volmdlr.Point2D(-0.2, 0.2)]
        contour = volmdlr.wires.ClosedPolygon2D(points)
        primitives = [volmdlr.edges.LineSegment2D(volmdlr.Point2D(-0.35, -0.1), volmdlr.Point2D(-0.1, 0)),
                      volmdlr.edges.LineSegment2D(volmdlr.Point2D(-0.1, 0), volmdlr.Point2D(0.2, 0.2)),
                      volmdlr.edges.LineSegment2D(volmdlr.Point2D(0.2, 0.2), volmdlr.Point2D(0.3, 0.3))]
        wire = volmdlr.wires.Wire2D(primitives)
        wire2 = volmdlr.wires.Wire2D(primitives[:-1])
        wire3 = volmdlr.wires.Wire2D(primitives[:-1] + [volmdlr.edges.LineSegment2D(
            volmdlr.Point2D(0.2, 0.2), volmdlr.Point2D(0.1, 0.0))])
        wire_crossings = contour.wire_crossings(wire)
        wire_crossings2 = contour.wire_crossings(wire2)
        wire_crossings3 = contour.wire_crossings(wire3)
        self.assertEqual(len(wire_crossings), 2)
        self.assertTrue(wire_crossings[0].is_close(volmdlr.Point2D(-0.26666666666666666, -0.0666666666666666)))
        self.assertTrue(wire_crossings[1].is_close(volmdlr.Point2D(0.2, 0.19999999999999996)))
        self.assertEqual(len(wire_crossings2), 1)
        self.assertTrue(wire_crossings2[0].is_close(volmdlr.Point2D(-0.26666666666666666, -0.0666666666666666)))
        self.assertEqual(len(wire_crossings3), 1)
        self.assertTrue(wire_crossings3[0].is_close(volmdlr.Point2D(-0.26666666666666666, -0.0666666666666666)))

    def test_split_with_two_points(self):
        expected_results = [(3, 0.3499999999999999, 3, 0.85),
                            (1, 0.1, 5, 1.0999999999999999),
                            (2, 0.575, 3, 0.625),
                            (4, 0.9249999999999998, 1, 0.275)]
        for i, (pt1, pt2) in enumerate(self.point_to_extract_with):
            inside_prims, outside_prims = self.contour_to_extract_from.split_with_two_points(pt1, pt2)
            expected_inside_ = expected_results[i][:2]
            expected_outside_ = expected_results[i][2:]
            self.assertEqual(len(inside_prims), expected_inside_[0])
            self.assertAlmostEqual(sum(prim.length() for prim in inside_prims), expected_inside_[1])
            self.assertEqual(len(outside_prims), expected_outside_[0])
            self.assertAlmostEqual(sum(prim.length() for prim in outside_prims), expected_outside_[1])

    def test_extract_with_points(self):
        list_expected_outside_params = [(3, 0.85), (5, 1.0999999999999999), (3, 0.625), (1, 0.275)]
        list_expected_inside_params = [(3, 0.3499999999999999), (1, 0.1), (2, 0.575), (4, 0.9249999999999998)]
        for i, (pt1, pt2) in enumerate(self.point_to_extract_with):
            inside_prims = self.contour_to_extract_from.extract_with_points(pt1, pt2, inside=True)
            outside_prims = self.contour_to_extract_from.extract_with_points(pt1, pt2, inside=False)
            self.assertEqual(len(inside_prims), list_expected_inside_params[i][0])
            self.assertAlmostEqual(sum(prim.length() for prim in inside_prims), list_expected_inside_params[i][1])
            self.assertEqual(len(outside_prims), list_expected_outside_params[i][0])
            self.assertAlmostEqual(sum(prim.length() for prim in outside_prims), list_expected_outside_params[i][1])

    def test_split_by_line(self):
        line = edges.Line2D(volmdlr.Point2D(volmdlr.TWO_PI, 0.1), volmdlr.Point2D(volmdlr.TWO_PI, -0.1))
        contour = wires.Contour2D.load_from_file("wires/contour_to_split.json")
        intersection = contour.line_intersections(line)[0][0]
        contour1, contour2 = contour.split_by_line(line)
        self.assertTrue(contour1.primitives[-1].end.is_close(intersection))
        self.assertTrue(contour2.primitives[0].start.is_close(intersection))

    def test_intersection_contour_with(self):
        vol = DessiaObject.load_from_file('wires/test_intersection_contour_with.json')
        contour2 = vol.primitives[0]
        contour3 = vol.primitives[1]
        intersection_contours1 = contour2.intersection_contour_with(contour3, abs_tol=1e-5)
        self.assertTrue(len(intersection_contours1), 1)
        self.assertAlmostEqual(intersection_contours1[0].length(), 0.1651419895544224)
        intersection_contours2 = self.contour2.intersection_contour_with(self.contour3, abs_tol=1e-6)
        self.assertTrue(len(intersection_contours1), 2)
        self.assertAlmostEqual(intersection_contours2[0].length(), 6.915893328290323)
        self.assertAlmostEqual(intersection_contours2[1].length(), 2.440629779494193)


if __name__ == '__main__':
    unittest.main()
