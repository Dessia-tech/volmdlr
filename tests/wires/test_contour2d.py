import unittest

from dessia_common.core import DessiaObject

import volmdlr
from volmdlr import edges, wires
from volmdlr.models.contours import contour2d_1, contour2d_2


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
        pass

    def test_offset(self):
        contour_to_offset = DessiaObject.load_from_file('wires/contour_to_offset.json')
        stringer_contour_offset = contour_to_offset.offset(4)
        self.assertEqual(len(stringer_contour_offset.primitives), 10)
        self.assertAlmostEqual(stringer_contour_offset.area(), 546.1486677646163)

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

    def test_closest_point_to_point2(self):
        point1 = volmdlr.Point2D(1.5, -1.5)
        point2 = volmdlr.Point2D(-1, -1)
        closest_point1 = self.contour2.closest_point_to_point2(point1)
        self.assertEqual(closest_point1, volmdlr.Point2D(1.0, -1.0))
        closest_point2 = self.contour2.closest_point_to_point2(point2)
        self.assertEqual(closest_point2, volmdlr.Point2D(-2.0, 0.7))

    def test_furthest_point_to_point2(self):
        point1 = volmdlr.Point2D(1.5, -1.5)
        point2 = volmdlr.Point2D(-1, -1)
        furthest_point1 = self.contour2.get_furthest_point_to_point2(point1)
        self.assertEqual(furthest_point1, volmdlr.Point2D(-2.0, 1.0))
        furthest_point2 = self.contour2.get_furthest_point_to_point2(point2)
        self.assertEqual(furthest_point2, volmdlr.Point2D(1.5, 1.0))


if __name__ == '__main__':
    unittest.main()
