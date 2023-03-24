import math
import unittest

from dessia_common.core import DessiaObject

import volmdlr
from volmdlr import edges, wires
from volmdlr.models.contours import contour2d_1, contour2d_2


class TestContour2D(unittest.TestCase):
    contour1 = wires.Contour2D([edges.FullArc2D(center=volmdlr.O2D, start_end=volmdlr.Point2D(0.029999999, 0))])
    not_ordered_contour = DessiaObject.load_from_file('wires/contour_not_ordered.json')
    ordered_contour = DessiaObject.load_from_file('wires/contour_ordered.json')

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

    def test_crossings(self):
        line_segment1 = edges.LineSegment2D(volmdlr.Point2D(1, -1), volmdlr.Point2D(1.5, 1))
        arc = edges.Arc2D(volmdlr.Point2D(1.5, 1), volmdlr.Point2D(1.3, 1.5), volmdlr.Point2D(0.5, 1.5))
        line_segment2 = edges.LineSegment2D(volmdlr.Point2D(0.5, 1.5), volmdlr.Point2D(-2, 1))
        line_segment3 = edges.LineSegment2D(volmdlr.Point2D(-2, 1), volmdlr.Point2D(-2, 0.7))
        lie_segment4 = edges.LineSegment2D(volmdlr.Point2D(-2, 0.7), volmdlr.Point2D(-1, 1))
        points2d = [volmdlr.Point2D(-1, 1),
                    volmdlr.Point2D(2, 2),
                    volmdlr.Point2D(-2, -2),
                    volmdlr.Point2D(1, -1)]
        bspline = edges.BSplineCurve2D(3, points2d, knot_multiplicities=[4, 4], knots=[0.0, 1.0])
        contour2 = wires.Contour2D([bspline, line_segment1, arc, line_segment2, line_segment3, lie_segment4])
        contour3 = contour2.rotation(volmdlr.Point2D(0.5, 0.5), math.pi / 1.5)
        contour3 = contour3.translation(volmdlr.Vector2D(-0.3, 0))
        contour_crossings = contour2.wire_crossings(contour3)
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


if __name__ == '__main__':
    unittest.main()
