import unittest

import dessia_common
import volmdlr
from volmdlr.models.contours import contour2d_1, contour2d_2
from volmdlr import wires, edges


class TestContour2D(unittest.TestCase):
    contour1 = wires.Contour2D([edges.FullArc2D(center=volmdlr.O2D, start_end=volmdlr.Point2D(0.029999999, 0))])
    not_ordered_contour = dessia_common.DessiaObject.load_from_file('wires/contour_not_ordered.json')
    ordered_contour = dessia_common.DessiaObject.load_from_file('wires/contour_ordered.json')

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
        for expected_primitive, primitve in zip(self.ordered_contour.primitives, ordered_contour.primitives):
            self.assertEqual(expected_primitive, primitve)

    def test_cut_by_wire(self):
        pass


if __name__ == '__main__':
    unittest.main()