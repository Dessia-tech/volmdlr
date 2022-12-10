import unittest

import dessia_common
import volmdlr
from volmdlr.models.contours import Contour2d_1
from volmdlr import wires, edges


class TestContour2D(unittest.TestCase):
    contour1 = wires.Contour2D([edges.FullArc2D(center=volmdlr.O2D, start_end=volmdlr.Point2D(0.029999999, 0))])

    def point_belongs(self):
        point1 = volmdlr.Point2D(0.0144822, 0.00595264)
        point2 = volmdlr.Point2D(0.02, 0.02)
        self.assertTrue(self.contour1.point_belongs(point1))
        self.assertFalse(self.contour1.point_belongs(point2))
        point3 = volmdlr.Point2D(0, 0.013)
        self.assertTrue(Contour2d_1.point_belongs(point3))
        self.assertFalse(Contour2d_1.point_belongs(point1))

    def test_cut_by_wire(self):
        pass

    def test_offset(self):
        contour_to_offset = dessia_common.DessiaObject.load_from_file('wires/contour_to_offset.json')
        stringer_contour_offset = contour_to_offset.offset(4)
        self.assertEqual(len(stringer_contour_offset.primitives), 10)
        self.assertAlmostEqual(stringer_contour_offset.area(), 546.1486677646163)


if __name__ == '__main__':
    unittest.main()
