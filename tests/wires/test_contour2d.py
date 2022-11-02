import unittest

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


if __name__ == '__main__':
    unittest.main()