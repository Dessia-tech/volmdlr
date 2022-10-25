"""
Unit tests for volmdlr.wires.Wire
"""
import unittest
import volmdlr
from volmdlr.models import contours


class TestWire(unittest.TestCase):

    def test_extract_primitives(self):
        point1 = volmdlr.Point2D(-0.007116025403784438, 0.010325317547305484)
        point2 = volmdlr.Point2D(-0.005383974596215561, 0.011325317547305485)

        prim_true = contours.Contour2d_1.extract_primitives(point1,
                                                            contours.Contour2d_1.primitives[3],
                                                            point2,
                                                            contours.Contour2d_1.primitives[3],
                                                            inside=True)

        prims_false = contours.Contour2d_1.extract_primitives(point1,
                                                              contours.Contour2d_1.primitives[3],
                                                              point2,
                                                              contours.Contour2d_1.primitives[3],
                                                              inside=False)

        self.assertAlmostEqual(volmdlr.wires.Wire2D(prim_true).length(),
                               0.002002125855992895)
        self.assertAlmostEqual(volmdlr.wires.Wire2D(prims_false).length(),
                               0.0797864911978587)
        self.assertEqual(len(prim_true), 1)
        self.assertEqual(len(prims_false), 7)

    def test_extract_without_primitives(self):
        point1 = volmdlr.Point2D(-0.007116025403784438, 0.010325317547305484)
        point2 = volmdlr.Point2D(-0.005383974596215561, 0.011325317547305485)

        prim_true = contours.Contour2d_1.extract_without_primitives(point1, point2, inside=True)
        prims_false = contours.Contour2d_1.extract_without_primitives(point1, point2, inside=False)

        self.assertAlmostEqual(volmdlr.wires.Wire2D(prim_true).length(),
                               0.002002125855992895)
        self.assertAlmostEqual(volmdlr.wires.Wire2D(prims_false).length(),
                               0.0797864911978587)
        self.assertEqual(len(prim_true), 1)
        self.assertEqual(len(prims_false), 7)


if __name__ == '__main__':
    unittest.main()
