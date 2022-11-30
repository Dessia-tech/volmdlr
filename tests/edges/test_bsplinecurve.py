"""
Unit tests for volmdlr.faces.BSplineCurve
"""
import unittest
from volmdlr.models import bspline_curves
import volmdlr
from volmdlr import edges
import dessia_common


class TestBSplineCurve(unittest.TestCase):

    def test_abscissa(self):
        bspline_curve2d = bspline_curves.bspline_curve2d_1
        point = volmdlr.Point2D(-0.31240117104573617,-2.8555856978321796)

        self.assertAlmostEqual(bspline_curve2d.abscissa(point), 7.747599410268476)

    def test_line_intersections(self):
        bspline_curve2d = dessia_common.DessiaObject.load_from_file('edges/bsplinecurve2d_1.json')
        line = edges.Line2D(volmdlr.Point2D(1.263163105753452, -0.002645572020392778),
                            volmdlr.Point2D(1.263163105753452, -0.001820963841291406))

        line_intersections = bspline_curve2d.line_intersections(line)
        self.assertEqual(len(line_intersections), 1)
        self.assertEqual(line_intersections[0], volmdlr.Point2D(1.2631631057526727, -0.0026450894385881708))


if __name__ == '__main__':
    unittest.main()

