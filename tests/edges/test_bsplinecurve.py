"""
Unit tests for volmdlr.faces.BSplineCurve
"""
import unittest
from volmdlr.models import bspline_curves
import volmdlr


class TestBSplineCurve(unittest.TestCase):

    def test_abscissa(self):
        bspline_curve2d = bspline_curves.bspline_curve2d_1
        point = volmdlr.Point2D(-0.31240117104573617,-2.8555856978321796)

        self.assertAlmostEqual(bspline_curve2d.abscissa(point), 7.747599410268476)


if __name__ == '__main__':
    unittest.main()
