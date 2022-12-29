"""
Unit tests for volmdlr.faces.BSplineCurve
"""
import unittest
from geomdl import BSpline
from geomdl import abstract

from volmdlr.models import bspline_curves
import volmdlr
import volmdlr.edges as vme


class TestBSplineCurve(unittest.TestCase):

    def test_abscissa(self):
        bspline_curve2d = bspline_curves.bspline_curve2d_1
        point = volmdlr.Point2D(-0.31240117104573617,-2.8555856978321796)

        self.assertAlmostEqual(bspline_curve2d.abscissa(point), 7.747599410268476)

    def test_discretization_points(self):
        bspline_curve2d = bspline_curves.bspline_curve2d_1

        curve = BSpline.Curve()
        curve.degree = 2
        curve.ctrlpts = [[1, 0, 0], [1, 1, 0], [0, 1, 0]]
        curve.knotvector = [0, 0, 0, 1, 1, 1]

        bspline_curve3d = vme.BSplineCurve3D.from_geomdl_curve(curve)
        # Test discretization with default number of points (20)
        points = bspline_curve3d.discretization_points()
        self.assertEqual(len(points), 20)

        # Test accuracy of first 5 discretized points
        expected_points = [volmdlr.Point3D(0.0, 0.0, 0.0),
                           volmdlr.Point3D(0.10526315789473684, 0.10526315789473684, 0.10526315789473684),
                           volmdlr.Point3D(0.21052631578947367, 0.21052631578947367, 0.21052631578947367),
                           volmdlr.Point3D(0.3157894736842105, 0.3157894736842105, 0.3157894736842105),
                           volmdlr.Point3D(0.42105263157894735, 0.42105263157894735, 0.42105263157894735)]
        for i in range(5):
            self.assertTrue(points[i], expected_points[i])

        # Test discretization with specified number of points
        points = bspline_curve2d.discretization_points(number_points=10)
        self.assertEqual(len(points), 10)

        # Test discretization with angle resolution
        points = bspline_curve2d.discretization_points(angle_resolution=10)
        self.assertEqual(len(points), 31)


if __name__ == '__main__':
    unittest.main()
