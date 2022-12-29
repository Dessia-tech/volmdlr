"""
Unit tests for volmdlr.faces.BSplineCurve
"""
import unittest
from geomdl import BSpline

from volmdlr.models import bspline_curves
import volmdlr
import volmdlr.edges as vme


class TestBSplineCurve(unittest.TestCase):

    def test_abscissa(self):
        bspline_curve2d = bspline_curves.bspline_curve2d_1
        point = volmdlr.Point2D(-0.31240117104573617, -2.8555856978321796)

        self.assertAlmostEqual(bspline_curve2d.abscissa(point), 7.747599410268476)

    def test_discretization_points(self):
        control_points_2d = [volmdlr.Point2D(1.5707963267948966, 2.3),
                             volmdlr.Point2D(1.680890866936472, 2.256043878001211),
                             volmdlr.Point2D(1.8428579918488803, 2.190912791233705),
                             volmdlr.Point2D(2.0551351923128847, 2.110710771857296),
                             volmdlr.Point2D(2.2068399827060317, 2.057538514554844),
                             volmdlr.Point2D(2.3561943231153806, 2.010935033351481),
                             volmdlr.Point2D(2.505548683644506, 1.9715519259143607),
                             volmdlr.Point2D(2.65725353031637, 1.940017133765504),
                             volmdlr.Point2D(2.8695307222689292, 1.908674758526091),
                             volmdlr.Point2D(3.031498051508191, 1.89997293414679),
                             volmdlr.Point2D(3.141592653589793, 1.9000000000000003)]
        bspline_curve2d = vme.BSplineCurve2D(3, control_points_2d, [4, 1, 1, 1, 1, 1, 1, 1, 4],
                                             [0.0, 0.2102659043588606, 0.30933566258662554, 0.40542083024287023,
                                              0.5000013075051806, 0.5945816603424732, 0.6906664654007513,
                                              0.7897356531977031, 1.0])

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
