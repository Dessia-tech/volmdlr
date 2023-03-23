"""
Unit tests for volmdlr.faces.BSplineCurve
"""
import unittest

from dessia_common.core import DessiaObject
from geomdl import BSpline, utilities

import volmdlr
import volmdlr.edges as vme
from volmdlr.models import bspline_curves


class TestBSplineCurve2D(unittest.TestCase):
    degree = 3
    points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1), volmdlr.Point2D(2, -1), volmdlr.Point2D(3, 0)]
    knotvector = utilities.generate_knot_vector(degree, len(points))
    knot_multiplicity = [1] * len(knotvector)
    bspline = vme.BSplineCurve2D(degree, points, knot_multiplicity, knotvector, None, False)

    def test_abscissa(self):
        bspline_curve2d = bspline_curves.bspline_curve2d_1
        point = volmdlr.Point2D(-0.31240117104573617, -2.8555856978321796)

        bspline = vme.BSplineCurve2D.load_from_file("edges/bg_bspline5_.json")
        point1 = bspline.points[25]
        point2 = bspline.points[75]

        abscissa1 = bspline.abscissa(point1)
        abscissa2 = bspline.abscissa(point2)

        test_point1 = bspline.point_at_abscissa(abscissa1)
        test_point2 = bspline.point_at_abscissa(abscissa2)

        self.assertTrue(point1.is_close(test_point1))
        self.assertTrue(point2.is_close(test_point2))

        abscissa3 = 0.00016294494116532595
        abscissa4 = 0.00017682955170114393

        point_at_abscissa3 = bspline.point_at_abscissa(abscissa3)
        point_at_abscissa4 = bspline.point_at_abscissa(abscissa4)

        test_abscissa3 = bspline.abscissa(point_at_abscissa3)
        test_abscissa4 = bspline.abscissa(point_at_abscissa4)

        self.assertAlmostEqual(abscissa3, test_abscissa3, 6)
        self.assertAlmostEqual(abscissa4, test_abscissa4, 6)

        self.assertAlmostEqual(bspline_curve2d.abscissa(point), 7.747599410268476)

    def test_line_intersections(self):
        bspline_curve2d = DessiaObject.load_from_file('edges/bsplinecurve2d_1.json')
        line = vme.Line2D(volmdlr.Point2D(1.263163105753452, -0.002645572020392778),
                          volmdlr.Point2D(1.263163105753452, -0.001820963841291406))

        line_intersections = bspline_curve2d.line_intersections(line)
        self.assertEqual(len(line_intersections), 1)
        self.assertTrue(line_intersections[0].is_close(volmdlr.Point2D(1.2631631057526727, -0.0026450894385881708)))

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
        self.assertEqual(len(points), 100)

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

    def test_offset(self):
        offseted_bspline = self.bspline.offset(-0.2)
        expected_distances = [0.2, 0.20000160183808904, 0.20053651951715856, 0.20372900125730523, 0.21044118400720574,
                              0.2192581584663399, 0.22774528008118392, 0.2340440381875313, 0.23739001591364056,
                              0.2379018126594174, 0.2362014374337063, 0.23307773295678147, 0.22924032294583793,
                              0.22517329538697972, 0.22109005047384114, 0.21697594011450796, 0.21267059325565962,
                              0.2079610665048543, 0.20299372351359257, 0.19999999999999987]
        for i, (point1, point2) in enumerate(zip(self.bspline.discretization_points(number_points=20),
                                                 offseted_bspline.discretization_points(number_points=20))):
            self.assertAlmostEqual(point1.point_distance(point2), expected_distances[i], 5)

    def test_point_distance(self):
        point = volmdlr.Point2D(1.5, 0.1)
        self.assertAlmostEqual(self.bspline.point_distance(point), 0.08945546033235202)
        point2 = self.bspline.point_at_abscissa(0.4)
        self.assertAlmostEqual(self.bspline.point_distance(point2), 0.0, 8)

    def test_point_belongs(self):
        point = volmdlr.Point2D(1.5, 0.1)
        self.assertFalse(self.bspline.point_belongs(point))
        point2 = self.bspline.point_at_abscissa(0.4)
        self.assertTrue(self.bspline.point_belongs(point2))

    def test_local_discretization(self):
        expected_points = [volmdlr.Point2D(0.22902909156524637, 0.17924444819399216),
                           volmdlr.Point2D(0.26974537451069974, 0.2013444443084787),
                           volmdlr.Point2D(0.3104616574561531, 0.22072505985054805),
                           volmdlr.Point2D(0.35117794040160644, 0.23747629494418182),
                           volmdlr.Point2D(0.3918942233470598, 0.25168814971336145),
                           volmdlr.Point2D(0.4326105062925132, 0.26345062428206867),
                           volmdlr.Point2D(0.4733267892379665, 0.2728537187742847),
                           volmdlr.Point2D(0.5140430721834197, 0.27998743331399134),
                           volmdlr.Point2D(0.5547593551288732, 0.28494176802517024),
                           volmdlr.Point2D(0.5954756380743265, 0.28780672303180266)]
        point1 = self.bspline.point_at_abscissa(0.25)
        point2 = self.bspline.point_at_abscissa(0.65)
        local_discretization = self.bspline.local_discretization(point1, point2, 10)
        for point1, point2 in zip(expected_points, local_discretization):
            self.assertTrue(point1.is_close(point2))


class TestBSplineCurve3D(unittest.TestCase):
    b_splinecurve3d = vme.BSplineCurve3D(degree=5, control_points=[
        volmdlr.Point3D(0.5334, 4.61e-10, -2.266), volmdlr.Point3D(0.5334, 0.236642912449, -2.26599999893),
        volmdlr.Point3D(0.5334, 0.473285829931, -2.23144925183),
        volmdlr.Point3D(0.5334, 0.70316976404, -2.16234807551),
        volmdlr.Point3D(0.5334, 1.13611540546, -1.95904362568), volmdlr.Point3D(0.5334, 1.49286052971, -1.64044168585),
        volmdlr.Point3D(0.5334, 1.64654439419, -1.45604332404), volmdlr.Point3D(0.5334, 1.77109261028, -1.25188280667),
        volmdlr.Point3D(0.5334, 1.86385510975, -1.03417888209)], knot_multiplicities=[6, 3, 6],
                                         knots=[0.0, 0.4999999725155696, 1.0])

    def test_line_intersections(self):
        line = vme.Line3D(volmdlr.Point3D(0.5334, -0.44659009801843536, 0.0),
                          volmdlr.Point3D(0.5334, 0.4342689853571558, -0.47337857496375274))
        bspline_line_intersections = self.b_splinecurve3d.line_intersections(line)
        self.assertEqual(bspline_line_intersections, [volmdlr.Point3D(0.5334, 1.784620481894723, -1.1990650295776075)])

    def test_linesegment_intersection(self):
        linesegment1 = vme.LineSegment3D(volmdlr.Point3D(0.5334, -0.44659009801843536, 0.0),
                                         volmdlr.Point3D(0.5334, 0.4342689853571558, -0.47337857496375274))
        linesegment2 = vme.LineSegment3D(volmdlr.Point3D(0.5334, -0.44659009801843536, 0.0),
                                         volmdlr.Point3D(0.5334, 2.1959871521083385, -1.4201357248912583))
        bspline_lineseg_intersections1 = self.b_splinecurve3d.linesegment_intersections(linesegment1)
        bspline_lineseg_intersections2 = self.b_splinecurve3d.linesegment_intersections(linesegment2)
        self.assertFalse(bspline_lineseg_intersections1)
        self.assertTrue(bspline_lineseg_intersections2[0].is_close(
            volmdlr.Point3D(0.5334, 1.784620481894723, -1.1990650295776075)))


if __name__ == '__main__':
    unittest.main()
