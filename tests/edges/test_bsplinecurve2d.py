import unittest
import os
from dessia_common.core import DessiaObject
import volmdlr
import volmdlr.nurbs.helpers as nurbs_helpers
import volmdlr.edges as vme
from volmdlr.models import bspline_curves
from volmdlr import curves
from geomdl import BSpline


DELTA = 0.001
folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'bsplinecurve_objects')


class TestBSplineCurve2D(unittest.TestCase):
    def test_bounding_rectangle(self):
        contour = DessiaObject.from_json(os.path.join(folder, "bounding_box_contour.json"))
        b_rec = contour.bounding_rectangle
        self.assertAlmostEqual(b_rec.area(), 0.48129011002687494, places=2)

    degree = 3
    points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1), volmdlr.Point2D(2, -1), volmdlr.Point2D(3, 0)]
    knotvector = nurbs_helpers.generate_knot_vector(degree, len(points))
    knot_multiplicity = [1] * len(knotvector)
    bspline1 = vme.BSplineCurve2D(degree, points, knot_multiplicity, knotvector, None)
    bspline2, bspline3 = bspline1.split(volmdlr.Point2D(1.5, 0.0))
    bspline4, bspline5 = bspline2.split(bspline2.point_at_abscissa(0.3 * bspline2.length()))
    bspline6 = bspline1.split(bspline1.point_at_abscissa(0.7 * bspline1.length()))[0]
    bspline7 = bspline1.split(bspline1.point_at_abscissa(0.3 * bspline1.length()))[1]
    degree = 3
    ctrlpts = [
        volmdlr.Point2D(5.0, 5.0),
        volmdlr.Point2D(10.0, 10.0),
        volmdlr.Point2D(20.0, 15.0),
        volmdlr.Point2D(35.0, 15.0),
        volmdlr.Point2D(45.0, 10.0),
        volmdlr.Point2D(50.0, 5.0),
    ]
    knot_multiplicities = [4, 1, 1, 4]
    knots = [0.0, 0.33, 0.66, 1.0]

    bspline2d = vme.BSplineCurve2D(degree, ctrlpts, knot_multiplicities, knots)
    weights = [0.5, 1.0, 0.75, 1.0, 0.25, 1.0]
    bspline2d_rational = vme.BSplineCurve2D(degree, ctrlpts, knot_multiplicities, knots, weights=weights)

    def test_evaluate_single(self):
        test_cases = [
            (0.0, (5.0, 5.0)),
            (0.3, (18.617, 13.377)),
            (0.5, (27.645, 14.691)),
            (0.6, (32.143, 14.328)),
            (1.0, (50.0, 5.0)),
        ]
        for param, res in test_cases:
            with self.subTest(param=param):
                evalpt = self.bspline2d.evaluate_single(param)
                self.assertAlmostEqual(evalpt[0], res[0], delta=DELTA)
                self.assertAlmostEqual(evalpt[1], res[1], delta=DELTA)

        test_cases = [
            (0.0, (5.0, 5.0)),
            (0.2, (13.8181, 11.5103)),
            (0.5, (28.1775, 14.7858)),
            (0.95, (48.7837, 6.0022)),
        ]
        for param, res in test_cases:
            with self.subTest(param=param):
                evalpt = self.bspline2d_rational.evaluate_single(param)

                self.assertAlmostEqual(evalpt[0], res[0], delta=DELTA)
                self.assertAlmostEqual(evalpt[1], res[1], delta=DELTA)

    def test_derivatives(self):
        derivatives = self.bspline2d.derivatives(u=0.35, order=2)
        expected_result = [[20.879272837543425, 13.96350686701158], [45.20015428165102, 9.987462558623653], [-1.334434093851499, -68.74685708529317]]
        for der, res in zip(derivatives, expected_result):
            self.assertAlmostEqual(der[0], res[0], delta=DELTA)
            self.assertAlmostEqual(der[1], res[1], delta=DELTA)

        test_cases = [
            (0.0, 1, ((5.0, 5.0), (90.9090, 90.9090))),
            (0.2, 2, ((13.8181, 11.5103), (40.0602, 17.3878), (104.4062, -29.3672))),
            (0.5, 3, ((28.1775, 14.7858), (39.7272, 2.2562), (-116.9254, -49.7367), (125.5276, 196.8865))),
            (0.95, 1, ((48.7837, 6.0022), (39.5178, -29.9962))),
        ]
        for param, order, res in test_cases:
            deriv = self.bspline2d_rational.derivatives(u=param, order=order)

            for computed, expected in zip(deriv, res):
                for c, e in zip(computed, expected):
                    self.assertAlmostEqual(c, e, delta=DELTA)

    def test_interpolate_curve(self):
        # The NURBS Book Ex9.1
        points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(3, 4), volmdlr.Point2D(-1, 4),
                  volmdlr.Point2D(-4, 0), volmdlr.Point2D(-4, -3)]
        degree = 3  # cubic curve

        # Do global curve interpolation
        curve = vme.BSplineCurve2D.from_points_interpolation(points, degree, centripetal=False)
        expected_ctrlpts = [
            [0.0, 0.0],
            [7.3169635171119936, 3.6867775257587367],
            [-2.958130565851424, 6.678276528176592],
            [-4.494953466891109, -0.6736915062424752],
            [-4.0, -3.0],
        ]
        for point, expected_point in zip(curve.control_points, expected_ctrlpts):
            self.assertAlmostEqual(point[0], expected_point[0], delta=DELTA)
            self.assertAlmostEqual(point[1], expected_point[1], delta=DELTA)

    def test_approximate_curve(self):
        # The NURBS Book Ex9.1
        points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(3, 4), volmdlr.Point2D(-1, 4),
                  volmdlr.Point2D(-4, 0), volmdlr.Point2D(-4, -3)]
        degree = 3  # cubic curve

        # Do global curve interpolation
        curve = vme.BSplineCurve2D.from_points_approximation(points, degree, centripetal=False)
        expected_ctrlpts = [
            [0.0, 0.0],
            [9.610024470158852, 8.200277881464892],
            [-8.160625855418692, 3.3820642030608417],
            [-4.0, -3.0],
        ]
        for point, expected_point in zip(curve.control_points, expected_ctrlpts):
            self.assertAlmostEqual(point[0], expected_point[0], delta=DELTA)
            self.assertAlmostEqual(point[1], expected_point[1], delta=DELTA)

        points2d = [volmdlr.Point2D(0, 0.1),
                    volmdlr.Point2D(0.2, 0.3),
                    volmdlr.Point2D(0.4, 0.4),
                    volmdlr.Point2D(0.5, 0.6),
                    volmdlr.Point2D(0.6, 0.7),
                    volmdlr.Point2D(0.8, 0.8),
                    volmdlr.Point2D(1, 0.9)]

        # %%% Approximation
        bspline_curve2d_approximated = vme.BSplineCurve2D.from_points_approximation(points2d, 3, ctrlpts_size=5)
        expected_ctrlpts = [volmdlr.Point2D(0.0, 0.1), volmdlr.Point2D(0.1686778402310228, 0.2366540266279785),
                            volmdlr.Point2D(0.466545895266623, 0.5077440536607246),
                            volmdlr.Point2D(0.7432185866086097, 0.852531277025759), volmdlr.Point2D(1.0, 0.9)]
        for point, expected_point in zip(bspline_curve2d_approximated.control_points, expected_ctrlpts):
            self.assertAlmostEqual(point[0], expected_point[0], delta=DELTA)
            self.assertAlmostEqual(point[1], expected_point[1], delta=DELTA)

    def test_length(self):
        total_length = self.bspline2d.length()
        self.assertAlmostEqual(total_length, 50.33433959792692, delta=DELTA)

    def test_abscissa(self):
        bspline_curve2d = bspline_curves.bspline_curve2d_1
        point = volmdlr.Point2D(-0.31240117104573617, -2.8555856978321796)

        bspline = vme.BSplineCurve2D.from_json(os.path.join(folder, "bg_bspline5_.json"))
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
        bspline_curve2d = DessiaObject.from_json(os.path.join(folder, "bsplinecurve2d_1.json"))
        line = curves.Line2D(volmdlr.Point2D(1.263163105753452, -0.002645572020392778),
                          volmdlr.Point2D(1.263163105753452, -0.001820963841291406))

        line_intersections = bspline_curve2d.line_intersections(line)
        self.assertEqual(len(line_intersections), 1)
        self.assertTrue(line_intersections[0].is_close(volmdlr.Point2D(1.263163105753452, -0.0026450893856384914)))

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
        offseted_bspline = self.bspline1.offset(-0.2)
        expected_distances = [0.2, 0.20000160183808904, 0.20053651951715856, 0.20372910969690097, 0.210441370708919,
                              0.2192581584663399, 0.22774528008118392, 0.23404460706854788, 0.23739001591364056,
                              0.2379018126594174, 0.2362014374337063, 0.23307773295678147, 0.22924032294583793,
                              0.22517329538697972, 0.22109005047384114, 0.21697594011450796, 0.21267059325565962,
                              0.2079610665048543, 0.20299372351359257, 0.19999999999999987]
        for i, (point1, point2) in enumerate(zip(self.bspline1.discretization_points(number_points=20),
                                                 offseted_bspline.discretization_points(number_points=20))):
            self.assertAlmostEqual(point1.point_distance(point2), expected_distances[i], 5)

    def test_point_distance(self):
        point = volmdlr.Point2D(1.5, 0.1)
        self.assertAlmostEqual(self.bspline1.point_distance(point), 0.08945546033235202)
        point2 = self.bspline1.point_at_abscissa(0.4)
        self.assertAlmostEqual(self.bspline1.point_distance(point2), 0.0, 7)

    def test_point_belongs(self):
        point = volmdlr.Point2D(1.5, 0.1)
        self.assertFalse(self.bspline1.point_belongs(point))
        point2 = self.bspline1.point_at_abscissa(0.4)
        self.assertTrue(self.bspline1.point_belongs(point2))

    def test_get_shared_primitives(self):
        shared_section1 = self.bspline1.get_shared_section(self.bspline2)
        self.assertEqual(len(shared_section1), 1)
        self.assertTrue(shared_section1[0].start.is_close(volmdlr.Point2D(0.0, 0.0)))
        self.assertTrue(shared_section1[0].end.is_close(volmdlr.Point2D(1.5, 0.0)))
        shared_section2 = self.bspline6.get_shared_section(self.bspline7)
        self.assertEqual(len(shared_section2), 1)
        self.assertTrue(shared_section2[0].start.is_close(volmdlr.Point2D(0.8999999, 0.252000000)))
        self.assertTrue(shared_section2[0].end.is_close(volmdlr.Point2D(2.09999999, -0.251999999)))
        self.assertAlmostEqual(shared_section2[0].length(), 1.3039875674329982, 6)
        shared_section3 = self.bspline1.get_shared_section(self.bspline5)
        self.assertEqual(shared_section3, [self.bspline5])
        shared_section4 = self.bspline5.get_shared_section(self.bspline1)
        self.assertEqual(shared_section4, [self.bspline5])
        self.assertFalse(self.bspline4.get_shared_section(self.bspline3))

    def test_delete_shared_primitives(self):
        remaining_section1 = self.bspline1.delete_shared_section(self.bspline2)
        self.assertEqual(len(remaining_section1), 1)
        self.assertTrue(remaining_section1[0].start.is_close(volmdlr.Point2D(1.5, 0.0)))
        self.assertTrue(remaining_section1[0].end.is_close(volmdlr.Point2D(3.0, 0.0)))
        self.assertAlmostEqual(remaining_section1[0].length(), 1.6373881438050524, 6)
        remaining_section2 = self.bspline6.delete_shared_section(self.bspline7)
        self.assertEqual(len(remaining_section2), 1)
        self.assertTrue(remaining_section2[0].start.is_close(volmdlr.Point2D(0.0, 0.0)))
        self.assertTrue(remaining_section2[0].end.is_close(volmdlr.Point2D(0.8999999997498065, 0.25200000006505024)))
        self.assertAlmostEqual(remaining_section2[0].length(), 0.9854029549808058, 6)
        remaining_section3 = self.bspline1.delete_shared_section(self.bspline5)
        self.assertEqual(len(remaining_section3), 2)
        self.assertTrue(remaining_section3[0].start.is_close(volmdlr.Point2D(0.0, 0.0)))
        self.assertTrue(remaining_section3[0].end.is_close(volmdlr.Point2D(0.44999999682593295, 0.26774999925409426)))
        self.assertAlmostEqual(remaining_section3[0].length(), 0.5305607215935024, 6)
        self.assertTrue(remaining_section3[1].start.is_close(volmdlr.Point2D(1.4999999878769186, 0.0)))
        self.assertTrue(remaining_section3[1].end.is_close(volmdlr.Point2D(3.0, 0.0)))
        self.assertAlmostEqual(remaining_section3[1].length(), 1.6373881438050524, 6)
        self.assertFalse(self.bspline5.delete_shared_section(self.bspline1))
        remaining_section4 = self.bspline4.delete_shared_section(self.bspline3)
        self.assertEqual(remaining_section4, [self.bspline4])

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
        point1 = self.bspline1.point_at_abscissa(0.25)
        point2 = self.bspline1.point_at_abscissa(0.65)
        local_discretization = self.bspline1.local_discretization(point1, point2, 10)
        for point1, point2 in zip(expected_points, local_discretization):
            self.assertTrue(point1.is_close(point2))

    def test_simplify(self):
        bsplinecurve = vme.BSplineCurve3D.from_json(os.path.join(folder, "bsplinecurve_fullarc.json"))
        fullarc = bsplinecurve.simplify
        self.assertTrue(isinstance(fullarc, vme.FullArc3D))

    def test_direction_independent_is_close(self):
        bsplinecurve1 = vme.BSplineCurve3D.from_json(os.path.join(folder, "bspline_curve1.json"))
        bsplinecurve2 = vme.BSplineCurve3D.from_json(os.path.join(folder, "bspline_curve2.json"))
        self.assertTrue(bsplinecurve1.direction_independent_is_close(bsplinecurve2))

    def test_split_curve(self):
        split_point = volmdlr.Point2D(28.1775252667145, 14.785855215217019)
        splitted_curves = self.bspline2d_rational.split(split_point)
        self.assertTrue(splitted_curves[0].start.is_close(volmdlr.Point2D(5.0, 5.0)))
        self.assertTrue(splitted_curves[0].end.is_close(split_point))
        self.assertTrue(splitted_curves[1].start.is_close(split_point))
        self.assertTrue(splitted_curves[1].end.is_close(volmdlr.Point2D(50.0, 5.0)))

        split_point = volmdlr.Point2D(0.04820589473987067, 0.011936395549382077)
        bsplinecurve = vme.BSplineCurve2D.from_json(os.path.join(folder, "bsplinecurve_split_bug.json"))
        splitted_curves = bsplinecurve.split(split_point)
        self.assertTrue(splitted_curves[0].start.is_close(volmdlr.Point2D(0.04873977000999985, 0.011815456390639745)))
        self.assertTrue(splitted_curves[0].end.is_close(split_point))
        self.assertTrue(splitted_curves[1].start.is_close(split_point))
        self.assertTrue(splitted_curves[1].end.is_close(volmdlr.Point2D(0.04793931370999993, 0.011891887758212483)))
        self.assertAlmostEqual(splitted_curves[0].length(), 0.0005535177002044544, 5)
        self.assertAlmostEqual(splitted_curves[1].length(), 0.0002710315376536523, 5)

    def test_merge_with(self):
        split_point = volmdlr.Point2D(27.64549230676716, 14.691702224146088)
        splitted_curves = self.bspline2d.split(split_point)
        merged_curve = splitted_curves[0].merge_with(splitted_curves[1])
        self.assertTrue(merged_curve.is_close(self.bspline2d))
        self.assertFalse(merged_curve.rational)

        split_point = volmdlr.Point2D(28.1775252667145, 14.785855215217019)
        splitted_curves = self.bspline2d_rational.split(split_point)
        merged_curve = splitted_curves[0].merge_with(splitted_curves[1])
        self.assertTrue(merged_curve.is_close(self.bspline2d_rational))

    def test_tangent(self):
        tangent = self.bspline1.tangent(0.5)
        tangent_rational = self.bspline2d_rational.tangent(0.5)
        self.assertTrue(tangent.is_close(volmdlr.Vector2D(0.8944271909999159, -0.4472135954999579)))
        self.assertTrue(tangent_rational.is_close(volmdlr.Vector2D(0.998391126712381, 0.05670236416572517)))


if __name__ == '__main__':
    unittest.main()
