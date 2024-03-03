"""
Unit tests for volmdlr.faces.BSplineSurface3D
"""
import unittest
import os
import numpy as np
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.faces as vmf
import volmdlr.grid
from volmdlr.models import bspline_surfaces
from volmdlr import surfaces


DELTA = 0.001
folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_bspline_test')


class TestBSplineSurface3D(unittest.TestCase):

    def setUp(self):
        """Create a B-spline surface instance as a fixture"""
        # Set degrees
        degree_u = 3
        degree_v = 3

        ctrlpts = [
            [-25.0, -25.0, -10.0], [-25.0, -15.0, -5.0], [-25.0, -5.0, 0.0], [-25.0, 5.0, 0.0],
            [-25.0, 15.0, -5.0], [-25.0, 25.0, -10.0], [-15.0, -25.0, -8.0], [-15.0, -15.0, -4.0],
            [-15.0, -5.0, -4.0], [-15.0, 5.0, -4.0], [-15.0, 15.0, -4.0], [-15.0, 25.0, -8.0],
            [-5.0, -25.0, -5.0], [-5.0, -15.0, -3.0], [-5.0, -5.0, -8.0], [-5.0, 5.0, -8.0],
            [-5.0, 15.0, -3.0], [-5.0, 25.0, -5.0], [5.0, -25.0, -3.0], [5.0, -15.0, -2.0],
            [5.0, -5.0, -8.0], [5.0, 5.0, -8.0], [5.0, 15.0, -2.0], [5.0, 25.0, -3.0],
            [15.0, -25.0, -8.0], [15.0, -15.0, -4.0], [15.0, -5.0, -4.0], [15.0, 5.0, -4.0],
            [15.0, 15.0, -4.0], [15.0, 25.0, -8.0], [25.0, -25.0, -10.0], [25.0, -15.0, -5.0],
            [25.0, -5.0, 2.0], [25.0, 5.0, 2.0], [25.0, 15.0, -5.0], [25.0, 25.0, -10.0],
        ]

        nb_u, nb_v = 6, 6

        # Set knot vectors
        knots_u = [0.0, 0.33, 0.66, 1.0]
        u_multiplicities = [4, 1, 1, 4]
        knots_v = [0.0, 0.33, 0.66, 1.0]
        v_multiplicities = [4, 1, 1, 4]
        # knotvector_u = [0.0, 0.0, 0.0, 0.0, 0.33, 0.66, 1.0, 1.0, 1.0, 1.0]
        # knotvector_v = [0.0, 0.0, 0.0, 0.0, 0.33, 0.66, 1.0, 1.0, 1.0, 1.0]

        self.spline_surf = surfaces.BSplineSurface3D(degree_u, degree_v, ctrlpts, nb_u, nb_v,
                                                     u_multiplicities, v_multiplicities, knots_u, knots_v)
        weights = [0.5, 1.0, 0.75, 1.0, 0.25, 1.0, 0.5, 1.0, 0.75, 1.0, 0.25, 1.0,
                   0.5, 1.0, 0.75, 1.0, 0.25, 1.0, 0.5, 1.0, 0.75, 1.0, 0.25, 1.0,
                   0.5, 1.0, 0.75, 1.0, 0.25, 1.0, 0.5, 1.0, 0.75, 1.0, 0.25, 1.0]
        self.nurbs_surf = surfaces.BSplineSurface3D(degree_u, degree_v, ctrlpts, nb_u, nb_v,
                                                     u_multiplicities, v_multiplicities, knots_u, knots_v, weights)

    def test_extract_curves(self):
        u = [0.0, 0.25, 0.5, 0.75, 1.0]
        v = [0.0, 0.25, 0.5, 0.75, 1.0]

        expected_results_u = [[volmdlr.Point3D(-25.0, -25.0, -10.0), volmdlr.Point3D(-25.0, 25.0, -10.0)],
                              [volmdlr.Point3D(-9.07716976931853, -25.0, -6.280643904610846),
                               volmdlr.Point3D(-9.07716976931853, 25.0, -6.280643904610846)],
                              [volmdlr.Point3D(0.11256465850708097, -25.0, -4.238138097577792),
                               volmdlr.Point3D(0.11256465850708097, 25.0, -4.238138097577792)],
                              [volmdlr.Point3D(9.309787132823084, -25.0, -5.579386459163334),
                               volmdlr.Point3D(9.309787132823084, 25.0, -5.579386459163334)],
                              [volmdlr.Point3D(25.0, -25.0, -10.0), volmdlr.Point3D(25.0, 25.0, -10.0)]
                              ]
        expected_results_v = [[volmdlr.Point3D(-25.0, -25.0, -10.0), volmdlr.Point3D(25.0, -25.0, -10.0)],
                              [volmdlr.Point3D(-25.0, -9.07716976931853, -2.397285527450817),
                               volmdlr.Point3D(25.0, -9.07716976931853, -1.3277054289450985)],
                              [volmdlr.Point3D(-25.0, 0.11256465850708097, -0.308297775853914),
                               volmdlr.Point3D(25.0, 0.11256465850708097, 1.5683831138045203)],
                              [volmdlr.Point3D(-25.0, 9.309787132823084, -2.497847912329014),
                               volmdlr.Point3D(25.0, 9.309787132823084, -1.4598916162388402)],
                              [volmdlr.Point3D(-25.0, 25.0, -10.0), volmdlr.Point3D(25.0, 25.0, -10.0)]
                              ]
        test = self.spline_surf.extract_curves(u, v)
        for curve, expected_result in zip(test["u"], expected_results_u):
            control_points = curve.control_points
            self.assertEqual(len(control_points), 6)
            self.assertTrue(curve.start.is_close(expected_result[0]))
            self.assertTrue(curve.end.is_close(expected_result[1]))

        for curve, expected_result in zip(test["v"], expected_results_v):
            control_points = curve.control_points
            self.assertEqual(len(control_points), 6)
            self.assertTrue(curve.start.is_close(expected_result[0]))
            self.assertTrue(curve.end.is_close(expected_result[1]))

    def test_extract_curves_rational(self):
        u = [0.0, 0.25, 0.5, 0.75, 1.0]
        v = [0.0, 0.25, 0.5, 0.75, 1.0]

        # points from surface evaluation
        expected_results_u = [[volmdlr.Point3D(-25.0, -25.0, -10.0),
                               volmdlr.Point3D(-25.0, -9.485960183093134, -2.651935673517623),
                               volmdlr.Point3D(-25.0, 0.3732766669129052, -0.21414478478298477),
                               volmdlr.Point3D(-25.000000000000004, 6.5923170800941175, -1.200587938543581),
                               volmdlr.Point3D(-25.0, 25.0, -10.0)],
                              [volmdlr.Point3D(-9.07716976931853, -25.000000000000004, -6.280643904610847),
                               volmdlr.Point3D(-9.07716976931853, -9.485960183093134, -4.7087289345950065),
                               volmdlr.Point3D(-9.07716976931853, 0.3732766669129052, -5.9676274770461415),
                               volmdlr.Point3D(-9.077169769318532, 6.5923170800941175, -5.601767665101728),
                               volmdlr.Point3D(-9.07716976931853, 25.000000000000004, -6.280643904610847)],
                              [volmdlr.Point3D(0.11256465850708161, -24.999999999999996, -4.238138097577793),
                               volmdlr.Point3D(0.1125646585070812, -9.485960183093134, -5.069196197307746),
                               volmdlr.Point3D(0.11256465850708078, 0.3732766669129052, -7.532144969741332),
                               volmdlr.Point3D(0.11256465850708117, 6.592317080094115, -6.7118731933170785),
                               volmdlr.Point3D(0.11256465850708161, 24.999999999999996, -4.238138097577793)],
                              [volmdlr.Point3D(9.309787132823084, -24.999999999999996, -5.579386459163333),
                               volmdlr.Point3D(9.309787132823084, -9.485960183093132, -4.462539699464645),
                               volmdlr.Point3D(9.309787132823082, 0.3732766669129051, -5.839760787576515),
                               volmdlr.Point3D(9.309787132823086, 6.592317080094118, -5.4233896471890475),
                               volmdlr.Point3D(9.309787132823084, 24.999999999999996, -5.579386459163333)],
                              [volmdlr.Point3D(25.0, -25.0, -10.0),
                               volmdlr.Point3D(25.0, -9.485960183093134, -1.6964667229125177),
                               volmdlr.Point3D(25.0, 0.3732766669129052, 1.7001973013038214),
                               volmdlr.Point3D(25.000000000000004, 6.5923170800941175, 0.37750338602002487),
                               volmdlr.Point3D(25.0, 25.0, -10.0)]
                              ]

        # points from surface evaluation
        expected_results_v = [[volmdlr.Point3D(-25.0, -25.0, -10.0),
                               volmdlr.Point3D(-9.07716976931853, -25.000000000000004, -6.280643904610847),
                               volmdlr.Point3D(0.11256465850708161, -24.999999999999996, -4.238138097577793),
                               volmdlr.Point3D(9.309787132823084, -24.999999999999996, -5.579386459163333),
                               volmdlr.Point3D(25.0, -25.0, -10.0)],
                              [volmdlr.Point3D(-25.0, -9.485960183093134, -2.651935673517623),
                               volmdlr.Point3D(-9.07716976931853, -9.485960183093134, -4.7087289345950065),
                               volmdlr.Point3D(0.1125646585070812, -9.485960183093134, -5.069196197307746),
                               volmdlr.Point3D(9.309787132823084, -9.485960183093132, -4.462539699464645),
                               volmdlr.Point3D(25.0, -9.485960183093134, -1.6964667229125177)],
                              [volmdlr.Point3D(-25.0, 0.3732766669129052, -0.21414478478298477),
                               volmdlr.Point3D(-9.07716976931853, 0.3732766669129052, -5.9676274770461415),
                               volmdlr.Point3D(0.11256465850708078, 0.3732766669129052, -7.532144969741332),
                               volmdlr.Point3D(9.309787132823082, 0.3732766669129051, -5.839760787576515),
                               volmdlr.Point3D(25.0, 0.3732766669129052, 1.7001973013038214)],
                              [volmdlr.Point3D(-25.000000000000004, 6.5923170800941175, -1.200587938543581),
                               volmdlr.Point3D(-9.077169769318532, 6.5923170800941175, -5.601767665101728),
                               volmdlr.Point3D(0.11256465850708117, 6.592317080094115, -6.7118731933170785),
                               volmdlr.Point3D(9.309787132823086, 6.592317080094118, -5.4233896471890475),
                               volmdlr.Point3D(25.000000000000004, 6.5923170800941175, 0.37750338602002487)],
                              [volmdlr.Point3D(-25.0, 25.0, -10.0),
                               volmdlr.Point3D(-9.07716976931853, 25.000000000000004, -6.280643904610847),
                               volmdlr.Point3D(0.11256465850708161, 24.999999999999996, -4.238138097577793),
                               volmdlr.Point3D(9.309787132823084, 24.999999999999996, -5.579386459163333),
                               volmdlr.Point3D(25.0, 25.0, -10.0)]]
        test = self.nurbs_surf.extract_curves(u, v)
        for curve, expected_result in zip(test["u"], expected_results_u):
            control_points = curve.control_points
            self.assertEqual(len(control_points), 6)
            for point in expected_result:
                self.assertTrue(curve.point_belongs(point))

        for curve, expected_result in zip(test["v"], expected_results_v):
            control_points = curve.control_points
            self.assertEqual(len(control_points), 6)
            for point in expected_result:
                self.assertTrue(curve.point_belongs(point))

    def test_point2d_to_3d(self):
        test_cases = [
            (volmdlr.Point2D(0.0, 0.0), (-25.0, -25.0, -10.0)),
            (volmdlr.Point2D(0.0, 0.2), (-25.0, -11.403, -3.385)),
            (volmdlr.Point2D(0.0, 1.0), (-25.0, 25.0, -10.0)),
            (volmdlr.Point2D(0.3, 0.0), (-7.006, -25.0, -5.725)),
            (volmdlr.Point2D(0.3, 0.4), (-7.006, -3.308, -6.265)),
            (volmdlr.Point2D(0.3, 1.0), (-7.006, 25.0, -5.725)),
            (volmdlr.Point2D(0.6, 0.0), (3.533, -25.0, -4.224)),
            (volmdlr.Point2D(0.6, 0.6), (3.533, 3.533, -6.801)),
            (volmdlr.Point2D(0.6, 1.0), (3.533, 25.0, -4.224)),
            (volmdlr.Point2D(1.0, 0.0), (25.0, -25.0, -10.0)),
            (volmdlr.Point2D(1.0, 0.8), (25.0, 11.636, -2.751)),
            (volmdlr.Point2D(1.0, 1.0), (25.0, 25.0, -10.0)),
        ]

        for param, res in test_cases:
            evalpt = self.spline_surf.point2d_to_3d(param)
            self.assertAlmostEqual(evalpt[0], res[0], delta=DELTA)
            self.assertAlmostEqual(evalpt[1], res[1], delta=DELTA)
            self.assertAlmostEqual(evalpt[2], res[2], delta=DELTA)

        test_data = [
            (volmdlr.Point2D(0.0, 0.0), (-25.0, -25.0, -10.0)),
            (volmdlr.Point2D(0.0, 0.2), (-25.0, -11.563, -3.489)),
            (volmdlr.Point2D(0.0, 1.0), (-25.0, 25.0, -10.0)),
            (volmdlr.Point2D(0.3, 0.0), (-7.006, -25.0, -5.725)),
            (volmdlr.Point2D(0.3, 0.4), (-7.006, -3.052, -6.196)),
            (volmdlr.Point2D(0.3, 1.0), (-7.006, 25.0, -5.725)),
            (volmdlr.Point2D(0.6, 0.0), (3.533, -25.0, -4.224)),
            (volmdlr.Point2D(0.6, 0.6), (3.533, 2.868, -7.257)),
            (volmdlr.Point2D(0.6, 1.0), (3.533, 25.0, -4.224)),
            (volmdlr.Point2D(1.0, 0.0), (25.0, -25.0, -10.0)),
            (volmdlr.Point2D(1.0, 0.8), (25.0, 9.425, -1.175)),
            (volmdlr.Point2D(1.0, 1.0), (25.0, 25.0, -10.0)),
        ]
        for param, res in test_data:
            evalpt = self.nurbs_surf.point2d_to_3d(param)
            self.assertAlmostEqual(evalpt[0], res[0], delta=DELTA)
            self.assertAlmostEqual(evalpt[1], res[1], delta=DELTA)
            self.assertAlmostEqual(evalpt[2], res[2], delta=DELTA)

    def test_parametric_points_to_3d(self):
        parametric_points = np.array([[0.0, 0.0], [0.0, 0.2], [0.0, 1.0], [0.3, 0.0],
                                      [0.3, 0.4], [0.3, 1.0], [0.6, 0.0], [0.6, 0.6],
                                      [0.6, 1.0], [1.0, 0.0], [1.0, 0.8], [1.0, 1.0],
                                      ])
        points3d = self.spline_surf.parametric_points_to_3d(parametric_points)
        expected_points = np.array([[-25.0, -25.0, -10.0], [-25.0, -11.40398, -3.3856], [-25.0, 25.0, -10.0],
                                    [-7.006, -25.0, -5.725], [-7.006, -3.308, -6.265], [-7.006, 25.0, -5.725],
                                    [3.533, -25.0, -4.224], [3.533, 3.533, -6.801], [3.533, 25.0, -4.224],
                                    [25.0, -25.0, -10.0], [25.0, 11.636, -2.751], [25.0, 25.0, -10.0]])
        for point, expected_point in zip(points3d, expected_points):
            self.assertAlmostEqual(np.linalg.norm(point - expected_point), 0.0, delta=DELTA)

        points3d = self.nurbs_surf.parametric_points_to_3d(parametric_points)
        expected_points = np.array([[-25.0, -25.0, -10.0], [-25.0, -11.563, -3.489], [-25.0, 25.0, -10.0],
                                    [-7.006, -25.0, -5.725], [-7.006, -3.052, -6.196], [-7.006, 25.0, -5.725],
                                    [3.533, -25.0, -4.224], [3.533, 2.868, -7.257], [3.533, 25.0, -4.224],
                                    [25.0, -25.0, -10.0], [25.0, 9.425, -1.175], [25.0, 25.0, -10.0]])
        for point, expected_point in zip(points3d, expected_points):
            self.assertAlmostEqual(np.linalg.norm(point - expected_point), 0.0, delta=DELTA)

    def test_derivatives(self):
        test_data = [
            (
                (0.0, 0.25),
                1,
                [
                    [[-25.0, -9.0771, -2.3972], [5.5511e-15, 43.6910, 17.5411]],
                    [[90.9090, 0.0, -15.0882], [-5.9750e-15, 0.0, -140.0367]],
                ],
            ),
            (
                (0.95, 0.75),
                2,
                [
                    [
                        [20.8948, 9.3097, -2.4845],
                        [-1.1347e-14, 43.7672, -15.0153],
                        [-5.0393e-30, 100.1022, -74.1165],
                    ],
                    [
                        [76.2308, -1.6965e-15, 18.0372],
                        [9.8212e-15, -5.9448e-15, -158.5462],
                        [4.3615e-30, -2.4356e-13, -284.3037],
                    ],
                    [
                        [224.5342, -5.6794e-14, 93.3843],
                        [4.9856e-14, -4.0400e-13, -542.6274],
                        [2.2140e-29, -1.88662e-12, -318.8808],
                    ],
                ],
            ),
        ]
        for param, order, res in test_data:
            deriv = self.spline_surf.derivatives(*param, order=order)
            for computed, expected in zip(deriv, res):
                for idx in range(order + 1):
                    for c, e in zip(computed[idx], expected[idx]):
                        self.assertAlmostEqual(c, e, delta=DELTA)

        test_data = [
            (
                (0.0, 0.25),
                1,
                [
                    [[-25.0, -9.4859, -2.6519], [9.1135e-15, 42.3855, 16.1635]],
                    [[90.9090, 0.0, -12.5504], [-3.6454e-14, 0.0, -135.9542]]
                ],
            ),
            (
                (0.95, 0.75),
                2,
                [
                    [
                        [20.8948, 6.5923, -1.4087],
                        [0.0, 42.0924, -14.6714],
                        [-4.4688e-14, 498.0982, -230.2790]
                    ],
                    [
                        [76.2308, -3.1813e-15, 31.9560],
                        [2.6692e-14, -1.6062e-14, -134.4819],
                        [-4.6035e-14, -4.5596e-13, -1646.1763]
                    ],
                    [
                        [224.5342, -1.9181e-14, 144.3825],
                        [-8.6754e-14, -4.3593e-13, -433.4414],
                        [1.2012e-12, -5.9424e-12, -4603.0856]
                    ]
                ],
            ),
        ]
        for param, order, res in test_data:
            deriv = self.nurbs_surf.derivatives(*param, order=order)
            for computed, expected in zip(deriv, res):
                for idx in range(order + 1):
                    for c, e in zip(computed[idx], expected[idx]):
                        self.assertAlmostEqual(c, e, delta=DELTA)

        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinesurface_derivatives_v_degree_1.json"))
        test_data = [((0.0, 0.41570203515189436), 2,
                     [[(2.686370456553301, -0.6157625276711683, 0.5584759391609816),
                       (0.04588611491000005, 0.02975096368800001, 0.05838959340700001),
                       (0.0, 0.0, 0.0)],
                      [(0.0, 0.11254942676189417, -0.05734675845308533),
                       (0.0, 0.06422355115999956, -0.03272351170500043),
                       (0.0, 0.0, 0.0)],
                      [(0.0, -0.06309176181313551, -0.12378660622984707),
                       (0.0, -0.035946467299996954, -0.07066413869999799), (0.0, 0.0, 0.0)]])]
        for param, order, res in test_data:
            deriv = surface.derivatives(*param, order=order)
            for computed, expected in zip(deriv, res):
                for idx in range(order + 1):
                    for c, e in zip(computed[idx], expected[idx]):
                        self.assertAlmostEqual(c, e, delta=DELTA)


    def test_interpolate_surface(self):
        points = [volmdlr.Point3D(1.0, 0.0, 0.0), volmdlr.Point3D(0.70710678, 0.70710678, 0.0),
                  volmdlr.Point3D(0.0, 1.0, 0.0), volmdlr.Point3D(-0.70710678, 0.70710678, 0.0),
                  volmdlr.Point3D(-1.0, 0.0, 0.0), volmdlr.Point3D(-0.70710678, -0.70710678, 0.0),
                  volmdlr.Point3D(0.0, -1.0, 0.0), volmdlr.Point3D(0.70710678, -0.70710678, 0.0),
                  volmdlr.Point3D(1.0, 0.0, 0.0), volmdlr.Point3D(1.0, 0.0, 1.0),
                  volmdlr.Point3D(0.70710678, 0.70710678, 1.0), volmdlr.Point3D(0.0, 1.0, 1.0),
                  volmdlr.Point3D(-0.70710678, 0.70710678, 1.0), volmdlr.Point3D(-1.0, 0.0, 1.0),
                  volmdlr.Point3D(-0.70710678, -0.70710678, 1.0), volmdlr.Point3D(0.0, -1.0, 1.0),
                  volmdlr.Point3D(0.70710678, -0.70710678, 1.0), volmdlr.Point3D(1.0, 0.0, 1.0)]

        degree_u = 1
        degree_v = 2
        size_u = 2
        size_v = 9
        surface = surfaces.BSplineSurface3D.from_points_interpolation(points, size_u, size_v,
                                                                                degree_u, degree_v)

        expected_points = [volmdlr.Point3D(1.0, 0.0, 0.0),
                           volmdlr.Point3D(0.9580995893491125, 0.6733882798117155, 0.0),
                           volmdlr.Point3D(-0.0005819501479292128, 1.0804111054393308, 0.0),
                           volmdlr.Point3D(-0.7628715805621301, 0.7627405224267781, 0.0),
                           volmdlr.Point3D(-1.0790428064792899, 0.0, 0.0),
                           volmdlr.Point3D(-0.7628715805621301, -0.7627405224267783, 0.0),
                           volmdlr.Point3D(-0.0005819501479290552, -1.0804111054393304, 0.0),
                           volmdlr.Point3D(0.9580995893491127, -0.6733882798117156, 0.0),
                           volmdlr.Point3D(1.0, 0.0, 0.0),
                           volmdlr.Point3D(1.0, 0.0, 1.0),
                           volmdlr.Point3D(0.9580995893491125, 0.6733882798117155, 1.0),
                           volmdlr.Point3D(-0.0005819501479292128, 1.0804111054393308, 1.0),
                           volmdlr.Point3D(-0.7628715805621301, 0.7627405224267781, 1.0),
                           volmdlr.Point3D(-1.0790428064792899, 0.0, 1.0),
                           volmdlr.Point3D(-0.7628715805621301, -0.7627405224267783, 1.0),
                           volmdlr.Point3D(-0.0005819501479290552, -1.0804111054393304, 1.0),
                           volmdlr.Point3D(0.9580995893491127, -0.6733882798117156, 1.0),
                           volmdlr.Point3D(1.0, 0.0, 1.0)]
        for point, expected_point in zip(surface.control_points, expected_points):
            self.assertTrue(point.is_close(expected_point))
        point1 = surface.point2d_to_3d(volmdlr.Point2D(0.0, 0.0))
        point2 = surface.point2d_to_3d(volmdlr.Point2D(0.25, 0.25))
        point3 = surface.point2d_to_3d(volmdlr.Point2D(0.5, 0.5))
        point4 = surface.point2d_to_3d(volmdlr.Point2D(0.75, 0.75))
        point5 = surface.point2d_to_3d(volmdlr.Point2D(1.0, 1.0))

        for point in [point1, point2, point3, point4, point5]:
            self.assertAlmostEqual(point.point_distance(volmdlr.Point3D(0.0, 0.0, point.z)), 1.0)

    def test_approximation_surface(self):
        construction_linesegment1 = vme.LineSegment3D(volmdlr.Point3D(0.5, -0.5, 0), volmdlr.Point3D(0.1, -0.5, 0))
        construction_linesegment2 = vme.LineSegment3D(volmdlr.Point3D(-0.5, -0.5, 0.5),
                                                      volmdlr.Point3D(-1.2, -0.5, 0.5))
        points = (construction_linesegment1.discretization_points(number_points=4) +
                  construction_linesegment2.discretization_points(number_points=5))
        for i in range(1, 4):
            temp_line1 = construction_linesegment1.translation(volmdlr.Vector3D(0, i*0.1, 0))
            temp_line2 = construction_linesegment2.translation(volmdlr.Vector3D(0, i * 0.1, 0))
            points.extend(temp_line1.discretization_points(number_points=4))
            points.extend(temp_line2.discretization_points(number_points=5))
        for i in range(7, 10):
            temp_line1 = construction_linesegment1.translation(volmdlr.Vector3D(0, i*0.1, -0.5))
            temp_line2 = construction_linesegment2.translation(volmdlr.Vector3D(0, i * 0.1, -0.5))
            points.extend(temp_line1.discretization_points(number_points=4))
            points.extend(temp_line2.discretization_points(number_points=5))

        degree_u = 2
        degree_v = 3
        size_u = 7
        size_v = 9
        surface = surfaces.BSplineSurface3D.from_points_approximation(points, size_u, size_v,
                                                                                    degree_u, degree_v,
                                                                                    ctrlpts_size_u=5, ctrlpts_size_v=6)

        expected_points = [
                           volmdlr.Point3D(0.5, -0.5, 0.0),
                           volmdlr.Point3D(0.4177877608688984, -0.49999999999999895, 0.018441240030724965),
                           volmdlr.Point3D(-0.08774883457241012, -0.5000000000000024, -0.12940652269371028),
                           volmdlr.Point3D(-0.4016938826094424, -0.4999999999999982, 0.7352898743278582),
                           volmdlr.Point3D(-1.046377148047132, -0.5000000000000008, 0.4409552496220567),
                           volmdlr.Point3D(-1.2, -0.5, 0.5),
                           volmdlr.Point3D(0.5000000000000002, -0.4376959060246239, -0.009076877140464948),
                           volmdlr.Point3D(0.41778776086889857, -0.4376959060246231, 0.009364362890260018),
                           volmdlr.Point3D(-0.08774883457241027, -0.4376959060246253, -0.13848339983417526),
                           volmdlr.Point3D(-0.4016938826094421, -0.43769590602462305, 0.7262129971873936),
                           volmdlr.Point3D(-1.046377148047133, -0.4376959060246242, 0.43187837248159167),
                           volmdlr.Point3D(-1.2000000000000002, -0.4376959060246239, 0.4909231228595352),
                           volmdlr.Point3D(0.49999999999999967, -0.04693124471297373, 0.04778074658303208),
                           volmdlr.Point3D(0.417787760868898, -0.04693124471297371, 0.06622198661375721),
                           volmdlr.Point3D(-0.08774883457240966, -0.04693124471297377, -0.08162577611067864),
                           volmdlr.Point3D(-0.4016938826094423, -0.04693124471297371, 0.7830706209108907),
                           volmdlr.Point3D(-1.0463771480471313, -0.046931244712973746, 0.4887359962050881),
                           volmdlr.Point3D(-1.2000000000000004, -0.04693124471297373, 0.5477807465830319),
                           volmdlr.Point3D(0.5, 0.1605350614381673, -0.5681325600258332),
                           volmdlr.Point3D(0.41778776086889857, 0.16053506143816723, -0.5496913199951079),
                           volmdlr.Point3D(-0.08774883457241042, 0.16053506143816734, -0.6975390827195442),
                           volmdlr.Point3D(-0.4016938826094424, 0.16053506143816715, 0.1671573143020255),
                           volmdlr.Point3D(-1.046377148047132, 0.16053506143816737, -0.12717731040377675),
                           volmdlr.Point3D(-1.1999999999999997, 0.1605350614381673, -0.06813256002583319),
                           volmdlr.Point3D(0.5, 0.4, -0.5),
                           volmdlr.Point3D(0.4177877608688984, 0.4000000000000003, -0.4815587599692746),
                           volmdlr.Point3D(-0.08774883457241012, 0.39999999999999913, -0.629406522693711),
                           volmdlr.Point3D(-0.4016938826094424, 0.4000000000000003, 0.2352898743278587),
                           volmdlr.Point3D(-1.046377148047132, 0.4000000000000001, -0.059044750377943586),
                           volmdlr.Point3D(-1.2, 0.4, 0.0)]

        for point, expected_point in zip(surface.control_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_contour2d_parametric_to_dimension(self):
        bspline_face = vmf.BSplineFace3D.from_surface_rectangular_cut(bspline_surfaces.bspline_surface_2, 0, 1, 0, 1)
        contour2d = bspline_surfaces.bspline_surface_2.contour3d_to_2d(bspline_face.outer_contour3d)
        grid2d = volmdlr.grid.Grid2D.from_properties((0, 1), (0, 1), (10, 10))
        contour2d_dim = bspline_surfaces.bspline_surface_2.contour2d_parametric_to_dimension(contour2d, grid2d)
        self.assertEqual(len(contour2d_dim.primitives), 4)
        self.assertAlmostEqual(contour2d_dim.area(), 16.657085821451233, places=2)
        self.assertAlmostEqual(contour2d_dim.length(), 16.823814079415172, places=2)

    def test_periodicity(self):
        bspline_suface = surfaces.BSplineSurface3D.from_json(os.path.join(folder, 'surface3d_8.json'))
        self.assertAlmostEqual(bspline_suface.x_periodicity,  0.8888888888888888)
        self.assertFalse(bspline_suface.y_periodicity)

    def test_bbox(self):
        surface = bspline_surfaces.bspline_surface_3
        bbox = surface.bounding_box
        volume = bbox.volume()

        # Check if the bounding box volume is correct
        self.assertAlmostEqual(volume, 3.97787, 2)

    def test_arc3d_to_2d(self):
        bspline_surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, 'BSplineSurface3D_with_Arc3D.json'))
        arc = vme.Arc3D.from_3_points(volmdlr.Point3D(-0.01, -0.013722146986970815, 0.026677756316261864),
                        volmdlr.Point3D(-0.01, 0.013517082603, 0.026782241839),
                        volmdlr.Point3D(-0.01, 0.029612430603, 0.004806657236))

        test = bspline_surface.arc3d_to_2d(arc3d=arc)[0]

        inv_prof = bspline_surface.linesegment2d_to_3d(test)[0]

        # Verifies the inversion operation
        self.assertTrue(inv_prof.start.is_close(arc.start))
        # self.assertTrue(inv_prof.interior.is_close(arc.interior))
        self.assertTrue(inv_prof.end.is_close(arc.end))

        # Strange case from step file
        bspline_surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, 'bsplinesurface_arc3d_to_2d_surface.json'))
        arc = vme.Arc3D.from_json(os.path.join(folder, "bsplinesurface_arc3d_to_2d_arc3d.json"))
        brep = bspline_surface.arc3d_to_2d(arc)[0]
        self.assertTrue(brep.start.is_close(volmdlr.Point2D(1, 0)))

    def test_bsplinecurve3d_to_2d(self):
        bspline_surface = bspline_surfaces.bspline_surface_4
        control_points = [
            volmdlr.Point3D(-0.012138106431296442, 0.11769707710908962, -0.10360094389690414),
            volmdlr.Point3D(-0.012153195391844274, 0.1177764571887428, -0.10360691055433219),
            volmdlr.Point3D(-0.01216612946601426, 0.11785649353385147, -0.10361063821784446),
            volmdlr.Point3D(-0.012176888504086755, 0.11793706145749239, -0.10361212108019317)
        ]
        weights = [1.0, 0.9994807070752826, 0.9994807070752826, 1.0]
        original_bspline = vme.BSplineCurve3D(3, control_points, [4, 4], [0, 1], weights, False)
        bspline_on_parametric_domain = bspline_surface.bsplinecurve3d_to_2d(original_bspline)[0]
        bspline_after_transfomation = bspline_surface.linesegment2d_to_3d(bspline_on_parametric_domain)[0]
        original_length = original_bspline.length()
        length_after_transformation = bspline_after_transfomation.length()
        point = original_bspline.point_at_abscissa(0.5 * original_length)
        point_test = bspline_after_transfomation.point_at_abscissa(0.5 * length_after_transformation)
        self.assertAlmostEqual(original_length, length_after_transformation, places=6)
        self.assertTrue(point.is_close(point_test, 1e-6))

        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinesurface_smallbsplinecurve.json"))
        bsplinecurve3d = vme.BSplineCurve3D.from_json(
            os.path.join(folder, "bsplinesurface_smallbsplinecurve_curve.json"))
        brep_primitive = surface.bsplinecurve3d_to_2d(bsplinecurve3d)[0]
        reversed_prof = surface.linesegment2d_to_3d(brep_primitive)[0]
        self.assertAlmostEqual(brep_primitive.length(), 0.0024101173639275997)
        self.assertAlmostEqual(bsplinecurve3d.length(), reversed_prof.length(), 5)

        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "periodic_surface_smallbsplinecurve3d.json"))
        bsplinecurve3d = vme.BSplineCurve3D.from_json(
            os.path.join(folder, "periodic_surface_smallbsplinecurve3d_curve.json"))
        brep_primitive = surface.bsplinecurve3d_to_2d(bsplinecurve3d)[0]
        reversed_prof = surface.linesegment2d_to_3d(brep_primitive)[0]
        self.assertAlmostEqual(brep_primitive.length(), 0.006419627118992597)
        self.assertTrue(bsplinecurve3d.start.is_close(reversed_prof.start))
        self.assertAlmostEqual(bsplinecurve3d.length(), reversed_prof.length(), 5)

        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinecurve3d_to_2d_vclosed_surface_test.json"))
        bsplinecurve3d = vme.BSplineCurve3D.from_json(
            os.path.join(folder, "bsplinecurve3d_to_2d_vclosed_surface_test_curve.json"))
        brep_primitive = surface.bsplinecurve3d_to_2d(bsplinecurve3d)[0]
        reversed_prof = surface.linesegment2d_to_3d(brep_primitive)[0]
        self.assertTrue(brep_primitive.end.is_close(volmdlr.Point2D(0.0, 1.0)))
        self.assertTrue(bsplinecurve3d.start.is_close(reversed_prof.start))

    def test_bsplinecurve2d_to_3d(self):
        surface = surfaces.BSplineSurface3D.from_json(os.path.join(folder, "bspline_surface_with_arcs.json"))
        contour3d = vmw.Contour3D.from_json(os.path.join(folder, "bspline_contour_with_arcs.json"))

        contour2d = surface.contour3d_to_2d(contour3d)
        bspline_1 = contour2d.primitives[0]
        arc3d = surface.bsplinecurve2d_to_3d(bspline_1)[0]
        self.assertTrue(isinstance(bspline_1, vme.BSplineCurve2D))
        self.assertTrue(isinstance(arc3d, vme.Arc3D))

    def test_arcellipse3d_to_2d(self):
        arcellipse = vme.ArcEllipse3D.from_json(os.path.join(folder, "arcellipse_on_bsplinesurface.json"))
        bsplinesurface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinesurface_with_arcellipse.json"))
        test = bsplinesurface.arcellipse3d_to_2d(arcellipse)[0]
        self.assertTrue(isinstance(test, vme.LineSegment2D))
        self.assertTrue(test.start.is_close(volmdlr.Point2D(0.5, 0.0), 1e-4))
        self.assertTrue(test.end.is_close(volmdlr.Point2D(0.5, 1), 1e-4))

        # todo: Uncomment this block when finish debugging contour2d healing
        # surface = surfaces.BSplineSurface3D.from_json(
        #     "surfaces/objects_bspline_test/bspline_surface_self_intersecting_contour.json")
        # contour3d = vmw.Contour3D.from_json(
        #     "surfaces/objects_bspline_test/bspline_contour_self_intersecting_contour.json")
        # face = surface.face_from_contours3d([contour3d])
        # self.assertTrue(face.surface2d.outer_contour.is_ordered())

    def test_fullarcellipse3d_to_2d(self):
        ellipse = vme.FullArcEllipse3D.from_json(
            os.path.join(folder, "bsplinesurface_with_fullarcellipse_fullarcellipse3d.json"))
        bsplinesurface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinesurface_with_fullarcellipse.json"))
        test = bsplinesurface.fullarcellipse3d_to_2d(ellipse)[0]
        self.assertAlmostEqual(test.length(), 1.0, 2)

    def test_contour3d_to_2d(self):
        surface = surfaces.BSplineSurface3D.from_json(os.path.join(folder, "periodicalsurface.json"))
        contour3d = vmw.Contour3D.from_json(os.path.join(folder, "periodicalsurface_contour.json"))
        contour2d = surface.contour3d_to_2d(contour3d)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 1/6, 5)

        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "contour3d_to_2d_small_primitives_surface.json"))
        contour3d = vmw.Contour3D.from_json(os.path.join(folder, "contour3d_to_2d_small_primitives_contour.json"))
        contour2d = surface.contour3d_to_2d(contour3d)
        self.assertTrue(contour2d.is_ordered(1e-2)) # 1e-2 is an acceptable value, because this is parametric dimension

        surface = surfaces.BSplineSurface3D.from_json(os.path.join(folder, "surface_with_singularity.json"))
        contour3d = vmw.Contour3D.from_json(os.path.join(folder, "surface_with_singularity_contour.json"))
        contour2d = surface.contour3d_to_2d(contour3d)
        self.assertTrue(contour2d.is_ordered())

        surface = surfaces.BSplineSurface3D.from_json(os.path.join(folder, "bsplinesurface_nan_bug.json"))
        contour3d = vmw.Contour3D.from_json(os.path.join(folder, "bsplinesurface_nan_bug_contour.json"))
        contour2d = surface.contour3d_to_2d(contour3d)
        self.assertTrue(contour2d.is_ordered())

        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinesurface_with_singularity_point3d_to_2d.json"))
        contour3d = vmw.Contour3D.from_json(
            os.path.join(folder, "bsplinesurface_with_singularity_point3d_to_2d_contour.json"))
        contour2d = surface.contour3d_to_2d(contour3d)
        self.assertIsNotNone(contour2d)

        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinesurface_with_singularity_linesegment3d_to_2d.json"))
        contour3d = vmw.Contour3D.from_json(
            os.path.join(folder, "bsplinesurface_with_singularity_linesegment3d_to_2d_contour.json"))
        contour2d = surface.contour3d_to_2d(contour3d)
        self.assertTrue(contour2d.is_ordered())



    def test_split_surface_u(self):
        surf1, surf2 = self.spline_surf.split_surface_u(0.33)
        expected_point_surf1 = [volmdlr.Point3D(-25.0, -25.0, -10.0), volmdlr.Point3D(-25.0, -15.0, -5.0),
                                volmdlr.Point3D(-25.0, -5.0, 0.0), volmdlr.Point3D(-25.0, 5.0, 0.0),
                                volmdlr.Point3D(-25.0, 15.0, -5.0), volmdlr.Point3D(-25.0, 25.0, -10.0),
                                volmdlr.Point3D(-15.0, -25.0, -8.0), volmdlr.Point3D(-15.0, -15.0, -4.0),
                                volmdlr.Point3D(-15.0, -5.0, -4.0), volmdlr.Point3D(-15.0, 5.0, -4.0),
                                volmdlr.Point3D(-15.0, 15.0, -4.0), volmdlr.Point3D(-15.0, 25.0, -8.0),
                                volmdlr.Point3D(-10.0, -25.0, -6.5), volmdlr.Point3D(-10.0, -15.0, -3.5),
                                volmdlr.Point3D(-10.0, -5.0, -6.0), volmdlr.Point3D(-10.0, 5.0, -6.0),
                                volmdlr.Point3D(-10.0, 15.0, -3.5), volmdlr.Point3D(-10.0, 25.0, -6.5),
                                volmdlr.Point3D(-5.85, -25.0, -5.42), volmdlr.Point3D(-5.85, -15.0, -3.085),
                                volmdlr.Point3D(-5.85, -5.0, -7.0), volmdlr.Point3D(-5.85, 5.0, -7.0),
                                volmdlr.Point3D(-5.85, 15.0, -3.085), volmdlr.Point3D(-5.85, 25.0, -5.42)]
        expected_point_surf2 = [
                volmdlr.Point3D(-5.85, -25.0, -5.42), volmdlr.Point3D(-5.85, -15.0, -3.085),
                volmdlr.Point3D(-5.85, -5.0, -7.0), volmdlr.Point3D(-5.85, 5.0, -7.0),
                volmdlr.Point3D(-5.85, 15.0, -3.085), volmdlr.Point3D(-5.85, 25.0, -5.42),
                volmdlr.Point3D(-1.6999999999999995, -25.0, -4.34), volmdlr.Point3D(-1.6999999999999995, -15.0, -2.67),
                volmdlr.Point3D(-1.6999999999999995, -5.0, -8.0), volmdlr.Point3D(-1.6999999999999995, 5.0, -8.0),
                volmdlr.Point3D(-1.6999999999999995, 15.0, -2.67), volmdlr.Point3D(-1.6999999999999995, 25.0, -4.34),
                volmdlr.Point3D(5.0, -25.0, -3.0), volmdlr.Point3D(5.0, -15.0, -2.0),
                volmdlr.Point3D(5.0, -5.0, -8.0), volmdlr.Point3D(5.0, 5.0, -8.0),
                volmdlr.Point3D(5.0, 15.0, -2.0), volmdlr.Point3D(5.0, 25.0, -3.0),
                volmdlr.Point3D(15.0, -25.0, -8.0), volmdlr.Point3D(15.0, -15.0, -4.0),
                volmdlr.Point3D(15.0, -5.0, -4.0), volmdlr.Point3D(15.0, 5.0, -4.0),
                volmdlr.Point3D(15.0, 15.0, -4.0), volmdlr.Point3D(15.0, 25.0, -8.0),
                volmdlr.Point3D(25.0, -25.0, -10.0), volmdlr.Point3D(25.0, -15.0, -5.0),
                volmdlr.Point3D(25.0, -5.0, 2.0), volmdlr.Point3D(25.0, 5.0, 2.0),
                volmdlr.Point3D(25.0, 15.0, -5.0), volmdlr.Point3D(25.0, 25.0, -10.0)]
        self.assertEqual(surf1.nb_u, 4)
        self.assertEqual(surf1.nb_v, 6)
        self.assertEqual(surf2.nb_u, 5)
        self.assertEqual(surf2.nb_v, 6)
        for point, expected_point in zip(surf1.control_points, expected_point_surf1):
            self.assertTrue(point.is_close(expected_point))
        for point, expected_point in zip(surf2.control_points, expected_point_surf2):
            self.assertTrue(point.is_close(expected_point))

    def test_split_surface_u_rational(self):
        surf1, surf2 = self.nurbs_surf.split_surface_u(0.5)
        expected_point_surf1 = [volmdlr.Point3D(-25.0, -25.0, -10.0), volmdlr.Point3D(-25.0, -15.0, -5.0),
                                volmdlr.Point3D(-25.0, -5.0, 0.0), volmdlr.Point3D(-25.0, 5.0, 0.0),
                                volmdlr.Point3D(-25.0, 15.0, -5.0), volmdlr.Point3D(-25.0, 25.0, -10.0),
                                volmdlr.Point3D(-15.0, -25.0, -8.0), volmdlr.Point3D(-15.0, -15.0, -4.0),
                                volmdlr.Point3D(-15.0, -5.0, -4.0), volmdlr.Point3D(-15.0, 5.0, -4.0),
                                volmdlr.Point3D(-15.0, 15.0, -4.0), volmdlr.Point3D(-15.0, 25.0, -8.0),
                                volmdlr.Point3D(-7.424242424242425, -25.0, -5.7272727272727275),
                                volmdlr.Point3D(-7.424242424242425, -15.0, -3.242424242424242),
                                volmdlr.Point3D(-7.424242424242425, -5.0, -7.03030303030303),
                                volmdlr.Point3D(-7.424242424242425, 5.0, -7.03030303030303),
                                volmdlr.Point3D(-7.424242424242425, 15.0, -3.242424242424242),
                                volmdlr.Point3D(-7.424242424242425, 25.0, -5.7272727272727275),
                                volmdlr.Point3D(-1.7998163452708908, -25.0, -4.418732782369146),
                                volmdlr.Point3D(-1.7998163452708908, -15.0, -2.6799816345270893),
                                volmdlr.Point3D(-1.7998163452708908, -5.0, -7.764921946740128),
                                volmdlr.Point3D(-1.7998163452708908, 5.0, -7.764921946740128),
                                volmdlr.Point3D(-1.7998163452708908, 15.0, -2.6799816345270893),
                                volmdlr.Point3D(-1.7998163452708908, 25.0, -4.418732782369146),
                                volmdlr.Point3D(0.11256465850708097, -25.0, -4.238138097577792),
                                volmdlr.Point3D(0.11256465850708097, -15.0, -2.588239271203505),
                                volmdlr.Point3D(0.1125646585070812, -5.0, -7.753361779316869),
                                volmdlr.Point3D(0.11256465850708097, 5.0, -7.753361779316869),
                                volmdlr.Point3D(0.11256465850708097, 15.0, -2.588239271203505),
                                volmdlr.Point3D(0.11256465850708097, 25.0, -4.238138097577792)]
        expected_point_surf2 = [volmdlr.Point3D(0.11256465850708097, -25.0, -4.238138097577792),
                                volmdlr.Point3D(0.11256465850708097, -15.0, -2.588239271203505),
                                volmdlr.Point3D(0.1125646585070812, -5.0, -7.753361779316869),
                                volmdlr.Point3D(0.11256465850708097, 5.0, -7.753361779316869),
                                volmdlr.Point3D(0.11256465850708097, 15.0, -2.588239271203505),
                                volmdlr.Point3D(0.11256465850708097, 25.0, -4.238138097577792),
                                volmdlr.Point3D(1.9124526620628202, -24.999999999999996, -4.068166629538872),
                                volmdlr.Point3D(1.9124526620628202, -15.0, -2.501893517487191),
                                volmdlr.Point3D(1.9124526620628204, -5.0, -7.742481621742036),
                                volmdlr.Point3D(1.9124526620628202, 5.0, -7.742481621742035),
                                volmdlr.Point3D(1.9124526620628202, 15.0, -2.501893517487191),
                                volmdlr.Point3D(1.9124526620628202, 24.999999999999996, -4.068166629538872),
                                volmdlr.Point3D(7.537313432835821, -24.999999999999996, -4.26865671641791),
                                volmdlr.Point3D(7.537313432835821, -15.0, -2.5074626865671643),
                                volmdlr.Point3D(7.537313432835821, -5.0, -6.985074626865671),
                                volmdlr.Point3D(7.537313432835821, 5.0, -6.985074626865671),
                                volmdlr.Point3D(7.537313432835821, 15.0, -2.5074626865671643),
                                volmdlr.Point3D(7.537313432835821, 24.999999999999996, -4.26865671641791),
                                volmdlr.Point3D(15.0, -25.0, -8.0), volmdlr.Point3D(15.0, -15.0, -4.0),
                                volmdlr.Point3D(15.0, -5.0, -4.0), volmdlr.Point3D(15.0, 5.0, -4.0),
                                volmdlr.Point3D(15.0, 15.0, -4.0), volmdlr.Point3D(15.0, 25.0, -8.0),
                                volmdlr.Point3D(25.0, -25.0, -10.0), volmdlr.Point3D(25.0, -15.0, -5.0),
                                volmdlr.Point3D(25.0, -5.0, 2.0), volmdlr.Point3D(25.0, 5.0, 2.0),
                                volmdlr.Point3D(25.0, 15.0, -5.0), volmdlr.Point3D(25.0, 25.0, -10.0)]
        self.assertEqual(surf1.nb_u, 5)
        self.assertEqual(surf1.nb_v, 6)
        self.assertEqual(surf2.nb_u, 5)
        self.assertEqual(surf2.nb_v, 6)
        for point, expected_point in zip(surf1.control_points, expected_point_surf1):
            self.assertTrue(point.is_close(expected_point))
        for point, expected_point in zip(surf2.control_points, expected_point_surf2):
            self.assertTrue(point.is_close(expected_point))

        expected_weights = [0.5, 1.0, 0.75, 1.0, 0.25, 1.0, 0.5, 1.0, 0.75, 1.0, 0.25, 1.0, 0.5, 1.0,
                            0.75, 1.0, 0.25, 1.0, 0.5, 1.0, 0.75, 1.0, 0.25, 1.0, 0.5, 1.0, 0.75, 1.0, 0.25, 1.0]

        for weight, expected_weight in zip(surf1.weights, expected_weights):
            self.assertAlmostEqual(weight, expected_weight)
        for weight, expected_weight in zip(surf2.weights, expected_weights):
            self.assertAlmostEqual(weight, expected_weight)

    def test_split_surface_v(self):
        surf1, surf2 = self.spline_surf.split_surface_v(0.66)
        expected_point_surf1 = [
            volmdlr.Point3D(-25.0, -25.0, -10.0), volmdlr.Point3D(-25.0, -15.0, -5.0),
            volmdlr.Point3D(-25.0, -5.0, 0.0), volmdlr.Point3D(-25.0, 1.6000000000000005, 0.0),
            volmdlr.Point3D(-25.0, 5.700556916908, -1.212965025618178), volmdlr.Point3D(-15.0, -25.0, -8.0),
            volmdlr.Point3D(-15.0, -15.0, -4.0), volmdlr.Point3D(-15.0, -5.0, -4.0),
            volmdlr.Point3D(-15.0, 1.6000000000000005, -4.0), volmdlr.Point3D(-15.0, 5.700556916908, -4.0),
            volmdlr.Point3D(-5.0, -25.0, -5.0), volmdlr.Point3D(-5.0, -15.0, -3.0),
            volmdlr.Point3D(-5.0, -5.0, -8.0), volmdlr.Point3D(-5.0, 1.6000000000000005, -8.0),
            volmdlr.Point3D(-5.0, 5.700556916908, -6.787034974381822), volmdlr.Point3D(5.0, -25.0, -3.0),
            volmdlr.Point3D(5.0, -15.0, -2.0), volmdlr.Point3D(5.0, -5.0, -8.0),
            volmdlr.Point3D(5.0, 1.6000000000000005, -8.0), volmdlr.Point3D(5.0, 5.700556916908, -6.544441969258186),
            volmdlr.Point3D(15.0, -25.0, -8.0), volmdlr.Point3D(15.0, -15.0, -4.0),
            volmdlr.Point3D(15.0, -5.0, -4.0), volmdlr.Point3D(15.0, 1.6000000000000005, -4.0),
            volmdlr.Point3D(15.0, 5.700556916908, -4.0), volmdlr.Point3D(25.0, -25.0, -10.0),
            volmdlr.Point3D(25.0, -15.0, -5.0), volmdlr.Point3D(25.0, -5.0, 2.0),
            volmdlr.Point3D(25.0, 1.6000000000000005, 2.0), volmdlr.Point3D(25.0, 5.700556916908, 0.30184896413455065)]

        expected_point_surf2 = [
            volmdlr.Point3D(-25.0, 5.700556916908, -1.212965025618178),
            volmdlr.Point3D(-25.0, 9.92537313432836, -2.4626865671641793),
            volmdlr.Point3D(-25.0, 15.0, -5.0), volmdlr.Point3D(-25.0, 25.0, -10.0),
            volmdlr.Point3D(-15.0, 5.700556916908, -4.0), volmdlr.Point3D(-15.0, 9.92537313432836, -4.0),
            volmdlr.Point3D(-15.0, 15.0, -4.0), volmdlr.Point3D(-15.0, 25.0, -8.0),
            volmdlr.Point3D(-5.0, 5.700556916908, -6.787034974381822),
            volmdlr.Point3D(-5.0, 9.92537313432836, -5.537313432835821),
            volmdlr.Point3D(-5.0, 15.0, -3.0), volmdlr.Point3D(-5.0, 25.0, -5.0),
            volmdlr.Point3D(5.0, 5.700556916908, -6.544441969258186),
            volmdlr.Point3D(5.0, 9.92537313432836, -5.044776119402984),
            volmdlr.Point3D(5.0, 15.0, -2.0), volmdlr.Point3D(5.0, 25.0, -3.0),
            volmdlr.Point3D(15.0, 5.700556916908, -4.0), volmdlr.Point3D(15.0, 9.92537313432836, -4.0),
            volmdlr.Point3D(15.0, 15.0, -4.0), volmdlr.Point3D(15.0, 25.0, -8.0),
            volmdlr.Point3D(25.0, 5.700556916908, 0.30184896413455065),
            volmdlr.Point3D(25.0, 9.92537313432836, -1.4477611940298512),
            volmdlr.Point3D(25.0, 15.0, -5.0), volmdlr.Point3D(25.0, 25.0, -10.0)]

        self.assertEqual(surf1.nb_u, 6)
        self.assertEqual(surf1.nb_v, 5)
        self.assertEqual(surf2.nb_u, 6)
        self.assertEqual(surf2.nb_v, 4)
        for point, expected_point in zip(surf1.control_points, expected_point_surf1):
            self.assertTrue(point.is_close(expected_point))
        for point, expected_point in zip(surf2.control_points, expected_point_surf2):
            self.assertTrue(point.is_close(expected_point))

    def test_split_surface_v_rational(self):
        surf1, surf2 = self.nurbs_surf.split_surface_v(0.33)
        expected_point_surf1 = [volmdlr.Point3D(-25.0, -25.0, -10.0), volmdlr.Point3D(-25.0, -15.0, -5.0),
                                volmdlr.Point3D(-25.0, -10.714285714285714, -2.857142857142857),
                                volmdlr.Point3D(-25.0, -5.995607613469985, -1.4641288433382138),
                                volmdlr.Point3D(-15.0, -25.0, -8.0), volmdlr.Point3D(-15.0, -15.0, -4.0),
                                volmdlr.Point3D(-15.0, -10.714285714285714, -4.0),
                                volmdlr.Point3D(-15.0, -5.995607613469985, -4.0),
                                volmdlr.Point3D(-5.0, -25.0, -5.0), volmdlr.Point3D(-5.0, -15.0, -3.0),
                                volmdlr.Point3D(-5.0, -10.714285714285714, -5.142857142857143),
                                volmdlr.Point3D(-5.0, -5.995607613469985, -6.535871156661786),
                                volmdlr.Point3D(5.0, -25.0, -3.0), volmdlr.Point3D(5.0, -15.0, -2.0),
                                volmdlr.Point3D(5.0, -10.714285714285714, -4.571428571428571),
                                volmdlr.Point3D(5.0, -5.995607613469985, -6.243045387994144),
                                volmdlr.Point3D(15.0, -25.0, -8.0), volmdlr.Point3D(15.0, -15.0, -4.0),
                                volmdlr.Point3D(15.0, -10.714285714285714, -4.0),
                                volmdlr.Point3D(15.0, -5.995607613469985, -4.0),
                                volmdlr.Point3D(25.0, -25.0, -10.0), volmdlr.Point3D(25.0, -15.0, -5.0),
                                volmdlr.Point3D(25.0, -10.714285714285714, -2.0),
                                volmdlr.Point3D(25.0, -5.995607613469985, -0.04978038067349925)]
        expected_point_surf2 = [
    volmdlr.Point3D(-25.0, -5.995607613469985, -1.4641288433382138), volmdlr.Point3D(-25.0, -1.0360360360360354, 0.0),
     volmdlr.Point3D(-25.0, 5.0, 0.0), volmdlr.Point3D(-25.0, 15.0, -5.0), volmdlr.Point3D(-25.0, 25.0, -10.0),
     volmdlr.Point3D(-15.0, -5.995607613469985, -4.0), volmdlr.Point3D(-15.0, -1.0360360360360354, -4.0),
     volmdlr.Point3D(-15.0, 5.0, -4.0), volmdlr.Point3D(-15.0, 15.0, -4.0), volmdlr.Point3D(-15.0, 25.0, -8.0),
     volmdlr.Point3D(-5.0, -5.995607613469985, -6.535871156661786),
     volmdlr.Point3D(-4.999999999999999, -1.0360360360360354, -8.0), volmdlr.Point3D(-5.0, 5.0, -8.0),
     volmdlr.Point3D(-5.0, 15.0, -3.0), volmdlr.Point3D(-5.0, 25.0, -5.0),
     volmdlr.Point3D(5.0, -5.995607613469985, -6.243045387994144),
     volmdlr.Point3D(4.999999999999999, -1.0360360360360354, -8.0),
     volmdlr.Point3D(5.0, 5.0, -8.0), volmdlr.Point3D(5.0, 15.0, -2.0), volmdlr.Point3D(5.0, 25.0, -3.0),
     volmdlr.Point3D(15.0, -5.995607613469985, -4.0), volmdlr.Point3D(15.0, -1.0360360360360354, -4.0),
     volmdlr.Point3D(15.0, 5.0, -4.0), volmdlr.Point3D(15.0, 15.0, -4.0), volmdlr.Point3D(15.0, 25.0, -8.0),
     volmdlr.Point3D(25.0, -5.995607613469985, -0.04978038067349925), volmdlr.Point3D(25.0, -1.0360360360360354, 2.0),
     volmdlr.Point3D(25.0, 5.0, 2.0), volmdlr.Point3D(25.0, 15.0, -5.0), volmdlr.Point3D(25.0, 25.0, -10.0)
        ]
        self.assertEqual(surf1.nb_u, 6)
        self.assertEqual(surf1.nb_v, 4)
        self.assertEqual(surf2.nb_u, 6)
        self.assertEqual(surf2.nb_v, 5)
        for point, expected_point in zip(surf1.control_points, expected_point_surf1):
            self.assertTrue(point.is_close(expected_point))
        for point, expected_point in zip(surf2.control_points, expected_point_surf2):
            self.assertTrue(point.is_close(expected_point))

        expected_weights_surf1 = [0.5, 1.0, 0.875, 0.85375, 0.5, 1.0, 0.875, 0.85375, 0.5, 1.0, 0.875, 0.85375,
                                  0.5, 1.0, 0.875, 0.85375, 0.5, 1.0, 0.875, 0.85375, 0.5, 1.0, 0.875, 0.85375]

        expected_weights_surf2 = [0.85375, 0.8325, 1.0, 0.25, 1.0, 0.85375, 0.8325, 1.0, 0.25, 1.0,
                                  0.85375, 0.8325, 1.0, 0.25, 1.0, 0.85375, 0.8325, 1.0, 0.25, 1.0,
                                  0.85375, 0.8325, 1.0, 0.25, 1.0, 0.85375, 0.8325, 1.0, 0.25, 1.0]

        for weight, expected_weight in zip(surf1.weights, expected_weights_surf1):
            self.assertAlmostEqual(weight, expected_weight)
        for weight, expected_weight in zip(surf2.weights, expected_weights_surf2):
            self.assertAlmostEqual(weight, expected_weight)

    def test_surface_curves(self):
        curves = self.spline_surf.surface_curves
        u_curves = curves["u"]
        v_curves = curves["v"]
        self.assertEqual(len(u_curves), 6)
        self.assertEqual(len(v_curves), 6)
        self.assertTrue(u_curves[0].start.is_close(volmdlr.Point3D(-25.0, -25.0, -10.0)))
        self.assertTrue(u_curves[0].end.is_close(volmdlr.Point3D(25.0, -25.0, -10.0)))
        self.assertTrue(u_curves[-1].start.is_close(volmdlr.Point3D(-25.0, 25.0, -10.0)))
        self.assertTrue(u_curves[-1].end.is_close(volmdlr.Point3D(25.0, 25.0, -10.0)))
        self.assertTrue(v_curves[0].start.is_close(volmdlr.Point3D(-25.0, -25.0, -10.0)))
        self.assertTrue(v_curves[0].end.is_close(volmdlr.Point3D(-25.0, 25.0, -10.0)))
        self.assertTrue(v_curves[-1].start.is_close(volmdlr.Point3D(25.0, -25.0, -10.0)))
        self.assertTrue(v_curves[-1].end.is_close(volmdlr.Point3D(25.0, 25.0, -10.0)))

        curves = self.nurbs_surf.surface_curves
        u_curves = curves["u"]
        v_curves = curves["v"]
        self.assertEqual(len(u_curves), 6)
        self.assertEqual(len(v_curves), 6)
        self.assertTrue(u_curves[0].start.is_close(volmdlr.Point3D(-25.0, -25.0, -10.0)))
        self.assertTrue(u_curves[0].end.is_close(volmdlr.Point3D(25.0, -25.0, -10.0)))
        self.assertTrue(u_curves[-1].start.is_close(volmdlr.Point3D(-25.0, 25.0, -10.0)))
        self.assertTrue(u_curves[-1].end.is_close(volmdlr.Point3D(25.0, 25.0, -10.0)))
        self.assertTrue(v_curves[0].start.is_close(volmdlr.Point3D(-25.0, -25.0, -10.0)))
        self.assertTrue(v_curves[0].end.is_close(volmdlr.Point3D(-25.0, 25.0, -10.0)))
        self.assertTrue(v_curves[-1].start.is_close(volmdlr.Point3D(25.0, -25.0, -10.0)))
        self.assertTrue(v_curves[-1].end.is_close(volmdlr.Point3D(25.0, 25.0, -10.0)))

        u0_expected_weights = [0.5, 0.5, 0.5, 0.5, 0.5, 0.5]
        v0_expected_weights = [0.5, 1.0, 0.75, 1.0, 0.25, 1.0]
        for curve, expected_weights in zip((u_curves[0], v_curves[0]), (u0_expected_weights, v0_expected_weights)):
            for weight, expected_weight in zip(curve.weights, expected_weights):
                self.assertAlmostEqual(weight, expected_weight)

    def test_plane_intersections(self):
        frame = volmdlr.Frame3D(volmdlr.O3D, volmdlr.Z3D, volmdlr.X3D, volmdlr.Y3D)
        plane = surfaces.Plane3D(frame)
        intersections = self.spline_surf.plane_intersections(plane)
        for point in intersections:
            self.assertTrue(plane.point_belongs(point))

    def test_decompose(self):
        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplineface_triangulation_problem_surface.json"))
        decompose_results = surface.decompose(return_params=True)
        self.assertEqual(len(decompose_results), 116)
        bezier_patches = [result[0] for result in decompose_results]
        params = [result[1] for result in decompose_results]
        self.assertEqual(len(bezier_patches), len(params))
        for patch, param in zip(bezier_patches, params):
            control_points = patch.control_points
            self.assertTrue(
                surface.point2d_to_3d(volmdlr.Point2D(param[0][0], param[1][0])).is_close(control_points[0]))
            self.assertTrue(
                surface.point2d_to_3d(volmdlr.Point2D(param[0][1], param[1][1])).is_close(control_points[-1]))
        decompose_results = surface.decompose(return_params=True, decompose_dir="u")
        bezier_patches = [result[0] for result in decompose_results]
        params = [result[1] for result in decompose_results]
        self.assertEqual(len(bezier_patches), 4)
        self.assertEqual(len(bezier_patches), len(params))
        for patch, param in zip(bezier_patches, params):
            control_points = patch.control_points
            self.assertTrue(
                surface.point2d_to_3d(volmdlr.Point2D(param[0][0], param[1][0])).is_close(control_points[0]))
            self.assertTrue(
                surface.point2d_to_3d(volmdlr.Point2D(param[0][1], param[1][1])).is_close(control_points[-1]))
        decompose_results = surface.decompose(return_params=True, decompose_dir="v")
        bezier_patches = [result[0] for result in decompose_results]
        params = [result[1] for result in decompose_results]
        self.assertEqual(len(bezier_patches), 29)
        self.assertEqual(len(bezier_patches), len(params))
        for patch, param in zip(bezier_patches, params):
            control_points = patch.control_points
            self.assertTrue(
                surface.point2d_to_3d(volmdlr.Point2D(param[0][0], param[1][0])).is_close(control_points[0]))
            self.assertTrue(
                surface.point2d_to_3d(volmdlr.Point2D(param[0][1], param[1][1])).is_close(control_points[-1]))

        bezier_patches = surface.decompose(decompose_dir="u")
        self.assertEqual(len(bezier_patches), 4)
        bezier_patches = surface.decompose(decompose_dir="v")
        self.assertEqual(len(bezier_patches), 29)
        bezier_patches = surface.decompose()
        self.assertEqual(len(bezier_patches), 116)

    def test_point_inversion_grid_search(self):
        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinesurface_point3d_to_2d_grid_search_1.json"))
        point = volmdlr.Point3D(-0.009668298046654873, 0.11887869426572631, -0.09560417062522625)
        _, distance = surface.point_inversion_grid_search(point, 5e-5, 2)
        self.assertLess(distance, 3e-5)
        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinesurface_point3d_to_2d_grid_search_2.json"))
        point = volmdlr.Point3D(0.001702815989525993, 0.003297577223278291, -0.026314554505063058)
        _, distance = surface.point_inversion_grid_search(point, 5e-5, 2)
        self.assertLess(distance, 2e-5)

        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinesurface_point3d_to_2d_grid_search_3.json"))
        point = volmdlr.Point3D(-0.008941313467488011, 0.01194521078356664, -0.000635664858372182)
        _, distance = surface.point_inversion_grid_search(point, 5e-5, 2)
        self.assertLess(distance, 5e-5)

        surface = surfaces.BSplineSurface3D.from_json(
            os.path.join(folder, "bsplinesurface_point3d_to_2d_grid_search_4.json"))
        point = volmdlr.Point3D(-0.02494082957642294, 0.03166087761892587, -0.07334489785111517)
        _, distance = surface.point_inversion_grid_search(point, 5e-5, 2)
        self.assertLess(distance, 5e-5)


if __name__ == '__main__':
    unittest.main(verbosity=0)
