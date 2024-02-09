"""
Unit tests for CylindriSurface3D
"""
import unittest
import math
import numpy as np
import os
import dessia_common.core
import volmdlr
from volmdlr import Point2D, Point3D, edges, wires, surfaces, curves
from volmdlr.models import cylindrical_surfaces


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_cylindrical_tests')


class TestCylindricalSurface3D(unittest.TestCase):
    cylindrical_surface = cylindrical_surfaces.cylindrical_surface
    cylindrical_surface2 = cylindrical_surfaces.cylindrical_surface2
    cylindrical_surface3 = cylindrical_surfaces.cylindrical_surface3
    cylindrical_surface4 = cylindrical_surfaces.cylindrical_surface4

    def test_line_intersections(self):
        line3d = curves.Line3D(volmdlr.O3D, volmdlr.Point3D(0.3, 0.3, 0.3))
        line_inters = self.cylindrical_surface.line_intersections(line3d)
        self.assertEqual(len(line_inters), 2)
        self.assertTrue(
            line_inters[0].is_close(volmdlr.Point3D(0.22627416, 0.22627416, 0.22627416))
        )
        self.assertTrue(
            line_inters[1].is_close(
                volmdlr.Point3D(-0.22627416, -0.22627416, -0.22627416)
            )
        )

    def test_plane_intersections(self):
        plane_surface = surfaces.Plane3D(volmdlr.OZXY)
        parallel_plane_secant_cylinder = plane_surface.frame_mapping(
            volmdlr.Frame3D(
                volmdlr.O3D,
                volmdlr.Point3D(0.5, 0, 0),
                volmdlr.Point3D(0, 0.5, 0),
                volmdlr.Point3D(0, 0, 0.2),
            ),
            "new",
        )
        cylinder_tanget_plane = plane_surface.frame_mapping(
            volmdlr.Frame3D(
                volmdlr.Point3D(0, 0.32 / 2, 0),
                volmdlr.Point3D(0.5, 0, 0),
                volmdlr.Point3D(0, 0.5, 0),
                volmdlr.Point3D(0, 0, 0.2),
            ),
            "new",
        )
        not_intersecting_cylinder_parallel_plane = plane_surface.frame_mapping(
            volmdlr.Frame3D(
                volmdlr.Point3D(0, 0.32, 0),
                volmdlr.Point3D(0.5, 0, 0),
                volmdlr.Point3D(0, 0.5, 0),
                volmdlr.Point3D(0, 0, 0.2),
            ),
            "new",
        )
        cylinder_concurrent_plane = plane_surface.rotation(
            volmdlr.O3D, volmdlr.X3D, math.pi / 4
        )
        cylinder_perpendicular_plane = plane_surface.rotation(
            volmdlr.O3D, volmdlr.X3D, math.pi / 2
        )

        cylinder_surface_secant_parallel_plane_intersec = (
            self.cylindrical_surface.plane_intersections(parallel_plane_secant_cylinder)
        )
        self.assertEqual(len(cylinder_surface_secant_parallel_plane_intersec), 2)
        self.assertTrue(
            isinstance(
                cylinder_surface_secant_parallel_plane_intersec[0], curves.Line3D
            )
        )
        self.assertTrue(
            isinstance(
                cylinder_surface_secant_parallel_plane_intersec[1], curves.Line3D
            )
        )

        cylinder_surface_tangent_plane = self.cylindrical_surface.plane_intersections(
            cylinder_tanget_plane
        )
        self.assertEqual(len(cylinder_surface_tangent_plane), 1)
        self.assertTrue(isinstance(cylinder_surface_tangent_plane[0], curves.Line3D))

        cylinder_surface_tangent_plane_not_intersecting = (
            self.cylindrical_surface.plane_intersections(
                not_intersecting_cylinder_parallel_plane
            )
        )
        self.assertEqual(len(cylinder_surface_tangent_plane_not_intersecting), 0)

        cylinder_surface_concurrent_plane_intersec = (
            self.cylindrical_surface.plane_intersections(cylinder_concurrent_plane)
        )
        self.assertTrue(
            isinstance(cylinder_surface_concurrent_plane_intersec[0], curves.Ellipse3D)
        )

        cylinder_surface_perpendicular_plane_intersec = (
            self.cylindrical_surface.plane_intersections(cylinder_perpendicular_plane)
        )
        self.assertTrue(
            isinstance(
                cylinder_surface_perpendicular_plane_intersec[0], curves.Circle3D
            )
        )

    def test_is_coincident(self):
        cyl_surface1 = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 1)
        cyl_surface2 = surfaces.CylindricalSurface3D(
            volmdlr.OXYZ.translation(volmdlr.Vector3D(0, 0, 1)), 1
        )
        plane_face = surfaces.Plane3D(volmdlr.OXYZ)
        self.assertTrue(cyl_surface1.is_coincident(cyl_surface2))
        self.assertFalse(cyl_surface1.is_coincident(plane_face))
        self.assertFalse(self.cylindrical_surface.is_coincident(cyl_surface1))

    def test_point_belongs(self):
        point = volmdlr.Point3D(0.32, 0, 1)
        point2 = volmdlr.Point3D(1, 1, 1)
        self.assertTrue(self.cylindrical_surface.point_belongs(point))
        self.assertFalse((self.cylindrical_surface.point_belongs(point2)))

    def test_linesegment2d_to_3d(self):
        point1 = volmdlr.Point2D(0.0, 0.0)
        u_linesegment = edges.LineSegment2D(point1, volmdlr.Point2D(0.0, 1.0))
        v_linesegment = edges.LineSegment2D(point1, volmdlr.Point2D(math.pi, 0.0))
        uv_linesegment = edges.LineSegment2D(point1, volmdlr.Point2D(math.pi, 1.0))

        linesegment3d = self.cylindrical_surface2.linesegment2d_to_3d(u_linesegment)[0]
        self.assertAlmostEqual(linesegment3d.length(), 1.0)
        self.assertTrue(linesegment3d.start.is_close(volmdlr.Point3D(1.0, 0.0, 0.0)))
        self.assertTrue(linesegment3d.end.is_close(volmdlr.Point3D(1.0, 0.0, 1.0)))

        arc3d = self.cylindrical_surface2.linesegment2d_to_3d(v_linesegment)[0]
        self.assertAlmostEqual(arc3d.radius, 1.0)
        self.assertTrue(arc3d.start.is_close(volmdlr.Point3D(1.0, 0.0, 0.0)))
        self.assertTrue(arc3d.end.is_close(volmdlr.Point3D(-1.0, 0.0, 0.0)))

        arc3d = self.cylindrical_surface2.linesegment2d_to_3d(v_linesegment.reverse())[0]
        self.assertAlmostEqual(arc3d.radius, 1.0)
        self.assertTrue(arc3d.start.is_close(volmdlr.Point3D(-1.0, 0.0, 0.0)))
        self.assertTrue(arc3d.end.is_close(volmdlr.Point3D(1.0, 0.0, 0.0)))
        self.assertAlmostEqual(self.cylindrical_surface2.frame.w.dot(arc3d.circle.normal), -1.0)

        bsplinecurve3d = self.cylindrical_surface2.linesegment2d_to_3d(uv_linesegment)[0]
        distances = [self.cylindrical_surface2.point_distance(point)
                     for point in bsplinecurve3d.discretization_points(number_points=100)]
        self.assertLess(max(distances), 1e-7)

    def test_parametric_points_to_3d(self):
        parametric_points = np.array([[0.0, 0.0], [0.5 * math.pi, 0.0], [math.pi, 0.0], [1.5 * math.pi, 0.0],
                                      [0.0, 1.0], [0.5 * math.pi, 1.0], [math.pi, 1.0], [1.5 * math.pi, 1.0]])
        points3d = self.cylindrical_surface.parametric_points_to_3d(parametric_points)
        expected_points = np.array([[0.32, 0.0, 0.0], [0.0, 0.32, 0.0], [-0.32, 0.0, 0.0], [0.0, -0.32, 0.0],
                                    [0.32, 0.0, 1.0], [0.0, 0.32, 1.0], [-0.32, 0.0, 1.0], [0.0, -0.32, 1.0]])
        for point, expected_point in zip(points3d, expected_points):
            self.assertAlmostEqual(np.linalg.norm(point - expected_point), 0.0)

    def test_arcellipse3d_to_2d(self):
        pass

    def test_linesegment3d_to_2d(self):
        surface = surfaces.CylindricalSurface3D.from_json(
            os.path.join(folder, "cylindricalsurface_with_linesegment3d.json"))
        linesegment3d = edges.LineSegment3D.from_json(
            os.path.join(folder, "cylindricalsurface_linesegment3d.json"))
        linesegment2d = surface.linesegment3d_to_2d(linesegment3d)[0]
        self.assertTrue(
            linesegment2d.start.is_close(
                volmdlr.Point2D(-0.021051754138835845, -0.0033749825505284136)
            )
        )
        self.assertTrue(
            linesegment2d.end.is_close(volmdlr.Point2D(0.0, -0.0033725697172752008))
        )

    def test_arc3d_to_2d(self):
        arc1 = edges.Arc3D.from_3_points(
            volmdlr.Point3D(1, 0, 0),
            volmdlr.Point3D(1 / math.sqrt(2), 1 / math.sqrt(2), 0),
            volmdlr.Point3D(0, 1, 0),
        )
        arc2 = edges.Arc3D.from_3_points(
            volmdlr.Point3D(1, 0, 0),
            volmdlr.Point3D(1 / math.sqrt(2), -1 / math.sqrt(2), 0),
            volmdlr.Point3D(0, -1, 0),
        )
        arc3 = edges.Arc3D.from_3_points(
            volmdlr.Point3D(-1 / math.sqrt(2), 1 / math.sqrt(2), 0),
            volmdlr.Point3D(-1, 0, 0),
            volmdlr.Point3D(-1 / math.sqrt(2), -1 / math.sqrt(2), 0),
        )
        arc4 = edges.Arc3D.from_3_points(
            volmdlr.Point3D(0, -1, 0),
            volmdlr.Point3D(-1 / math.sqrt(2), 1 / math.sqrt(2), 0),
            volmdlr.Point3D(1, 0, 0),
        )
        test1 = self.cylindrical_surface2.arc3d_to_2d(arc3d=arc1)[0]
        test2 = self.cylindrical_surface2.arc3d_to_2d(arc3d=arc2)[0]
        test3 = self.cylindrical_surface2.arc3d_to_2d(arc3d=arc3)[0]
        test4 = self.cylindrical_surface2.arc3d_to_2d(arc3d=arc4)[0]

        inv_prof = self.cylindrical_surface2.linesegment2d_to_3d(test4)[0]

        # Assert that the returned object is an edges.LineSegment2D
        self.assertIsInstance(test1, edges.LineSegment2D)
        self.assertIsInstance(test2, edges.LineSegment2D)
        self.assertIsInstance(test3, edges.LineSegment2D)
        self.assertIsInstance(test4, edges.LineSegment2D)

        # Assert that the returned object is right on the parametric domain (take into account periodicity)
        self.assertTrue(test1.start.is_close(volmdlr.Point2D(0, 0)))
        self.assertTrue(test1.end.is_close(volmdlr.Point2D(0.5 * math.pi, 0)))
        self.assertTrue(test2.start.is_close(volmdlr.Point2D(0, 0)))
        self.assertTrue(test2.end.is_close(volmdlr.Point2D(-0.5 * math.pi, 0)))
        self.assertTrue(test3.start.is_close(volmdlr.Point2D(0.75 * math.pi, 0)))
        self.assertTrue(test3.end.is_close(volmdlr.Point2D(1.25 * math.pi, 0)))
        self.assertTrue(test4.start.is_close(volmdlr.Point2D(-0.5 * math.pi, 0)))
        self.assertTrue(test4.end.is_close(volmdlr.Point2D(-2 * math.pi, 0)))

        # Verifies the inversion operation
        self.assertIsInstance(inv_prof, edges.Arc3D)
        # self.assertEqual(inv_prof, arc4)
        self.assertTrue(inv_prof.start.is_close(arc4.start))
        self.assertTrue(inv_prof.middle_point().is_close(arc4.middle_point()))
        self.assertTrue(inv_prof.end.is_close(arc4.end))

    def test_contour3d_to_2d(self):
        center1 = Point3D(0, 0, 0.013)
        start_end1 = Point3D(0.03, 0, 0.013)
        circle1 = curves.Circle3D(
            volmdlr.Frame3D(center1, volmdlr.X3D, -volmdlr.Y3D, -volmdlr.Z3D),
            center1.point_distance(start_end1),
        )
        center2 = Point3D(0, 0, 0.003)
        start_end2 = Point3D(0.03, 0, 0.003)
        circle2 = curves.Circle3D(
            volmdlr.Frame3D(center2, volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D),
            center2.point_distance(start_end2),
        )
        primitives_cylinder = [
            edges.LineSegment3D(Point3D(0.03, 0, 0.003), Point3D(0.03, 0, 0.013)),
            edges.FullArc3D(circle1, start_end1),
            edges.LineSegment3D(Point3D(0.03, 0, 0.013), Point3D(0.03, 0, 0.003)),
            edges.FullArc3D(circle2, start_end2),
        ]
        contour_cylinder = wires.Contour3D(primitives_cylinder)

        contour2d_cylinder = self.cylindrical_surface4.contour3d_to_2d(contour_cylinder)
        # ax = contour2d_cylinder.plot()
        # ax.set_aspect("auto")
        area = contour2d_cylinder.area()
        fullarc2d = contour2d_cylinder.primitives[3]
        linesegment2d = contour2d_cylinder.primitives[0]

        self.assertEqual(area, 0.02 * math.pi)
        self.assertEqual(fullarc2d.start, Point2D(-volmdlr.TWO_PI, 0.003))
        self.assertEqual(fullarc2d.end, Point2D(0, 0.003))
        self.assertEqual(linesegment2d.start, Point2D(0, 0.003))
        self.assertEqual(linesegment2d.end, Point2D(0, 0.013))

        surface = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, "cylindrical_surface_bspline_openned_contour.json"))
        contour = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder,"cylindrical_contour_bspline_openned_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 2)
        self.assertFalse(contour2d.is_ordered())

        surface = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, "test_contour3d_to_2d_surface.json"
        ))
        contour = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, "test_contour3d_to_2d_contour.json"
        ))

        surface = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, "test_contour3d_to_2d_surface.json"))
        contour = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, "test_contour3d_to_2d_contour.json"))

        contour2d = surface.contour3d_to_2d(contour)
        self.assertAlmostEqual(contour2d.area(), 0.29361767646954695, 2)
        self.assertTrue(contour2d.is_ordered())

        surface = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, "cylindricalsurface_small_periodic_bsplinecurve.json"))
        contour = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, "cylindricalsurface_small_periodic_bsplinecurve_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertAlmostEqual(contour2d.area(), 0.0, 6)
        self.assertTrue(contour2d.is_ordered())


    def test_bsplinecurve3d_to_2d(self):
        surface = dessia_common.core.DessiaObject.from_json(os.path.join(folder, "cylindrical_surf_bug.json"))
        bsplinecurve3d = dessia_common.core.DessiaObject.from_json(os.path.join(folder, "bsplinecurve3d_bug.json"))
        primitive2d = surface.bsplinecurve3d_to_2d(bsplinecurve3d)[0]
        self.assertTrue(
            primitive2d.start.is_close(
                volmdlr.Point2D(-0.001540582016168617, -0.0006229082591074433)
            )
        )
        self.assertTrue(
            primitive2d.end.is_close(
                volmdlr.Point2D(0.004940216577284154, -0.000847814405768888)
            )
        )

        # Test to _fix_angle_discontinuity_on_discretization_points
        z = np.linspace(0, 2 * math.pi, 50)
        theta = math.pi + 0.5 * math.pi * np.cos(z)
        points_2d = [volmdlr.Point2D(x, y / (2 * math.pi)) for x, y in zip(theta, z)]
        cylinder = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 1)
        points_3d = [cylinder.point2d_to_3d(point) for point in points_2d]
        bspline = edges.BSplineCurve3D.from_points_interpolation(points_3d, 3)
        result = cylinder.bsplinecurve3d_to_2d(bspline)[0]
        self.assertTrue(all(point.x < 0 for point in result.points))

    def test_point_projection(self):
        test_points = [
            Point3D(-2.0, -2.0, 0.0),
            Point3D(0.0, -2.0, 0.0),
            Point3D(2.0, -2.0, 0.0),
            Point3D(2.0, 0.0, 0.0),
            Point3D(2.0, 2.0, 0.0),
            Point3D(0.0, 2.0, 0.0),
            Point3D(-2.0, 2.0, 0.0),
            Point3D(-2.0, 0.0, 0.0),
        ]
        expected_points = [
            volmdlr.Point3D(-0.5 * math.sqrt(2), -0.5 * math.sqrt(2), 0.0),
            volmdlr.Point3D(0.0, -1.0, 0.0),
            volmdlr.Point3D(0.5 * math.sqrt(2), -0.5 * math.sqrt(2), 0.0),
            volmdlr.Point3D(1.0, 0.0, 0.0),
            volmdlr.Point3D(0.5 * math.sqrt(2), 0.5 * math.sqrt(2), 0.0),
            volmdlr.Point3D(0.0, 1.0, 0.0),
            volmdlr.Point3D(-0.5 * math.sqrt(2), 0.5 * math.sqrt(2), 0.0),
            volmdlr.Point3D(-1.0, 0.0, 0.0),
        ]

        for i, point in enumerate(test_points):
            self.assertTrue(
                self.cylindrical_surface2.point_projection(point).is_close(
                    expected_points[i]
                )
            )

    def test_plot(self):
        ax = self.cylindrical_surface.plot()
        self.assertTrue(ax)

    def test_ellipse_intersections(self):
        cyl_surface = surfaces.CylindricalSurface3D(
            volmdlr.Frame3D(
                origin=volmdlr.Point3D(1.0, 1.0, 0.0),
                u=volmdlr.Vector3D(-5.551115123125783e-17, 0.0, 0.9999999999999998),
                v=volmdlr.Vector3D(0.0, 0.9999999999999998, 0.0),
                w=volmdlr.Vector3D(-0.9999999999999998, 0.0, -5.551115123125783e-17),
            ),
            1.5,
        )

        frame = volmdlr.Frame3D(
            origin=volmdlr.Point3D(0.0, 0.0, 0.0),
            u=volmdlr.Vector3D(
                0.5773502691896258, 0.5773502691896258, 0.5773502691896258
            ),
            v=volmdlr.Vector3D(
                0.8164965809277258, -0.40824829046386313, -0.40824829046386313
            ),
            w=volmdlr.Vector3D(0.0, 0.7071067811865476, -0.7071067811865476),
        )
        ellipse = curves.Ellipse3D(2, 1, frame)
        ellipse_intersections = cyl_surface.ellipse_intersections(ellipse)
        self.assertEqual(len(ellipse_intersections), 2)
        self.assertTrue(
            ellipse_intersections[0],
            volmdlr.Point3D(
                -1.3695411442140175, -0.4354143466934857, -0.4354143466934857
            ),
        )
        self.assertTrue(
            ellipse_intersections[1],
            volmdlr.Point3D(
                0.7889886819560367, -0.4354143466934856, -0.4354143466934856
            ),
        )

    def test_coinicalsurface_intersections(self):
        expected_slutions = [
            [3.710032833168665],
            [2.754671034122705, 0.7935213452250598],
            [2.07512669883945, 2.075126698839449],
            [2.569944707187624, 2.569944707187624],
        ]
        conical_surface = surfaces.ConicalSurface3D(volmdlr.OXYZ, math.pi / 6)
        cylindrical_surface1 = surfaces.CylindricalSurface3D(
            volmdlr.Frame3D(
                volmdlr.Point3D(0.3, 0.3, 0.8), volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D
            ),
            0.3,
        )
        cylindrical_surface2 = surfaces.CylindricalSurface3D(
            volmdlr.Frame3D(
                volmdlr.Point3D(0, 0, 0.5), volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D
            ),
            0.3,
        )
        cylindrical_surface3 = surfaces.CylindricalSurface3D(
            volmdlr.Frame3D(
                volmdlr.Point3D(0, 0, 1), volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D
            ),
            0.3,
        )
        cylindrical_surface4 = surfaces.CylindricalSurface3D(
            volmdlr.Frame3D(
                volmdlr.Point3D(0.0, 0.41068360252295905, 1.2886751345948129),
                volmdlr.Y3D,
                volmdlr.Z3D,
                volmdlr.X3D,
            ),
            math.tan(conical_surface.semi_angle) / 2,
        )
        for i, cylindrical_surface in enumerate(
            [
                cylindrical_surface1,
                cylindrical_surface2,
                cylindrical_surface3,
                cylindrical_surface4,
            ]
        ):
            list_curves = cylindrical_surface.conicalsurface_intersections(
                conical_surface
            )
            for curve, expected_length in zip(list_curves, expected_slutions[i]):
                self.assertAlmostEqual(curve.length(), expected_length)

    def test_circle_intersections(self):
        cylindrical_surface = surfaces.CylindricalSurface3D(
            volmdlr.Frame3D(
                origin=volmdlr.Point3D(1.0, 1.0, 0.0),
                u=volmdlr.Vector3D(0.7071067811865476, 0.0, -0.7071067811865475),
                v=volmdlr.Vector3D(0.0, 1.0, 0.0),
                w=volmdlr.Vector3D(0.7071067811865475, 0.0, 0.7071067811865476),
            ),
            1.5,
        )
        circle = curves.Circle3D(
            volmdlr.Frame3D(
                origin=volmdlr.Point3D(1.9960534568565431, 0.12558103905862675, 0.0),
                u=volmdlr.Vector3D(0.9980267284282716, 0.06279051952931337, 0.0),
                v=volmdlr.Vector3D(0.0, 0.0, 1.0),
                w=volmdlr.Vector3D(0.06279051952931337, -0.9980267284282716, 0.0),
            ),
            1,
        )
        circle_intersections = cylindrical_surface.circle_intersections(circle)
        self.assertTrue(
            circle_intersections[0].is_close(
                volmdlr.Point3D(2.975410303031, 0.187196949158, 0.19251871881)
            )
        )
        self.assertTrue(
            circle_intersections[1].is_close(
                volmdlr.Point3D(1.740409914248, 0.10949731064, -0.966637220568)
            )
        )

    def test_sphericalsurface_intersections(self):
        spherical_surface = surfaces.SphericalSurface3D(
            volmdlr.OXYZ.translation(volmdlr.Vector3D(0.5, 0.5, 0)), 2
        )

        # test 1
        cylindrical_surface = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 1)
        inters = spherical_surface.surface_intersections(cylindrical_surface)
        self.assertEqual(len(inters), 2)
        self.assertAlmostEqual(inters[0].length(), 6.613411150146185)
        self.assertAlmostEqual(inters[1].length(), 6.613411150146188)
      
        # test2
        cylindrical_surface = surfaces.CylindricalSurface3D(
            volmdlr.OXYZ.translation(volmdlr.X3D * 1.5), 1
        )
        inters = spherical_surface.surface_intersections(cylindrical_surface)
        self.assertEqual(len(inters), 1)
        self.assertAlmostEqual(inters[0].length(), 13.535573972180183)

        # test 3 - sphere inside the cylinder, sphere radius < cylinder radius.
        spherical_surface = surfaces.SphericalSurface3D(
            volmdlr.OXYZ.translation(volmdlr.Vector3D(0.5, 0.5, 0)), 0.5
        )
        cylindrical_surface = surfaces.CylindricalSurface3D(
            volmdlr.OXYZ.translation(volmdlr.X3D * 1.5), 2
        )
        inters3 = spherical_surface.surface_intersections(cylindrical_surface)
        self.assertFalse(inters3)

        # test 4 - sphere outside the cylinder
        spherical_surface = surfaces.SphericalSurface3D(
            volmdlr.OXYZ.translation(volmdlr.Vector3D(0.5, 3, 0)), 0.5
        )
        self.assertFalse(spherical_surface.surface_intersections(cylindrical_surface))

        # test 5 - Sphere on the cylinder axis, sphere radius  == cylinder radius -> One circle 3D at z = 0.
        spherical_surface = surfaces.SphericalSurface3D(volmdlr.OXYZ, 2)
        # spherical_surface = surfaces.SphericalSurface3D(volmdlr.OXYZ, 3)
        cylindrical_surface = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 2)

        inters5 = spherical_surface.surface_intersections(cylindrical_surface)
        self.assertTrue(len(inters5), 1)
        self.assertTrue(inters5[0].center.is_close(volmdlr.O3D))
        self.assertEqual(inters5[0].radius, spherical_surface.radius)

        # test 6 - Sphere on cylinder axis, sphere R > cylinder r -> Two circle 3D at z = -+ sqrt(R**2 - r**2).
        spherical_surface = surfaces.SphericalSurface3D(volmdlr.OXYZ, 3)
        cylindrical_surface = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 2)

        inters6 = spherical_surface.surface_intersections(cylindrical_surface)
        self.assertTrue(len(inters5), 2)
        self.assertTrue(
            inters6[0].center.is_close(volmdlr.Point3D(0.0, 0.0, 2.23606797749979))
        )
        self.assertAlmostEqual(inters6[0].radius, 1.9999999999999998)
        self.assertTrue(
            inters6[1].center.is_close(volmdlr.Point3D(0.0, 0.0, -2.23606797749979))
        )
        self.assertAlmostEqual(inters6[1].radius, 1.9999999999999998)

        # test 7 sphere tangent to cylinder.
        spherical_surface = surfaces.SphericalSurface3D(volmdlr.OXYZ, 3)
        cylindrical_surface = surfaces.CylindricalSurface3D(
            volmdlr.OXYZ.translation(volmdlr.X3D), 2
        )
        inters7 = spherical_surface.surface_intersections(cylindrical_surface)
        self.assertAlmostEqual(inters7[0].length(), 7.006825334082508)
        self.assertAlmostEqual(inters7[1].length(), 7.006825334082508)
        self.assertAlmostEqual(inters7[2].length(), 7.006825334082508)
        self.assertAlmostEqual(inters7[3].length(), 7.006825334082508)

    def test_cylindricalsurface_intersections(self):
        cylindrical_surface1 = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 2)
        cylindrical_surface2 = surfaces.CylindricalSurface3D(volmdlr.OYZX, 1)
        inters = cylindrical_surface1.surface_intersections(cylindrical_surface2)
        expected_lengths1 = [6.393300778078848, 6.393300265079942]
        for intersection, expected_length in zip(inters, expected_lengths1):
            self.assertAlmostEqual(intersection.length(), expected_length)
        cylindrical_surface2 = surfaces.CylindricalSurface3D(
            volmdlr.OYZX.rotation(volmdlr.O3D, volmdlr.Y3D, math.pi / 4), 1)

        # test 2
        inters = cylindrical_surface1.surface_intersections(cylindrical_surface2)
        expected_lengths2 = [7.767042217039914, 7.767042239472898]
        for intersection, expected_length in zip(inters, expected_lengths2):
            self.assertAlmostEqual(intersection.length(), expected_length, 6)

        # test 3
        cylindrical_surface2 = surfaces.CylindricalSurface3D(volmdlr.OXYZ.translation(volmdlr.X3D * .5), 2)
        inters = cylindrical_surface1.surface_intersections(cylindrical_surface2)
        self.assertTrue(len(inters), 2)
        self.assertTrue(inters[0].point1.is_close(volmdlr.Point3D(0.25, -1.984313483298, 0.0)))
        self.assertTrue(inters[0].point2.is_close(volmdlr.Point3D(0.25, -1.984313483298, 1.0)))
        self.assertTrue(inters[1].point1.is_close(volmdlr.Point3D(0.25, 1.984313483298, 0.0)))
        self.assertTrue(inters[1].point2.is_close(volmdlr.Point3D(0.25, 1.984313483298, 1.0)))

        # test 4
        cylindrical_surface2 = surfaces.CylindricalSurface3D(volmdlr.OXYZ.translation(volmdlr.X3D * .5), .8)
        inters = cylindrical_surface1.surface_intersections(cylindrical_surface2)
        self.assertFalse(inters)


if __name__ == "__main__":
    unittest.main()
