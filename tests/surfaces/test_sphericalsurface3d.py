import unittest
import math
import os
import numpy as np
import volmdlr
from volmdlr import surfaces, wires, edges, curves


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_spherical_tests')


class TestSphericalSurface3D(unittest.TestCase):
    surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)

    def test_parametric_points_to_3d(self):
        parametric_points = np.array([[0.0, 0.0], [0.0, 0.5 * math.pi], [0.0, math.pi], [0.0, 1.5 * math.pi],
                                      [0.5 * math.pi, 0.0], [0.5 * math.pi, 0.5 * math.pi],
                                      [0.5 * math.pi, math.pi], [0.5 * math.pi, 1.5 * math.pi],
                                      [math.pi, 0.0], [math.pi, 0.5 * math.pi],
                                      [math.pi, math.pi], [math.pi, 1.5 * math.pi],
                                      [1.5 * math.pi, 0.0], [1.5 * math.pi, 0.5 * math.pi],
                                      [1.5 * math.pi, math.pi], [1.5 * math.pi, 1.5 * math.pi]])
        points3d = self.surface3d.parametric_points_to_3d(parametric_points)
        expected_points = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [-1.0, 0.0, 0.0], [0.0, 0.0, -1.0],
                                    [0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0], [0.0, 0.0, -1.0],
                                    [-1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 0.0, -1.0],
                                    [0.0, -1.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [0.0, 0.0, -1.0]])
        for point, expected_point in zip(points3d, expected_points):
            self.assertAlmostEqual(np.linalg.norm(point - expected_point), 0.0)

    def test_arc3d_to_2d(self):
        arc_with_two_singularities = edges.Arc3D.from_json(
            os.path.join(folder, "arc_with_two_singularities.json"))
        test_arc = self.surface3d.arc3d_to_2d(arc_with_two_singularities)
        self.assertEqual(len(test_arc), 5)
        self.assertTrue(test_arc[0].start.is_close(volmdlr.Point2D(0, 0)))
        self.assertTrue(test_arc[-1].end.is_close(volmdlr.Point2D(0.0, 0.8975979010256554)))

    def test_contour3d_to_2d(self):
        surface = surfaces.SphericalSurface3D.from_json(os.path.join(folder, "sphericalsurface1.json"))
        contour = wires.Contour3D.from_json(os.path.join(folder, "spericalsurface1_contour0.json"))
        contour2d, primitives_mapping = surface.contour3d_to_2d(contour, return_primitives_mapping=True)
        self.assertEqual(len(contour2d.primitives), 6)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 4.107527949001648, 2)
        self.assertEqual(len(primitives_mapping), len(contour.primitives))
        self.assertIsNone(primitives_mapping.get(contour2d.primitives[3]))
        self.assertEqual(contour.primitives[0], primitives_mapping.get(contour2d.primitives[0]))
        self.assertEqual(contour.primitives[-1], primitives_mapping.get(contour2d.primitives[-1]))

        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "spherical_surface_arc3d_to_2d.json"))
        contour = wires.Contour3D.from_json(
            os.path.join(folder, "spherical_surface_arc3d_to_2d_contour3d.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 4)
        self.assertTrue(contour2d.is_ordered(1e-2))
        self.assertAlmostEqual(contour2d.area(), 1.7779412219307336, 2)

        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "buggy_contour3d_to_2d_surface.json"))
        contour = wires.Contour3D.from_json(
            os.path.join(folder, "buggy_contour3d_to_2d_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)

        self.assertTrue(contour2d.is_ordered(1e-2))
        self.assertAlmostEqual(contour2d.area(), 0.028684788284169843, 2)

        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "contour3d_to_2d_surface_bspline_with_singularity.json"))
        contour = wires.Contour3D.from_json(
            os.path.join(folder, "contour3d_to_2d_contour_bspline_with_singularity.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 4)
        self.assertTrue(contour2d.is_ordered(1e-2))
        self.assertAlmostEqual(contour2d.area(), 1.1836145679685492, 2)

        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "repair_primitives2d_periodicity_surface.json"))
        contour = wires.Contour3D.from_json(
            os.path.join(folder, "repair_primitives2d_periodicity_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 4)
        self.assertTrue(contour2d.is_ordered(1e-2))
        self.assertAlmostEqual(contour2d.area(), 0.6254993351001795, 2)

        contour_left_side = wires.Contour3D.from_json(os.path.join(folder, "contour_left_side.json"))
        test = self.surface3d.contour3d_to_2d(contour_left_side)
        theta_min, theta_max, _, _ = test.bounding_rectangle.bounds()
        self.assertEqual(theta_min, -math.pi)
        self.assertEqual(theta_max, 0)
        contour_rigth_side = wires.Contour3D.from_json(os.path.join(folder, "contour_rigth_side.json"))
        test = self.surface3d.contour3d_to_2d(contour_rigth_side)
        theta_min, theta_max, _, _ = test.bounding_rectangle.bounds()
        self.assertEqual(theta_min, 0)
        self.assertEqual(theta_max, math.pi)
        contour_left_side = wires.Contour3D.from_json(os.path.join(folder, "contour_upper_side.json"))
        test = self.surface3d.contour3d_to_2d(contour_left_side)
        _, _, phi_min, phi_max = test.bounding_rectangle.bounds()
        self.assertEqual(phi_min, 0)
        self.assertEqual(phi_max, 0.5 * math.pi)
        contour_rigth_side = wires.Contour3D.from_json(os.path.join(folder, "contour_lower_side.json"))
        test = self.surface3d.contour3d_to_2d(contour_rigth_side)
        _, _, phi_min, phi_max = test.bounding_rectangle.bounds()
        self.assertEqual(phi_min, -0.5 * math.pi)
        self.assertEqual(phi_max, 0)

        contour_any_direction_upper = wires.Contour3D.from_json(
            os.path.join(folder, "contour_any_direction_upper_side.json"))
        test = self.surface3d.contour3d_to_2d(contour_any_direction_upper)
        theta_min, theta_max, phi_min, phi_max = test.bounding_rectangle.bounds()
        self.assertAlmostEqual(theta_min, 0, 6)
        self.assertAlmostEqual(theta_max, math.pi, 4)
        self.assertAlmostEqual(phi_max, 0.5 * math.pi, 3)

        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "test_sphericalsurface_repair_periodicity_surface.json"))
        contour = wires.Contour3D.from_json(
            os.path.join(folder, "test_sphericalsurface_repair_periodicity_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 6)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 6.129921072323977, 2)

        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "test_2_sphericalsurface_repair_periodicity_surface.json"))
        contour = wires.Contour3D.from_json(
            os.path.join(folder, "test_2_sphericalsurface_repair_periodicity_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 8)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 2.1665348983853794, 2)

        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "sphericalsurface_contour3d_to_2d_positive_singularity.json"))
        contour = wires.Contour3D.from_json(
            os.path.join(folder, "sphericalsurface_contour3d_to_2d_positive_singularity_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 4)
        self.assertTrue(contour2d.is_ordered(0.0001))
        self.assertAlmostEqual(contour2d.area(), 0.10067063484647809, 2)

    def test_plane_intersections(self):
        frame = volmdlr.Frame3D.from_point_and_vector(volmdlr.Point3D(0, 0.5, 0.5),
                                                      volmdlr.Vector3D(0, 1 / math.sqrt(2), 1 / math.sqrt(2)),
                                                      volmdlr.Z3D)
        plane = surfaces.Plane3D(frame)
        fullarc = self.surface3d.plane_intersections(plane)[0]
        self.assertAlmostEqual(fullarc.circle.radius, 1/math.sqrt(2))
        self.assertTrue(fullarc.circle.normal.is_colinear_to(frame.w))

    def test_line_intersections(self):
        spherical_surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)
        line = curves.Line3D(volmdlr.Point3D(-.01, -.01, -.01), volmdlr.Point3D(0.8, 0.8, 0.8))
        line_intersections = spherical_surface3d.line_intersections(line)
        self.assertEqual(len(line_intersections), 2)
        self.assertTrue(line_intersections[0].is_close(
            volmdlr.Point3D(0.5773502691896257, 0.5773502691896257, 0.5773502691896257)))
        self.assertTrue(line_intersections[1].is_close(
            volmdlr.Point3D(-0.5773502691896257, -0.5773502691896257, -0.5773502691896257)))

        line2 = curves.Line3D(volmdlr.Point3D(0, 1, -0.5), volmdlr.Point3D(0, 1, 0.5))
        line_intersections2 = spherical_surface3d.line_intersections(line2)
        self.assertEqual(len(line_intersections2), 1)
        self.assertTrue(line_intersections2[0].is_close(volmdlr.Point3D(0.0, 1.0, 0.0)))

        line3 = curves.Line3D(volmdlr.Point3D(0, 2, -0.5), volmdlr.Point3D(0, 2, 0.5))
        self.assertFalse(spherical_surface3d.line_intersections(line3))

    def test_linesegment_intersections(self):
        spherical_surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)
        linesegment = edges.LineSegment3D(volmdlr.Point3D(-.01, -.01, -.01), volmdlr.Point3D(0.5, -0.4, 0.8))
        linesegment_intersections = spherical_surface3d.linesegment_intersections(linesegment)
        self.assertEqual(len(linesegment_intersections), 1)
        self.assertTrue(linesegment_intersections[0].is_close(
            volmdlr.Point3D(0.4878134615813934, -0.3906808823857714, 0.7806449095704483)))

        linesegment2 = edges.LineSegment3D(volmdlr.Point3D(-0.8, -0.8, -0.8), volmdlr.Point3D(0.8, 0.8, 0.8))
        linesegment2_intersections = spherical_surface3d.linesegment_intersections(linesegment2)
        self.assertEqual(len(linesegment2_intersections), 2)
        self.assertTrue(linesegment2_intersections[0].is_close(
            volmdlr.Point3D(0.5773502691896257, 0.5773502691896257, 0.5773502691896257)))
        self.assertTrue(linesegment2_intersections[1].is_close(
            volmdlr.Point3D(-0.5773502691896257, -0.5773502691896257, -0.5773502691896257)))

        linesegment3 = edges.LineSegment3D(volmdlr.Point3D(-.01, -.01, -.01), volmdlr.Point3D(0.3, -0.4, 0.5))
        self.assertFalse(spherical_surface3d.linesegment_intersections(linesegment3))

        linesegment4 = edges.LineSegment3D(volmdlr.Point3D(-.01, -1.5, -.01), volmdlr.Point3D(0.3, -1.3, 0.5))
        self.assertFalse(spherical_surface3d.linesegment_intersections(linesegment4))

    def test_circle_intersections(self):
        spherical_surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)

        # test1
        circle = curves.Circle3D(volmdlr.OXYZ.translation(volmdlr.Vector3D(.5, .5, 2)), .5)
        circle_intersections = spherical_surface3d.circle_intersections(circle)
        self.assertFalse(circle_intersections)

        # test2
        circle = curves.Circle3D(volmdlr.OXYZ.translation(volmdlr.Vector3D(.5, .5, .5)), .5)
        circle_intersections = spherical_surface3d.circle_intersections(circle)
        self.assertEqual(len(circle_intersections), 2)
        self.assertTrue(circle_intersections[0], volmdlr.Point3D(0.853553390593, 0.146446609407, 0.5))
        self.assertTrue(circle_intersections[1], volmdlr.Point3D(0.146446609407, 0.853553390593, 0.5))

    def test_arc_intersections(self):
        # test1
        spherical_surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        vector3 = vector1.cross(vector2)
        frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)
        circle = curves.Circle3D(frame.translation(volmdlr.Vector3D(.5, .5, .5)), .5)

        point1 = circle.point_at_abscissa(0.2)
        point2 = circle.point_at_abscissa(1.5)
        arc = edges.Arc3D(circle, point1, point2)
        arc_intersections = spherical_surface3d.arc_intersections(arc)
        self.assertEqual(len(arc_intersections), 1)
        self.assertTrue(arc_intersections[0].is_close(volmdlr.Point3D(0.908248290464, 0.295875854768, 0.295875854768)))
        # test 2:
        point2 = circle.point_at_abscissa(2.5)
        arc = edges.Arc3D(circle, point1, point2)
        arc_intersections = spherical_surface3d.arc_intersections(arc)
        self.assertEqual(len(arc_intersections), 2)
        self.assertTrue(arc_intersections[0].is_close(volmdlr.Point3D(0.091751709536, 0.704124145232, 0.704124145232)))
        self.assertTrue(arc_intersections[1].is_close(volmdlr.Point3D(0.908248290464, 0.295875854768, 0.295875854768)))

    def test_arcellipse_intersections(self):
        spherical_surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        vector3 = vector1.cross(vector2)
        frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)
        ellipse = curves.Ellipse3D(1, .5, frame.translation(volmdlr.Vector3D(.5, .5, .5)))

        point1 = ellipse.point_at_abscissa(0.2)
        # test 2
        point2 = ellipse.point_at_abscissa(2.5)
        arcellipse = edges.ArcEllipse3D(ellipse, point1, point2)
        arcellipse_intersections = spherical_surface3d.arcellipse_intersections(arcellipse)
        self.assertEqual(len(arcellipse_intersections), 1)
        self.assertTrue(arcellipse_intersections[0].is_close(
            volmdlr.Point3D(0.908248233153, 0.295875797457, 0.295875797457)))

        # test 3
        point2 = ellipse.point_at_abscissa(4.5)
        arcellipse = edges.ArcEllipse3D(ellipse, point1, point2)
        arcellipse_intersections = spherical_surface3d.arcellipse_intersections(arcellipse)
        self.assertEqual(len(arcellipse_intersections), 2)
        self.assertTrue(arcellipse_intersections[0].is_close(
            volmdlr.Point3D(0.908248233153, 0.295875797457, 0.295875797457)))
        self.assertTrue(arcellipse_intersections[1].is_close(
            volmdlr.Point3D(0.091751652225, 0.704124087921, 0.704124087921)))

    def test_sphericalsurface_intersections(self):
        #test1
        spherical_surface1 = surfaces.SphericalSurface3D(
            volmdlr.OXYZ.translation(volmdlr.Vector3D(0.5, 0.5, 0)), 2)

        spherical_surface2 = spherical_surface1.translation(volmdlr.Vector3D(1, 1, 1))

        inters = spherical_surface1.surface_intersections(spherical_surface2)
        self.assertEqual(len(inters), 1)
        self.assertAlmostEqual(inters[0].length(), 11.327173398039175)

        #test2
        spherical_surface2 = surfaces.SphericalSurface3D(
            volmdlr.OXYZ.translation(volmdlr.Vector3D(0.5, 0.5, 0)), 1)
        spherical_surface2 = spherical_surface2.translation(volmdlr.Vector3D(1, 1, 1))
        inters = spherical_surface1.surface_intersections(spherical_surface2)
        self.assertEqual(len(inters), 1)
        self.assertAlmostEqual(inters[0].length(), 6.283185306688713)

    def test_normal_at_point(self):
        spherical_surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)
        point = volmdlr.Point3D(0.908248233153, 0.295875797457, 0.295875797457)
        normal = spherical_surface3d.normal_at_point(point)
        self.assertTrue(normal.is_close(volmdlr.Vector3D(0.9082483187228496, 0.2958758253326914, 0.30974418137711035)))


if __name__ == '__main__':
    unittest.main()
