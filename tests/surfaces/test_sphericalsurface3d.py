import unittest
import math
import volmdlr
from volmdlr import surfaces, wires, edges, curves


class TestSphericalSurface3D(unittest.TestCase):
    surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)

    def test_arc3d_to_2d(self):
        arc_with_two_singularities = edges.Arc3D.load_from_file(
            "surfaces/objects_spherical_tests/arc_with_two_singularities.json")
        test_arc = self.surface3d.arc3d_to_2d(arc_with_two_singularities)
        self.assertEqual(len(test_arc), 5)
        self.assertTrue(test_arc[0].start.is_close(volmdlr.Point2D(0, 0)))
        self.assertTrue(test_arc[-1].end.is_close(volmdlr.Point2D(0.0, 0.8975979010256554)))

    def test_contour3d_to_2d(self):
        surface = surfaces.SphericalSurface3D.load_from_file("surfaces/objects_spherical_tests/sphericalsurface1.json")
        contour = wires.Contour3D.load_from_file("surfaces/objects_spherical_tests/spericalsurface1_contour0.json")
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 6)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 4.107527949001648, 2)

        surface = surfaces.SphericalSurface3D.load_from_file(
            "surfaces/objects_spherical_tests/spherical_surface_arc3d_to_2d.json")
        contour = wires.Contour3D.load_from_file(
            "surfaces/objects_spherical_tests/spherical_surface_arc3d_to_2d_contour3d.json")
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 4)
        self.assertTrue(contour2d.is_ordered(1e-2))
        self.assertAlmostEqual(contour2d.area(), 1.7779412219307336, 2)

        surface = surfaces.SphericalSurface3D.load_from_file(
            "surfaces/objects_spherical_tests/buggy_contour3d_to_2d_surface.json")
        contour = wires.Contour3D.load_from_file(
            "surfaces/objects_spherical_tests/buggy_contour3d_to_2d_contour.json")
        contour2d = surface.contour3d_to_2d(contour)

        self.assertTrue(contour2d.is_ordered(1e-2))
        self.assertAlmostEqual(contour2d.area(), 0.028684788284169843, 2)

        surface = surfaces.SphericalSurface3D.load_from_file(
            "surfaces/objects_spherical_tests/contour3d_to_2d_surface_bspline_with_singularity.json")
        contour = wires.Contour3D.load_from_file(
            "surfaces/objects_spherical_tests/contour3d_to_2d_contour_bspline_with_singularity.json")
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 4)
        self.assertTrue(contour2d.is_ordered(1e-2))
        self.assertAlmostEqual(contour2d.area(), 1.1836145679685492, 2)

        surface = surfaces.SphericalSurface3D.load_from_file(
            "surfaces/objects_spherical_tests/repair_primitives2d_periodicity_surface.json")
        contour = wires.Contour3D.load_from_file(
            "surfaces/objects_spherical_tests/repair_primitives2d_periodicity_contour.json")
        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 4)
        self.assertTrue(contour2d.is_ordered(1e-2))
        self.assertAlmostEqual(contour2d.area(), 0.6254993351001795, 2)

        contour_left_side = wires.Contour3D.load_from_file("surfaces/objects_spherical_tests/contour_left_side.json")
        test = self.surface3d.contour3d_to_2d(contour_left_side)
        theta_min, theta_max, _, _ = test.bounding_rectangle.bounds()
        self.assertEqual(theta_min, -math.pi)
        self.assertEqual(theta_max, 0)
        contour_rigth_side = wires.Contour3D.load_from_file("surfaces/objects_spherical_tests/contour_rigth_side.json")
        test = self.surface3d.contour3d_to_2d(contour_rigth_side)
        theta_min, theta_max, _, _ = test.bounding_rectangle.bounds()
        self.assertEqual(theta_min, 0)
        self.assertEqual(theta_max, math.pi)
        contour_left_side = wires.Contour3D.load_from_file("surfaces/objects_spherical_tests/contour_upper_side.json")
        test = self.surface3d.contour3d_to_2d(contour_left_side)
        _, _, phi_min, phi_max = test.bounding_rectangle.bounds()
        self.assertEqual(phi_min, 0)
        self.assertEqual(phi_max, 0.5 * math.pi)
        contour_rigth_side = wires.Contour3D.load_from_file("surfaces/objects_spherical_tests/contour_lower_side.json")
        test = self.surface3d.contour3d_to_2d(contour_rigth_side)
        _, _, phi_min, phi_max = test.bounding_rectangle.bounds()
        self.assertEqual(phi_min, -0.5 * math.pi)
        self.assertEqual(phi_max, 0)

        contour_any_direction_upper = wires.Contour3D.load_from_file(
            "surfaces/objects_spherical_tests/contour_any_direction_upper_side.json")
        test = self.surface3d.contour3d_to_2d(contour_any_direction_upper)
        theta_min, theta_max, phi_min, phi_max = test.bounding_rectangle.bounds()
        self.assertAlmostEqual(theta_min, 0, 6)
        self.assertAlmostEqual(theta_max, math.pi, 4)
        self.assertAlmostEqual(phi_max, 0.5 * math.pi, 3)

    def test_plane_intersection(self):
        frame = volmdlr.Frame3D.from_point_and_vector(volmdlr.Point3D(0, 0.5, 0.5),
                                                      volmdlr.Vector3D(0, 1 / math.sqrt(2), 1 / math.sqrt(2)),
                                                      volmdlr.Z3D)
        plane = surfaces.Plane3D(frame)
        fullarc = self.surface3d.plane_intersection(plane)[0]
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


if __name__ == '__main__':
    unittest.main()
