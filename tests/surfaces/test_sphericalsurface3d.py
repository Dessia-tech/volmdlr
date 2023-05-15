import unittest
import math
import volmdlr
from volmdlr import surfaces, wires, edges


class TestSphericalSurface3D(unittest.TestCase):
    surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)


    def arc3d_to_2d(self):
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
        self.assertAlmostEqual(theta_max, math.pi, 6)
        self.assertAlmostEqual(phi_max, 0.5 * math.pi, 3)

    def test_plane_intersection(self):
        frame = volmdlr.Frame3D.from_point_and_vector(volmdlr.Point3D(0, 0.5, 0.5),
                                                      volmdlr.Vector3D(0, 1 / math.sqrt(2), 1 / math.sqrt(2)),
                                                      volmdlr.Z3D)
        plane = surfaces.Plane3D(frame)
        circle = self.surface3d.plane_intersection(plane)[0]
        self.assertAlmostEqual(circle.radius, 1/math.sqrt(2))
        self.assertTrue(circle.normal.is_colinear_to(frame.w))


if __name__ == '__main__':
    unittest.main()
