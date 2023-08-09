import unittest
import math

import volmdlr
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.faces
from volmdlr import surfaces


class TestRevolutionSurface3D(unittest.TestCase):
    linesegment = vme.LineSegment3D(volmdlr.Point3D(0.5, 0, 0), volmdlr.Point3D(0.5, 0, 0.5))
    arc = vme.Arc3D.from_3_points(volmdlr.Point3D(0.5, 0, 0.5),
                    volmdlr.Point3D(0.3 + 0.2 * math.cos(math.pi / 6), 0, 0.5 + 0.2 * math.sin(math.pi / 6)),
                    volmdlr.Point3D(0.3 + 0.2 * math.cos(math.pi / 3), 0, 0.5 + 0.2 * math.sin(math.pi / 3)))

    axis_point = volmdlr.O3D
    axis = volmdlr.Z3D
    surface = surfaces.RevolutionSurface3D(arc, axis_point, axis)

    def test_point2d_to_3d(self):
        surface = surfaces.RevolutionSurface3D(self.arc, self.axis_point, self.axis)

        point2d = volmdlr.Point2D(math.pi, 0.2)
        point3d = surface.point2d_to_3d(point2d)
        expected_point3d = volmdlr.Point3D(-0.4080604, 0, 0.66829419)

        self.assertTrue(point3d.is_close(expected_point3d))

    def test_point3d_to_2d(self):
        surface = surfaces.RevolutionSurface3D(self.arc, self.axis_point, self.axis)

        point3d = volmdlr.Point3D(-0.4080604, 0, 0.66829419)
        point2d = surface.point3d_to_2d(point3d)
        expected_point2d = volmdlr.Point2D(math.pi, 0.2)

        self.assertTrue(point2d.is_close(expected_point2d))

    def test_rectangular_cut(self):
        surface = surfaces.RevolutionSurface3D(edge=self.arc, axis_point=self.axis_point, axis=self.axis)
        rectangular_cut = volmdlr.faces.RevolutionFace3D.from_surface_rectangular_cut(
            surface, 0, volmdlr.TWO_PI, 0, 1)
        self.assertEqual(rectangular_cut.surface2d.area(), volmdlr.TWO_PI)

    def arc3d_to_2d(self):
        surface = surfaces.RevolutionSurface3D.load_from_file(
            "surfaces/objects_revolution_tests/revolution_surface_bug_0.json")
        contour = vmw.Contour3D.load_from_file("surfaces/objects_revolution_tests/revolution_contour_bug_0.json")
        arc = contour.primitives[3]
        linesegment2d = surface.arc3d_to_2d(arc)[0]
        self.assertTrue(linesegment2d.start.is_close(volmdlr.Point2D(-2.6665730021726306, 0.0003141532401719152)))
        self.assertTrue(linesegment2d.end.is_close(volmdlr.Point2D(-2.6665730021726324, 0)))

        surface = surfaces.RevolutionSurface3D.load_from_file(
            "surfaces/objects_revolution_tests/revolution_surface_bug_1.json")
        contour = vmw.Contour3D.load_from_file("surfaces/objects_revolution_tests/revolution_contour_bug_1.json")
        arc = contour.primitives[3]
        linesegment2d = surface.arc3d_to_2d(arc)[0]
        self.assertTrue(linesegment2d.start.is_close(volmdlr.Point2D(1.0582040468439544, 0.0)))
        self.assertTrue(linesegment2d.end.is_close(volmdlr.Point2D(6.0956438652626925, 0.0)))

    def test_frame_mapping(self):
        surface = self.surface
        new_frame = volmdlr.Frame3D(volmdlr.Point3D(0, 1, 0), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        new_surface = surface.frame_mapping(new_frame, "old")
        self.assertEqual(new_surface.edge.start.y, 1)
        self.assertTrue(new_surface.frame.origin.is_close(volmdlr.Point3D(0, 1, 0)))

    def test_simplify(self):
        rev1 = surfaces.RevolutionSurface3D.load_from_file(
            "surfaces/objects_revolution_tests/revolutionsurface_simplify_spherical.json")
        rev2 = surfaces.RevolutionSurface3D.load_from_file(
            "surfaces/objects_revolution_tests/revolutionsurface_simplify_conical.json")
        rev3 = surfaces.RevolutionSurface3D.load_from_file(
            "surfaces/objects_revolution_tests/revolutionsurface_simplify_cylindrical.json")

        sphere = rev1.simplify()
        self.assertTrue(isinstance(sphere, surfaces.SphericalSurface3D))

        cone = rev2.simplify()
        self.assertTrue(isinstance(cone, surfaces.ConicalSurface3D))

        cylinder = rev3.simplify()
        self.assertTrue(isinstance(cylinder, surfaces.CylindricalSurface3D))

    def test_linesegment2d_to_3d(self):
        surface = surfaces.RevolutionSurface3D.load_from_file(
            "surfaces/objects_revolution_tests/revolutionsurface_linesegment2d_to_3d.json")
        linesegment1 = vme.LineSegment2D.load_from_file("surfaces/objects_revolution_tests/linesegment2d_arc3d.json")
        arc = surface.linesegment2d_to_3d(linesegment1)[0]
        self.assertAlmostEqual(arc.circle.radius, 0.02404221842799788)

        linesegment2 = vme.LineSegment2D.load_from_file(
            "surfaces/objects_revolution_tests/linesegment2d_rotated_primitive.json")
        arc = surface.linesegment2d_to_3d(linesegment2)[0]
        self.assertAlmostEqual(arc.circle.radius, 0.022500000035448893)
        self.assertAlmostEqual(arc.angle, 0.7195087615152496, 5)

        linesegment3 = vme.LineSegment2D.load_from_file(
            "surfaces/objects_revolution_tests/linesegment2d_split_primitive.json")
        arc = surface.linesegment2d_to_3d(linesegment3)[0]
        self.assertAlmostEqual(arc.circle.radius, 0.022500000035448893)
        self.assertAlmostEqual(arc.angle, 0.15581712793343738)

    def test_contour3d_to_2d(self):
        surface = surfaces.RevolutionSurface3D.load_from_file(
            "surfaces/objects_revolution_tests/revolutionface_surface.json")
        contour = vmw.Contour3D.load_from_file("surfaces/objects_revolution_tests/revolutionface_contour.json")
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.031887165433924704 * 0.5 * math.pi, 2)

        surface = surfaces.RevolutionSurface3D.load_from_file(
            "surfaces/objects_revolution_tests/revolutionface_surface_1.json")
        contour = vmw.Contour3D.load_from_file("surfaces/objects_revolution_tests/revolutionface_contour_1.json")
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())

        surface = surfaces.RevolutionSurface3D.load_from_file(
            "surfaces/objects_revolution_tests/revolutionface_surface_2.json")
        contour = vmw.Contour3D.load_from_file("surfaces/objects_revolution_tests/revolutionface_contour_2.json")
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.00031415327300491437 * math.pi, 2)


if __name__ == '__main__':
    unittest.main()
