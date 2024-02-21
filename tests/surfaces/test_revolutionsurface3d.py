import unittest
import math
import os
import numpy as np
import volmdlr
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.faces
from volmdlr import surfaces


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_revolution_tests')


class TestRevolutionSurface3D(unittest.TestCase):
    linesegment = vme.LineSegment3D(volmdlr.Point3D(0.5, 0, 0), volmdlr.Point3D(0.5, 0, 0.5))
    arc = vme.Arc3D.from_3_points(volmdlr.Point3D(0.5, 0, 0.5),
                    volmdlr.Point3D(0.3 + 0.2 * math.cos(math.pi / 6), 0, 0.5 + 0.2 * math.sin(math.pi / 6)),
                    volmdlr.Point3D(0.3 + 0.2 * math.cos(math.pi / 3), 0, 0.5 + 0.2 * math.sin(math.pi / 3)))

    axis_point = volmdlr.O3D
    axis = volmdlr.Z3D
    surface = surfaces.RevolutionSurface3D(arc, axis_point, axis)

    def test_parametric_points_to_3d(self):
        parametric_points = np.array([[0.0, 0.0], [0.5 * math.pi, 0.0], [math.pi, 0.0], [1.5 * math.pi, 0.0],
                                      [0.0, 0.2], [0.5 * math.pi, 0.2], [math.pi, 0.2], [1.5 * math.pi, 0.2]])
        points3d = self.surface.parametric_points_to_3d(parametric_points)
        expected_points = np.array([[0.5, 0.0, 0.5], [0.0, 0.5, 0.5], [-0.5, 0.0, 0.5], [0.0, -0.5, 0.5],
                                    [0.40806046117362793, 0.0, 0.6682941969615793],
                                    [0.0, 0.40806046117362793, 0.6682941969615792],
                                    [-0.40806046117362793, 0.0, 0.6682941969615793], [0.0, -0.40806046117362793, 0.6682941969615793]])
        for point, expected_point in zip(points3d, expected_points):
            self.assertAlmostEqual(np.linalg.norm(point - expected_point), 0.0)

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
        surface = surfaces.RevolutionSurface3D.from_json(os.path.join(folder, "revolution_surface_bug_0.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "revolution_contour_bug_0.json"))
        arc = contour.primitives[3]
        linesegment2d = surface.arc3d_to_2d(arc)[0]
        self.assertTrue(linesegment2d.start.is_close(volmdlr.Point2D(-2.6665730021726306, 0.0003141532401719152)))
        self.assertTrue(linesegment2d.end.is_close(volmdlr.Point2D(-2.6665730021726324, 0)))

        surface = surfaces.RevolutionSurface3D.from_json(os.path.join(folder, "revolution_surface_bug_1.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "revolution_contour_bug_1.json"))
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
        rev1 = surfaces.RevolutionSurface3D.from_json(
            os.path.join(folder, "revolutionsurface_simplify_spherical.json"))
        rev2 = surfaces.RevolutionSurface3D.from_json(
            os.path.join(folder, "revolutionsurface_simplify_conical.json"))
        rev3 = surfaces.RevolutionSurface3D.from_json(
            os.path.join(folder, "revolutionsurface_simplify_cylindrical.json"))

        sphere = rev1.simplify()
        self.assertTrue(isinstance(sphere, surfaces.SphericalSurface3D))

        cone = rev2.simplify()
        self.assertTrue(isinstance(cone, surfaces.ConicalSurface3D))

        cylinder = rev3.simplify()
        self.assertTrue(isinstance(cylinder, surfaces.CylindricalSurface3D))

    def test_linesegment2d_to_3d(self):
        surface = surfaces.RevolutionSurface3D.from_json(
            os.path.join(folder, "revolutionsurface_linesegment2d_to_3d.json"))
        linesegment1 = vme.LineSegment2D.from_json(os.path.join(folder, "linesegment2d_arc3d.json"))
        arc = surface.linesegment2d_to_3d(linesegment1)[0]
        self.assertAlmostEqual(arc.circle.radius, 0.02404221842799788)

        linesegment2 = vme.LineSegment2D.from_json(
            os.path.join(folder, "linesegment2d_rotated_primitive.json"))
        arc = surface.linesegment2d_to_3d(linesegment2)[0]
        self.assertAlmostEqual(arc.circle.radius, 0.022500000035448893)
        self.assertAlmostEqual(arc.angle, 0.7195087615152496, 5)

        linesegment3 = vme.LineSegment2D.from_json(
            os.path.join(folder, "linesegment2d_split_primitive.json"))
        arc = surface.linesegment2d_to_3d(linesegment3)[0]
        self.assertAlmostEqual(arc.circle.radius, 0.022500000035448893)
        self.assertAlmostEqual(arc.angle, 0.15581712793343738)

        surface = surfaces.RevolutionSurface3D.from_json(
            os.path.join(folder, "revolutionsurface_periodical_linesegment2d_to_3d.json"))
        linesegment = vme.LineSegment2D.from_json(
            os.path.join(folder, "revolutionsurface_periodical_linesegment2d_to_3d_linesegment2d.json"))
        arc = surface.linesegment2d_to_3d(linesegment)[0]
        self.assertAlmostEqual(arc.radius, 0.017000000000019)
        self.assertTrue(arc.center.is_close(volmdlr.Point3D(0.0, 0.007299999999984744, -8.104628079745562e-19)))

    def test_contour3d_to_2d(self):
        surface = surfaces.RevolutionSurface3D.from_json(os.path.join(folder, "revolutionface_surface.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "revolutionface_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.031887165433924704 * 0.5 * math.pi, 2)

        surface = surfaces.RevolutionSurface3D.from_json(os.path.join(folder, "revolutionface_surface_1.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "revolutionface_contour_1.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())

        surface = surfaces.RevolutionSurface3D.from_json(os.path.join(folder, "revolutionface_surface_2.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "revolutionface_contour_2.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.00031415327300491437 * math.pi, 2)

        surface = surfaces.RevolutionSurface3D.from_json(os.path.join(folder,
                                                                           "revolutionsurface_with_singularity.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "revolutionsurface_with_singularity_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), surface.edge.length() * math.pi, 2)

        surface = surfaces.RevolutionSurface3D.from_json(os.path.join(folder,
                                                                        "revolutionsurface_with_singularity_1.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder,
                                                            "revolutionsurface_with_singularity_contour_1.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), surface.edge.length() * math.pi, 2)


    def test_arc3d_to_2d(self):
        surface = surfaces.RevolutionSurface3D.from_json(os.path.join(folder, "arc3d_to_2d_surface.json"))
        arc3d = vme.Arc3D.from_json(os.path.join(folder, "arc3d_to_2d_arc3d.json"))
        brep = surface.arc3d_to_2d(arc3d)[0]
        self.assertTrue(brep.start.is_close(volmdlr.Point2D(-math.pi, 0.0)))
        self.assertTrue(brep.end.is_close(volmdlr.Point2D(-math.pi, 0.0038322109949349634)))

        surface = surfaces.RevolutionSurface3D.from_json(os.path.join(folder, "arc3d_to_2d_surface_2.json"))
        arc3d = vme.Arc3D.from_json(os.path.join(folder, "arc3d_to_2d_arc3d_2.json"))
        brep = surface.arc3d_to_2d(arc3d)[0]
        self.assertTrue(brep.start.is_close(volmdlr.Point2D(-math.pi, 0.0038322109949349634)))
        self.assertTrue(brep.end.is_close(volmdlr.Point2D(-math.pi, 0.0)))

    def test_v_iso(self):
        surface = surfaces.RevolutionSurface3D.from_json(
            os.path.join(folder, "revolutionsurface_periodical_linesegment2d_to_3d.json"))
        v = 0.023550776716126855
        arc = surface.v_iso(v)
        self.assertAlmostEqual(arc.radius, 0.017000000000019)
        self.assertTrue(arc.center.is_close(volmdlr.Point3D(0.0, 0.007299999999984744, -8.104628079745562e-19)))




if __name__ == '__main__':
    unittest.main()
