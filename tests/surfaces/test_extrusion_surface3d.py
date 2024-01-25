import unittest
import os
import numpy as np
import volmdlr
import volmdlr.edges as vme
import volmdlr.faces as vmf
import volmdlr.step as vms
import volmdlr.wires as vmw
from volmdlr import surfaces


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_extrusion_tests')


class TestExtrusionSurface3D(unittest.TestCase):
    control_points = [
        volmdlr.Point3D(-0.025917292, 0.002544355, 0.0),
        volmdlr.Point3D(-0.005449685, -0.007265721, 0.0),
        volmdlr.Point3D(0.0, 0.0, 0.0),
        volmdlr.Point3D(0.014457705000000001, -0.002636091, 0.0),
        volmdlr.Point3D(0.013503079, -0.014007147, 0.0)]
    edge = vme.BSplineCurve3D(3, control_points, [4, 1, 4], [0.0, 0.5, 1.0])
    surface = surfaces.ExtrusionSurface3D(edge, -volmdlr.Z3D)

    def test_parametric_points_to_3d(self):
        parametric_points = np.array([[0.0, 0.0], [0.25 * self.edge.length(), 0.0], [0.5 * self.edge.length(), 0.0],
                                      [0.75 * self.edge.length(), 0.0], [self.edge.length(), 0.0],
                                      [0.0, 1.0], [0.25 * self.edge.length(), 1.0], [0.5 * self.edge.length(), 1.0],
                                      [0.75 * self.edge.length(), 1.0], [self.edge.length(), 1.0],
                                      [0.0, -1.0], [0.25 * self.edge.length(), -1.0], [0.5 * self.edge.length(), -1.0],
                                      [0.75 * self.edge.length(), -1.0], [self.edge.length(), -1.0]])
        points3d = self.surface.parametric_points_to_3d(parametric_points)
        expected_points = np.array([[-0.025917292, 0.002544355, 0.0], [-0.006023608687500001, -0.0040783553125, 0.0],
                                    [0.0022520050000000005, -0.002475453, 0.0],
                                    [0.010101844562500002, -0.0035431261875, 0.0], [0.013503079, -0.014007147, 0.0],
                                    [-0.025917292, 0.002544355, -1.0], [-0.006023608687500001, -0.0040783553125, -1.0],
                                    [0.0022520050000000005, -0.002475453, -1.0],
                                    [0.010101844562500002, -0.0035431261875, -1.0], [0.013503079, -0.014007147, -1.0],
                                    [-0.025917292, 0.002544355, 1.0], [-0.006023608687500001, -0.0040783553125, 1.0],
                                    [0.0022520050000000005, -0.002475453, 1.0],
                                    [0.010101844562500002, -0.0035431261875, 1.0], [0.013503079, -0.014007147, 1.0]
                                    ])
        for point, expected_point in zip(points3d, expected_points):
            self.assertAlmostEqual(np.linalg.norm(point - expected_point), 0.0)

    def test_point2d_to_3d(self):
        point3d = self.surface.point2d_to_3d(volmdlr.Point2D(0.5 * self.edge.length(), 0.5))
        self.assertTrue(point3d.is_close(volmdlr.Point3D(0.002252005, -0.002475453, -0.5)))

    def test_point3d_to_2d(self):
        point2d_1 = self.surface.point3d_to_2d(self.edge.start)
        self.assertEqual(point2d_1, volmdlr.Point2D(0, 0))
        point2d_2 = self.surface.point3d_to_2d(self.edge.end)
        self.assertEqual(point2d_2, volmdlr.Point2D(self.edge.length(), 0))

        surface = surfaces.ExtrusionSurface3D.from_json(os.path.join(folder, "point3d_to_2d_surface.json"))
        point3d = volmdlr.Point3D(-0.1322515585788849, -0.14157776310991646, 0.0)
        point2d_3 = surface.point3d_to_2d(point3d)
        self.assertTrue(point2d_3.is_close(volmdlr.Point2D(0.8082617614489124, 0.0), 1e-3))

    def test_rectangular_cut(self):
        face = vmf.ExtrusionFace3D.from_surface_rectangular_cut(self.surface, 0, self.edge.length(), 0, 2)
        self.assertEqual(face.surface2d.area(), 2 * self.edge.length())

    def test_from_step(self):
        step = vms.Step.from_file(os.path.join(folder, "bspline_extruded_simple.step"))
        model = step.to_volume_model()
        extrusion_surface = model.primitives[0].primitives[0].surface3d
        self.assertEqual(extrusion_surface.direction, -volmdlr.Z3D)
        self.assertEqual(extrusion_surface.edge.degree, 3)
        self.assertEqual(extrusion_surface.edge.knot_multiplicities.tolist(), [4, 1, 4])

    def test_linesegment2d_to_3d(self):
        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusion_surface_undefined_direction_linesegment.json"))
        point1 = volmdlr.Point2D(0.9020984833336293, -0.08534036750789999)
        point2 = volmdlr.Point2D(0.9286913444016728, -0.07799341694)
        linesegment2d = vme.LineSegment2D(point1, point2)
        start3d = surface.point2d_to_3d(point1)
        end3d = surface.point2d_to_3d(point2)
        result = surface.linesegment2d_to_3d(linesegment2d)[0]
        self.assertTrue(result.start.is_close(start3d))
        self.assertTrue(result.end.is_close(end3d))
        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusion_surface_test_linesegment2d_to_3d.json"))
        linesegment2d = vme.LineSegment2D.from_json(os.path.join(folder, "linesegment2d_to_linesegment3d.json"))
        result = surface.linesegment2d_to_3d(linesegment2d)[0]
        self.assertIsInstance(result, vme.LineSegment3D)

        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusion_surface_test_linesegment2d_to_3d.json"))
        linesegment2d = vme.LineSegment2D.from_json(os.path.join(folder, "linesegment2d_to_self_edge.json"))
        result = surface.linesegment2d_to_3d(linesegment2d)[0]
        self.assertEqual(result, surface.edge)

    def test_arc3d_to_2d(self):
        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusion_surface_test_arc3d_to_2d.json"))
        contour3d = vmw.Contour3D.from_json(os.path.join(folder, "extrusion_contour_test_arc3d_to_2d.json"))
        arc3d = contour3d.primitives[2]
        result = surface.arc3d_to_2d(arc3d)[0]
        self.assertTrue(result.start.is_close(volmdlr.Point2D(0.0034138143320201525, 0.0032000000499998738)))
        self.assertTrue(result.end.is_close(volmdlr.Point2D(0.00046275860846800896, 0.0032000000499998738)))

    def test_bsplinecurve3d_to_2d(self):
        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "periodical_extrusionsurface.json"))
        bsplinecurve3d = vme.BSplineCurve3D.from_json(
            os.path.join(folder, "periodical_extrusionsurface_bsplinecurve.json"))
        result = surface.bsplinecurve3d_to_2d(bsplinecurve3d)[0]
        inverse_prof = surface.linesegment2d_to_3d(result)[0]
        self.assertTrue(result.start.is_close(volmdlr.Point2D(5.475029217377275, 0.02477709130796299)))
        self.assertTrue(result.end.is_close(volmdlr.Point2D(5.616985468553595, 0.024781167649391873)))
        self.assertTrue(inverse_prof.is_close(inverse_prof))

    def test_fullarcellipse3d_to_2d(self):
        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusionsurface_fullarcellipse3d_to_2d.json"))
        ellipse = vme.FullArcEllipse3D.from_json(
            os.path.join(folder, "extrusionsurface_fullarcellipse3d_to_2d_fullarcellipse3d.json"))
        result = surface.fullarcellipse3d_to_2d(ellipse)[0]
        self.assertTrue(result.start.is_close(volmdlr.Point2D(0.0, 0.01)))
        self.assertTrue(result.end.is_close(volmdlr.Point2D(0.025526998862788763, 0.01)))

    def test_frame_mapping(self):
        surface = self.surface
        new_frame = volmdlr.Frame3D(volmdlr.Point3D(0, 0, 1), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        new_surface = surface.frame_mapping(new_frame, "old")
        self.assertEqual(new_surface.edge.start.z, 1)
        self.assertTrue(new_surface.frame.origin.is_close(volmdlr.Point3D(-0.025917292, 0.002544355, 1.0)))

    def test_contour3d_to_2d(self):
        surface = surfaces.ExtrusionSurface3D.from_json(os.path.join(folder, "contour3d_to_2d_surface.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "contour3d_to_2d_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.00032168769592775094, 6)

        surface = surfaces.ExtrusionSurface3D.from_json(os.path.join(folder, "contour3d_to_2d_surface_2.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "contour3d_to_2d_contour_2.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.05992365409316021, 6)

        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusionsurface_edge_not_in_normal_plane.json"))
        contour = vmw.Contour3D.from_json(
            os.path.join(folder, "extrusionsurface_edge_not_in_normal_plane_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.00019036534467768707, 6)

        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusionsurface_with_small_edge.json"))
        contour = vmw.Contour3D.from_json(
            os.path.join(folder, "extrusionsurface_with_small_edge_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 3.74649557711703e-09, 10)

        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusionsurface_fullarcellipse.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "extrusionsurface_fullarcellipse_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.012120134592666365, 6)

        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusionsurface_fullarc.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "extrusionsurface_fullarc_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 2.0719721732132054e-06, 8)

        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusionsurface_fullarcellipse.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "extrusionsurface_fullarcellipse_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.012120134592666365, 2)

        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusionsurface_periodic.json"))
        contour = vmw.Contour3D.from_json(os.path.join(folder, "extrusionsurface_periodic_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 2.009851332304794e-06, 8)

        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "periodical_extrusionsurface_linesegment3d_to_2d.json"))
        contour = vmw.Contour3D.from_json(
            os.path.join(folder, "periodical_extrusionsurface_linesegment3d_to_2d_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)
        self.assertTrue(contour2d.is_ordered(1e-5))
        self.assertAlmostEqual(contour2d.area(), 0.007376809172328507, 2)


if __name__ == '__main__':
    unittest.main()
