import unittest
import volmdlr
import volmdlr.edges as vme
import volmdlr.faces as vmf
import volmdlr.step as vms
import volmdlr.wires as vmw
from volmdlr import surfaces


class TestExtrusionSurface3D(unittest.TestCase):
    control_points = [
        volmdlr.Point3D(-0.025917292, 0.002544355, 0.0),
        volmdlr.Point3D(-0.005449685, -0.007265721, 0.0),
        volmdlr.Point3D(0.0, 0.0, 0.0),
        volmdlr.Point3D(0.014457705000000001, -0.002636091, 0.0),
        volmdlr.Point3D(0.013503079, -0.014007147, 0.0)]
    edge = vme.BSplineCurve3D(3, control_points, [4, 1, 4], [0.0, 0.5, 1.0])
    surface = surfaces.ExtrusionSurface3D(edge, -volmdlr.Z3D)

    def test_point2d_to_3d(self):
        point3d = self.surface.point2d_to_3d(volmdlr.Point2D(0.5, 0.5))
        self.assertTrue(point3d.is_close(volmdlr.Point3D(0.002252005, -0.002475453, -0.5)))

    def test_point3d_to_2d(self):
        point2d_1 = self.surface.point3d_to_2d(self.edge.start)
        self.assertEqual(point2d_1, volmdlr.Point2D(0, 0))
        point2d_2 = self.surface.point3d_to_2d(self.edge.end)
        self.assertEqual(point2d_2, volmdlr.Point2D(1.0, 0))

    def test_rectangular_cut(self):
        face = vmf.ExtrusionFace3D.from_surface_rectangular_cut(self.surface, 0, 1, 0, 2)
        self.assertEqual(face.surface2d.area(), 2)

    def test_from_step(self):
        step = vms.Step.from_file("surfaces/objects_extrusion_tests/bspline_extruded_simple.step")
        model = step.to_volume_model()
        extrusion_surface = model.primitives[0].primitives[0].surface3d
        self.assertEqual(extrusion_surface.direction, -volmdlr.Z3D)
        self.assertEqual(extrusion_surface.edge.degree, 3)
        self.assertEqual(extrusion_surface.edge.knot_multiplicities, [4, 1, 4])

    def test_linesegment2d_to_3d(self):
        surface = surfaces.ExtrusionSurface3D.load_from_file(
            "surfaces/objects_extrusion_tests/extrusion_surface_undefined_direction_linesegment.json")
        point1 = volmdlr.Point2D(0.9020984833336293, -0.08534036750789999)
        point2 = volmdlr.Point2D(0.9286913444016728, -0.07799341694)
        linesegment2d = vme.LineSegment2D(point1, point2)
        start3d = surface.point2d_to_3d(point1)
        end3d = surface.point2d_to_3d(point2)
        result = surface.linesegment2d_to_3d(linesegment2d)[0]
        self.assertTrue(result.start.is_close(start3d))
        self.assertTrue(result.end.is_close(end3d))
        surface = surfaces.ExtrusionSurface3D.load_from_file(
            "surfaces/objects_extrusion_tests/extrusion_surface_test_linesegment2d_to_3d.json")
        linesegment2d = vme.LineSegment2D.load_from_file(
            "surfaces/objects_extrusion_tests/linesegment2d_to_linesegment3d.json")
        result = surface.linesegment2d_to_3d(linesegment2d)[0]
        self.assertIsInstance(result, vme.LineSegment3D)

        surface = surfaces.ExtrusionSurface3D.load_from_file(
            "surfaces/objects_extrusion_tests/extrusion_surface_test_linesegment2d_to_3d.json")
        linesegment2d = vme.LineSegment2D.load_from_file(
            "surfaces/objects_extrusion_tests/linesegment2d_to_self_edge.json")
        result = surface.linesegment2d_to_3d(linesegment2d)[0]
        self.assertEqual(result, surface.edge)

    def test_arc3d_to_2d(self):
        surface = surfaces.ExtrusionSurface3D.load_from_file(
            "surfaces/objects_extrusion_tests/extrusion_surface_test_arc3d_to_2d.json")
        contour3d = vmw.Contour3D.load_from_file(
            "surfaces/objects_extrusion_tests/extrusion_contour_test_arc3d_to_2d.json")
        arc3d = contour3d.primitives[2]
        result = surface.arc3d_to_2d(arc3d)[0]
        self.assertTrue(result.start.is_close(volmdlr.Point2D(1.0, 0.0032000000499998738)))
        self.assertTrue(result.end.is_close(volmdlr.Point2D(0.13555464614559587, 0.0032000000499998738)))

    def test_frame_mapping(self):
        surface = self.surface
        new_frame = volmdlr.Frame3D(volmdlr.Point3D(0, 0, 1), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        new_surface = surface.frame_mapping(new_frame, "old")
        self.assertEqual(new_surface.edge.start.z, 1)
        self.assertTrue(new_surface.frame.origin.is_close(volmdlr.Point3D(-0.025917292, 0.002544355, 1.0)))



if __name__ == '__main__':
    unittest.main()
