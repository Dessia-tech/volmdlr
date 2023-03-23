import unittest
import volmdlr
import volmdlr.edges as vme
import volmdlr.faces as vmf
import volmdlr.step as vms


class TestExtrusionSurface3D(unittest.TestCase):
    control_points = [
        volmdlr.Point3D(-0.025917292, 0.002544355, 0.0),
        volmdlr.Point3D(-0.005449685, -0.007265721, 0.0),
        volmdlr.Point3D(0.0, 0.0, 0.0),
        volmdlr.Point3D(0.014457705000000001, -0.002636091, 0.0),
        volmdlr.Point3D(0.013503079, -0.014007147, 0.0)]
    edge = vme.BSplineCurve3D(3, control_points, [4, 1, 4], [0.0, 0.5, 1.0])
    surface = vmf.ExtrusionSurface3D(edge, -volmdlr.Z3D)

    def test_point2d_to_3d(self):
        point3d = self.surface.point2d_to_3d(volmdlr.Point2D(0.5, 0.5))
        self.assertTrue(point3d.is_close(volmdlr.Point3D(0.002252005, -0.002475453, -0.5)))

    def test_point3d_to_2d(self):
        point2d_1 = self.surface.point3d_to_2d(self.edge.start)
        self.assertEqual(point2d_1, volmdlr.Point2D(0, 0))
        point2d_2 = self.surface.point3d_to_2d(self.edge.end)
        self.assertEqual(point2d_2, volmdlr.Point2D(1.0, 0))

    def test_rectangular_cut(self):
        face = self.surface.rectangular_cut(0, 1, 0, 2)
        self.assertEqual(face.surface2d.area(), 2)

    def test_from_step(self):
        step = vms.Step.from_file("faces/objects_extrusion_tests/bspline_extruded_simple.step")
        model = step.to_volume_model()
        extrusion_surface = model.primitives[0].primitives[0].surface3d
        self.assertEqual(extrusion_surface.direction, -volmdlr.Z3D)
        self.assertEqual(extrusion_surface.edge.degree, 3)
        self.assertEqual(extrusion_surface.edge.knot_multiplicities, [4, 1, 4])


if __name__ == '__main__':
    unittest.main()
