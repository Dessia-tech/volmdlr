import unittest
import math

import volmdlr
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.faces
from volmdlr import surfaces


class TestRevolutionSurface3D(unittest.TestCase):
    linesegment = vme.LineSegment3D(volmdlr.Point3D(0.5, 0, 0), volmdlr.Point3D(0.5, 0, 0.5))
    arc = vme.Arc3D(volmdlr.Point3D(0.5, 0, 0.5),
                    volmdlr.Point3D(0.3 + 0.2 * math.cos(math.pi / 6), 0, 0.5 + 0.2 * math.sin(math.pi / 6)),
                    volmdlr.Point3D(0.3 + 0.2 * math.cos(math.pi / 3), 0, 0.5 + 0.2 * math.sin(math.pi / 3)))

    wire = vmw.Wire3D([linesegment, arc])
    axis_point = volmdlr.O3D
    axis = volmdlr.Z3D
    surface = surfaces.RevolutionSurface3D(wire, axis_point, axis)

    def test_point2d_to_3d(self):
        surface = surfaces.RevolutionSurface3D(self.wire, self.axis_point, self.axis)

        point2d = volmdlr.Point2D(math.pi, 0.5)
        point3d = surface.point2d_to_3d(point2d)
        expected_point3d = volmdlr.Point3D(-0.5, 0, 0.5)

        self.assertTrue(point3d.is_close(expected_point3d))

    def test_point3d_to_2d(self):
        surface = surfaces.RevolutionSurface3D(self.wire, self.axis_point, self.axis)

        point3d = volmdlr.Point3D(-0.5, 0, 0.5)
        point2d = surface.point3d_to_2d(point3d)
        expected_point2d = volmdlr.Point2D(math.pi, 0.5)

        self.assertTrue(point2d.is_close(expected_point2d))

    def test_rectangular_cut(self):
        surface = surfaces.RevolutionSurface3D(wire=self.wire, axis_point=self.axis_point, axis=self.axis)
        rectangular_cut = volmdlr.faces.RevolutionFace3D.from_surface_rectangular_cut(
            surface, 0, volmdlr.TWO_PI, 0, 1)
        self.assertEqual(rectangular_cut.surface2d.area(), volmdlr.TWO_PI)

    def test_frame_mapping(self):
        surface = self.surface
        new_frame = volmdlr.Frame3D(volmdlr.Point3D(0, 1, 0), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        new_surface = surface.frame_mapping(new_frame, "old")
        self.assertEqual(new_surface.wire.primitives[0].start.y, 1)
        self.assertTrue(new_surface.frame.origin.is_close(volmdlr.Point3D(0, 1, 0)))


if __name__ == '__main__':
    unittest.main()
