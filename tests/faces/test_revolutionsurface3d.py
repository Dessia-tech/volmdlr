import unittest
import math

import volmdlr
import volmdlr.edges as vme
import volmdlr.wires as vmw
import volmdlr.faces


class TestRevolutionSurface3D(unittest.TestCase):
    linesegment = vme.LineSegment3D(volmdlr.Point3D(0.5, 0, 0), volmdlr.Point3D(0.5, 0, 0.5))
    arc = vme.Arc3D(
        volmdlr.Point3D(0.5, 0, 0.5),
        volmdlr.Point3D(0.3 + 0.2 * math.cos(math.pi / 6), 0, 0.5 + 0.2 * math.sin(math.pi / 6)),
        volmdlr.Point3D(0.3 + 0.2 * math.cos(math.pi / 3), 0, 0.5 + 0.2 * math.sin(math.pi / 3)),
    )

    wire = vmw.Wire3D([linesegment, arc])
    axis_point = volmdlr.O3D
    axis = volmdlr.Z3D

    def test_init(self):

        surface = volmdlr.faces.RevolutionSurface3D(self.wire, self.axis_point, self.axis)

        self.assertEqual(surface.x_periodicity, volmdlr.TWO_PI)
        self.assertEqual(surface.y_periodicity, None)
        self.assertEqual(surface.frame.origin, self.axis_point)
        self.assertEqual(surface.axis_point, self.axis_point)
        self.assertEqual(surface.axis, self.axis)

    def test_point2d_to_3d(self):
        surface = volmdlr.faces.RevolutionSurface3D(self.wire, self.axis_point, self.axis)

        point2d = volmdlr.Point2D(math.pi, 0.7047817224492219)
        point3d = surface.point2d_to_3d(point2d)
        expected_point3d = volmdlr.Point3D(-0.5, 0, 0.5)

        self.assertEqual(point3d, expected_point3d)

    def test_point3d_to_2d(self):
        surface = volmdlr.faces.RevolutionSurface3D(self.wire, self.axis_point, self.axis)

        point3d = volmdlr.Point3D(-0.5, 0, 0.5)
        point2d = surface.point3d_to_2d(point3d)
        expected_point2d = volmdlr.Point2D(math.pi, 0.7047817224492219)

        self.assertEqual(point2d, expected_point2d)

    def test_rectangular_cut(self):
        surface = volmdlr.faces.RevolutionSurface3D(wire=self.wire, axis_point=self.axis_point, axis=self.axis)
        rectangular_cut = surface.rectangular_cut(0, volmdlr.TWO_PI, 0, 1)
        self.assertEqual(rectangular_cut.surface2d.area(), volmdlr.TWO_PI)


if __name__ == "__main__":
    unittest.main()
