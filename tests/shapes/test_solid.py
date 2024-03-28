import math
import unittest

import volmdlr
from volmdlr import shapes, wires, curves


class TestSolid(unittest.TestCase):

    def test_make_extrusion(self):
        length, width, height, radius = 0.4, 0.3, 0.08, 0.1
        outer_contour2d = wires.Contour2D.rectangle_from_center_and_sides(volmdlr.O2D, x_length=length, y_length=width,
                                                                          is_trigo=True)
        inner_contours2d = [wires.Contour2D.from_circle(
            circle=curves.Circle2D.from_center_and_radius(volmdlr.O2D, radius, is_trigo=False))]
        solid = shapes.Solid.make_extrusion_from_frame_and_wires(volmdlr.OXYZ, outer_contour2d,
                                                                 inner_contours2d, height)
        self.assertAlmostEqual(solid.volume(), (length * width - math.pi * (radius ** 2)) * height)

    def test_make_wedge(self):
        dx, dy, dz = 1, 2, 1
        solid = shapes.Solid.make_wedge(dx=dx, dy=dy, dz=dz, xmin=dx / 2, xmax=dx / 2, zmin=dz / 2, zmax=dz / 2,
                                        local_frame_origin=volmdlr.Point3D(-0.5, 0.5, 0.0),
                                        local_frame_direction=-volmdlr.Y3D,
                                        local_frame_x_direction=volmdlr.X3D)

        self.assertAlmostEqual(solid.volume(), (1 / 3) * dy)

        solid = shapes.Solid.make_wedge(dx=dx, dy=dy, dz=dz, xmin=dx / 4, xmax=3 * dx / 4,
                                        zmin=dz / 4, zmax=3 * dz / 4,
                                        local_frame_origin=volmdlr.Point3D(-0.5, 0.5, 0.0),
                                        local_frame_direction=-volmdlr.Y3D,
                                        local_frame_x_direction=volmdlr.X3D)

        self.assertAlmostEqual(solid.volume(), (1 / 3) * dy * (1 + 0.5**2 + 0.5))


if __name__ == '__main__':
    unittest.main()
