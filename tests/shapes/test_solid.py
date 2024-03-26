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
        solid = shapes.Solid.make_extrusion_from_frame_and_wires(volmdlr.OXYZ, outer_contour2d, inner_contours2d, height)
        self.assertAlmostEqual(solid.volume(), (length * width - math.pi * (radius ** 2)) * height)


if __name__ == '__main__':
    unittest.main()
