import math
import unittest

import volmdlr
from volmdlr import shapes, wires, curves
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp


class TestShell(unittest.TestCase):

    def test_make_extrusion(self):
        length, width, height = 0.4, 0.3, 0.08
        contour = wires.Contour2D.rectangle_from_center_and_sides(volmdlr.O2D, x_length=length, y_length=width,
                                                                  is_trigo=True).to_3d(volmdlr.O3D, volmdlr.X3D,
                                                                                       volmdlr.Y3D)
        shell = shapes.Shell.make_extrusion(contour, extrusion_direction=volmdlr.Z3D, extrusion_length=height)
        gprops = GProp_GProps()
        BRepGProp.VolumeProperties_s(shell.wrapped, gprops)
        volume = gprops.Mass()
        self.assertAlmostEqual(volume, length * width * height)


if __name__ == '__main__':
    unittest.main()
