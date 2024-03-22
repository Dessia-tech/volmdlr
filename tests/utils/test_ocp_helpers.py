import unittest
import math
from OCP.BRepBuilderAPI import BRepBuilderAPI_MakeFace
from OCP.gp import gp_Ax3, gp_Cylinder
from volmdlr.utils import ocp_helpers


class TestOCPHelpers(unittest.TestCase):
    def test_plot_face(self):
        frame = gp_Ax3()
        cylinder = gp_Cylinder(frame, 8)
        cylindrical_face = BRepBuilderAPI_MakeFace(cylinder, 0, math.pi, 0, 15).Face()
        ax = ocp_helpers.plot_face(cylindrical_face)
        self.assertEqual(len(ax.lines), 4)


if __name__ == '__main__':
    unittest.main()
