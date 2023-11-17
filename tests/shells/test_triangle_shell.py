import os
import unittest
from dessia_common.core import DessiaObject
from volmdlr.shells import (
    OpenTriangleShell3D,
    ClosedTriangleShell3D,
    DisplayTriangleShell3D,
)
from volmdlr.primitives3d import Block
from volmdlr import OXYZ

SHOW_BABYLONJS = True


folder = os.path.dirname(os.path.realpath(__file__))


class TestTrianglShell3D(unittest.TestCase):
    def setUp(self) -> None:
        self.faces = Block(OXYZ).triangulation().faces

    def test_open_triangle_shell(self):
        open_triangle_shell = OpenTriangleShell3D(self.faces)
        self.assertEqual(12, len(open_triangle_shell.faces))

        if SHOW_BABYLONJS:
            open_triangle_shell.babylonjs()

    def test_closed_triangle_shell(self):
        closed_triangle_shell = ClosedTriangleShell3D(self.faces)
        self.assertEqual(12, len(closed_triangle_shell.faces))
        self.assertEqual(1.0, closed_triangle_shell.volume())

        if SHOW_BABYLONJS:
            closed_triangle_shell.babylonjs()

    def test_display_triangle_shell(self):
        display_triangle_shell = OpenTriangleShell3D(
            self.faces
        ).to_display_triangle_shell()
        self.assertEqual(12, len(display_triangle_shell.indices))
        self.assertEqual(8, len(display_triangle_shell.positions))

        self.assertEqual(
            display_triangle_shell,
            DisplayTriangleShell3D.dict_to_object(display_triangle_shell.to_dict()),
        )

        if SHOW_BABYLONJS:
            display_triangle_shell.babylonjs()

    def test_turn_normals_outwards(self):
        closed_shell = DessiaObject.load_from_file(
            os.path.join(folder, "closedtriangleshell3d.json")
        )
        self.assertFalse(closed_shell.are_normals_pointing_outwards())
        new_closed_shell = closed_shell.turn_normals_outwards()
        self.assertTrue(new_closed_shell.are_normals_pointing_outwards())

    def test_turn_normals_inwards(self):
        closed_shell = DessiaObject.load_from_file(
            os.path.join(folder, "closedtriangleshell3d.json")
        )
        self.assertFalse(closed_shell.are_normals_pointing_inwards())
        new_closed_shell = closed_shell.turn_normals_inwards()
        self.assertTrue(new_closed_shell.are_normals_pointing_inwards())

    def test_closedtriagleshell3d_subtraction(self):
        shell1 = DessiaObject.load_from_file(os.path.join(folder, "shell1(1).json"))
        shell2 = DessiaObject.load_from_file(os.path.join(folder, "shell2(1).json"))
        new_shell = shell2.subtract_to_closed_shell(shell1)[0]
        self.assertEqual(len(new_shell.faces), 76)


if __name__ == "__main__":
    unittest.main()
