import unittest
from volmdlr.shells import OpenTriangleShell3D, ClosedTriangleShell3D, DisplayTriangleShell3D
from volmdlr.primitives3d import Block
from volmdlr import OXYZ

SHOW_BABYLONJS = True


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
        display_triangle_shell = OpenTriangleShell3D(self.faces).to_display_triangle_shell()
        self.assertEqual(12, len(display_triangle_shell.indices))
        self.assertEqual(8, len(display_triangle_shell.positions))

        self.assertEqual(display_triangle_shell, DisplayTriangleShell3D.dict_to_object(display_triangle_shell.to_dict()))

        if SHOW_BABYLONJS:
            display_triangle_shell.babylonjs()
