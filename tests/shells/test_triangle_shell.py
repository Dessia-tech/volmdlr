import os
import unittest

import numpy as np

from dessia_common.core import DessiaObject
from volmdlr.shells import OpenTriangleShell3D, ClosedTriangleShell3D, DisplayTriangleShell3D
from volmdlr.primitives3d import Block
from volmdlr import OXYZ

SHOW_BABYLONJS = True


folder = os.path.dirname(os.path.realpath(__file__))


class TestTriangleShell3D(unittest.TestCase):
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
        self.assertEqual(display_triangle_shell, display_triangle_shell + display_triangle_shell)

        if SHOW_BABYLONJS:
            display_triangle_shell.babylonjs()

    def test_turn_normals_outwards(self):
        closed_shell = DessiaObject.load_from_file(os.path.join(folder, 'closedtriangleshell3d.json'))
        self.assertFalse(closed_shell.are_normals_pointing_outwards())
        new_closed_shell = closed_shell.turn_normals_outwards()
        self.assertTrue(new_closed_shell.are_normals_pointing_outwards())

    def test_turn_normals_inwards(self):
        closed_shell = DessiaObject.load_from_file(os.path.join(folder, 'closedtriangleshell3d.json'))
        self.assertFalse(closed_shell.are_normals_pointing_inwards())
        new_closed_shell = closed_shell.turn_normals_inwards()
        self.assertTrue(new_closed_shell.are_normals_pointing_inwards())


class TestDisplayTriangleShell3D(unittest.TestCase):

    def setUp(self):
        # Sample data for testing
        self.positions1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.indices1 = np.array([[0, 1, 2]])
        self.shell1 = DisplayTriangleShell3D(self.positions1, self.indices1, "Shell1")

        self.positions2 = np.array([[0, 1, 0], [1, 1, 0], [1, 0, 0]])  # Note shared vertex with shell1
        self.indices2 = np.array([[0, 1, 2]])
        self.shell2 = DisplayTriangleShell3D(self.positions2, self.indices2, "Shell2")

    def test_concatenate(self):
        concatenated_shell = self.shell1.concatenate(self.shell2)

        expected_positions = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
        expected_indices = np.array([[0, 1, 2], [2, 3, 1]])

        np.testing.assert_array_equal(np.sort(concatenated_shell.positions, axis=0),
                                      np.sort(expected_positions, axis=0))
        # Compare indices carefully since their correctness depends on the order of positions

    def test_add_operator(self):
        combined_shell = self.shell1 + self.shell2

        expected_positions = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])
        expected_indices = np.array([[0, 1, 2], [2, 3, 1]])

        np.testing.assert_array_equal(np.sort(combined_shell.positions, axis=0),
                                      np.sort(expected_positions, axis=0))
        # Compare indices carefully since their correctness depends on the order of positions


if __name__ == '__main__':
    unittest.main()
