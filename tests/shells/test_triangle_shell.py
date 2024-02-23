import os
import unittest

import numpy as np

from dessia_common.core import DessiaObject
from volmdlr.shells import (
    OpenTriangleShell3D,
    ClosedTriangleShell3D,
    DisplayTriangleShell3D,
)
from volmdlr.primitives3d import Block
from volmdlr import OXYZ

SHOW_BABYLONJS = False


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
        display_triangle_shell = OpenTriangleShell3D(
            self.faces
        ).to_display_triangle_shell()
        self.assertEqual(12, len(display_triangle_shell.indices))
        self.assertEqual(8, len(display_triangle_shell.positions))
        self.assertEqual(display_triangle_shell, DisplayTriangleShell3D.dict_to_object(display_triangle_shell.to_dict()))
        self.assertEqual(len(display_triangle_shell.indices), len((display_triangle_shell + display_triangle_shell).indices))

        if SHOW_BABYLONJS:
            display_triangle_shell.babylonjs()

    def test_turn_normals_outwards(self):
        closed_shell = DessiaObject.from_json(
            os.path.join(folder, "closedtriangleshell3d.json")
        )
        self.assertFalse(closed_shell.are_normals_pointing_outwards())
        new_closed_shell = closed_shell.turn_normals_outwards()
        self.assertTrue(new_closed_shell.are_normals_pointing_outwards())

    def test_turn_normals_inwards(self):
        closed_shell = DessiaObject.from_json(
            os.path.join(folder, "closedtriangleshell3d.json")
        )
        self.assertFalse(closed_shell.are_normals_pointing_inwards())
        new_closed_shell = closed_shell.turn_normals_inwards()
        self.assertTrue(new_closed_shell.are_normals_pointing_inwards())

    def test_closedtriagleshell3d_subtraction(self):
        shell1 = DessiaObject.from_json(os.path.join(folder, "shell1(1).json"))
        shell2 = DessiaObject.from_json(os.path.join(folder, "shell2(1).json"))
        new_shell = shell2.subtract_to_closed_shell(shell1)[0]
        self.assertEqual(len(new_shell.faces), 76)


class TestDisplayTriangleShell3D(unittest.TestCase):

    def setUp(self):
        # Sample data for testing
        self.positions1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.indices1 = np.array([[0, 1, 2]])
        self.shell1 = DisplayTriangleShell3D(self.positions1, self.indices1, "Shell1")

        self.positions2 = np.array([[0, 1, 0], [1, 1, 0], [1, 0, 0]])  # Note shared vertex with shell1
        self.indices2 = np.array([[0, 1, 2]])
        self.shell2 = DisplayTriangleShell3D(self.positions2, self.indices2, "Shell2")

        self.positions3 = np.array(
            [
                [0.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                [0.0, 1.0, 0.0],
                [0.0, 1.0, 1.0],
                [1.0, 0.0, 0.0],
                [1.0, 0.0, 1.0],
                [1.0, 1.0, 0.0],
                [1.0, 1.0, 1.0],
            ]
        )
        self.indices3 = np.array(
            [
                [2, 6, 7],
                [0, 4, 5],
                [1, 7, 5],
                [0, 2, 6],
                [4, 6, 7],
                [1, 3, 7],
                [0, 2, 3],
                [2, 7, 3],
                [0, 6, 4],
                [4, 7, 5],
                [0, 5, 1],
                [0, 3, 1],
            ]
        )
        self.shell3 = DisplayTriangleShell3D(self.positions3, self.indices3, "Shell3")

        self.positions4 = np.array(
            [
                [0.0, 0.0, 1.0],
                [0.0, 0.0, 2.0],
                [0.0, 1.0, 1.0],
                [0.0, 1.0, 2.0],
                [1.0, 0.0, 1.0],
                [1.0, 0.0, 2.0],
                [1.0, 1.0, 1.0],
                [1.0, 1.0, 2.0],
            ]
        )
        self.indices4 = np.array(
            [
                [2, 7, 3],
                [1, 7, 5],
                [0, 6, 4],
                [4, 7, 5],
                [0, 3, 1],
                [0, 2, 6],
                [4, 6, 7],
                [2, 6, 7],
                [0, 4, 5],
                [1, 3, 7],
                [0, 2, 3],
                [0, 5, 1],
            ]
        )
        self.shell4 = DisplayTriangleShell3D(self.positions4, self.indices4, "Shell4")

    def test_concatenate(self):
        concatenated_shell = self.shell1.concatenate(self.shell2)

        expected_positions = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])

        np.testing.assert_array_equal(np.sort(concatenated_shell.positions, axis=0),
                                      np.sort(expected_positions, axis=0))
        # Compare indices carefully since their correctness depends on the order of positions

    def test_add_operator(self):
        combined_shell = self.shell1 + self.shell2

        expected_positions = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])

        np.testing.assert_array_equal(np.sort(combined_shell.positions, axis=0),
                                      np.sort(expected_positions, axis=0))
        # Compare indices carefully since their correctness depends on the order of positions

    def test_concatenate_cube(self):
        combined_shell = self.shell3 + self.shell4

        self.assertEqual(22, len(combined_shell.indices))
        self.assertEqual(12, len(combined_shell.positions))

    def test_equality(self):
        self.assertNotEqual(self.shell1, self.shell2)
        self.assertFalse(self.shell1._data_eq(self.shell2))

    def test_concatenate_empty(self):
        empty_shell = DisplayTriangleShell3D(np.array([]), np.array([]))

        self.assertEqual(self.shell1, self.shell1 + empty_shell)
        self.assertEqual(self.shell1, empty_shell + self.shell1)
        self.assertNotEqual(self.shell1, self.shell1 + self.shell2)


if __name__ == '__main__':
    unittest.main()
