import unittest

import numpy as np

from volmdlr.display import Mesh3D

SHOW_BABYLONJS = True


class TestMesh3D(unittest.TestCase):

    def setUp(self):
        # Sample data for testing
        self.positions1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.indices1 = np.array([[0, 1, 2]])
        self.shell1 = Mesh3D(self.positions1, self.indices1, "Shell1")

        self.positions2 = np.array([[0, 1, 0], [1, 1, 0], [1, 0, 0]])  # Note shared vertex with shell1
        self.indices2 = np.array([[0, 1, 2]])
        self.shell2 = Mesh3D(self.positions2, self.indices2, "Shell2")

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
        self.shell3 = Mesh3D(self.positions3, self.indices3, "Shell3")

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
        self.shell4 = Mesh3D(self.positions4, self.indices4, "Shell4")

    def test_concatenate(self):
        concatenated_shell = self.shell1.concatenate(self.shell2)

        expected_positions = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])

        np.testing.assert_array_equal(np.sort(concatenated_shell.positions, axis=0),
                                      np.sort(expected_positions, axis=0))
        # Compare indices carefully since their correctness depends on the order of positions

    def test_add_operator(self):
        combined_shell = self.shell1 + self.shell2

        expected_positions = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])

        np.testing.assert_array_equal(np.sort(combined_shell.vertices, axis=0),
                                      np.sort(expected_positions, axis=0))
        # Compare indices carefully since their correctness depends on the order of positions

    def test_concatenate_cube(self):
        combined_shell = self.shell3 + self.shell4

        self.assertEqual(22, len(combined_shell.triangles))
        self.assertEqual(12, len(combined_shell.vertices))

    def test_equality(self):
        self.assertNotEqual(self.shell1, self.shell2)
        self.assertFalse(self.shell1._data_eq(self.shell2))

    def test_concatenate_empty(self):
        empty_shell = Mesh3D(np.array([]), np.array([]))

        self.assertEqual(self.shell1, self.shell1 + empty_shell)
        self.assertEqual(self.shell1, empty_shell + self.shell1)
        self.assertNotEqual(self.shell1, self.shell1 + self.shell2)


if __name__ == '__main__':
    unittest.main()
