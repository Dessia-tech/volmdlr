import unittest

import numpy as np

from volmdlr.display import Mesh3D

SHOW_BABYLONJS = True


class TestMesh3D(unittest.TestCase):

    def setUp(self):
        # Sample data for testing
        self.positions1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.indices1 = np.array([[0, 1, 2]])
        self.mesh1 = Mesh3D(self.positions1, self.indices1, "Mesh1")

        self.positions2 = np.array([[0, 1, 0], [1, 1, 0], [1, 0, 0]])  # Note shared vertex with mesh1
        self.indices2 = np.array([[0, 1, 2]])
        self.mesh2 = Mesh3D(self.positions2, self.indices2, "Mesh2")

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
        self.mesh3 = Mesh3D(self.positions3, self.indices3, "Mesh3")

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
        self.mesh4 = Mesh3D(self.positions4, self.indices4, "Mesh4")

    def test_concatenate(self):
        concatenated_mesh = self.mesh1.concatenate(self.mesh2)

        expected_positions = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])

        np.testing.assert_array_equal(np.sort(concatenated_mesh.positions, axis=0),
                                      np.sort(expected_positions, axis=0))
        # Compare indices carefully since their correctness depends on the order of positions

    def test_add_operator(self):
        combined_mesh = self.mesh1 + self.mesh2

        expected_positions = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])

        np.testing.assert_array_equal(np.sort(combined_mesh.vertices, axis=0),
                                      np.sort(expected_positions, axis=0))
        # Compare indices carefully since their correctness depends on the order of positions

    def test_concatenate_cube(self):
        combined_mesh = self.mesh3 + self.mesh4

        self.assertEqual(22, len(combined_mesh.triangles))
        self.assertEqual(12, len(combined_mesh.vertices))

    def test_merge_cube(self):
        merged_mesh_1 = self.mesh3.merge(self.mesh4, mutualize=False)

        self.assertEqual(24, len(merged_mesh_1.triangles))
        self.assertEqual(16, len(merged_mesh_1.vertices))

        merged_mesh_2 = self.mesh3.merge(self.mesh4, mutualize=True)

        self.assertEqual(22, len(merged_mesh_2.triangles))
        self.assertEqual(12, len(merged_mesh_2.vertices))

    def test_equality(self):
        self.assertNotEqual(self.mesh1, self.mesh2)
        self.assertEqual(self.mesh1, self.mesh1)
        self.assertFalse(self.mesh1._data_eq(self.mesh2))
        self.assertTrue(self.mesh1._data_eq(self.mesh1))

    def test_concatenate_empty(self):
        empty_mesh = Mesh3D(np.array([]), np.array([]))

        self.assertEqual(self.mesh1, self.mesh1 + empty_mesh)
        self.assertEqual(self.mesh1, empty_mesh + self.mesh1)
        self.assertNotEqual(self.mesh1, self.mesh1 + self.mesh2)


if __name__ == '__main__':
    unittest.main()
