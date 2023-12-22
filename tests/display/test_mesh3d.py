import unittest

import numpy as np

from volmdlr.display import Mesh3D

SHOW_BABYLONJS = True


class TestMesh3D(unittest.TestCase):
    def setUp(self):
        # Sample data for testing
        self.vertices1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0]])
        self.triangles1 = np.array([[0, 1, 2]])
        self.mesh1 = Mesh3D(self.vertices1, self.triangles1, "Mesh1")

        self.vertices2 = np.array([[0, 1, 0], [1, 1, 0], [1, 0, 0]])  # Note shared vertex with mesh1
        self.triangles2 = np.array([[0, 1, 2]])
        self.mesh2 = Mesh3D(self.vertices2, self.triangles2, "Mesh2")

        self.vertices3 = np.array(
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
        self.triangles3 = np.array(
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
        self.mesh3 = Mesh3D(self.vertices3, self.triangles3, "Mesh3")

        self.vertices4 = np.array(
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
        self.triangles4 = np.array(
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
        self.mesh4 = Mesh3D(self.vertices4, self.triangles4, "Mesh4")

    def test_merge_without_mutualization(self):
        merged_meshes = self.mesh1.merge(self.mesh2, False, False)
        expected_vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0], [0, 1, 0], [1, 1, 0]])

        np.testing.assert_array_equal(np.sort(merged_meshes.vertices, axis=0), np.sort(expected_vertices, axis=0))

    def test_merge_with_mutualization(self):
        merged_meshes = self.mesh1.merge(self.mesh2, True, True)
        expected_vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])

        np.testing.assert_array_equal(np.sort(merged_meshes.vertices, axis=0), np.sort(expected_vertices, axis=0))

    def test_add_operator(self):
        # No mutualization
        added_mesh = self.mesh1 + self.mesh2
        expected_vertices = np.array([[0, 0, 0], [1, 0, 0], [1, 0, 0], [0, 1, 0], [0, 1, 0], [1, 1, 0]])

        np.testing.assert_array_equal(np.sort(added_mesh.vertices, axis=0), np.sort(expected_vertices, axis=0))

    def test_union_operator(self):
        # Mutualization
        added_mesh = self.mesh1 | self.mesh2
        expected_vertices = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [1, 1, 0]])

        np.testing.assert_array_equal(np.sort(added_mesh.vertices, axis=0), np.sort(expected_vertices, axis=0))

    def test_merge_cube_without_mutualization(self):
        merged_mesh_1 = self.mesh3.merge(self.mesh4, mutualize_vertices=False, mutualize_triangles=False)

        self.assertEqual(24, len(merged_mesh_1.triangles))
        self.assertEqual(16, len(merged_mesh_1.vertices))

    def test_merge_cube_with_mutualization(self):
        merged_mesh_2 = self.mesh3.merge(self.mesh4, mutualize_vertices=True, mutualize_triangles=True)

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


if __name__ == "__main__":
    unittest.main()
