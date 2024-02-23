"""
Unit testing of volmdlr.display.Mesh3D class.
"""
import math
import os
import tempfile
import unittest

import numpy as np
import trimesh
from dessia_common.serialization import BinaryFile
from volmdlr import Point3D
from volmdlr.display import Mesh3D
from volmdlr.faces import Triangle3D
from volmdlr.shells import ClosedTriangleShell3D, OpenTriangleShell3D

SHOW_BABYLONJS = False

FOLDER = os.path.dirname(os.path.realpath(__file__))


class TestMesh3D(unittest.TestCase):
    def setUp(self):
        # Sample data for testing
        self.vertices1 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        self.triangles1 = np.array([[0, 1, 2]])
        self.mesh1 = Mesh3D(self.vertices1, self.triangles1, name="Mesh1")

        self.vertices2 = np.array([[0, 1, 0], [1, 1, 0], [1, 0, 0]])  # Note shared vertex with mesh1
        self.triangles2 = np.array([[0, 1, 2]])
        self.mesh2 = Mesh3D(self.vertices2, self.triangles2, name="Mesh2")

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
        self.mesh3 = Mesh3D(self.vertices3, self.triangles3, name="Mesh3")

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
        self.mesh4 = Mesh3D(self.vertices4, self.triangles4, name="Mesh4")

    def test_resize(self):
        exepected_vertices = np.array([[0.0, 0.0, 0.0], [0.001, 0.0, 0.0], [0.0, 0.001, 0.0]])
        expected_mesh = Mesh3D(exepected_vertices, self.triangles1)

        self.assertEqual(expected_mesh, self.mesh1.resize(0.001))

    def test_round_vertices(self):
        self.assertEqual(self.mesh4, self.mesh4.round_vertices())

        vertices_1_bis = np.array([[0.0, 0.0, 0.0], [1.000000000001, 0.0, 0.0], [0.0, 0.99999999999999, 0.0]])
        mesh_1_bis = Mesh3D(vertices_1_bis, self.triangles1)

        self.assertNotEqual(self.mesh1, mesh_1_bis)
        self.assertEqual(self.mesh1, mesh_1_bis.round_vertices(6))

    def test_remove_degenerated_triangles(self):
        triangles = np.array([[0, 1, 2], [0, 1, 1]])
        degenerated_mesh = Mesh3D(self.vertices1, triangles)

        self.assertNotEqual(self.mesh1, degenerated_mesh)
        self.assertEqual(self.mesh1, degenerated_mesh.remove_degenerate_triangles())

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
        merged_mesh_1 = self.mesh3.merge(self.mesh4, merge_vertices=False, merge_triangles=False)

        self.assertEqual(24, len(merged_mesh_1.triangles))
        self.assertEqual(16, len(merged_mesh_1.vertices))

    def test_merge_cube_with_mutualization(self):
        merged_mesh_2 = self.mesh3.merge(self.mesh4, merge_vertices=True, merge_triangles=True)

        self.assertEqual(22, len(merged_mesh_2.triangles))
        self.assertEqual(12, len(merged_mesh_2.vertices))

    def test_split_shared_vertices(self):
        merged_mesh_2 = self.mesh3.merge(self.mesh4, merge_vertices=True, merge_triangles=True)

        self.assertEqual(3 * len(merged_mesh_2.triangles), len(merged_mesh_2.split_shared_vertices().vertices))

        self.assertEqual(3 * len(self.mesh1.triangles), len(self.mesh1.split_shared_vertices().vertices))
        self.assertEqual(3 * len(self.mesh2.triangles), len(self.mesh2.split_shared_vertices().vertices))
        self.assertEqual(3 * len(self.mesh3.triangles), len(self.mesh3.split_shared_vertices().vertices))
        self.assertEqual(3 * len(self.mesh4.triangles), len(self.mesh4.split_shared_vertices().vertices))

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

    def test_serialization(self):
        dict_ = self.mesh4.to_dict()
        mesh = Mesh3D.dict_to_object(dict_)

        self.assertEqual(self.mesh4, mesh)

    def test_to_triangles3d(self):
        expected_triangles3d_1 = [Triangle3D(Point3D(0, 0, 0), Point3D(1, 0, 0), Point3D(0, 1, 0))]
        self.assertEqual(expected_triangles3d_1, self.mesh1.to_triangles3d())

        expected_triangles3d_4 = []
        for triangle in self.triangles4:
            expected_triangles3d_4.append(
                Triangle3D(
                    Point3D(*self.vertices4[triangle[0]]),
                    Point3D(*self.vertices4[triangle[1]]),
                    Point3D(*self.vertices4[triangle[2]]),
                )
            )
        self.assertEqual(expected_triangles3d_4, self.mesh4.to_triangles3d())

    def test_check_concistency(self):
        self.assertTrue(self.mesh1.check_consistency())
        self.assertTrue(self.mesh2.check_consistency())
        self.assertTrue(self.mesh3.check_consistency())
        self.assertTrue(self.mesh4.check_consistency())

        unconsistent_mesh = Mesh3D(self.vertices1, np.array([[0, 1, 3]]))
        self.assertFalse(unconsistent_mesh.check_consistency())

    def test_area(self):
        self.assertEqual(0.5, self.mesh1.area())

    def test_get_edges_triangles(self):
        # Test the get_edges_triangles method
        expected_edges = np.array([[[0, 1], [0, 2], [1, 2]]])
        np.testing.assert_array_equal(expected_edges, self.mesh1.get_edges_triangles())

    def test_compute_len_edges(self):
        # Test the compute_len_edges method
        expected_lengths = np.array([[1.0, 1.0, math.sqrt(2)]])

        lengths, _ = self.mesh1.compute_len_edges()
        np.testing.assert_array_almost_equal(lengths, expected_lengths, decimal=6)

    def test_get_mesh_border(self):
        # Test the get_mesh_border method
        border_edges, all_edges = self.mesh1.get_mesh_border()
        expected_border_edges = np.array([[0, 1], [0, 2], [1, 2]])
        expected_all_edges = np.array([[0, 1], [0, 2], [1, 2]])

        np.testing.assert_array_equal(border_edges, expected_border_edges)
        np.testing.assert_array_equal(all_edges, expected_all_edges)

    def test_remove_large_triangles(self):
        # Test the remove_large_triangles method
        threshold_edge_length = 2.0
        expected_triangles = np.array([[0, 1, 2]])
        new_mesh = self.mesh1.remove_large_triangles(threshold_edge_length)
        np.testing.assert_array_equal(new_mesh.triangles, expected_triangles)

        new_mesh = self.mesh1 + self.mesh1.resize(10)
        self.assertEqual(2, new_mesh.n_triangles)
        new_mesh = new_mesh.remove_large_triangles(threshold_edge_length)
        self.assertEqual(1, new_mesh.n_triangles)

    def test_minimum_distance(self):
        self.assertEqual(0.0, self.mesh1.minimum_distance(self.mesh2))
        self.assertEqual(0.0, self.mesh1.minimum_distance(self.mesh3))
        self.assertEqual(0.0, self.mesh3.minimum_distance(self.mesh1))
        self.assertEqual(1.0, self.mesh1.minimum_distance(self.mesh4))
        self.assertEqual(1.0, self.mesh4.minimum_distance(self.mesh1))


class TestMesh3DImport(unittest.TestCase):
    def setUp(self) -> None:
        self.stl_file_path = os.path.join(FOLDER, "models", "simple.stl")

    def test_from_trimesh(self):
        mesh = Mesh3D.from_trimesh(trimesh.Trimesh(vertices=[[0, 0, 0], [0, 0, 1], [0, 1, 0]], faces=[[0, 1, 2]]))

        if SHOW_BABYLONJS:
            mesh.babylonjs()

    def test_from_stl_file(self):
        mesh = Mesh3D.from_stl_file(self.stl_file_path)

        if SHOW_BABYLONJS:
            mesh.babylonjs()

    def test_from_stl_stream(self):
        with open(self.stl_file_path, "rb") as file:
            binary_content = BinaryFile()
            binary_content.write(file.read())

        mesh = Mesh3D.from_stl_stream(binary_content)

        if SHOW_BABYLONJS:
            mesh.babylonjs()


class TestMesh3DExport(unittest.TestCase):
    def setUp(self) -> None:
        self.mesh = Mesh3D(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, 1.0, 0.0]]), np.array([[0, 1, 2]]))
        self.degenerated_mesh = Mesh3D(np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 1.0]]), np.array([[0, 1, 1]]))

        self.triangle3d = Triangle3D(Point3D(*[0, 0, 0]), Point3D(*[0, 0, 1]), Point3D(*[0, 1, 0]))

    def test_to_triangles3d(self):
        self.assertEqual([self.triangle3d], self.mesh.to_triangles3d())
        self.assertEqual([], self.degenerated_mesh.to_triangles3d())

    def test_plot(self):
        self.mesh.plot()
        self.degenerated_mesh.plot()

    def test_to_closed_shell(self):
        self.assertEqual(ClosedTriangleShell3D([self.triangle3d]), self.mesh.to_closed_shell())
        self.assertEqual(ClosedTriangleShell3D([]), self.degenerated_mesh.to_closed_shell())

    def test_to_open_shell(self):
        self.assertEqual(OpenTriangleShell3D([self.triangle3d]), self.mesh.to_open_shell())
        self.assertEqual(OpenTriangleShell3D([]), self.degenerated_mesh.to_open_shell())

    def test_to_trimesh(self):
        trimesh_ = trimesh.Trimesh(vertices=[[0, 0, 0], [0, 0, 1], [0, 1, 0]], faces=[[0, 1, 2]])

        np.testing.assert_array_equal(trimesh_.vertices, self.mesh.to_trimesh().vertices)
        np.testing.assert_array_equal(trimesh_.faces, self.mesh.to_trimesh().faces)

    def test_save_to_stl_file(self):
        # Create a temporary STL file for testing
        with tempfile.NamedTemporaryFile(suffix=".stl", delete=False) as temp_stl_file:
            temp_stl_filename = temp_stl_file.name

            # Call the save_to_stl_file method to write to the temporary file
            self.mesh.save_to_stl_file(temp_stl_filename)

            # Call the from_stl_file method to read from the temporary file
            new_mesh = Mesh3D.from_stl_file(temp_stl_filename)

            # Assert that the original and new meshes are equal
            self.assertEqual(self.mesh, new_mesh)

        # Clean up the temporary file after the test
        os.remove(temp_stl_filename)

    def test_save_to_stl_stream(self):
        stream = BinaryFile()
        self.mesh.save_to_stl_stream(stream)
        mesh_from_stream = Mesh3D.from_stl_stream(stream)

        self.assertEqual(self.mesh, mesh_from_stream)

    def test_save_to_obj_file(self):
        # Create a temporary OBJ file for testing
        with tempfile.NamedTemporaryFile(suffix=".obj", delete=False) as temp_obj_file:
            temp_obj_filename = temp_obj_file.name

            # Call the save_to_obj_file method to write to the temporary file
            self.mesh.save_to_obj_file(temp_obj_filename)

            # Call the from_obj_file method to read from the temporary file
            new_mesh = Mesh3D.from_obj_file(temp_obj_filename)

            # Assert that the original and new meshes are equal
            self.assertEqual(self.mesh, new_mesh)

        # Clean up the temporary file after the test
        os.remove(temp_obj_filename)

    def test_save_to_obj_stream(self):
        stream = BinaryFile()
        self.mesh.save_to_obj_stream(stream)
        mesh_from_stream = Mesh3D.from_obj_stream(stream)

        self.assertEqual(self.mesh, mesh_from_stream)

    def test_save_to_ply_file(self):
        # Create a temporary PLY file for testing
        with tempfile.NamedTemporaryFile(suffix=".ply", delete=False) as temp_ply_file:
            temp_ply_filename = temp_ply_file.name

            # Call the save_to_ply_file method to write to the temporary file
            self.mesh.save_to_ply_file(temp_ply_filename)

            # Call the from_ply_file method to read from the temporary file
            new_mesh = Mesh3D.from_ply_file(temp_ply_filename)

            # Assert that the original and new meshes are equal
            self.assertEqual(self.mesh, new_mesh)

        # Clean up the temporary file after the test
        os.remove(temp_ply_filename)

    def test_save_to_ply_stream(self):
        stream = BinaryFile()
        self.mesh.save_to_ply_stream(stream)
        mesh_from_stream = Mesh3D.from_ply_stream(stream)

        self.assertEqual(self.mesh, mesh_from_stream)

    def test_save_to_off_file(self):
        # Create a temporary OFF file for testing
        with tempfile.NamedTemporaryFile(suffix=".off", delete=False) as temp_off_file:
            temp_off_filename = temp_off_file.name

            # Call the save_to_off_file method to write to the temporary file
            self.mesh.save_to_off_file(temp_off_filename)

            # Call the from_off_file method to read from the temporary file
            new_mesh = Mesh3D.from_off_file(temp_off_filename)

            # Assert that the original and new meshes are equal
            self.assertEqual(self.mesh, new_mesh)

        # Clean up the temporary file after the test
        os.remove(temp_off_filename)

    def test_save_to_off_stream(self):
        stream = BinaryFile()
        self.mesh.save_to_off_stream(stream)
        mesh_from_stream = Mesh3D.from_off_stream(stream)

        self.assertEqual(self.mesh, mesh_from_stream)

    def test_save_to_3mf_file(self):
        # Create a temporary 3MF file for testing
        with tempfile.NamedTemporaryFile(suffix=".3mf", delete=False) as temp_3mf_file:
            temp_3mf_filename = temp_3mf_file.name

            # Call the save_to_3mf_file method to write to the temporary file
            self.mesh.save_to_3mf_file(temp_3mf_filename)

            # Call the from_3mf_file method to read from the temporary file
            new_mesh = Mesh3D.from_3mf_file(temp_3mf_filename)

            # Assert that the original and new meshes are equal
            self.assertEqual(self.mesh, new_mesh)

        # Clean up the temporary file after the test
        os.remove(temp_3mf_filename)

    def test_save_to_3mf_stream(self):
        stream = BinaryFile()
        self.mesh.save_to_3mf_stream(stream)
        mesh_from_stream = Mesh3D.from_3mf_stream(stream)

        self.assertEqual(self.mesh, mesh_from_stream)


if __name__ == "__main__":
    unittest.main()
