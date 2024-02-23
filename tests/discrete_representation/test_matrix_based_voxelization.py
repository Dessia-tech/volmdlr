"""
Unit testing of 'MatrixBasedVoxelization' class.
"""
import os
import unittest

import volmdlr
from volmdlr.core import BoundingBox, VolumeModel
from volmdlr.discrete_representation import MatrixBasedVoxelization
from volmdlr.primitives3d import Block, Cylinder, Sphere
from volmdlr.shells import ClosedTriangleShell3D, DisplayTriangleShell3D

SHOW_BABYLONJS = False


class TestMatrixBasedVoxelizationCreation(unittest.TestCase):
    """
    Unit testing of voxelization creation methods.
    """

    def setUp(self):
        self.block = Block(frame=volmdlr.OXYZ, name="block")
        self.sphere = Sphere(center=volmdlr.Point3D(0.0, 0.0, 0.1), radius=0.1, name="sphere")
        self.cylinder = Cylinder(frame=volmdlr.OXYZ, radius=0.1, length=0.2, name="cylinder")
        self.volume_model = VolumeModel(primitives=[self.sphere, self.cylinder], name="volume model")

    def test_voxelize_block(self):
        block_voxelization = MatrixBasedVoxelization.from_shell(self.block, 0.2, name="voxelization")
        self.assertEqual(152, len(block_voxelization))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel([self.block, block_voxelization.to_closed_triangle_shell()])
            volume_model.babylonjs()

    def test_voxelize_translated_block(self):
        translated_block = self.block.translation(volmdlr.Vector3D(11, 1.8, 4.8))
        translated_block_voxelization = MatrixBasedVoxelization.from_shell(translated_block, 0.1, name="voxelization")
        self.assertEqual(1216, len(translated_block_voxelization))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel([translated_block, translated_block_voxelization.to_closed_triangle_shell()])
            volume_model.babylonjs()

    def test_voxelize_sphere(self):
        sphere_voxelization = MatrixBasedVoxelization.from_shell(self.sphere, 0.01, name="sphere voxelization")
        self.assertEqual(1900, len(sphere_voxelization))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel([self.sphere, sphere_voxelization.to_closed_triangle_shell()])
            volume_model.babylonjs()

    def test_voxelize_cylinder(self):
        cylinder_voxelization = MatrixBasedVoxelization.from_shell(self.cylinder, 0.01, name="cylinder voxelization")
        self.assertEqual(2788, len(cylinder_voxelization))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel([self.cylinder, cylinder_voxelization.to_closed_triangle_shell()])
            volume_model.babylonjs()

    def test_from_volume_model(self):
        volume_model_voxelization = MatrixBasedVoxelization.from_volume_model(
            self.volume_model, 0.01, name="voxelization"
        )
        self.assertEqual(4288, len(volume_model_voxelization))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel(
                self.volume_model.primitives + [volume_model_voxelization.to_closed_triangle_shell()]
            )
            volume_model.babylonjs()

    def test_from_mesh_data(self):
        file_path = os.path.join(os.path.dirname(os.path.realpath(__file__)), "stanford_bunny.json")
        stanford_bunny = DisplayTriangleShell3D.from_json(file_path)
        voxelization1 = MatrixBasedVoxelization.from_mesh_data(stanford_bunny.positions, stanford_bunny.indices, 0.005)
        voxelization2 = MatrixBasedVoxelization.from_shell(stanford_bunny, 0.005)

        self.assertEqual(1735, len(voxelization1))
        self.assertEqual(1735, len(voxelization2))


class TestMatrixBasedVoxelizationBooleanOperation(unittest.TestCase):
    """
    Unit testing of voxelization boolean operation methods.
    """

    def setUp(self):
        self.sphere = Sphere(center=volmdlr.Point3D(0.0, 0.0, 0.1), radius=0.1, name="sphere")
        self.cylinder = Cylinder(frame=volmdlr.OXYZ, radius=0.1, length=0.2, name="cylinder")
        self.volume_model = VolumeModel(primitives=[self.sphere, self.cylinder], name="volume model")

        self.sphere_voxelization = MatrixBasedVoxelization.from_shell(self.sphere, 0.01, name="sphere voxelization")
        self.cylinder_voxelization = MatrixBasedVoxelization.from_shell(
            self.cylinder, 0.01, name="cylinder voxelization"
        )
        self.volume_model_voxelization = MatrixBasedVoxelization.from_volume_model(
            self.volume_model, 0.01, name="volume model voxelization"
        )

    def test_union_1(self):
        union_1 = self.sphere_voxelization.union(self.cylinder_voxelization)
        union_2 = self.cylinder_voxelization.union(self.sphere_voxelization)

        self.assertEqual(len(union_1), len(union_2))
        self.assertEqual(4288, len(union_1))
        self.assertEqual(self.volume_model_voxelization, union_1)

        if SHOW_BABYLONJS:
            volume_model = VolumeModel([self.sphere, self.cylinder, union_1.to_closed_triangle_shell()])
            volume_model.babylonjs()

    def test_union_2(self):
        block_1 = Block(frame=volmdlr.OXYZ)
        block_2 = block_1.translation(volmdlr.Vector3D(0.5, 0.5, 0.5))

        voxelization_1 = MatrixBasedVoxelization.from_shell(block_1, 0.1)
        voxelization_2 = MatrixBasedVoxelization.from_shell(block_2, 0.1)

        union_3 = voxelization_1.union(voxelization_2)
        self.assertEqual(2312, len(union_3))

        if SHOW_BABYLONJS:
            VolumeModel(
                [voxelization_1.to_closed_triangle_shell(), voxelization_2.to_closed_triangle_shell()]
            ).babylonjs()
            union_3.babylonjs()

    def test_intersection(self):
        intersection_1 = self.sphere_voxelization.intersection(self.cylinder_voxelization)
        intersection_2 = self.cylinder_voxelization.intersection(self.sphere_voxelization)

        self.assertEqual(len(intersection_1), len(intersection_2))
        self.assertEqual(400, len(intersection_1))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel([self.sphere, self.cylinder, intersection_1.to_closed_triangle_shell()])
            volume_model.babylonjs()

    def test_difference(self):
        difference_1 = self.sphere_voxelization.difference(self.cylinder_voxelization)
        difference_2 = self.cylinder_voxelization.difference(self.sphere_voxelization)

        self.assertNotEqual(len(difference_1), len(difference_2))
        self.assertEqual(1500, len(difference_1))
        self.assertEqual(2388, len(difference_2))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel(
                [
                    self.sphere,
                    self.cylinder,
                    difference_1.to_closed_triangle_shell(),
                    difference_2.to_closed_triangle_shell(),
                ]
            )
            volume_model.babylonjs()

    def test_symmetric_difference(self):
        symmetric_difference_1 = self.sphere_voxelization.symmetric_difference(self.cylinder_voxelization)
        symmetric_difference_2 = self.cylinder_voxelization.symmetric_difference(self.sphere_voxelization)

        self.assertEqual(len(symmetric_difference_1), len(symmetric_difference_2))
        self.assertEqual(3888, len(symmetric_difference_1))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel([self.sphere, self.cylinder, symmetric_difference_1.to_closed_triangle_shell()])
            volume_model.babylonjs()


class TestMatrixBasedVoxelizationManipulation(unittest.TestCase):
    """
    Unit testing of voxelization manipulation methods.
    """

    def setUp(self):
        self.cylinder = Cylinder(frame=volmdlr.OXYZ, radius=0.1, length=0.2, name="cylinder")
        self.cylinder_voxelization = MatrixBasedVoxelization.from_shell(
            self.cylinder, 0.01, name="cylinder voxelization"
        )

    def test_inverse(self):
        inverse_cylinder_voxelization = self.cylinder_voxelization.inverse()

        self.assertEqual(6452, len(inverse_cylinder_voxelization))
        self.assertEqual(0, len(self.cylinder_voxelization.intersection(inverse_cylinder_voxelization)))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel([self.cylinder, inverse_cylinder_voxelization.to_closed_triangle_shell()])
            volume_model.babylonjs()

    def test_fill_outer_voxels(self):
        outer_filled_voxelization = self.cylinder_voxelization.fill_outer_voxels()

        self.assertEqual(4416, len(outer_filled_voxelization))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel([self.cylinder, outer_filled_voxelization.to_closed_triangle_shell()])
            volume_model.babylonjs()

    def test_fill_enclosed_voxels(self):
        inner_filled_voxelization = self.cylinder_voxelization.fill_enclosed_voxels()

        self.assertEqual(7612, len(inner_filled_voxelization))

        if SHOW_BABYLONJS:
            volume_model = VolumeModel([self.cylinder, inner_filled_voxelization.to_closed_triangle_shell()])
            volume_model.babylonjs()


class TestMatrixBasedVoxelizationExport(unittest.TestCase):
    """
    Unit testing of voxelization export methods.
    """

    def setUp(self):
        self.sphere = Sphere(center=volmdlr.O3D, radius=0.1, name="sphere")
        self.sphere_voxelization = MatrixBasedVoxelization.from_shell(self.sphere, 0.01, name="sphere voxelization")

    def test_min_grid_center(self):
        self.assertEqual((-0.105, -0.105, -0.105), self.sphere_voxelization.min_grid_center)

    def test_max_grid_center(self):
        self.assertEqual((0.095, 0.105, 0.105), self.sphere_voxelization.max_grid_center)

    def test_bounding_box(self):
        self.assertEqual(
            BoundingBox(-0.115, 0.105, -0.115, 0.115, -0.115, 0.115), self.sphere_voxelization.bounding_box
        )

    def test_to_triangles(self):
        self.assertEqual(7504, len(self.sphere_voxelization.to_triangles()))

    def test_to_closed_triangle_shell(self):
        self.assertIsInstance(self.sphere_voxelization.to_closed_triangle_shell(), ClosedTriangleShell3D)

    def test_to_point_voxelization(self):
        point_voxelization = self.sphere_voxelization.to_point_based_voxelization()

        self.assertEqual(len(self.sphere_voxelization), len(point_voxelization))
        self.assertEqual(self.sphere_voxelization.min_grid_center, point_voxelization.min_grid_center)

    def test_serialization(self):
        self.assertEqual(
            self.sphere_voxelization, MatrixBasedVoxelization.dict_to_object(self.sphere_voxelization.to_dict())
        )


if __name__ == "__main__":
    unittest.main()
