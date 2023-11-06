"""
Shell intersection detection.
"""
import math

import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.step import Step
from volmdlr.discrete_representation import OctreeBasedVoxelization

# Load
shell_1 = Step.from_file("../../step/block.step").to_volume_model().get_shells()[0]
shell_1.color = (0, 1, 0)
# shell2 = shell1.translation(volmdlr.Vector3D(0.01, 0.01, 0.01))
shell_2 = shell_1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi / 4)
shell_2.color = (0, 0, 1)

volume_model = VolumeModel([shell_1, shell_2])
volume_model.babylonjs()

# face_combinations = shell1.intersecting_faces_combinations(shell2)

VOXEL_SIZE = 0.001
voxelization_1 = OctreeBasedVoxelization.from_shell(shell_1, VOXEL_SIZE)
voxelization_2 = OctreeBasedVoxelization.from_shell(shell_2, VOXEL_SIZE)
intersection = voxelization_1.intersection(voxelization_2)
