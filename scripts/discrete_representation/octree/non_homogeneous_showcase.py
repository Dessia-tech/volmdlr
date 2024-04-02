"""
Showcase of non homogeneous voxelization.
"""
import random

import volmdlr
from volmdlr.model import VolumeModel
from volmdlr.discrete_representation import MatrixBasedVoxelization
from volmdlr.primitives3d import Cylinder

VOXEL_SIZE = 0.01

# Create the geometry
cylinder = Cylinder(volmdlr.OXYZ, 0.1, 0.1)

# Create the inner grwoing voxelization
non_homogeneous_voxelization = (
    MatrixBasedVoxelization.from_shell(cylinder, VOXEL_SIZE)
    .fill_enclosed_voxels()
    .to_octree_based_voxelization()
    .to_non_homogeneous_point_based_voxelizations()
)


primitives = []
for vox in non_homogeneous_voxelization:
    print(vox.voxel_size)
    shell = vox.to_closed_triangle_shell()
    shell.color = (random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1))
    shell.alpha = 0.5
    shell.name = str(vox.voxel_size)
    primitives.append(shell)

# Display
volume_model = VolumeModel(primitives)
volume_model.babylonjs()
