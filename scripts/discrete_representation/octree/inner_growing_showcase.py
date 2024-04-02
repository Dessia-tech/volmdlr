"""
Showcase of inner growing voxelization.
"""
import random

import volmdlr
from volmdlr.model import VolumeModel
from volmdlr.discrete_representation import MatrixBasedVoxelization
from volmdlr.primitives3d import Block

VOXEL_SIZE = 0.05
LAYERS_MINIMAL_THICKNESS = 2 * VOXEL_SIZE


# Create the geometry
block = Block(volmdlr.OXYZ)

# Create the inner grwoing voxelization
inner_growing_voxelization = (
    MatrixBasedVoxelization.from_shell(block, VOXEL_SIZE)
    .fill_enclosed_voxels()
    .to_inner_growing_voxelizations(LAYERS_MINIMAL_THICKNESS)
)


primitives = []
for vox in inner_growing_voxelization:
    print(vox.voxel_size)
    shell = vox.to_closed_triangle_shell()
    shell.color = (random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1))
    shell.alpha = 0.5
    shell.name = str(vox.voxel_size)
    primitives.append(shell)

# Display
volume_model = VolumeModel(primitives)
volume_model.babylonjs()
