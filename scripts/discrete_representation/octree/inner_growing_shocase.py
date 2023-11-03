"""
Showcase of inner growing voxelization.
"""
import random

from volmdlr.discrete_representation import OctreeBasedVoxelization
from volmdlr.core import VolumeModel

LAYERS_MINIMAL_THICKNESS = 5

octree = [[[[[[1] for _ in range(8)] for _ in range(8)] for _ in range(8)] for _ in range(8)] for _ in range(8)]
# octree[0][0][0][0][0] = []


# Create the voxelization
octree_based_voxelization = OctreeBasedVoxelization(
    octree=octree,
    root_center=(0, 0, 0),
    octree_depth=5,
    voxel_size=0.1,
    triangles=[],
)

# Get inner growing point based voxelization
inner_growing_voxelization = octree_based_voxelization.to_inner_growing_point_based_voxelizations(
    LAYERS_MINIMAL_THICKNESS
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
