"""
Showcase of computing the inverse of a voxelization.
"""
import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.discrete_representation import MatrixBasedVoxelization
from volmdlr.primitives3d import Cylinder

VOXEL_SIZE = 0.01

# Create a volume model
cylinder = Cylinder(volmdlr.OXYZ, 0.1, 0.2, name="Cylinder")
volume_model = VolumeModel([cylinder])

# Voxelize the volume model (it uses the triangulated model to create the voxelization)
voxelization = MatrixBasedVoxelization.from_volume_model(volume_model, VOXEL_SIZE, name="Voxelization")

# Compute the inverse of the voxelization
voxelization_inverse = ~voxelization  # equivalent to voxelization.inverse()
voxelization_inverse.name = "Invert voxelization"

# Display the result
voxelization_cs = voxelization.to_mesh().split_shared_vertices()
voxelization_cs.color = (1, 0, 0)

voxelization_inverse_cs = voxelization_inverse.to_mesh().split_shared_vertices()
voxelization_inverse_cs.color = (0, 0, 1)

volume_model.primitives.extend([voxelization_cs, voxelization_inverse_cs])
volume_model.babylonjs()
