"""
Showcase of filling a voxelization.
"""
import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.discrete_representation import MatrixBasedVoxelization
from volmdlr.primitives3d import Cylinder

VOXEL_SIZE = 0.005

# Create a volume model
cylinder = Cylinder(volmdlr.OXYZ, 0.1, 0.2, name="Cylinder")
volume_model = VolumeModel([cylinder])

# Voxelize the volume model (it uses the triangulated model to create the voxelization)
voxelization = MatrixBasedVoxelization.from_volume_model(volume_model, VOXEL_SIZE, name="Voxelization")

# Fill the voxelization
filled_voxelization = voxelization.fill_enclosed_voxels()  # filling the enclosed voxels of the voxelization
filled_voxelization.name = "Filled voxelization"

# Print the volumes
# The filled voxelization volume will be always over-estimated,
# but become more precise as we reduce the voxelization size

print(f"Cylinder volume: {cylinder.volume()}")
print(f"Voxelization volume: {voxelization.volume}")
print(f"Filled voxelization volume: {filled_voxelization.volume}")

# Display the result
voxelization_cs = voxelization.to_closed_triangle_shell()
voxelization_cs.color = (1, 0, 0)
voxelization_cs.alpha = 0.5

filled_voxelization_cs = filled_voxelization.to_mesh().split_shared_vertices()
filled_voxelization_cs.color = (0, 0, 1)
filled_voxelization_cs.alpha = 0.5

volume_model = VolumeModel([voxelization_cs, filled_voxelization_cs])
volume_model.babylonjs()
