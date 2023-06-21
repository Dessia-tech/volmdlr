"""
Showcase of computing the inverse of a voxelization.
"""
import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.primitives3d import Cylinder
from volmdlr.voxelization import Voxelization

VOXEL_SIZE = 0.01

# Create a volume model
cylinder = Cylinder(volmdlr.Point3D(0.0, 0.0, 0.1), volmdlr.Z3D, 0.1, 0.2, name="Cylinder")
volume_model = VolumeModel([cylinder])

# Voxelize the volume model (it uses the triangulated model to create the voxelization)
voxelization = Voxelization.from_volume_model(volume_model, VOXEL_SIZE, name="Voxelization")

# Compute the inverse of the voxelization
voxelization_inverse = ~voxelization  # equivalent to voxelization.inverse()
voxelization_inverse.name = "Inversed voxelization"

# Display the result
voxelization_cs = voxelization.to_closed_triangle_shell()
voxelization_cs.color = (1, 0, 0)
# voxelization_cs.alpha = 0.9

voxelization_inverse_cs = voxelization_inverse.to_closed_triangle_shell()
voxelization_inverse_cs.color = (0, 1, 0)
# voxelization_inverse_cs.alpha = 0.9

volume_model.primitives.extend([voxelization_cs, voxelization_inverse_cs])
volume_model.babylonjs()
