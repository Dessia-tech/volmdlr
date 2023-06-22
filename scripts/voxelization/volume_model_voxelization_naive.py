"""
Voxelization of a volume model using "naive" method.
"""
import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.primitives3d import Cylinder, Sphere
from volmdlr.voxelization import Voxelization

VOXEL_SIZE = 0.005

# Create a volume model
sphere = Sphere(volmdlr.O3D, 0.1, name="Sphere")
cylinder = Cylinder(volmdlr.Point3D(0.0, 0.0, 0.1), volmdlr.Z3D, 0.1, 0.2, name="Cylinder")
volume_model = VolumeModel([sphere, cylinder])

# Voxelize the volume model (it uses the triangulated model to create the voxelization)
voxelization = Voxelization.from_volume_model(volume_model, VOXEL_SIZE, method="naive", name="Voxelization")

# Display the result
volume_model.primitives.append(voxelization.to_closed_shell())
volume_model.babylonjs()
