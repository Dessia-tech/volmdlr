"""
Voxelization of a closed shell.
"""
from volmdlr.voxelization import Voxelization
from volmdlr.primitives3d import Sphere
from volmdlr.core import VolumeModel
import volmdlr

# Create a volume model
sphere = Sphere(volmdlr.O3D, 0.1, name="Sphere")
volume_model = VolumeModel([sphere])

# Voxelize the volume model (it uses the triangulated model to create the voxelization)
voxelization = Voxelization.from_volume_model(volume_model, 0.01)

# Display the result
volume_model.primitives.extend(voxelization.volmdlr_primitives())
volume_model.babylonjs()
