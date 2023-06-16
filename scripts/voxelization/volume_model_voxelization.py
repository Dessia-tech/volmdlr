"""
Voxelization of a closed shell.
"""
from volmdlr.voxelization import Voxelization
from volmdlr.primitives3d import Sphere, Cylinder
from volmdlr.core import VolumeModel
import volmdlr

# Create a volume model
sphere = Sphere(volmdlr.O3D, 0.1, name="Sphere")
cylinder = Cylinder(volmdlr.Point3D(0.0, 0.0, 0.1), volmdlr.Z3D, 0.1, 0.2, name="Cylinder")
volume_model = VolumeModel([sphere, cylinder])

# Voxelize the volume model (it uses the triangulated model to create the voxelization)
voxelization = Voxelization.from_volume_model(volume_model, 0.005, name="Voxels")

# Display the result
volume_model.primitives.extend(voxelization.volmdlr_primitives())
volume_model.babylonjs()

# Compare triangulation methods
voxelization.to_closed_shell().babylonjs()
voxelization.to_closed_triangle_shell().babylonjs()
