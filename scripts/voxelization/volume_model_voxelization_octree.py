"""
Voxelization of a volume model using "octree" method.
"""
from volmdlr.voxelization import Voxelization
from volmdlr.primitives3d import Sphere, Cylinder
from volmdlr.core import VolumeModel
import volmdlr

VOXEL_SIZE = 0.01

# Create a volume model
sphere = Sphere(volmdlr.O3D, 0.1, name="Sphere")
cylinder = Cylinder(volmdlr.Point3D(0.0, 0.0, 0.1), volmdlr.Z3D, 0.1, 0.2, name="Cylinder")
volume_model = VolumeModel([sphere, cylinder])

# Voxelize the volume model (it uses the triangulated model to create the voxelization)
voxelization = Voxelization.from_volume_model(volume_model, VOXEL_SIZE, method="octree", name="Voxels")

# Display the result
volume_model.primitives.append(voxelization.to_closed_triangle_shell())
volume_model.babylonjs()
