"""
Comparison of "iterative" and "octree" voxelization methods.
"""
import time

import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.primitives3d import Cylinder, Sphere
from volmdlr.voxelization import Voxelization

VOXEL_SIZE = 0.01

# Create a volume model
sphere = Sphere(volmdlr.O3D, 0.1, name="Sphere")
cylinder = Cylinder(volmdlr.OXYZ, 0.1, 0.2, name="Cylinder")
volume_model = VolumeModel([sphere, cylinder])

# Voxelize the volume model with "iterative" method
start = time.perf_counter()
voxelization_iterative = Voxelization.from_volume_model(
    volume_model, VOXEL_SIZE, method="iterative", name="Iterative voxelization"
)
print(f"Iterative voxelization computed in {round((time.perf_counter() - start) * 1000)}ms")

# Voxelize the volume model with "octree" method
start = time.perf_counter()
voxelization_octree = Voxelization.from_volume_model(
    volume_model, VOXEL_SIZE, method="octree", name="Octree voxelization"
)
print(f"Octree voxelization computed in {round((time.perf_counter() - start) * 1000)}ms")

# Check for difference between voxelization
print(f"Both voxelization are equal: {voxelization_iterative == voxelization_octree}")

# Display the difference if there is one
if voxelization_iterative != voxelization_octree:
    symmetric_difference = voxelization_iterative.symmetric_difference(voxelization_octree)
    symmetric_difference.babylonjs()
