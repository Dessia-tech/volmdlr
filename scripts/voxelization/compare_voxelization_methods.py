"""
Comparaison of "naive" and "octree" voxelization methods.
"""
from volmdlr.voxelization import Voxelization
from volmdlr.primitives3d import Sphere, Cylinder
from volmdlr.core import VolumeModel
import volmdlr
import time

VOXEL_SIZE = 0.1

# Create a volume model
sphere = Sphere(volmdlr.O3D, 0.1, name="Sphere")
cylinder = Cylinder(volmdlr.Point3D(0.0, 0.0, 0.1), volmdlr.Z3D, 0.1, 0.2, name="Cylinder")
volume_model = VolumeModel([sphere, cylinder])

# Voxelize the volume model with "naive" method
start = time.perf_counter()
voxelization_naive = Voxelization.from_volume_model(volume_model, VOXEL_SIZE, method="naive", name="Voxels")
print(f"Naive voxelization computed in {round((time.perf_counter() - start) * 1000)}ms")

# Voxelize the volume model with "octree" method
start = time.perf_counter()
voxelization_octree = Voxelization.from_volume_model(volume_model, VOXEL_SIZE, method="octree", name="Voxels")
print(f"Octree voxelization computed in {round((time.perf_counter() - start) * 1000)}ms")

print(f"Both voxelization are equal: {voxelization_naive == voxelization_octree}")
