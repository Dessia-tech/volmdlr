"""
Demo script for interference computing using voxelization
"""
import math
import time

import volmdlr
from volmdlr.voxelization import Voxelization
from volmdlr.core import VolumeModel
from volmdlr.stl import Stl

VOXEL_SIZE = 0.05

# Define the volume model from file
volume_model = volmdlr.stl.Stl.load_from_file("model.stl").to_volume_model()
volume_model.primitives[0].color = (0, 1, 0)

# Move the volume model
moved_volume_model = volume_model.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 2)
moved_volume_model = moved_volume_model.translation(volmdlr.Vector3D(0.0, 0.5, 0.3))
moved_volume_model.primitives[0].color = (0, 0, 1)

# Voxelize the volume model
voxelization = Voxelization.from_volume_model(volume_model, VOXEL_SIZE, method="iterative")

# Move the voxelization
start = time.perf_counter()
moved_voxelization = voxelization.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 2)
print(f"Voxels rotation computing time: {(time.perf_counter() - start)*1000}ms")

start = time.perf_counter()
moved_voxelization = moved_voxelization.translation(volmdlr.Vector3D(0.0, 0.5, 0.3))
print(f"Voxels translation computing time: {(time.perf_counter() - start)*1000}ms")

# Compute the intersection
start = time.perf_counter()
voxelization_intersection = voxelization.intersection(moved_voxelization)
print(f"Intersection computing time: {(time.perf_counter() - start)*1000}ms")

# Display the result
voxelization_intersection_cs = voxelization_intersection.to_closed_triangle_shell()
voxelization_intersection_cs.color = (1, 0, 0)

volume_model = VolumeModel(volume_model.primitives + moved_volume_model.primitives + [voxelization_intersection_cs])
volume_model.babylonjs()

print(f"Interfence percentage: {voxelization.interference(moved_voxelization)*100}%")
