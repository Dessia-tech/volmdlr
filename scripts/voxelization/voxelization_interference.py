"""
Demo script for interference computing using voxelization
"""
import math
import time

import volmdlr
from volmdlr.voxelization import Voxelization
from volmdlr.core import VolumeModel
from volmdlr.stl import Stl

STL_MODEL_FILE_PATH = "../stl/simple.stl"
VOXEL_SIZE = 0.002
TRANSLATION_VECTOR = volmdlr.Vector3D(-0.0015, 0.005, 0.01)
ROTATION_ANGLE = math.pi / 2

# Define the volume model from file
volume_model = Stl.load_from_file(STL_MODEL_FILE_PATH).to_volume_model()
volume_model.primitives[0].color = (0, 1, 0)

# Move the volume model
moved_volume_model = volume_model.rotation(volmdlr.O3D, volmdlr.X3D, ROTATION_ANGLE)
moved_volume_model = moved_volume_model.translation(TRANSLATION_VECTOR)
moved_volume_model.primitives[0].color = (0, 0, 1)

# Voxelize the volume model
voxelization = Voxelization.from_volume_model(volume_model, VOXEL_SIZE, method="iterative")

# Move the voxelization
start = time.perf_counter()
moved_voxelization = voxelization.rotation(volmdlr.O3D, volmdlr.X3D, ROTATION_ANGLE)
print(f"Voxels rotation computing time: {(time.perf_counter() - start)*1000}ms")

start = time.perf_counter()
moved_voxelization = moved_voxelization.translation(TRANSLATION_VECTOR)
print(f"Voxels translation computing time: {(time.perf_counter() - start)*1000}ms")

# Compute the intersection
start = time.perf_counter()
voxelization_intersection = voxelization.intersection(moved_voxelization)
print(f"\nIntersection computing time: {(time.perf_counter() - start)*1000}ms")

# Check if the voxelization are intersecting
start = time.perf_counter()
print(f"Voxelizations are intersecting: {voxelization.is_intersecting(moved_voxelization)}")
print(f"Voxelization.is_intersecting computing time: {(time.perf_counter() - start)*1000}ms")

print(f"Voxelization interference percentage: {voxelization.interference(moved_voxelization)*100}%")

# ClosedTriangleShell3D is_intersecting_with check
start = time.perf_counter()
closed_triangle_shell_intersection = volume_model.primitives[0].is_intersecting_with(moved_volume_model.primitives[0])
print(f"\nClosedTriangleShell3D are intersecting: {closed_triangle_shell_intersection}")
print(f"ClosedTriangleShell3D.is_intersecting_with computing time: {(time.perf_counter() - start)*1000}ms")

# Display the result
voxelization_intersection_cs = voxelization_intersection.to_closed_triangle_shell()
voxelization_intersection_cs.color = (1, 0, 0)

volume_model = VolumeModel(volume_model.primitives + moved_volume_model.primitives + [voxelization_intersection_cs])
volume_model.babylonjs()
