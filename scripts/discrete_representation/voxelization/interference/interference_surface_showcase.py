"""
Showcase of surface interference computing using voxelization.
"""
import math
import time

import volmdlr
from volmdlr.model import VolumeModel
from volmdlr.discrete_representation import PointBasedVoxelization
from volmdlr.step import Step

STEP_MODEL_FILE_PATH = "../../../step/tore1.step"
VOXEL_SIZE = 0.0001
TRANSLATION_VECTOR = volmdlr.Vector3D(-0.0015, 0.005, 0.01)
ROTATION_ANGLE = math.pi / 2

# Define the volume model from file
volume_model = Step.from_file(STEP_MODEL_FILE_PATH).to_volume_model()
volume_model.primitives[0].color = (0, 1, 0)

# Move the volume model
moved_volume_model = volume_model.rotation(volmdlr.O3D, volmdlr.X3D, ROTATION_ANGLE)
moved_volume_model = moved_volume_model.translation(TRANSLATION_VECTOR)
moved_volume_model.primitives[0].color = (0, 0, 1)

# Voxelize the volume models
start = time.perf_counter()
voxelization = PointBasedVoxelization.from_volume_model(volume_model, VOXEL_SIZE)
voxelization_moved = PointBasedVoxelization.from_volume_model(moved_volume_model, VOXEL_SIZE)
print(f"\nBoth voxelizations computing time: {(time.perf_counter() - start)*1000}ms")

# Compute the intersection
start = time.perf_counter()
voxelization_intersection = voxelization.intersection(voxelization_moved)
print(f"\nIntersection computing time: {(time.perf_counter() - start)*1000}ms")

# Check if the voxelization are intersecting
start = time.perf_counter()
print(f"Voxelizations are intersecting: {voxelization.is_intersecting(voxelization_moved)}")
print(f"Voxelization.is_intersecting computing time: {(time.perf_counter() - start)*1000}ms")
print(f"Voxelization interference percentage: {voxelization.interference(voxelization_moved)*100}%")

# Display the result
voxelization_intersection_cs = voxelization_intersection.to_closed_triangle_shell()
voxelization_intersection_cs.color = (1, 0, 0)

volume_model = VolumeModel(volume_model.primitives + moved_volume_model.primitives + [voxelization_intersection_cs])
volume_model.babylonjs()
