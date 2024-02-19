"""
Demo script for moving voxelization
"""
import math
import time

import volmdlr
from volmdlr.discrete_representation import PointBasedVoxelization
from volmdlr.stl import Stl

STL_MODEL_FILE_PATH = "../../stl/simple.stl"
VOXEL_SIZE = 0.0025
TRANSLATION_VECTOR = volmdlr.Vector3D(0.01, 0.02, -0.01)
ROTATION_ANGLE = math.pi / 2

# Define the volume model from file
volume_model = Stl.load_from_file(STL_MODEL_FILE_PATH).to_volume_model()
volume_model.color = (1, 0, 0)

# Move the volume model
start = time.perf_counter()
moved_volume_model = volume_model.rotation(volmdlr.O3D, volmdlr.X3D, ROTATION_ANGLE)
print(f"Volume model rotation computing time: {(time.perf_counter() - start)*1000}ms")
start = time.perf_counter()
moved_volume_model = moved_volume_model.translation(TRANSLATION_VECTOR)
print(f"Volume model translation computing time: {(time.perf_counter() - start)*1000}ms")
moved_volume_model.color = (0, 0, 1)

# Define a voxelization from the none-moved volume_model and move it
voxelization = PointBasedVoxelization.from_volume_model(volume_model, VOXEL_SIZE)

start = time.perf_counter()
moved_voxelization = voxelization.rotation(volmdlr.O3D, volmdlr.X3D, ROTATION_ANGLE)
print(f"Voxels rotation computing time: {(time.perf_counter() - start)*1000}ms")

start = time.perf_counter()
moved_voxelization = moved_voxelization.translation(TRANSLATION_VECTOR)
print(f"Voxels translation computing time: {(time.perf_counter() - start)*1000}ms")

# Define a voxelization from the moved volume_model to compare it
voxelization_from_moved_volume_model = PointBasedVoxelization.from_volume_model(
    moved_volume_model,
    VOXEL_SIZE,
)

# Display the result
moved_volume_model.primitives.append(moved_voxelization.to_mesh().split_shared_vertices())
moved_volume_model.primitives[-1].alpha = 0.5
moved_volume_model.babylonjs()

# Display the difference between the voxelization from moved volume model and the moved voxelization
equality = voxelization_from_moved_volume_model == moved_voxelization
print(f"Moved volume model voxelization and moved voxelization are equal: {equality}")
if not equality:
    print("Displaying the symmetric difference")
    symmetric_diff = moved_voxelization ^ voxelization_from_moved_volume_model
    # equivalent to moved_voxelization.symmetric_difference(voxelization_from_moved_volume_model)

    symmetric_diff.babylonjs()
