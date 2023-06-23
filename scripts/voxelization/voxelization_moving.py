"""
Demo script for moving voxelization
"""
import math
import time

import volmdlr
from volmdlr.voxelization import Voxelization
from volmdlr.stl import Stl

VOXEL_SIZE = 0.02

# Define the volume model from file
volume_model = Stl.load_from_file("model.stl").to_volume_model()
volume_model.color = (1, 0, 0)

# Move the volume model
moved_volume_model = volume_model.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 2.1)
moved_volume_model = moved_volume_model.translation(volmdlr.Vector3D(0.0, 0.5, 0.3))
moved_volume_model.color = (0, 0, 1)

# Define a voxelization from the none-moved volume_model and move it
voxelization = Voxelization.from_volume_model(volume_model, VOXEL_SIZE, method="iterative")

start = time.perf_counter()
moved_voxelization = voxelization.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 2.1)
print(f"Voxels rotation computing time: {(time.perf_counter() - start)*1000}ms")

start = time.perf_counter()
moved_voxelization = moved_voxelization.translation(volmdlr.Vector3D(0.0, 0.5, 0.3))
print(f"Voxels translation computing time: {(time.perf_counter() - start)*1000}ms")

# Define a voxelization from the moved volume_model to compare it
voxelization_from_moved_volume_model = Voxelization.from_volume_model(
    moved_volume_model, VOXEL_SIZE, method="iterative"
)

# Display the result
moved_volume_model.primitives.append(moved_voxelization.to_closed_triangle_shell())
moved_volume_model.babylonjs()

# Display the difference between the voxelization from moved volume model and the moved voxelization
equality = voxelization_from_moved_volume_model == moved_voxelization
print(f"Moved volume model voxelization and moved voxelization are equal: {equality}")
if not equality:
    print("Displaying the symmetric difference")
    symmetric_diff = moved_voxelization ^ voxelization_from_moved_volume_model
    # equivalent to moved_voxelization.symmetric_difference(voxelization_from_moved_volume_model)

    symmetric_diff.babylonjs()
