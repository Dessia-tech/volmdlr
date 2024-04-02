"""
Performance comparaison on 'ìs_intersecting' method.
"""
import math
import time

import volmdlr
from volmdlr.model import VolumeModel
from volmdlr.discrete_representation import MatrixBasedVoxelization, OctreeBasedVoxelization, PointBasedVoxelization
from volmdlr.step import Step

STEP_MODEL_FILE_PATH = "../../../step/tore1.step"
VOXEL_SIZE = 0.00005
TRANSLATION_VECTOR = volmdlr.Vector3D(-0.0015, 0.005, 0.01)
ROTATION_ANGLE = math.pi / 2

# Define the volume model from file
volume_model_1 = Step.from_file(STEP_MODEL_FILE_PATH).to_volume_model()
volume_model_1.primitives[0].color = (1, 0, 0)

# Define the second volume model
volume_model_2 = volume_model_1.rotation(volmdlr.O3D, volmdlr.X3D, ROTATION_ANGLE)
volume_model_2 = volume_model_2.translation(TRANSLATION_VECTOR)
volume_model_2.primitives[0].color = (0, 1, 0)

# Define the third volume model
volume_model_3 = volume_model_2.rotation(volmdlr.O3D, volmdlr.X3D, ROTATION_ANGLE)
volume_model_3 = volume_model_3.translation(TRANSLATION_VECTOR)
volume_model_3.primitives[0].color = (0, 0, 1)

# Display
volume_model = VolumeModel(volume_model_1.primitives + volume_model_2.primitives + volume_model_3.primitives)
# volume_model.babylonjs()

print("\nPointBasedVoxelization")
t0 = time.perf_counter()
voxelization_1 = PointBasedVoxelization.from_volume_model(volume_model_1, VOXEL_SIZE)
voxelization_2 = PointBasedVoxelization.from_volume_model(volume_model_2, VOXEL_SIZE)
voxelization_3 = PointBasedVoxelization.from_volume_model(volume_model_3, VOXEL_SIZE)
t1 = time.perf_counter()
print(f"Test 1/2: {voxelization_1.is_intersecting(voxelization_2)}")
t2 = time.perf_counter()
print(f"Test 1/3: {voxelization_1.is_intersecting(voxelization_3)}")
t3 = time.perf_counter()
print(f"Test 2/3: {voxelization_2.is_intersecting(voxelization_3)}")
t4 = time.perf_counter()
print(
    f"Voxelization: {(t1-t0)*1000:.3f}ms, 1/2: {(t2-t1)*1000:.3f}ms, 1/3: {(t3-t2)*1000:.3f}ms, 2/3: {(t4-t3)*1000:.3f}ms"
)

# print("\nMatrixBasedVoxelization")
# t0 = time.perf_counter()
# voxelization_1 = MatrixBasedVoxelization.from_volume_model(volume_model_1, VOXEL_SIZE)
# voxelization_2 = MatrixBasedVoxelization.from_volume_model(volume_model_2, VOXEL_SIZE)
# voxelization_3 = MatrixBasedVoxelization.from_volume_model(volume_model_3, VOXEL_SIZE)
# t1 = time.perf_counter()
# print(f"Test 1/2: {voxelization_1.is_intersecting(voxelization_2)}")
# t2 = time.perf_counter()
# print(f"Test 1/3: {voxelization_1.is_intersecting(voxelization_3)}")
# t3 = time.perf_counter()
# print(f"Test 2/3: {voxelization_2.is_intersecting(voxelization_3)}")
# t4 = time.perf_counter()
# print(
#     f"Voxelization: {(t1-t0)*1000:.3f}ms, 1/2: {(t2-t1)*1000:.3f}ms, 1/3: {(t3-t2)*1000:.3f}ms, 2/3: {(t4-t3)*1000:.3f}ms"
# )

print("\nOctreeBasedVoxelization")
t0 = time.perf_counter()
voxelization_1 = OctreeBasedVoxelization.from_volume_model(volume_model_1, VOXEL_SIZE)
voxelization_2 = OctreeBasedVoxelization.from_volume_model(volume_model_2, VOXEL_SIZE)
voxelization_3 = OctreeBasedVoxelization.from_volume_model(volume_model_3, VOXEL_SIZE)
t1 = time.perf_counter()
print(f"Test 1/2: {voxelization_1.is_intersecting(voxelization_2)}")
t2 = time.perf_counter()
print(f"Test 1/3: {voxelization_1.is_intersecting(voxelization_3)}")
t3 = time.perf_counter()
print(f"Test 2/3: {voxelization_2.is_intersecting(voxelization_3)}")
t4 = time.perf_counter()
print(
    f"Voxelization: {(t1-t0)*1000:.3f}ms, 1/2: {(t2-t1)*1000:.3f}ms, 1/3: {(t3-t2)*1000:.3f}ms, 2/3: {(t4-t3)*1000:.3f}ms"
)

# print("\nNative")
# t1 = time.perf_counter()
# print(f"Test 1/2: {volume_model_1.primitives[0].is_intersecting_with(volume_model_2.primitives[0])}")
# t2 = time.perf_counter()
# print(f"Test 1/3: {voxelization_1.primitives[0].is_intersecting_with(volume_model_3.primitives[0])}")
# t3 = time.perf_counter()
# print(f"Test 2/3: {volume_model_2.primitives[0].is_intersecting_with(volume_model_3.primitives[0])}")
# t4 = time.perf_counter()
# print(
#     f"1/2: {(t2-t1)*1000:.3f}ms, 1/3: {(t3-t2)*1000:.3f}ms, 2/3: {(t4-t3)*1000:.3f}ms"
# )
