"""
Performance comparaison on 'union' method.
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

print("\nOctreeBasedVoxelization")
t0 = time.perf_counter()
voxelization_1 = OctreeBasedVoxelization.from_volume_model(volume_model_1, VOXEL_SIZE)
voxelization_2 = OctreeBasedVoxelization.from_volume_model(volume_model_2, VOXEL_SIZE)
voxelization_3 = OctreeBasedVoxelization.from_volume_model(volume_model_3, VOXEL_SIZE)
t1 = time.perf_counter()
union_21 = voxelization_1.union(voxelization_2)
t2 = time.perf_counter()
union_22 = voxelization_1.union(voxelization_3)
t3 = time.perf_counter()
union_23 = voxelization_2.union(voxelization_3)
t4 = time.perf_counter()
print(
    f"Voxelization: {(t1-t0)*1000:.3f}ms, 1/2: {(t2-t1)*1000:.3f}ms, 1/3: {(t3-t2)*1000:.3f}ms, 2/3: {(t4-t3)*1000:.3f}ms"
)

print("\nPointBasedVoxelization")
t0 = time.perf_counter()
voxelization_1 = voxelization_1.to_point_based_voxelization()
voxelization_2 = voxelization_2.to_point_based_voxelization()
voxelization_3 = voxelization_3.to_point_based_voxelization()
t1 = time.perf_counter()
union_11 = voxelization_1.union(voxelization_2)
t2 = time.perf_counter()
union_12 = voxelization_1.union(voxelization_3)
t3 = time.perf_counter()
union_13 = voxelization_2.union(voxelization_3)
t4 = time.perf_counter()
print(
    f"Voxelization: {(t1-t0)*1000:.3f}ms, 1/2: {(t2-t1)*1000:.3f}ms, 1/3: {(t3-t2)*1000:.3f}ms, 2/3: {(t4-t3)*1000:.3f}ms"
)


print(union_11 == union_21.to_point_based_voxelization())
print(union_12 == union_22.to_point_based_voxelization())
print(union_13 == union_23.to_point_based_voxelization())
