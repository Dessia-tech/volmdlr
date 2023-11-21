"""
Compare memory usage of differents voxelization.
"""
import math
import time

import volmdlr
from pympler.asizeof import asizeof
from volmdlr.discrete_representation import MatrixBasedVoxelization, OctreeBasedVoxelization, PointBasedVoxelization
from volmdlr.step import Step

STEP_MODEL_FILE_PATH = "../../step/tore1.step"
VOXEL_SIZE = 0.00001
TRANSLATION_VECTOR = volmdlr.Vector3D(-0.0015, 0.005, 0.01)
ROTATION_ANGLE = math.pi / 2

# Define the volume model from file
volume_model = Step.from_file(STEP_MODEL_FILE_PATH).to_volume_model()

matrix_based = MatrixBasedVoxelization.from_volume_model(volume_model, VOXEL_SIZE)
print(asizeof(matrix_based) + matrix_based.matrix.nbytes)

# point_based = PointBasedVoxelization.from_volume_model(volume_model, VOXEL_SIZE)
# print(asizeof(point_based))

octree_based = OctreeBasedVoxelization.from_volume_model(volume_model, VOXEL_SIZE)
print(asizeof(octree_based))
