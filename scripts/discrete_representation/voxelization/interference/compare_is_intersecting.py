"""
Performance comparaison on 'ìs_intersecting' method.
"""
import math
import time

import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.discrete_representation import PointBasedVoxelization
from volmdlr.step import Step

STEP_MODEL_FILE_PATH = "../../../step/tore1.step"
VOXEL_SIZE = 0.0001
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
volume_model.babylonjs()
