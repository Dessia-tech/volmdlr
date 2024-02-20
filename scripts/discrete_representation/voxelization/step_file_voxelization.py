"""
Example of voxelization from a STEP file.
"""
from volmdlr.discrete_representation import MatrixBasedVoxelization
from volmdlr.step import Step

VOXEL_SIZE = 0.001
STEP_MODEL_FILE_PATH = "../../step/tore1.step"

# Load and convert the SETP
volume_model = Step.from_file(STEP_MODEL_FILE_PATH).to_volume_model()

# Voxelize the model
voxelization = MatrixBasedVoxelization.from_volume_model(volume_model, VOXEL_SIZE)

# Display the result
voxelization_primitive = voxelization.to_mesh().split_shared_vertices()
voxelization_primitive.alpha = 0.5
voxelization_primitive.color = (1, 0, 0)

volume_model.primitives.append(voxelization_primitive)
volume_model.babylonjs()
