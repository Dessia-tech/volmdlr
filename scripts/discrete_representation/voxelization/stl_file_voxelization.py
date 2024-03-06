"""
Example of voxelization from an STL file.
"""
from volmdlr.discrete_representation import MatrixBasedVoxelization
from volmdlr.display import Mesh3D
from volmdlr.core import VolumeModel

VOXEL_SIZE = 0.0015
STL_MODEL_FILE_PATH = "../../stl/simple.stl"

# Load and convert the STL
volume_model = VolumeModel([Mesh3D.from_stl_file(STL_MODEL_FILE_PATH)])

# Voxelize the model
voxelization = MatrixBasedVoxelization.from_volume_model(volume_model, VOXEL_SIZE, name="Voxelization")

# Display the result
voxelization_primitive = voxelization.to_mesh().split_shared_vertices()
voxelization_primitive.alpha = 0.5
voxelization_primitive.color = (1, 0, 0)

volume_model.primitives.append(voxelization_primitive)
volume_model.babylonjs()
