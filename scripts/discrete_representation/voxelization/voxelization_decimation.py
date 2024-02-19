"""
Example of voxelization decimation.
"""
from volmdlr.discrete_representation import MatrixBasedVoxelization
from volmdlr.stl import Stl

VOXEL_SIZE = 0.0015
STL_MODEL_FILE_PATH = "../../stl/simple.stl"

# Load and convert the STL
volume_model = Stl.load_from_file(STL_MODEL_FILE_PATH).to_volume_model()

# Voxelize the model
voxelization = MatrixBasedVoxelization.from_volume_model(volume_model, VOXEL_SIZE, name="Voxelization")

# Display the result
voxelization_primitive = voxelization.to_mesh()
voxelization_primitive.alpha = 0.5
voxelization_primitive.color = (1, 0, 0)

volume_model.primitives.append(voxelization_primitive)
volume_model.babylonjs()

voxelization_primitive_decimated = voxelization_primitive.decimate(500, preserve_border=True)
voxelization_primitive_decimated.alpha = 0.5
voxelization_primitive_decimated.color = (1, 0, 0)

volume_model.primitives[-1] = voxelization_primitive_decimated
volume_model.babylonjs()

print(f"Number of triangles before decimation: {voxelization_primitive.n_triangles}")
print(f"Number of triangles after decimation: {voxelization_primitive_decimated.n_triangles}")
