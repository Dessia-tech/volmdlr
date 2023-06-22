"""
Example of voxelization from a STL file.
"""
from volmdlr.voxelization import Voxelization
from volmdlr.stl import Stl

VOXEL_SIZE = 0.01

# Load and convert the STL
volume_model = Stl.load_from_file("model.stl").to_volume_model()

# Voxelize the model
voxelization = Voxelization.from_volume_model(volume_model, VOXEL_SIZE, method="iterative")

# Display the result
volume_model.primitives.append(voxelization.to_closed_shell())
volume_model.babylonjs()
