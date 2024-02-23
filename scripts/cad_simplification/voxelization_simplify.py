"""
Showcase of the 'VoxelizationSimplify' class.
"""
import time

from volmdlr.cad_simplification import VoxelizationSimplify
from volmdlr.step import Step

VOXEL_SIZE = 0.002

# Load
volume_model = Step.from_file("../step/engine.step").to_volume_model()

# Simplify
start = time.perf_counter()
simplifier = VoxelizationSimplify(volume_model=volume_model)
simplified_volume_model = simplifier.simplify(voxel_size=VOXEL_SIZE)

print(f"Simplification took {time.perf_counter() - start:.6f} seconds\n")

# Count faces triangles:
n_faces = sum(len(shell.faces) for shell in volume_model.get_shells())
n_triangles = sum(len(shell.triangulation().triangles) for shell in volume_model.get_shells())
print(f"Given model has {n_faces} faces and {n_triangles} triangles when triangulated for display.")

n_faces = sum(len(shell.faces) for shell in simplified_volume_model.get_shells())
n_triangles = sum(len(shell.triangulation().faces) for shell in simplified_volume_model.get_shells())
print(f"Simplified model has {n_faces} faces and {n_triangles} triangles when triangulated for display.")

# Display
for simplified_shell in simplified_volume_model.get_shells():
    simplified_shell.color = (1.0, 0.0, 0.0)
    simplified_shell.alpha = 0.6
    volume_model.primitives.append(simplified_shell)

volume_model.babylonjs()
