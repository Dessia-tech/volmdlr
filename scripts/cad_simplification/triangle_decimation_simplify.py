"""
Showcase of the 'TriangleDecimationSimplify'
"""
import time

from volmdlr.cad_simplification import TriangleDecimationSimplify
from volmdlr.step import Step

VOXEL_SIZE = 0.005

# Load
volume_model = Step.from_file("engine.step").to_volume_model()

# Simplify
start = time.perf_counter()
octre_block_simplify = TriangleDecimationSimplify(volume_model=volume_model)
simplified_volume_model = octre_block_simplify.simplify(target_ratio=0.2, preserve_border=True)

print(f"Simplification took {time.perf_counter() - start:.6f} seconds\n")

# Count faces triangles:
n_faces = sum(len(shell.faces) for shell in volume_model.get_shells())
n_triangles = sum(len(shell.triangulation().faces) for shell in volume_model.get_shells())
print(f"Given model has {n_faces} faces and {n_triangles} triangles when triangulated for display.")

n_faces = sum(len(shell.faces) for shell in simplified_volume_model.get_shells())
n_triangles = sum(len(shell.triangulation().faces) for shell in simplified_volume_model.get_shells())
print(f"Simplified model has {n_faces} faces and {n_triangles} triangles when triangulated for display.")

# Display
for simplified_shell in simplified_volume_model.get_shells():
    simplified_shell.color = (1.0, 0.0, 0.0)
    simplified_shell.alpha = 0.6
    volume_model.primitives.append(simplified_shell)

volume_model.babylonjs(use_cdn=False)
