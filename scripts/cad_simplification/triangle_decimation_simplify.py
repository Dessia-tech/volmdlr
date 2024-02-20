"""
Showcase of the 'TriangleDecimationSimplify' class.
"""
import time

from volmdlr.cad_simplification import TriangleDecimationSimplify
from volmdlr.step import Step

# Load
volume_model = Step.from_file("../step/engine.step").to_volume_model()

# Simplify
start = time.perf_counter()
simplifier = TriangleDecimationSimplify(volume_model=volume_model)
simplified_volume_model = simplifier.simplify(target_ratio=0.6, preserve_border=True, preserve_shells=False)

print(f"Simplification took {time.perf_counter() - start:.6f} seconds\n")

# Count faces triangles:
n_faces = sum(len(shell.faces) for shell in volume_model.get_shells())
n_triangles = sum(len(shell.triangulation().triangles) for shell in volume_model.get_shells())
print(f"Given model has {n_faces} faces and {n_triangles} triangles when triangulated for display.")

n_triangles = sum(len(shell.triangles) for shell in simplified_volume_model.get_shells())
print(f"Simplified model mesh has {n_triangles} triangles.")


# Display
for simplified_shell in simplified_volume_model.get_shells():
    simplified_shell.color = (1.0, 0.0, 0.0)
    simplified_shell.alpha = 0.6
    volume_model.primitives.append(simplified_shell)

volume_model.babylonjs()
