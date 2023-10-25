"""
Showcase of the 'TriangleDecimationSimplify' class.
"""
import math
import time

import volmdlr
from volmdlr.cad_simplification import TriangleDecimationSimplify
from volmdlr.core import VolumeModel
from volmdlr.primitives3d import HollowCylinder

VOXEL_SIZE = 0.005

# Create volume model
hollow_cylinder = HollowCylinder.from_center_point_and_axis(volmdlr.O3D, volmdlr.Z3D, 0.1, 0.2, 1.0)
volume_model = VolumeModel(
    [
        hollow_cylinder.rotation(volmdlr.O3D, volmdlr.X3D, -math.pi / 4),
        hollow_cylinder.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 4),
    ]
)

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

volume_model.babylonjs()