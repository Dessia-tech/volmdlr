"""
Showcase of the 'OctreBlockSimplify'
"""
import time

from volmdlr.cad_simplification import OctreeBlockSimplify
from volmdlr.step import Step

# Load
volume_model = Step.from_file("engine.step").to_volume_model()

# Simplify
start = time.perf_counter()
octre_block_simplify = OctreeBlockSimplify(volume_model=volume_model)
simplified_volume_model = octre_block_simplify.simplify(depth=3)

print(f"Simplification took {time.perf_counter() - start:.6f} seconds")

# Display
simplified_closed_shell = simplified_volume_model.get_shells()[0]
simplified_closed_shell.color = (1.0, 0.0, 0.0)
simplified_closed_shell.alpha = 0.6

volume_model.primitives.append(simplified_closed_shell)
volume_model.babylonjs()
