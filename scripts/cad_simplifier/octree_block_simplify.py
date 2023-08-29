"""
Showcase of the 'OctreBlockSimplify'
"""
import time

from volmdlr.cad_simplification import OctreeBlockSimplify
from volmdlr.step import Step
from volmdlr.core import VolumeModel

# Load
volume_model = Step.from_file("engine.step").to_volume_model()
closed_shell = volume_model.get_shells()[0]

# Simplify
start = time.perf_counter()
octre_block_simplify = OctreeBlockSimplify(closed_shell=closed_shell)
simplified_closed_shell = octre_block_simplify.simplify(depth=3)

print(f"Simplification took {time.perf_counter() - start:.6f} seconds")

# Display
simplified_closed_shell.color = (1.0, 0.0, 0.0)
simplified_closed_shell.alpha = 0.6

volume_model = VolumeModel([closed_shell, simplified_closed_shell])
volume_model.babylonjs()
