"""
Showcase of the 'OctreBlockSimplify'
"""
import time

from volmdlr.shells import ClosedShell3D
from volmdlr.cad_simplification import OctreeBlockSimplify
from volmdlr.core import VolumeModel

closed_shell = ClosedShell3D.load_from_file("model.json")

start = time.perf_counter()
octre_block_simplify = OctreeBlockSimplify(closed_shell=closed_shell)
simplified_closed_shell = octre_block_simplify.simplify(depth=3)

print(f"Simplification took {time.perf_counter() - start:.6f} seconds")

simplified_closed_shell.alpha = 0.6
volume_model = VolumeModel([closed_shell, simplified_closed_shell])
volume_model.babylonjs()
