import math
import time
from itertools import product

import matplotlib.pyplot as plt

import dessia_common.core
import volmdlr
from volmdlr import edges, primitives3d, wires, cloud, cad_simplification
from volmdlr.core import EdgeStyle

obj = dessia_common.core.DessiaObject.load_from_file('test_cad_simplification.json')

octre_block_simplify = cad_simplification.OctreeBlockSimplify(obj)
time_s = time.time()
simplified_shell = octre_block_simplify.simplify(3)
simplified_shell.alpha = 0.6
volmdlr.core.VolumeModel([obj, simplified_shell]).babylonjs()
simplified_model2 = cad_simplification.TrippleExtrusionSimplify(volmdlr.core.VolumeModel([obj])).simplify()
simplified_model2.alpha = 0.6
volmdlr.core.VolumeModel([obj, simplified_model2]).babylonjs()
time_e = time.time()
print(f'took {time_e - time_s} seconds')
