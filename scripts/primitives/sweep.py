#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A sweep example.
"""

import random
import matplotlib.pyplot as plt
import volmdlr

import volmdlr as vm
import volmdlr.primitives3d as primitives3d
import volmdlr.wires as wires
from volmdlr import curves
from volmdlr.models.open_rounded_line_segments import open_rounded_line_segements

# contour = wires.Circle2D(vm.O2D, 0.008)
contour = wires.ClosedPolygon2D([volmdlr.Point2D(-0.004, -0.004), volmdlr.Point2D(0.004, -0.004),
                                 volmdlr.Point2D(0.004, 0.004), volmdlr.Point2D(-0.004, 0.004)])

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for prim in open_rounded_line_segements.primitives:
    prim.plot(ax=ax)


r1 = open_rounded_line_segements.to_dict()
r2 = primitives3d.OpenRoundedLineSegments3D.dict_to_object(r1)
c1 = contour.to_dict()
c2 = curves.Circle2D.dict_to_object(c1)

sweep = primitives3d.Sweep(contour, open_rounded_line_segements, name='Random pipe')

model = vm.core.VolumeModel([sweep])
model._check_platform()
model.babylonjs()

model.to_step('sweep.step')
