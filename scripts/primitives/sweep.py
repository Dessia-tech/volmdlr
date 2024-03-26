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
r1 = open_rounded_line_segements.to_dict()
r2 = primitives3d.OpenRoundedLineSegments3D.dict_to_object(r1)
c1 = contour.to_dict()
c2 = curves.Circle2D.dict_to_object(c1)
origin = open_rounded_line_segements.primitives[0].start
w = open_rounded_line_segements.primitives[0].unit_direction_vector(0.)
u = open_rounded_line_segements.primitives[0].unit_normal_vector(0.)
if not u:
    u = w.deterministic_unit_normal_vector()
v = w.cross(u)
frame = volmdlr.Frame3D(origin, u, v, w)
frame.plot(ax, ratio=0.01)
open_rounded_line_segements.primitives[0].start.plot(ax)
for prim in open_rounded_line_segements.primitives:
    prim.plot(ax=ax)
    frame = prim.move_frame_along(frame)
    frame.plot(ax, ratio=0.025)



sweep = primitives3d.Sweep(contour, open_rounded_line_segements, name='Random pipe')

model = vm.model.VolumeModel([sweep])
model._check_platform()
model.babylonjs()

model.to_step('sweep.step')

contour = wires.ClosedPolygon2D([volmdlr.Point2D(-0.008, -0.004), volmdlr.Point2D(0.008, -0.004),
                                 volmdlr.Point2D(0.008, 0.004), volmdlr.Point2D(-0.008, 0.004)])

sweep = primitives3d.Sweep(contour, open_rounded_line_segements, name='Random pipe')
model = vm.model.VolumeModel([sweep])
model._check_platform()
model.babylonjs()

model.to_step('sweep.step')
