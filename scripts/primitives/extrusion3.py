#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math

import volmdlr

import volmdlr as vm
import volmdlr.core as vmc
import volmdlr.edges as vme
import volmdlr.wires as vmw
from volmdlr.primitives3d import ExtrudedProfile

number_holes = 5

outer_circle = vmw.Circle2D(vm.O2D, 0.06)


circles = []
delta_angle = 2*math.pi/number_holes
inner_circle = vmw.Circle2D(vm.O2D, 0.04)
first_circle = vmw.Circle2D(vm.Point2D(0, 0.05), 0.005)
circles = [inner_circle, first_circle]
for i in range(1, number_holes):
    circles.append(first_circle.rotation(vm.O2D, i*delta_angle))

frame = volmdlr.Frame3D(vm.O3D, vm.Y3D, vm.Z3D, vm.X3D)

extrusion = ExtrudedProfile(frame, outer_circle, circles, length=0.1)

print(extrusion.volume())
model = vmc.VolumeModel([extrusion])
model.babylonjs()