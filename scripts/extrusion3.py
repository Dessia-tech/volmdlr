#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
import volmdlr as vm
from volmdlr.primitives3D import ExtrudedProfile

number_holes = 5

outer_circle = vm.Circle2D(vm.O2D, 0.06)


circles = []
delta_angle = 2*math.pi/number_holes
inner_circle = vm.Circle2D(vm.O2D, 0.04)
first_circle = vm.Circle2D(vm.Point2D((0, 0.05)), 0.005)
circles = [inner_circle, first_circle]
for i in range(1, number_holes):
    circles.append(first_circle.Rotation(vm.O2D, i*delta_angle))


extrusion = ExtrudedProfile(vm.O3D, vm.Y3D, vm.Z3D, outer_circle,
                                            circles, vm.X3D*0.1)

model = vm.VolumeModel([extrusion])
model.babylonjs()