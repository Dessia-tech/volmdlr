#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
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

extrusion_vector = vm.X3D*0.1

for i in range(1, number_holes):
    circles.append(first_circle.rotation(vm.O2D, i*delta_angle))


extrusion = ExtrudedProfile(vm.O3D, vm.Y3D, vm.Z3D, outer_circle,
                                            circles, extrusion_vector)

inner_circles_area = sum([c.area() for c in circles])
assert math.isclose(extrusion.volume(), (outer_circle.area() - inner_circles_area)*extrusion_vector.norm())

model = vmc.VolumeModel([extrusion])
model.babylonjs()
model._check_platform()