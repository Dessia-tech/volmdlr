#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math

import volmdlr as vm
import volmdlr.core as vmc
import volmdlr.edges
import volmdlr.edges as vme
import volmdlr.wires as vmw
from volmdlr.primitives3d import ExtrudedProfile
from volmdlr import curves
number_holes = 5

outer_circle = vmw.Contour2D([volmdlr.edges.FullArc2D.from_curve(curves.Circle2D(vm.O2D, 0.06))])


delta_angle = 2*math.pi/number_holes
inner_circle = vmw.Contour2D([volmdlr.edges.FullArc2D.from_curve(curves.Circle2D(vm.O2D, 0.04))])
first_circle = vmw.Contour2D([volmdlr.edges.FullArc2D.from_curve(curves.Circle2D(vm.Point2D(0, 0.05), 0.005))])
circles = [inner_circle, first_circle]

extrusion_length = 0.1

for i in range(1, number_holes):
    circles.append(first_circle.rotation(vm.O2D, i*delta_angle))


extrusion = ExtrudedProfile(vm.OYZX, outer_circle, circles, extrusion_length)

inner_circles_area = sum([c.area() for c in circles])
assert math.isclose(extrusion.volume(), (outer_circle.area() - inner_circles_area)*extrusion_length)

model = vmc.VolumeModel([extrusion])
model.babylonjs()
model._check_platform()
