#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 09:56:29 2017

@author: steven
"""

import numpy as npy
import volmdlr

import volmdlr as vm
import volmdlr.edges as edges
import volmdlr.wires as wires
import volmdlr.primitives3d as primitives3d
import volmdlr.primitives2d as primitives2d

p1 = vm.Point2D(0, 0)
p2 = vm.Point2D(0.1, 0.)
p3 = vm.Point2D(0.1, 0.2)
p4=vm.Point2D(0.05, 0.1)
p5=vm.Point2D(0.,0.21)
p6=vm.Point2D(0.05, 0.05)

p7 = vm.Point2D(0.06, 0.05)
p8 = vm.Point2D(0.04, 0.07)

radius = {0: 0.01, 2: 0.01, 3: 0.015}

outer_profile = primitives2d.ClosedRoundedLineSegments2D([p1, p2, p3, p4, p5], radius)
l1 = edges.LineSegment2D(p6, p7)
l2 = edges.LineSegment2D(p7, p8)
l3 = edges.LineSegment2D(p8, p6)
c2 = wires.Contour2D([l1, l2, l3])

frame = volmdlr.Frame3D(vm.O3D, vm.Y3D, vm.Z3D, vm.X3D)
profile=primitives3d.ExtrudedProfile(frame, outer_profile, [c2], 0.1, name='extrusion')

model=vm.core.VolumeModel([profile])

model.babylonjs()

model.to_step('extrusion')