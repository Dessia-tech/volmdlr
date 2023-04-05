#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 13:03:54 2023

@author: masfaraud
"""

import math
import volmdlr as vm
from volmdlr import core as vmc
from volmdlr import wires as vmw
from volmdlr import primitives2d as p2d
from volmdlr import primitives3d as p3d

GAP = 8
SNAKE_WIDTH = 61
TAIL_OFFSET = 72
EYE_RADIUS = 12

radius = 0.7*SNAKE_WIDTH

center = vm.Point3D(0., -SNAKE_WIDTH - 0.5*GAP, 0)
eye_center = vm.Point2D(-36, 36)

p1 = vm.Point2D(0., 0.)
p2 = vm.Point2D(0., GAP)
p3 = vm.Point2D(-SNAKE_WIDTH, GAP)
p4 = vm.Point2D(-SNAKE_WIDTH, GAP+SNAKE_WIDTH)
p5 = vm.Point2D(SNAKE_WIDTH, GAP+SNAKE_WIDTH)
p6 = vm.Point2D(SNAKE_WIDTH, -SNAKE_WIDTH)
p7 = vm.Point2D(-TAIL_OFFSET, -SNAKE_WIDTH)
p8 = vm.Point2D(-TAIL_OFFSET, -SNAKE_WIDTH - TAIL_OFFSET)
p9 = vm.Point2D(-TAIL_OFFSET-SNAKE_WIDTH, -SNAKE_WIDTH - TAIL_OFFSET)
p10 = vm.Point2D(-TAIL_OFFSET-SNAKE_WIDTH, 0.)

snake1_eye = vmw.Circle2D(eye_center, EYE_RADIUS)

snake1_contour = p2d.ClosedRoundedLineSegments2D([p1, p2, p3, p4, p5, p6, p7, p8, p9, p10],
                                                 {3: radius, 4: radius, 5: radius, 6: radius, 8: radius, 9: radius})

snake1_contour.plot()

snake1 = p3d.ExtrudedProfile(vm.O3D, vm.X3D, vm.Y3D, snake1_contour, [snake1_eye], GAP*vm.Z3D,
                             color=(0.188, 0.412, 0.596), name='Blue snake')
snake2 = snake1.rotation(center, vm.Z3D, math.pi)
snake2.color = (1, 0.867, 0.329)
snake2.name = 'Yellow snake'

model = vmc.VolumeModel([snake1, snake2])
model.babylonjs()
