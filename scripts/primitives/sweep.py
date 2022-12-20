#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 15:16:33 2018

@author: steven
"""

import random
import matplotlib.pyplot as plt

import volmdlr as vm
import volmdlr.primitives3d as primitives3d
import volmdlr.wires as wires

p1 = vm.Point3D(0, 0, 0)
p2 = vm.Point3D(-0.150, 0, 0)
p3 = vm.Point3D(-0.150, 0.215, 0)
p4 = vm.Point3D(-0.150, 0.215, -0.058)
p5 = vm.Point3D(-0.175, 0.186, -0.042)

points = [p1, p2, p3, p4, p5]
radius = {1: 0.015, 2: 0.020, 3: 0.005}

current_point = p5

for i in range(6):
    current_point += vm.Point3D.random(-0.1, 0.3, -0.1, 0.3, -0.1, 0.3)
    points.append(current_point)
    radius[4+i] = 0.01 + 0.03 * random.random()

contour = wires.Circle2D(vm.O2D, 0.008)

rl = primitives3d.OpenRoundedLineSegments3D(points, radius, adapt_radius=True, name='wire')

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for prim in rl.primitives:
    prim.plot(ax=ax)


r1 = rl.to_dict()
r2 = primitives3d.OpenRoundedLineSegments3D.dict_to_object(r1)
c1 = contour.to_dict()
c2 = vm.wires.Circle2D.dict_to_object(c1)

sweep = primitives3d.Sweep(contour, rl, name='Random pipe')

model = vm.core.VolumeModel([sweep])
model.babylonjs()


model.to_step('sweep.step')
