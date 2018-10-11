#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 11 15:16:33 2018

@author: steven
"""

import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import numpy as npy
import random
p1 = vm.Point3D((0, 0, 0))
p2 = vm.Point3D((-0.150, 0, 0))
p3 = vm.Point3D((-0.150, 0.215, 0))
p4 = vm.Point3D((-0.150, 0.215, -0.058))
p5 = vm.Point3D((-0.175, 0.186, -0.042))

points = [p1, p2, p3, p4, p5]
radius = {1: 0.015, 2: 0.020, 3: 0.005}

current_point = p5.vector
#points = [p1, p2]
#radius = {1: 0.010}
for i in range(14):
    current_point += 0.300 * (npy.random.random(3) -0.5)
    points.append(vm.Point3D(current_point))
    radius[4+i] = 0.01 + 0.03 * random.random()
#print(radius)
c = vm.Circle3D(p1, 0.008, p2-p1)

rl = primitives3D.RoundedLineSegments3D(points, radius, closed=False, adapt_radius=True, name='wire')
contour = vm.Contour3D([c])

sweep = primitives3D.Sweep(contour, rl, name = 'Random pipe')

m = vm.VolumeModel([('Random Pipe', [sweep])])

m.FreeCADExport('sweep')

