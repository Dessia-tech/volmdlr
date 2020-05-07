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
for i in range(4):
    current_point += 0.300 * (npy.random.random(3) -0.5)
    points.append(vm.Point3D(current_point))
    radius[4+i] = 0.01 + 0.03 * random.random()
#print(radius)
# c = vm.Circle3D(p1, 0.008, p2-p1)
c = vm.Circle2D(vm.Point2D((0,0)), 0.008)

rl = primitives3D.OpenedRoundedLineSegments3D(points, radius, adapt_radius=True, name='wire')
contour = vm.Contour2D([c])

r1 = rl.to_dict()
r2 = primitives3D.OpenedRoundedLineSegments3D.dict_to_object(r1)
c1 = c.to_dict()
c2 = vm.Circle2D.dict_to_object(c1)

c1 = contour.to_dict()
c2 = vm.Contour2D.dict_to_object(c1)


sweep = primitives3D.Sweep(contour, rl, name = 'Random pipe')
sweepy = sweep.copy()
v1 = vm.Vector3D((1,1,1))
v1.Normalize()
v2 = v1.deterministic_unit_normal_vector()
v3 = v1.Cross(v2)
frame0 = vm.Frame3D(vm.Point3D((0,0,0)), v1, v2, v3)

frame_mapping = sweepy.frame_mapping(frame0, 'new', False)

m = vm.VolumeModel([sweepy])
# m = vm.VolumeModel([frame_mapping])
m.babylonjs()
# m.FreeCADExport('sweep')


