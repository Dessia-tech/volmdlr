#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 17:02:36 2017

@author: steven
"""

import volmdlr as vm
import volmdlr.wires as wires
import volmdlr.primitives2d as primitives2d
import time
p1 = vm.Point2D(0, 0)
p2 = vm.Point2D(1, 0)
p3 = vm.Point2D(2, 1)
p4 = vm.Point2D(1, 0.5)
p5 = vm.Point2D(-0.5, 1)

polygon = vm.wires.ClosedPolygon2D([p1, p2, p3, p4, p5])
a = polygon.self_intersects()
print(a)
#print(c.Area())
#print(c.SecondMomentArea(p1))

import numpy as npy
npy.seterr(divide='raise')

points_inside=[]
points_outside=[]

for i in range(100):
    pt=vm.Point2D.random(-0.5, 2, -0.1, 1.1)
#    print(p.PointDistance(pt))
    if polygon.point_belongs(pt):
        points_inside.append(pt)
    else:
        points_outside.append(pt)
 
    
#polygon.MPLPlot()
#points_inside.MPLPlot()
#c1=vm.CompositePrimitive2D([polygon, *points_inside])
#c1.MPLPlot()
a = polygon.plot()
for point in points_inside:
    point.plot(a, color='b')
for point in points_outside:
    point.plot(a, color = 'r')
#
#c2=vm.CompositePrimitive2D([polygon, *points_outside])
#c2.MPLPlot()
#
#cog_p = polygon.CenterOfMass()
#c3 = vm.CompositePrimitive2D([polygon, cog_p])
#c3.MPLPlot()

# Speed test
t = time.time()
n = 100000
for i in range(n):
    pt=vm.Point2D.random(-0.3, 0.7, -0.3, 0.7)
#    print(p.PointDistance(pt))
    polygon.point_belongs(pt)
t= time.time() - t 
print('time spent: {}s, {}s/eval'.format(t, t/n))   