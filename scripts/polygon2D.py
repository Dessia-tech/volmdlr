#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 17:02:36 2017

@author: steven
"""

import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import time
p1=vm.Point2D((0, 0))
p2=vm.Point2D((1, 0))
p3=vm.Point2D((2, 1))
p4=vm.Point2D((1, 0.5))
p5=vm.Point2D((-1, 0.1))

polygon = vm.Polygon2D([p1, p2, p3, p4, p5])

#print(c.Area())
#print(c.SecondMomentArea(p1))

import numpy as npy

points_inside=[]
points_outside=[]

for i in range(100):
    pt=vm.Point2D(2*npy.random.random(2) - 0.3)
#    print(p.PointDistance(pt))
    if polygon.PointBelongs(pt):
        points_inside.append(pt)
    else:
        points_outside.append(pt)
 
    
#polygon.MPLPlot()
#points_inside.MPLPlot()
#c1=vm.CompositePrimitive2D([polygon, *points_inside])
#c1.MPLPlot()
f, a = polygon.MPLPlot()
for point in points_inside:
    point.MPLPlot(a, style = 'ob')
for point in points_outside:
    point.MPLPlot(a, style = 'or')
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
    pt=vm.Point2D(2*npy.random.random(2) - 0.3)
#    print(p.PointDistance(pt))
    polygon.PointBelongs(pt)
t= time.time() - t 
print('time spent: {}s, {}s/eval'.format(t, t/n))   