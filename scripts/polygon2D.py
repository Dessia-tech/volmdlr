#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 17:02:36 2017

@author: steven
"""

import volmdlr as vm
import volmdlr.primitives2D as primitives2D

p1=vm.Point2D((0,0))
p2=vm.Point2D((1,0))
p3=vm.Point2D((2,1))
p4=vm.Point2D((1,0.5))
p5=vm.Point2D((-1,0.1))

p=vm.Polygon2D([p1,p2,p3,p4,p5])

#print(c.Area())
#print(c.SecondMomentArea(p1))

import numpy as npy

points_inside=[]
points_outside=[]
for i in range(100):
    pt=vm.Point2D(2*npy.random.random(2)-0.3)
#    print(p.PointDistance(pt))
    if p.PointBelongs(pt):
        points_inside.append(pt)
    else:
        points_outside.append(pt)
        
p.MPLPlot()
#points_inside.MPLPlot()
c1=vm.CompositePrimitive2D([p,*points_inside])
c1.MPLPlot()

c2=vm.CompositePrimitive2D([p,*points_outside])
c2.MPLPlot()

cog_p=p.CenterOfMass()
cog_p.MPLPlot('ro')