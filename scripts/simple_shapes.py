#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 14 15:32:37 2018

@author: Steven Masfaraud masfaraud@dessia.tech
"""
import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import numpy as npy

#for i in range(20):
triangle_points=[vm.Point2D(npy.random.random(2)) for i in range(3)]
triangle=vm.Polygon2D(triangle_points)


cog_triangle=triangle.CenterOfMass()
c1 = vm.CompositePrimitive2D([triangle, cog_triangle])
c1.MPLPlot()

print(triangle.Area())

p0=vm.Point2D((-0.5,0))
p1=vm.Point2D((0,0.5))
p2=vm.Point2D((0.5,0))

a=vm.Arc2D(p0,p1,p2)
l=vm.LineSegment2D(p2,p0)

c=vm.Contour2D([a,l])

c2 = vm.CompositePrimitive2D([c])
c2.MPLPlot()
print(c.Area())