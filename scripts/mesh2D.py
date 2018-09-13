#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 12 10:04:09 2017

@author: steven
"""

import volmdlr as vm

l=0.1

p1=vm.Point2D((0,0))
p2=vm.Point2D((l,0))
p3=vm.Point2D((l,l))
p4=vm.Point2D((0,l))



l1=vm.LineSegment2D(p1,p2)
l2=vm.LineSegment2D(p2,p3)
l3=vm.LineSegment2D(p3,p4)
l4=vm.LineSegment2D(p4,p1)

p5=vm.Point2D((l/2,l/2))
c1=vm.Circle2D(p5,l/5)

ct1=vm.Contour2D([l4,l3,l2,l1])
ct2=vm.Contour2D([c1])
mesh=vm.Mesh2D([ct1,ct2],{},0.01)

print(mesh.GeoScript('mesh2D.geo'))