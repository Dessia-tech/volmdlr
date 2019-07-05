#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 17:32:34 2019

@author: ringhausen
"""

import volmdlr as vm
import numpy as npy

#nb = 5
#p0 = vm.Point2D((0,0))
#
#pts = []
#for i in range(nb):    
#    pts.append(vm.Point2D(npy.random.random(2)))
#    
#polygone = vm.Polygon2D(pts)
#polygone = vm.Polygon2D([vm.Point2D((0.,0.)),vm.Point2D((1.,0.)),vm.Point2D((1.,1.)),vm.Point2D((1.,0.)),vm.Point2D((0.,1.))])
##polygone = vm.Polygon2D([vm.Point2D((0.,0.)),vm.Point2D((1.,0.)),vm.Point2D((0.,0.)),vm.Point2D((1.,0.))])
#polygone.Rotation(p0 , 0.5)
#polygone.MPLPlot()
#a,b,c = polygone.SelfIntersect()
#print(a)
#if b is not None and c is not None:
#    print(b.points)
#    print(c.points)


#while polygone.SelfIntersect():
#    pts = []
#    for i in range(nb):    
#        pts.append(vm.Point2D(npy.random.random(2)))
#    polygone = vm.Polygon2D(pts)
#    #    print(polygone.SelfIntersect(), '\n\n')
#    a,b,c = polygone.SelfIntersect()
#    
#    if polygone.SelfIntersect():    
#        polygone.MPLPlot(style='r')
#    else:
#        polygone.MPLPlot()


test = [[0., 0.035], [0., 0.045], [0.03, 0.04], [0.03, 0.035], [0., 0.035], [0., 0.055], [0., 0.065], [0.03, 0.065]]
test = [[0,0], [0,1], [0,0], [1,1]]
test = [[0.,    0.035], [0.,    0.045], [0. ,   0.055], [0. ,   0.065], [0.03,  0.065], [0.03, 0.04], [0.03,  0.035]] #, [0.,    0.035]]
pts = []
for i in test:
    pts.append(vm.Point2D(i))
    
polygon = vm.Polygon2D(pts)
a,b,c = polygon.SelfIntersect()
print(a)
print(b.points,c.points)