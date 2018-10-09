#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 22:01:35 2017

@author: steven


"""

#import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
#import math

p1=vm.Point2D((0, 0))
p2=vm.Point2D((0.1, 0.))
p3=vm.Point2D((0.1, 0.2))
p4=vm.Point2D((0, 0.1))
p5=vm.Point2D((-0.01, 0.05))

#p6=vm.Point2D((0.1,0.3))

l1=primitives2D.RoundedLineSegments2D([p1, p2, p3, p4], {2: 0.01}, False)
l2=vm.Arc2D(p4, p5, p1)
c1=vm.Contour2D([l1, l2])
c2=c1.Rotation(vm.Point2D((0,0)),npy.pi,True)
c1.MPLPlot()
c2.MPLPlot()
c3=vm.Contour2D([c1, c2])
c3.MPLPlot()

po = vm.Point3D((0, 0, 0))
xp = vm.Vector3D((1, 0, 0))
yp = vm.Vector3D((0, 1, 0))


profile = primitives3D.ExtrudedProfile(po, xp, yp, c1, [], vm.Vector3D((0, 0.1, 0.2)))

model = vm.VolumeModel([('', [profile])])

#profile.MPLPlot((0,0,0),(1,0,0),(0,1,0))

#model.MPLPlot()

model.FreeCADExport('extrusion2')