#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  3 09:56:29 2017

@author: steven
"""

import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

p1=vm.Point2D((0,0))
p2=vm.Point2D((0.1,0.))
p3=vm.Point2D((0.1,0.2))
p4=vm.Point2D((0,0.1))

p5=vm.Point2D((0.1,0.3))

l1=primitives2D.RoundedLines2D([p1,p2,p3,p4],{0:0.01,2:0.01})
#l2=vm.Circle2D(p5,0.01)
c1=vm.Contour2D([l1])

po=vm.Point3D((0,0,0))
xp=vm.Vector3D((1,0,0))
yp=vm.Vector3D((0,1,0))
c1.MPLPlot()

profile=primitives3D.ExtrudedProfile(po,xp,yp,[c1],(0,0,0.2))

model=vm.VolumeModel([profile])

#profile.MPLPlot((0,0,0),(1,0,0),(0,1,0))

#model.MPLPlot()

model.FreeCADExport('python','extrusion','/usr/lib/freecad/lib/')