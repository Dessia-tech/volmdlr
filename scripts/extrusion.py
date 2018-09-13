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

p1=vm.Point2D((0, 0))
p2=vm.Point2D((0.1, 0.))
p3=vm.Point2D((0.1, 0.2))
p4=vm.Point2D((0.05, 0.1))
p5=vm.Point2D((0.,0.21))
p6=vm.Point2D((0.05, 0.05))


radius = {0: 0.01, 2: 0.01, 3: 0.015}

outer_profile = primitives2D.RoundedLineSegments2D([p1, p2, p3, p4, p5], radius, closed = True)
hole = vm.Circle2D(p6, 0.01)
c1 = vm.Contour2D([outer_profile])
c2 = vm.Contour2D([hole])

po=vm.Point3D((0, 0, 0))
xp=vm.Vector3D((1, 0, 0))
yp=vm.Vector3D((0, 1, 0))
c1.MPLPlot()
c2.MPLPlot()

profile=primitives3D.ExtrudedProfile(po, xp, yp, c1, [c2], (0,0, 0.2), name = 'extrusion')

model=vm.VolumeModel([('profile', [profile])])

#profile.MPLPlot((0,0,0),(1,0,0),(0,1,0))

#model.MPLPlot()

model.FreeCADExport('extrusion')