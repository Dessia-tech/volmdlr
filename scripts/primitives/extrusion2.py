#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 22:01:35 2017

@author: steven


"""

import math

#import numpy as npy
import volmdlr as vm
import volmdlr.core
import volmdlr.core as vmc
import volmdlr.edges as vme
import volmdlr.primitives2d as primitives2d
import volmdlr.primitives3d as primitives3d
import volmdlr.wires as vmw

p1=vm.Point2D(0, 0)
p2=vm.Point2D(0.1, 0.)
p3=vm.Point2D(0.1, 0.2)
p4=vm.Point2D(0, 0.1)
p5=vm.Point2D(-0.01, 0.05)

#p6=vm.Point2D((0.1,0.3))

l1 = primitives2d.OpenedRoundedLineSegments2D([p1, p2, p3, p4], {2: 0.01})
l2 = vme.Arc2D.from_3_points(p4, p5, p1)
c1 = vmw.Contour2D([l1, l2])
c2 = c1.rotation(vm.Point2D(0,0), math.pi)
ax = c1.plot()
c2.plot(ax=ax, edge_style=volmdlr.core.EdgeStyle(color='r'))
#c3 = vm.Contour2D([c1, c2])
#c3.MPLPlot()




profile = primitives3d.ExtrudedProfile(vm.OYZX, c1, [], 0.1)

model = vmc.VolumeModel([profile])
model.babylonjs()

#profile.MPLPlot((0,0,0),(1,0,0),(0,1,0))

#model.MPLPlot()

#model.FreeCADExport('extrusion2')
