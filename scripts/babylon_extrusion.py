#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  6 10:25:22 2019

@author: ringhausen
"""

import volmdlr as vm
import volmdlr.primitives3D as p3d
import volmdlr.primitives2D as p2d
import math

p1=vm.Point2D((0, 0))
p2=vm.Point2D((0, 2))
p3=vm.Point2D((2, 4))
p4=vm.Point2D((4, 4))
p5=vm.Point2D((4, 3))
p6=vm.Point2D((3, 2))
p7=vm.Point2D((3, 0))

l1=p2d.RoundedLineSegments2D([p7, p1, p2], {}, False)
l2=vm.Arc2D(p2, vm.Point2D((math.sqrt(2)/2, 3+math.sqrt(2)/2)), p3)
print('===')
print(l2.is_trigo)
print('===')
l3=p2d.RoundedLineSegments2D([p3, p4, p5, p6], {}, False, adapt_radius=True)
l4=vm.Arc2D(p6, vm.Point2D((4, 1)), p7)
#l4=vm.Arc2D(p7, vm.Point2D((4, 1)), p6)
print('===')
print(l4.is_trigo)
print('===')
c1=vm.Contour2D([l1, l2, l3, l4])

#ax = l1.MPLPlot()
#l2.MPLPlot(ax=ax)
#l3.MPLPlot(ax=ax)
#l4.MPLPlot(ax=ax)

profile = p3d.ExtrudedProfile(vm.o3D, vm.x3D, vm.y3D, c1, [], vm.Vector3D((0,0,1)))

model = vm.VolumeModel([profile])
model.BabylonShow()

#%%

p1=vm.Point2D((0, 0))
p2=vm.Point2D((0, 2))
p3=vm.Point2D((2, 2))
p4=vm.Point2D((2, 0))

l1 = p2d.RoundedLineSegments2D([p1, p2, p3, p4], {}, True)
c1 = vm.Contour2D([l1])

profile = p3d.ExtrudedProfile(vm.o3D, vm.x3D, vm.y3D, c1, [], vm.Vector3D((0,0,1)))

model = vm.VolumeModel([profile])
model.BabylonShow()
