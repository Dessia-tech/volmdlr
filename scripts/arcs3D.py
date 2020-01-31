#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 18:06:45 2018

@author: steven
"""

import volmdlr as vm




i = vm.Point3D((1,0,0))
e = i.Rotation(vm.O3D, vm.Z3D, 1)
s = i.Rotation(vm.O3D, vm.Z3D, -3.5)

a = vm.Arc3D(s, i, e)
assert a.angle == 4.5


# Random arc
i = vm.Point3D.random(-1,1,-1,1,-1,1)
e = vm.Point3D.random(-1,1,-1,1,-1,1)
s = vm.Point3D.random(-1,1,-1,1,-1,1)

a = vm.Arc3D(s, i, e)
fig, ax = a.MPLPlot()

for p in a.tessellation_points():
    p.MPLPlot(ax=ax)
    
s.MPLPlot(ax=ax, color='r')
e.MPLPlot(ax=ax, color='g')
i.MPLPlot(ax=ax, color='b')