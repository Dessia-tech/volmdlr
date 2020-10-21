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
ax = a.MPLPlot()

for p in a.tessellation_points():
    p.MPLPlot(ax=ax)
    
s.MPLPlot(ax=ax, color='r')
e.MPLPlot(ax=ax, color='g')
i.MPLPlot(ax=ax, color='b')


arc1 = vm.Arc3D(vm.Point3D([-0.030962035803739997, 0.0011626900994054661, -0.02]),
                vm.Point3D([-0.031209642286239472, -0.00040063570451895954, -0.02]),
                vm.Point3D([-0.026119083, 0.0, -0.02]),
                vm.Vector3D([0.0, 0.0, 0.001]))

ax = arc1.MPLPlot()
for p in arc1.tessellation_points():
    p.MPLPlot(ax=ax)


arc1.start.MPLPlot(ax=ax, color='r')
arc1.end.MPLPlot(ax=ax, color='m')
arc1.interior.MPLPlot(ax=ax, color='b')
arc1.center.MPLPlot(ax=ax, color='g')


print(arc1.center)
print(arc1.center-vm.Point3D([-0.030962035803739997, 0.0011626900994054661, -0.02]))
print(arc1.center-vm.Point3D([-0.031209642286239472, -0.00040063570451895954, -0.02]))
print(arc1.center-vm.Point3D([-0.026119083, 0.0, -0.02]))