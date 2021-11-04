#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 17:46:36 2021

@author: steven
"""

import volmdlr as vm
import volmdlr.faces as vmf

base_points1 = [vm.Point2D(0, 0),
               vm.Point2D(0.1, 0.05),
               vm.Point2D(0.2, 0.1),
               vm.Point2D(0.2, 0.3)]
      
grid1 = []
for z in [0, 0.1, 0.2, 0.3]:
    for x, y in base_points1:
        grid1.append(vm.Point3D(x, y, z))
        

surface1 = vmf.BSplineSurface3D.from_pointgrid(grid1, 4, 4, 2, 2)
face1 = surface1.rectangular_cut(0, 1, 0, 1)


base_points2 = [vm.Point2D(0, 0),
               vm.Point2D(-0.1, 0.05),
               vm.Point2D(-0.2, 0.13),
               vm.Point2D(-0.2, 0.25)]
      
grid2 = []
for z in [0.05, 0.12, 0.21, 0.28]:
    for x, y in base_points2:
        grid2.append(vm.Point3D(x, y, z))
        

surface2 = vmf.BSplineSurface3D.from_pointgrid(grid2, 4, 4, 2, 2)
face2 = surface2.rectangular_cut(0, 1, 0, 1)

shell = vmf.OpenShell3D([face1, face2])
shell.babylonjs()
# ax = surface.plot()
# for p in grid:
#     p.plot(ax=ax, color='r')
# face1.babylonjs()
