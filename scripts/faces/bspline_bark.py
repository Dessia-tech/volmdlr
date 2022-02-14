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
        
size_u, size_v, degree_u, degree_v = 4, 4, 2, 2
surface1 = vm.faces.BSplineSurface3D.points_fitting_into_bspline_surface(grid1, size_u, size_v, degree_u, degree_v)
# surface1 = vmf.BSplineSurface3D.from_pointgrid(grid1, 4, 4, 2, 2)

p1 = surface1.point2d_to_3d(vm.Point2D(0, 0))
p2 = surface1.point2d_to_3d(vm.Point2D(0.5, 0.5))

assert p1.point_distance(p2) < surface1.geodesic_distance(p1, p2)

face1 = surface1.rectangular_cut(0, 1, 0, 1)


base_points2 = [vm.Point2D(0, 0),
               vm.Point2D(-0.1, 0.05),
               vm.Point2D(-0.2, 0.13),
               vm.Point2D(-0.2, 0.25)]
      
grid2 = []
for z in [0.05, 0.12, 0.21, 0.28]:
    for x, y in base_points2:
        grid2.append(vm.Point3D(x, y, z))
        
size_u, size_v, degree_u, degree_v = 4, 4, 2, 2
surface2 = vm.faces.BSplineSurface3D.points_fitting_into_bspline_surface(grid2, size_u, size_v, degree_u, degree_v)
# surface2 = vmf.BSplineSurface3D.from_pointgrid(grid2, 4, 4, 2, 2)
face2 = surface2.rectangular_cut(0, 1, 0, 1)

shell = vmf.OpenShell3D([face1, face2])
shell._check_platform()
shell.babylonjs()
# ax = surface.plot()
# for p in grid:
#     p.plot(ax=ax, color='r')
# face1.babylonjs()
