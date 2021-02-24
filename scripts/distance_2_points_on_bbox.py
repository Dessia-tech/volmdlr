#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:35:43 2019

@author: ringhausen
"""

import volmdlr as vm
import volmdlr.primitives3d as p3d

#bbox = vm.BoundingBox(-1, 1, -1, 1, -1, 1)
##point1 = vm.Point3D((-0.8, 1, -0.8))
##point2 = vm.Point3D((1, 0.5, 0.5))
#point1 = vm.Point3D((-1, 0.8, 0.8))
#point2 = vm.Point3D((1, -0.7, -0.7))
#
#d = bbox.distance_between_two_points_on_bbox(point1, point2)
#print(d)


origin = vm.Point3D((-0.5,0.5,0))
u = vm.Vector3D((0.7,0,0))
v = vm.Vector3D((0,0.7,0))
w = vm.Vector3D((0,0,0.7))
frame = vm.Frame3D(origin, u, v, w)
block_obstacle1 = p3d.Block(frame, 'test', (0.5,0,0))
block_obstacle1.Translation(vm.Vector3D((-5,2,3)), False)


bbox = block_obstacle1.bounding_box
point1 = vm.Point3D((-0.85, 0.44736842, 0.)) 
point2 = vm.Point3D((-0.285, 0.15, 0.))
point1.Translation(vm.Vector3D((-5,2,3)), False)
point2.Translation(vm.Vector3D((-5,2,3)), False)

d = bbox.distance_between_two_points_on_bbox(point1, point2)
print(d)

#point2 = vm.Point3D((-1, 0.5, 0.5))
#point1 = vm.Point3D((1, 0.5, 0.5))
#d = bbox.distance_between_two_points_on_bbox(point1, point2)
#print(d)
#
#point1 = vm.Point3D((0.5, -1, 0.5))
#point2 = vm.Point3D((0.5, 1, 0.5))
#d = bbox.distance_between_two_points_on_bbox(point1, point2)
#print(d)
#
#point2 = vm.Point3D((0.5, -1, 0.5))
#point1 = vm.Point3D((0.5, 1, 0.5))
#d = bbox.distance_between_two_points_on_bbox(point1, point2)
#print(d)
#
#point1 = vm.Point3D((0.5, 0.5, -1))
#point2 = vm.Point3D((0.5, 0.5, 1))
#d = bbox.distance_between_two_points_on_bbox(point1, point2)
#print(d)
#
#point2 = vm.Point3D((0.5, 0.5, -1))
#point1 = vm.Point3D((0.5, 0.5, 1))
#d = bbox.distance_between_two_points_on_bbox(point1, point2)
#print(d)

ax = bbox.plot()
point1.MPLPlot(ax)
point2.MPLPlot(ax)

primitives = [bbox, point1, point2]
primitives.extend(d)
volumemodel = vm.VolumeModel(primitives)
volumemodel.BabylonShow()
