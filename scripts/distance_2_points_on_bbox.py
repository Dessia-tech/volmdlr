#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 14:35:43 2019

@author: ringhausen
"""

import volmdlr as vm

bbox = vm.BoundingBox(-1, 1, -1, 1, -1, 1)
point1 = vm.Point3D((-1, 0.5, 0.5))
point2 = vm.Point3D((1, 0.5, 0.5))
d = bbox.distance_between_two_points_on_bbox(point1, point2)
print(d)

point2 = vm.Point3D((-1, 0.5, 0.5))
point1 = vm.Point3D((1, 0.5, 0.5))
d = bbox.distance_between_two_points_on_bbox(point1, point2)
print(d)

point1 = vm.Point3D((0.5, -1, 0.5))
point2 = vm.Point3D((0.5, 1, 0.5))
d = bbox.distance_between_two_points_on_bbox(point1, point2)
print(d)

point2 = vm.Point3D((0.5, -1, 0.5))
point1 = vm.Point3D((0.5, 1, 0.5))
d = bbox.distance_between_two_points_on_bbox(point1, point2)
print(d)

point1 = vm.Point3D((0.5, 0.5, -1))
point2 = vm.Point3D((0.5, 0.5, 1))
d = bbox.distance_between_two_points_on_bbox(point1, point2)
print(d)

point2 = vm.Point3D((0.5, 0.5, -1))
point1 = vm.Point3D((0.5, 0.5, 1))
d = bbox.distance_between_two_points_on_bbox(point1, point2)
print(d)

ax = bbox.plot()
point1.MPLPlot(ax)
point2.MPLPlot(ax)