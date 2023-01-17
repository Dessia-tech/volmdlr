#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 22 14:20:01 2021

@author: dasilva
"""


import math
import time

import matplotlib.pyplot as plt

import volmdlr as vm
import volmdlr.primitives2d as primitives2d
import volmdlr.wires as vmw

n = 1000

points = []
for i in range(n):
    points.append(vm.Point2D.random(-0.3, 0.7, -0.3, 0.7))
    
vec1 = vm.Point2D(1, 0.5)
vec2 = vm.Point2D(0.5, -1)
new_points = []
for point in points:
    new_points.append(point.translation(vec1))
    new_points.append(point.translation(vec2))
    
# fig = plt.figure()
# ax = fig.add_subplot(111)
# for pt in points + new_points:
#     pt.plot(ax=ax)
# polygon = vm.wires.ClosedPolygon2D.points_convex_hull(points+new_points)
# polygon = vm.wires.ClosedPolygon2D.convex_hull_points(points+new_points)
# polygon = vm.wires.ClosedPolygon2D.hull(points+new_points, 0.06)
polygon = vm.wires.ClosedPolygon2D.convex_hull(points+new_points)
# polygon,nearby_points = vm.wires.ClosedPolygon2D.concave_hull(new_points, concavity=-1, scale_factor=0.005)

fig = plt.figure()
ax = fig.add_subplot(111)
for pt in points + new_points:
    pt.plot(ax=ax)
for point in polygon.points:
    point.plot(ax= ax, color = 'g')
for line in polygon.line_segments:
    line.plot(ax= ax, color = 'r')
# for point in nearby_points:
#     point.plot(ax=ax, color = 'r')
