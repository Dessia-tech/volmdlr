# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 10:51:54 2020

@author: Mack Pro
"""


import volmdlr as vm
import volmdlr.step
import volmdlr.primitives3d as primitives3d
import math
import matplotlib.pyplot as plt

radius = 5e-3 #Choose the radius
center = vm.Point3D(0,0,0) #Choose the coordinate of the center
normal = vm.Vector3D(0,0,1) #Choose the normal
cylinder = primitives3d.Cylinder(center, normal, radius, length=0.1, name='Cylinder')

h = 10e-3 #Height of cylinder
angle = 3*math.pi/2 #Arc's angle 

#You have to create a cutting pattern in 2D

# center2d = center.to_2d(center, plane.vectors[0], plane.vectors[1])
# segbh = vm.LineSegment2D(center2d, center2d + vm.Point2D((0,h)))
# circlestart = vm.LineSegment2D(segbh.points[1], segbh.points[1]+vm.Point2D((angle,0)))
# seghb = vm.LineSegment2D(circlestart.points[1],circlestart.points[1]-segbh.points[1])
# circlend = vm.LineSegment2D(seghb.points[1],segbh.points[0])
# edges = [segbh, circlestart, seghb, circlend]
# points = edges[0].points
# contours =  [vm.Contour2D(edges)]
#
# cylinder = vm.CylindricalFace3D(contours, cylindersurface3d, points)
#
# pts1, tangle1 = cylinder.triangulation(resolution=12)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# [pt.MPLPlot(ax=ax) for pt in pts1]
# pt1 = vm.Point3D((radius*math.cos(2*math.pi/3),
#                radius*math.sin(2*math.pi/3),
#                h/4))
# p1 = frame.OldCoordinates(pt1)
# p1.MPLPlot(ax=ax, color='r')

# shell = vm.Shell3D([cylinder])
model = vm.core.VolumeModel([cylinder], name='cylinder model')

print(model.to_step('cylinder.step'))

# model.babylonjs()

# Reading own step
step = volmdlr.step.Step('cylinder.step')
model2 = step.to_volume_model()
model2.babylonjs()