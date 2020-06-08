# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 10:51:54 2020

@author: Mack Pro
"""


import volmdlr as vm
import math

radius = 5e-3 #Choose the radius
center = vm.Point3D([0,0,0]) #Choose the coordinate of the center
normal = vm.Vector3D([1,1,1]) #Choose the normal
normal.Normalize() #Normalize the normal if it is not the case
plane = vm.Plane3D.from_normal(center, normal) #Create a plane to give us two others vector

frame = vm.Frame3D(center, plane.vectors[0], plane.vectors[1], normal) #Frame in the center of the cylinder
cylindersurface3d = vm.CylindricalSurface3D(frame, radius*1000) #*1000 because cylsurf3d /1000

h = 10e-3 #Height of cylinder
angle = 3*math.pi/2 #Arc's angle 

#You have to create a cutting pattern in 2D

center2d = center.To2D(center, plane.vectors[0], plane.vectors[1])
segbh = vm.LineSegment2D(center2d, center2d + vm.Point2D((0,h))) 
circlestart = vm.LineSegment2D(segbh.points[1], segbh.points[1]+vm.Point2D((angle*radius,0)))
seghb = vm.LineSegment2D(circlestart.points[1],circlestart.points[1]-segbh.points[1])
circlend = vm.LineSegment2D(seghb.points[1],segbh.points[0])
edges = [segbh, circlestart, seghb, circlend]
points = edges[0].points 
contours =  [vm.Contour2D(edges)]

cylinder = vm.CylindricalFace3D(contours, cylindersurface3d, points)

shell = vm.Shell3D([cylinder])
m = vm.VolumeModel([shell])
m.babylonjs()