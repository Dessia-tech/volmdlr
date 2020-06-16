# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 14:15:50 2020

@author: Mack Pro
"""

import numpy as npy
import volmdlr as volmdlr
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
import random
import math

rmin, rmax = 100, 1000
posmin, posmax = -100, 100
x1, y1, z1 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100
x2, y2, z2 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100


R1, R2 = random.randrange(rmin, rmax, 1)/1000, random.randrange(rmin, rmax, 1)/1000 #Radius of the generative arc3D
r1, r2 = random.randrange(rmin/10, rmax/10, 1)/1000, random.randrange(rmin/10, rmax/10, 1)/1000 #Radius of the arc3d generated

c1, c2 = volmdlr.Point3D([x1,y1,z1]), volmdlr.Point3D([x2,y2,z2]) #Choose the coordinate of the center

x3, y3, z3 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100
x4, y4, z4 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100

n1, n2 = volmdlr.Vector3D([x3,y3,z3]), volmdlr.Vector3D([x4,y4,z4]) #Choose the normal
n1.Normalize() #Normalize the normal if it is not the case
n2.Normalize()
plane1, plane2 = volmdlr.Plane3D.from_normal(c1, n1), volmdlr.Plane3D.from_normal(c2, n2) #Create a plane to give us two others vector

frame1 = volmdlr.Frame3D(c1, plane1.vectors[0], plane1.vectors[1], n1) #Frame in the center of the Tore
frame2 = volmdlr.Frame3D(c2, plane2.vectors[0], plane2.vectors[1], n2)
toresurface1 = volmdlr.ToroidalSurface3D(frame1, R1*1000, r1*1000) #*1000 because torsurf3d /1000
toresurface2 = volmdlr.ToroidalSurface3D(frame2, R2*1000, r2*1000)


angle_min, angle_max = 0, 2*3.14*100

theta1 = random.randrange(angle_min, angle_max, 20)/100 #Tore's length
phi1 = 2*math.pi #angle of circle 
offset_theta1 = random.randrange(angle_min, angle_max, 20)/100 #Theta's offset if you want to turn it with normal's reference
offset_phi1 = random.randrange(angle_min, angle_max, 20)/100 #Idem but with circle's normal

print('param1', phi1, theta1, offset_phi1, offset_theta1)

#You have to create a cutting pattern in 2D

pt1, pt2, pt3, pt4 = volmdlr.Point2D((offset_theta1, offset_phi1)), volmdlr.Point2D((offset_theta1, offset_phi1+phi1)), volmdlr.Point2D((offset_theta1+theta1, offset_phi1+phi1)), volmdlr.Point2D((offset_theta1+theta1, offset_phi1))
seg1, seg2, seg3, seg4 = volmdlr.LineSegment2D(pt1, pt2), volmdlr.LineSegment2D(pt2, pt3), volmdlr.LineSegment2D(pt3, pt4), volmdlr.LineSegment2D(pt4, pt1) 
edges = [seg1, seg2, seg3, seg4]
contours2d =  [volmdlr.Contour2D(edges)]
points = [theta1, phi1] 



theta2 = random.randrange(angle_min, angle_max, 20)/100 #Tore's length
phi2 = random.randrange(angle_min, angle_max, 20)/100 #angle of circle 
offset_theta2 = random.randrange(angle_min, angle_max, 20)/100 #Theta's offset if you want to turn it with normal's reference
offset_phi2 = random.randrange(angle_min, angle_max, 20)/100 #Idem but with circle's normal

print('param2', phi2, theta2, offset_phi2, offset_theta2)

#You have to create a cutting pattern in 2D

pt1_2, pt2_2, pt3_2, pt4_2 = volmdlr.Point2D((offset_theta2, offset_phi2)), volmdlr.Point2D((offset_theta2, offset_phi2+phi2)), volmdlr.Point2D((offset_theta2+theta2, offset_phi2+phi2)), volmdlr.Point2D((offset_theta2+theta2, offset_phi2))
seg1_2, seg2_2, seg3_2, seg4_2 = volmdlr.LineSegment2D(pt1_2, pt2_2), volmdlr.LineSegment2D(pt2_2, pt3_2), volmdlr.LineSegment2D(pt3_2, pt4_2), volmdlr.LineSegment2D(pt4_2, pt1_2) 
edges_2 = [seg1_2, seg2_2, seg3_2, seg4_2]
contours2d_2 =  [volmdlr.Contour2D(edges_2)]
points_2 = [theta2, phi2]

toroidalface1 = volmdlr.ToroidalFace3D(contours2d, toresurface1, points)
toroidalface2 = volmdlr.ToroidalFace3D(contours2d_2, toresurface2, points_2)

pts1, tangle1 = toroidalface1.triangulation(resolution=10)
pts2, tangle2 = toroidalface2.triangulation(resolution=10)

p1, p2 = toroidalface1.minimum_distance_points_tore(toroidalface2)
print('p1, p2', p1,p2)
print(p1.point_distance(p2))

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# [pt.MPLPlot(ax=ax) for pt in pts1]
# [pt.MPLPlot(ax=ax) for pt in pts2]
# p1.MPLPlot(ax=ax, color='r')
# p2.MPLPlot(ax=ax, color='b')
# toroidalface1.start.MPLPlot(ax=ax, color='m')
# toroidalface2.start.MPLPlot(ax=ax, color='g')

# LS = volmdlr.LineSegment3D(p1, p2)

shell = volmdlr.Shell3D([toroidalface1,toroidalface2])
vol = volmdlr.VolumeModel([shell, p1, p2])
vol.babylonjs_from_script()
# m = volmdlr.VolumeModel([shell])
# m.babylonjs()
