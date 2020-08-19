# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 14:56:13 2020

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


R1 = random.randrange(rmin, rmax, 1)/1000 #Radius of the generative arc3D
r1, r2 = random.randrange(rmin/10, rmax/10, 1)/1000, random.randrange(rmin/2, rmax/2, 1)/1000 #Radius of the arc3d generated

c1, c2 = volmdlr.Point3D([x1,y1,z1]), volmdlr.Point3D([x2,y2,z2]) #Choose the coordinate of the center

x3, y3, z3 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100
x4, y4, z4 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100

n1, n2 = volmdlr.Vector3D([x3,y3,z3]), volmdlr.Vector3D([x4,y4,z4]) #Choose the normal
n1.Normalize() #Normalize the normal if it is not the case
n2.Normalize()
plane1, plane2 = volmdlr.Plane3D.from_normal(c1, n1), volmdlr.Plane3D.from_normal(c2, n2) #Create a plane to give us two others vector

frame1 = volmdlr.Frame3D(c1, plane1.vectors[0], plane1.vectors[1], n1) #Frame in the center of the Tore
frame2 = volmdlr.Frame3D(c2, plane2.vectors[0], plane2.vectors[1], n2)
toresurface1 = volmdlr.ToroidalSurface3D(frame1, R1, r1) 
cylsurface2 = volmdlr.CylindricalSurface3D(frame2, r2)


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

#Cylinder
hmin, hmax = -50, 50

h2 = random.randrange(hmin, hmax, 5)/100 #Height of cylinder
angle_cyl = random.randrange(angle_min, angle_max, 20)/100

center2d2 = c2.To2D(c2, plane2.vectors[0], plane2.vectors[1])
segbh2 = volmdlr.LineSegment2D(center2d2, center2d2 + volmdlr.Point2D((0,h2)) + volmdlr.Point2D((angle_cyl/3,0))) #### Minus Pt2D because of Step adaptation
circlestart2 = volmdlr.LineSegment2D(segbh2.points[1], segbh2.points[1]+volmdlr.Point2D((angle_cyl,0)) - volmdlr.Point2D((0,h2/10))) #You can change 2*pi by an other angle
seghb2 = volmdlr.LineSegment2D(circlestart2.points[1],circlestart2.points[1]-segbh2.points[1] + volmdlr.Point2D((angle_cyl/3,0)))
circlend2 = volmdlr.LineSegment2D(seghb2.points[1],segbh2.points[0])
edges2 = [segbh2, circlestart2, seghb2, circlend2]
points2 = edges2[0].points 
contours2 =  [volmdlr.Contour2D(edges2)]


toroidalface1 = volmdlr.ToroidalFace3D(contours2d, toresurface1, points)
cyl2 = volmdlr.CylindricalFace3D(contours2, cylsurface2, points2)

pts1, tangle1 = toroidalface1.triangulation(resolution=10)
pts2, tangle2 = cyl2.triangulation(resolution=12)

distance, p1, p2 = toroidalface1.minimum_distance(cyl2, return_points=True)
print(distance)

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# [pt.MPLPlot(ax=ax) for pt in pts1]
# [pt.MPLPlot(ax=ax) for pt in pts2]
# p1.MPLPlot(ax=ax, color='r')
# p2.MPLPlot(ax=ax, color='b')
# toroidalface1.start.MPLPlot(ax=ax, color='m')
# toroidalface2.start.MPLPlot(ax=ax, color='g')

# LS = volmdlr.LineSegment3D(p1, p2)

shell = volmdlr.Shell3D([toroidalface1,cyl2])
vol = volmdlr.VolumeModel([shell, p1, p2])
vol.babylonjs_from_script()
# m = volmdlr.VolumeModel([shell])
# m.babylonjs()