# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:19:52 2020

@author: Mack Pro
"""


import numpy as npy
import volmdlr as volmdlr
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
import random

#### Cyl Cyl
rmin, rmax = 10, 100
posmin, posmax = -100, 100
x1, y1, z1 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100
x2, y2, z2 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100

r1, r2 = random.randrange(rmin, rmax, 1)/1000, random.randrange(rmin, rmax, 1)/1000 #Choose the radius
c1, c2 = volmdlr.Point3D([x1,y1,z1]), volmdlr.Point3D([x2,y2,z2]) #Choose the coordinate of the center

x3, y3, z3 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100
x4, y4, z4 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100

n1, n2 = volmdlr.Vector3D([x3,y3,z3]), volmdlr.Vector3D([x4,y4,z4]) #Choose the normal
n1.Normalize() #Normalize the normal if it is not the case
n2.Normalize()
plane1, plane2 = volmdlr.Plane3D.from_normal(c1, n1), volmdlr.Plane3D.from_normal(c2, n2) #Create a plane to give us two others vector

frame1 = volmdlr.Frame3D(c1, plane1.vectors[0], plane1.vectors[1], n1) #Frame in the center of the cylinder
frame2 = volmdlr.Frame3D(c2, plane2.vectors[0], plane2.vectors[1], n2)
cylsurface1 = volmdlr.CylindricalSurface3D(frame1, r1*1000) #*1000 because cylsurf3d /1000
cylsurface2 = volmdlr.CylindricalSurface3D(frame2, r2*1000)

hmin, hmax = -50, 50

h1, h2 = random.randrange(hmin, hmax, 5)/100, random.randrange(hmin, hmax, 5)/100 #Height of cylinder

center2d = c1.To2D(c1, plane1.vectors[0], plane1.vectors[1])
segbh = volmdlr.LineSegment2D(center2d, center2d + volmdlr.Point2D((0,h1))) #### Minus Pt2D because of Step adaptation
circlestart = volmdlr.LineSegment2D(segbh.points[1], segbh.points[1]+volmdlr.Point2D((2*math.pi*r1*3/4,0))) #You can change 2*pi by an other angle
seghb = volmdlr.LineSegment2D(circlestart.points[1],circlestart.points[1]-segbh.points[1])
circlend = volmdlr.LineSegment2D(seghb.points[1],segbh.points[0])
edges = [segbh, circlestart, seghb, circlend]
points = edges[0].points 
contours =  [volmdlr.Contour2D(edges)]

center2d2 = c2.To2D(c2, plane2.vectors[0], plane2.vectors[1])
segbh2 = volmdlr.LineSegment2D(center2d2, center2d2 + volmdlr.Point2D((0,h2))) #### Minus Pt2D because of Step adaptation
circlestart2 = volmdlr.LineSegment2D(segbh2.points[1], segbh2.points[1]+volmdlr.Point2D((2*math.pi*r2,0))) #You can change 2*pi by an other angle
seghb2 = volmdlr.LineSegment2D(circlestart2.points[1],circlestart2.points[1]-segbh2.points[1])
circlend2 = volmdlr.LineSegment2D(seghb2.points[1],segbh2.points[0])
edges2 = [segbh2, circlestart2, seghb2, circlend2]
points2 = edges2[0].points 
contours2 =  [volmdlr.Contour2D(edges2)]


cyl1 = volmdlr.CylindricalFace3D(contours, cylsurface1, points)
cyl2 = volmdlr.CylindricalFace3D(contours2, cylsurface2, points2)

pts1, tangle1 = cyl1.triangulation(resolution=12)
pts2, tangle2 = cyl2.triangulation(resolution=12)

p1, p2 = cyl1.minimum_distance_points_cyl(cyl2)
print(p1.point_distance(p2))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
[pt.MPLPlot(ax=ax) for pt in pts1]
[pt.MPLPlot(ax=ax) for pt in pts2]
p1.MPLPlot(ax=ax, color='r')
p2.MPLPlot(ax=ax, color='b')


shell = volmdlr.Shell3D([cyl1,cyl2])
m = volmdlr.VolumeModel([shell])
m.babylonjs()