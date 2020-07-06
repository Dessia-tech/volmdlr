# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 12:11:37 2020

@author: Mack Pro
"""


import numpy as npy
import volmdlr as volmdlr
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
import random
import math

rmin, rmax = 10, 100
posmin, posmax = -100, 100
x1, y1, z1 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100

r1 = random.randrange(rmin, rmax, 1)/1000
c1 = volmdlr.Point3D([x1,y1,z1])

x3, y3, z3 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100

n1 = volmdlr.Vector3D([x3,y3,z3])
n1.Normalize() #Normalize the normal if it is not the case
plane1 = volmdlr.Plane3D.from_normal(c1, n1)

frame1 = volmdlr.Frame3D(c1, plane1.vectors[0], plane1.vectors[1], n1) #Frame in the center of the cylinder
cylsurface1 = volmdlr.CylindricalSurface3D(frame1, r1)

hmin, hmax = -50, 50

h1 = random.randrange(hmin, hmax, 5)/100

center2d = c1.To2D(c1, plane1.vectors[0], plane1.vectors[1])
#Classic Contour
segbh = volmdlr.LineSegment2D(center2d, center2d + volmdlr.Point2D((0,h1))) #### Minus Pt2D because of Step adaptation
circlestart = volmdlr.LineSegment2D(segbh.points[1], segbh.points[1]+volmdlr.Point2D((2*math.pi*r1*3/4,0))) #You can change 2*pi by an other angle
seghb = volmdlr.LineSegment2D(circlestart.points[1],circlestart.points[1]-segbh.points[1])
circlend = volmdlr.LineSegment2D(seghb.points[1],segbh.points[0])
edges = [segbh, circlestart, seghb, circlend]
points = edges[0].points 
contours =  [volmdlr.Contour2D(edges)]
# segbh = volmdlr.LineSegment2D(center2d, center2d + volmdlr.Point2D((0,h1/2))) #### Minus Pt2D because of Step adaptation
# circlestart = volmdlr.LineSegment2D(segbh.points[1], segbh.points[1]+volmdlr.Point2D((2*math.pi*r1*3/4,-h1/10))) #You can change 2*pi by an other angle
# seghb = volmdlr.LineSegment2D(circlestart.points[1],circlestart.points[1]-segbh.points[1])
# circlend = volmdlr.LineSegment2D(seghb.points[1],segbh.points[0])
# edges = [segbh, circlestart, seghb, circlend]
# points = edges[0].points 
# contours =  [volmdlr.Contour2D(edges)]

cyl1 = volmdlr.CylindricalFace3D(contours, cylsurface1, points)

pts1, tangle1 = cyl1.triangulation(resolution=12)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
[pt.MPLPlot(ax=ax) for pt in pts1]





# number_holes = 5

# outer_circle = volmdlr.Circle2D(volmdlr.O2D, 0.06)

# circles = []
# delta_angle = 2*math.pi/number_holes
# inner_circle = volmdlr.Circle2D(volmdlr.O2D, 0.04)
# first_circle = volmdlr.Circle2D(volmdlr.Point2D((0, 0.05)), 0.005)
# circles = [inner_circle, first_circle]
# for i in range(1, number_holes):
#     circles.append(first_circle.Rotation(volmdlr.O2D, i*delta_angle))
    
# xn, yn, zn = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100
# xc, yc, zc = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100

# n = volmdlr.Vector3D([xn,yn,zn])
# n.Normalize()
# c = volmdlr.Point3D((xc,yc,zc))
# plane = volmdlr.Plane3D.from_normal(c, n)
# contours = [volmdlr.Contour2D([outer_circle])]#, inner_circle])]#, volmdlr.Contour2D(circles)]
# planeface = volmdlr.PlaneFace3D(contours,  plane)


# p1, p2 = cyl1.minimum_distance_points_plane(planeface)
# # print(p1.point_distance(p2))
# # print(p1)
# # print(p2)
# pts2, t2 = planeface.triangulation()
# [pt.MPLPlot(ax=ax) for pt in pts2]
# p1.MPLPlot(ax=ax, color='r')
# p2.MPLPlot(ax=ax, color='b')


##### extrusion 1 

p1=volmdlr.Point2D((0, 0))
p2=volmdlr.Point2D((0.1, 0.))
p3=volmdlr.Point2D((0.1, 0.2))
p4=volmdlr.Point2D((0.05, 0.1))
p5=volmdlr.Point2D((0.,0.21))
p6=volmdlr.Point2D((0.05, 0.05))

p7 = volmdlr.Point2D((0.06, 0.05))
p8 = volmdlr.Point2D((0.04, 0.07))

radius = {0: 0.01, 2: 0.01, 3: 0.015}

outer_profile = primitives2D.ClosedRoundedLineSegments2D([p1, p2, p3, p4, p5], radius)
#hole = volmdlr.Circle2D(p6, 0.01)
#inner_profile = primitives2D.RoundedLineSegments2D([p6, p7, p8], {0: 0.5}, closed = True)
l1 = volmdlr.LineSegment2D(p6, p7)
l2 = volmdlr.LineSegment2D(p7, p8)
l3 = volmdlr.LineSegment2D(p8, p6)
c2 = volmdlr.Contour2D([l1,l2,l3])

#c1 = volmdlr.Contour2D([outer_profile])
#c2 = volmdlr.Contour2D([inner_profile])

# f, a = outer_profile.MPLPlot()
# c2.MPLPlot(a)

profile=primitives3D.ExtrudedProfile(volmdlr.O3D, volmdlr.Y3D, volmdlr.Z3D, outer_profile, [c2], volmdlr.X3D*0.1, name = 'extrusion')
dmin, p1 ,p2 = cyl1.minimum_distance(profile.faces[0], return_points=True)
face_min = profile.faces[0]
# print('dmin', dmin)
for face in profile.faces[1:] :
    dtest, ptest1, ptest2 = cyl1.minimum_distance(face, return_points=True)
    # print('dtest', dtest)
    if dtest < dmin :
        p1, p2 = ptest1, ptest2
        dmin = dtest
        face_min = face


# print('>>>>>>>>> distance minimale', dmin)

shell = volmdlr.Shell3D([cyl1])
vol = volmdlr.VolumeModel([shell, profile, p1, p2])
vol.babylonjs_from_script()
# m = volmdlr.VolumeModel([shell, profile])
# m.babylonjs()