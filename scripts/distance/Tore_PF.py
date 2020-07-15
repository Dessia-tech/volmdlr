# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 17:08:56 2020

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
r1 = random.randrange(rmin/10, rmax/10, 1)/1000 #Radius of the arc3d generated

c1 = volmdlr.Point3D([x1,y1,z1]) #Choose the coordinate of the center

x3, y3, z3 = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100

n1 = volmdlr.Vector3D([x3,y3,z3]) #Choose the normal
n1.Normalize() #Normalize the normal if it is not the case
plane1 = volmdlr.Plane3D.from_normal(c1, n1) #Create a plane to give us two others vector

frame1 = volmdlr.Frame3D(c1, plane1.vectors[0], plane1.vectors[1], n1) #Frame in the center of the Tore
toresurface1 = volmdlr.ToroidalSurface3D(frame1, R1*1000, r1*1000) #*1000 because torsurf3d /1000

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

toroidalface1 = volmdlr.ToroidalFace3D(contours2d, toresurface1, points)

pts1, tangle1 = toroidalface1.triangulation(resolution=10)

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
dmin, p1 ,p2 = toroidalface1.minimum_distance(profile.faces[0], return_points=True)
face_min = profile.faces[0]
# print('dmin', dmin)
for face in profile.faces[1:] :
    dtest, ptest1, ptest2 = toroidalface1.minimum_distance(face, return_points=True)
    # print('dtest', dtest)
    if dtest < dmin :
        p1, p2 = ptest1, ptest2
        dmin = dtest
        face_min = face


# print('>>>>>>>>> distance minimale', dmin)

shell = volmdlr.Shell3D([toroidalface1])
vol = volmdlr.VolumeModel([shell, profile, p1, p2])
vol.babylonjs_from_script()
# m = volmdlr.VolumeModel([shell, profile])
# m.babylonjs()