# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 14:15:50 2020

@author: Mack Pro
"""

import volmdlr as volmdlr
import volmdlr.primitives3d as p3d
import volmdlr.faces as vmf
import volmdlr.edges as vme
import volmdlr.wires as vmw
import random
import math

rmin, rmax = 100, 1000
mini, maxi = -1, 1

R1, R2 = random.randrange(rmin, rmax, 1)/1000, random.randrange(rmin, rmax, 1)/1000 #Radius of the generative arc3D
r1, r2 = random.randrange(rmin/10, rmax/10, 1)/1000, random.randrange(rmin/10, rmax/10, 1)/1000 #Radius of the arc3d generated

c1 = volmdlr.Point3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center
c2 = volmdlr.Point3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center


n1 = volmdlr.Vector3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center
n2 = volmdlr.Vector3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center

n1.normalize() #Normalize the normal if it is not the case
n2.normalize()
plane1, plane2 = vmf.Plane3D.from_normal(c1, n1), vmf.Plane3D.from_normal(c2, n2) #Create a plane to give us two others vector

toresurface1 = vmf.ToroidalSurface3D(plane1.frame, R1, r1) 
toresurface2 = vmf.ToroidalSurface3D(plane2.frame, R2, r2)


angle_min, angle_max = 100, 2*3.14*100

theta1 = random.randrange(angle_min, angle_max, 20)/100 #Tore's length
phi1 = 2*math.pi #angle of circle 
offset_theta1 = random.randrange(angle_min, angle_max, 20)/100 #Theta's offset if you want to turn it with normal's reference
offset_phi1 = random.randrange(angle_min, angle_max, 20)/100 #Idem but with circle's normal

print('param1', phi1, theta1, offset_phi1, offset_theta1)

#You have to create a cutting pattern in 2D

pt1, pt2 = volmdlr.Point2D(offset_theta1, offset_phi1), volmdlr.Point2D(offset_theta1, offset_phi1+phi1)
pt3, pt4 = volmdlr.Point2D(offset_theta1+theta1, offset_phi1+phi1), volmdlr.Point2D(offset_theta1+theta1, offset_phi1)
seg1, seg2 = vme.LineSegment2D(pt1, pt2), vme.LineSegment2D(pt2, pt3)
seg3, seg4 = vme.LineSegment2D(pt3, pt4), vme.LineSegment2D(pt4, pt1)
edges = [seg1, seg2, seg3, seg4]
surface2d_1 = vmf.Surface2D(outer_contour = vmw.Contour2D(edges),
                            inner_contours = [])



theta2 = random.randrange(angle_min, angle_max, 20)/100 #Tore's length
phi2 = random.randrange(angle_min, angle_max, 20)/100 #angle of circle 
offset_theta2 = random.randrange(angle_min, angle_max, 20)/100 #Theta's offset if you want to turn it with normal's reference
offset_phi2 = random.randrange(angle_min, angle_max, 20)/100 #Idem but with circle's normal

print('param2', phi2, theta2, offset_phi2, offset_theta2)

#You have to create a cutting pattern in 2D

pt1_2, pt2_2 = volmdlr.Point2D(offset_theta2, offset_phi2), volmdlr.Point2D(offset_theta2, offset_phi2+phi2)
pt3_2, pt4_2 = volmdlr.Point2D(offset_theta2+theta2, offset_phi2+phi2), volmdlr.Point2D(offset_theta2+theta2, offset_phi2)
seg1_2, seg2_2 = vme.LineSegment2D(pt1_2, pt2_2), vme.LineSegment2D(pt2_2, pt3_2)
seg3_2, seg4_2 = vme.LineSegment2D(pt3_2, pt4_2), vme.LineSegment2D(pt4_2, pt1_2)
edges_2 = [seg1_2, seg2_2, seg3_2, seg4_2]
surface2d_2 = vmf.Surface2D(outer_contour = vmw.Contour2D(edges_2),
                            inner_contours = [])

toroidalface1 = vmf.ToroidalFace3D(toresurface1, surface2d_1)
toroidalface2 = vmf.ToroidalFace3D(toresurface2, surface2d_2)

p1, p2 = toroidalface1.minimum_distance_points_tore(toroidalface2)
print('distance', p1.point_distance(p2))

# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# [pt.MPLPlot(ax=ax) for pt in pts1]
# [pt.MPLPlot(ax=ax) for pt in pts2]
# p1.MPLPlot(ax=ax, color='r')
# p2.MPLPlot(ax=ax, color='b')
# toroidalface1.start.MPLPlot(ax=ax, color='m')
# toroidalface2.start.MPLPlot(ax=ax, color='g')

# LS = volmdlr.LineSegment3D(p1, p2)

spheres = [p3d.Sphere(p1, 5e-3, color = (250,0,0)),
           p3d.Sphere(p2, 5e-3, color = (0,0,250))]

vol = volmdlr.core.VolumeModel([toroidalface1,toroidalface2] + spheres)
vol.babylonjs()
# m = volmdlr.VolumeModel([shell])
# m.babylonjs()
