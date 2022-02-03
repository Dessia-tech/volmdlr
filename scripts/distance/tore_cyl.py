# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 14:56:13 2020

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
mini, maxi = -0.5, 0.5

R1 = random.randrange(rmin, rmax, 1)/1000 #Radius of the generative arc3D
r1, r2 = random.randrange(rmin/10, rmax/10, 1)/1000, random.randrange(rmin/2, rmax/2, 1)/1000 #Radius of the arc3d generated

c1 = volmdlr.Point3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center
c2 = volmdlr.Point3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center

n1 = volmdlr.Vector3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center
n2 = volmdlr.Vector3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the coordinate of the center

n1.normalize() #Normalize the normal if it is not the case
n2.normalize()
plane1, plane2 = vmf.Plane3D.from_normal(c1, n1), vmf.Plane3D.from_normal(c2, n2) #Create a plane to give us two others vector

toresurface1 = vmf.ToroidalSurface3D(plane1.frame, R1, r1) 
cylsurface2 = vmf.CylindricalSurface3D(plane2.frame, r2)


angle_min, angle_max = 0, 2*3.14*100

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

#Cylinder
hmin, hmax = -50, 50

h2 = random.randrange(hmin, hmax, 5)/100 #Height of cylinder
if h2 == 0 :
    h2 = 1/100
angle_cyl = random.randrange(angle_min, angle_max, 20)/100

center2d2 = c2.to_2d(c2, plane2.frame.u, plane2.frame.v)
segbh2 = vme.LineSegment2D(center2d2, center2d2 + volmdlr.Point2D(0,h2) + volmdlr.Point2D(angle_cyl/3,0)) #### Minus Pt2D because of Step adaptation
circlestart2 = vme.LineSegment2D(segbh2.end, segbh2.end+volmdlr.Point2D(angle_cyl,0) - volmdlr.Point2D(0,h2/10)) #You can change 2*pi by an other angle
seghb2 = vme.LineSegment2D(circlestart2.end,circlestart2.end-segbh2.end + volmdlr.Point2D(angle_cyl/3,0))
circlend2 = vme.LineSegment2D(seghb2.end,segbh2.start)
edges2 = [segbh2, circlestart2, seghb2, circlend2]
surface2d_2 = vmf.Surface2D(outer_contour = vmw.Contour2D(edges2),
                            inner_contours = [])


toroidalface1 = vmf.ToroidalFace3D(toresurface1, surface2d_1)
cyl2 = vmf.CylindricalFace3D(cylsurface2, surface2d_2)

distance, p1, p2 = toroidalface1.minimum_distance(cyl2, return_points=True)
print('distance ', distance)

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

vol = volmdlr.core.VolumeModel([toroidalface1,cyl2] + spheres)
vol.babylonjs()
# m = volmdlr.VolumeModel([shell])
# m.babylonjs()