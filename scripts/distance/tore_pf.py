# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 17:08:56 2020

@author: Mack Pro
"""


import numpy as npy
import volmdlr as volmdlr
import volmdlr.primitives3d as p3d
import volmdlr.primitives2d as p2d
import volmdlr.faces as vmf
import volmdlr.edges as vme
import volmdlr.wires as vmw
import matplotlib.pyplot as plt
import random
import math

rmin, rmax = 100, 1000
mini, maxi = -1, 1

R1 = random.randrange(rmin, rmax, 1)/1000 #Radius of the generative arc3D
r1 = random.randrange(rmin/10, rmax/10, 1)/1000 #Radius of the arc3d generated

c1 = volmdlr.Point3D.random(mini, maxi, mini, maxi, mini, maxi)#Choose the coordinate of the center


n1 = volmdlr.Vector3D.random(mini, maxi, mini, maxi, mini, maxi) #Choose the normal
n1.normalize() #Normalize the normal if it is not the case
plane1 = vmf.Plane3D.from_normal(c1, n1) #Create a plane to give us two others vector

frame1 = volmdlr.Frame3D(c1, plane1.frame.u, plane1.frame.v, n1) #Frame in the center of the Tore
toresurface1 = vmf.ToroidalSurface3D(frame1, R1, r1)

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

toroidalface1 = vmf.ToroidalFace3D(toresurface1, surface2d_1)

##### extrusion 1 

p1=volmdlr.Point2D(0, 0)
p2=volmdlr.Point2D(0.1, 0.)
p3=volmdlr.Point2D(0.1, 0.2)
p4=volmdlr.Point2D(0.05, 0.1)
p5=volmdlr.Point2D(0.,0.21)
p6=volmdlr.Point2D(0.05, 0.05)

p7 = volmdlr.Point2D(0.06, 0.05)
p8 = volmdlr.Point2D(0.04, 0.07)

radius = {0: 0.01, 2: 0.01, 3: 0.015}

outer_profile = p2d.ClosedRoundedLineSegments2D([p1, p2, p3, p4, p5], radius)
l1 = vme.LineSegment2D(p6, p7)
l2 = vme.LineSegment2D(p7, p8)
l3 = vme.LineSegment2D(p8, p6)
c2 = vmw.Contour2D([l1,l2,l3])


profile=p3d.ExtrudedProfile(volmdlr.O3D, volmdlr.Y3D, volmdlr.Z3D, outer_profile, [c2], volmdlr.X3D*0.1, name = 'extrusion')
dmin, p1 ,p2 = toroidalface1.minimum_distance(profile.faces[0], return_points=True)
face_min = profile.faces[0]
for face in profile.faces[1:] :
    dtest, ptest1, ptest2 = toroidalface1.minimum_distance(face, return_points=True)
    if dtest < dmin :
        p1, p2 = ptest1, ptest2
        dmin = dtest
        face_min = face
 
print('distance ', dmin)    
 
spheres = [p3d.Sphere(p1, 5e-3, color = (250,0,0)),
           p3d.Sphere(p2, 5e-3, color = (0,0,250))]

vol = volmdlr.core.VolumeModel([toroidalface1, profile, face_min] + spheres)
vol.babylonjs()
