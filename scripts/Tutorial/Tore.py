# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 15:25:47 2020

@author: Mack Pro
"""


import volmdlr as vm
import math
# import matplotlib.pyplot as plt

r1 = 10e-3 #Radius of the generative arc3D
r2 = 3e-3 #Radius of the arc3d generated

center = vm.Point3D([0,0,0]) #Choose the coordinate of the center
normal1 = vm.Vector3D([0,0,1]) #Choose the normal of the generative
normal1.Normalize() #Normalize the normal if it is not the case
vec1 = normal1.deterministic_unit_normal_vector()


frame = vm.Frame3D(center, vec1, normal1.Cross(vec1), normal1) #Frame in the center of the generative arc3D
toroidalsurface3d = vm.ToroidalSurface3D(frame, r1, r2)

theta = 4*math.pi/3 #Tore's length
phi = 2*math.pi #angle of circle 
offset_theta = math.pi/4 #Theta's offset if you want to turn it with normal's reference
offset_phi = math.pi #Idem but with circle's normal

#You have to create a cutting pattern in 2D

pt1, pt2, pt3, pt4 = vm.Point2D((offset_theta, offset_phi)), vm.Point2D((offset_theta, offset_phi+phi)), vm.Point2D((offset_theta+theta, offset_phi+phi)), vm.Point2D((offset_theta+theta, offset_phi))
seg1, seg2, seg3, seg4 = vm.LineSegment2D(pt1, pt2), vm.LineSegment2D(pt2, pt3), vm.LineSegment2D(pt3, pt4), vm.LineSegment2D(pt4, pt1) 
edges = [seg1, seg2, seg3, seg4]
contours2d =  [vm.Contour2D(edges)]
points = [theta, phi] 



toroidalface = vm.ToroidalFace3D(contours2d, toroidalsurface3d, points)

shell = vm.Shell3D([toroidalface])
m = vm.VolumeModel([shell])
m.babylonjs(debug=True)