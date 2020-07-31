# -*- coding: utf-8 -*-
"""
Created on Thu Apr  2 14:09:52 2020

@author: Mack Pro
"""


import volmdlr as vm
import math

Gradius = 8e-3 #Grand radius
Sradius = 3e-3 #Small radius
h = 20e-3 #Height of the Ellipse

center = vm.Point3D([0,0,0]) #Choose the coordinate of the center
normal = vm.Vector3D([1,1,1]) #Choose the normal
normal.Normalize() #Normalize the normal if it is not the case
plane = vm.Plane3D.from_normal(center, normal) #Create a plane to give us two others vector
Gdir = plane.vectors[0] #Gradius direction
Sdir = plane.vectors[1] #Sradius direction 

ellipse = vm.Ellipse3D(Gradius, Sradius, center, normal, Gdir)
bsplinextru = vm.BSplineExtrusion(ellipse, -normal) #Perhaps the normal needs to be the opposite

angle = 5*math.pi/4 #Angle of pint if start=end

# position of point in ellipse with an angle
# Gradius*vm.Point3D(Gdir.vector)*math.cos(angle)+Sradius*vm.Point3D(Sdir.vector)*math.sin(angle)+center


pointellipse = center + Gradius*Gdir #Point on Ellipse
pint = Gradius*vm.Point3D(Gdir.vector)*math.cos(angle)+Sradius*vm.Point3D(Sdir.vector)*math.sin(angle)+center
extra1 = Gradius*vm.Point3D(Gdir.vector)*math.cos(math.pi/2)+Sradius*vm.Point3D(Sdir.vector)*math.sin(math.pi/2)+center

segbh = vm.LineSegment3D(pointellipse, pointellipse + vm.Point3D([i*h for i in normal.vector])) #point on the ellipse not on the center

#IF you want to do a complete ellipse, you need to add an 'extra' point ---> see angle

ellipse1 = vm.ArcEllipse3D(pointellipse, pint, pointellipse, center,Gdir, normal, 'ellipse1', extra = extra1)
seghb = vm.LineSegment3D(segbh.points[1], segbh.points[0])

center2 = center + vm.Point3D([i*h for i in normal.vector])
pointellipse2 = center2 + Gradius*Gdir
pint2 = pint + vm.Point3D([i*h for i in normal.vector])
extra2 = extra1 + vm.Point3D([i*h for i in normal.vector])
ellipse2 = vm.ArcEllipse3D(pointellipse2, pint2, pointellipse2, center2, Gdir, normal, 'ellipse2',extra = extra2)

edges = [segbh, ellipse1, seghb, ellipse2]
points = segbh.points+ellipse1.points+seghb.points+ellipse2.points
contours = [vm.Contour3D(edges)]

EllipseFace = vm.BSplineFace3D(contours, bsplinextru, points)

shell = vm.Shell3D([EllipseFace])
m = vm.VolumeModel([shell])
m.babylonjs()