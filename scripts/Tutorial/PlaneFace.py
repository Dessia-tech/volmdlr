# -*- coding: utf-8 -*-
"""
Created on Tue May 26 09:45:01 2020

@author: Mack Pro
"""

import volmdlr as volmdlr
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import random
import math

posmin, posmax = -100, 100

number_holes = 5

outer_circle = volmdlr.Circle2D(volmdlr.O2D, 0.06)

circles = []
delta_angle = 2*math.pi/number_holes
inner_circle = volmdlr.Circle2D(volmdlr.O2D, 0.04)
first_circle = volmdlr.Circle2D(volmdlr.Point2D((0, 0.05)), 0.005)
circles = [inner_circle, first_circle]
for i in range(1, number_holes):
    circles.append(first_circle.Rotation(volmdlr.O2D, i*delta_angle))
    
xn, yn, zn = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100
xc, yc, zc = random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100, random.randrange(posmin, posmax, 1)/100

n = volmdlr.Vector3D([xn,yn,zn])
n.Normalize()
c = volmdlr.Point3D((xc,yc,zc))

#Plane to place your PlaneFace3D
plane = volmdlr.Plane3D.from_normal(c, n)

#To do a PlaneFace3D, you need to use contours2D
#The first is always the basis, the others are used to cut the basis

contours = [volmdlr.Contour2D([outer_circle]), volmdlr.Contour2D([circles[0]]),
            volmdlr.Contour2D([circles[1]]), volmdlr.Contour2D([circles[2]]),
            volmdlr.Contour2D([circles[3]]),volmdlr.Contour2D([circles[4]]),
            volmdlr.Contour2D([circles[5]])]

planeface = volmdlr.PlaneFace3D(contours,  plane)

shell = volmdlr.Shell3D([planeface])
m = volmdlr.VolumeModel([shell])
m.babylonjs()