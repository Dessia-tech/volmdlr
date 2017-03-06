#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  6 14:51:21 2017

@author: steven
"""

import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D

import math

r1=0.032
r2=0.04
Z=12
theta=2*math.pi/Z
theta1=0.4*theta
theta2=0.4*theta
theta3=theta-theta1-theta2
r3=0.001
e=0.030

pc=vm.Point2D((0,0))

p1=vm.Point2D((0,r2))
p2=p1.Rotation(pc,0.5*theta1)
p3=p1.Rotation(pc,theta1)
a1=vm.Arc2D(p1,p2,p3)

p4=vm.Point2D((0,r1))
p4.Rotation(pc,theta1+theta3,False)
p5=p4.Rotation(pc,theta2)
p6=p1.Rotation(pc,theta)
l1=primitives2D.RoundedLines2D([p3,p4,p5,p6],{1:r3,2:r3})

#l1=primitives2D.RoundedLines2D([p1,p2,p3,p4],{0:0.01,2:0.01})
#l2=vm.Circle2D(p5,0.01)
L=[a1,l1]
for i in range(Z-1):
    thetar=(i+1)*theta
    L.append(a1.Rotation(pc,thetar))
    L.append(l1.Rotation(pc,thetar))
c1=vm.Contour2D(L)

po=vm.Point3D((0,0,0))
xp=vm.Vector3D((1,0,0))
yp=vm.Vector3D((0,1,0))
c1.MPLPlot()

profile=primitives3D.ExtrudedProfile(po,xp,yp,[c1],(0,0,e))

model=vm.VolumeModel([profile])


model.FreeCADExport('python','gear','/usr/lib/freecad/lib/')