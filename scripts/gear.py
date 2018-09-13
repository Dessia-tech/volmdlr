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
r4=0.015

pc=vm.Point2D((0,0))

p1=vm.Point2D((0,r2))
p2=p1.Rotation(pc,0.5*theta1)
p3=p1.Rotation(pc,theta1)
a1=vm.Arc2D(p1,p2,p3)

p4=vm.Point2D((0,r1))
p4.Rotation(pc,theta1+theta3,False)
p5=p4.Rotation(pc,theta2)
p6=p1.Rotation(pc,theta)

l1=primitives2D.RoundedLineSegments2D([p3,p4,p5,p6],{1:r3,2:r3})

#l1=primitives2D.RoundedLines2D([p1,p2,p3,p4],{0:0.01,2:0.01})
#l2=vm.Circle2D(p5,0.01)
L=[a1,l1]
for i in range(Z-1):
    thetar=(i+1)*theta
    L.append(a1.Rotation(pc,thetar,True))
    L.append(l1.Rotation(pc,thetar,True))
#p7=vm.Point2D((0,r4))
l2=vm.Circle2D(pc,r4)

c1=vm.Contour2D(L)
c2=vm.Contour2D([l2])

po=vm.Point3D((0,0,0))
xp=vm.Vector3D((1,0,0))
yp=vm.Vector3D((0,1,0))



#c1.MPLPlot()
#extr_vect=vm.Vector3D((0,0,e))

profile_straight = primitives3D.ExtrudedProfile(po,xp,yp, c1, [c2], (0,0,e),
                                                name='straight')
#
#model_straight=vm.VolumeModel([profile_straight])

profile_helical = primitives3D.HelicalExtrudedProfile(po, xp, yp, 
                                                    vm.Vector3D((0,0,0)),
                                                    vm.Vector3D((0,0,e)),
                                                    28*3.14/180, c1, [c2],
                                                    name='helical')

model=vm.VolumeModel([('helical', [profile_helical]), ('straight', [profile_straight])])


#resp=model_straight.FreeCADExport('python','gear-straight','/usr/lib/freecad/lib/',['stl','fcstd'])
#print(resp)

resp=model.FreeCADExport('gear')
print(resp)