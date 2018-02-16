#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:14:49 2017

@author: steven
"""

import volmdlr as vm
import volmdlr.primitives3D as primitives3D
import math

D=55*1e-3
D1=47.6*1e-3
d=30*1e-3
d1=38.3*1e-3
B=13*1e-3
n_balls=13
D_balls=8*1e-3


y_ball0=(D+d)/4
theta_b=math.acos((y_ball0-d1/2)/D_balls/2)

x=vm.Vector3D((1,0,0))
y=vm.Vector3D((0,1,0))

primitives=[]
center=vm.Point3D((0,0,0))
center_ball1=vm.Point3D((0,y_ball0,0))

for i in range(n_balls):    
    angle=i*2*math.pi/n_balls
    center_ball=center_ball1.Rotation(center,x,angle,True)
#    print(center_ball.vector)
    primitives.append(primitives3D.Sphere(center_ball,D_balls/2,'Ball{}'.format(i+1)))

# inner
pi1=vm.Point2D((-B/2,d/2))
pi2=vm.Point2D((-B/2,d1/2))
pi3=vm.Point2D((-0.5*D_balls*math.sin(theta_b),d1/2))
pi4=vm.Point2D((0,y_ball0-D_balls/2))
pi5=vm.Point2D((0.5*D_balls*math.sin(theta_b),d1/2))
pi6=vm.Point2D((B/2,d1/2))
pi7=vm.Point2D((B/2,d/2))

li1=vm.Line2D(pi1,pi2)
li2=vm.Line2D(pi2,pi3)
ci3=vm.Arc2D(pi3,pi4,pi5)
li4=vm.Line2D(pi5,pi6)
li5=vm.Line2D(pi6,pi7)
li6=vm.Line2D(pi7,pi1)

ci=vm.Contour2D([li1,li2,ci3,li4,li5,li6])
inner=primitives3D.RevolvedProfile(center,x,y,[ci],center,x,math.pi*2,'innerring')
primitives.append(inner)

# outter
po1=vm.Point2D((-B/2,D/2))
po2=vm.Point2D((-B/2,D1/2))
po3=vm.Point2D((-0.5*D_balls*math.sin(theta_b),D1/2))
po4=vm.Point2D((0,y_ball0+D_balls/2))
po5=vm.Point2D((0.5*D_balls*math.sin(theta_b),D1/2))
po6=vm.Point2D((B/2,D1/2))
po7=vm.Point2D((B/2,D/2))

lo1=vm.Line2D(po1,po2)
lo2=vm.Line2D(po2,po3)
co3=vm.Arc2D(po3,po4,po5)
lo4=vm.Line2D(po5,po6)
lo5=vm.Line2D(po6,po7)
lo6=vm.Line2D(po7,po1)

co=vm.Contour2D([lo1,lo2,co3,lo4,lo5,lo6])
outter=primitives3D.RevolvedProfile(center,x,y,[co],center,x,math.pi*2)

primitives.append(outter)


model=vm.VolumeModel(primitives)
resp=model.FreeCADExport('/usr/bin/python','bearing','/usr/lib/freecad/lib/',['stl','fcstd'])
print(resp)