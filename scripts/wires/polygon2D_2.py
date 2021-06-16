#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 21 16:38:31 2017

@author: steven
"""

import volmdlr as vm
import volmdlr.wires as vmw
import math

r1 = 1.78*0.5
r2 = r1+0.3

theta1 = 12*2*math.pi/360
theta2 = 33*2*math.pi/360


pm1=vm.Point2D(0, -r1)
pm2=vm.Point2D(0, -r2)
pc=vm.Point2D(0, 0)
p1=pm1.rotation(pc,-theta2)
p2=pm1.rotation(pc,-theta1)
p3=pm1.rotation(pc,theta1)
p4=pm1.rotation(pc,theta2)

p8=pm2.rotation(pc,-theta2)
p7=pm2.rotation(pc,-theta1)
p6=pm2.rotation(pc,theta1)
p5=pm2.rotation(pc,theta2)

border = vmw.ClosedPolygon2D([p1,p2,p3,p4,p5,p6,p7,p8], name='border')

ptest=vm.Point2D(-0.08366,-0.91306)

projections = []
for line in border.primitives:
    projections.append(line.point_projection(ptest))



#print('test: ',border.PointDistance(ptest)+(ptest.vector[1]-p2.vector[1]))
cog_p = border.center_of_mass()

ax = border.plot()
ptest.plot(ax=ax)
cog_p.plot(ax=ax, color='r')


mesh = border.triangulation()
mesh.plot()
