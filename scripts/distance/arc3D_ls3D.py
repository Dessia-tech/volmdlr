# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:17:00 2020

@author: Mack Pro
"""


import volmdlr as vm
import volmdlr.edges as vme
import matplotlib.pyplot as plt
import random

### Cas arc/LS

mini, maxi = -5, 5

pt1 = vm.Point3D.random(mini, maxi, mini, maxi, mini, maxi)
pt2 = vm.Point3D.random(mini, maxi, mini, maxi, mini, maxi)
ptmid = ( pt1 + pt2 )/2
pt_midmid = pt1 + (pt2-pt1)/4
pt_midmid2 = pt2 + (pt1-pt2)/4
LS1 = vme.LineSegment3D(pt1, pt2)

pt = vm.Point3D.random(mini, maxi, mini, maxi, mini, maxi)
radius = 2
start, interior, end = pt, pt + vm.Point3D(0,-radius,radius),pt + vm.Point3D(0,-radius,-radius)
arc = vme.Arc3D(start, interior, end)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
pt1.plot(ax=ax)
pt2.plot(ax=ax, color='r')
LS1.plot(ax=ax)
start.plot(ax=ax,color='g')
interior.plot(ax=ax,color='b')
end.plot(ax=ax,color='y')
arc.plot(ax=ax)
ptmid.plot(ax=ax)

pta1, pta2 = arc.minimum_distance_points_line(LS1)
pta1.plot(ax=ax, color='m')
pta2.plot(ax=ax, color='m')

print('int',(interior-pt2).norm(), (interior-pt1).norm())
print('start',(start-pt2).norm(), (start-pt1).norm())
print('end',(end-pt2).norm(), (end-pt1).norm())

d_min = LS1.minimum_distance(arc)
print('d_min',d_min)