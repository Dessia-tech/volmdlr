# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:17:00 2020

@author: Mack Pro
"""


import numpy as npy
import volmdlr as volmdlr
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
import random

### Cas arc/LS

mini, maxi = -5, 5

pt1 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
pt2 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
ptmid = ( pt1 + pt2 )/2
pt_midmid = pt1 + (pt2-pt1)/4
pt_midmid2 = pt2 + (pt1-pt2)/4
LS1 = volmdlr.LineSegment3D(pt1, pt2)

pt = volmdlr.Point3D((random.randint(2*mini, 2*maxi),random.randint(2*mini, 2*maxi),random.randint(2*mini, 2*maxi)))
radius = 2
start, interior, end = pt, pt + volmdlr.Point3D((0,-radius,radius)),pt + volmdlr.Point3D((0,-radius,-radius))
arc = volmdlr.Arc3D(start, interior, end)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
pt1.MPLPlot(ax=ax)
pt2.MPLPlot(ax=ax, color='r')
LS1.MPLPlot(ax=ax)
start.MPLPlot(ax=ax,color='g')
interior.MPLPlot(ax=ax,color='b')
end.MPLPlot(ax=ax,color='y')
arc.MPLPlot(ax=ax)
ptmid.MPLPlot(ax=ax)

pta1, pta2 = arc.minimum_distance_points_line(LS1)
pta1.MPLPlot(ax=ax, color='m')
pta2.MPLPlot(ax=ax, color='m')

print('int',(interior-pt2).Norm(), (interior-pt1).Norm())
print('start',(start-pt2).Norm(), (start-pt1).Norm())
print('end',(end-pt2).Norm(), (end-pt1).Norm())

d_min = LS1.minimum_distance(arc)
print('d_min',d_min)