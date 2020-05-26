# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:17:51 2020

@author: Mack Pro
"""

import numpy as npy
import volmdlr as volmdlr
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
import random

##### cas 9 arc/arc

mini, maxi = -5, 5
rad_min, rad_max = -2, 2

pt1 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
rad1 = random.randint(rad_min, rad_max)
start1, interior1, end1 = pt1, pt1 + volmdlr.Point3D((0,-rad1,rad1)), pt1 + volmdlr.Point3D((0,-2*rad1,0))
arc1 = volmdlr.Arc3D(start1, interior1, end1)

pt2 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
rad2 = random.randint(rad_min, rad_max)
start2, interior2, end2 = pt2, pt2 + volmdlr.Point3D((-2*rad2,0,0)), pt2 + volmdlr.Point3D((-rad2,rad2,0))
arc2 = volmdlr.Arc3D(start2, interior2, end2)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
pt1.MPLPlot(ax=ax)
arc1.MPLPlot(ax=ax)
end1.MPLPlot(ax=ax)
interior1.MPLPlot(ax=ax)

pt2.MPLPlot(ax=ax, color='r')
interior2.MPLPlot(ax=ax, color='r')
end2.MPLPlot(ax=ax,color='r')
arc2.MPLPlot(ax=ax)

pta1, pta2 = arc1.minimum_distance_points_arc(arc2)
pta1.MPLPlot(ax=ax, color='m')
pta2.MPLPlot(ax=ax, color='m')

d_min = arc1.minimum_distance(arc2)
print(d_min)