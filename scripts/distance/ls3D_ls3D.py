# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:14:46 2020

@author: Mack Pro
"""

import numpy as npy
import volmdlr as volmdlr
import volmdlr.primitives3D as primitives3D
import volmdlr.primitives2D as primitives2D
import matplotlib.pyplot as plt
import random

### Cas random
mini, maxi = -5, 5

pt1 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
pt2 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
ptmid = ( pt1 + pt2 )/2
pt3 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
pt4 = volmdlr.Point3D((random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi)))
ptmid2 = (pt3 + pt4)/2

LS1 = volmdlr.LineSegment3D(pt1, pt2)
LS2 = volmdlr.LineSegment3D(pt3, pt4)


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
pt1.MPLPlot(ax=ax)
pt2.MPLPlot(ax=ax, color='r')
LS1.MPLPlot(ax=ax)

pt3.MPLPlot(ax=ax, color='g')
pt4.MPLPlot(ax=ax, color='b')
LS2.MPLPlot(ax=ax)
ptmid.MPLPlot(ax=ax)
ptmid2.MPLPlot(ax=ax)

d_min = LS1.minimum_distance(LS2)
print(d_min)