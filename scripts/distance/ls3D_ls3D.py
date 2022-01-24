# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:14:46 2020

@author: Mack Pro
"""

import numpy as npy
import volmdlr as volmdlr
import volmdlr.primitives3d as primitives3D
import volmdlr.primitives2d as primitives2D
import volmdlr.edges as vme
import matplotlib.pyplot as plt
import random

### Cas random
mini, maxi = -5, 5

pt1 = volmdlr.Point3D(random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi))
pt2 = volmdlr.Point3D(random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi))
ptmid = ( pt1 + pt2 )/2
pt3 = volmdlr.Point3D(random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi))
pt4 = volmdlr.Point3D(random.randint(mini, maxi),random.randint(mini, maxi),random.randint(mini, maxi))
ptmid2 = (pt3 + pt4)/2

LS1 = vme.LineSegment3D(pt1, pt2)
LS2 = vme.LineSegment3D(pt3, pt4)


ax = pt1.plot()
pt2.plot(ax=ax, color='r')
LS1.plot(ax=ax)

pt3.plot(ax=ax, color='g')
pt4.plot(ax=ax, color='b')
LS2.plot(ax=ax)
ptmid.plot(ax=ax)
ptmid2.plot(ax=ax)

d_min = LS1.minimum_distance(LS2)
print(d_min)