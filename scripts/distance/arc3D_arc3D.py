# -*- coding: utf-8 -*-
"""
Created on Fri May 22 11:17:51 2020

@author: Mack Pro
"""

import volmdlr as vm
import volmdlr.edges as vme
import matplotlib.pyplot as plt
import random

##### cas 9 arc/arc

mini, maxi = -5, 5
rad_min, rad_max = 1, 3

pt1 = vm.Point3D.random(mini, maxi, mini, maxi, mini, maxi)
rad1 = random.randint(rad_min, rad_max)
start1, interior1, end1 = pt1, pt1 + vm.Point3D(0,-rad1,rad1), pt1 + vm.Point3D(0,-2*rad1,0)
arc1 = vme.Arc3D(start1, interior1, end1)

pt2 = vm.Point3D.random(mini, maxi, mini, maxi, mini, maxi)
rad2 = random.randint(rad_min, rad_max)
start2, interior2, end2 = pt2, pt2 + vm.Point3D(-2*rad2,0,0), pt2 + vm.Point3D(-rad2,rad2,0)

# TODO testcase to make robust
try:
    arc2 = vme.Arc3D(start2, interior2, end2)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    pt1.plot(ax=ax)
    arc1.plot(ax=ax)
    end1.plot(ax=ax)
    interior1.plot(ax=ax)
    
    pt2.plot(ax=ax, color='r')
    interior2.plot(ax=ax, color='r')
    end2.plot(ax=ax,color='r')
    arc2.plot(ax=ax)
    
    pta1, pta2 = arc1.minimum_distance_points_arc(arc2)
    pta1.plot(ax=ax, color='m')
    pta2.plot(ax=ax, color='m')
    
    d_min = arc1.minimum_distance(arc2)
    print(d_min)
    
except ZeroDivisionError:
    pass