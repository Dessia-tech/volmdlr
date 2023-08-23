#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  7 16:38:35 2023

@author: masfaraud
"""

import random
import volmdlr
import volmdlr.cloud

random.seed(561)
points = [volmdlr.Point3D.random(0.1, 0.3, -0.2, 0.37, -0.4, -0.05) for _ in range(500)]
# points.extend([volmdlr.Point3D.random(0.15, 0.25, -0.17, 0.05, -0.1, 0.05) for _ in range(150)])


cloud3d = volmdlr.cloud.PointCloud3D(points)
cloud3d.plot()

wrapping = cloud3d.to_shell()

# wrapping.babylonjs()

ax = None
for tri in wrapping.faces[1:40]:
    ax = tri.plot(ax=ax)
    normal = 0.03*tri.normal()
    
    # print('NORMAL', tri.normal())
    normal.plot(ax=ax, starting_point=tri.middle())
    
