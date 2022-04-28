#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 18:06:45 2018

@author: steven
"""

import volmdlr
import volmdlr.edges


# Random arc
i = volmdlr.Point2D.random(-1,1,-1,1)
e = volmdlr.Point2D.random(-1,1,-1,1)
s = volmdlr.Point2D.random(-1,1,-1,1)


a = volmdlr.edges.Arc2D(s, i, e)
ax = a.plot()

for p in a.polygon_points(10):
    p.plot(ax=ax)
    
s.plot(ax=ax, color='r')
e.plot(ax=ax, color='g')
i.plot(ax=ax, color='b')