#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 18:06:45 2018

@author: steven
"""

import volmdlr as vm






# Random arc
i = vm.Point2D.random(-1,1,-1,1)
e = vm.Point2D.random(-1,1,-1,1)
s = vm.Point2D.random(-1,1,-1,1)

# s = vm.Point2D((-0.9519892596753585, -0.3892760733608891))
# e = vm.Point2D((-0.6728670157346539, -0.7397140760705987))
# i = vm.Point2D((-0.5843044181756767, 0.47828817571533633))


a = vm.Arc2D(s, i, e)
ax = a.MPLPlot()



for p in a.tessellation_points():
    p.MPLPlot(ax=ax)
    
s.MPLPlot(ax=ax, color='r')
e.MPLPlot(ax=ax, color='g')
i.MPLPlot(ax=ax, color='b')