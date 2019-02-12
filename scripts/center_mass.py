#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 12 19:23:27 2019

@author: steven
"""

import volmdlr as vm
import math

r = 0.23

p1 = vm.Point2D((0, 0))
p2 = vm.Point2D((r, r))
p3 = vm.Point2D((2*r, 0))

a = vm.Arc2D(p1, p2, p3)
l = vm.LineSegment2D(p3, p1)

c = vm.Contour2D([a, l])
f, ax = c.MPLPlot()
cog = c.CenterOfMass()
cog.MPLPlot(ax)

assert math.isclose(a.radius, r)
assert math.isclose(c.Area(), math.pi*r**2/2.)
assert math.isclose(cog[0], r)
print(cog[1], 4*r/3.*math.pi)
assert math.isclose(cog[1], 4*r/3.*math.pi)