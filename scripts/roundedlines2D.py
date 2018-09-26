#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 19:30:39 2018

Script checking offset and Curvilinear absissa of roundedline2D
"""

import volmdlr as vm
import volmdlr.primitives2D as primitives2D

p1 = vm.Point2D((0, 0))
p2 = vm.Point2D((0.3, 0))
p3 = vm.Point2D((0.33, 0.22))
p4 = vm.Point2D((0.25, 0.08))
p5 = vm.Point2D((0.05, 0.20))

c1 = primitives2D.RoundedLineSegments2D([p1, p2, p3, p4, p5], {0: 0.01, 1: 0.01, 2: 0.015, 3: 0.010, 4:0.008}, True)
c2 = c1.Offset(0.010)
f,a = c1.MPLPlot()
c2.MPLPlot(a)

l = c1.Length()
n = 250
for i in range(n):
    s = i * l/n
#    print(s)
    p = c1.PointAtCurvilinearAbscissa(s)
    p.MPLPlot(a)