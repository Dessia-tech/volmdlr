#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  8 10:58:19 2017

@author: steven
"""

import volmdlr as vm

p1=vm.Point2D((1,1.45))
p2=vm.Point2D((0.4,0.1))
l1=vm.Line2D(p1,p2)

p3=vm.Point2D((-1,0))
p4=vm.Point2D((1,-0.5))
l2=vm.Line2D(p3,p4)

p5=vm.Point2D.LinesIntersection(l1,l2)
p6=vm.Point2D.MiddlePoint(p1,p3)

c=vm.CompositePrimitive2D([l1,l2,p5,p6])
c.MPLPlot()

