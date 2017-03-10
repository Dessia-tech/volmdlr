#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 17:02:36 2017

@author: steven
"""

import volmdlr as vm
import volmdlr.primitives2D as primitives2D

p1=vm.Point2D((0,0))
p2=vm.Point2D((1,0))
p3=vm.Point2D((2,1))
p4=vm.Point2D((1,0.5))
p5=vm.Point2D((-1,0.1))

p=primitives2D.RoundedLines2D([p1,p2,p3,p4,p5],{2:0.1},True)

c=vm.Contour2D([p])
print(c.Area())
p.MPLPlot()