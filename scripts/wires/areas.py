#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 13:55:59 2017

@author: steven
"""

# Debug of area and inertia

import volmdlr as vm

r=0.05

p1=vm.Point2D((0,0))
p2=vm.Point2D((0,-r))
p3=vm.Point2D((r,0))
p4=vm.Point2D((0,r))

c1=vm.Arc2D(p2,p3,p4)

print(c1.SecondMomentArea(p1))

c2=vm.Circle2D(p1,r)

print(c2.SecondMomentArea(p1))
