#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Minimum distance line between two lines
"""

import volmdlr as vm
import numpy as npy

p11 = vm.Point3D(npy.random.random(3))
p12 = vm.Point3D(npy.random.random(3))
p21 = vm.Point3D(npy.random.random(3))
p22 = vm.Point3D(npy.random.random(3))

l1 = vm.Line3D(p11, p12)
l2 = vm.Line3D(p21, p22)

pmd1, pmd2 = l1.MinimumDistancePoints(l2)

u = p12 - p11 # vector of line1
v = p22 - p21 # vector of line2
w = pmd2 - pmd1

print(u.Dot(w), v.Dot(w))

m = vm.VolumeModel([('', [l1, l2, pmd1, pmd2])])
m.MPLPlot()