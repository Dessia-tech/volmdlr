#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Minimum distance line between two lines
"""

import volmdlr as vm
import volmdlr.primitives3d as primitives3d
#import numpy as npy

p11 = vm.Point3D.random(-1,1,-1,1,-1,1)
p12 = vm.Point3D.random(-1,1,-1,1,-1,1)
p21 = vm.Point3D.random(-1,1,-1,1,-1,1)
p22 = vm.Point3D.random(-1,1,-1,1,-1,1)

l1 = vm.edges.Line3D(p11, p12)
l2 = vm.edges.Line3D(p21, p22)

pmd1, pmd2 = l1.minimum_distance_points(l2)

u = p12 - p11 # vector of line1
v = p22 - p21 # vector of line2
w = pmd2 - pmd1

print(u.dot(w), v.dot(w))

m = vm.core.VolumeModel([('', [l1, l2, pmd1, pmd2])])

m.MPLPlot() #Problem

#m.mpl_plot() ?
#m.babylonjs() ? 

