#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 18:06:45 2018

@author: steven
"""

import volmdlr as vm
import numpy as npy
i = vm.Point3D((1,0,0))
e = i.Rotation(vm.o3D, vm.z3D, 1)
s = i.Rotation(vm.o3D, vm.z3D, -3.5)

a = vm.Arc3D(s, i, e)
assert a.angle == 4.5


# Random arc
i = vm.Point3D(npy.random.random(3))
e = vm.Point3D(npy.random.random(3))
s = vm.Point3D(npy.random.random(3))

a = vm.Arc3D(s, i, e)
a.MPLPlot()