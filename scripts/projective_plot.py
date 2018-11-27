#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Drawing 3D primitives in 2D
"""

import volmdlr as vm
import volmdlr.primitives3D as primitives3D



b = primitives3D.Block(vm.Frame3D(vm.Point3D((1, 2.3, 4)), vm.Vector3D((0.5,0.3, -0.1)), vm.y3D, vm.z3D))

b.MPLPlot(vm.y3D, vm.z3D)