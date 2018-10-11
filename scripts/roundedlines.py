#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 10:50:18 2018

@author: steven
"""

import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import volmdlr.primitives3D as primitives3D


# =============================================================================
#  2D version
# =============================================================================

p1 = vm.Point2D((0, 0))
p2 = vm.Point2D((1, 0))
p3 = vm.Point2D((1, 1))
p4 = vm.Point2D((2.1, 3.23))
p5 = vm.Point2D((0, 1.23))


rl2D = primitives2D.RoundedLineSegments2D([p1, p2, p3, p4, p5],
                                        {0: 0.7, 1:0.3, 2:0.1, 3:2},
                                        closed=True, adapt_radius=True)
rl2D.MPLPlot()


# =============================================================================
#  3D Version
# =============================================================================

p1 = vm.Point3D((0, 0, 0))
p2 = vm.Point3D((1, 0, 0))
p3 = vm.Point3D((1, 1, 1.2))
p4 = vm.Point3D((2.1, 3.23, 0.3))
p5 = vm.Point3D((0, 1.23,4.6))


rl3D = primitives3D.RoundedLineSegments3D([p1, p2, p3, p4, p5],
                                        {2:4},
                                        closed=True, adapt_radius=True)
rl3D.MPLPlot()
m = vm.VolumeModel([('rl3D', [rl3D])])
m.FreeCADExport('rl3D')