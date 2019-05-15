#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 10:50:18 2018

@author: steven
"""

import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import volmdlr.primitives3D as primitives3D

# Ajout commentaire juste pour tester
# =============================================================================
#  2D version
# =============================================================================

p1 = vm.Point2D((0, 0))
p11 = vm.Point2D((0.5, 0))
p2 = vm.Point2D((1, 0))
p3 = vm.Point2D((1, 1))
p4 = vm.Point2D((2.1, 3.23))
p5 = vm.Point2D((0, 1.23))
#p4 = vm.Point2D((0,1))
#p5 = vm.Point2D((0,2))



rl2D_o = primitives2D.RoundedLineSegments2D([p1, p11, p2, p3, p4, p5],
                                        {2:0.3, 4:0.1, 3:0.1},
                                        closed=False, adapt_radius=True)



rl2D_o2 = rl2D_o.Offset(0.1)
f, ax= rl2D_o.MPLPlot()
rl2D_o2.MPLPlot(ax=ax, style='r')


rl2D_c = primitives2D.RoundedLineSegments2D([p1, p2, p3, p4, p5],
                                        {0: 1, 1:0.05, 2:0.05, 3:1},
                                        closed=True, adapt_radius=True)
rl2D_c2 = rl2D_c.Offset(0.2)
f2, ax2 = rl2D_c.MPLPlot()
rl2D_c2.MPLPlot(ax=ax2, style='r')


# =============================================================================
#  3D Version
# =============================================================================
#
#p1 = vm.Point3D((0, 0, 0))
#p2 = vm.Point3D((1, 0, 0))
#p3 = vm.Point3D((1, 1, 1.2))
#p4 = vm.Point3D((2.1, 3.23, 0.3))
#p5 = vm.Point3D((0, 1.23,4.6))
#
#
#rl3D = primitives3D.RoundedLineSegments3D([p1, p2, p3, p4, p5],
#                                        {2:4},
#                                        closed=True, adapt_radius=True)
#rl3D.MPLPlot()
#m = vm.VolumeModel([('rl3D', [rl3D])])
#m.FreeCADExport('rl3D')

#{1:0.3, 2:0.1, 3:0.5}