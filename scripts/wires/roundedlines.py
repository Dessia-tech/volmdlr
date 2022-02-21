#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 10:50:18 2018

@author: steven
"""

import math
import numpy as np
import volmdlr as vm
import volmdlr.primitives2d as primitives2d
import volmdlr.primitives3d as primitives3d


# Ajout commentaire juste pour tester
# =============================================================================
#  2D version
# =============================================================================

#p1 = vm.Point2D((0, 0))
#p11 = vm.Point2D((0.5, 0))
#p2 = vm.Point2D((1, 0))
#p3 = vm.Point2D((1, 1))
#p4 = vm.Point2D((2.1, 3.23))
#p5 = vm.Point2D((0, 1.23))
#p4 = vm.Point2D((0,1))
#p5 = vm.Point2D((0,2))

p0 = vm.Point2D(0.2, -0.5)

p1 = vm.Point2D(0  , 0  )
p2 = vm.Point2D(0.5, 0  )
p3 = vm.Point2D(1  , 0  )
p4 = vm.Point2D(1.5, 0  )

p5 = vm.Point2D(2  , 0.4)
p6 = vm.Point2D(2.5, 0.4)
p7 = vm.Point2D(3  , 0.4)
p8 = vm.Point2D(3.5, 0.2)





rl2D_o = primitives2d.OpenedRoundedLineSegments2D(
    [p0, p1, p2, p3, p4, p5, p6, p7, p8], {},
#                                        {2:0.3, 4:0.1, 3:0.1},
                                        adapt_radius=True)



rl2D_o2 = rl2D_o.offset_lines([2], -1.25)
ax= rl2D_o.plot()
rl2D_o2.plot(ax=ax)


rl2D_c = primitives2d.ClosedRoundedLineSegments2D([p0, p1, p2, p3, p4, p5, p6, p7, p8], {},
                                       # {0:1, 1:0.05, 2:0.05, 3:1},
                                        adapt_radius=True)
rl2D_c2 = rl2D_c.offset_lines([2], 0.2)
ax2 = rl2D_c.plot()
rl2D_c2.plot(ax=ax2)


com = rl2D_c2.center_of_mass()
cut_line = vm.edges.Line2D(com, com+ vm.Point2D.random(-1, 1, 0, 1))
ax3 = rl2D_c2.plot()
cut_line.plot(ax=ax3, color='red')
com.plot(color='b', ax=ax3)

cutted_contours = rl2D_c2.cut_by_line(cut_line)
list_centers_of_mass = [contour.center_of_mass() for contour in cutted_contours]

u1 = cut_line.point_projection(
    list_centers_of_mass[0])[0] - list_centers_of_mass[0]
u1.normalize()
sign = -1
for contour in cutted_contours:
    r, g, b = np.random.random(), np.random.random(), np.random.random()
    color = (r, g, b)
    contour.translation(sign*0.05*u1).plot(ax=ax3, color=color)
    sign *= -1

# assert math.isclose(cutted_contours[0].area()+cutted_contours[1].area(),
#                     rl2D_c2.area(), abs_tol=1e-12)


# ax4 = rl2D_c2.plot()
# mesh = rl2D_c2.grid_triangulation(25, 10)
# mesh.plot(ax=ax4)


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