#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Surface2D examples
"""

import volmdlr
import volmdlr.primitives2d as p2d
import volmdlr.faces as faces
import volmdlr.wires as wires

l = 0.12
e1 = 0.015
e2 = 0.015
h1 = 0.035
h2 = 0.052
h3 = 0.018
r = 0.5*e2


p1 = volmdlr.Point2D(0, 0)
p2 = volmdlr.Point2D(l, 0)
p3 = volmdlr.Point2D(l, h1)
p4 = volmdlr.Point2D(e1, h1)
p5 = volmdlr.Point2D(0., h2+h3)
p6 = volmdlr.Point2D(-2.5*e2, h2+h3)

pc = volmdlr.Point2D(-e2, h2)

contour = volmdlr.primitives2d.ClosedRoundedLineSegments2D([p1, p2,
                                                            p3, p4,
                                                            p5, p6],
                                                            {0: r,
                                                             1: r,
                                                             4: r,
                                                             5: r},
                                                           adapt_radius=True)

# contour.plot()
hole = wires.Circle2D(pc, 0.5*r)
surface = faces.Surface2D(contour, [hole])
surface.plot()
contours = contour.split_regularly(10)
ax =  contours[0].plot()
for c in contours[1:]:
    c.plot(ax=ax)
mesh = surface.triangulation()
ax = surface.plot()
mesh.plot(ax=ax)