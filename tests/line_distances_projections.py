#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Distance to lines and projection debug

"""

import numpy as npy

import volmdlr as vm



p1 = vm.Point2D(2*npy.random.random(2)-1)
p2 = vm.Point2D(2*npy.random.random(2)-1)

line = vm.Line2D(p1, p2)

line_segment = vm.LineSegment2D(p1, p2)

for i in range(100):

    point = vm.Point2D(4*npy.random.random(2)-2)
    
    point_projection_line = line.PointProjection(point)
    point_projection_line_segment = line_segment.PointProjection(point)

    assert point_projection_line.PointDistance(point) == line.PointDistance(point)
    assert point_projection_line_segment.PointDistance(point) == line_segment.PointDistance(point)

#fig, ax = line.MPLPlot()
#point.MPLPlot(ax=ax)
#point_projection_line.MPLPlot(ax=ax)
#point_projection_line_segment.MPLPlot(ax=ax)