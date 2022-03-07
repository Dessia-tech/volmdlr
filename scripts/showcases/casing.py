#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script checking offset and Curvilinear absissa of roundedline2D
"""

import plot_data
from volmdlr.models import casing


casing._check_platform()

bottom, sides, belt = casing.primitives

ax = belt.outer_contour2d.plot()
# outer_contour.plot(a)

ax = belt.outer_contour3d.plot()
l = belt.outer_contour3d.length()
for i in range(100):
    p = belt.outer_contour3d.point_at_abscissa(i*l/100)
    p.plot(ax=ax)

ax = belt.outer_contour3d.plot()


casing.babylonjs()
casing.to_step('casing')
casing.to_stl('casing')


contour = belt.outer_contour2d.plot_data()
primitive_group = plot_data.PrimitiveGroup(primitives=[contour])
plot_data.plot_canvas(plot_data_object=primitive_group, debug_mode=True)


