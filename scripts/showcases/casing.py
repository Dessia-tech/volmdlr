#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to show volmdlr capabilities on extrusions
"""

import math

import plot_data

from volmdlr.models import casing
from volmdlr.primitives3d import Block

casing._check_platform()

bottom, sides, belt = casing.primitives

bbox = casing.bounding_box

assert bbox.zmin == -0.005
assert math.isclose(bbox.xmax, 0.34067, abs_tol=1e-5)

box = Block.from_bounding_box(bbox)
box.alpha = 0.3
casing.primitives.append(box)
casing.babylonjs()

casing.to_step('casing')
casing.to_stl('casing')


contour = belt.outer_contour2d.plot_data()
primitive_group = plot_data.PrimitiveGroup(primitives=[contour])
plot_data.plot_canvas(plot_data_object=primitive_group, local=True)
