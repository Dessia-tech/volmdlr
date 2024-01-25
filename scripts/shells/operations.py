#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Goal: call shells methods (point_belongs, volume)
"""

import time

import volmdlr as vm
from volmdlr.models import casing

bottom, sides, belt = casing.primitives

t1 = time.time()
assert not sides.point_inside(vm.Point3D(0.1, 0.1, 0.1))
t1 = time.time() - t1
