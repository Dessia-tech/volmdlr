#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Testing core module functions
"""

import math
import numpy as npy
import volmdlr as vm

v2D_1 = vm.Vector2D(3*npy.random.random(2)-1.5)
v2D_2 = vm.Vector2D(3*npy.random.random(2)-1.5)

p2D_1 = vm.Point2D(3*npy.random.random(2)-1.5)
p2D_2 = vm.Point2D(3*npy.random.random(2)-1.5)


v3D_1 = vm.Vector3D(3*npy.random.random(3)-1.5)
v3D_2 = vm.Vector3D(3*npy.random.random(3)-1.5)

p3D_1 = vm.Point3D(3*npy.random.random(3)-1.5)
p3D_2 = vm.Point3D(3*npy.random.random(3)-1.5)


# Testing if normalized vector has norm ==1 and is still colinear to original vector
v2D_1_normalized = v2D_1.copy()
v2D_1_normalized.Normalize()
assert math.isclose(v2D_1_normalized.Norm(), 1, abs_tol=1e-9)
assert math.isclose(v2D_1_normalized.Dot(v2D_1),v2D_1.Norm(), abs_tol=1e-9)

# Testing normal vector
normal_v2D_1 = v2D_1.NormalVector()
assert  math.isclose(normal_v2D_1.Dot(v2D_1), 0, abs_tol=1e-9)
normal_unit_v2D_1 = v2D_1.NormalVector(unit=True)
assert math.isclose(normal_unit_v2D_1.Norm(), 1., abs_tol=1e-9)

# Testing Cross
