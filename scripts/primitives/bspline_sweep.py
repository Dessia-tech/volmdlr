#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A sweep using a bspline.
"""

import matplotlib.pyplot as plt

import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.wires as vmw
from volmdlr import primitives3d, curves
from volmdlr.core import EdgeStyle

degree = 5
control_points = [vm.Point3D(0, 0, 0),
                  vm.Point3D(0.3, 0.2, 0.1),
                  vm.Point3D(0.5, -0.1, 0.4),
                  vm.Point3D(0.5, -0.4, 0.0),
                  vm.Point3D(-0.1, -0.2, -0.3),
                  vm.Point3D(-0.3, 0.4, 0.1)]
knots = [0.0, 1.0]
knot_multiplicities = [6, 6]
weights = None  # [1, 2, 1, 2, 1, 2]
bspline_curve3d = vme.BSplineCurve3D(degree=degree,
                                     control_points=control_points,
                                     knot_multiplicities=knot_multiplicities,
                                     knots=knots,
                                     weights=weights,
                                     periodic=False,
                                     name='B Spline Curve 3D 1')

circle = curves.Circle2D(vm.O2D, 0.015)
contour = vmw.Contour2D(circle.split_at_abscissa(circle.length()*.5))

# rl = primitives3d.OpenRoundedLineSegments3D(points, radius, adapt_radius=True, name='wire')


sweep = primitives3d.Sweep(contour, vmw.Wire3D([bspline_curve3d]), name='Random pipe')

model = vm.core.VolumeModel([sweep])
model._check_platform()
model.babylonjs()

model.to_step('bspline_sweep.step')
