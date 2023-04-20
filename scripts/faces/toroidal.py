#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 12:40:31 2020

@author: masfaraud
"""

import volmdlr as vm
from volmdlr import faces, surfaces

# u = vm.Vector3D.random(0, 1, 0, 1, 0, 1)
# u.Normalize()
# v = u.RandomUnitNormalVector()
# w = u.Cross(v)

R = 0.2
r = 0.03
ts = surfaces.ToroidalSurface3D(vm.OXYZ, R, r)

tf = faces.ToroidalFace3D.from_surface_rectangular_cut(ts, 0, 0.6, 0., 1.3)
# tf.outer_contour2d.MPLPlot()
tf.babylonjs(debug=True, use_cdn=False)
