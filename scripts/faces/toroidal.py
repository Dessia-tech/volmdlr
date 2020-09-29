#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 12:40:31 2020

@author: masfaraud
"""

import volmdlr as vm

# u = vm.Vector3D.random(0, 1, 0, 1, 0, 1)
# u.Normalize()
# v = u.RandomUnitNormalVector()
# w = u.Cross(v)

R = 0.5
r = 0.03
ts = vm.ToroidalSurface3D(vm.OXYZ, R, r)

tf = ts.rectangular_cut(0, 0.3, 0., 0.8)
tf.outer_contour2d.MPLPlot()
tf.babylonjs(debug=True, use_cdn=False)