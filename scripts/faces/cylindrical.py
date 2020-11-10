#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 12:40:31 2020

@author: masfaraud
"""

import volmdlr
import volmdlr.faces

R = 0.32

ts = volmdlr.faces.CylindricalSurface3D(volmdlr.OXYZ, R)

tf = ts.rectangular_cut(-0.01, 1.3, -0.1, 0.3)
# tf.outer_contour2d.plot()
tf.babylonjs(debug=True, use_cdn=False)