#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import volmdlr as vm

R = 0.32
alpha = 0.3
cs = vm.ConicalSurface3D(vm.OXYZ, alpha)

cf = cs.rectangular_cut(-0.01, 1.3, 0., 0.3)
cf.outer_contour2d.MPLPlot()
cf.babylonjs(debug=True, use_cdn=False)