#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import volmdlr as vm
import volmdlr.faces as faces

R = 0.32
alpha = 0.3
cs = vm.faces.ConicalSurface3D(vm.OZXY, alpha)

cf = cs.rectangular_cut(-0.01, 1.3, 0., 0.3)
cf.surface2d.plot()
cf.babylonjs(debug=True, use_cdn=False)