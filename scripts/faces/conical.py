#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import volmdlr as vm
import volmdlr.faces as faces

R = 0.32
alpha = 0.2
cs = vm.faces.ConicalSurface3D(vm.Frame3D(vm.Point3D.random(-0.1, 0.1, -0.1, 0.2, -0.2, 0.1),
                                          vm.X3D, vm.Y3D, vm.Z3D), alpha)

cf = cs.rectangular_cut(-0.01, 1.3, 0., 0.3)

cf.surface2d.plot()
cf.babylonjs(debug=True, use_cdn=False)