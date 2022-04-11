#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm
import volmdlr.faces as vmf

# 2D
control_points = [vm.Point3D(0, 0, 0),
                  vm.Point3D(0.1, 0.02, 0),
                  vm.Point3D(0.2, 0.02, 0),
                  vm.Point3D(0, 0, 0.15),
                  vm.Point3D(0.1, 0.02, 0.15),
                  vm.Point3D(0.2, 0.02, 0.15),
                  vm.Point3D(0, 0, 0.3),
                  vm.Point3D(0.1, 0.021, 0.3),
                  vm.Point3D(0.2, 0.022, 0.3),
                  ]


surface = vmf.BSplineSurface3D(degree_u=2,
                               degree_v=2,
                               control_points=control_points,
                               nb_u=3,
                               nb_v=3,
                               u_multiplicities=[1, 2, 2, 1],
                               v_multiplicities=[1, 2, 2, 1],
                               u_knots=[0.1, 0.3, 0.5, 0.7],
                               v_knots=[0.1, 0.3, 0.5, 0.7])

face = surface.rectangular_cut(0, 1, 0, 1)
face._check_platform()

# face.babylonjs()
                 # degree=2,
                 #               control_points=points2d,
                 #               knot_multiplicities=[1, 2, 2, 1],
                 #               knots=[0.1, 0.3, 0.5, 0.7])
