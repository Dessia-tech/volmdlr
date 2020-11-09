#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Demo of cylinders

"""

import volmdlr as vm
import volmdlr.primitives3d as primitives3d

cylinder1 = primitives3d.Cylinder(vm.O3D, vm.X3D,
                                0.03, 0.02, 'cylinder1')
cylinder2 = primitives3d.HollowCylinder(vm.Y3D, vm.X3D,
                                      0.02, 0.06, 0.03,'cylinder2')

model=vm.core.VolumeModel([cylinder1, cylinder2])


model.babylonjs(debug=True)
