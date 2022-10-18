#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Demo of cylinders

"""

import volmdlr as vm
import volmdlr.primitives3d as primitives3d

cylinder1 = primitives3d.Cylinder(vm.O3D, vm.X3D,
                                0.03, 0.02, name='cylinder1')
cyl2_center = vm.Point3D(0.1, 0, 0)
cylinder2 = primitives3d.HollowCylinder(cyl2_center, vm.X3D,
                                      0.02, 0.06, 0.03, name='cylinder2')

model=vm.core.VolumeModel([cylinder1, cylinder2])


model.babylonjs(debug=True)

point1 = vm.Point3D(0.1, 0.3, -1.3)
point2 = vm.Point3D(-0.1, 0.2, 1)
cyl3 = primitives3d.Cylinder.from_extremal_points(point1, point2, 0.2)
cyl4 = primitives3d.HollowCylinder.from_extremal_points(point1, point2, 0.1, 0.2)

# Uncomment when Cylindrical surface implements linesegment intersections
# assert not cylinder2.point_belongs(cyl2_center)