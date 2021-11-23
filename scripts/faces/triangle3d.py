#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""


"""

import volmdlr as vm
import volmdlr.core as vmc
import volmdlr.faces as faces

p1 = vm.Point3D.random(0, 0.1, 0, 0.1, -0.1, 0)
p2 = vm.Point3D.random(0, 0.1, 0, 0.1, -0.1, 0)
p3 = vm.Point3D.random(0, 0.1, 0, 0.1, -0.1, 0)


tri = faces.Triangle3D(p1, p2, p3)
tri.babylonjs()
# tri.subdescription()

# vmc.VolumeModel([faces.OpenShell3D([tri])])