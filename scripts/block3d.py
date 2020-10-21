#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 23:39:42 2020

@author: Pierrem
"""
import math
import volmdlr as vm
import volmdlr.primitives3d as primitives3d

bx0 = primitives3d.Block(
    vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.3, 0, 0), vm.Vector3D(0, 0.3, 0), vm.Vector3D(0, 0, 0.3)),
    color=(0.2, 1, 0.1), alpha=0.6)
vol = vm.core.VolumeModel([bx0])
# vol.babylonjs()
# vol1 = vol.copy()
vol1 = vol.frame_mapping(vm.Frame3D(vm.Point3D(0, 1, 0), vm.Vector3D(1, 0, 0),
                         vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'new', copy=True)
vol2 = vol.frame_mapping(vm.Frame3D(vm.Point3D(0, 0.5, 0), vm.Vector3D(1, 0, 0),
                         vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'old', copy=True)

vol1.primitives.extend(vol.primitives)
vol1.primitives.extend(vol2.primitives)
# vol1.babylonjs()

print(vol.primitives[0].faces[0]==vol.primitives[0].faces[0])