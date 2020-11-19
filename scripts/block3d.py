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
vol1 = vol.frame_mapping(vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.9, 0, 0),
                         vm.Vector3D(0, 0.9, 0), vm.Vector3D(0, 0, 0.9)), 'new', copy=True)
vol2 = vol.frame_mapping(vm.Frame3D(vm.Point3D(0, 0.8, 0), vm.Vector3D(1, 0, 0),
                         vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'new', copy=True)
vol3 = vol.frame_mapping(vm.Frame3D(vm.Point3D(0, 0.7, 0), vm.Vector3D(1, 0, 0),
                         vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'old', copy=True)

vol1.primitives.extend(vol.primitives)
vol1.primitives.extend(vol2.primitives)
vol1.babylonjs()

print(vol.primitives[0].faces[0]==vol.primitives[0].faces[0])
print(vol.primitives[0].is_inside_shell(vol1.primitives[0]))

print(vol.primitives[0].distance_to_shell(vol1.primitives[0]))
print(vol2.primitives[0].shell_intersection(vol3.primitives[0]))
print(vol2.primitives[0].intersection_internal_aabb_volume(vol3.primitives[0]))
print(vol2.primitives[0].intersection_external_aabb_volume(vol3.primitives[0]))
