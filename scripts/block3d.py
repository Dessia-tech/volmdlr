#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 11 23:39:42 2020

@author: Pierrem
"""
import math
import volmdlr as vm
import volmdlr.step as vm_step
import volmdlr.primitives3d as primitives3d
resolution = 0.0010

box = primitives3d.Block(
    vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.3, 0, 0),
               vm.Vector3D(0, 0.3, 0), vm.Vector3D(0, 0, 0.3)),
    alpha=0.6)

box_red = primitives3d.Block(
    vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.4, 0, 0),
               vm.Vector3D(0, 0.4, 0), vm.Vector3D(0, 0, 0.4)),
    color=(0.2, 1, 0.4), alpha=0.6)

p1_ray = vm.Point3D(-0.15, -0.15, -0.15)
p2_ray = vm.Point3D(0.009855980224206917, 0.6250574317556334, -0.1407142090413507)
# p1_ray = vm.Point3D(-0.15, -0.12999999999999992, 0.15)
# p2_ray = vm.Point3D(0.09377883804318171, 0.17764785706502192, 0.19256693676483136)
ray = vm.edges.LineSegment3D(p1_ray, p2_ray)


ax = ray.plot(color='b')
p1_ray.plot(ax=ax, color='b')
p2_ray.plot(ax=ax, color='b')
box_red.plot(ax=ax, color='r')
for face, inter_points in box_red.linesegment_intersections(ray):
    # print('ip', inter_point)
    face.plot(ax=ax, color='b')
    for inter_point in inter_points:
        inter_point.plot(ax=ax, color='r')


box_red.color = (1, 0.1, 0.1)
box_red.name = 'box_red'

box_green = box.frame_mapping(vm.Frame3D(vm.Point3D(0, 0.8, 0), vm.Vector3D(1, 0, 0),
                         vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'new', copy=True)

box_green.color = (0.1, 1, 0.1)
box_green.name = 'box_green'


box_blue = box.frame_mapping(vm.Frame3D(vm.Point3D(0, 0.2, 0), vm.Vector3D(1, 0, 0),
                         vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'old', copy=True)
box_blue.color = (0.1, 0.1, 1)
box_blue.name = 'box_blue'

assert box.faces[0] == box.faces[0]
print(box.distance_to_shell(box_red, resolution))
print(box_green.shell_intersection(box_blue, resolution))
print(box_green.intersection_internal_aabb_volume(box_blue, resolution))
print(box_green.intersection_external_aabb_volume(box_blue, resolution))
model = vm.core.VolumeModel([box, box_red, box_green, box_blue])
model.babylonjs(debug=True)

assert box.is_inside_shell(box_red, resolution) == True
assert box_red.is_inside_shell(box, resolution) == False

assert box.is_inside_shell(box_green, resolution) == False
assert box_green.is_inside_shell(box, resolution) == False

assert box.is_inside_shell(box_blue, resolution) == False
assert box_blue.is_inside_shell(box, resolution) == False

model = vm.core.VolumeModel([box_red])
model.to_step('block.step')

step = vm_step.Step('block.step')
model2 = step.to_volume_model()

# box = primitives3d.Block(
#     vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.3, 0, 0),
#                vm.Vector3D(0, 0.3, 0), vm.Vector3D(0, 0, 0.3)),
#     alpha=0.6)
# box_red = primitives3d.Block(
#     vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.1, 0, 0),
#                vm.Vector3D(0, 0.1, 0), vm.Vector3D(0, 0, 0.1)),
#     color=(0.2, 1, 0.4), alpha=0.6)
#
# for i in range(30):
#     print('----NEW STEP----', box_red.is_inside_shell(box, resolution))
#     print('distance_to_shell', box.distance_to_shell(box_red, resolution))
#     # print('shell_intersection', box.shell_intersection(box_red, resolution))
#     print('volume', box_red.bounding_box.volume(), box.bounding_box.volume())
#     print('intersection_internal_aabb_volume', box.intersection_internal_aabb_volume(box_red, resolution), box_red.intersection_internal_aabb_volume(box, resolution))
#     print('intersection_external_aabb_volume', box.intersection_external_aabb_volume(box_red, resolution), box_red.intersection_external_aabb_volume(box, resolution))
#     # print('is_inside_shell', box.is_inside_shell(box_red, resolution), box_red.is_inside_shell(box, resolution))
#     if not box_red.is_inside_shell(box, resolution):
#         model = vm.core.VolumeModel([box, box_red])
#         model.babylonjs(debug=True)
#         raise
#     box_red = box_red.translation(vm.Vector3D(0.01, 0, 0), copy=True)
#
#
# model = vm.core.VolumeModel([box, box_red])
# model.babylonjs(debug=True)
