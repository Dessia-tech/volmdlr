#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:34:09 2021

@author: dasilva
"""
import volmdlr as vm
import volmdlr.edges as edges
import matplotlib.pyplot as plt
import volmdlr.primitives3d as primitives3D
# p1, p2, p3, p4 = vm.Point3D(-0.086, 0.11399999999999999, 0.21),vm.Point3D(-0.086, 0.114, 0.17), vm.Point3D(-0.086, -0.126, 0.2042857142857143), vm.Point3D(-0.086, 0.124, 0.2042857142857143)
# l1=edges.Line3D(p1,p2)
# l2=edges.Line3D(p3,p4)
# print(l1.intersection(l2))
# fig = plt.figure()
# ax = fig.add_subplot(111,projection='3d')
# l1.plot(ax=ax)
# l2.plot(ax=ax)
# p1.plot(ax=ax)
# p2.plot(ax=ax)
# p3.plot(ax=ax)
# p4.plot(ax=ax)

# cyl1 = primitives3D.Cylinder(vm.Point3D(0, 0,0), vm.Vector3D(1,0,0), 0.055, 0.15)
# cyl2 = cyl1.copy()
# # cyl3 = cyl2.frame_mapping(vm.Frame3D(vm.Point3D(0.3,0,0), vm.Vector3D(1,0,0), vm.Vector3D(0,1,0), vm.Vector3D(0,0,1)), 'old')
# cyl2.frame_mapping(vm.Frame3D(vm.Point3D(0.3,0,0), vm.Vector3D(1,0,0), vm.Vector3D(0,1,0), vm.Vector3D(0,0,1)), 'old', False)
# vol = vm.core.VolumeModel([cyl1, cyl2])
# vol.babylonjs(debug = True)


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
new_box = vm.faces.ClosedShell3D.unions(box_red, box_blue)
new_box.color = (1, 0.1, 0.1)
new_box.alpha = 0.6
# for face in new_box.face:
#     face.color = 
vm.core.VolumeModel([new_box]).babylonjs()

