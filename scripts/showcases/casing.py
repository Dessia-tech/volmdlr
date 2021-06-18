#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script checking offset and Curvilinear absissa of roundedline2D
"""

import volmdlr as vm
import volmdlr.primitives2d as primitives2d
import volmdlr.primitives3d as primitives3d
import matplotlib.pyplot as plt
import plot_data

thickness = 0.005
height = 0.080
screw_holes_diameter = 0.006
screw_holes_clearance = 0.003
n_screws = 25

p1 = vm.Point2D(0, 0)
p2 = vm.Point2D(0.3, 0)
p3 = vm.Point2D(0.33, 0.22)
p4 = vm.Point2D(0.2, 0.08)
p5 = vm.Point2D(0.16, 0.18)
p6 = vm.Point2D(0.05, 0.20)

inner_contour = primitives2d.ClosedRoundedLineSegments2D([p1, p2, p3, p4, p5, p6],
                                                         {0: 0.01, 1: 0.01, 2: 0.015, 3: 0.020, 4: 0.012, 5:0.008},
                                                         adapt_radius=True)

outer_contour = inner_contour.offset(-thickness)

a = inner_contour.plot()
outer_contour.plot(a)



sides = primitives3d.ExtrudedProfile(vm.O3D, vm.X3D, vm.Y3D,
                                      outer_contour, [inner_contour],
                                      (height-2*thickness) * vm.Z3D, name='sides')

bottom = primitives3d.ExtrudedProfile(vm.O3D, vm.X3D, vm.Y3D, outer_contour, [],
                                      -thickness * vm.Z3D, name='bottom')

screw_holes_rl = inner_contour.offset(-(thickness+screw_holes_clearance + 0.5 * screw_holes_diameter))
screw_holes = []
l = screw_holes_rl.length()
for i in range(n_screws):
    s = i * l/n_screws
    p = screw_holes_rl.point_at_abscissa(s)
    screw_holes.append(vm.wires.Circle2D(p, screw_holes_diameter*0.5))
# ###
# fig, ax = plt.subplots()
# [prim.MPLPlot(ax=ax) for prim in inner_contour.primitives]
# ###
belt_outer_contour = inner_contour.offset(-(2*screw_holes_clearance + screw_holes_diameter+thickness))
belt = primitives3d.ExtrudedProfile(vm.Z3D*(height - 2*thickness), vm.X3D, vm.Y3D,
                                      belt_outer_contour,
                                      [inner_contour]+screw_holes,
                                      thickness * vm.Z3D, name='belt')

ax = inner_contour.plot()
belt_outer_contour.plot(ax=ax)


ax = belt.outer_contour3d.plot()
l = belt.outer_contour3d.length()
for i in range(100):
    p = belt.outer_contour3d.point_at_abscissa(i*l/100)
    p.plot(ax=ax)

ax = belt.outer_contour3d.plot()
# l = belt.outer_contour3d.Length()
# for i in range(100):



model = vm.core.VolumeModel([bottom, sides, belt], name='Casing')
model.babylonjs()
model.to_step('casing')
model.to_stl('casing')

#
# model = vm.VolumeModel([sides])
# model.babylonjs('sides')
#
# model = vm.VolumeModel([belt])
# model.babylonjs('belt')

# dict_ = model.to_dict()
# model2 = vm.VolumeModel.dict_to_object(dict_)

#casing = vm.primitives3d.Fuse([bottom, sides, belt], 'Lower Casing')
#model = vm.VolumeModel([casing])
##model.FreeCADExport('casing')
#model.BabylonShow()

contour = screw_holes_rl.plot_data()
primitive_group = plot_data.PrimitiveGroup(primitives=[contour])
plot_data.plot_canvas(plot_data_object=primitive_group, debug_mode=True)

