#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 19:30:39 2018

Script checking offset and Curvilinear absissa of roundedline2D
"""

import volmdlr as vm
import volmdlr.primitives2D as primitives2D
import volmdlr.primitives3D as primitives3D

thickness = 0.005
height = 0.080
screw_holes_diameter = 0.006
screw_holes_clearance = 0.003
n_screws = 25

p1 = vm.Point2D((0, 0))
p2 = vm.Point2D((0.3, 0))
p3 = vm.Point2D((0.33, 0.22))
p4 = vm.Point2D((0.2, 0.08))
p5 = vm.Point2D((0.16, 0.18))
p6 = vm.Point2D((0.05, 0.20))

inner_contour = primitives2D.ClosedRoundedLineSegments2D([p1, p2, p3, p4, p5, p6], {0: 0.01, 1: 0.01, 2: 0.015, 3: 0.010, 4: 0.012, 5:0.008}, True)
outer_contour = inner_contour.Offset(-thickness)
#inner_contour = vm.Contour2D([inner_rl])
#outer_contour = vm.Contour2D([outer_rl])

f, a = inner_contour.MPLPlot()
outer_contour.MPLPlot(a)

    

sides = primitives3D.ExtrudedProfile(vm.o3D, vm.x3D, vm.y3D,
                                     outer_contour, [inner_contour],
                                     (height-2*thickness) * vm.z3D, 'sides')
bottom = primitives3D.ExtrudedProfile(vm.o3D, vm.x3D, vm.y3D, outer_contour, [], -thickness * vm.z3D, 'bottom')

screw_holes_rl = inner_contour.Offset(-(thickness+screw_holes_clearance + 0.5 * screw_holes_diameter))
screw_holes = []
l = screw_holes_rl.Length()
for i in range(n_screws):
    s = i * l/n_screws
    p = screw_holes_rl.PointAtCurvilinearAbscissa(s)
    screw_holes.append(vm.Contour2D([vm.Circle2D(p, screw_holes_diameter*0.5)]))

belt_outer_contour = inner_contour.Offset(-(2*screw_holes_clearance + screw_holes_diameter+thickness))
#belt_outer_contour = vm.Contour2D([belt_outer_rl])
belt = primitives3D.ExtrudedProfile(vm.z3D*(height - 2*thickness), vm.x3D, vm.y3D,
                                      belt_outer_contour, [inner_contour]+screw_holes, thickness * vm.z3D, 'belt')

model = vm.VolumeModel([bottom, sides, belt])
model.BabylonShow('bottom')

#model = vm.VolumeModel([sides])
#model.BabylonShow('sides')
#
#model = vm.VolumeModel([belt])
#model.BabylonShow('belt')

model.to_dict()

#casing = vm.primitives3D.Fuse([bottom, sides, belt], 'Lower Casing')
#model = vm.VolumeModel([casing])
##model.FreeCADExport('casing')
#model.BabylonShow()