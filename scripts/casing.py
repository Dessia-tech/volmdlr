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

inner_contour = primitives2D.RoundedLineSegments2D([p1, p2, p3, p4, p5, p6], {0: 0.01, 1: 0.01, 2: 0.015, 3: 0.010, 4: 0.012, 5:0.008}, True)
outer_contour = inner_contour.Offset(-thickness)
f, a = inner_contour.MPLPlot()
outer_contour.MPLPlot(a)

    

sides = primitives3D.ExtrudedProfile(vm.o3D, vm.x3D, vm.y3D,
                                     outer_contour, [inner_contour],
                                     (height-2*thickness) * vm.z3D, 'sides')
bottom = primitives3D.ExtrudedProfile(vm.o3D, vm.x3D, vm.y3D, outer_contour, [], -thickness * vm.z3D, 'bottom')

screw_holes_contour = inner_contour.Offset(-(thickness+screw_holes_clearance + 0.5 * screw_holes_diameter))
screw_holes = []
l = screw_holes_contour.Length()
for i in range(n_screws):
    s = i * l/n_screws
    p = screw_holes_contour.PointAtCurvilinearAbscissa(s)
    screw_holes.append(vm.Contour2D([vm.Circle2D(p, screw_holes_diameter*0.5)]))

belt_outer_contour = inner_contour.Offset(-(2*screw_holes_clearance + screw_holes_diameter+thickness))
belt = primitives3D.ExtrudedProfile(vm.z3D*(height - 2*thickness), vm.x3D, vm.y3D,
                                      belt_outer_contour, [inner_contour]+screw_holes, thickness * vm.z3D, 'belt')

casing = vm.primitives3D.Fuse([bottom, sides, belt], 'Lower Casing')
model = vm.VolumeModel([('casing', [casing])])
model.FreeCADExport('casing')