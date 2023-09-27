#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test case representing a casing. Use extrusion features.
"""

import volmdlr as vm
import volmdlr.wires
from volmdlr import primitives2d, primitives3d, curves, edges

THICKNESS = 0.005
HEIGHT = 0.080
SCREW_HOLES_DIAMETER = 0.006
SCREW_HOLES_CLEARANCE = 0.003
N_SCREWS = 25

p1 = vm.Point2D(0, 0)
p2 = vm.Point2D(0.3, 0)
p3 = vm.Point2D(0.33, 0.22)
p4 = vm.Point2D(0.2, 0.08)
p5 = vm.Point2D(0.16, 0.18)
p6 = vm.Point2D(0.05, 0.20)

inner_contour = primitives2d.ClosedRoundedLineSegments2D([p1, p2, p3, p4, p5, p6],
                                                         {0: 0.01, 1: 0.01, 2: 0.015, 3: 0.020, 4: 0.012, 5: 0.008},
                                                         adapt_radius=True)

outer_contour = inner_contour.offset(-THICKNESS)


sides = primitives3d.ExtrudedProfile(vm.OXYZ,
                                     outer_contour, [inner_contour],
                                     HEIGHT-2*THICKNESS, name='sides')

bottom = primitives3d.ExtrudedProfile(vm.OXYZ, outer_contour, [],
                                      -THICKNESS, name='bottom')

screw_holes_rl = inner_contour.offset(-(THICKNESS+SCREW_HOLES_CLEARANCE + 0.5 * SCREW_HOLES_DIAMETER))
screw_holes = []
length = screw_holes_rl.length()
for i in range(N_SCREWS):
    s = i * length / N_SCREWS
    p = screw_holes_rl.point_at_abscissa(s)
    circle = curves.Circle2D(p, SCREW_HOLES_DIAMETER*0.5)
    screw_holes.append(volmdlr.wires.Contour2D([edges.FullArc2D.from_curve(circle)]))

belt_outer_contour = inner_contour.offset(-(2*SCREW_HOLES_CLEARANCE + SCREW_HOLES_DIAMETER+THICKNESS))
belt = primitives3d.ExtrudedProfile(volmdlr.Frame3D(vm.Point3D(0, 0, 1) * (HEIGHT - 2*THICKNESS),
                                                    vm.X3D, vm.Y3D, vm.Z3D),
                                    belt_outer_contour,
                                    [inner_contour]+screw_holes,
                                    THICKNESS, name='belt')

casing = vm.core.VolumeModel([bottom, sides, belt], name='Casing')
