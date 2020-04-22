#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr.core as vm
from volmdlr.primitives2D import ClosedRoundedLineSegments2D
from volmdlr.primitives3D import RevolvedProfile

p1 = vm.Point2D((-0.08, 0.1755))
p2 = vm.Point2D((0., 0.1755))
p3 = vm.Point2D((0.0, 0.04))
p4 = vm.Point2D((0.01, 0.04))
p5 = vm.Point2D((0.01, 0.))
p6 = vm.Point2D((-0.06, 0.0))
p7 = vm.Point2D((-0.06, 0.04))
p8 = vm.Point2D((-0.08, 0.04))

points1 = [p1, p2, p3, p4, p5, p6, p7, p8]

c1 = vm.Polygon2D(points1)
c2 = ClosedRoundedLineSegments2D(points1, {0: 0.002, 1: 0.002, 3: 0.002}) 

# c1.MPLPlot(plot_points=True)
# c2.MPLPlot(plot_points=True)

p21 = vm.Point2D((0.1, 0.))
p22 = vm.Point2D((0.1, 0.1))
p23 = vm.Point2D((0.08, 0.1))
p24 = vm.Point2D((0.08, 0.2))
p25 = vm.Point2D((0., 0.2))
p26 = vm.Point2D((0., 0.05))
p27 = vm.Point2D((0.04, 0.05))
p28 = vm.Point2D((0.04, 0.))

points2 = [p21, p22, p23, p24, p25, p26, p27, p28]
c3 = vm.Polygon2D(points2)
c4 = ClosedRoundedLineSegments2D(points2, {0: 0.002, 1: 0.002, 3: 0.002})

# c3.MPLPlot(plot_points=True)
# c4.MPLPlot(plot_points=True)

profile1 = RevolvedProfile(vm.O3D, vm.Y3D, vm.Z3D, c1, vm.O3D, vm.Y3D)
profile2 = RevolvedProfile(0.5*vm.X3D, vm.Y3D, vm.Z3D, c2, 0.5*vm.X3D, vm.Y3D)
profile3 = RevolvedProfile(-0.5*vm.X3D, vm.Y3D, vm.Z3D, c1, -0.5*vm.X3D, vm.Y3D)
profile4 = RevolvedProfile(vm.X3D, vm.Y3D, vm.Z3D, c2, vm.X3D, vm.Y3D)

contour_dict = {'name': '',
 'package_version': None,
 'object_class': 'volmdlr.core.Contour2D',
 'primitives': [{'name': '',
   'package_version': None,
   'object_class': 'volmdlr.core.Arc2D',
   'interior': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0, 0.031875]},
   'start': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.0041097537330863606, 0.0345]},
   'end': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.0041097537330863606, 0.0345]},
   'center': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.0, 0.03640465728506993]},
   'radius': 0.004529657285069946,
   'is_trigo': True,
   'angle1': -2.707611789919429,
   'angle2': -0.4339808636703642,
   'angle': 2.273630926249065},
  {'name': '',
   'point1': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.0041097537330863606, 0.0345]},
   'point2': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.004999999999999999, 0.0345]},
   'object_class': 'volmdlr.core.LineSegment2D'},
  {'name': '',
   'package_version': None,
   'object_class': 'volmdlr.core.Arc2D',
   'interior': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.005353553390593273, 0.034353553390593275]},
   'start': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.004999999999999999, 0.0345]},
   'end': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.0055, 0.034]},
   'center': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.005000000000000158, 0.034000000000000294]},
   'radius': 0.000499999999999709,
   'is_trigo': False,
   'angle1': -5.828670879283917e-13,
   'angle2': 1.570796326795214,
   'angle': 1.570796326795797},
  {'name': '',
   'point1': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.0055, 0.034]},
   'point2': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.0055, 0.028]},
   'object_class': 'volmdlr.core.LineSegment2D'},
  {'name': '',
   'package_version': None,
   'object_class': 'volmdlr.core.Arc2D',
   'interior': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.005353553390593273, 0.027646446609406728]},
   'start': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.0055, 0.028]},
   'end': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.004999999999999999, 0.0275]},
   'center': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.005000000000000029, 0.028000000000000094]},
   'radius': 0.000499999999999971,
   'is_trigo': False,
   'angle1': -1.5707963267949556,
   'angle2': -1.8735013540550605e-13,
   'angle': 1.5707963267947682},
  {'name': '',
   'point1': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [0.004999999999999999, 0.0275]},
   'point2': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.004999999999999999, 0.0275]},
   'object_class': 'volmdlr.core.LineSegment2D'},
  {'name': '',
   'package_version': None,
   'object_class': 'volmdlr.core.Arc2D',
   'interior': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.005353553390593273, 0.027646446609406728]},
   'start': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.004999999999999999, 0.0275]},
   'end': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.0055, 0.028]},
   'center': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.005000000000000029, 0.028000000000000094]},
   'radius': 0.0005000000000000941,
   'is_trigo': False,
   'angle1': -3.1415926535896057,
   'angle2': -1.5707963267948377,
   'angle': 1.570796326794768},
  {'name': '',
   'point1': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.0055, 0.028]},
   'point2': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.0055, 0.034]},
   'object_class': 'volmdlr.core.LineSegment2D'},
  {'name': '',
   'package_version': None,
   'object_class': 'volmdlr.core.Arc2D',
   'interior': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.005353553390593273, 0.034353553390593275]},
   'start': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.0055, 0.034]},
   'end': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.004999999999999999, 0.0345]},
   'center': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.005000000000000158, 0.034000000000000294]},
   'radius': 0.0004999999999998417,
   'is_trigo': False,
   'angle1': 1.5707963267945793,
   'angle2': -3.1415926535892105,
   'angle': 1.5707963267957965},
  {'name': '',
   'point1': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.004999999999999999, 0.0345]},
   'point2': {'name': '',
    'package_version': None,
    'object_class': 'volmdlr.core_compiled.Point2D',
    'vector': [-0.0041097537330863606, 0.0345]},
   'object_class': 'volmdlr.core.LineSegment2D'}],
 '_utd_analysis': False}

contour = vm.Contour2D.dict_to_object(contour_dict)
profile5 = RevolvedProfile(0.15*vm.Y3D, vm.X3D, vm.Z3D, contour, 0.15*vm.Y3D, vm.X3D)

model = vm.VolumeModel([profile1, profile2, profile3, profile4, profile5])
model.babylonjs()


