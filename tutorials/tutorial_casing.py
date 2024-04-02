#!/usr/bin/env python
# coding: utf-8

# [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/Dessia-tech/tutorials/blob/notebook/notebook/power_transmission.ipynb)

# # Initialization (only for Colab users)



# Reload the web page (quit) and execute this cell
import sys
print("User Current Version:-", sys.version)
if sys.version[:3] != "3.9":
  raise SystemError("Try to reload the web page (quit)")


# # Python imports



import math
import numpy as npy
from scipy.optimize import minimize
from typing import List, Tuple
from dessia_common.core import PhysicalObject, DessiaObject
import volmdlr.primitives3d as primitives3d
import volmdlr.primitives2d as primitives2d


# # Casing

p1 = volmdlr.Point2D(0, 0)
p2 = volmdlr.Point2D(0.3, 0)
p3 = volmdlr.Point2D(0.33, 0.22)
p4 = volmdlr.Point2D(0.2, 0.08)
p5 = volmdlr.Point2D(0.16, 0.18)
p6 = volmdlr.Point2D(0.05, 0.20)

inner_contour = primitives2d.ClosedRoundedLineSegments2D(
    points=[p1, p2, p3, p4, p5, p6],
    radius={0: 0.01, 1: 0.01, 2: 0.015, 3: 0.020, 4: 0.012, 5:0.008},
    adapt_radius=True)

thickness = 0.005
outer_contour = inner_contour.offset(-thickness)



height = 0.080

sides = primitives3d.ExtrudedProfile(
    frame=vm.OXYZ,
    outer_contour2d=outer_contour,
    inner_contours2d=[inner_contour],
    extrusion_length=(height - 2 * thickness),
    name='sides')

# Jupyter Notebook usage
sides.save_babylonjs_to_file(filename='/tmp/sides')


bottom = primitives3d.ExtrudedProfile(
    frame=vm.OXYZ,
    outer_contour2d=outer_contour,
    inner_contours2d=[],
    extrusion_length=-thickness,
    name='bottom')

# Jupyter Notebook usage
bottom.save_babylonjs_to_file(filename='/tmp/bottom')


belt_width = 0.011
belt_outer_contour = inner_contour.offset(-belt_width)

belt = primitives3d.ExtrudedProfile(
    frame=vm.Frame3D((vm.Z3D*(height - 2*thickness)).to_point(), vm.X3D, vm.Y3D, vm.Z3D),
    outer_contour2d=belt_outer_contour,
    inner_contours2d=[inner_contour],
    extrusion_length=thickness,
    name='belt')

# Jupyter Notebook usage
belt.save_babylonjs_to_file(filename='/tmp/belt')


casing = volmdlr.model.VolumeModel([bottom, sides, belt], name='Casing')

# Jupyter Notebook usage
casing.save_babylonjs_to_file(filename='/tmp/casing')


lid = primitives3d.ExtrudedProfile(
    frame=vm.Frame3D((vm.Z3D*(height - thickness)).to_point(), vm.X3D, vm.Y3D, vm.Z3D),
    outer_contour2d=belt_outer_contour,
    inner_contours2d=[],
    extrusion_length=thickness,
    name='lid')

casing = volmdlr.model.VolumeModel([bottom, sides, belt, lid], name='Casing')

# Jupyter Notebook usage
casing.save_babylonjs_to_file(filename='/tmp/casing')


casing.to_step('/tmp/casing')
