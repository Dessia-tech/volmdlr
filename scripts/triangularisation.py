#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 2 2022

@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import volmdlr.wires
import volmdlr.faces

# import matplotlib.pyplot as plt

# %% Contours (inners & outer)


p = [vm.Point2D(-0.3, -0.2), vm.Point2D(0.3, -0.2),
      vm.Point2D(0.2, 0.2), vm.Point2D(0, 0.3), vm.Point2D(-0.2, 0.2)]
contour = vm.wires.Contour2D(vm.wires.ClosedPolygon2D(p).primitives)


pts = vm.wires.Circle2D(vm.O2D, radius=0.05).tessellation_points(resolution=40)
inner_contour = vm.wires.ClosedPolygon2D(pts)
pts2 = vm.wires.Circle2D(vm.Point2D(-0.13, 0.15), radius=0.02).tessellation_points(resolution=40)
inner_contour = [vm.wires.ClosedPolygon2D(pts), vm.wires.ClosedPolygon2D(pts2)]

# %% 1/ Surface2D + inners

surface = volmdlr.faces.Surface2D(outer_contour=contour, inner_contours=inner_contour)

# %%% Triangularisation

x_density, y_density = 5, 5

triangularisation = surface.triangularisation_2(x_density, y_density)

ax =surface.plot()
for element in triangularisation:
    element.plot(ax, 'b')

# %% 2/ Surface2D + not inners

surface = volmdlr.faces.Surface2D(outer_contour=contour, inner_contours=[])

# %%% Triangularisation

triangularisation = surface.triangularisation_2(x_density, y_density)

ax =surface.plot()
for element in triangularisation:
    element.plot(ax, 'b')
