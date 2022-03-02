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

# %% Surface2D + inners

surface = volmdlr.faces.Surface2D(outer_contour=contour, inner_contours=[inner_contour])

# %% Triangularisation

x_density, y_density = 5, 5

triangularisation = surface.triangularisation_2(x_density, y_density)

ax =surface.plot()
for element in triangularisation:
    element.plot(ax, 'b')


primitives = [vm.edges.LineSegment2D(vm.Point2D(-0.35, -0.1), vm.Point2D(-0.1, 0)),
              vm.edges.LineSegment2D(vm.Point2D(-0.1, 0), vm.Point2D(0.2, 0.2)),
              vm.edges.LineSegment2D(vm.Point2D(0.2, 0.2), vm.Point2D(0.3, 0.3))]

wire = vm.wires.Wire2D(primitives)


# %% Surface2D + not inners

surface = volmdlr.faces.Surface2D(outer_contour=contour, inner_contours=[])

# %% Triangularisation

triangularisation = surface.triangularisation_2(x_density, y_density)

ax =surface.plot()
for element in triangularisation:
    element.plot(ax, 'b')
