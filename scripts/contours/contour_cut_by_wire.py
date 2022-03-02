#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 2 2022

@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import volmdlr.wires
import matplotlib.pyplot as plt

# %% Contour2d

p = [vm.Point2D(-0.3, -0.2), vm.Point2D(0.3, -0.2),
      vm.Point2D(0.2, 0.2), vm.Point2D(0, 0.3), vm.Point2D(-0.2, 0.2)]

contour = vm.wires.Contour2D(vm.wires.ClosedPolygon2D(p).primitives)

# %% Wire2d

primitives = [vm.edges.LineSegment2D(vm.Point2D(-0.35, -0.1), vm.Point2D(-0.1, 0)),
              vm.edges.LineSegment2D(vm.Point2D(-0.1, 0), vm.Point2D(0.2, 0.2)),
              vm.edges.LineSegment2D(vm.Point2D(0.2, 0.2), vm.Point2D(0.3, 0.3))]

wire = vm.wires.Wire2D(primitives)

# %% Cut_by_wire

contours = contour.cut_by_wire(wire)

# %% Plots

fig, axs = plt.subplots(1, 3)

contour.plot(ax=axs[0])
for prim in wire.primitives:
    prim.plot(ax=axs[0], width=2, color='r')
axs[0].set_title("Initial Contour2d + Wire2d")


contour.plot(ax=axs[1])
wire.plot(ax=axs[1])
for prim in contours[0].primitives:
    prim.plot(ax=axs[1], width=2, color='g')
axs[1].set_title("1st Cutted Contour2d 'green'")


contour.plot(ax=axs[2])
wire.plot(ax=axs[2])
for prim in contours[1].primitives:
    prim.plot(ax=axs[2], width=2, color='b')
axs[2].set_title("2nd Cutted Contour2d 'blue'")
