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

titles = ["Initial Contour2d + Wire2d", "1st Cutted Contour2d 'green'", "2nd Cutted Contour2d 'blue'"]
colors = ['g', 'b']
for i in range(len(axs)):
    contour.plot(ax=axs[i])
    for prim in wire.primitives:
        prim.plot(ax=axs[i], width=2, color='r')
    axs[i].set_title(titles[i])
    if i !=0:
        for prim in contours[i-1].primitives:
            prim.plot(ax=axs[i], width=2, color=colors[i-1])
