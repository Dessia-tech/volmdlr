#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 
@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import volmdlr.edges
import volmdlr.wires
import matplotlib.pyplot as plt

# %% Initial Data

points = [vm.Point2D(6, 0), vm.Point2D(6, 2),
          vm.Point2D(3, 2), vm.Point2D(3, 1),
          vm.Point2D(1, 1), vm.Point2D(1, 2), 
          vm.Point2D(0, 2), vm.Point2D(0, 0)]

contour = vm.wires.Contour2D.from_points(points)
line = vm.edges.Line2D(vm.Point2D(-1, 3), vm.Point2D(7, 3))

# %% Symmetry

axial_contour = contour.axial_symmetry(line)

fig, ax = plt.subplots()
ax.set_aspect('equal')

line.plot(ax)

contour.plot(ax=ax, color='r')
axial_contour.plot(ax=ax, color='g')
