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

points = [vm.Point2D(0, 1),
          vm.Point2D(0.1, 0.8),
          vm.Point2D(0.2, 1.5)]

line = vm.edges.Line2D(vm.Point2D(-0.5, 1), vm.Point2D(-0.5, 8))

# %% Symmetry

axial_points = [point.axial_symmetry(line) for point in points]

fig, ax = plt.subplots()
ax.set_aspect('equal')
line.plot(ax)
[p.plot(ax=ax, color='r') for p in points]
[p.plot(ax=ax, color='g') for p in axial_points]
