#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 2 2022

@author: s.bendjebla
"""

# %% Librairies

import volmdlr
import volmdlr.grid
import matplotlib.pyplot as plt

# %% Grid2D

x_limits = (0,1)
y_limits = (0,1)
points_nbr = (5,5)
direction = ['+x','+y']

grid2d = volmdlr.grid.Grid2D.from_properties(x_limits, y_limits, points_nbr, direction)

patterns = grid2d.grid_pattern()

# %% Displays

print('x limits = ', grid2d.limits_xy[0])
print('y limits = ', grid2d.limits_xy[1])
print('x points = ', grid2d.points_xy[0])
print('y points = ', grid2d.points_xy[1])


fig, axs = plt.subplots(1, 2)
titles = ['Grid2d', 'Patterns']

for i, ax in enumerate(axs):
    for points in grid2d.list_points:
        for point in points:
            point.plot(ax=ax)
    ax.set_title(titles[i])
    if i == 1:
        for pattern in patterns:
            pattern.plot(ax=ax, color='b')    
