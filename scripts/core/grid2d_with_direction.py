#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: s.bendjebla
"""

# %% Librairies

import volmdlr.grid
import matplotlib.pyplot as plt

# %% Grid2d definition 'with direction'

fig, axs = plt.subplots(2, 4)

points_x, points_y, xmin, xmax, ymin, ymax = 5, 5, 0, 1, 0, 1

grid2d_directions = [[['+x','+y'], ['-x','+y'], ['+y','+x'], ['-y','+x']],
                     [['+x','-y'], ['-x','-y'], ['+y','-x'], ['-y','-x']]]

for i in range(0, len(grid2d_directions)):
    for j in range(0, len(grid2d_directions[i])):
        grid2d = volmdlr.grid.Grid2D.from_properties((xmin, xmax), (ymin, ymax), (points_x, points_y), grid2d_directions[i][j]).points
    
        for k, p in enumerate(grid2d):
            if k<points_x:
                p.plot(ax=axs[i][j], color='r')
                x, y = p.x, p.y
                axs[i][j].annotate(k, (x, y))
            else:
                p.plot(ax=axs[i][j])

        axs[i][j].set_title(grid2d_directions[i][j])
