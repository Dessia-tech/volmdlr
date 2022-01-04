#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
 
@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import volmdlr.wires
import matplotlib.pyplot as plt

# %% Contours2d

p1 = [vm.Point2D(4, 0), vm.Point2D(4, 2), 
      vm.Point2D(3, 2), vm.Point2D(3, 1),
      vm.Point2D(1, 1), vm.Point2D(1, 2), 
      vm.Point2D(0, 2), vm.Point2D(0, 0)]

p2 = [vm.Point2D(0, 2), vm.Point2D(1, 2),
      vm.Point2D(1, 3), vm.Point2D(3, 3),
      vm.Point2D(3, 2), vm.Point2D(4, 2),
      vm.Point2D(4, 4), vm.Point2D(0, 4)]

contour1 = vm.wires.ClosedPolygon2D(p1)
contour2 = vm.wires.ClosedPolygon2D(p2)

to_plot, title, colors = [], [], []

# %% is_sharing_primitives_with

is_sharing = contour1.is_sharing_primitives_with(contour2)
if is_sharing:
    print('Contour1 is sharing primitives with Contour2')
else:
    print('The contours are not adjacent. They dont share primitives')


# %% shared_primitives_extremities

shared_primitives_extremities = contour1.shared_primitives_extremities(contour2)

to_plot.append(shared_primitives_extremities)
title.append('shared_primitives_extremities')
colors.append('k')

# %% shared_primitives_with

shared_primitives = contour1.shared_primitives_with(contour2)

to_plot.append(shared_primitives)
title.append('shared_primitives')
colors.append('r')

# %% merge_primitives_with

merged_primitives = contour1.merge_primitives_with(contour2)

to_plot.append(merged_primitives)
title.append('primitives_to_be_merged')
colors.append('g')

# %% merge_with

merged_contour = contour1.merge_with(contour2)

to_plot.append([merged_contour])
title.append('merged_contour')
colors.append('b')

# %% Display

fig, axs = plt.subplots(1, 4)

for i, plot_ in enumerate(to_plot):
    if i == 0: 
        contour1.plot(ax=axs[i], color='b')
        contour2.plot(ax=axs[i], color='m')
    if i != 3 and i != 0: 
        contour1.plot(ax=axs[i])
        contour2.plot(ax=axs[i])
    for p in plot_:
        if isinstance(p, list):
            for element in p:
                element.plot(ax=axs[i], color=colors[i], width=3)
        else:
            if isinstance(p, volmdlr.edges.LineSegment2D):
                p.plot(ax=axs[i], color=colors[i], width=3)
            else:
                p.plot(ax=axs[i], color=colors[i])
    axs[i].set_title(title[i])
    