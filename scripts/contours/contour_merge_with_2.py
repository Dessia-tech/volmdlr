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

p1 = [vm.Point2D(6, 0), vm.Point2D(6, 2),
      vm.Point2D(3, 2), vm.Point2D(3, 1),
      vm.Point2D(1, 1), vm.Point2D(1, 2), 
      vm.Point2D(0, 2), vm.Point2D(0, 0)]

p2 = [vm.Point2D(-1, 2), vm.Point2D(1, 2),
      vm.Point2D(1, 3), vm.Point2D(3, 3),
      vm.Point2D(3, 2), vm.Point2D(4, 2),
      vm.Point2D(4, 4), vm.Point2D(-1, 4)]

contour1 = volmdlr.wires.Contour2D((vm.wires.ClosedPolygon2D(p1).discretized_primitives(5))).order_contour()
contour2 = volmdlr.wires.Contour2D((vm.wires.ClosedPolygon2D(p2).discretized_primitives(3))).order_contour()


fig, axs = plt.subplots(2, 4)
for i in range(0, 2):
    for j in range(0, 4):
        contour1.plot(ax=axs[i][j], color='k')
        contour2.plot(ax=axs[i][j], color='k')
        if (i == 0 and j == 0):
            for p in contour1.primitives:
                p.plot(ax=axs[i][j], color='b', width=4)
            axs[i][j].set_title('contour1-shared_primitives_extremities')

        if (i == 1 and j == 0):
            for p in contour2.primitives:
                p.plot(ax=axs[i][j], color='m', width=4)
            axs[i][j].set_title('contour2-shared_primitives_extremities')


# %% is_sharing_primitives_with

is_sharing = contour1.is_sharing_primitives_with(contour2)
if is_sharing:
    print('Contour1 is sharing primitives with Contour2')
else:
    print('The contours are not adjacent. They dont share primitives')


# %% shared_primitives_extremities

shared_primitives_extremities = contour1.shared_primitives_extremities(contour2)

j=0
for i in range(0, 2):
    for p in shared_primitives_extremities:
        p.plot(ax=axs[i][j])
    axs[i][j].set_title('shared_primitives_extremities')


# %% shared_primitives_with

shared_primitives = contour1.shared_primitives_with(contour2)

title = ['contour1-shared_primitives_extremities', 'contour2-shared_primitives_extremities']

j=j+1
for i in range(0, 2):
    for prim in shared_primitives:
        for p in prim:
            p.plot(ax=axs[i][j], color='r', width=4)
    axs[i][j].set_title(title[i])


# %% merge_primitives_with

merged_primitives = contour1.merge_primitives_with(contour2)

title = ['contour1-primitives_to_be_merged', 'contour2-primitives_to_be_merged']

j=j+1
for i in range(0, 2):
    for p in merged_primitives:
        if contour1.point_over_contour(p.start) and contour1.point_over_contour(p.end):
            p.plot(ax=axs[0][j], color='g', width=4)
        else:
            p.plot(ax=axs[1][j], color='g', width=4)

    axs[i][j].set_title(title[i])


# %% merge_with

merged_contour = contour1.merge_with(contour2)

title = ['merged_contour_outer', 'merged_contour_inner']
colors = ['b', 'm']

j=j+1
for i in range(0, 2):
    for p in merged_contour[i].primitives:
        p.plot(ax=axs[i][j], color=colors[i], width=4)
    axs[i][j].set_title(title[i])
