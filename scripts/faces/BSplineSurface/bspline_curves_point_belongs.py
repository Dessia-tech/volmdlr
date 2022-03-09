#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import volmdlr.edges
import matplotlib.patches as mpatches

#%%  BSpline-curve 2D definition 

# Parameters
knots = [0.0, 0.29927104501765267, 0.46883825500950527, 0.6072437200108631, 0.7456491850122211, 1.0]
knot_multiplicities = [3, 1, 1, 1, 1, 3]
degree = 2
control_points = [vm.Point2D(0.0, 0.10000000000000002),
                  vm.Point2D(0.10497581068453007, 0.28289448655392097),
                  vm.Point2D(0.4204851473531051, 0.37983356450729255),
                  vm.Point2D(0.47891562865564064, 0.5900187216711663),
                  vm.Point2D(0.6122988905036818, 0.718784316931106),
                  vm.Point2D(0.852910655633803, 0.8204050240455999),
                  vm.Point2D(1.0, 0.9)]

# Bspline_curve
bspline_curve2d = vm.edges.BSplineCurve2D(degree,
                                          control_points,
                                          knot_multiplicities,
                                          knots)


# %% Point-belongs

points2d = [vm.Point2D(0.18337173875665347, 0.29040613195152876),
            vm.Point2D(0.4, 0.5),
            vm.Point2D(0.5, 0.9),
            vm.Point2D(0.7632838172192403, 0.7833495785460081),
            vm.Point2D(1, 0.7),
            vm.Point2D(0.3393332354977111, 0.2220637856052557)]

colors = {True:'g', False:'r'}

ax = bspline_curve2d.plot()
for i,p in enumerate(points2d):
    belongs = bspline_curve2d.point_belongs(p)
    p.plot(ax=ax, color=colors[belongs])

ax.legend(handles=[mpatches.Patch(color='green', label='Belongs'),
                   mpatches.Patch(color='red', label='Not')])

