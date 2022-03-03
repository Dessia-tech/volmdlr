#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import matplotlib.patches as mpatches
from volmdlr.models import bspline_surfaces

#%%  BSpline-surface definition 

# script = 'bspline_surface_definition.py'
# print('\n## Executing script {}'.format(script))
# exec(open(script).read())

bspline_surface = bspline_surfaces.bspline_surface_1

# %% BSpline-curve 2D definition

points2d = [vm.Point2D(0, 0.1),
            vm.Point2D(0.2, 0.3),
            vm.Point2D(0.4, 0.4),
            vm.Point2D(0.5, 0.6),
            vm.Point2D(0.6, 0.7),
            vm.Point2D(0.8, 0.8),
            vm.Point2D(1,0.9)]

# %%% Approximation

bspline_curve2d_approximated = vm.edges.BSplineCurve2D.from_points_approximation(points2d, 3, ctrlpts_size = 5)

# %%% Interpolation

bspline_curve2d_interpolated = vm.edges.BSplineCurve2D.from_points_interpolation(points2d, 2)

# %%% Display

ax = bspline_surfaces.bspline_face_1.surface2d.plot()
for p in points2d: 
    p.plot(ax=ax, color ='k')
    
bspline_curve2d_approximated.plot(ax=ax, color='g')
bspline_curve2d_interpolated.plot(ax=ax, color='r')

ax.legend(handles=[mpatches.Patch(color='green', label='Approximation'),
                   mpatches.Patch(color='red', label='Interpolation')])

# %% BSpline-curve 3D definition

bspline_curve3d_approximated = bspline_surface.bsplinecurve2d_to_3d(bspline_curve2d_approximated)[0]
bspline_curve3d_approximated.color = 'g'

bspline_curve3d_interpolated = bspline_surface.bsplinecurve2d_to_3d(bspline_curve2d_interpolated)[0]
bspline_curve3d_interpolated.color = 'r'

# %%% Display

ax = bspline_surfaces.bspline_face_1.plot()
bspline_curve3d_approximated.plot(ax=ax, color=bspline_curve3d_approximated.color)
bspline_curve3d_interpolated.plot(ax=ax, color=bspline_curve3d_interpolated.color)

ax.legend(handles=[mpatches.Patch(color='green', label='Approximation'),
                   mpatches.Patch(color='red', label='Interpolation')])

# primitives = [bspline_face, bspline_curve3d_approximated, bspline_curve3d_interpolated]
# vm.core.VolumeModel(primitives).babylonjs()
