#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import numpy as npy

#%%  BSpline-surface definition 

script = 'bspline_surface_definition.py'
print('\n## Executing script {}'.format(script))
exec(open(script).read())

#%% (1) BSpline-surface split - u_split 
# split the bspline_surface at the input parametric coordinate (u) on the u-direction, into 2 surfaces

u = 0.2
splitted_surfaces = bspline_surface.split_surface_u(u)

random_colors = []
splitted_faces = []
for i,s in enumerate(splitted_surfaces):
    splitted_faces.append(s.rectangular_cut(0,1,0,1))
    random_colors.append([list(npy.random.choice(range(255), size=1))[0] / 256, 
                          list(npy.random.choice(range(255), size=1))[0] / 256, 
                          list(npy.random.choice(range(255), size=1))[0] / 256])
    splitted_faces[i].color = random_colors[i]
    
# %%% Display

ax = bspline_surface.rectangular_cut(0, 1, 0, 1).plot()
for f in splitted_faces:
    f.plot(ax=ax, color=f.color)

# vm.core.VolumeModel(splitted_faces).babylonjs()

# %% (2) BSpline-surface split - v_split 
# split the bspline_surface at the input parametric coordinate (v) on the v-direction, into 2 surfaces

v = 0.4
splitted_surfaces = bspline_surface.split_surface_v(v)

splitted_faces = []
for i,s in enumerate(splitted_surfaces):
    splitted_faces.append(s.rectangular_cut(0,1,0,1))
    splitted_faces[i].color = random_colors[i]
    
# %%% Display

ax = bspline_surface.rectangular_cut(0, 1, 0, 1).plot()
for f in splitted_faces:
    f.plot(ax=ax, color=f.color)

# vm.core.VolumeModel(splitted_faces).babylonjs()

# %% (3) BSpline-surface split - bspline_curve_split 
# split the bspline_surface, into 2 surfaces, using a bspline curve

# %%% Bspline-curve definition

points2d = [vm.Point2D(0, 0.1),
            vm.Point2D(0.2, 0.3),
            vm.Point2D(0.4, 0.4),
            vm.Point2D(0.5, 0.6),
            vm.Point2D(0.6, 0.7),
            vm.Point2D(0.8, 0.8),
            vm.Point2D(1,0.9)]

bspline_curve3d = bspline_surface.bsplinecurve2d_to_3d(vm.edges.BSplineCurve2D.from_points_interpolation(points2d, 2))[0]

# %%% Split surface

splitted_surfaces = bspline_surface.split_surface_with_bspline_curve(bspline_curve3d)

splitted_faces = []
for i,s in enumerate(splitted_surfaces):
    splitted_faces.append(s.rectangular_cut(0,1,0,1))
    splitted_faces[i].color = random_colors[i]
    
# %%% Display

ax = bspline_surface.rectangular_cut(0, 1, 0, 1).plot()
for f in splitted_faces:
    f.plot(ax=ax, color=f.color)

# vm.core.VolumeModel(splitted_faces).babylonjs()
