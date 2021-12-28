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

splitted_faces = []
for i,s in enumerate(splitted_surfaces):
    splitted_faces.append(s.rectangular_cut(0,1,0,1))
    splitted_faces[i].color = [list(npy.random.choice(range(255), size=1))[0] / 256, 
                               list(npy.random.choice(range(255), size=1))[0] / 256, 
                               list(npy.random.choice(range(255), size=1))[0] / 256]
    
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
    splitted_faces[i].color = [list(npy.random.choice(range(255), size=1))[0] / 256, 
                               list(npy.random.choice(range(255), size=1))[0] / 256, 
                               list(npy.random.choice(range(255), size=1))[0] / 256]
    
# %%% Display

ax = bspline_surface.rectangular_cut(0, 1, 0, 1).plot()
for f in splitted_faces:
    f.plot(ax=ax, color=f.color)

# vm.core.VolumeModel(splitted_faces).babylonjs()


# %% (3) BSpline-surface split - bspline_curve_split 
# split the bspline_surface, into 2 surfaces, using a bspline curve

bspline_curve3d = vm.edges.BSplineCurve3D.load_from_file('bspline_curve3d_interpolated.json')
splitted_surfaces = bspline_surface.split_surface_with_bspline_curve(bspline_curve3d)

splitted_faces = []
for i,s in enumerate(splitted_surfaces):
    splitted_faces.append(s.rectangular_cut(0,1,0,1))
    splitted_faces[i].color = [list(npy.random.choice(range(255), size=1))[0] / 256, 
                               list(npy.random.choice(range(255), size=1))[0] / 256, 
                               list(npy.random.choice(range(255), size=1))[0] / 256]

# %%% Display

ax = bspline_surface.rectangular_cut(0, 1, 0, 1).plot()
for f in splitted_faces:
    f.plot(ax=ax, color=f.color)

# vm.core.VolumeModel(splitted_faces).babylonjs()
