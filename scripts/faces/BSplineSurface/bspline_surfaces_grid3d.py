#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import volmdlr.step as vms

# %% Read Step file

files_path = ['../scripts/faces/BSplineSurface/bspline_surface_1.step', '../scripts/faces/BSplineSurface/bspline_surface_2.step']
bspline_faces = []

for file_path in files_path: 
    step_file = vms.Step.from_file(file_path)
    
    model = step_file.to_volume_model()
    primitives = model.primitives
    
    faces = []
    for primitive in primitives:
        faces.extend(primitive.faces)
    
    bspline_faces.append(faces[0])


# %% Bspline surfaces

bspline_surfaces = [bspline_faces[0].surface3d, bspline_faces[1].surface3d]


# %% Grdi3d initial

points_x, points_y, xmin, xmax, ymin, ymax = 5, 5, 0, 1, 0, 1

points3d = []
for i, bspline in enumerate(bspline_surfaces):
    grid2d = vm.Point2D.grid2d(points_x, points_y, xmin, xmax, ymin, ymax)
    grid3d = []
    for p in grid2d:
        grid3d.append(bspline.point2d_to_3d(p))

    points3d.append(grid3d)


# %%% Display

ax1 = points3d[0][0].plot()
for i, points in enumerate(points3d):
    for k,p in enumerate(points):
        if k<points_x:
            if k == 0:
                p.plot(ax=ax1, color='g')  
            else:
                p.plot(ax=ax1, color='r')
        else:
            p.plot(ax=ax1)


# %% Grdi3d with directions

corresponding_directions, grid2d_direction = bspline_faces[0].pair_with(bspline_faces[1])
points_x, points_y, xmin, xmax, ymin, ymax = 5, 5, 0, 1, 0, 1

points3d = []
for i, bspline in enumerate(bspline_surfaces):
    grid2d = vm.Point2D.grid2d_with_direction(points_x, points_y, xmin, xmax, ymin, ymax, grid2d_direction[i])[0]
    grid3d = []
    for p in grid2d:
        grid3d.append(bspline.point2d_to_3d(p))

    points3d.append(grid3d)


# %%% Display

ax1 = points3d[0][0].plot()
for i, points in enumerate(points3d):
    for k, p in enumerate(points):
        if k < points_x:
            if k == 0:
                p.plot(ax=ax1, color='g')
            else:
                p.plot(ax=ax1, color='r')
        else:
            p.plot(ax=ax1)
