#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: s.bendjebla
"""

# %% Libraries

# import volmdlr.step as vms
import matplotlib.pyplot as plt

import volmdlr.grid
from volmdlr.models import bspline_surfaces
from volmdlr import faces
# %% Read Step file

# file_path = 'bspline_surface_2.step'

# step_file = vms.Step.from_file(file_path)

# model = step_file.to_volume_model()
# primitives = model.primitives

# faces = []
# for primitive in primitives:
#     faces.extend(primitive.faces)

# %% Bspline face/surface/contour

bspline_face = faces.BSplineFace3D.from_surface_rectangular_cut(bspline_surfaces.bspline_surface_2, 0, 1, 0, 1)

bspline_surface = bspline_surfaces.bspline_surface_2

contour3d = bspline_face.outer_contour3d

# %% To 2d
# %%% Parametric frame

contour2d = bspline_surface.contour3d_to_2d(contour3d)


# %%% Dimensionned frame

grid2d = volmdlr.grid.Grid2D.from_properties((0,1),(0,1),(10,10))
contour2d_dim = bspline_surface.contour2d_parametric_to_dimension(contour2d, grid2d)

# %%% Display

fig, (ax1, ax2) = plt.subplots(1, 2)
contour2d.plot(ax=ax1)
contour2d_dim.plot(ax=ax2)

ax1.set_title('Parametric frame')
ax2.set_title('Dimensionned frame')
