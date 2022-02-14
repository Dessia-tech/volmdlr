#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: s.bendjebla
"""

# %% Librairies

import volmdlr.step as vms
import matplotlib.pyplot as plt

# %% Read Step file

file_path = 'bspline_surface_2.step'

step_file = vms.Step.from_file(file_path)

model = step_file.to_volume_model()
primitives = model.primitives

faces = []
for primitive in primitives:
    faces.extend(primitive.faces)

# %% Bspline face/surface/contour

bspline_face = faces[0]

bspline_surface = bspline_face.surface3d

contour3d = bspline_face.outer_contour3d

# %% To 2d
# %%% Parametric frame

contour2d = bspline_surface.contour3d_to_2d(contour3d)


# %%% Dimensionned frame

contour2d_dim = bspline_surface.contour3d_to_2d_with_dimension(contour3d, 10, 10)

# %%% Display

fig, (ax1, ax2) = plt.subplots(1, 2)
contour2d.plot(ax=ax1)
contour2d_dim.plot(ax=ax2)

ax1.set_title('Parametric frame')
ax2.set_title('Dimensionned frame')

