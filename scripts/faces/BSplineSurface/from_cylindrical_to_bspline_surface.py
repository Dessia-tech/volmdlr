#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import volmdlr.faces as vmf
import volmdlr.step as vms
import numpy as npy

# %% Read Step file

file_path = 'cylindrical_surface_1.step'

# Chargement des fichiers step
step_file = vms.Step(file_path)

# Extraction des primitives et faces
model = step_file.to_volume_model()
primitives = model.primitives

faces = []
for primitive in primitives:
    faces.extend(primitive.faces)

# %% Cylindrical face 

cylindrical_face = faces[0]

# %% Bspline-surface/face 

degree_u, degree_v = 3, 3
bspline_surface = vmf.BSplineSurface3D.from_cylindrical_face(cylindrical_face, degree_u, degree_v)

bspline_face = bspline_surface.rectangular_cut(0, 1, 0, 1)


# %% Display

cylindrical_face.color = [list(npy.random.choice(range(255), size=1))[0] / 256, 
                          list(npy.random.choice(range(255), size=1))[0] / 256, 
                          list(npy.random.choice(range(255), size=1))[0] / 256]

bspline_face.color = [list(npy.random.choice(range(255), size=1))[0] / 256, 
                      list(npy.random.choice(range(255), size=1))[0] / 256, 
                      list(npy.random.choice(range(255), size=1))[0] / 256]

vm.core.VolumeModel([cylindrical_face, bspline_face]).babylonjs()