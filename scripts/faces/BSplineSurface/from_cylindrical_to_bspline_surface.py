#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: s.bendjebla
"""

# %% Librairies

import volmdlr as vm
import volmdlr.faces as vmf
import numpy as npy

# %% Cylindrical-surface/face 

cylindrical_surface = vmf.CylindricalSurface3D.load_from_file('cylindrical_surface_1.json')

cylindrical_face = cylindrical_surface.rectangular_cut(10, 11, 0, 5)

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