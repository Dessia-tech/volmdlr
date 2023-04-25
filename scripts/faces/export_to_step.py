#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 2022

@author: s.bendjebla
"""

# %% Libraries

import volmdlr.core
import volmdlr.step
from volmdlr import faces, shells
from volmdlr.models import bspline_surfaces

# %% BsplineFaces3D

bspline_faces = [faces.BSplineFace3D.from_surface_rectangular_cut(bspline_surfaces.bspline_surface_2, 0,1,0,1)]

# %% Export

model = volmdlr.core.VolumeModel([shells.OpenShell3D(bspline_faces)])

model.to_step('model_to_step.stp')

# %% Check

filepath = 'model_to_step.stp'

step = volmdlr.step.Step.from_file(filepath)
model_imported = step.to_volume_model()

# assert model_imported == model
