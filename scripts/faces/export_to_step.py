#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 2022

@author: s.bendjebla
"""

# %% Librairies

from volmdlr.models import bspline_surfaces
import volmdlr.core
import volmdlr.step

# %% BsplineFaces3D

bspline_faces = [bspline_surfaces.bspline_surface_2.rectangular_cut(0,1,0,1)]

# %% Export

model = volmdlr.core.VolumeModel(bspline_faces)

model.to_step('model_to_step.stp')

# %% Check

filepath = 'model_to_step.stp'

step = volmdlr.step.Step.from_file(filepath)
model_imported = step.to_volume_model()

assert model_imported == model
