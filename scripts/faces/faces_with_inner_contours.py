#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 2022

@author: s.bendjebla
"""

import os
import volmdlr.step

for step_file in [
    'bsplineface_with_inner_contours.stp',
    'planeface_with_inner_contour.step']:
    
    print('Reading step file: ', step_file)
    
    filepath = os.path.join('../faces', step_file)
    step = volmdlr.step.Step.from_file(filepath=filepath)
    model = step.to_volume_model()    
    primitives = model.primitives

    faces = []
    for primitive in primitives:
        faces.extend(primitive.faces)

    surfaces_3d, outer_contours_3d, inner_contours_3d = [], [], []
    for face in faces:
        surfaces_3d.append(face.surface3d)
        outer_contours_3d.append(face.outer_contour3d)
        inner_contours_3d.extend(face.inner_contours3d)
    
    assert len(inner_contours_3d) > 0
