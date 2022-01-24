#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import io

import volmdlr as vm
import volmdlr.step
import volmdlr.cloud as vmcd

for step_file in [
    'tore1.step',
    'cone1.step',
    'cone2.step',
    'cylinder.step',
    'block.step',
    # '2_bspline_faces.stp'# Uncomment when bug of delta fixed!
  ]:
    print('filename: ', step_file)
    step = volmdlr.step.Step('step/'+step_file)
    model = step.to_volume_model()
    assert len(model.primitives) > 0.
    model.to_step(step_file+'_reexport')
    model.babylonjs()

    file_io = io.FileIO('step/' + step_file, 'r')
    step = volmdlr.step.Step(stream=file_io)
    model = step.to_volume_model()
    assert len(model.primitives) > 0.
    model.to_step(step_file + '_reexport')

    model2 = model.copy()
    
    # model2 = model.copy()
    # assert model == model2

    model._check_platform()