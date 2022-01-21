#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import volmdlr as vm
import volmdlr.step
import volmdlr.cloud as vmcd

for step_file in [
    'tore1.step',
    'cone1.step',
    'cone2.step',
    'cylinder.step',
    'block.step',
    # '4_bspline_faces.step',
    '2_bspline_faces.stp'
    # 'STEP_test1.stp',
    # 'iso4162M16x55.step',
    # 'aircraft_engine.step'
  ]:
    print('filename: ', step_file)
    step = volmdlr.step.Step('step/'+step_file)
    model = step.to_volume_model()
    assert len(model.primitives) > 0.
    model.to_step(step_file+'_reexport')

    model.babylonjs()
    
model2 = model.copy()

assert model == model2
