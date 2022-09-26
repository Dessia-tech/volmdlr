#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import io
import os
# import volmdlr as vm
import volmdlr.step

# import volmdlr.cloud as vmcd

for step_file in [
    'tore1.step',
    # 'iter8finaldesign.stp',
    'cone1.step',
    'cone2.step',
    'cylinder.step',
    'block.step',
    'bracket2.step',
    # '2_bspline_faces.stp'# Uncomment when bug of delta fixed!
]:
    print('Reading step file: ', step_file)
    # filepath = os.path.join('step', step_file)
    step = volmdlr.step.Step.from_file(filepath=step_file)
    model = step.to_volume_model()
    assert len(model.primitives) > 0.
    model.to_step(step_file + '_reexport')
    model.babylonjs()

    file_io = io.FileIO(step_file, 'r')
    step = volmdlr.step.Step.from_stream(stream=file_io)
    model = step.to_volume_model()
    assert len(model.primitives) > 0.
    model.to_step(step_file + '_reexport')

    model2 = model.copy()

    # model2 = model.copy()
    # assert model == model2

    model._check_platform()