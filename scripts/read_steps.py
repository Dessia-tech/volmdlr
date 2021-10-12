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

    # 'STEP_test1.stp',
    # File "/home/axel/Documents/git/volmdlr/volmdlr/step.py",
    # line 406, in instanciate
    #     volmdlr_object = object_dict[arguments[1]]
    # KeyError: 1176

    #  'iso4162M16x55.step',
    # File "/home/axel/Documents/git/volmdlr/volmdlr/wires.py",
    # line 1153, in cut_by_line
    # NotImplementedError: 3 intersections not supported yet

    # 'aircraft_engine.step'
    # File "/home/axel/Documents/git/volmdlr/volmdlr/wires.py",
    # line 2662, in from_step
    # NotImplementedError: ('Edges of contour not follwing each
    # other', 'delta = 0.0001325146419990375, 0.00021337880199912007,
    # 0.06408522304249997, 1.520937009991119e-05,
    # 3.126516413497999, 4.779537000842993e-06')
  ]:
    print('filename: ', step_file)
    step = volmdlr.step.Step('step/'+step_file)
    model = step.to_volume_model()
    model.to_step(step_file+'_reexport')

    model.babylonjs()
    
model2 = model.copy()

assert model == model2
