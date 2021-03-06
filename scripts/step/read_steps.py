#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm
import volmdlr.step

for step_file in [#'tore1.step',
                   'cone1.step',
                  #'cone2.step',
                  #'iso4162M16x55.step',
                   'cylindre.step',
                  'block.step',
                  ]:
    step = volmdlr.step.Step(step_file)
    model = step.to_volume_model()
    print(model.primitives)
    
    model.babylonjs()
    
model2 = model.copy()

assert model == model2