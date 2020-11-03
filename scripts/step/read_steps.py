#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm
import volmdlr.step

for step_file in [#'tore1.step',
                  'cylindre.step',
                  'block.step',
                  ]:
    step = volmdlr.step.Step(step_file)
    model = step.to_volume_model()
    print(model.primitives)
    
    model.babylonjs()
    del model
