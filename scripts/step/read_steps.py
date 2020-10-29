#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm
import volmdlr.step

for step_file in ['tore1.step', 'cylindre.step']:
    step = volmdlr.step.Step(step_file)
    model = step.to_volume_model()
    model.babylonjs()