#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm
import volmdlr.step

step1 = volmdlr.step.Step('tore1.step')
step2 = volmdlr.step.Step('cylindre.step')
step3 = volmdlr.step.Step('block.step')
model = step3.to_volume_model()

model.babylonjs()