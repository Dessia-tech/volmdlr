#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm
import volmdlr.step

step = volmdlr.step.Step('tore1.step')
step = volmdlr.step.Step('cylindre.step')
model = step.to_volume_model()

model.babylonjs()