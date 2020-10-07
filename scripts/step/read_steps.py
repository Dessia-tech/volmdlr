#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm
import volmdlr.step

step = volmdlr.step.Step('tore1.step')
model = vm.VolumeModel.from_step('tore1.step')

model.babylonjs()