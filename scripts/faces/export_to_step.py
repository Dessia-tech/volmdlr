#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 14 2022

@author: s.bendjebla
"""

import volmdlr.core
import volmdlr.step


#Export
filepath = '../faces/volume_model.json'
model = volmdlr.core.VolumeModel.load_from_file(filepath)

#Check
filepath = '/home/bendjebla/Téléchargements/step.step'
step = volmdlr.step.Step.from_file(filepath=filepath)
model = step.to_volume_model()    
model.babylonjs()
