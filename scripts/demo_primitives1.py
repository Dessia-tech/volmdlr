#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 15:33:01 2017

@author: steven
"""

import volmdlr as vm
import volmdlr.primitives as primitives

cylinder1=primitives.Cylinder((0,0,0),(1,0,0),0.05,0.02)
cylinder2=primitives.HollowCylinder((0,0.1,0),(1,0,0),0.02,0.06,0.03)
profile=primitives.ExtrudedProfile([(0,0,0),(0.1,0.02,0),(0.2,0.3,0),(0.1,0.4,0),(-0.1,0.2,0),(-0.2,0.1,0)],[],[],(0,0,0.2))

model=vm.VolumeModel([cylinder1,cylinder2,profile])
model.FreeCADExport('python','cylinders','/usr/lib/freecad-daily/lib/')