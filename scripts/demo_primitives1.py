#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 15:33:01 2017

@author: steven
"""
import numpy as npy
import volmdlr as vm
import volmdlr.primitives3D as primitives3D

cylinder1 = primitives3D.Cylinder(vm.Vector3D((0.,0.,0.)), vm.Vector3D((1.,0.,0.)),
                                0.03, 0.02, 'cylinder1')
cylinder2 = primitives3D.HollowCylinder(vm.Vector3D((0,0.1,0.)), vm.Vector3D((1.,0.,0.)),
                                      0.02, 0.06, 0.03,'cylinder2')
#profile=primitives3D.ExtrudedProfile((0,0,0),(1,0,0),(0,1,0),[(0,0),(0.1,0.),(0.15,0.4),(0.,0.3)],{0:0.05,2:0.01},(0,0,0.2))

model=vm.VolumeModel([('Cylinder 1', [cylinder1]), ('Cylinder 2', [cylinder2])])

#profile.MPLPlot((0,0,0),(1,0,0),(0,1,0))

model.FreeCADExport('cylinders')

#print(model.BabylonScript())
#model.BabylonShow()