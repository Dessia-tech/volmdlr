#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 10:54:45 2020

@author: masfaraud
"""

import volmdlr as vm

import volmdlr.primitives3D as primitives3D


block1 = primitives3D.Block(vm.Frame3D(vm.O3D.copy(),
                                      vm.X3D.copy(), vm.Y3D.copy(), vm.Z3D.copy()),
                           color=(0.8,0.1,0.1),
                           alpha=0.5,
                           name='Block 1')
      
block2 = block1.copy()
block2.color = [0.1,0.1,0.8]
block2.name='Block 2'

                     
f1 = vm.OXYZ.copy()

f2 = f1.Rotation(vm.Z3D, 0.1, copy=True)
f2.Translation(0.1*vm.X3D, copy=False)

f3 = f2.Rotation(vm.Z3D, 0.1, copy=True)
f3.Translation(0.1*vm.X3D, copy=False)

f4 = f3.Rotation(vm.Z3D, 0.1, copy=True)
f4.Translation(0.1*vm.X3D, copy=False)

f5 = f4.Rotation(vm.Z3D, 0.1, copy=True)
f5.Translation(0.1*vm.X3D, copy=False)


model = vm.MovingVolumeModel([block1, block2], [[f1, -f1], [f2, -f2], [f3, -f3], [f4, -f4], [f5, -f5]])
model.babylonjs()

# Fetching baybylon data to put custom labels
babylon_data = model.babylon_data()

for i, d in enumerate(babylon_data['steps']):
    d['label'] = 'custom label {}'.format(i+1)
    
model.babylonjs_from_babylon_data(babylon_data, debug=True)