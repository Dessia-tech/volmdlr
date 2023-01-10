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

f2 = f1.rotation(vm.Z3D, 0.1)
f2.translation_inplace(0.1*vm.X3D)

f3 = f2.rotation(vm.Z3D, 0.1)
f3.translation_inplace(0.1*vm.X3D)

f4 = f3.rotation(vm.Z3D, 0.1)
f4.translation_inplace(0.1*vm.X3D)

f5 = f4.rotation(vm.Z3D, 0.1)
f5.translation_inplace(0.1*vm.X3D)


model = vm.MovingVolumeModel([block1, block2], [[f1, -f1], [f2, -f2], [f3, -f3], [f4, -f4], [f5, -f5]])
model.babylonjs()

# Fetching baybylon data to put custom labels
babylon_data = model.babylon_data()

for i, d in enumerate(babylon_data['steps']):
    d['label'] = 'custom label {}'.format(i+1)
    
model.babylonjs_from_babylon_data(babylon_data)