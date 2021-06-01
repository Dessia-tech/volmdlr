#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr as vm
import volmdlr.step
import volmdlr.cloud as vmcd

faces = []
for step_file in [
                  # 'tore1.step',
                  # 'cone1.step',
                  # 'cone2.step',
                  # 'cylinder.step',
                    'STEP_test1.stp',
                  # 'block.step',
                   # 'iso4162M16x55.step',
                   # 'aircraft_engine.step'
                   # 'OCIM0BDME8.stp',
                  # 'boite_vitesse.stp'
                  'engine_body.stp'
                  ]:
    print('filename: ', step_file)
    # step = volmdlr.step.Step(step_file)
    # model = step.to_volume_model()
    # model.to_step(step_file+'_reexport')
    # print(model.primitives)

    # model.babylonjs()
    
    cloud = vmcd.PointCloud3D.from_step(step_file)
    to_add = cloud.subdescription_2d()
    faces.extend(to_add)
volum = volmdlr.core.VolumeModel(faces)

volum.babylonjs()
    
# model2 = model.copy()

# assert model == model2