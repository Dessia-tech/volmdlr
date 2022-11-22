#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import io
import os
# import volmdlr as vm
import volmdlr.step
import volmdlr.faces as vmf

# import volmdlr.cloud as vmcd


for step_file in [
    # 'cylinder-test.step',
    # 'bracket2.step',
    # 'read_test1.step',
    # 'read_test2.step',
    # 'read_test3.step',
    # 'read_test4.step',
    # 'read_test6.step',
    # 'read_test7.step',
    # 'iso4162M16x55.step',
    # 'spanners.step',
    # 'tormach_wrench.step',
    # 'water_tank.step',
    # 'angle_bar.step',
    # 'tore1.step',
    # 'block.step',
    # 'Hollow_Loft.step',
    # 'spherical_surface_body.step',
    # 'bracket2_cut3.step',
    # 'demi_sphere.step',
    # 'Video_Version.step',
    # 'bracket3-cut1.step',
    # 'quart_cone.step',
    # 'demi_cone.step',
    # 'cone1.step',
    # 'cone2.step',
    # 'HRG_BOT.stp',
    # 'bracket3.step',
    # 'strange_gasket.,
    # '2_bspline_faces.stp',
    'bracket1.step'

]:

    print('Reading step file: ', step_file)
    # filepath = os.path.join('step', step_file)
    step = volmdlr.step.Step.from_file(filepath=step_file)

    model = step.to_volume_model()


    # assert len(model.primitives) > 0.
    # model.to_step(step_file+'_reexport')
    model.primitives[0].alpha = 1
    model.primitives[0].color = (1, 0.1, 0.1)

    model.babylonjs()
    # faces = []
    # error = {}
    # contours = {}
    # for closedshell in model.primitives:
    #     for i, face in enumerate(closedshell.faces):
    #         try:
    #             face.triangulation()
    #         except Exception:
    #             error[i] = face
    #             contours[i] = face.surface2d
    #         else:
    #             faces.append(face)
    # model2 = volmdlr.core.VolumeModel(faces)
    # model2.babylonjs()
    # assert len(model.primitives) > 0.
    # model.to_step(step_file + '_reexport')
    # model.babylonjs()
    # file_io = io.FileIO(step_file, 'r')
    # step = volmdlr.step.Step.from_stream(stream=file_io)
    # model = step.to_volume_model()
    # assert len(model.primitives) > 0.
    # model.to_step(step_file + '_reexport')
    #
    # model2 = model.copy()
    #
    # # model2 = model.copy()
    # # assert model == model2
    #
    # model._check_platform()

    # closedshell = model.primitives[0]
    # ax0 = closedshell.faces[2].plot()
    # ax = closedshell.faces[2].surface2d.plot()
    # ax1 = closedshell.faces[2].triangulation().plot()