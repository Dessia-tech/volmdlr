#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import io
import os
# import volmdlr as vm
import volmdlr.step

# import volmdlr.cloud as vmcd


for step_file in [
    #CylindricalSurfaces
    # 'cylinder-test.step',
    # 'read_test1.step'
    # 'read_test2.step',
    # 'read_test3.step',
    # 'read_test4.step',
    # 'read_test8.step',
    # 'read_test10.step',
    # 'angle_bar.step',
    # 'pipe.step',
    'bracket2_cut3.step',
    # 'bracket2_cut4.step',
    # 'bracket1.step'
    # 'Video_Version.step'
    #ToroidalSurface
    # 'read_test9.step',
    # 'read_test15.step',
    # 'tore1.step',
    # 'iter8finaldesign.stp',
    # 'cone1.step',
    # 'cone2.step',
    #bspline
    # 'Hollow_Loft.step'
    # # 'block.step',
    # # 'read_test6.step',
    # 'dimond.step'
    # 'bottle.step' # not implemented OFFSET_SURFACE

    # '2_bspline_faces.stp'# Uncomment when bug of delta fixed!
]:
    print('Reading step file: ', step_file)
    # filepath = os.path.join('step', step_file)
    step = volmdlr.step.Step.from_file(filepath=step_file)

    model = step.to_volume_model()
    faces = []
    closedshell = model.primitives[0]
    error = {}
    contours = {}
    for i, face in enumerate(closedshell.faces):
        try:
            face.triangulation()
        except Exception:
            error[i] = face
            contours[i] = face.surface2d.outer_contour
        else:
            faces.append(face)
    # for primitive in model.primitives:
    #     # primitive.color = (1, 0.2, 0.1)
    #     primitive.alpha = 0.98
    model2 = volmdlr.core.VolumeModel(faces)
    model2.babylonjs()
    # model.babylonjs()
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