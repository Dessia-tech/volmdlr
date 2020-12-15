#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 11:26:19 2020

@author: masfaraud
"""

import dessia_common.workflow 
import volmdlr.core
import volmdlr

packer = dessia_common.workflow.Sequence(3, float)

point3D_instanciation = dessia_common.workflow.InstanciateModel(
    volmdlr.Point3D)
pipe_vector = dessia_common.workflow.Pipe(packer.outputs[0],
                                          point3D_instanciation.inputs[0])
point3D_instanciator = dessia_common.workflow.Workflow(
    [packer, point3D_instanciation],
    [pipe_vector],
    point3D_instanciation.outputs[0],
    name='Point3D instanciator')
point3D_instanciator.inputs[0].name = 'X'
point3D_instanciator.inputs[1].name = 'Y'
point3D_instanciator.inputs[2].name = 'Z'
