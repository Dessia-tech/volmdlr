#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 17:26:59 2020

@author: mehdigasmi
"""

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as npy
import volmdlr.edges as edges
import volmdlr.wires as wires
import volmdlr.faces as faces
import volmdlr.core_compiled as vm
import volmdlr.mesh as vmmesh


rotor_external=wires.Contour2D.load_from_file('/Users/mehdigasmi/Desktop/eletric_machine/rotor_external_contour.json') 
rotor_internal=wires.Contour2D.load_from_file('/Users/mehdigasmi/Desktop/eletric_machine/rotor_internal_contour.json')
rotor_magnet=wires.Contour2D.load_from_file('/Users/mehdigasmi/Desktop/eletric_machine/rotor_magnet_contour.json')
stator_external=wires.Contour2D.load_from_file('/Users/mehdigasmi/Desktop/eletric_machine/stator_external_contour.json')
stator_internal=wires.Contour2D.load_from_file('/Users/mehdigasmi/Desktop/eletric_machine/stator_internal_contour.json')
rotor_internal=rotor_internal.primitives[0]
stator_external=stator_external.primitives[0]


all_rotor_mesh=vmmesh.Mesher([rotor_internal],[rotor_magnet,rotor_external],60)

surface=faces.Surface2D(stator_external,[stator_internal])
pattern=surface.get_pattern_single_inner()
all_patterns=surface.contour_from_pattern()
stator_mesh=vmmesh.Mesher([],all_patterns,60)
