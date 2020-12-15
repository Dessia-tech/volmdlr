#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 11:46:37 2020

@author: milan gasmi
"""


from volmdlr.models import stator_mesher,rotor_mesher,pattern_mesher


rotor_mesh=rotor_mesher
rotor_triangle_elements=rotor_mesh.generate_mesh(4)
print(rotor_mesh.check_conformity(rotor_triangle_elements))

pattern_mesh=pattern_mesher
pattern_triangle_elements=pattern_mesh.generate_mesh(4)

stator_mesh=stator_mesher
stator_triangle_elements=stator_mesh.generate_mesh(4)