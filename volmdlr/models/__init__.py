#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 11:03:00 2020

@author: masfaraud
"""

from .electric_machine import all_rotor_mesh
from .electric_machine import stator_mesh

all_rotor_triangle_elements=all_rotor_mesh.generate_mesh(4)

stator_triangle_elements=stator_mesh.generate_mesh(4)