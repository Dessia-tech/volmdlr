#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 2022

@author: s.bendjebla
"""

import volmdlr
import volmdlr.gmsh

# %% 2D

file_path = 'model.msh'

gmsh_parser = volmdlr.gmsh.GmshParser.from_file(file_path)

mesh = gmsh_parser.define_triangular_element_mesh()

mesh.plot()

# %% 3D

file_path = 'block.msh'

gmsh_parser = volmdlr.gmsh.GmshParser.from_file(file_path)

mesh = gmsh_parser.define_tetrahedron_element_mesh()

# mesh.plot()
