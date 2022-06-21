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

gmsh = volmdlr.gmsh.Gmsh.from_file(file_path)

mesh = gmsh.define_triangular_element_mesh()

mesh.plot()

# %% 3D

file_path = 'block.msh'

gmsh = volmdlr.gmsh.Gmsh.from_file(file_path)

mesh = gmsh.define_tetrahedron_element_mesh()
