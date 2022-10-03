#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 2022

@author: s.bendjebla
"""

import volmdlr
import volmdlr.gmsh_vm

# %% 2D

file_path = 'model.msh'

mesh_gmsh = volmdlr.gmsh_vm.Gmsh.from_file(file_path)

mesh = mesh_gmsh.define_triangular_element_mesh()

mesh.plot()

# %% 3D

file_path = 'block.msh'

mesh_gmsh = volmdlr.gmsh_vm.Gmsh.from_file(file_path)

mesh = mesh_gmsh.define_tetrahedron_element_mesh()

# mesh.plot()
