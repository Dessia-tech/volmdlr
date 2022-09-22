#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 21 2022

@author: S.Bendjebla
"""

import volmdlr as vm
import volmdlr.wires
import volmdlr.faces
import matplotlib.pyplot as plt
import volmdlr.gmsh

# %% Contour2D

p = [vm.Point2D(-0.3, -0.2), vm.Point2D(0.3, -0.2),
      vm.Point2D(0.2, 0.2), vm.Point2D(0, 0.3), vm.Point2D(-0.2, 0.2)]

contour = vm.wires.Contour2D.from_points(p)

# %% Surface2D

surface = vm.faces.Surface2D(contour, [])


# %% Create .msh file

file_name, mesh_dimension, factor = 'mesh_surface', 2, 0.5

surface.to_msh(file_name, mesh_dimension, factor) #mesh_dimension = 2 (triangles), factor: mesh size between 0, 1

# %% Mesh

gmsh = volmdlr.gmsh.Gmsh.from_file(file_name+'.msh')

mesh = gmsh.define_triangular_element_mesh()

ax = mesh.plot()