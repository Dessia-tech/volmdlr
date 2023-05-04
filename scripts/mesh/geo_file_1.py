#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 2022

@author: s.bendjebla
"""

import volmdlr
import volmdlr as vm
import volmdlr.primitives3d as primitives3d

# import gmsh

# %% Extrusion

points = [vm.Point2D(0, 0),
          vm.Point2D(0.1, 0.),
          vm.Point2D(0.1, 0.2),
          vm.Point2D(0.03, 0.15),
          vm.Point2D(0.,0.21)]

outer_profile = vm.wires.Contour2D.from_points(points)

profile=primitives3d.ExtrudedProfile(vm.O3D, vm.Y3D, vm.Z3D, outer_profile, [], vm.X3D*0.1, name = 'extrusion')

model=vm.core.VolumeModel([profile])
model.to_geo(file_name = 'model_1_geo',
             factor = 0.5,
             curvature_mesh_size = 0,
             min_points = None,
             initial_mesh_size = 5)

# %% gmsh file generation

# gmsh.initialize()
# gmsh.open("model.geo")

# gmsh.model.geo.synchronize()
# gmsh.model.mesh.generate(2)

# gmsh.write("model.msh")

# gmsh.finalize()

# %% DIRECT: gmsh file generation

# %%% 2D & Order 1

model.to_msh(file_name = 'mesh_2d_order_1',
             mesh_dimension = 2,
             factor = 1,
             curvature_mesh_size = 0,
             min_points = None,
             initial_mesh_size = 5)

# %%% 3D & Order 1

model.to_msh(file_name = 'mesh_3d_order_1',
             mesh_dimension = 3,
             mesh_order = 1,
             factor = 1,
             curvature_mesh_size = 0,
             min_points = None,
             initial_mesh_size = 5)

# %%% 3D & Order 2

model.to_msh(file_name = 'mesh_3d_order_2',
             mesh_dimension = 3,
             mesh_order = 2,
             factor = 1,
             curvature_mesh_size = 0,
             min_points = None,
             initial_mesh_size = 5)
