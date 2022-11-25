#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 2022

@author: s.bendjebla
"""

import volmdlr
import volmdlr as vm
import volmdlr.primitives3d as primitives3d
import volmdlr.edges
# import gmsh

# %% Extrusion

edges = [volmdlr.edges.LineSegment2D(vm.Point2D(0, 0), vm.Point2D(0.1, 0.)),
         volmdlr.edges.LineSegment2D(vm.Point2D(0.1, 0.), vm.Point2D(0.1, 0.2)),
         volmdlr.edges.Arc2D(start=vm.Point2D(0.1, 0.2),
                             end=vm.Point2D(0.,0.2),
                             interior=vm.Point2D(0.05,0.25)),
         volmdlr.edges.LineSegment2D(vm.Point2D(0.,0.2), vm.Point2D(0, 0))]

outer_profile = vm.wires.Contour2D(edges)


profile=primitives3d.ExtrudedProfile(vm.O3D, vm.Y3D, vm.Z3D, outer_profile, [], vm.X3D*0.1, name = 'extrusion')

model=vm.core.VolumeModel([profile])

model.to_geo(file_name = 'model_2_geo',
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

model.to_msh(file_name = 'model_2',
             mesh_dimension = 2,
             factor = 1,
             curvature_mesh_size = 0,
             min_points = None,
             initial_mesh_size = 5)
