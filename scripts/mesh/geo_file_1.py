#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 2022

@author: s.bendjebla
"""

import volmdlr
import volmdlr as vm
import volmdlr.primitives3d as primitives3d
import gmsh

# %% Extrusion

points = [vm.Point2D(0, 0),
          vm.Point2D(0.1, 0.),
          vm.Point2D(0.1, 0.2),
          vm.Point2D(0.03, 0.15),
          vm.Point2D(0.,0.21)]

outer_profile = vm.wires.Contour2D.from_points(points)

profile=primitives3d.ExtrudedProfile(vm.O3D, vm.Y3D, vm.Z3D, outer_profile, [], vm.X3D*0.1, name = 'extrusion')

model=vm.core.VolumeModel([profile])
model.to_geo('model')

# %% gmsh file generation

# gmsh.initialize()
# gmsh.open("model.geo")

# gmsh.model.geo.synchronize()
# gmsh.model.mesh.generate(2)

# gmsh.write("model.msh")

# gmsh.finalize()

# %% DIRECT: gmsh file generation

model.to_msh('model_msh')
