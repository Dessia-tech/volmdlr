"""
Test code for casing module
Generate random cylinders and create the casing for them
"""

from volmdlr.primitives3d import Cylinder
import volmdlr as vm

cylinders = [Cylinder(position=vm.O3D, axis=vm.X3D, radius=.01, length=.1),
             Cylinder(position=vm.O3D, axis=vm.Y3D, radius=.01, length=.1),
             Cylinder(position=vm.Point3D(0, 0.1, 0), axis=vm.X3D, radius=.01, length=.1)]

volume_model = vm.core.VolumeModel(cylinders)
volume_model.babylonjs()
