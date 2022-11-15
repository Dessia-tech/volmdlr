"""
Test for method Cylinder.point_belongs and Cylinder.random_point_inside
"""
import volmdlr

from volmdlr.primitives3d import Cylinder

cylinder = Cylinder(position=volmdlr.O3D, axis=volmdlr.Z3D, radius=0.02, length=0.1)

p = cylinder.random_point_inside()

if not cylinder.point_belongs(p):
    raise ValueError("Point should belong to cylinder")
