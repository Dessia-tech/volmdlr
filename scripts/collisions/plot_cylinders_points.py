"""
Test script for cylinder discretization
"""

import volmdlr as vm
from volmdlr.primitives3d import Cylinder

c = Cylinder(
        position=vm.Point3D(0, 0.1, 0),
        axis=vm.Vector3D(1, 0, 0),
        radius=0.01,
        length=0.1,
        color=(1, 0, 0),
)

ay = c.random_point_inside().plot()
for _ in range(5000):
    c.random_point_inside().plot(ax=ay)
