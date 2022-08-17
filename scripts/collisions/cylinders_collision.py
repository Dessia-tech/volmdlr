"""
Test code for intersections between cylinders
Generate random cylinders and create the casing for them
"""

from time import time
from volmdlr.primitives3d import Cylinder
import volmdlr.core
import volmdlr as vm


cylinders = [
    Cylinder(
        position=vm.Point3D(0, 0.1, 0),
        axis=vm.Vector3D(1, 0, 0),
        radius=0.01,
        length=0.1,
    ),
    Cylinder(
        position=vm.Point3D(0, 0.05, 0),
        axis=vm.Vector3D(1, 1, 0),
        radius=0.005,
        length=0.1,
    ),
    Cylinder(
        position=vm.Point3D(0, 0.159, 0),
        axis=vm.Vector3D(1, 1, 0),
        radius=0.02,
        length=0.1,
    ),
]

cylinders[0].color = (1, 0, 0)
cylinders[1].color = (0, 1, 0)
cylinders[2].color = (0, 0, 1)

volume_model = vm.core.VolumeModel(cylinders)
volume_model.babylonjs()

print("Collision detection methods")
start = time()

print("Intersecting red & green:", cylinders[0].is_intersecting_other_cylinder(cylinders[1]))
print("Intersecting red & blue:", cylinders[0].is_intersecting_other_cylinder(cylinders[2]))
print("Intersecting green & blue:", cylinders[1].is_intersecting_other_cylinder(cylinders[2]))

print(f"\nIntersection computing duration: {time() - start}s")
