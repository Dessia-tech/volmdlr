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

print("Collision detection methods\n")
start = time()

print(
    f"""\nRed & green:
min distance is {cylinders[0].min_distance_to_other_cylinder(cylinders[1])}m, 
collision: {cylinders[0].is_intersecting_other_cylinder(cylinders[1])}"""
)

print(
    f"""\nRed & blue:
min distance is {cylinders[0].min_distance_to_other_cylinder(cylinders[2])}m, 
collision: {cylinders[0].is_intersecting_other_cylinder(cylinders[2])}"""
)

print(
    f"""\nGreen & blue:
min distance is {cylinders[1].min_distance_to_other_cylinder(cylinders[2])}m, 
collision: {cylinders[1].is_intersecting_other_cylinder(cylinders[2])}"""
)

print(f"\nCollision detection computing time: {time() - start}s")
