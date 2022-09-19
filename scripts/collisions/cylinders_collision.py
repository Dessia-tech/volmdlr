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
        color=(1, 0, 0),
    ),
    Cylinder(
        position=vm.Point3D(0, 0.05, 0),
        axis=vm.Vector3D(1, 1, 0),
        radius=0.005,
        length=0.1,
        color=(0, 1, 0),
    ),
    Cylinder(
        position=vm.Point3D(0, 0.159, 0),
        axis=vm.Vector3D(1, 1, 0),
        radius=0.02,
        length=0.1,
        color=(0, 0, 1),
    ),
    Cylinder(
        position=vm.Point3D(0, 0.1, 0.016),
        axis=vm.Vector3D(0, 1, 0),
        radius=0.01,
        length=0.1,
        color=(1, 0, 1),
    ),
]

volume_model = vm.core.VolumeModel(cylinders)
volume_model.babylonjs()

print("Collision detection methods")
start = time()

print(
    f"""\nRed & green:
min distance computed is {cylinders[0].min_distance_to_other_cylinder(cylinders[1])}m, 
collision: {cylinders[0].is_intersecting_other_cylinder(cylinders[1])}"""
)

print(
    f"""\nRed & blue:
min distance computed is {cylinders[0].min_distance_to_other_cylinder(cylinders[2])}m, 
collision: {cylinders[0].is_intersecting_other_cylinder(cylinders[2])}"""
)

print(
    f"""\nGreen & blue:
min distance computed is {cylinders[1].min_distance_to_other_cylinder(cylinders[2])}m, 
collision: {cylinders[1].is_intersecting_other_cylinder(cylinders[2])}"""
)

print(
    f"""\nRed & purple:
min distance computed is {cylinders[0].min_distance_to_other_cylinder(cylinders[3])}m, 
collision: {cylinders[0].is_intersecting_other_cylinder(cylinders[3])}"""
)

print(
    f"""\nGreen & purple:
min distance computed is {cylinders[1].min_distance_to_other_cylinder(cylinders[3])}m, 
collision: {cylinders[1].is_intersecting_other_cylinder(cylinders[3])}"""
)

print(
    f"""\nPurple & blue:
min distance computed is {cylinders[3].min_distance_to_other_cylinder(cylinders[2])}m, 
collision: {cylinders[3].is_intersecting_other_cylinder(cylinders[2])}"""
)

end = time()
print(f"\nTotal collision detection computation time: {end - start}s")
print(f"Time per collision: {(end - start)/6}s")
