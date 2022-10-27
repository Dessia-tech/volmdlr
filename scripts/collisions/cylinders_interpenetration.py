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
        position=vm.Point3D(0, 0, 0),
        axis=vm.Vector3D(1, 0, 0),
        radius=0.02,
        length=0.1,
        color=(0, 0, 1),
    ),
    Cylinder(
        position=vm.Point3D(0.05, 0, 0),
        axis=vm.Vector3D(0, 1, 0),
        radius=0.002,
        length=0.01,
        color=(1, 0, 1),
    ),
]

volume_model = vm.core.VolumeModel(cylinders)
volume_model.babylonjs()

print(
    f"""
Purple & blue:
min distance computed is {cylinders[0].min_distance_to_other_cylinder(cylinders[1])}m, 
collision: {cylinders[0].is_intersecting_other_cylinder(cylinders[0])}"""
)


print("Interpenetration")
start = time()

for _ in range(100):
    print(
        f"interpenetration volume: {cylinders[1].interference_volume_with_other_cylinder(cylinders[0], n_points=1000)}"
    )


end = time()
print(f"\nTotal interpenetration computation time: {end - start}s")
print(f"Time per calculus: {(end - start)/100}s")
