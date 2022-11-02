"""
Test code for intersections between cylinders
Generate random cylinders and create the casing for them
"""

from time import time
from volmdlr.primitives3d import Cylinder
import volmdlr.core
import volmdlr as vm
import numpy


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
        radius=0.005,
        length=0.01,
        color=(1, 0, 1),
    ),
]

volume_model = vm.core.VolumeModel(cylinders)
# volume_model.babylonjs()

print(
    f"""
Purple & blue:
min distance computed is {cylinders[0].min_distance_to_other_cylinder(cylinders[1])}m, 
collision: {cylinders[0].is_intersecting_other_cylinder(cylinders[0])}"""
)


print("Interpenetration")
start = time()

volumes = []
n = 100

for _ in range(n):
    volumes.append(cylinders[1].interference_volume_with_other_cylinder(cylinders[0], n_points=300))
    print(
        f"interpenetration volume: {volumes[-1]}"
    )

end = time()
print(f"\nTotal interpenetration computation time: {end - start}s")
print(f"Time per calculus: {(end - start)/100}s")
print(f"Mean volume : {sum(volumes) / n} m³")
print(f'Standard deviation : {numpy.std(volumes)}')
print(f'Variation coefficient : {(numpy.std(volumes) / (sum(volumes) / n))*100} %')
