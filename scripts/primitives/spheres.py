"""
Demo script of spheres displayed on the skin of a defined sphere.
"""
import math
import random

import volmdlr as vm
import volmdlr.core as vmc
import volmdlr.primitives3d as p3d

center = vm.Point3D.random(0, 0.1, 0, 0.1, -0.1, 0)
radius = random.randrange(5, 50, 5)
radius *= 1e-3


sphere = p3d.Sphere(center=center, radius=radius, color=(0.25, 0.33, 0.6), alpha=0.5)

resolution = max([math.pi * radius / 10, 1e-3])

skin_points = sphere.skin_points(resolution=resolution)

spheres = [sphere]
for pt in skin_points:
    spheres.append(p3d.Sphere(center=pt, radius=resolution / 4, color=(1, 0, 0)))


points_in = sphere.inner_points(resolution=resolution)

for pt in points_in:
    spheres.append(p3d.Sphere(center=pt, radius=resolution / 4, color=(0, 1, 0)))

# all_points = point_skin + point_in
# print(len(all_points))

# vol = vmc.VolumeModel(spheres)
# vol.babylonjs()
