# -*- coding: utf-8 -*-


import math
import random

import volmdlr as vm
import volmdlr.core as vmc
import volmdlr.primitives3d as p3d

center = vm.Point3D.random(0, 0.1, 0, 0.1, -0.1, 0)
radius = random.randrange(5,500,5)
radius *= 1e-3


sphere = p3d.Sphere(center = center, 
                    radius = radius,
                    color = (0.25,0.33,0.6), 
                    alpha = 0.5)

resolution = max([math.pi*radius/10, 1e-3])


point_skin = sphere.to_point_skin(resolution = resolution)

spheres = [sphere]
for pt in point_skin :
    spheres.append(p3d.Sphere(center = pt, 
                              radius = resolution/4,
                              color = (1,0,0)))


point_in = sphere.to_point_in(resolution = resolution)

for pt in point_in :
    spheres.append(p3d.Sphere(center = pt, 
                              radius = resolution/4,
                              color = (0,1,0)))

# all_points = point_skin + point_in
# print(len(all_points))

# vol = vmc.VolumeModel(spheres)
# vol.babylonjs()
