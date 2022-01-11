# -*- coding: utf-8 -*-


import volmdlr.primitives3d as p3d
import volmdlr as vm

import volmdlr.core as vmc
import random
import math

center = vm.Point3D.random(0, 0.1, 0, 0.1, -0.1, 0)
radius = random.randrange(5,500,5)
radius *= 1e-3

sphere = p3d.Sphere(center = center, 
                    radius = radius,
                    color = (0.25,0.33,0.6), 
                    alpha = 0.5)

resolution = max([math.pi*radius/10, 1e-3])


pointcloud_skin = sphere.to_pointcloud3d_skin(resolution = resolution)

spheres = [sphere]
for pt in pointcloud_skin :
    spheres.append(p3d.Sphere(center = pt, 
                              radius = resolution/4,
                              color = (1,0,0)))


pointcloud_in = sphere.to_pointcloud3d_in(resolution = resolution)

for pt in pointcloud_in :
    spheres.append(p3d.Sphere(center = pt, 
                              radius = resolution/4,
                              color = (0,1,0)))


# vol = vmc.VolumeModel(spheres)
# vol.babylonjs()