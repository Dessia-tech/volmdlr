#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun  4 15:15:37 2021

@author: dasilva
"""


import volmdlr.cloud
import volmdlr.core
import volmdlr as vm
import volmdlr.wires as vmw
import volmdlr.faces as vmf
import volmdlr.edges as vme
import matplotlib.pyplot as plt

points_poly1 = [vm.Point3D(1,0,0),
            vm.Point3D(0.8660,0.5,0),
            vm.Point3D(0.7071,0.7071,0),
            vm.Point3D(0.5,0.8660,0),
            vm.Point3D(0,1,0),
            vm.Point3D(-0.5,0.8660,0),
            vm.Point3D(-0.7071,0.7071,0),
            vm.Point3D(-0.8660,0.5,0),
            vm.Point3D(-1,0,0),
            vm.Point3D(-0.8660,-0.5,0),
            vm.Point3D(-0.7071,-0.7071,0),
            vm.Point3D(-0.5,-0.8660,0),
            vm.Point3D(0,-1,0),
            vm.Point3D(0.5,-0.8660,0),
            vm.Point3D(0.7071,-0.7071,0),
            vm.Point3D(0.8660,-0.5,0)]
points_poly2=[vm.Point3D(-0.25,0.33,0.5),
           vm.Point3D(-0.33,-0.25,0.5),
           vm.Point3D(0.3,0.1,0.5)]

polygon1 = vmw.ClosedPolygon3D(points_poly1)
polygon2 = vmw.ClosedPolygon3D(points_poly2)


vec1, vec2 = vm.X3D, vm.Y3D
normal = vm.Z3D
faces=[]
coords = polygon1.sewing_with(polygon2, vec1, vec2, normal)
# coords = polygon1.sewing(polygon2, vec1, vec2, normal)
# for trio in coords :
#     faces.append(vmf.Triangle3D(trio[0], trio[1], trio[2]))

# volum = volmdlr.core.VolumeModel(faces)
# volum.babylonjs()  

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
polygon1.plot(ax=ax, color='g')
polygon2.plot(ax=ax, color='r')
    # volum = volmdlr.core.VolumeModel(cloud_faces)
    # faces.extend(cloud_faces)
    # volum.save_to_file(file)