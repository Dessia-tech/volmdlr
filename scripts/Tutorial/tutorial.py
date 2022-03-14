# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 11:22:27 2020

@author: Mack Pro
"""

import numpy as npy
import volmdlr as vm
import volmdlr.primitives3d as primitives3D
import volmdlr.primitives2d as primitives2D
import matplotlib.pyplot as plt

p0 = vm.Point2D(0,0)
p1 = vm.Point2D(0.02, 0) 
p2 = vm.Point2D(0.02, -0.025)
p3 = vm.Point2D(0.230, -0.025)
p4 = vm.Point2D(0.230, 0)
p5 = vm.Point2D(0.285, 0.055)
p6 = vm.Point2D(0.285, 0.065)
p7 = vm.Point2D(0.275, 0.075)
p8 = vm.Point2D(0.225, 0.175)
p9 = vm.Point2D(0.225, 0.3)
p10 = vm.Point2D(0.09, 0.3)
p11 = vm.Point2D(0.09, 0.245)
p12 = vm.Point2D(0.08, 0.235)
p13 = vm.Point2D(0.08, 0.205)
p14 = vm.Point2D(0.02, 0.150)
p15 = vm.Point2D(0.03, 0.05)

points1 = [p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15]
c1 = primitives2D.ClosedRoundedLineSegments2D(points1, {2: 0.005, 3: 0.005, 4: 0.005, 8: 0.005, 9: 0.005, 10: 0.005, 14: 0.005}) 



circle0 = vm.wires.Circle2D(vm.Point2D(0.12, 0.25), 0.02)
circle1 = vm.wires.Circle2D(vm.Point2D(0.19, 0.25), 0.02)
circle2 = vm.wires.Circle2D(vm.Point2D(0.15, 0.08), 0.08)


c2 = vm.wires.Contour2D([circle0])
c3 = vm.wires.Contour2D([circle1])
c4 = vm.wires.Contour2D([circle2])


e0 = vm.Point2D(0, 0.00)
e1 = vm.Point2D(0, 0.03)
e2 = vm.Point2D(0.11, 0.03)
e3 = vm.Point2D(0.11, 0.04)
e4 = vm.Point2D(0.13, 0.04)
e5 = vm.Point2D(0.13, 0)

points2 = [e0,e1,e2,e3,e4,e5]
c5 = primitives2D.ClosedRoundedLineSegments2D(points2, {})




profile=primitives3D.ExtrudedProfile(vm.O3D, vm.Y3D, vm.Z3D, c1, [], vm.X3D*0.40, name = 'extrusion')
profile1=primitives3D.Cylinder(position = vm.Point3D(-0.0075, 0.14, 0.25), axis = vm.X3D, 
                               radius = 0.02, length = 0.015,
                               name = 'profile1')

profile2=primitives3D.Cylinder(position = vm.Point3D(-0.0075, 0.19, 0.25), axis = vm.X3D, 
                               radius = 0.02, length = 0.015,
                               name = 'profile2')

profile3=primitives3D.Cylinder(position = vm.Point3D(-0.0075, 0.15, 0.08), axis = vm.X3D, 
                               radius = 0.08, length = 0.015,
                               name = 'profile3')



y = vm.X3D.random_unit_normal_vector()
z = vm.X3D.cross(y)
profile4=primitives3D.RevolvedProfile(vm.Z3D*0.08-0.015*vm.Y3D, vm.X3D, z, c5, vm.Z3D*0.08-0.015*vm.Y3D , vm.X3D)
profile5=primitives3D.RevolvedProfile(vm.Z3D*0.22+0.035*vm.Y3D, vm.X3D, z, c5, vm.Z3D*0.22+0.035*vm.Y3D , vm.X3D)

c = vm.wires.Circle2D(vm.Point2D(0,0), 0.008)
pt0 = vm.Point3D(0.01, 0.04, 0.16)
pt1 = vm.Point3D(0.03, 0, 0.2)
pt2 = vm.Point3D(0.45, 0.01, 0.1)
pt3 = vm.Point3D(0.45, 0, -0.1)
pt4 = vm.Point3D(0.3, 0.04, -0.02)
pts = [pt0, pt1, pt2, pt3, pt4]
radius = {1: 0.03, 2: 0.01, 3: 0.07}
rl = primitives3D.OpenRoundedLineSegments3D(pts, radius, adapt_radius=True, name='wire')
sweep = primitives3D.Sweep(c, rl, name = 'pipe')


pt10 = vm.Point3D(0.02, 0.22, 0.25)
pt11 = vm.Point3D(0.02, 0.24, 0.25)
pt12 = vm.Point3D(0.6, 0.24, 0.20)
pt13 = vm.Point3D(0.40, 0.17, 0.13)
pts1 = [pt10, pt11, pt12, pt13]
radius1 = {1: 0.01, 2: 0.05}
rl1 = primitives3D.OpenRoundedLineSegments3D(pts1, radius1, adapt_radius=True, name='wire1')
sweep1 = primitives3D.Sweep(contour2d = c, 
                            wire3d = rl1, name = 'pipe1')



model=vm.core.VolumeModel([profile, profile1, profile2, profile3, profile4, profile5, sweep, sweep1])

model.babylonjs()

# from dessia_api_client import Client
# import json

# C = Client()
# C.api_url = 'https://api.platform.dessia.tech'

# d = model.to_dict()
# j = json.dumps(d)
# d2 = json.loads(j)
# v2 = vm.VolumeModel.dict_to_object(d2)

# r = C.create_object_from_python_object(model)