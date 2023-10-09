 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import math

import volmdlr as vm
import volmdlr.wires as vmw
from volmdlr.primitives2d import (ClosedRoundedLineSegments2D,
                                  OpenedRoundedLineSegments2D)
from volmdlr.primitives3d import RevolvedProfile

p1 = vm.Point2D(-0.08, 0.1755)
p2 = vm.Point2D(0., 0.1755)
p3 = vm.Point2D(0.0, 0.04)
p4 = vm.Point2D(0.01, 0.04)
p5 = vm.Point2D(0.01, 0.)
p6 = vm.Point2D(-0.06, 0.0)
# p7 = vm.Point2D(-0.06, 0.04)
p8 = vm.Point2D(-0.08, 0.06)

points1 = [p1, p2, p3, p4, p5, p6, p8]

c1 = vm.wires.ClosedPolygon2D(points1)
c2 = ClosedRoundedLineSegments2D(points1, {0: 0.002, 1: 0.002, 3: 0.002}) 

# c1.MPLPlot(plot_points=True)

# c2.MPLPlot(plot_points=True)

p21 = vm.Point2D(0.1, 0.)
p22 = vm.Point2D(0.1, 0.1)
p23 = vm.Point2D(0.08, 0.1)
p24 = vm.Point2D(0.08, 0.2)
p25 = vm.Point2D(0., 0.2)
p26 = vm.Point2D(0., 0.05)
p27 = vm.Point2D(0.04, 0.05)
p28 = vm.Point2D(0.04, 0.)

points2 = [p21, p22, p23, p24, p25, p26, p27, p28]
c3 = vm.wires.ClosedPolygon2D(points2)
c4 = ClosedRoundedLineSegments2D(points2, {0: 0.002, 1: 0.002, 3: 0.002})

# c3.MPLPlot(plot_points=True)
# c4.MPLPlot(plot_points=True)
frame = vm.Frame3D(vm.Point3D(0, 1, 0), vm.Y3D, vm.Z3D, -vm.X3D)
profile1 = RevolvedProfile(frame, c1, vm.O3D, vm.Y3D, 3)
# profile2 = RevolvedProfile(0.5*vm.X3D, vm.Y3D, vm.Z3D, c2, 0.5*vm.X3D, vm.Y3D)
# profile3 = RevolvedProfile(-0.5*vm.X3D, vm.Y3D, vm.Z3D, c1, -0.5*vm.X3D, vm.Y3D)
# profile4 = RevolvedProfile(vm.X3D, vm.Y3D, vm.Z3D, c2, vm.X3D, vm.Y3D)


# B = 0.020
# d1 = 0.015
# h = 0.005
# radius = 0.002
# F = 0.025
# d = 0.010
B = 0.057
d1 = 0.45200000000000007
h = 0.007778409372698711
radius = 0.005 #with 0.005 it's better to debug
F = 0.42500000000000004
d = 0.38


# Internal ring contour
pbi2 = vm.Point2D(-B/2., d1/2.)
pbi1 = pbi2.translation(vm.Vector2D(h, 0))
pbi3 = vm.Point2D(-B/2., d/2.)
pbi4 = vm.Point2D(B/2., d/2.)
pbi5 = vm.Point2D(B/2., d1/2.)
pbi6 = pbi5.translation(vm.Vector2D(-h, 0))
bi1 = OpenedRoundedLineSegments2D([pbi6, pbi5, pbi4, pbi3, pbi2, pbi1],
                                          {1: radius,
                                          2: radius,
                                          3: radius,
                                          4: radius},
                                          adapt_radius=True)
cbi1 = vm.edges.Arc2D.from_3_points(pbi1, vm.Point2D(0, F/2), pbi6)
c5 = vm.wires.Contour2D([cbi1] + bi1.primitives)
# c5.MPLPlot(plot_points=True)

y = vm.X3D.random_unit_normal_vector()
z = vm.X3D.cross(y)
frame2 = vm.Frame3D(0.15*vm.Y3D.to_point(), vm.X3D, z, vm.X3D.cross(z))
profile5 = RevolvedProfile(frame2, c5, 0.15*vm.Y3D.to_point(), vm.X3D,
                           angle=0.7, name='strange part')

# model = vm.VolumeModel([profile1, profile2, profile3, profile4, profile5])

r = 0.15
R = 0.2
Rb = 0.25
w = 0.2
wb = 0.02
th = 0.008

rim_contour = vmw.ClosedPolygon2D([vm.Point2D(-0.5*w, Rb),
                                   vm.Point2D(-0.5*w+wb, Rb),
                                   vm.Point2D(-0.5*w+wb, R),
                                   vm.Point2D(-0.05 * w, r),
                                   vm.Point2D(0, r),
                                   vm.Point2D(0, r-th),
                                   vm.Point2D(-0.05 * w, r - th),
                                   vm.Point2D(-0.5*w, R-th),
                                   ])

rim_contour.plot()
frame2 = vm.Frame3D(0.5*vm.X3D.to_point(), vm.X3D, z, vm.X3D.cross(z))
conical_rim = RevolvedProfile(frame2, rim_contour, 0.5*vm.X3D.to_point(), vm.X3D,
                              angle=3.14, name='conical rim')


model = vm.core.VolumeModel([profile1, profile5, conical_rim])

translated_model = model.translation(vm.Point3D.random(0, 1, 0, 1, 0, 1))
turned_model = model.rotation(0.3*vm.X3D.to_point(), vm.Z3D, 0.4)

# model.primitives.extend(translated_model.primitives+turned_model.primitives)
model.babylonjs()
# translated_model.babylonjs()
# turned_model.babylonjs()
model.to_step('revolved_profile.step')

# Checking face triangulation
planar_face = profile1.faces[1]

planar_face_triangulation = planar_face.surface2d.triangulation
assert math.isclose(planar_face_triangulation().area(), planar_face.surface2d.area(), abs_tol=1e-3)

model_copy = model.copy()
assert model_copy == model
