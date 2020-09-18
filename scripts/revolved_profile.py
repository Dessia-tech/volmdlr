 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr.core as vm
from volmdlr.primitives2D import ClosedRoundedLineSegments2D, OpenedRoundedLineSegments2D
from volmdlr.primitives3D import RevolvedProfile

p1 = vm.Point2D((-0.08, 0.1755))
p2 = vm.Point2D((0., 0.1755))
p3 = vm.Point2D((0.0, 0.04))
p4 = vm.Point2D((0.01, 0.04))
p5 = vm.Point2D((0.01, 0.))
p6 = vm.Point2D((-0.06, 0.0))
p7 = vm.Point2D((-0.06, 0.04))
p8 = vm.Point2D((-0.08, 0.04))

points1 = [p1, p2, p3, p4, p5, p6, p7, p8]

c1 = vm.Polygon2D(points1)
c2 = ClosedRoundedLineSegments2D(points1, {0: 0.002, 1: 0.002, 3: 0.002}) 

# c1.MPLPlot(plot_points=True)

# c2.MPLPlot(plot_points=True)

p21 = vm.Point2D((0.1, 0.))
p22 = vm.Point2D((0.1, 0.1))
p23 = vm.Point2D((0.08, 0.1))
p24 = vm.Point2D((0.08, 0.2))
p25 = vm.Point2D((0., 0.2))
p26 = vm.Point2D((0., 0.05))
p27 = vm.Point2D((0.04, 0.05))
p28 = vm.Point2D((0.04, 0.))

points2 = [p21, p22, p23, p24, p25, p26, p27, p28]
c3 = vm.Polygon2D(points2)
c4 = ClosedRoundedLineSegments2D(points2, {0: 0.002, 1: 0.002, 3: 0.002})

# c3.MPLPlot(plot_points=True)
# c4.MPLPlot(plot_points=True)

# profile1 = RevolvedProfile(vm.O3D, vm.Y3D, vm.Z3D, c1, vm.O3D, vm.Y3D)
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
pbi2 = vm.Point2D((-B/2., d1/2.))
pbi1 = pbi2.Translation(vm.Vector2D((h, 0)))
pbi3 = vm.Point2D((-B/2., d/2.))
pbi4 = vm.Point2D((B/2., d/2.))
pbi5 = vm.Point2D((B/2., d1/2.))
pbi6 = pbi5.Translation(vm.Vector2D((-h, 0)))
bi1 = OpenedRoundedLineSegments2D([pbi6, pbi5, pbi4, pbi3, pbi2, pbi1],
                                          {1: radius,
                                          2: radius,
                                          3: radius,
                                          4: radius},
                                          adapt_radius=True)
cbi1 = vm.Arc2D(pbi1, vm.Point2D((0, F/2)), pbi6)
c5 = vm.Contour2D([cbi1] + bi1.primitives)
# c5.MPLPlot(plot_points=True)

y = vm.X3D.RandomUnitNormalVector()
z = vm.X3D.Cross(y)
profile5 = RevolvedProfile(0.15*vm.Y3D, vm.X3D, z, c5, 0.15*vm.Y3D, vm.X3D)

# model = vm.VolumeModel([profile1, profile2, profile3, profile4, profile5])
model = vm.VolumeModel([profile5])
model.babylonjs()


