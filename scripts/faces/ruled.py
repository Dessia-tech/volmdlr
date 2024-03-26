#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import volmdlr
import volmdlr.edges
import volmdlr.primitives3d
import volmdlr.wires
from volmdlr import faces, surfaces

n = 5

# First line
radius1 = {}
current_point = volmdlr.Point3D(0,0,0)
points1 = []
for i in range(n):
    points1.append(current_point)
    delta = volmdlr.Point3D.random(-0.01, 0.01, 0.02, 0.06, -0.01, 0.01)
    r = 0.01*delta.norm()
    radius1[i] = r
    current_point += delta
del radius1[n-1]

wire1 = volmdlr.primitives3d.OpenRoundedLineSegments3D(points1, radius1, adapt_radius=True)

# First line
radius2 = {}
points2 = []
current_point = volmdlr.Point3D(0.1,0,0)
for i in range(n):
    points2.append(current_point)
    delta = volmdlr.Point3D.random(-0.01, 0.01, 0.02, 0.06, -0.01, 0.01)
    r = 0.01*delta.norm()
    radius2[i] = r
    current_point += delta
del radius2[n-1]

wire2 = volmdlr.primitives3d.OpenRoundedLineSegments3D(points2, radius2, adapt_radius=True)

ruled_surface = surfaces.RuledSurface3D(wire1, wire2)
# ruled_surface.babylonjs()
ax = wire1.plot()
wire2.plot(ax=ax)

face = volmdlr.faces.RuledFace3D.from_surface_rectangular_cut(ruled_surface, 0, 1,0, 1)
face.babylonjs()

circle1 = volmdlr.curves.Circle3D(volmdlr.OXYZ, 0.1)
circle2 = volmdlr.curves.Circle3D(volmdlr.Frame3D(0.1*volmdlr.Z3D,
                                                 volmdlr.X3D,
                                                 volmdlr.Y3D,
                                                 volmdlr.Z3D),
                                 0.12, volmdlr.Z3D)
ruled_surface = surfaces.RuledSurface3D(circle1, circle2)
face2 = faces.RuledFace3D.from_surface_rectangular_cut(ruled_surface, 0, 1,0, 1)
face2.babylonjs()
# model = volmdlr.model.VolumeModel([face]).babylonjs(debug=True)
