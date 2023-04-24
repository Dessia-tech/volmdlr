#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 18 12:40:31 2020

@author: masfaraud
"""

import volmdlr
import volmdlr.faces
from volmdlr import faces, surfaces
R = 0.32

surface = surfaces.CylindricalSurface3D(volmdlr.OXYZ, R)

face = faces.CylindricalFace3D.from_surface_rectangular_cut(surface, -0.01, 1.3, -0.1, 0.3)
face.babylonjs(debug=True, use_cdn=False)

lines_x, lines_y = face.triangulation_lines()
ax = face.surface2d.plot()
for line in lines_x+lines_y:
    line.plot(ax=ax, color='r')

ax2 = face.surface2d.plot()
for surface in face.surface2d.split_by_lines(lines_x):
    surface.plot(ax=ax2, color='b')
