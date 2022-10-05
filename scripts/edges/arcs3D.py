#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 23 18:06:45 2018

@author: steven
"""

import volmdlr
import volmdlr.edges




i = volmdlr.X3D
e = i.rotation(volmdlr.O3D, volmdlr.Z3D, 1)
s = i.rotation(volmdlr.O3D, volmdlr.Z3D, -3.5)

a = volmdlr.edges.Arc3D(s, i, e)
assert a.angle == 4.5


# Random arc
i = volmdlr.Point3D.random(-1,1,-1,1,-1,1)
e = volmdlr.Point3D.random(-1,1,-1,1,-1,1)
s = volmdlr.Point3D.random(-1,1,-1,1,-1,1)

a = volmdlr.edges.Arc3D(s, i, e)
ax = a.plot()

# for p in a.polygon_points():
#     p.plot(ax=ax)
    
s.plot(ax=ax, color='r')
e.plot(ax=ax, color='g')
i.plot(ax=ax, color='b')


arc1 = volmdlr.edges.Arc3D(volmdlr.Point3D(-0.03096, 0.001162, -0.02),
                volmdlr.Point3D(-0.03120, -0.000400635, -0.02),
                volmdlr.Point3D(-0.026119083, 0.0, -0.02),
                volmdlr.Vector3D(0.0, 0.0, 0.001))

ax = arc1.plot(ax=ax)
# for p in arc1.polygon_points():
#     p.plot(ax=ax)


arc1.start.plot(ax=ax, color='r')
arc1.end.plot(ax=ax, color='g')
arc1.interior.plot(ax=ax, color='b')
arc1.center.plot(ax=ax, color='m')


print(arc1.center)
print(arc1.center-volmdlr.Point3D(-0.030962035803739997, 0.0011626900994054661, -0.02))
print(arc1.center-volmdlr.Point3D(-0.031209642286239472, -0.00040063570451895954, -0.02))
print(arc1.center-volmdlr.Point3D(-0.026119083, 0.0, -0.02))