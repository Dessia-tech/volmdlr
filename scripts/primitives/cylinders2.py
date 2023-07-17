"""
Test script to define a cylinder based on a ClosedShell instead of a RevolvedProfile
"""
import math

import volmdlr
import volmdlr.primitives3d
import volmdlr.edges
import volmdlr.faces

POSITION = volmdlr.O3D
AXIS = volmdlr.X3D
RADIUS = 0.1
LENGTH = 0.5
frame = volmdlr.Frame3D.from_point_and_normal(POSITION, AXIS)
cylinder = volmdlr.primitives3d.Cylinder(frame, RADIUS, LENGTH, alpha=0.1)
cylinder.babylonjs()


def shell_faces(cyl: volmdlr.primitives3d.Cylinder):
    normal_vector = cyl.axis.random_unit_normal_vector()
    p1 = cyl.position + 1/2 * cyl.length * cyl.axis + cyl.radius * normal_vector
    p2 = cyl.position - 1/2 * cyl.length * cyl.axis + cyl.radius * normal_vector
    ls = volmdlr.edges.LineSegment3D(start=p1, end=p2)

    return ls.revolution(axis_point=cyl.position, axis=cyl.axis, angle=volmdlr.TWO_PI)


cs = volmdlr.shells.ClosedShell3D(faces=shell_faces(cylinder))
cs.babylonjs()

# cs2 = cs.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi/2)
#
# cs2.babylonjs()
#
# cs3 = cs.intersection(cs2)[0]
# cs3.babylonjs()
