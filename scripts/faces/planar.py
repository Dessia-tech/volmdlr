
import volmdlr
import volmdlr.edges
import volmdlr.wires
import volmdlr.faces

p1 = volmdlr.Point3D(0.15, 0.48, 0.5)
p2 = volmdlr.Point3D(0.15, 0.1, 0.5)

p1s = volmdlr.Point2D(0, 0)
p2s = volmdlr.Point2D(0.1, 0)
p3s = volmdlr.Point2D(0.2, 0.1)
p4s = volmdlr.Point2D(-0.01, 0.05)
surface2d = volmdlr.faces.Surface2D(volmdlr.wires.ClosedPolygon2D([p1s, p2s, p3s, p4s]), [])

u = volmdlr.Vector3D(0.1, 0.7, -0.5)
u.normalize()
v = u.deterministic_unit_normal_vector()
w = u.cross(v)
plane = volmdlr.faces.Plane3D(frame=volmdlr.Frame3D(0.1*volmdlr.X3D, u, v, w))
face = volmdlr.faces.PlaneFace3D(plane, surface2d)

ax = face.plot()
p1.plot(ax=ax, color='b')
p2.plot(ax=ax, color='g')

l1 = volmdlr.edges.LineSegment3D(p1, p1+w)
l2 = volmdlr.edges.LineSegment3D(p2, p2+w)

l1.plot(ax=ax, color='b')
l2.plot(ax=ax, color='g')

i1 = face.linesegment_intersections(l1)
if i1:
    i1[0].plot(ax=ax, color='r')

i2 = face.linesegment_intersections(l2)
if i2:
    i2[0].plot(ax=ax, color='r')

plane_inter_1 = plane.linesegment_intersections(l1)
if plane_inter_1:
    plane_inter_1[0].plot(ax=ax, color='b')
plane_inter_2 = plane.linesegment_intersections(l2)
if plane_inter_2:
    plane_inter_2[0].plot(ax=ax, color='g')

plane_inter_1_2d = plane.point3d_to_2d(plane_inter_1[0])
plane_inter_2_2d = plane.point3d_to_2d(plane_inter_2[0])

ax2 = face.surface2d.plot()
plane_inter_1_2d.plot(ax=ax2, color='b')
plane_inter_2_2d.plot(ax=ax2, color='g')

assert surface2d.point_belongs(plane_inter_1_2d) == True
assert surface2d.point_belongs(plane_inter_2_2d) == False

p1_2dto3d = plane.point2d_to_3d(plane_inter_1_2d)
p1_2dto3d.plot(ax=ax, color='b')
assert p1_2dto3d == plane_inter_1[0]