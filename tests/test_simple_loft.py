import volmdlr
import volmdlr.wires as vmw
import volmdlr.primitives3d as p3d

points_polygon1 = [volmdlr.Point3D(0.5, 0.5, 0), volmdlr.Point3D(0.5, -0.5, 0), volmdlr.Point3D(-0.5, -0.5, 0),
                   volmdlr.Point3D(-0.5, 0.5, 0)]

points_polygon2 = [volmdlr.Point3D(0.25, 0.25, 0.5), volmdlr.Point3D(0.25, -0.25, 0.5), volmdlr.Point3D(0.0, -0.4, 0.5),
                   volmdlr.Point3D(-0.25, -0.25, 0.5), volmdlr.Point3D(-0.25, 0.25, 0.5), volmdlr.Point3D(0.0, 0.4, 0.5)
                   ]

polygon1 = vmw.ClosedPolygon3D(points_polygon1)

polygon2 = vmw.ClosedPolygon3D(points_polygon2)

lists_profiles = [polygon1, polygon2]

loft = p3d.Loft(lists_profiles, color=(1, 0.5, 0.2))
loft.babylonjs()