import volmdlr
import volmdlr.wires as vmw
import volmdlr.faces as vmf
import volmdlr.primitives3d as p3d

points_polygon1 = [volmdlr.Point3D(0.5, 0.5, 0), volmdlr.Point3D(0.5, -0.5, 0), volmdlr.Point3D(-0.5, -0.5, 0),
                   volmdlr.Point3D(-0.5, 0.5, 0)]

# points_polygon2 = [volmdlr.Point3D(0.25, 0.25, 0.5), volmdlr.Point3D(0.25, -0.25, 0.5), volmdlr.Point3D(0.0, -0.4, 0.5),
#                    volmdlr.Point3D(-0.25, -0.25, 0.5), volmdlr.Point3D(-0.25, 0.25, 0.5), volmdlr.Point3D(0.0, 0.4, 0.5)]

# points_polygon2 = [volmdlr.Point3D(0.25, 0.25, 0.5), volmdlr.Point3D(0.0, 0.4, 0.5), volmdlr.Point3D(-0.25, 0.25, 0.5),
#                    volmdlr.Point3D(-0.25, -0.25, 0.5), volmdlr.Point3D(0.0, -0.4, 0.5), volmdlr.Point3D(0.25, -0.25, 0.5)]

points_polygon2 = [volmdlr.Point3D(0.25, 0.15, 0.5), volmdlr.Point3D(0.25, -0.15, 0.5), volmdlr.Point3D(0.0, -0.3, 0.5),
                   volmdlr.Point3D(-0.25, -0.15, 0.5), volmdlr.Point3D(-0.25, 0.15, 0.5), volmdlr.Point3D(0.0, 0.3, 0.5)]



polygon1 = vmw.ClosedPolygon3D(points_polygon1)

polygon2 = vmw.ClosedPolygon3D(points_polygon2)

frame = volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0.5), volmdlr.Vector3D(1, 0, 0), volmdlr.Vector3D(1, 0, 0), volmdlr.Vector3D(1, 0, 0))
circle1 = vmw.Circle3D(frame, 0.3)


lists_profiles = [polygon1, polygon2]

loft = p3d.Loft(lists_profiles, color=(1, 0.5, 0.2))
loft.babylonjs()

# list_triangles_points = polygon1.sewing(polygon2, volmdlr.X3D, volmdlr.Y3D)
# faces = [vmf.Triangle3D(*triangle_points, alpha=0.9,
#                         color=(1, 0.1, 0.1))
#          for triangle_points in list_triangles_points]
# volmdlr.core.VolumeModel(faces).babylonjs()
#
# loft2 = p3d.Loft([polygon1, circle], color=(1, 0.5, 0.2))
# loft2.babylonjs()

# ruled_surface = vmf.RuledSurface3D(polygon1, polygon2, name='Test')
# ruled_face = ruled_surface.rectangular_cut(0, 1, 0, 1, name=ruled_surface.name)
# ruled_face.babylonjs()

