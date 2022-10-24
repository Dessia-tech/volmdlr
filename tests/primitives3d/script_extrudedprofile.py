from volmdlr import primitives3d as p3d
import volmdlr as vm
import volmdlr.core
import volmdlr.wires as vmw
import volmdlr.edges as vme

points = [vm.Point2D(1, 0.5), vm.Point2D(-1, 0.5), vm.Point2D(-1, -0.5), vm.Point2D(1, -0.5)]
line_segments = [vme.LineSegment2D(points[0], points[1]), vme.LineSegment2D(points[1], points[2]),
                 vme.LineSegment2D(points[2], points[3]), vme.LineSegment2D(points[3], points[0])]

circle = vmw.Circle2D(vm.O2D, 0.2)
contour = vmw.Contour2D(line_segments)
polygon = vmw.ClosedPolygon2D(points, "square")

frame = volmdlr.Frame3D(vm.O3D, vm.X3D, vm.Y3D, vm.Z3D)
profiles = [circle, contour, polygon]

frame3 = volmdlr.Frame3D(vm.O3D, vm.X3D, vm.Z3D, vm.Y3D)
ext4 = p3d.ExtrudedProfile(frame3, polygon, [circle], 2)
ext4.babylonjs()

# for profil in profiles:
#     ext = p3d.ExtrudedProfile(frame, profil, [], length=2, color=(0.8, 0.2, 0.1), alpha=0.8, name="block")
#     ext.babylonjs()