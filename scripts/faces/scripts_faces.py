import volmdlr as vm
import volmdlr.edges as edges
import matplotlib.pyplot as plt
import volmdlr.primitives3d as primitives3D
import volmdlr.wires as vmw
import volmdlr.faces as vmf
import volmdlr.edges as vme

import math
import volmdlr.step as vm_step

"""########################################### UNION 1 #################################### """
# resolution = 0.0010
#
# box_red = primitives3D.Block(
#     vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.4, 0, 0),
#                 vm.Vector3D(0, 0.4, 0), vm.Vector3D(0, 0, 0.4)),
#     color=(0.2, 1, 0.4), alpha=0.6)
#
# box_red.color = (1, 0.1, 0.1)
# box_red.name = 'box_red'


poly1_vol1 = vmw.ClosedPolygon3D([vm.Point3D(-0.1, -0.05, 0),
                                  vm.Point3D(-0.15, 0.1, 0),
                                  vm.Point3D(0.05, 0.2, 0),
                                  vm.Point3D(0.12, 0.15, 0),
                                  vm.Point3D(0.1, -0.02, 0)])

poly2_vol1 = poly1_vol1.rotation(vm.O3D, vm.Z3D, math.pi).translation(
    0.2 * vm.Z3D)
poly3_vol1 = poly2_vol1.rotation(vm.O3D, vm.Z3D, math.pi / 8).translation(
    0.1 * (vm.Z3D + vm.X3D + vm.Y3D))

point_triangles = poly1_vol1.sewing(poly2_vol1, vm.X3D,
                                    vm.Y3D) + poly2_vol1.sewing(poly3_vol1,
                                                                vm.X3D, vm.Y3D)
faces = [vmf.Triangle3D(trio[0], trio[1], trio[2]) for trio in point_triangles]

plane3d_1 = vmf.Plane3D.from_plane_vectors(vm.O3D, vm.X3D, vm.Y3D)
surf2d_1 = vmf.Surface2D(poly1_vol1.to_2d(vm.O3D, vm.X3D, vm.Y3D), [])

plane3d_2 = vmf.Plane3D.from_plane_vectors(0.3 * vm.Z3D, vm.X3D, vm.Y3D)
surf2d_2 = vmf.Surface2D(poly3_vol1.to_2d(vm.O3D, vm.X3D, vm.Y3D), [])
faces += [vmf.PlaneFace3D(plane3d_1, surf2d_1),
          vmf.PlaneFace3D(plane3d_2, surf2d_2)]

shell = vmf.ClosedShell3D(faces)
shell.color = (0.1, 1, 0.1)
shell.alpha = 0.4
ax = faces[0].outer_contour3d.plot()
for face in faces[1:]:
    face.outer_contour3d.plot(ax=ax)

point = vm.Point3D(0.6, -0.05, 0.15)
point.plot(ax=ax, color="r")
dist_point = shell.minimum_distance_point(point)
line = edges.LineSegment3D(point, dist_point)
line.plot(ax=ax, color="y")
print("minimal distance to shell from point is :",
      point.point_distance(dist_point))

shell2 = shell.rotation(vm.O3D, vm.Z3D, math.pi).translation(0.7 * vm.Z3D - 0.1 * vm.Y3D)
for face in shell2.faces:
    face.outer_contour3d.plot(ax=ax, color='b')
distance_to_shell = shell.distance_to_shell(shell2, resolution=10)
print("distance from shell to shell2 :", distance_to_shell)
