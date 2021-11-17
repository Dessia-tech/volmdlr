import volmdlr as vm
import volmdlr.primitives3d as primitives3D
import volmdlr.faces as vmf
import math
import volmdlr.wires as vmw

poly1_vol1 = vmw.ClosedPolygon3D([vm.Point3D(-0.1, -0.05, 0),
                                  vm.Point3D(-0.15, 0.1, 0),
                                  vm.Point3D(0.05, 0.2, 0),
                                  vm.Point3D(0.12, 0.15, 0),
                                  vm.Point3D(0.1, -0.02, 0)])

poly2_vol1 = poly1_vol1.rotation(vm.O3D, vm.Z3D, math.pi).translation(0.2*vm.Z3D)
poly3_vol1 = poly2_vol1.rotation(vm.O3D, vm.Z3D, math.pi/8).translation(0.1*(vm.Z3D+vm.X3D+vm.Y3D))

point_triangles = poly1_vol1.sewing(poly2_vol1, vm.X3D, vm.Y3D) + poly2_vol1.sewing(poly3_vol1, vm.X3D, vm.Y3D)
faces = [vmf.Triangle3D(trio[0], trio[1], trio[2]) for trio in point_triangles]

plane3d_1 = vmf.Plane3D.from_plane_vectors(vm.O3D, vm.X3D, vm.Y3D)
surf2d_1 = vmf.Surface2D(poly1_vol1.to_2d(vm.O3D, vm.X3D, vm.Y3D),[])

plane3d_2 = vmf.Plane3D.from_plane_vectors(0.3*vm.Z3D, vm.X3D, vm.Y3D)
surf2d_2 = vmf.Surface2D(poly3_vol1.to_2d(vm.O3D, vm.X3D, vm.Y3D),[])
faces += [vmf.PlaneFace3D(plane3d_1, surf2d_1), vmf.PlaneFace3D(plane3d_2, surf2d_2)]

shell1 = vmf.ClosedShell3D(faces)
shell1.color = (0.1, 1, 0.1)
shell1.alpha = 0.4

poly1_vol2 = vmw.ClosedPolygon3D([vm.Point3D(-0.1, -0.1, -0.2),
                                  vm.Point3D(-0.15, -0.1, -0.05),
                                  vm.Point3D(0.05, -0.1, 0.2),
                                  vm.Point3D(0.12, -0.1, 0.05),
                                  vm.Point3D(0.1, -0.1, -0.02)])


poly2_vol2 = poly1_vol2.rotation(vm.O3D, vm.Y3D, math.pi/2).translation(0.02*vm.Y3D)
poly3_vol2 = poly2_vol2.rotation(vm.O3D, vm.Y3D, math.pi/8).translation(0.1*(vm.Z3D+vm.X3D+vm.Y3D))
poly4_vol2 = poly3_vol2.rotation(vm.O3D, vm.Y3D, math.pi/4).translation(0.05*vm.Y3D)
poly5_vol2 = poly4_vol2.rotation(vm.O3D, vm.Y3D, math.pi/10).translation(0.2*vm.Y3D)

point_triangles_2 = poly1_vol2.sewing(poly2_vol2, vm.X3D, vm.Z3D) + poly2_vol2.sewing(poly3_vol2, vm.X3D, vm.Z3D) +\
                    poly3_vol2.sewing(poly4_vol2, vm.X3D, vm.Z3D) + poly4_vol2.sewing(poly5_vol2, vm.X3D, vm.Z3D)

faces_2 = [vmf.Triangle3D(trio[0], trio[1], trio[2]) for trio in point_triangles_2]

plane3d_3 = vmf.Plane3D.from_plane_vectors(-0.1*vm.Y3D, vm.X3D, vm.Z3D)
surf2d_3 = vmf.Surface2D(poly1_vol2.to_2d(vm.O3D, vm.X3D, vm.Z3D),[])

plane3d_4 = vmf.Plane3D.from_plane_vectors(0.27*vm.Y3D, vm.X3D, vm.Z3D)
surf2d_4 = vmf.Surface2D(poly5_vol2.to_2d(vm.O3D, vm.X3D, vm.Z3D),[])
faces_2 += [vmf.PlaneFace3D(plane3d_3, surf2d_3), vmf.PlaneFace3D(plane3d_4, surf2d_4)]


shell2 = vmf.ClosedShell3D(faces_2)
union_box = shell1.union(shell2)
subtraction_box = shell1.subtract(shell2)
intersection_box = shell1.intersection(shell2)

for new_box in [union_box, subtraction_box, intersection_box]:
    for shell in new_box:
        shell.color = (1, 0.1, 0.1)
        shell.alpha = 0.6
    vm.core.VolumeModel(new_box + [shell1, shell2]).babylonjs()