#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:34:09 2021

@author: dasilva
"""
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
resolution = 0.0010

box = primitives3D.Block(
    vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.4, 0, 0),
                vm.Vector3D(0, 0.1, 0), vm.Vector3D(0, 0, 0.8)),
    alpha=0.6)

box_red = primitives3D.Block(
    vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(0.4, 0, 0),
                vm.Vector3D(0, 0.4, 0), vm.Vector3D(0, 0, 0.4)),
    color=(0.2, 1, 0.4), alpha=0.6)

box_red.color = (1, 0.1, 0.1)
box_red.name = 'box_red'

box_green = box.frame_mapping(vm.Frame3D(vm.Point3D(0, 0.8, 0), vm.Vector3D(1, 0, 0),
                          vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'new', copy=True)

box_green.color = (0.1, 1, 0.1)
box_green.name = 'box_green'


box_blue = box.frame_mapping(vm.Frame3D(vm.Point3D(0, 0, 0), vm.Vector3D(1, 0, 0),
                          vm.Vector3D(0, 1, 0), vm.Vector3D(0, 0, 1)), 'old', copy=True)
box_blue.color = (0.1, 0.1, 1)
box_blue.name = 'box_blue'


new_box =box_red.union(box_blue)
# new_box =box_red.subtract(box_blue)
for shell in new_box:
    shell.color = (1, 0.1, 0.1)
    shell.alpha = 0.6
vm.core.VolumeModel(new_box).babylonjs()

# box_blue = box_blue.translation(vm.Point3D(0,0,0.3))
# vm.core.VolumeModel([box_blue, box_red]).babylonjs()
# new_box = box_red.union(box_blue)
# new_box =box_red.subtract(box_blue)
# for shell in new_box:
#     shell.color = (1, 0.1, 0.1)
#     shell.alpha = 0.6
# vm.core.VolumeModel(new_box).babylonjs()


"""########################################### UNION 2 #################################### """





# poly1_vol1 = vmw.ClosedPolygon3D([vm.Point3D(-0.1, -0.05, 0),
#                                   vm.Point3D(-0.15, 0.1, 0),
#                                   vm.Point3D(0.05, 0.2, 0),
#                                   vm.Point3D(0.12, 0.15, 0),
#                                   vm.Point3D(0.1, -0.02, 0)])

# poly2_vol1 = poly1_vol1.rotation(vm.O3D, vm.Z3D, math.pi).translation(0.2*vm.Z3D)
# poly3_vol1 = poly2_vol1.rotation(vm.O3D, vm.Z3D, math.pi/8).translation(0.1*(vm.Z3D+vm.X3D+vm.Y3D))

# point_triangles = poly1_vol1.sewing(poly2_vol1, vm.X3D, vm.Y3D) + poly2_vol1.sewing(poly3_vol1, vm.X3D, vm.Y3D)
# faces = [vmf.Triangle3D(trio[0], trio[1], trio[2]) for trio in point_triangles]

# plane3d_1 = vmf.Plane3D.from_plane_vectors(vm.O3D, vm.X3D, vm.Y3D)
# surf2d_1 = vmf.Surface2D(poly1_vol1.to_2d(vm.O3D, vm.X3D, vm.Y3D),[])

# plane3d_2 = vmf.Plane3D.from_plane_vectors(0.3*vm.Z3D, vm.X3D, vm.Y3D)
# surf2d_2 = vmf.Surface2D(poly3_vol1.to_2d(vm.O3D, vm.X3D, vm.Y3D),[])
# faces += [vmf.PlaneFace3D(plane3d_1, surf2d_1), vmf.PlaneFace3D(plane3d_2, surf2d_2)]

# shell1 = vmf.ClosedShell3D(faces)
# shell1.color = (0.1, 1, 0.1)
# shell1.alpha = 0.4

# poly1_vol2 = vmw.ClosedPolygon3D([vm.Point3D(-0.1, -0.1, -0.2),
#                                   vm.Point3D(-0.15, -0.1, -0.05),
#                                   vm.Point3D(0.05, -0.1, 0.2),
#                                   vm.Point3D(0.12, -0.1, 0.05),
#                                   vm.Point3D(0.1, -0.1, -0.02)])


# poly2_vol2 = poly1_vol2.rotation(vm.O3D, vm.Y3D, math.pi/2).translation(0.02*vm.Y3D)
# poly3_vol2 = poly2_vol2.rotation(vm.O3D, vm.Y3D, math.pi/8).translation(0.1*(vm.Z3D+vm.X3D+vm.Y3D))
# poly4_vol2 = poly3_vol2.rotation(vm.O3D, vm.Y3D, math.pi/4).translation(0.05*vm.Y3D)
# poly5_vol2 = poly4_vol2.rotation(vm.O3D, vm.Y3D, math.pi/10).translation(0.2*vm.Y3D)

# point_triangles_2 = poly1_vol2.sewing(poly2_vol2, vm.X3D, vm.Z3D) + poly2_vol2.sewing(poly3_vol2, vm.X3D, vm.Z3D) +\
#                     poly3_vol2.sewing(poly4_vol2, vm.X3D, vm.Z3D) + poly4_vol2.sewing(poly5_vol2, vm.X3D, vm.Z3D)

# faces_2 = [vmf.Triangle3D(trio[0], trio[1], trio[2]) for trio in point_triangles_2]

# plane3d_3 = vmf.Plane3D.from_plane_vectors(-0.1*vm.Y3D, vm.X3D, vm.Z3D)
# surf2d_3 = vmf.Surface2D(poly1_vol2.to_2d(vm.O3D, vm.X3D, vm.Z3D),[])

# plane3d_4 = vmf.Plane3D.from_plane_vectors(0.27*vm.Y3D, vm.X3D, vm.Z3D)
# surf2d_4 = vmf.Surface2D(poly5_vol2.to_2d(vm.O3D, vm.X3D, vm.Z3D),[])
# faces_2 += [vmf.PlaneFace3D(plane3d_3, surf2d_3), vmf.PlaneFace3D(plane3d_4, surf2d_4)]


# shell2 = vmf.ClosedShell3D(faces_2)
# # new_box = shell1.union(shell2)
# new_box = shell1.subtract(shell2)
# for shell in new_box:
#     shell.color = (1, 0.1, 0.1)
#     shell.alpha = 0.6
# vm.core.VolumeModel(new_box).babylonjs()

# """# =============================================================================
# # USECASE 3
# # =============================================================================

# # Union between shell1 and shell2
# # Union between shell1 and shell3
# # Union between shell2 and shell3
# # Union between shell1, shell2 and shell3
# """

# number_points = 50
#
# poly_1 = vmw.ClosedPolygon3D([vm.Point3D(-0.3, 0.05, -0.20),
#                                 vm.Point3D(0, 0.25, -0.20),
#                                 vm.Point3D(0.25, 0.1, -0.20),
#                                 vm.Point3D(0.2, -0.15, -0.20),
#                                 vm.Point3D(-0.2, -0.12, -0.20)])
#
# length_poly_11 = poly_1.length()
#
# points_poly_11 = [poly_1.point_at_abscissa(k*length_poly_11/(number_points)) for k in range(number_points)]
#
# new_poly_11 = vmw.ClosedPolygon3D(points_poly_11)
#
# new_poly_12 = new_poly_11.translation(0.3*vm.Z3D).rotation(vm.O3D, vm.Z3D, math.pi/2)
#
# new_poly_13 = new_poly_12.translation(0.05*vm.Z3D)
#
# new_poly_14 = new_poly_13.translation(0.2*vm.Z3D).rotation(vm.O3D, vm.Z3D, math.pi/4)
#
# points_triangles_1 = new_poly_11.sewing(new_poly_12, vm.X3D, vm.Y3D) + new_poly_12.sewing(new_poly_13, vm.X3D, vm.Y3D) + new_poly_13.sewing(new_poly_14, vm.X3D, vm.Y3D)
# faces1 = [vmf.Triangle3D(trio[0], trio[1], trio[2]) for trio in points_triangles_1]
#
# plane3d_1 = vmf.Plane3D.from_plane_vectors(-0.2*vm.Z3D, vm.X3D, vm.Y3D)
# surf2d_1 = vmf.Surface2D(new_poly_11.to_2d(vm.O3D, vm.X3D, vm.Y3D),[])
#
# plane3d_2 = vmf.Plane3D.from_plane_vectors(0.35*vm.Z3D, vm.X3D, vm.Y3D)
# surf2d_2 = vmf.Surface2D(new_poly_14.to_2d(vm.O3D, vm.X3D, vm.Y3D),[])
# faces1 += [vmf.PlaneFace3D(plane3d_1, surf2d_1), vmf.PlaneFace3D(plane3d_2, surf2d_2)]
#
# shell1 = vmf.ClosedShell3D(faces1)
#
# poly_2 = vmw.ClosedPolygon3D([vm.Point3D(-0.10, 0.05, 0),
#                                 vm.Point3D(-0.07, 0.05, 0.05),
#                                 vm.Point3D(0, 0.05, 0.10),
#                                 vm.Point3D(0.05, 0.05, 0.07),
#                                 vm.Point3D(0.10, 0.05, 0)])
#
# length_poly_2 = poly_2.length()
#
# points_poly_2 = [poly_2.point_at_abscissa(k*length_poly_2/(number_points)) for k in range(number_points)]
#
# new_poly_21 = vmw.ClosedPolygon3D(points_poly_2)
# new_poly_22 = new_poly_21.translation(0.1*vm.Y3D).rotation(vm.O3D, vm.Y3D, math.pi/2)
# new_poly_23 = new_poly_22.translation(0.05*vm.Y3D)
# new_poly_24 = new_poly_23.translation(0.2*vm.Y3D).rotation(vm.O3D, vm.Y3D, math.pi/4)
# points_triangles_2 = new_poly_21.sewing(new_poly_22, vm.X3D, vm.Z3D) + new_poly_23.sewing(new_poly_22, vm.X3D, vm.Z3D) + new_poly_23.sewing(new_poly_24, vm.X3D, vm.Z3D)
#
# faces2 = [vmf.Triangle3D(trio[0], trio[1], trio[2]) for trio in points_triangles_2]
#
# plane3d_3 = vmf.Plane3D.from_plane_vectors(0.05*vm.Y3D, vm.Z3D, vm.X3D)
# surf2d_3 = vmf.Surface2D(new_poly_21.to_2d(vm.O3D, vm.Z3D, vm.X3D),[])
#
# plane3d_4 = vmf.Plane3D.from_plane_vectors(0.4*vm.Y3D, vm.Z3D, vm.X3D)
# surf2d_4 = vmf.Surface2D(new_poly_24.to_2d(vm.O3D, vm.Z3D, vm.X3D),[])
# faces2 += [vmf.PlaneFace3D(plane3d_3, surf2d_3), vmf.PlaneFace3D(plane3d_4, surf2d_4)]
#
# shell2 = vmf.ClosedShell3D(faces2)
# new_box = shell1.union(shell2)
# # for shell in new_box:
# #     shell.color = (1, 0.1, 0.1)
# #     shell.alpha = 0.6
# # vm.core.VolumeModel(new_box).babylonjs()
#
#
# shell3 = shell2.rotation(vm.O3D, vm.Z3D, math.pi).translation(0.3*vm.Z3D-0.1*vm.Y3D)
# # # new_box = shell1.union(shell3)
# # # for shell in new_box:
# # #     shell.color = (1, 0.1, 0.1)
# # #     shell.alpha = 0.6
# # # vm.core.VolumeModel(new_box).babylonjs()
#
# new_box = new_box[0].union(shell3)
# for shell in new_box:
#     shell.color = (1, 0.1, 0.1)
#     shell.alpha = 0.6
# vm.core.VolumeModel(new_box).babylonjs()