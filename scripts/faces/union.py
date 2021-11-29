#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:34:09 2021

@author: dasilva
"""
import volmdlr as vm
import volmdlr.wires as vmw
import volmdlr.faces as vmf

import math

number_points = 50

poly_1 = vmw.ClosedPolygon3D([vm.Point3D(-0.3, 0.05, -0.20),
                                vm.Point3D(0, 0.25, -0.20),
                                vm.Point3D(0.25, 0.1, -0.20),
                                vm.Point3D(0.2, -0.15, -0.20),
                                vm.Point3D(-0.2, -0.12, -0.20)])

length_poly_11 = poly_1.length()

points_poly_11 = [poly_1.point_at_abscissa(k*length_poly_11/(number_points)) for k in range(number_points)]

new_poly_11 = vmw.ClosedPolygon3D(points_poly_11)

new_poly_12 = new_poly_11.translation(0.3*vm.Z3D).rotation(vm.O3D, vm.Z3D, math.pi/2)

new_poly_13 = new_poly_12.translation(0.05*vm.Z3D)

new_poly_14 = new_poly_13.translation(0.2*vm.Z3D).rotation(vm.O3D, vm.Z3D, math.pi/4)

faces1 = new_poly_11.sewing(new_poly_12, vm.X3D, vm.Y3D) + new_poly_12.sewing(new_poly_13, vm.X3D, vm.Y3D) + new_poly_13.sewing(new_poly_14, vm.X3D, vm.Y3D)
# faces1 = [vmf.Triangle3D(trio[0], trio[1], trio[2]) for trio in points_triangles_1]

plane3d_1 = vmf.Plane3D.from_plane_vectors(-0.2*vm.Z3D, vm.X3D, vm.Y3D)
surf2d_1 = vmf.Surface2D(new_poly_11.to_2d(vm.O3D, vm.X3D, vm.Y3D),[])

plane3d_2 = vmf.Plane3D.from_plane_vectors(0.35*vm.Z3D, vm.X3D, vm.Y3D)
surf2d_2 = vmf.Surface2D(new_poly_14.to_2d(vm.O3D, vm.X3D, vm.Y3D),[])
faces1 += [vmf.PlaneFace3D(plane3d_1, surf2d_1), vmf.PlaneFace3D(plane3d_2, surf2d_2)]

shell1 = vmf.ClosedShell3D(faces1)

poly_2 = vmw.ClosedPolygon3D([vm.Point3D(-0.10, 0.05, 0),
                                vm.Point3D(-0.07, 0.05, 0.05),
                                vm.Point3D(0, 0.05, 0.10),
                                vm.Point3D(0.05, 0.05, 0.07),
                                vm.Point3D(0.10, 0.05, 0)])

length_poly_2 = poly_2.length()

points_poly_2 = [poly_2.point_at_abscissa(k*length_poly_2/(number_points)) for k in range(number_points)]

new_poly_21 = vmw.ClosedPolygon3D(points_poly_2)
new_poly_22 = new_poly_21.translation(0.1*vm.Y3D).rotation(vm.O3D, vm.Y3D, math.pi/2)
new_poly_23 = new_poly_22.translation(0.05*vm.Y3D)
new_poly_24 = new_poly_23.translation(0.2*vm.Y3D).rotation(vm.O3D, vm.Y3D, math.pi/4)
faces2 = new_poly_21.sewing(new_poly_22, vm.X3D, vm.Z3D) + new_poly_23.sewing(new_poly_22, vm.X3D, vm.Z3D) + new_poly_23.sewing(new_poly_24, vm.X3D, vm.Z3D)

# faces2 = [vmf.Triangle3D(trio[0], trio[1], trio[2]) for trio in points_triangles_2]

plane3d_3 = vmf.Plane3D.from_plane_vectors(0.05*vm.Y3D, vm.Z3D, vm.X3D)
surf2d_3 = vmf.Surface2D(new_poly_21.to_2d(vm.O3D, vm.Z3D, vm.X3D),[])

plane3d_4 = vmf.Plane3D.from_plane_vectors(0.4*vm.Y3D, vm.Z3D, vm.X3D)
surf2d_4 = vmf.Surface2D(new_poly_24.to_2d(vm.O3D, vm.Z3D, vm.X3D),[])
faces2 += [vmf.PlaneFace3D(plane3d_3, surf2d_3), vmf.PlaneFace3D(plane3d_4, surf2d_4)]

shell2 = vmf.ClosedShell3D(faces2)
new_box = shell1.union(shell2)
subtract_to_closed_shell = shell1.subtract_to_closed_shell(shell2)
# new_box = shell1.intersection(shell2)
for shell in [new_box, subtract_to_closed_shell]:
    shell[0].color = (1, 0.1, 0.1)
    shell[0].alpha = 0.6
    vm.core.VolumeModel(shell).babylonjs()


shell3 = shell2.rotation(vm.O3D, vm.Z3D, math.pi).translation(0.3*vm.Z3D-0.1*vm.Y3D)
# # new_box = shell1.union(shell3)
# # for shell in new_box:
# #     shell.color = (1, 0.1, 0.1)
# #     shell.alpha = 0.6
# # vm.core.VolumeModel(new_box).babylonjs()

new_box = new_box[0].union(shell3)
subtract_to_closed_shell = subtract_to_closed_shell[0].subtract_to_closed_shell(shell3)
# new_box = new_box[0].intersection(shell3)
for shell in [new_box, subtract_to_closed_shell]:
    shell[0].color = (1, 0.1, 0.1)
    shell[0].alpha = 0.6
    vm.core.VolumeModel(shell).babylonjs()
