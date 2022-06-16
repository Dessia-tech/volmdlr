#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 2022

@author: s.bendjebla
"""
import volmdlr
import volmdlr as vm
import volmdlr.edges as edges
import volmdlr.wires as wires
import volmdlr.primitives3d as primitives3d
import volmdlr.primitives2d as primitives2d
import gmsh

# %% Extrusion

points = [vm.Point2D(0, 0),
          vm.Point2D(0.1, 0.),
          vm.Point2D(0.1, 0.2),
          vm.Point2D(0.03, 0.15),
          vm.Point2D(0.,0.21)]

outer_profile = vm.wires.Contour2D.from_points(points)

profile=primitives3d.ExtrudedProfile(vm.O3D, vm.Y3D, vm.Z3D, outer_profile, [], vm.X3D*0.1, name = 'extrusion')

model=vm.core.VolumeModel([profile])
model.to_geo('model')
# model.babylonjs()

# %% .geo file lines

# from itertools import chain

# faces = profile.faces
# primitives = []
# points = set()
# for face in faces:
#     for c, contour in enumerate(list(chain(*[[face.outer_contour3d], face.inner_contours3d]))):
#         if isinstance(contour, volmdlr.wires.Circle2D):
#             points.add(volmdlr.Point3D(contour.radius, contour.center.y, 0))
#             points.add(volmdlr.Point3D(contour.center.x, contour.center.y, 0))
#             points.add(volmdlr.Point3D(-contour.radius, contour.center.y, 0))

#         else:
#             for p, primitive in enumerate(contour.primitives):
#                 if isinstance(primitive, volmdlr.edges.LineSegment):
#                     points.add(primitive.start)
#                     points.add(primitive.end)
#                 if ((primitive not in primitives)
#                     and (primitive.reverse() not in primitives)):
#                     primitives.append(primitive)

# indices_check = len(primitives)*[None]

# points = list(points)
# lines = []
# for p, point in enumerate(points):
#     lines.append(point.get_geo_lines(tag=p+1, mesh_size=1))

# line_account = 1
# line_surface = []
# lines_tags = []  # points_tags = []
# line_loop_account = 1

# for face in faces:
#     line_surface = []
#     for c, contour in enumerate(list(chain(*[[face.outer_contour3d], face.inner_contours3d]))):
#         lines_tags = []
#         if isinstance(contour, volmdlr.wires.Circle2D):
#             pass
#         else:
#             for p, primitive in enumerate(contour.primitives):
#                 if isinstance(primitive, volmdlr.edges.LineSegment):
#                     try:
#                         index = primitives.index(primitive)
#                     except ValueError:
#                         index = primitives.index(primitive.reverse())
#                     if indices_check[index]:
#                         lines_tags.append(indices_check[index])
#                     else:
#                         start_point_tag = points.index(primitive.start)+1
#                         end_point_tag = points.index(primitive.end)+1
#                         lines.append(primitive.get_geo_lines(tag=line_account,
#                                                              start_point_tag=start_point_tag,
#                                                              end_point_tag=end_point_tag))

#                         lines_tags.append(line_account)
#                         indices_check[index] = line_account
#                         line_account += 1

#             lines.append('Line Loop(' + str(line_loop_account) + ') = {' + str(lines_tags)[1:-1] + '};')
#             line_surface.append(line_loop_account)
#             line_loop_account += 1
#             lines_tags = []

#             lines.append('Plane Surface(' + str(1) + ') = {' + str(line_surface)[1:-1] + '};')
#     line_surface = []
    


# # lines = face.get_geo_lines()

# # with open('face_geo.geo', 'w') as f:
# #     for line in lines:
# #         f.write(line)
# #         f.write('\n')

# # # %% gmsh file generation

# # gmsh.initialize()
# # gmsh.open("face_geo.geo")

# # gmsh.model.geo.synchronize()
# # gmsh.model.mesh.generate(2)

# # gmsh.write("face_geo.msh")

# # gmsh.finalize()

