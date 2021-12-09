#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  1 14:53:49 2019

@author: dumouchel
"""
from dessia_common import DessiaObject
import volmdlr as vm
import volmdlr.primitives3d as primitives3d
import volmdlr.stl
import volmdlr.cloud as vmc
# import mechanical_components.optimization.wires_protected as wires_opt
from dectree import DecisionTree
# import layout_3d.connector as conn
import plot_data

from random import random
import math
import os
from typing import Union, List, Tuple, Dict, BinaryIO
import numpy as np
import time
import datetime
import networkx as nx

from dessia_common.vectored_objects import Catalog, Objective, ParetoSettings, \
    ObjectiveSettings
from dessia_common import InstanceOf
from itertools import combinations, product, permutations
import inspect

class Boxes(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = ['volume_model']
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _eq_is_data_eq = True

    def _proximity(self, box1, box2):
        for face1, face2 in product(box1.faces, box2.faces):
            if face1 == face2:
                return True
        bb1 = box1.bounding_box
        bb2 = box2.bounding_box
        if bb1.distance_to_bbox(bb2) == 0:
            return True
        return False

    def _connected_boxes(self, boxes: List[vm.primitives3d.Block]):
        graph = nx.Graph()
        add_edges_from = []
        for box1, box2 in product(boxes, repeat=2):
            if math.isclose((box1.bounding_box + box2.bounding_box).volume(), box1.bounding_box.volume() + box2.bounding_box.volume()):
                add_edges_from.extend([(box1, box2)])
        graph.add_edges_from(add_edges_from)
        closed_shells = []
        for g in nx.connected_components(graph):
            elem_boxes = list(g)
            closed_shells.append(elem_boxes)
        return closed_shells

    def _merge_boxes(self, boxes: List[vm.primitives3d.Block]):
        graph = nx.Graph()
        for b1, b2 in product(boxes, repeat=2):
            if b1 != b2:
                if self._proximity(b1, b2):
                    graph.add_edges_from([(b1, b2)])
        node0 = list(graph.nodes())[0]

        dfs_tree = nx.dfs_tree(graph, node0)
        list_dfs_tree = list(dfs_tree.edges())
        boxe_cumul = list_dfs_tree[0][0]
        for i, (node1, node2) in enumerate(list_dfs_tree):
            output = boxe_cumul.union(node2, tol=1e-6)
            boxe_cumul = output[0]
            print(i)
        return boxe_cumul

    def define_block(self, adresse: Tuple[int, int, int]):
        origin = vm.Point3D(self.xmin + adresse[0] * self.cx + self.cx / 2,
                            self.ymin + adresse[1] * self.cy + self.cy / 2,
                            self.zmin + adresse[2] * self.cz + self.cz / 2)
        frame = vm.Frame3D(origin, self.cx * vm.X3D, self.cy * vm.Y3D,
                           self.cz * vm.Z3D)
        return vm.primitives3d.Block(frame)

    def define_adresses(self, box: vm.primitives3d.Block):
        bb = box.bounding_box
        nx = self.nx - int(round((self.xmax - bb.xmax) / self.cx, 0))
        ny = self.ny - int(round((self.ymax - bb.ymax) / self.cy, 0))
        nz = self.nz - int(round((self.zmax - bb.zmax) / self.cz, 0))
        return (nx - 1, ny - 1, nz - 1)

    def define_adresses_with_frame(self, frame: vm.Frame3D):
        xmax = frame.origin[0] + frame.u.norm() / 2
        ymax = frame.origin[1] + frame.v.norm() / 2
        zmax = frame.origin[2] + frame.w.norm() / 2
        nx = self.nx - int(round((self.xmax - xmax) / self.cx, 0))
        ny = self.ny - int(round((self.ymax - ymax) / self.cy, 0))
        nz = self.nz - int(round((self.zmax - zmax) / self.cz, 0))
        return (nx - 1, ny - 1, nz - 1)

    def _size_boxes(self):
        nb1 = math.ceil((self.xmax - self.xmin) / self.size_boxe)
        nb2 = math.ceil((self.ymax - self.ymin) / self.size_boxe)
        nb3 = math.ceil((self.zmax - self.zmin) / self.size_boxe)
        round_value = 1e6
        c1 = math.ceil((self.xmax - self.xmin) / (nb1) * round_value) / (round_value)
        c2 = math.ceil((self.ymax - self.ymin) / (nb2) * round_value) / (round_value)
        c3 = math.ceil((self.zmax - self.zmin) / (nb3) * round_value) / (round_value)
        return nb1, nb2, nb3, c1, c2, c3

    def _search_all_bounding_box(self, bounding_box):
        bb = bounding_box
        nx1 = int((bb.xmin - self.xmin) / self.cx)
        nx2 = self.nx - int((self.xmax - bb.xmax) / self.cx)
        ny1 = int((bb.ymin - self.ymin) / self.cy)
        ny2 = self.ny - int((self.ymax - bb.ymax) / self.cy)
        nz1 = int((bb.zmin - self.zmin) / self.cz)
        nz2 = self.nz - int((self.zmax - bb.zmax) / self.cz)

        adresses = []
        for i, j, k in product(range(nx1, nx2), range(ny1, ny2), range(nz1, nz2)):
            adresses.append((i, j, k))
        return adresses

    def smart_box_discretization(self):
        adresses = {}
        for j, primitive in enumerate(self.volume_model.primitives):
            for i, face in enumerate(primitive.faces):
                adds = self._search_all_bounding_box(face.bounding_box)
                for add in adds:
                    if add in adresses:
                        adresses[add].append((j, i))
                    else:
                        adresses[add] = [(j, i)]

        return adresses

class GiftWrap(Boxes, DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = ['volume_model']
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _eq_is_data_eq = True

    def __init__(self, volume_model: vm.core.VolumeModel = None,
                 name: str = ''):
        DessiaObject.__init__(self, name=name)
        self.volume_model = volume_model

    def contour_intersection(self, plane: vm.faces.PlaneFace3D):
        primitives = []
        for face in self.volume_model.primitives[0].faces:
            sol = plane.face_intersections(face)
            if sol:
                if isinstance(sol[0], volmdlr.wires.Wire3D):
                    primitives.extend([p for s in sol for p in s.primitives])
        if primitives == []:
            return None
        return primitives

    def _analysis_stl(self):
        graph = nx.Graph()
        add_edges_from = []
        for face in self.volume_model.primitives[0].faces:
            add_edges_from.extend([(face.point1, face.point2),
                                   (face.point1, face.point3),
                                   (face.point3, face.point2)])
        graph.add_edges_from(add_edges_from)
        closed_shells = []
        for g in nx.connected_components(graph):
            nodes = list(g)
            faces = []
            for face in self.volume_model.primitives[0].faces:
                if face.point1 not in nodes:
                    continue
                if face.point2 not in nodes:
                    continue
                if face.point3 not in nodes:
                    continue
                faces.append(face)
            if faces:
                closed_shells.append(vm.faces.ClosedShell3D(faces))
        return closed_shells

    def _draw_2d(self, primitives: List[vm.edges.LineSegment3D],
                 plane: vm.faces.PlaneFace3D):
        normal = plane.surface3d.frame.w
        vect1 = plane.surface3d.frame.u
        vect2 = plane.surface3d.frame.v
        primitive2d = []
        for p in primitives:
            point1 = p.start.dot(vect1) * vm.Vector2D(1, 0) + p.start.dot(vect2) * vm.Vector2D(0, 1)
            point2 = p.end.dot(vect1) * vm.Vector2D(1, 0) + p.end.dot(vect2) * vm.Vector2D(0, 1)
            if point1 != point2:
                primitive2d.append(vm.edges.LineSegment2D(point1, point2))
        contour = vm.wires.Contour2D(primitive2d)
        contour.plot()

    def _draw_3d(self, plane: vm.faces.PlaneFace3D):
        vol = vm.core.VolumeModel(self.volume_model.primitives + [plane])
        vol.babylonjs()

    def _contour_order(self, primitives: List[vm.edges.LineSegment3D],
                       plane: vm.faces.PlaneFace3D):
        graph = nx.Graph()
        for primitive in primitives:
            graph.add_edges_from([(primitive.start,
                                   primitive.end)])
        areas = []
        polys = []
        for g in nx.cycle_basis(graph):
            poly = vm.wires.ClosedPolygon3D(g)
            poly2d = poly.to_2d(g[0], plane.surface3d.frame.u, plane.surface3d.frame.v)
            areas.append(poly2d.area())
            polys.append(poly)
        print('areas', areas)
        areas_sort = np.argsort(np.array(areas))
        if polys != [] and areas[areas_sort[-1]] > 1e-4:
            return polys[areas_sort[-1]]
        else:
            return None

    def _update_orientation(self, polygon: vm.wires.ClosedPolygon3D,
                            plane: vm.faces.PlaneFace3D):
        center = polygon.average_center_point()
        point1, point2 = polygon.points[0], polygon.points[int(len(polygon.points) / 4)]
        vect1 = point1 - center
        vect1.normalize()
        vect2 = point2 - point1
        vect2.normalize()
        vect3 = vect1.cross(vect2)
        if plane.surface3d.frame.w.dot(vect3) > 0:
            return vm.wires.ClosedPolygon3D(polygon.points[::-1])
        return polygon

    def _discretize_polygon(self, polygon: vm.wires.ClosedPolygon3D):
        points = polygon.points[-1:] + polygon.points + polygon.points[0:1]
        new_points = []
        for point1, point2, point3 in zip(points[0:-2], points[1:-1], points[2:]):
            l1, l2, l3 = (point2 - point1).norm(), (point2 - point3).norm(), (point3 - point1).norm()
            cos = (l1 ** 2 + l2 ** 2 - l3 ** 2) / (2 * l1 * l2)
            if abs(cos) < 0.9:
                new_points.append(point2)
        return new_points

    def _graph(self, points1: List[vm.Point3D], points2: List[vm.Point3D]):
        graph = nx.Graph()
        previous_connection_points = []
        for pt1 in points1:
            all_dist = []
            for pt2 in points2:
                all_dist.append(pt1.point_distance(pt2))
            sort_all_dist = np.argsort(np.array(all_dist))

            next_available_indice = list(range(len(points2)))
            if previous_connection_points:
                all_indice = [points2.index(e) for e in previous_connection_points]
                print('all_indice', all_indice)
                previous_indice = points2.index(previous_connection_points[-1])
                next_available_indice = []
                if all_indice[-1] > all_indice[0] or len(all_indice) == 1:
                    next_available_indice += list(range(previous_indice, len(points2)))
                print(1, next_available_indice)
                if len(all_indice) > 1:
                    if all_indice[-1] > all_indice[0]:
                        next_available_indice += list(range(all_indice[0]))
                    else:
                        next_available_indice += list(range(all_indice[-1], all_indice[0]))
                if len(all_indice) == 1:
                    next_available_indice += list(range(all_indice[0]))
                print(2, next_available_indice)
                for i in all_indice[0:-1]:
                    if i in next_available_indice:
                        next_available_indice.remove(i)
            print('next_available_indice', next_available_indice)
            if not next_available_indice:
                break

            check = False
            for indice in sort_all_dist:
                if indice not in next_available_indice:
                    continue
                if points2[indice] not in previous_connection_points:
                    graph.add_edges_from([(pt1, points2[indice])])
                    check = True
                    break
                if previous_connection_points:
                    if points2[indice] == previous_connection_points[-1]:
                        graph.add_edges_from([(pt1, points2[indice])])
                        check = True
                        break
            if check and points2[indice] not in previous_connection_points:
                previous_connection_points.append(points2[indice])

        for p1, p2 in zip(points1, points1[1:] + points1[0:1]):
            graph.add_edges_from([(p1, p2)])

        dict_pos = nx.kamada_kawai_layout(graph)

        points_without_connection = []
        for pt2 in points2:
            if pt2 in graph.nodes():
                continue
            points_without_connection.append(pt2)
        print([points2.index(i) for i in points_without_connection])

        graph_points2 = nx.Graph()
        for p1, p2 in zip(points2, points2[1:] + points2[0:1]):
            graph_points2.add_edges_from([(p1, p2)])

        group_without_connection = []
        check_graph = nx.Graph()
        if points_without_connection:
            last_node = points_without_connection[0]
            check_graph.add_node(points_without_connection[0])
            if len(points_without_connection) > 1:
                for i in points_without_connection[1:]:
                    if i in list(graph_points2[last_node]):
                        check_graph.add_edges_from([(last_node, i)])
                    else:
                        check_graph.add_node(i)
                    last_node = i
                    # group_without_connection.append(last_group)
                    # last_group = [i]
            if points2[0] in check_graph.nodes() and \
                    points2[-1] in check_graph.nodes():
                check_graph.add_edges_from([(points2[0], points2[-1])])
            for g in nx.connected_components(check_graph):
                all_indice = list(g)
                print([points2.index(i) for i in all_indice])
                if len(all_indice) > 1:
                    for indice in all_indice:
                        print(check_graph.degree[indice])
                        if check_graph.degree[indice] == 1:
                            extreme = indice
                    oriented_graph = nx.dfs_tree(check_graph, source=extreme)
                    list_edge = list(oriented_graph.edges())
                    new_list = [list_edge[0][0]]
                    new_list += [e[1] for e in oriented_graph.edges()]
                else:
                    new_list = all_indice
                group_without_connection.append(new_list)
        print([[points2.index(e) for e in pts] for pts in group_without_connection])

        for g in group_without_connection:
            print(g)
            extremes = []
            for e in g:
                successor = list(graph_points2[e])
                if successor[0] not in g:
                    extremes.append(successor[0])
                if successor[1] not in g:
                    extremes.append(successor[1])
            first_point = extremes[0]
            last_point = extremes[1]
            print('first_point', first_point)
            print('last_point', last_point)
            path = nx.shortest_path(graph, source=first_point, target=last_point)

            length_path = len(path) - 2
            print(len(path))

            if length_path == 2:
                indice = 1
                for i, e in enumerate(g):
                    if i == 0:
                        graph.add_edges_from([(e, path[indice])])
                        indice += 1
                    else:
                        graph.add_edges_from([(e, path[indice])])
            elif length_path == 1:
                # nx.draw(graph)
                # raise
                indice = 1
                for i, e in enumerate(g):
                    print('e', e)
                    print('path', path[indice])
                    graph.add_edges_from([(e, path[indice])])

            else:
                raise KeyError('not implemented')

        for p1, p2 in zip(points2, points2[1:] + points2[0:1]):
            graph.add_edges_from([(p1, p2)])

        pos_x = 0.05
        for pt in points2:
            if pt not in dict_pos:
                dict_pos[pt] = np.array([pos_x, 0])
                pos_x += 0.1
                print('add point', pt)
        # dict_pos = nx.kamada_kawai_layout(graph)
        # nx.draw(graph, pos=dict_pos)
        return graph

    def _search_cycles(self, polygon1, polygon2):
        points1 = polygon1.points
        points2 = polygon2.points
        graph = self._graph(points1, points2)

        elementary_edges = []
        for pt1 in points1:
            all_cycles = nx.cycle_basis(graph, root=pt1)
            for cycle in all_cycles:
                if cycle in elementary_edges or cycle[::-1] in elementary_edges:
                    continue
                if len(cycle) == 3:
                    elementary_edges.append(cycle)
                elif len(cycle) == 4:
                    elementary_edges.append(cycle)

        print('we obtain {} different cycle to analyze'.format(len(elementary_edges)))
        return elementary_edges

    def _map_with_faces(self, elementary_edges: List[List[vm.Point3D]]):

        def if_face_plane(point1, point2, point3, point4):
            vect1 = point2 - point1
            vect2 = point3 - point1
            vect3 = vect1.cross(vect2)
            dot = vect3.dot(point4 - point1)
            if abs(dot) > 1e-20:
                return False
            return True

        def triangles_with_4_points(point1, point2, point3, point4):
            triangles = []
            triangles.append(vm.faces.Triangle3D(point1, point2, point3))
            triangles.append(vm.faces.Triangle3D(point1, point3, point4))
            return triangles

        def face_plane_with_4_points(point1, point2, point3, point4):
            vect1 = point2 - point1
            vect1.normalize()
            vect2 = point4 - point1
            vect2.normalize()
            vect3 = vect1.cross(vect2)
            vect4 = vect1.cross(vect3)
            pt1 = vm.Point2D(vect1.dot(point1 - point1), vect4.dot(point1 - point1))
            pt2 = vm.Point2D(vect1.dot(point2 - point1), vect4.dot(point2 - point1))
            pt3 = vm.Point2D(vect1.dot(point3 - point1), vect4.dot(point3 - point1))
            pt4 = vm.Point2D(vect1.dot(point4 - point1), vect4.dot(point4 - point1))
            plane = vm.faces.Plane3D(vm.Frame3D(point1, vect1, vect4, vect3))
            contour = vm.wires.ClosedPolygon2D([pt1, pt2, pt3, pt4])
            # contour.plot()
            plane_face = vm.faces.PlaneFace3D(plane, vm.faces.Surface2D(contour, []))
            return plane_face

        faces = []
        for coords in elementary_edges:
            # for p in coords:
            #     faces.append(vm.primitives3d.Sphere(p, 0.001))
            if len(coords) == 3:
                faces.append(vm.faces.Triangle3D(*coords))
            elif len(coords) == 4:
                if if_face_plane(*coords):
                    # new_coords = order_points(*coords)
                    faces.extend(triangles_with_4_points(*coords))
                else:
                    faces.extend(triangles_with_4_points(*coords))
        return faces

    def _def_plane_face(self, center, vect1, vect2, normal, dict_orientation):
        plane = vm.faces.Plane3D(vm.Frame3D(center, vect1, vect2, normal))
        contour = vm.wires.ClosedPolygon2D([vm.Point2D(dict_orientation[vect1][0], dict_orientation[vect2][0]),
                                            vm.Point2D(dict_orientation[vect1][0], dict_orientation[vect2][1]),
                                            vm.Point2D(dict_orientation[vect1][1], dict_orientation[vect2][1]),
                                            vm.Point2D(dict_orientation[vect1][1], dict_orientation[vect2][0])])
        contour_offset = contour.offset(-0.01)
        plane_face = vm.faces.PlaneFace3D(plane, vm.faces.Surface2D(contour_offset, []))
        return plane_face

    def _initial_face(self, polygon, position, normal, vect1, vect2, plane_face):
        for p1, p2 in zip(polygon.points[0:-1], polygon.points[1:]):
            if p1.point_distance(p2) < 1e-6:
                raise
        surf2d = vm.faces.Surface2D(polygon.to_2d(position * normal, vect1, vect2), [])
        return vm.faces.PlaneFace3D(plane_face, surf2d)

    def generate_polygon(self, normal: vm.Vector3D,
                         number_contour: int = 5,
                         number_points: int = 10):
        bb = self.volume_model.bounding_box
        xmin, xmax, ymin, ymax, zmin, zmax = bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax
        dict_orientation = {vm.X3D: [xmin, xmax],
                            vm.Y3D: [ymin, ymax],
                            vm.Z3D: [zmin, zmax]}
        all_vect = list(dict_orientation.keys())
        indice_normal = all_vect.index(normal)
        new_axis = list(range(3))
        new_axis.remove(indice_normal)
        vect1 = all_vect[new_axis[0]]
        vect2 = all_vect[new_axis[1]]

        pas_trans = (dict_orientation[normal][1] - dict_orientation[normal][0]) / number_contour
        center = vm.Point3D(0, 0, 0)
        center.translation(dict_orientation[normal][0] * normal, copy=False)
        polygons = []

        for p in [dict_orientation[normal][0] + pas_trans * i for i in range(number_contour)]:
            center.translation(pas_trans * normal, copy=False)
            plane_face = self._def_plane_face(center, vect1, vect2, normal, dict_orientation)
            primitives = self.contour_intersection(plane_face)

            if primitives is None:
                continue
            # self._draw_2d(primitives, plane_face)
            # self._draw_3d(plane_face)
            polygon = self._contour_order(primitives, plane_face)
            if polygon is None:
                continue
            poly = self._update_orientation(polygon, plane_face)
            # poly = poly.simplify()
            # points = poly.points
            # poly.plot()
            points = self._discretize_polygon(poly)
            new_poly = vm.wires.ClosedPolygon3D(points)
            # new_poly.plot()
            # length = poly.length()
            # pas = length / number_points
            # points = [poly.point_at_abscissa(i * pas) for i in range(number_points)]
            polygons.append(vm.wires.ClosedPolygon3D(points))
            # polygons[-1].plot()

        polygons_simplified = polygons
        return polygons_simplified

    def generate_shell(self, polygons: List[vm.wires.ClosedPolygon3D],
                       normal: vm.Vector3D,
                       number_contour: int = 5,
                       number_points: int = 10):
        bb = self.volume_model.bounding_box
        xmin, xmax, ymin, ymax, zmin, zmax = bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax
        dict_orientation = {vm.X3D: [xmin, xmax],
                            vm.Y3D: [ymin, ymax],
                            vm.Z3D: [zmin, zmax]}
        all_vect = list(dict_orientation.keys())
        indice_normal = all_vect.index(normal)
        new_axis = list(range(3))
        new_axis.remove(indice_normal)
        vect1 = all_vect[new_axis[0]]
        vect2 = all_vect[new_axis[1]]

        pas_trans = (dict_orientation[normal][1] - dict_orientation[normal][0]) / number_contour
        center = vm.Point3D(0, 0, 0)
        pos_initial_plane = polygons[0].primitives[0].start.dot(normal)
        center.translation(pos_initial_plane * normal, copy=False)
        # center.translation(pas_trans * normal, copy=False)
        plane = vm.faces.Plane3D(vm.Frame3D(center, vect1, vect2, normal))

        initial_face = self._initial_face(polygons[0], pos_initial_plane, normal, vect1, vect2, plane)
        faces = [initial_face]

        for i, (poly1, poly2) in enumerate(zip(polygons[0:-1], polygons[1:])):
            output = self._search_cycles(poly1, poly2)
            faces.extend(self._map_with_faces(output))

        center = vm.Point3D(0, 0, 0)
        pos_final_plane = polygons[-1].primitives[0].start.dot(normal)
        center.translation(pos_final_plane * normal, copy=False)
        # center.translation(len(polygons)*pas_trans * normal, copy=False)
        plane = vm.faces.Plane3D(vm.Frame3D(center, vect1, vect2, normal))
        final_face = self._initial_face(polygons[-1],
                                        pos_final_plane,
                                        normal, vect1, vect2, plane)
        faces.append(final_face)

        return vm.faces.ClosedShell3D(faces)
        # return vmc.PointCloud3D.generate_shell(polygons, normal, vect1, vect2)

    def _extrude_one_direction(self, resolution: int, normal: vm.Vector3D):
        points = self.upload_stl.stl.extract_points()
        pointcloud3d = vmc.PointCloud3D(points)
        polygons, vec1, vec2, position_plane = pointcloud3d.to_polygon(resolution=resolution, normal=normal)
        polygon2d = polygons[0].to_2d(position_plane[0] * normal, vec1, vec2)
        return vm.primitives3d.ExtrudedProfile(-1 * normal, vec1, vec2,
                                               polygon2d, [], 2 * normal)

    def _generate_middle_polygon(self, points: List[vm.Point3D], resolution: int,
                                 normal: vm.Vector3D, vect1: vm.Vector3D, vect2: vm.Vector3D):
        pt2d = [p.to_2d(normal, vect1, vect2) for p in points]
        new_pt2d = []
        for p in pt2d:
            check = True
            for pp in new_pt2d:
                if p.point_distance(pp) < 1e-4:
                    check = False
                    break
            if check:
                new_pt2d.append(p)

        # polygon = vm.wires.ClosedPolygon2D.concave_hull(new_pt2d, -0.1, 0.00005)
        polygon = vm.wires.ClosedPolygon2D.convex_hull_points(new_pt2d)

        for p1, p2 in zip(polygon.points[0:-1], polygon.points[1:]):
            if p1.point_distance(p2) < 1e-6:
                raise

        pt_normal = [p.dot(normal) for p in points]
        position_plane = (min(pt_normal) + max(pt_normal)) / 2
        polygon3d = polygon.to_3d(position_plane * normal, vect1, vect2)
        for p1, p2 in zip(polygon3d.points[0:-1], polygon3d.points[1:]):
            if p1.point_distance(p2) < 1e-6:
                raise
        return polygon3d


class BoxesWrap(GiftWrap):
    _standalone_in_db = True
    _non_serializable_attributes = ['volume_model', '_filter']
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _eq_is_data_eq = True

    def __init__(self, volume_model: vm.core.VolumeModel = None,
                 size_boxe: float = None,
                 name: str = ''):
        GiftWrap.__init__(self, volume_model=volume_model, name=name)

        self.bounding_box = self.volume_model.bounding_box
        self.xmin, self.xmax, self.ymin, self.ymax, self.zmin, self.zmax = self.bounding_box.xmin, self.bounding_box.xmax, \
                                                                           self.bounding_box.ymin, self.bounding_box.ymax, \
                                                                           self.bounding_box.zmin, self.bounding_box.zmax
        if size_boxe is None:
            min_dist = min([self.xmax - self.xmin, self.ymax - self.ymin, self.zmax - self.zmin])
            self.size_boxe = min_dist / 10
        else:
            self.size_boxe = size_boxe
        self.dict_orientation = {vm.X3D: [self.xmin, self.xmax],
                                 vm.Y3D: [self.ymin, self.ymax],
                                 vm.Z3D: [self.zmin, self.zmax]}
        self.nx, self.ny, self.nz, self.cx, self.cy, self.cz = self._size_boxes()
        self.dict_size = {vm.X3D: self.cx, vm.Y3D: self.cy, vm.Z3D: self.cz}
        self.dict_number = {vm.X3D: self.nx, vm.Y3D: self.ny, vm.Z3D: self.nz}
        self._filter = {}
        self.adresses = self.smart_box_discretization()

    def _smart_glue_boxes(self, propagation_dir: vm.Vector3D,
                          boxes: List[vm.primitives3d.Block]):
        bb = self.volume_model.bounding_box
        xmin, xmax, ymin, ymax, zmin, zmax = bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax
        dict_orientation = {vm.X3D: [xmin, xmax],
                            vm.Y3D: [ymin, ymax],
                            vm.Z3D: [zmin, zmax]}

        all_vect = [vm.X3D, vm.Y3D, vm.Z3D]
        for v in all_vect:
            if v.dot(propagation_dir) == 1:
                direction = 1
                all_v = all_vect[:]
                all_v.remove(v)
                face_dir = all_v[0]
                line_dir = all_v[1]
                propagation_min_max = self.dict_orientation[v]
                propagation_number = self.dict_number[v]
            elif v.dot(propagation_dir) == -1:
                direction = -1
                all_v = all_vect[:]
                all_v.remove(v)
                face_dir = all_v[0]
                line_dir = all_v[1]
                propagation_min_max = self.dict_orientation[v][::-1]
                propagation_number = self.dict_number[v]
        face_number = self.dict_number[face_dir]
        face_min_max = self.dict_orientation[face_dir]
        line_number = self.dict_number[line_dir]
        line_min_max = self.dict_orientation[line_dir]
        propagation_step = abs(propagation_min_max[1] - propagation_min_max[0]) / propagation_number
        face_step = abs(face_min_max[1] - face_min_max[0]) / face_number
        line_step = abs(line_min_max[1] - line_min_max[0]) / line_number

        boxes_localisation = []
        for i, box in enumerate(boxes):
            origin = box.frame.origin
            ind_prop = int(
                (abs(origin.dot(propagation_dir) - direction * propagation_min_max[
                    0]) - propagation_step / 2) / propagation_step)
            ind_face = int(
                (abs(origin.dot(face_dir) - face_min_max[0]) - face_step / 2) / face_step)
            ind_line = int(
                (abs(origin.dot(line_dir) - line_min_max[0]) - line_step / 2) / line_step)
            boxes_localisation.append([ind_prop, ind_face, ind_line, i])

        boxes_localisation = np.array(boxes_localisation)
        dict_boxes = {}
        polygons = []
        graph = nx.DiGraph()
        for i, k in enumerate(set(boxes_localisation[:, 0])):
            print('analyze ', i)
            select_boxes = np.array(
                [list(boxes_localisation[j]) for j, i in enumerate(boxes_localisation[:, 0]) if i == k])
            all_boxes = [boxes[int(s)] for s in select_boxes[:, 3:4]]

            print('start _connected_boxes')
            elem_all_boxes = self._connected_boxes(all_boxes)
            print(len(elem_all_boxes))
            polygon_temp = []
            for e_all_boxes in elem_all_boxes:
                vol_all_boxes = vm.core.VolumeModel(e_all_boxes)
                polygon = self._glue_boxes_with_shape(e_all_boxes, propagation_dir, face_dir, line_dir)
                plane_face = self._def_plane_face(vol_all_boxes.bounding_box.center, face_dir, line_dir,
                                                  propagation_dir, dict_orientation)
                polygon = self._update_orientation(polygon, plane_face)
                polygon_temp.append(polygon)
                connected_polygons = self._smart_proximity_for_polygon(dict_boxes, e_all_boxes)
                for p in connected_polygons:
                    graph.add_edges_from([(p, polygon)])
                dict_boxes[polygon] = e_all_boxes
            polygons.append(polygon_temp)

        print('start order')
        all_edges = self._define_order_edges(graph, polygons)

        print('start generate shell')
        primitives = []
        for edge in all_edges:
            poly1, poly2 = edge
            if graph.degree[poly1] == 2 and graph.degree[poly2] == 2:
                primitives.append(self.generate_shell([poly1, poly2], propagation_dir))
            else:
                origin = poly1.points[0]
                contour = poly1.to_2d(origin, face_dir, line_dir)
                primitives.append(vm.primitives3d.ExtrudedProfile(origin, face_dir, line_dir,
                                                                  contour, [], propagation_step * propagation_dir))

        primitive_merge = primitives[0]
        for i, s in enumerate(primitives[1:]):
            print(i)
            output = primitive_merge.union(s, tol=1e-6)
            print('number volume', len(output))
            primitive_merge = output[0]
        return primitive_merge

    def _define_order_edges(self, graph, polygons):
        all_edges = []
        for polys in polygons:
            for poly in polys:
                edges = nx.dfs_edges(graph, poly)
                for edge in edges:
                    if edge not in all_edges:
                        all_edges.append(edge)
        return all_edges

    def _smart_proximity_for_polygon(self, past_dict_boxes: Dict[vm.wires.ClosedPolygon3D, List[vm.primitives3d.Block]],
                                     current_boxes: List[vm.primitives3d.Block]):
        group_boxes = []
        polygons = []
        for polygon, boxes in past_dict_boxes.items():
            group_boxes.append(boxes)
            polygons.append(polygon)

        nb_group_boxes = len(group_boxes)
        connected_polygons = []
        for ind_group_box in range(nb_group_boxes):
            group_box = group_boxes[ind_group_box]
            polygon = polygons[ind_group_box]
            for box1, box2 in product(current_boxes, group_box):
                if self._proximity(box1, box2):
                    if polygon not in connected_polygons:
                        connected_polygons.append(polygon)
        return connected_polygons

    def _glue_boxes_with_shape(self, boxes: List[vm.primitives3d.Block], normal: vm.Vector3D,
                                vect1: vm.Vector3D, vect2: vm.Vector3D):
        points = []
        for box in boxes:
            for face in box.faces:
                points.extend(face.outer_contour3d.points)
        polygons = self._generate_middle_polygon(points, 2, normal, vect1, vect2)
        return polygons

    def _analyze_profile(self, indices: List[int], propagation_step: float,
                         transversal_step: float):
        evolution_vector = []
        for i, (ni, nm, np) in enumerate(zip([None] + indices[0:-1], indices, indices[1:] + [None])):
            if ni is None and np is not None:
                if nm - np == 1:
                    evolution_vector.append(vm.Vector2D(0, transversal_step))
                elif nm - np == -1:
                    evolution_vector.append(vm.Vector2D(propagation_step, transversal_step))
                else:
                    evolution_vector.append(vm.Vector2D(0, transversal_step))
                    if np != nm:
                        evolution_vector.append(vm.Vector2D((np - nm) * propagation_step, 0))
            elif ni is not None and np is None:
                if ni - nm == 1:
                    evolution_vector.append(vm.Vector2D(-propagation_step, transversal_step))
                elif ni - nm == -1:
                    evolution_vector.append(vm.Vector2D(0, transversal_step))
                else:
                    evolution_vector.append(vm.Vector2D(0, transversal_step))
            elif ni is None and np is None:
                evolution_vector.append(vm.Vector2D(0, transversal_step))
            else:
                if ni - nm == 1 and nm - np == 1:
                    evolution_vector.append(vm.Vector2D(-propagation_step, transversal_step))
                elif ni - nm == -1 and nm - np == -1:
                    evolution_vector.append(vm.Vector2D(propagation_step, transversal_step))
                elif ni - nm != 1 and nm - np == 1:
                    evolution_vector.append(vm.Vector2D(0, transversal_step))
                elif ni <= nm and nm - np == -1:
                    evolution_vector.append(vm.Vector2D(propagation_step, transversal_step))
                elif ni - nm == 1 and nm >= np:
                    evolution_vector.append(vm.Vector2D(-propagation_step, transversal_step))
                    if np != nm:
                        evolution_vector.append(vm.Vector2D((np - nm) * propagation_step, 0))
                elif ni - nm == -1 and nm - np != -1:
                    evolution_vector.append(vm.Vector2D(0, transversal_step))
                    if np != nm:
                        evolution_vector.append(vm.Vector2D((np - nm) * propagation_step, 0))
                else:
                    evolution_vector.append(vm.Vector2D(0, transversal_step))
                    if np != nm:
                        evolution_vector.append(vm.Vector2D((np - nm) * propagation_step, 0))

        return evolution_vector

    def _glue_boxes_with_line(self, boxes: List[vm.primitives3d.Block]):
        glue_dir = vm.X3D
        other_dir1 = vm.Y3D
        other_dir2 = vm.Z3D
        dict_boxes = {}
        for box in boxes:
            origin = box.frame.origin
            dir1 = origin.dot(other_dir1)
            dir2 = origin.dot(other_dir2)
            pt2 = vm.Point2D(dir1, dir2)
            if pt2 in dict_boxes:
                dict_boxes[pt2].append(box)
            else:
                dict_boxes[pt2] = [box]

        all_boxes = []
        for pt, boxes in dict_boxes.items():
            dist = []
            for box in boxes:
                origin = box.frame.origin
                dist.append(origin.dot(glue_dir))
            sort_dist = np.argsort(np.array(dist))
            last_indice = sort_dist[0]
            groups_boxes = []
            boxes_temp = [boxes[last_indice]]
            if len(boxes) > 1:
                for s in sort_dist[1:]:
                    if self._proximity(boxes[last_indice], boxes[s]):
                        boxes_temp.append(boxes[s])
                    else:
                        all_boxes.append(self._simplify_boxes(boxes_temp))
                        boxes_temp = [boxes[s]]
                    last_indice = s
            all_boxes.append(self._simplify_boxes(boxes_temp))
        return all_boxes

    def _simplify_boxes(self, boxes: List[vm.primitives3d.Block]):
        bb = boxes[0].bounding_box
        for box in boxes[1:]:
            bb += box.bounding_box
        frame = bb.to_frame()
        return vm.primitives3d.Block(frame)

    def _generate_inside_boxes(self):
        all_boxes = self._generate_boxes_one_direction(vm.X3D, False)
        all_boxes_outside = self._generate_boxes()
        all_boxes_inside = []
        for box in all_boxes:
            if box not in all_boxes_outside:
                all_boxes_inside.append(box)
        return all_boxes_inside

    def _generate_boxes(self):
        all_vect = [vm.X3D, -vm.X3D, vm.Y3D, -vm.Y3D, vm.Z3D, -vm.Z3D]
        all_boxes = []
        for vect in all_vect:
            boxes = self._generate_boxes_one_direction(vect)
            for box in boxes:
                if box not in all_boxes:
                    all_boxes.append(box)
        return all_boxes

    def _generate_boxes_one_direction(self, propagation_dir: vm.Vector3D,
                                      check_primitives: bool = True):
        all_vect = [vm.X3D, vm.Y3D, vm.Z3D]
        for v in all_vect:
            if v.dot(propagation_dir) == 1:
                all_v = all_vect[:]
                all_v.remove(v)
                face_dir = all_v[0]
                line_dir = all_v[1]
            elif v.dot(propagation_dir) == -1:
                all_v = all_vect[:]
                all_v.remove(v)
                face_dir = all_v[0]
                line_dir = all_v[1]
        for v in all_vect:
            if v.dot(face_dir) == 1:
                number_face = self.dict_number[v]
            elif v.dot(face_dir) == -1:
                number_face = self.dict_number[v]

        all_boxes = []
        for n in range(number_face):
            boxes = self._generate_face(propagation_dir, face_dir, line_dir, n, check_primitives)
            all_boxes.extend(boxes)
        return all_boxes

    def _generate_face(self, propagation_dir: vm.Vector3D,
                       face_dir: vm.Vector3D,
                       line_dir: vm.Vector3D,
                       face_indice: int,
                       check_primitives: bool = True):

        all_vect = [vm.X3D, vm.Y3D, vm.Z3D]
        for v in all_vect:
            if v.dot(propagation_dir) == 1:
                direction_propagation = 1
                propagation_min_max = self.dict_orientation[v]
                number_propagation = self.dict_number[v]
                size_propagation = self.dict_size[v]
            elif v.dot(propagation_dir) == -1:
                direction_propagation = -1
                propagation_min_max = self.dict_orientation[v][::-1]
                number_propagation = self.dict_number[v]
                size_propagation = self.dict_size[v]
            if v.dot(face_dir) == 1:
                direction_face = 1
                face_min_max = self.dict_orientation[v]
                number_face = self.dict_number[v]
                size_face = self.dict_size[v]
            elif v.dot(face_dir) == -1:
                direction_face = -1
                face_min_max = self.dict_orientation[v][::-1]
                number_face = self.dict_number[v]
                size_face = self.dict_size[v]

        face_position = direction_face * face_dir * (
                face_min_max[0] + direction_face * size_face / 2 + direction_face * face_indice * size_face)
        all_boxes = []
        acceptable_boxes = None
        for n in range(number_propagation):
            propagation_position = direction_propagation * propagation_dir * (propagation_min_max[
                                                                                  0] + direction_propagation * size_propagation / 2 + direction_propagation * n * size_propagation)
            boxes, acceptable_boxes = self._generate_line(line_dir, propagation_dir, face_position,
                                                          propagation_position,
                                                          acceptable_boxes, check_primitives)
            all_boxes.extend(boxes)
        return all_boxes

    def _generate_line(self, line_dir: vm.Vector3D, propagation_dir: vm.Vector3D,
                       face_position: vm.Vector3D,
                       propagation_position: vm.Vector3D,
                       acceptable_boxes: List[int] = None,
                       check_primitives: bool = True):

        all_vect = [vm.X3D, vm.Y3D, vm.Z3D]
        for v in all_vect:
            if v.dot(line_dir) == 1:
                direction = 1
                min_max = self.dict_orientation[v]
                number = self.dict_number[v]
                size = self.dict_size[v]
            elif v.dot(line_dir) == -1:
                direction = -1
                min_max = self.dict_orientation[v][::-1]
                number = self.dict_number[v]
                size = self.dict_size[v]

        boxes = []
        center_init = (min_max[0] + direction * size / 2) * line_dir * direction + face_position + propagation_position
        frame_init = vm.Frame3D(center_init, self.cx * vm.X3D, self.cy * vm.Y3D, self.cz * vm.Z3D)

        if acceptable_boxes is None:
            range_number = range(number)
        else:
            range_number = acceptable_boxes

        next_acceptable_boxes = []
        for n in range_number:
            frame = frame_init.translation(n * size * line_dir, copy=True)
            adresse = self.define_adresses_with_frame(frame)

            if check_primitives:
                if not self._analyze_boxe(adresse, propagation_dir, frame):
                    boxes.append(adresse)
                    next_acceptable_boxes.append(n)
            else:
                boxes.append(adresse)

        if check_primitives:
            return boxes, next_acceptable_boxes
        else:
            return boxes, None

    def _analyze_boxe(self, adresse: Tuple[int, int, int], direction: vm.Vector3D, frame: vm.Frame3D):
        # add = self.define_adresses(box)
        if adresse in self.adresses:
            print(self.adresses[adresse])
            all_indices = self.adresses[adresse]
            all_faces = []
            for (primitive_indice, face_indice) in all_indices:
                all_faces.append(self.volume_model.primitives[primitive_indice].faces[face_indice])
            box = vm.primitives3d.Block(frame)

            for face in all_faces:
                for plane in box.faces:
                    sol = plane.face_intersections(face)
                    if sol:
                        if isinstance(sol[0], volmdlr.wires.Wire3D):
                            return True
            return False
        else:
            return False

    def generate_polygon(self, normal: vm.Vector3D,
                         number_contour: int = 5,
                         number_points: int = 10):
        center = vm.Point3D(0, 0, 0)
        center.translation(self.dict_orientation[normal][0] * normal, copy=False)
        polygons = []

        all_primitives = []
        for p in [self.dict_orientation[normal][0] + pas_trans * i for i in range(number_contour)]:
            center.translation(pas_trans * normal, copy=False)
            plane_face = self._def_plane_face(center, self.vect1, self.vect2, normal, self.dict_orientation)
            primitives = self.contour_intersection(plane_face)
            all_primitives.append(primitives)
            wires = vm.wires.Contour3D(primitives)
            wires.plot()

        return all_primitives

class IntersectionWrap(GiftWrap):
    _standalone_in_db = True
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _eq_is_data_eq = True

    def __init__(self, volume_model: vm.core.VolumeModel = None,
                 name: str = ''):
        DessiaObject.__init__(self, name=name)
        self.volume_model = volume_model

    def _generate_extrude(self):
        primitives = []
        primitives.append(self._extrude_one_direction(2, vm.X3D))
        primitives.append(self._extrude_one_direction(2, vm.Y3D))
        primitives.append(self._extrude_one_direction(2, vm.Z3D))
        output_primitive = primitives[0]
        for primitive in primitives[1:]:
            output = output_primitive.intersection(primitive)
            output_primitive = output[0]
        return output_primitive
