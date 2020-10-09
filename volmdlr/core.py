#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from packaging import version
import warnings
import math
import numpy as npy

npy.seterr(divide='raise')

# from geomdl import NURBS
from geomdl import BSpline
from geomdl import utilities
import matplotlib.pyplot as plt
import mpl_toolkits
from matplotlib.patches import Arc, FancyArrowPatch
from mpl_toolkits.mplot3d import Axes3D
# from mpl_toolkits.mplot3d import proj3d
from matplotlib import __version__ as _mpl_version

import networkx as nx

from volmdlr.core_compiled import (Vector2D, Vector3D, Point2D, Point3D,
                            O2D, X2D, Y2D, OXY,
                            Basis2D, Basis3D, Frame2D, Frame3D,
                            O3D, X3D, Y3D, Z3D,
                            LineSegment2DPointDistance,
                            PolygonPointBelongs, Matrix22
                            )

from scipy.linalg import solve
from scipy.spatial import Delaunay

import volmdlr.geometry as geometry
from volmdlr import plot_data
# from volmdlr import triangulation as tri
import triangle  # doc : https://rufat.be/triangle/

import dessia_common as dc
from jinja2 import Environment, PackageLoader, select_autoescape

import webbrowser
import os
import tempfile
import subprocess
import random

import scipy as scp
import scipy.optimize

from typing import List, Tuple

# import volmdlr.faces3d
# import volmdlr.surfaces3d as surfaces3d
# import volmdlr.primitives3D

two_pi = 2*math.pi

def standardize_knot_vector(knot_vector):
    u0 = knot_vector[0]
    u1 = knot_vector[-1]
    standard_u_knots = []
    if u0 != 0 or u1 != 1:
        x = 1 / (u1 - u0)
        y = u0 / (u0 - u1)
        for u in knot_vector:
            standard_u_knots.append(u * x + y)
        return standard_u_knots
    else:
        return knot_vector

def find_and_replace(string, find, replace):
    """
    Finds a string in a string and replace it
    """
    index = string.find(find)
    if index != -1:
        try:
            # verifie si il reste pas des chiffre apres
            int(string[index + len(find)])
        except (ValueError, IndexError):
            # on remplace
            return string[:index] + replace + string[index + len(find):]
        else:
            return string[:index] + find_and_replace(string[index + len(find)],
                                                     find, replace)
    return string



def set_to_list(step_set):
    char_list = step_set.split(',')
    char_list[0] = char_list[0][1:]
    char_list[-1] = char_list[-1][:-1]
    return [elem for elem in char_list]


def delete_node_and_predecessors(graph, node):
    predecessors = list(graph.predecessors(node))
    graph.remove_node(node)
    for predecessor in predecessors:
        delete_node_and_predecessors(graph, predecessor)


def delete_node_and_successors(graph, node):
    successors = list(graph.successors(node))
    graph.remove_node(node)
    for successor in successors:
        delete_node_and_successors(graph, successor)


def clockwise_angle(vector1, vector2):
    """
    Return the clockwise angle in radians between vector1 and vector2.
    """
    vector0 = Vector2D((0, 0))
    if vector0 in (vector1, vector2):
        return 0

    dot = vector1.Dot(vector2)
    norm_vec_1 = vector1.Norm()
    norm_vec_2 = vector2.Norm()
    sol = dot / (norm_vec_1 * norm_vec_2)
    cross = vector1.Cross(vector2)
    if math.isclose(sol, 1, abs_tol=1e-6):
        inner_angle = 0
    elif math.isclose(sol, -1, abs_tol=1e-6):
        inner_angle = math.pi
    else:
        inner_angle = math.acos(sol)

    if cross < 0:
        return inner_angle

    return two_pi - inner_angle


def vectors3d_angle(vector1, vector2):
    dot = vector1.Dot(vector2)
    theta = math.acos(dot / (vector1.Norm() * vector2.Norm()))

    return theta


def sin_cos_angle(u1, u2):
    # You have to give cos(theta)=u1, sin(theta)=u2
    # return an angle between 0 and 2pi
    if u1 < -1:
        u1 = -1
    elif u1 > 1:
        u1 = 1
    if u2 < -1:
        u2 = -1
    elif u2 > 1:
        u2 = 1

    if u1 > 0:
        if u2 >= 0:
            theta = math.acos(u1)
        else:
            theta = two_pi + math.asin(u2)
    else:
        if u2 >= 0:
            theta = math.acos(u1)
        else:
            theta = two_pi - math.acos(u1)
    return theta


def delete_double_pos(points, triangles):
    points_to_indexes = {}

    for index, point in enumerate(points):
        if point not in points_to_indexes:
            points_to_indexes[point] = [index]
        else:
            points_to_indexes[point].append(index)

    new_points = []
    index_to_modified_index = {}
    for i, (point, indexes) in enumerate(points_to_indexes.items()):
        new_points.append(point)
        index_to_modified_index[indexes[0]] = i

    index_to_new_index = {}

    for indexes in points_to_indexes.values():
        for index in indexes[1:]:
            index_to_new_index[index] = indexes[0]

    new_triangles = []
    for face_triangles in triangles:
        if face_triangles is None:
            continue
        # print('face_triangles', face_triangles)
        new_face_triangles = []
        for triangle_ in face_triangles:
            # print('triangle', triangle)
            new_triangle = []
            for index in triangle_:
                if index in index_to_new_index:
                    modified_index = index_to_modified_index[
                        index_to_new_index[index]]
                else:
                    modified_index = index_to_modified_index[index]
                new_triangle.append(modified_index)
            new_face_triangles.append(tuple(new_triangle))
        new_triangles.append(new_face_triangles)

    return new_points, new_triangles


def determinant(vec1, vec2, vec3):
    a = npy.array((vec1.vector, vec2.vector, vec3.vector))
    return npy.linalg.det(a)


def delete_double_point(list_point):
    points = []
    for pt in list_point:
        if pt not in points:
            points.append(pt)
        else:
            continue
    return points


def max_pos(list_of_float):
    pos_max, max_float = 0, list_of_float[0]
    for pos, fl in enumerate(list_of_float):
        if pos == 0:
            continue
        else:
            if fl > max_float:
                max_float = fl
                pos_max = pos
    return max_float, pos_max


def min_pos(list_of_float):
    pos_min, min_float = 0, list_of_float[0]
    for pos, fl in enumerate(list_of_float):
        if pos == 0:
            continue
        else:
            if fl < min_float:
                min_float = fl
                pos_min = pos
    return min_float, pos_min


def check_singularity(all_points):
    plus_pi, moins_pi = [], []
    for enum, pt in enumerate(all_points):
        if pt.vector[0] > math.pi * 1.01:
            plus_pi.append(enum)
        elif pt.vector[0] < math.pi * 0.99:
            moins_pi.append(enum)

    if len(moins_pi) <= 2 and len(all_points) > 4:
        for pos in moins_pi:
            new_pt = all_points[pos].copy() + Point2D((two_pi, 0))
            if new_pt.vector[0] > two_pi:
                new_pt.vector[0] = two_pi
            all_points[pos] = new_pt
    elif len(plus_pi) <= 2 and len(all_points) > 4:
        for pos in plus_pi:
            new_pt = all_points[pos].copy() - Point2D((two_pi, 0))
            if new_pt.vector[0] < 0:
                new_pt.vector[0] = 0
            all_points[pos] = new_pt
    if 3 * len(moins_pi) <= len(plus_pi) and len(all_points) > 4:
        for pos in moins_pi:
            new_pt = all_points[pos].copy() + Point2D((two_pi, 0))
            if new_pt.vector[0] > two_pi:
                new_pt.vector[0] = two_pi
            all_points[pos] = new_pt
    elif 3 * len(plus_pi) <= len(moins_pi) and len(all_points) > 4:
        for pos in plus_pi:
            new_pt = all_points[pos].copy() - Point2D((two_pi, 0))
            if new_pt.vector[0] < 0:
                new_pt.vector[0] = 0
            all_points[pos] = new_pt

    return all_points


def posangle_arc(start, end, radius, frame=None):
    if frame is None:
        p1_new, p2_new = start, end
    else:
        p1_new, p2_new = frame.NewCoordinates(start), frame.NewCoordinates(end)
    # Angle pour le p1
    u1, u2 = p1_new.vector[0] / radius, p1_new.vector[1] / radius
    theta1 = sin_cos_angle(u1, u2)
    # Angle pour le p2
    u3, u4 = p2_new.vector[0] / radius, p2_new.vector[1] / radius
    theta2 = sin_cos_angle(u3, u4)

    if math.isclose(theta1, theta2, abs_tol=1e-6):
        if math.isclose(theta2, 0, abs_tol=1e-6):
            theta2 += two_pi
        elif math.isclose(theta1, two_pi, abs_tol=1e-6):
            theta1 -= two_pi

    return theta1, theta2


def offset_angle(trigo, angle_start, angle_end):
    if trigo:
        offset = angle_start
    else:
        offset = angle_end
    if angle_start > angle_end:
        angle = angle_start - angle_end
    else:
        angle = angle_end - angle_start
    return offset, angle


def angle_principal_measure(angle, min_angle=-math.pi):
    """
    returns angle between O and 2 pi
    """
    max_angle = min_angle + two_pi
    angle = angle % (two_pi)

    if math.isclose(angle, min_angle, abs_tol=1e-9):
        return min_angle
    if math.isclose(angle, max_angle, abs_tol=1e-9):
        return max_angle

    return angle


class Primitive2D(dc.DessiaObject):
    def __init__(self, name=''):
        self.name = name

        dc.DessiaObject.__init__(self, name=name)


class CompositePrimitive2D(Primitive2D):
    """
    A collection of simple primitives
    """

    def __init__(self, primitives, name=''):
        Primitive2D.__init__(self, name)
        self.primitives = primitives
        self.UpdateBasisPrimitives()

    def UpdateBasisPrimitives(self):
        basis_primitives = []
        for primitive in self.primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.basis_primitives)
            else:
                basis_primitives.append(primitive)

        self.basis_primitives = basis_primitives

    def Rotation(self, center, angle, copy=True):
        if copy:
            return self.__class__([p.Rotation(center, angle, copy=True) \
                                   for p in self.primitives])
        else:
            for p in self.primitives:
                p.Rotation(center, angle, copy=False)
            self.UpdateBasisPrimitives()

    def Translation(self, offset, copy=True):
        if copy:
            return self.__class__([p.Translation(offset, copy=True) \
                                   for p in self.primitives])
        else:
            for p in self.primitives:
                p.Translation(offset, copy=False)
            self.UpdateBasisPrimitives()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            return self.__class__([p.frame_mapping(frame, side, copy=True) \
                                   for p in self.primitives])
        else:
            for p in self.primitives:
                p.frame_mapping(frame, side, copy=False)
            self.UpdateBasisPrimitives()

    def To3D(self, plane_origin, x, y, name=None):
        if name is None:
            name = '3D of {}'.format(self.name)
        primitives3D = [p.To3D(plane_origin, x, y) for p in self.primitives]
        return CompositePrimitive3D(primitives3D, name)

    def MPLPlot(self, ax=None, color='k', arrow=False, width=None,
                plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        # else:
        #     fig = ax.figure

        for element in self.primitives:
            if element.__class__.__name__ == 'LineSegment2D':
                element.MPLPlot(ax, color, arrow, width,
                                plot_points=plot_points)
            else:

                element.MPLPlot(ax, color=color)

        ax.margins(0.1)
        plt.show()

        return ax

    def plot_data(self, name, fill=None, color='black', stroke_width=0.2,
                  opacity=1):
        plot_data = {}
        plot_data['fill'] = fill
        plot_data['name'] = name
        plot_data['type'] = 'wire'
        plot_data['plot_data'] = []
        for item in self.primitives:
            plot_data['plot_data'].append(item.plot_data(color=color,
                                                         stroke_width=stroke_width,
                                                         opacity=opacity))
        return plot_data







class Primitive3D(dc.DessiaObject):
    def __init__(self, color=None, alpha=0.5, name=''):
        self.color = color
        self.alpha = alpha

        dc.DessiaObject.__init__(self, name=name)

    def volmdlr_primitives(self):
        return [self]

    def babylon_meshes(self):
        mesh = self.triangulation()
        positions, indices = mesh.to_babylon()

        babylon_mesh = {'positions': positions,
                        'indices': indices,
                        'alpha': self.alpha,
                        'name': self.name
                        }

        if self.color is None:
            babylon_mesh['color'] = [0.8, 0.8, 0.8]
        else:
            babylon_mesh['color'] = list(self.color)

        return [babylon_mesh]



class CompositePrimitive3D(Primitive3D):
    _standalone_in_db = True
    _generic_eq = True
    _non_serializable_attributes = ['basis_primitives']
    _non_eq_attributes = ['name', 'basis_primitives']
    _non_hash_attributes = []
    """
    A collection of simple primitives3D
    """

    def __init__(self, primitives, name=''):
        self.primitives = primitives

        Primitive3D.__init__(self, name=name)



    def UpdateBasisPrimitives(self):
        # TODO: This is a copy/paste from CompositePrimitive2D, in the future make a Common abstract class
        basis_primitives = []
        for primitive in self.primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.primitives)
            else:
                basis_primitives.append(primitive)

        self.primitives = basis_primitives

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        for primitive in self.primitives:
            primitive.MPLPlot(ax)

        return ax




class Shell3D(CompositePrimitive3D):
    _standalone_in_db = True
    _generic_eq = True
    _non_serializable_attributes = ['bounding_box']
    _non_eq_attributes = ['name', 'color', 'alpha' 'bounding_box', 'contours']
    _non_hash_attributes = []

    def __init__(self, faces, color=None, alpha=1., name=''):
        self.faces = faces
        self.name = name
        self.color = color
        self.alpha = alpha
        self.bounding_box = self._bounding_box()

    def __hash__(self):
        return sum([hash(f) for f in self.faces])

    def __eq__(self, other_):
        if self.__class__ != other_.__class__:
            return False
        equal = True
        for face, other_face in zip(self.faces, other_.faces):
            equal = (equal and face == other_face)
        return equal

    @classmethod
    def from_step(cls, arguments, object_dict):
        faces = []
        for face in arguments[1]:
            faces.append(object_dict[int(face[1:])])
        return cls(faces, name=arguments[0][1:-1])

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_faces = [face.Rotation(center, axis, angle, copy=True) for face
                         in self.faces]
            return Shell3D(new_faces, name=self.name)
        else:
            for face in self.faces:
                face.Rotation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def Translation(self, offset, copy=True):
        if copy:
            new_faces = [face.Translation(offset, copy=True) for face in
                         self.faces]
            return Shell3D(new_faces, name=self.name)
        else:
            for face in self.faces:
                face.Translation(offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_faces = [face.frame_mapping(frame, side, copy=True) for face in
                         self.faces]
            return Shell3D(new_faces, name=self.name)
        else:
            for face in self.faces:
                face.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()

    def copy(self):
        new_faces = [face.copy() for face in self.faces]
        return Shell3D(new_faces, name=self.name)

    def union(self, shell2):
        new_faces = [face for face in self.faces + shell2.faces]
        new_name = self.name + ' union ' + shell2.name
        new_color = self.color
        return Shell3D(new_faces, name=new_name, color=new_color)

    def volume(self):
        """
        Does not consider holes
        """
        volume = 0
        for i, face in enumerate(self.faces):
            points_3D, triangles_indexes = face.triangulation()
            for triangle_indexes in triangles_indexes[0]:
                point1 = points_3D[triangle_indexes[0]]
                point2 = points_3D[triangle_indexes[1]]
                point3 = points_3D[triangle_indexes[2]]

                v321 = point3[0] * point2[1] * point1[2]
                v231 = point2[0] * point3[1] * point1[2]
                v312 = point3[0] * point1[1] * point2[2]
                v132 = point1[0] * point3[1] * point2[2]
                v213 = point2[0] * point1[1] * point3[2]
                v123 = point1[0] * point2[1] * point3[2]
                volume_tetraedre = 1 / 6 * (
                            -v321 + v231 + v312 - v132 - v213 + v123)

                volume += volume_tetraedre

        return abs(volume)

    def _bounding_box(self):
        """
        Returns the boundary box
        """
        points = []
        bbox = self.faces[0]._bounding_box()

        for face in self.faces:
            bbox += face._bounding_box()

        return bbox

    def point_belongs(self, point, nb_rays=1):
        """
        Ray Casting algorithm
        Returns True if the point is inside the Shell, False otherwise
        """
        epsilon = 10

        bbox = self.bounding_box
        if point[0] < bbox.xmin or point[0] > bbox.xmax:
            return False
        if point[1] < bbox.ymin or point[1] > bbox.ymax:
            return False
        if point[2] < bbox.zmin or point[2] > bbox.zmax:
            return False

        rays = []
        for k in range(0, nb_rays):
            rays.append(LineSegment3D(point, Point3D((random.uniform(0,
                                                                     1) * epsilon,
                                                      random.uniform(0,
                                                                     1) * epsilon,
                                                      random.uniform(0,
                                                                     1) * epsilon))))

        rays = sorted(rays, key=lambda ray: ray.Length())

        rays_intersections = []
        tests = []

        # for ray in rays[:3]:
        for ray in rays[:nb_rays]:
            count = 0
            ray_intersection = []
            is_inside = True
            for face in self.faces:
                intersection_point = face.linesegment_intersection(ray)
                if intersection_point is not None:
                    ray_intersection.append(intersection_point)
                    count += 1
            if count % 2 == 0:
                is_inside = False
            tests.append(is_inside)
            rays_intersections.append(ray_intersection)

        for test1, test2 in zip(tests[:-1], tests[1:]):
            if test1 != test2:
                raise ValueError

        return tests[0]

    def is_inside_shell(self, shell2):
        """
        Returns True if all the points of self are inside shell2 and no face \
        are intersecting
        """
        bbox1 = self.bounding_box
        bbox2 = shell2.bounding_box
        if not bbox1.is_inside_bbox(bbox2):
            return False

        points = []
        for face in self.faces:
            points.extend(face.contours3d[0].tessel_points)

        for point in points:
            if not shell2.point_belongs(point):
                return False

        # Check if any faces are intersecting
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    #                    print('Two faces are intersecting :', face1, face2)
                    return False

        return True

    def shell_intersection(self, shell2):
        """
        Return None if disjointed
        Return (1, 0) or (0, 1) if one is inside the other
        Return (n1, n2) if intersection

        4 cases :
            (n1, n2) with face intersection             => (n1, n2)
            (0, 0) with face intersection               => (0, 0)
            (0, 0) with no face intersection            => None
            (1, 0) or (0, 1) with no face intersection  => 1
        """
        # Check if boundary boxes don't intersect
        bbox1 = self.bounding_box
        bbox2 = shell2.bounding_box
        if not bbox1.bbox_intersection(bbox2):
            #            print("No intersection of shells' BBox")
            return None

        # Check if any point of the first shell is in the second shell
        points1 = []
        for face in self.faces:
            points1.extend(face.contours3d[0].tessel_points)
        points2 = []
        for face in shell2.faces:
            points2.extend(face.contours3d[0].tessel_points)

        nb_pts1 = len(points1)
        nb_pts2 = len(points2)
        compteur1 = 0
        compteur2 = 0
        for point1 in points1:
            if shell2.point_belongs(point1):
                compteur1 += 1
        for point2 in points2:
            if self.point_belongs(point2):
                compteur2 += 1

        inter1 = compteur1 / nb_pts1
        inter2 = compteur2 / nb_pts2
        #        print('shell intersection')
        #        print('shell1 intersecte shell2 à', inter1*100, '%')
        #        print('shell2 intersecte shell1 à', inter2*100, '%')

        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    #                    print('Two faces are intersecting :', face1, face2)
                    #                    ax = face1.plot()
                    #                    face2.plot(ax)
                    return inter1, inter2
        if (inter1, inter2) == (0, 0):
            return None
        return 1

    def distance_to_shell(self, shell2, add_to_volumemodel=None):
        """
        Returns a Mesure object if the distance is not zero, otherwise returns None
        """

        if self.shell_intersection(
                shell2) is not None and self.shell_intersection(shell2) != 1:
            return None

        # distance_min, point1_min, point2_min = self.faces[0].distance_to_face(shell2.faces[0], return_points=True)
        distance_min, point1_min, point2_min = self.faces[0].minimum_distance(
            shell2.faces[0], return_points=True)
        for face1 in self.faces:
            bbox1 = face1.bounding_box
            for face2 in shell2.faces:
                bbox2 = face2.bounding_box
                bbox_distance = bbox1.distance_to_bbox(bbox2)
                if bbox_distance < distance_min:
                    # distance, point1, point2 = face1.distance_to_face(face2, return_points=True)
                    distance, point1, point2 = face1.minimum_distance(face2,
                                                                      return_points=True)
                    if distance == 0:
                        return None
                    elif distance < distance_min:
                        distance_min, point1_min, point2_min = distance, point1, point2

        mesure = Measure3D(point1_min, point2_min)

        if add_to_volumemodel is not None:
            add_to_volumemodel.primitives.append(mesure)

        return mesure

    def distance_to_point(self, point, add_to_volumemodel=None):
        """
        Computes the distance of a point to a Shell3D, whether it is inside or outside the Shell3D
        """
        distance_min, point1_min = self.faces[0].distance_to_point(point,
                                                                   return_other_point=True)
        for face in self.faces[1:]:
            bbox_distance = self.bounding_box.distance_to_point(point)
            if bbox_distance < distance_min:
                distance, point1 = face.distance_to_point(point,
                                                          return_other_point=True)
                if distance < distance_min:
                    distance_min, point1_min = distance, point1

        mesure = Measure3D(point, point1_min)

        if add_to_volumemodel is not None:
            add_to_volumemodel.primitives.append(mesure)

        return mesure

    def intersection_internal_aabb_volume(self, shell2):
        """
        aabb made of the intersection points and the points of self internal to shell2
        """
        intersections_points = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    intersections_points.extend(intersection_points)

        shell1_points_inside_shell2 = []
        for face in self.faces:
            for point in face.contours3d[0].tessel_points:
                if shell2.point_belongs(point):
                    shell1_points_inside_shell2.append(point)

        if len(intersections_points + shell1_points_inside_shell2) == 0:
            return 0
        bbox = BoundingBox.from_points(
            intersections_points + shell1_points_inside_shell2)
        return bbox.volume()

    def intersection_external_aabb_volume(self, shell2):
        """
        aabb made of the intersection points and the points of self external to shell2
        """
        intersections_points = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    intersections_points.extend(intersection_points)

        shell1_points_outside_shell2 = []
        for face in self.faces:
            for point in face.contours3d[0].tessel_points:
                if not shell2.point_belongs(point):
                    shell1_points_outside_shell2.append(point)

        if len(intersections_points + shell1_points_outside_shell2) == 0:
            return 0
        bbox = BoundingBox.from_points(
            intersections_points + shell1_points_outside_shell2)
        return bbox.volume()

    def triangulation(self):
        positions = []
        indices = []

        nb_points = 0
        mesh = DisplayMesh3D([], [])
        for i, face in enumerate(self.faces):
            face_mesh = face.triangulation()
            mesh += face_mesh
        return mesh

    def babylon_script(self, name='primitive_mesh'):
        s = 'var {} = new BABYLON.Mesh("{}", scene);\n'.format(name, name)

        mesh = self.babylon_meshes()[0]

        s += 'var positions = {};\n'.format(mesh['positions'])
        s += 'var indices = {};\n'.format(mesh['indices'])
        s += 'var normals = [];\n'
        s += 'var vertexData = new BABYLON.VertexData();\n'
        s += 'BABYLON.VertexData.ComputeNormals(positions, indices, normals);\n'
        s += 'vertexData.positions = positions;\n'
        s += 'vertexData.indices = indices;\n'
        s += 'vertexData.normals = normals;\n'
        s += 'vertexData.applyToMesh({});\n'.format(name)
        s += '{}.enableEdgesRendering(0.9);\n'.format(name)
        s += '{}.edgesWidth = 0.1;\n'.format(name)
        s += '{}.edgesColor = new BABYLON.Color4(0, 0, 0, 0.6);\n'.format(name)
        s += 'var mat = new BABYLON.StandardMaterial("mat", scene);\n'
        #        s += 'mat.diffuseColor = BABYLON.Color3.Green();\n'
        #        s += 'mat.specularColor = new BABYLON.Color3(0.5, 0.6, 0.87);\n'
        #        s += 'mat.emissiveColor = new BABYLON.Color3(1, 1, 1);\n'
        #        s += 'mat.ambientColor = new BABYLON.Color3(0.23, 0.98, 0.53);\n'
        s += 'mat.backFaceCulling = false;\n'
        s += 'mat.alpha = {};\n'.format(self.alpha)
        s += '{}.material = mat;\n'.format(name)
        if self.color is not None:
            s += 'mat.diffuseColor = new BABYLON.Color3({}, {}, {});\n'.format(
                *self.color)
        return s


class BoundingBox(dc.DessiaObject):
    """
    An axis aligned boundary box
    """

    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax, name=''):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.points = [Point3D((self.xmin, self.ymin, self.zmin)), \
                       Point3D((self.xmax, self.ymin, self.zmin)), \
                       Point3D((self.xmax, self.ymax, self.zmin)), \
                       Point3D((self.xmin, self.ymax, self.zmin)), \
                       Point3D((self.xmin, self.ymin, self.zmax)), \
                       Point3D((self.xmax, self.ymin, self.zmax)), \
                       Point3D((self.xmax, self.ymax, self.zmax)), \
                       Point3D((self.xmin, self.ymax, self.zmax))]
        self.center = (self.points[0] + self.points[-2]) / 2
        self.name = name

    def __hash__(self):
        return sum([hash(p) for p in self.points])

    def __add__(self, other_bbox):
        return BoundingBox(min(self.xmin, other_bbox.xmin),
                           max(self.xmax, other_bbox.xmax),
                           min(self.ymin, other_bbox.ymin),
                           max(self.ymax, other_bbox.ymax),
                           min(self.zmin, other_bbox.zmin),
                           max(self.zmax, other_bbox.zmax))

    def plot(self, ax=None, color=''):
        fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')

        bbox_edges = [[self.points[0], self.points[1]],
                      [self.points[0], self.points[3]],
                      [self.points[0], self.points[4]],
                      [self.points[1], self.points[2]],
                      [self.points[1], self.points[5]],
                      [self.points[2], self.points[3]],
                      [self.points[2], self.points[6]],
                      [self.points[3], self.points[7]],
                      [self.points[4], self.points[5]],
                      [self.points[5], self.points[6]],
                      [self.points[6], self.points[7]],
                      [self.points[7], self.points[4]]]

        x = [p[0] for p in self.points]
        y = [p[1] for p in self.points]
        z = [p[2] for p in self.points]
        ax.scatter(x, y, z)
        for edge in bbox_edges:
            ax.plot3D([edge[0][0], edge[1][0]],
                      [edge[0][1], edge[1][1]],
                      [edge[0][2], edge[1][2]],
                      'gray')
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.show()
        return ax

    @classmethod
    def from_points(cls, points):
        # if len(points) == 0:
        #     return (0, 0, 0, 0, 0, 0)
        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])
        return cls(xmin, xmax, ymin, ymax, zmin, zmax)

    def volume(self):
        return (self.xmax - self.xmin) * (self.ymax - self.ymin) * (
                    self.zmax - self.zmin)

    def bbox_intersection(self, bbox2):
        return (self.xmin < bbox2.xmax and self.xmax > bbox2.xmin \
                and self.ymin < bbox2.ymax and self.ymax > bbox2.ymin \
                and self.zmin < bbox2.zmax and self.zmax > bbox2.zmin)

    def is_inside_bbox(self, bbox2):
        return (self.xmin > bbox2.xmin and self.xmax < bbox2.xmax \
                and self.ymin > bbox2.ymin and self.ymax < bbox2.ymax \
                and self.zmin > bbox2.zmin and self.zmax < bbox2.zmax)

    def intersection_volume(self, bbox2):
        if not self.bbox_intersection(bbox2):
            return 0
        if self.is_inside_bbox(bbox2) or bbox2.is_inside_bbox(self):
            return min(self.volume(), bbox2.volume())

        lx = min(self.xmax, bbox2.xmax) - max(self.xmin, bbox2.xmin)
        ly = min(self.ymax, bbox2.ymax) - max(self.ymin, bbox2.ymin)
        lz = min(self.zmax, bbox2.zmax) - max(self.zmin, bbox2.zmin)

        return lx * ly * lz

    # def intersection_volume(self, bbox2):
    #     if not self.bbox_intersection(bbox2):
    #         return 0
    #
    #     permute_bbox1 = self
    #     permute_bbox2 = bbox2
    #
    #     if permute_bbox2.xmin < permute_bbox1.xmin:
    #         permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
    #     lx = permute_bbox1.xmax - permute_bbox2.xmin
    #
    #     if permute_bbox2.ymin < permute_bbox1.ymin:
    #         permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
    #     ly = permute_bbox1.ymax - permute_bbox2.ymin
    #
    #     if permute_bbox2.zmin < permute_bbox1.zmin:
    #         permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
    #     lz = permute_bbox1.zmax - permute_bbox2.zmin
    #
    #     return lx*ly*lz

    def distance_to_bbox(self, bbox2):
        if self.bbox_intersection(bbox2):
            return 0

        permute_bbox1 = self
        permute_bbox2 = bbox2

        if permute_bbox2.xmin < permute_bbox1.xmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dx = permute_bbox2.xmin - permute_bbox1.xmax
        if dx < 0:
            dx = 0

        if permute_bbox2.ymin < permute_bbox1.ymin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dy = permute_bbox2.ymin - permute_bbox1.ymax
        if dy < 0:
            dy = 0

        if permute_bbox2.zmin < permute_bbox1.zmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dz = permute_bbox2.zmin - permute_bbox1.zmax
        if dz < 0:
            dz = 0

        return (dx ** 2 + dy ** 2 + dz ** 2) ** 0.5

    def point_belongs(self, point):
        return self.xmin < point[0] and point[0] < self.xmax \
               and self.ymin < point[1] and point[1] < self.ymax \
               and self.zmin < point[2] and point[2] < self.zmax

    def distance_to_point(self, point):
        if self.point_belongs(point):
            return min([self.xmax - point[0], point[0] - self.xmin,
                        self.ymax - point[1], point[1] - self.ymin,
                        self.zmax - point[2], point[2] - self.zmin])
        else:
            if point[0] < self.xmin:
                dx = self.xmin - point[0]
            elif self.xmax < point[0]:
                dx = point[0] - self.xmax
            else:
                dx = 0

            if point[1] < self.ymin:
                dy = self.ymin - point[1]
            elif self.ymax < point[1]:
                dy = point[1] - self.ymax
            else:
                dy = 0

            if point[2] < self.zmin:
                dz = self.zmin - point[2]
            elif self.zmax < point[2]:
                dz = point[2] - self.zmax
            else:
                dz = 0
        return (dx ** 2 + dy ** 2 + dz ** 2) ** 0.5

    def distance_between_two_points_on_bbox(self, point1, point2):

        if math.isclose(point1[0], self.xmin, abs_tol=1e-8):
            face_point1 = 5
        elif math.isclose(point1[0], self.xmax, abs_tol=1e-8):
            face_point1 = 3
        elif math.isclose(point1[1], self.ymin, abs_tol=1e-8):
            face_point1 = 4
        elif math.isclose(point1[1], self.ymax, abs_tol=1e-8):
            face_point1 = 2
        elif math.isclose(point1[2], self.zmin, abs_tol=1e-8):
            face_point1 = 6
        elif math.isclose(point1[2], self.zmax, abs_tol=1e-8):
            face_point1 = 1
        else:
            raise NotImplementedError

        if math.isclose(point2[0], self.xmin, abs_tol=1e-8):
            face_point2 = 5
        elif math.isclose(point2[0], self.xmax, abs_tol=1e-8):
            face_point2 = 3
        elif math.isclose(point2[1], self.ymin, abs_tol=1e-8):
            face_point2 = 4
        elif math.isclose(point2[1], self.ymax, abs_tol=1e-8):
            face_point2 = 2
        elif math.isclose(point2[2], self.zmin, abs_tol=1e-8):
            face_point2 = 6
        elif math.isclose(point2[2], self.zmax, abs_tol=1e-8):
            face_point2 = 1
        else:
            raise NotImplementedError

        point1_copy = point1.copy()
        point2_copy = point2.copy()
        if face_point1 > face_point2:
            point1, point2 = point2, point1
            face_point1, face_point2 = face_point2, face_point1

        # The points are on the same face
        if face_point1 == face_point2:
            return point1.point_distance(point2)

        deltax = self.xmax - self.xmin
        deltay = self.ymax - self.ymin
        deltaz = self.zmax - self.zmin

        point1_2d_coordinate_dict = {1: Point2D((point1[
                                                     0] - self.xmin - deltax / 2,
                                                 point1[
                                                     1] - self.ymin - deltay / 2)),
                                     2: Point2D((point1[
                                                     2] - self.zmin - deltaz / 2,
                                                 point1[
                                                     0] - self.xmin - deltax / 2)),
                                     3: Point2D((point1[
                                                     1] - self.ymin - deltay / 2,
                                                 point1[
                                                     2] - self.zmin - deltaz / 2)),
                                     4: Point2D((point1[
                                                     0] - self.xmin - deltax / 2,
                                                 point1[
                                                     2] - self.zmin - deltaz / 2)),
                                     5: Point2D((point1[
                                                     2] - self.zmin - deltaz / 2,
                                                 point1[
                                                     1] - self.ymin - deltay / 2)),
                                     6: Point2D((point1[
                                                     1] - self.ymin - deltay / 2,
                                                 point1[
                                                     0] - self.xmin - deltax / 2))}

        point2_2d_coordinate_dict = {1: Point2D((point2[
                                                     0] - self.xmin - deltax / 2,
                                                 point2[
                                                     1] - self.ymin - deltay / 2)),
                                     2: Point2D((point2[
                                                     2] - self.zmin - deltaz / 2,
                                                 point2[
                                                     0] - self.xmin - deltax / 2)),
                                     3: Point2D((point2[
                                                     1] - self.ymin - deltay / 2,
                                                 point2[
                                                     2] - self.zmin - deltaz / 2)),
                                     4: Point2D((point2[
                                                     0] - self.xmin - deltax / 2,
                                                 point2[
                                                     2] - self.zmin - deltaz / 2)),
                                     5: Point2D((point2[
                                                     2] - self.zmin - deltaz / 2,
                                                 point2[
                                                     1] - self.ymin - deltay / 2)),
                                     6: Point2D((point2[
                                                     1] - self.ymin - deltay / 2,
                                                 point2[
                                                     0] - self.xmin - deltax / 2))}

        vertex_2d_coordinate_dict = {1: [Point2D((
                                                 self.xmin - self.xmin - deltax / 2,
                                                 self.ymin - self.ymin - deltay / 2)),
                                         Point2D((
                                                 self.xmin - self.xmin - deltax / 2,
                                                 self.ymax - self.ymin - deltay / 2)),
                                         Point2D((
                                                 self.xmax - self.xmin - deltax / 2,
                                                 self.ymax - self.ymin - deltay / 2)),
                                         Point2D((
                                                 self.xmax - self.xmin - deltax / 2,
                                                 self.ymin - self.ymin - deltay / 2))],
                                     2: [Point2D((
                                                 self.zmin - self.zmin - deltaz / 2,
                                                 self.xmin - self.xmin - deltax / 2)),
                                         Point2D((
                                                 self.zmin - self.zmin - deltaz / 2,
                                                 self.xmax - self.xmin - deltax / 2)),
                                         Point2D((
                                                 self.zmax - self.zmin - deltaz / 2,
                                                 self.xmax - self.xmin - deltax / 2)),
                                         Point2D((
                                                 self.zmax - self.zmin - deltaz / 2,
                                                 self.xmin - self.xmin - deltax / 2))],
                                     3: [Point2D((
                                                 self.ymin - self.ymin - deltay / 2,
                                                 self.zmin - self.zmin - deltaz / 2)),
                                         Point2D((
                                                 self.ymin - self.ymin - deltay / 2,
                                                 self.zmax - self.zmin - deltaz / 2)),
                                         Point2D((
                                                 self.ymax - self.ymin - deltay / 2,
                                                 self.zmax - self.zmin - deltaz / 2)),
                                         Point2D((
                                                 self.ymax - self.ymin - deltay / 2,
                                                 self.zmin - self.zmin - deltaz / 2))],
                                     4: [Point2D((
                                                 self.xmin - self.xmin - deltax / 2,
                                                 self.zmin - self.zmin - deltaz / 2)),
                                         Point2D((
                                                 self.xmin - self.xmin - deltax / 2,
                                                 self.zmax - self.zmin - deltaz / 2)),
                                         Point2D((
                                                 self.xmax - self.xmin - deltax / 2,
                                                 self.zmax - self.zmin - deltaz / 2)),
                                         Point2D((
                                                 self.xmax - self.xmin - deltax / 2,
                                                 self.zmin - self.zmin - deltaz / 2))],
                                     5: [Point2D((
                                                 self.zmin - self.zmin - deltaz / 2,
                                                 self.ymin - self.ymin - deltay / 2)),
                                         Point2D((
                                                 self.zmin - self.zmin - deltaz / 2,
                                                 self.ymax - self.ymin - deltay / 2)),
                                         Point2D((
                                                 self.zmax - self.zmin - deltaz / 2,
                                                 self.ymax - self.ymin - deltay / 2)),
                                         Point2D((
                                                 self.zmax - self.zmin - deltaz / 2,
                                                 self.ymin - self.ymin - deltay / 2))],
                                     6: [Point2D((
                                                 self.ymin - self.ymin - deltay / 2,
                                                 self.xmin - self.xmin - deltax / 2)),
                                         Point2D((
                                                 self.ymin - self.ymin - deltay / 2,
                                                 self.xmax - self.xmin - deltax / 2)),
                                         Point2D((
                                                 self.ymax - self.ymin - deltay / 2,
                                                 self.xmax - self.xmin - deltax / 2)),
                                         Point2D((
                                                 self.ymax - self.ymin - deltay / 2,
                                                 self.xmin - self.xmin - deltax / 2))], }

        vertex_to_3d_dict = {1: (2, self.zmax, 0, 1),
                             2: (1, self.ymax, 2, 0),
                             3: (0, self.xmax, 1, 2),
                             4: (1, self.ymin, 0, 2),
                             5: (0, self.xmin, 2, 1),
                             6: (2, self.zmin, 1, 0)}

        offset_dict = {0: self.xmin + deltax / 2,
                       1: self.ymin + deltay / 2,
                       2: self.zmin + deltaz / 2}

        opposite_face_dict = {1: 6, 2: 4, 3: 5, 4: 2, 5: 3, 6: 1}

        combination_dict = {
            (1, 2): Frame2D(Point2D((0, deltay / 2 + deltaz / 2)),
                            Vector2D((0, -1)), Vector2D((1, 0))),
            (2, 1): Frame2D(Point2D((deltay / 2 + deltaz / 2, 0)),
                            Vector2D((0, 1)), Vector2D((-1, 0))),
            (1, 3): Frame2D(Point2D((deltax / 2 + deltaz / 2, 0)),
                            Vector2D((0, 1)), Vector2D((-1, 0))),
            (3, 1): Frame2D(Point2D((0, deltax / 2 + deltaz / 2)),
                            Vector2D((0, -1)), Vector2D((1, 0))),
            (1, 4): Frame2D(Point2D((0, -deltay / 2 - deltaz / 2)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (4, 1): Frame2D(Point2D((-deltay / 2 - deltaz / 2, 0)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (1, 5): Frame2D(Point2D((-deltax / 2 - deltaz / 2, 0)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (5, 1): Frame2D(Point2D((0, -deltax / 2 - deltaz / 2)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (2, 3): Frame2D(Point2D((0, deltax / 2 + deltay / 2)),
                            Vector2D((0, -1)), Vector2D((1, 0))),
            (3, 2): Frame2D(Point2D((deltax / 2 + deltay / 2, 0)),
                            Vector2D((0, 1)), Vector2D((-1, 0))),
            (2, 5): Frame2D(Point2D((0, -deltax / 2 - deltay / 2)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (5, 2): Frame2D(Point2D((-deltax / 2 - deltay / 2, 0)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (2, 6): Frame2D(Point2D((-deltaz / 2 - deltay / 2, 0)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (6, 2): Frame2D(Point2D((0, -deltaz / 2 - deltay / 2)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (3, 4): Frame2D(Point2D((-deltay / 2 - deltax / 2, 0)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (4, 3): Frame2D(Point2D((0, -deltay / 2 - deltax / 2)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (3, 6): Frame2D(Point2D((0, -deltaz / 2 - deltax / 2)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (6, 3): Frame2D(Point2D((-deltaz / 2 - deltax / 2, 0)),
                            Vector2D((1, 0)), Vector2D((0, 1))),
            (4, 5): Frame2D(Point2D((-deltax / 2 - deltay / 2, 0)),
                            Vector2D((0, 1)), Vector2D((-1, 0))),
            (5, 4): Frame2D(Point2D((0, -deltax / 2 - deltay / 2)),
                            Vector2D((0, -1)), Vector2D((1, 0))),
            (4, 6): Frame2D(Point2D((0, -deltaz / 2 - deltay / 2)),
                            Vector2D((0, -1)), Vector2D((1, 0))),
            (6, 4): Frame2D(Point2D((-deltaz / 2 - deltay / 2, 0)),
                            Vector2D((0, 1)), Vector2D((-1, 0))),
            (5, 6): Frame2D(Point2D((-deltaz / 2 - deltax / 2, 0)),
                            Vector2D((0, 1)), Vector2D((-1, 0))),
            (6, 5): Frame2D(Point2D((0, -deltaz / 2 - deltax / 2)),
                            Vector2D((0, -1)), Vector2D((1, 0)))}

        point1_2d = point1_2d_coordinate_dict[face_point1]
        point2_2d = point2_2d_coordinate_dict[face_point2]

        # The points are on adjacent faces
        if opposite_face_dict[face_point1] != face_point2:
            frame = combination_dict[(face_point1, face_point2)]
            net_point2 = frame.OldCoordinates(point2_2d)

            # Computes the 3D intersection between the net_line and the edges of the face_point1
            net_line = LineSegment2D(point1_2d, net_point2)
            vertex_points = vertex_2d_coordinate_dict[face_point1]
            edge_lines = [LineSegment2D(p1, p2) for p1, p2 in
                          zip(vertex_points,
                              vertex_points[1:] + [vertex_points[0]])]
            for line in edge_lines:
                edge_intersection_point, a, b = Point2D.LinesIntersection(
                    net_line, line, curvilinear_abscissa=True)
                if edge_intersection_point is not None \
                        and a > 0 and a < 1 and b > 0 and b < 1:
                    break
            offset_indice, offset, indice1, indice2 = vertex_to_3d_dict[
                face_point1]
            disordered_coordinate = [
                (indice1, edge_intersection_point[0] + offset_dict[indice1]),
                (indice2, edge_intersection_point[1] + offset_dict[indice2]),
                (offset_indice, offset)]
            disordered_coordinate = sorted(disordered_coordinate,
                                           key=lambda a: a[0])
            intersection_point_3d = Point3D(
                tuple([p[1] for p in disordered_coordinate]))

            mesures = [Measure3D(point1_copy, intersection_point_3d),
                       Measure3D(intersection_point_3d, point2_copy)]

            return mesures

        # The points are on opposite faces
        else:
            net_points2_and_frame = []

            faces_number = [1, 2, 3, 4, 5, 6]
            faces_number.remove(face_point1)
            faces_number.remove(face_point2)
            pathes = []
            for face_nb in faces_number:
                path = [(face_point1, face_nb), (face_nb, face_point2)]
                pathes.append(path)

            for path in pathes:
                frame1 = combination_dict[(path[0][0], path[0][1])]
                frame2 = combination_dict[(path[1][0], path[1][1])]
                frame = frame1 + frame2
                net_points2_and_frame.append(
                    (Point2D(frame.OldCoordinates(point2_2d).vector), frame))
            net_point2, frame = min(net_points2_and_frame,
                                    key=lambda pt: pt[0].point_distance(
                                        point1_2d))
            net_line = LineSegment2D(point1_2d, net_point2)

            # Computes the 3D intersection between the net_line and the edges of the face_point1
            vertex_points = vertex_2d_coordinate_dict[face_point1]
            edge_lines = [LineSegment2D(p1, p2) for p1, p2 in
                          zip(vertex_points,
                              vertex_points[1:] + [vertex_points[0]])]
            for line in edge_lines:
                edge_intersection_point1, a, b = Point2D.LinesIntersection(
                    net_line, line, curvilinear_abscissa=True)
                if edge_intersection_point1 is not None \
                        and a > 0 and a < 1 and b > 0 and b < 1:
                    break
            offset_indice, offset, indice1, indice2 = vertex_to_3d_dict[
                face_point1]
            disordered_coordinate = [
                (indice1, edge_intersection_point1[0] + offset_dict[indice1]),
                (indice2, edge_intersection_point1[1] + offset_dict[indice2]),
                (offset_indice, offset)]
            disordered_coordinate = sorted(disordered_coordinate,
                                           key=lambda a: a[0])
            intersection_point1_3d = Point3D(
                tuple([p[1] for p in disordered_coordinate]))

            # Computes the 3D intersection between the net_line and the edges of the face_point2
            vertex_points = [frame.OldCoordinates(p) for p in
                             vertex_2d_coordinate_dict[face_point2]]
            edge_lines = [LineSegment2D(p1, p2) for p1, p2 in
                          zip(vertex_points,
                              vertex_points[1:] + [vertex_points[0]])]
            for line in edge_lines:
                edge_intersection_point2, a, b = Point2D.LinesIntersection(
                    net_line, line, curvilinear_abscissa=True)
                if edge_intersection_point2 is not None \
                        and a > 0 and a < 1 and b > 0 and b < 1:
                    break
            edge_intersection_point2 = Point2D(
                frame.NewCoordinates(edge_intersection_point2))
            offset_indice, offset, indice1, indice2 = vertex_to_3d_dict[
                face_point2]
            disordered_coordinate = [
                (indice1, edge_intersection_point2[0] + offset_dict[indice1]),
                (indice2, edge_intersection_point2[1] + offset_dict[indice2]),
                (offset_indice, offset)]
            disordered_coordinate = sorted(disordered_coordinate,
                                           key=lambda a: a[0])
            intersection_point2_3d = Point3D(
                tuple([p[1] for p in disordered_coordinate]))

            if point1 == point1_copy:
                mesures = [Measure3D(point1, intersection_point1_3d),
                           Measure3D(intersection_point1_3d,
                                     intersection_point2_3d),
                           Measure3D(intersection_point2_3d, point2)]
            else:
                mesures = [Measure3D(point2, intersection_point2_3d),
                           Measure3D(intersection_point2_3d,
                                     intersection_point1_3d),
                           Measure3D(intersection_point1_3d, point1)]
            return mesures

    def babylon_script(self):
        height = self.ymax - self.ymin
        width = self.xmax - self.xmin
        depth = self.zmax - self.zmin
        s = 'var box = BABYLON.MeshBuilder.CreateBox("box", {{height: {}, width: {}, depth: {}}}, scene);\n'.format(
            height, width, depth)
        s += 'box.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n'.format(
            self.center[0], self.center[1], self.center[2])
        s += 'var bboxmat = new BABYLON.StandardMaterial("bboxmat", scene);\n'
        s += 'bboxmat.alpha = 0.4;\n'
        s += 'var DTWidth = {};\n'.format(width * 60)
        s += 'var DTHeight = {};\n'.format(height * 60)
        s += 'var font_type = "Arial";\n'
        s += 'var text = "{}";\n'.format(self.name)
        s += 'var dynamicTexture = new BABYLON.DynamicTexture("DynamicTexture", {width:DTWidth, height:DTHeight}, scene);\n'
        s += 'var ctx = dynamicTexture.getContext();\n'
        s += 'var size = 0.8;\n'
        s += 'ctx.font = size + "px " + font_type;\n'
        s += 'var textWidth = ctx.measureText(text).width;\n'
        s += 'var ratio = textWidth/size;\n'
        s += 'var font_size = Math.floor(DTWidth / ratio);\n'
        s += 'var font = font_size + "px " + font_type;\n'
        s += 'dynamicTexture.drawText(text, null, null, font, "#000000", "#ffffff", false);\n'
        s += 'bboxmat.diffuseTexture = dynamicTexture;\n'
        s += 'box.material = bboxmat;\n'
        return s




class VolumeModel(dc.DessiaObject):
    _standalone_in_db = True
    _generic_eq = True
    _non_serializable_attributes = ['shells', 'bounding_box']
    _non_eq_attributes = ['name', 'shells', 'bounding_box', 'contours',
                          'faces']
    _non_hash_attributes = []
    """
    :param groups: A list of two element tuple. The first element is a string naming the group and the second element is a list of primitives of the group
    """

    def __init__(self, primitives, name=''):
        self.primitives = primitives
        self.name = name
        self.shells = []
        if self.primitives:
            self.shells = self._extract_shells()
        if self.shells:
            self.bounding_box = self._bounding_box()
        else:
            self.bounding_box = BoundingBox(-1, 1, -1, 1, -1, 1)

    # def __hash__(self):
    #     return sum([hash(p) for p in self.primitives])

    def _extract_shells(self):
        shells = []
        for primitive in self.primitives:
            if isinstance(primitive, Shell3D):
                shells.append(primitive)
        return shells

    def Volume(self):
        volume = 0
        for primitive in self.primitives:
            volume += primitive.Volume()
        return volume

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_primitives = [
                primitive.Rotation(center, axis, angle, copy=True) for
                primitive in self.primitives]
            return VolumeModel(new_primitives, self.name)
        else:
            for primitives in self.primitives:
                primitives.Translation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def Translation(self, offset, copy=True):
        if copy:
            new_primitives = [primitive.Translation(offset, copy=True) for
                              primitive in self.primitives]
            return VolumeModel(new_primitives, self.name)
        else:
            for primitives in self.primitives:
                primitives.Translation(offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_primitives = [primitive.frame_mapping(frame, side, copy=True)
                              for primitive in self.primitives]
            return VolumeModel(new_primitives, self.name)
        else:
            for primitives in self.primitives:
                primitives.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()

    def copy(self):
        new_primitives = [primitive.copy() for primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def _bounding_box(self):
        bboxes = []
        points = []
        for primitive in self.primitives:
            if hasattr(primitive, 'bounding_box'):
                bboxes.append(primitive.bounding_box)
            else:
                if primitive.__class__.__name__ == 'Point3D':
                    points.append(primitive)
        if bboxes:
            xmin = min([box.xmin for box in bboxes])
            xmax = max([box.xmax for box in bboxes])
            ymin = min([box.ymin for box in bboxes])
            ymax = max([box.ymax for box in bboxes])
            zmin = min([box.zmin for box in bboxes])
            zmax = max([box.zmax for box in bboxes])
        elif points:
            xmin = min([p[0] for p in points])
            xmax = max([p[0] for p in points])
            ymin = min([p[1] for p in points])
            ymax = max([p[1] for p in points])
            zmin = min([p[2] for p in points])
            zmax = max([p[2] for p in points])
        else:
            # raise ValueError('Bounding box cant be determined')
            return BoundingBox(-1, 1, -1, 1, 1 - 1, 1)
        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def plot(self, ax=None, color=None):
        fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')

        for i, shell in enumerate(self.shells):
            bbox = shell.bbox()
            bbox.plot(ax, color[i])

        return ax

    def MPLPlot(self):
        """
        Matplotlib plot of model.
        To use for debug.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', adjustable='box')
        #        ax.set_aspect('equal')
        for primitive in self.primitives:
            primitive.MPLPlot(ax)
        ax.set_aspect('equal')
        ax.margins(0.1)
        return ax

    def FreeCADScript(self, fcstd_filepath,
                      freecad_lib_path='/usr/lib/freecad/lib',
                      export_types=('fcstd',),
                      save_to='',
                      tolerance=0.0001):
        """
        Generate python a FreeCAD definition of model
        :param fcstd_filename: a filename without extension to give the name at the fcstd part written in python code
        :type fcstd_filename:str
        """
        fcstd_filepath = os.path.abspath(fcstd_filepath)
        fcstd_filepath = fcstd_filepath.replace('\\', '\\\\')
        freecad_lib_path = freecad_lib_path.replace('\\', '\\\\')

        s = '# -*- coding: utf-8 -*-\n'
        if freecad_lib_path != '':
            s += "import sys\nsys.path.append('" + freecad_lib_path + "')\n"

        s += "import math\nimport FreeCAD as fc\nimport Part\n\ndoc=fc.newDocument('doc')\n\n"
        for ip, primitive in enumerate(self.primitives):
            if primitive.name == '':
                primitive_name = 'Primitive_{}'.format(ip)
            else:
                primitive_name = 'Primitive_{}_{}'.format(ip, primitive.name)
            s += "part = doc.addObject('App::Part','{}')\n".format(
                primitive_name)
            if hasattr(primitive, 'FreeCADExport'):
                sp = primitive.FreeCADExport(ip)
                if sp != '':
                    #                        s += (sp+'\n')
                    s += (sp)
                    s += 'shapeobj = doc.addObject("Part::Feature","{}")\n'.format(
                        primitive_name)
                    if isinstance(primitive, BSplineCurve3D) \
                            or isinstance(primitive, BSplineSurface3D) \
                            or isinstance(primitive, Circle3D) \
                            or isinstance(primitive, LineSegment3D) \
                            or isinstance(primitive, Ellipse3D):
                        #                            print(primitive)
                        #                            s += 'S = Part.Shape([primitive{}])\n'.format(ip)
                        #                            s += 'shapeobj.Shape = S\n'
                        s += 'shapeobj.Shape = primitive{}.toShape()\n'.format(
                            ip)
                    else:
                        s += "shapeobj.Shape = primitive{}\n".format(ip)
                    s += 'part.addObject(shapeobj)\n\n'.format(ip,
                                                               primitive.name)
            # --------------------DEBUG-------------------
        #                else:
        #                    raise NotImplementedError
        # ---------------------------------------------

        s += 'doc.recompute()\n'
        if 'fcstd' in export_types:
            s += "doc.saveAs('" + fcstd_filepath + ".fcstd')\n\n"
        if 'stl' in export_types:
            s += "import Mesh\nMesh.export(doc.Objects,'{}.stl', tolerance={})\n".format(
                fcstd_filepath, tolerance)
        if 'step' in export_types:
            s += "Part.export(doc.Objects,'{}.step')\n".format(fcstd_filepath)

        if save_to != '':
            with open(os.path.abspath(save_to), 'w') as file:
                file.write(s)
        return s

    def FreeCADExport(self, fcstd_filepath,
                      python_path='python',
                      freecad_lib_path='/usr/lib/freecad/lib',
                      export_types=('fcstd',),
                      tolerance=0.0001):
        """
        Export model to .fcstd FreeCAD standard

        :param python_path: path of python binded to freecad

            * on windows: something like C:\\\\Program Files\\\\FreeCAD X.XX\\\\bin\\\\python
            * on linux: python if installed by a dstribution package
        :param filepath: path of fcstd file (without extension)
        :param freecad_lib_path: FreeCAD.so lib path (/usr/lib/freecad/lib in general)
        :param tolerance: the tolerance of tesselation for mesh exports

        """
        fcstd_filepath = os.path.abspath(fcstd_filepath)
        s = self.FreeCADScript(fcstd_filepath,
                               freecad_lib_path=freecad_lib_path,
                               export_types=export_types,
                               tolerance=tolerance)
        with tempfile.NamedTemporaryFile(suffix=".py", delete=False) as f:
            f.write(bytes(s, 'utf8'))

        arg = f.name
        output = subprocess.call([python_path, arg])

        f.close()
        os.remove(f.name)
        return output

    def babylon_script(self, use_cdn=True, debug=False):

        env = Environment(loader=PackageLoader('volmdlr', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))

        template = env.get_template('babylon.html')

        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])

        primitives_strings = []
        for primitive in self.primitives:
            if hasattr(primitive, 'babylon_script'):
                primitives_strings.append(primitive.babylon_script())

        return template.render(name=self.name,
                               center=tuple(center),
                               length=2 * max_length,
                               primitives_strings=primitives_strings,
                               use_cdn=use_cdn,
                               debug=debug)

    def babylonjs_from_script(self, page_name=None, use_cdn=True, debug=False):
        script = self.babylon_script(use_cdn=use_cdn, debug=debug)

        if page_name is None:
            with tempfile.NamedTemporaryFile(suffix=".html",
                                             delete=False) as file:
                file.write(bytes(script, 'utf8'))
            page_name = file.name
        else:
            page_name += '.html'
            with open(page_name, 'w')  as file:
                file.write(script)

        webbrowser.open('file://' + os.path.realpath(page_name))

    def babylon_data(self):
        meshes = []
        for primitive in self.primitives:
            if hasattr(primitive, 'babylon_meshes'):
                meshes.extend(primitive.babylon_meshes())
        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])

        babylon_data = {'meshes': meshes,
                        'max_length': max_length,
                        'center': list(center)}
        return babylon_data

    @classmethod
    def babylonjs_from_babylon_data(cls, babylon_data, page_name=None,
                                    use_cdn=True, debug=False):
        env = Environment(loader=PackageLoader('volmdlr', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))

        template = env.get_template('babylon_unpacker.html')

        script = template.render(babylon_data=babylon_data,
                                 use_cdn=use_cdn,
                                 debug=debug
                                 )
        if page_name is None:
            with tempfile.NamedTemporaryFile(suffix=".html",
                                             delete=False) as file:
                file.write(bytes(script, 'utf8'))
            page_name = file.name
        else:
            page_name += '.html'
            with open(page_name, 'w')  as file:
                file.write(script)

        webbrowser.open('file://' + os.path.realpath(page_name))

    def babylonjs(self, page_name=None, use_cdn=True, debug=False):
        # print('self.primitives', self.primitives)
        babylon_data = self.babylon_data()
        self.babylonjs_from_babylon_data(babylon_data, page_name=page_name,
                                         use_cdn=use_cdn, debug=debug)


class MovingVolumeModel(VolumeModel):
    def __init__(self, primitives, step_frames, name=''):
        VolumeModel.__init__(self, primitives=primitives, name=name)
        self.step_frames = step_frames

        if not self.is_consistent():
            raise dc.ConsistencyError

    def is_consistent(self):
        n_primitives = len(self.primitives)
        for frames in self.step_frames:
            if len(frames) != n_primitives:
                return False
        return True

    def step_volume_model(self, istep):
        primitives = []
        for primitive, frame in zip(self.primitives, self.step_frames[istep]):
            primitives.append(
                primitive.frame_mapping(frame, side='old', copy=True))
        return VolumeModel(primitives)

    def babylon_script(self, use_cdn=True, debug=False):

        env = Environment(loader=PackageLoader('volmdlr', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))

        template = env.get_template('babylon.html')

        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])

        primitives_strings = []
        for primitive in self.primitives:
            if hasattr(primitive, 'babylon_script'):
                primitives_strings.append(primitive.babylon_script())

        positions = []
        orientations = []
        for step in self.step_frames:
            step_positions = []
            step_orientations = []

            for frame in step:
                step_positions.append(list(frame.origin))
                step_orientations.append([list(frame.u),
                                          list(frame.v),
                                          list(frame.w)])

            positions.append(step_positions)
            orientations.append(step_orientations)

        return template.render(name=self.name,
                               center=tuple(center),
                               length=2 * max_length,
                               primitives_strings=primitives_strings,
                               positions=positions,
                               orientations=orientations,
                               use_cdn=use_cdn,
                               debug=debug)

    def babylon_data(self):
        meshes = []
        primitives_to_meshes = []
        for ip, primitive in enumerate(self.primitives):
            if hasattr(primitive, 'babylon_meshes'):
                meshes.extend(primitive.babylon_meshes())
                primitives_to_meshes.append(ip)

        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])

        steps = []
        for istep, frames in enumerate(self.step_frames):

            # step_positions = []
            # step_orientations = []
            step = {'time': istep}
            for iframe, frame in enumerate(frames):
                if iframe in primitives_to_meshes:
                    imesh = primitives_to_meshes.index(iframe)
                    step[imesh] = {}
                    step[imesh]['position'] = list(round(frame.origin, 6))
                    step[imesh]['orientations'] = [list(round(frame.u, 6)),
                                                   list(round(frame.v, 6)),
                                                   list(round(frame.w, 6))]

            steps.append(step)

        babylon_data = {'meshes': meshes,
                        'max_length': max_length,
                        'center': list(center),
                        'steps': steps}
        return babylon_data


class Routing:
    def __init__(self, point1, point2, volumemodel):
        self.points = [point1, point2]
        self.volumemodel = volumemodel

    def straight_line(self):
        """
        Returns 2 distances :
            - no collision distance
            - collision distance
        """
        line = LineSegment3D(self.points[0], self.points[1])

        intersection_points = []
        abscissea_list = []
        for shell in self.volumemodel.shells:
            for face in shell.faces:
                intersection_point, intersection_abscissea = face.linesegment_intersection(
                    line, abscissea=True)
                if intersection_point is not None and intersection_abscissea != 0 and intersection_abscissea != 1:
                    not_in_abscissea_list = True
                    for abscissea in abscissea_list:
                        if math.isclose(abscissea, intersection_abscissea,
                                        abs_tol=1e-8):
                            not_in_abscissea_list = False
                    if not_in_abscissea_list:
                        intersection_points.append(
                            (intersection_point, intersection_abscissea))
                        abscissea_list.append(intersection_abscissea)

        if len(intersection_points) % 2 != 0:
            raise NotImplementedError

        intersection_points = sorted(intersection_points,
                                     key=lambda abscissea: abscissea[1])
        all_points_abscissea = [(self.points[0], 0)] + intersection_points[
                                                       :] + [
                                   (self.points[1], 1)]
        all_points = [p[0] for p in all_points_abscissea]

        no_collision_mesures = []
        collision_mesures = []
        i = 0
        for pt1, pt2 in zip(all_points[:-1], all_points[1:]):
            if i % 2 == 0:
                no_collision_mesures.append(
                    Measure3D(pt1, pt2, color=(0, 0, 1)))
            else:
                collision_mesures.append(Measure3D(pt1, pt2, color=(1, 0, 0)))
            i += 1

        return no_collision_mesures, collision_mesures

    def straight_line2(self):
        """
        Returns the distance of the line going around each shell's bbox encountered along the path
        """
        line = LineSegment3D(self.points[0], self.points[1])

        all_mesures_abscissea = []
        intersection_points = []
        for shell in self.volumemodel.shells:
            shell_intersection_points = []
            bbox = shell.bounding_box
            for face in shell.faces:
                intersection_point, intersection_abscissea = face.linesegment_intersection(
                    line, abscissea=True)
                if intersection_point is not None and intersection_abscissea != 0 and intersection_abscissea != 1:
                    intersection_points.append(
                        (intersection_point, intersection_abscissea))
                    shell_intersection_points.append(
                        (intersection_point, intersection_abscissea))

            if len(shell_intersection_points) == 2:
                shell_intersection_points = sorted(shell_intersection_points,
                                                   key=lambda abscissea:
                                                   abscissea[1])
                abscissea1 = shell_intersection_points[0][1]
                abscissea2 = shell_intersection_points[1][1]
                shell_intersection_points = [p[0] for p in
                                             shell_intersection_points]
                around_bbox_mesures = bbox.distance_between_two_points_on_bbox(
                    shell_intersection_points[0], shell_intersection_points[1])
                all_mesures_abscissea.append(
                    (around_bbox_mesures, abscissea1, abscissea2))
            elif len(shell_intersection_points) > 2 or len(
                    shell_intersection_points) == 1:
                raise NotImplementedError

        intersection_points = sorted(intersection_points,
                                     key=lambda abscissea: abscissea[1])
        all_mesures_abscissea = sorted(all_mesures_abscissea,
                                       key=lambda abscissea: abscissea[1])
        all_points_abscissea = [(self.points[0], 0)] + intersection_points[
                                                       :] + [
                                   (self.points[1], 1)]
        all_points = [p[0] for p in all_points_abscissea]

        no_collision_mesures = []
        i = 0
        for pt1, pt2 in zip(all_points[:-1], all_points[1:]):
            if i % 2 == 0:
                no_collision_mesures.append(
                    Measure3D(pt1, pt2, color=(0, 0, 1)))
            else:
                no_collision_mesures.extend(all_mesures_abscissea[i // 2][0])
            i += 1

        return no_collision_mesures


class ViewIso:  # TODO: rename this in IsoView
    def __init__(self, component, frame, size):
        self.component = component
        self.frame = frame
        self.size = size
        self.plot_datas = self.plot_data()

    def plot_data(self, detail=True):
        wide = min(self.size) / 2
        plot_datas = []
        plot_datas.extend(self.component.plot_data(self.frame, detail=detail))
        plot_datas.extend(self.component.plot_data(Frame3D(
            self.frame.origin + Point3D(
                (0, self.size[1] / 2 + self.size[2] / 2 + wide, 0)),
            self.frame.u, self.frame.w, self.frame.v), detail=detail))
        plot_datas.extend(self.component.plot_data(Frame3D(
            self.frame.origin + Point3D(
                (self.size[0] / 2 + self.size[2] / 2 + wide, 0, 0)),
            self.frame.w, self.frame.v, self.frame.u), detail=detail))
        return plot_datas

    def plot(self):
        plot_data.plot(self.plot_datas)



