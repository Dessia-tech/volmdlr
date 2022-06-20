#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
import os
import subprocess
import tempfile
import webbrowser
from datetime import datetime
# import volmdlr.stl as vmstl
from typing import List

import matplotlib.pyplot as plt
import numpy as npy

import dessia_common as dc
import dessia_common.files as dcf
import volmdlr
import volmdlr.templates

# from mpl_toolkits.mplot3d import Axes3D

npy.seterr(divide='raise')


# TODO: put voldmlr metadata in this freecad header
STEP_HEADER = '''ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('{name}'),'2;1');
FILE_NAME('{filename}','{timestamp}',('Author'),(''),'Volmdlr v{version}','','Unknown');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN {{ 1 0 10303 214 1 1 1 1 }}'));
ENDSEC;
DATA;
#1 = APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#2);
#2 = APPLICATION_CONTEXT('core data for automotive mechanical design processes');
#3 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );
#4 = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );
#5 = ( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() );
#6 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#3,'distance_accuracy_value','confusion accuracy');
#7 = ( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#6)) GLOBAL_UNIT_ASSIGNED_CONTEXT ((#3,#4,#5)) REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY') );
'''

STEP_FOOTER = '''ENDSEC;
END-ISO-10303-21;
'''


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
    return list(char_list)


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
    vector0 = volmdlr.O2D
    if vector0 in (vector1, vector2):
        return 0

    dot = vector1.dot(vector2)
    norm_vec_1 = vector1.norm()
    norm_vec_2 = vector2.norm()
    sol = dot / (norm_vec_1 * norm_vec_2)
    cross = vector1.cross(vector2)
    if math.isclose(sol, 1, abs_tol=1e-6):
        inner_angle = 0
    elif math.isclose(sol, -1, abs_tol=1e-6):
        inner_angle = math.pi
    else:
        inner_angle = math.acos(sol)

    if cross < 0:
        return inner_angle

    return volmdlr.TWO_PI - inner_angle


def vectors3d_angle(vector1, vector2):
    dot = vector1.dot(vector2)
    theta = math.acos(dot / (vector1.norm() * vector2.norm()))

    return theta


def sin_cos_angle(u1, u2):
    """
    cos(theta)=u1, sin(theta)=u2
    returns an angle between 0 and 2pi
    """
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
            theta = volmdlr.TWO_PI + math.asin(u2)
    else:
        if u2 >= 0:
            theta = math.acos(u1)
        else:
            theta = volmdlr.TWO_PI - math.acos(u1)
    if math.isclose(theta, volmdlr.TWO_PI, abs_tol=1e-9):
        return 0.
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
        new_face_triangles = []
        for triangle_ in face_triangles:
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
            new_pt = all_points[pos].copy() + volmdlr.Point2D((volmdlr.TWO_PI, 0))
            if new_pt.vector[0] > volmdlr.TWO_PI:
                new_pt.vector[0] = volmdlr.TWO_PI
            all_points[pos] = new_pt
    elif len(plus_pi) <= 2 and len(all_points) > 4:
        for pos in plus_pi:
            new_pt = all_points[pos].copy() - volmdlr.Point2D((volmdlr.TWO_PI, 0))
            if new_pt.vector[0] < 0:
                new_pt.vector[0] = 0
            all_points[pos] = new_pt
    if 3 * len(moins_pi) <= len(plus_pi) and len(all_points) > 4:
        for pos in moins_pi:
            new_pt = all_points[pos].copy() + volmdlr.Point2D((volmdlr.TWO_PI, 0))
            if new_pt.vector[0] > volmdlr.TWO_PI:
                new_pt.vector[0] = volmdlr.TWO_PI
            all_points[pos] = new_pt
    elif 3 * len(plus_pi) <= len(moins_pi) and len(all_points) > 4:
        for pos in plus_pi:
            new_pt = all_points[pos].copy() - volmdlr.Point2D((volmdlr.TWO_PI, 0))
            if new_pt.vector[0] < 0:
                new_pt.vector[0] = 0
            all_points[pos] = new_pt

    return all_points


def posangle_arc(start, end, radius, frame=None):
    if frame is None:
        p1_new, p2_new = start, end
    else:
        p1_new, p2_new = frame.new_coordinates(start), frame.new_coordinates(end)
    # Angle pour le p1
    u1, u2 = p1_new.x / radius, p1_new.y / radius
    theta1 = sin_cos_angle(u1, u2)
    # Angle pour le p2
    u3, u4 = p2_new.x / radius, p2_new.y / radius
    theta2 = sin_cos_angle(u3, u4)

    if math.isclose(theta1, theta2, abs_tol=1e-6):
        if math.isclose(theta2, 0, abs_tol=1e-6):
            theta2 += volmdlr.TWO_PI
        elif math.isclose(theta1, volmdlr.TWO_PI, abs_tol=1e-6):
            theta1 -= volmdlr.TWO_PI

    return theta1, theta2


def clockwise_interior_from_circle3d(start, end, circle):
    """
    Returns the clockwise interior point between start and end on the circle
    """
    start2d = start.to_2d(plane_origin=circle.frame.origin,
                          x=circle.frame.u, y=circle.frame.v)
    end2d = end.to_2d(plane_origin=circle.frame.origin,
                      x=circle.frame.u, y=circle.frame.v)

    # Angle pour le p1
    u1, u2 = start2d.x / circle.radius, start2d.y / circle.radius
    theta1 = sin_cos_angle(u1, u2)
    # Angle pour le p2
    u3, u4 = end2d.x / circle.radius, end2d.y / circle.radius
    theta2 = sin_cos_angle(u3, u4)

    if theta1 > theta2:
        theta3 = (theta1 + theta2) / 2
    elif theta2 > theta1:
        theta3 = (theta1 + theta2) / 2 + volmdlr.TWO_PI / 2
    else:
        raise NotImplementedError

    if theta3 > volmdlr.TWO_PI:
        theta3 -= volmdlr.TWO_PI

    interior2d = volmdlr.Point2D(circle.radius * math.cos(theta3),
                                 circle.radius * math.sin(theta3))
    interior3d = interior2d.to_3d(plane_origin=circle.frame.origin,
                                  vx=circle.frame.u, vy=circle.frame.v)
    return interior3d


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
    max_angle = min_angle + volmdlr.TWO_PI
    angle = angle % (volmdlr.TWO_PI)

    if math.isclose(angle, min_angle, abs_tol=1e-9):
        return min_angle
    if math.isclose(angle, max_angle, abs_tol=1e-9):
        return max_angle

    return angle


def step_ids_to_str(ids):
    return ','.join(['#{}'.format(i) for i in ids])


class CompositePrimitive(dc.DessiaObject):
    def __init__(self, name=''):
        self.name = name
        self._primitives_to_index = None
        self._utd_primitives_to_index = False
        self.basis_primitives = []

        dc.DessiaObject.__init__(self, name=name)

    def primitive_to_index(self, primitive):
        if not self._utd_primitives_to_index:
            self._primitives_to_index = {p: ip for ip, p in enumerate(self.primitives)}
            self._utd_primitives_to_index = True
        return self._primitives_to_index[primitive]

    def update_basis_primitives(self):
        basis_primitives = []
        for primitive in self.primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.basis_primitives)
            else:
                basis_primitives.append(primitive)

        self.basis_primitives = basis_primitives


class Primitive2D(CompositePrimitive):
    def __init__(self, name=''):
        self.name = name

        dc.DessiaObject.__init__(self, name=name)


class CompositePrimitive2D(Primitive2D):
    """
    A collection of simple primitives
    """
    _non_serializable_attributes = ['name', '_utd_primitives_to_index',
                                    '_primitives_to_index']
    _non_data_hash_attributes = ['name', '_utd_primitives_to_index',
                                 '_primitives_to_index']

    def __init__(self, primitives, name=''):
        Primitive2D.__init__(self, name)
        self.primitives = primitives
        self.update_basis_primitives()

        self._utd_primitives_to_index = False

    # def primitive_to_index(self, primitive):
    #     if not self._utd_primitives_to_index:
    #         self._primitives_to_index = {p: ip for ip, p in enumerate(self.primitives)}
    #         self._utd_primitives_to_index = True
    #     return self._primitives_to_index[primitive]

    # def update_basis_primitives(self):
    #     basis_primitives = []
    #     for primitive in self.primitives:
    #         if hasattr(primitive, 'basis_primitives'):
    #             basis_primitives.extend(primitive.basis_primitives)
    #         else:
    #             basis_primitives.append(primitive)

    #     self.basis_primitives = basis_primitives

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        CompositePrimitive2D rotation
        :param center: rotation center
        :param angle: angle rotation
        :return: a new rotated CompositePrimitive2D
        """
        return self.__class__([point.rotation(center, angle)
                               for point in self.primitives])

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        CompositePrimitive2D rotation. Object is updated inplace
        :param center: rotation center
        :param angle: rotation angle
        """
        for point in self.primitives:
            point.rotation_inplace(center, angle)
        self.update_basis_primitives()

    def translation(self, offset: volmdlr.Vector2D):
        """
        CompositePrimitive2D translation
        :param offset: translation vector
        :return: A new translated CompositePrimitive2D
        """
        return self.__class__([point.translation(offset)
                               for point in self.primitives])

    def translation_inplace(self, offset: volmdlr.Vector2D):
        """
        CompositePrimitive2D translation. Object is updated inplace
        :param offset: translation vector
        """
        for point in self.primitives:
            point.translation_inplace(offset)
        self.update_basis_primitives()

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and return a new CompositePrimitive2D
        side = 'old' or 'new'
        """
        return self.__class__([primitive.frame_mapping(frame, side)
                               for primitive in self.primitives])

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and the object is updated inplace
        side = 'old' or 'new'
        """
        for primitive in self.primitives:
            primitive.frame_mapping_inplace(frame, side)
        self.update_basis_primitives()

    def plot(self, ax=None, color='k', alpha=1,
             plot_points=False, equal_aspect=True):

        if ax is None:
            fig, ax = plt.subplots()

        if equal_aspect:
            ax.set_aspect('equal')

        for element in self.primitives:
            element.plot(ax=ax, color=color, alpha=alpha)  # , plot_points=plot_points)

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


class Primitive3D(dc.PhysicalObject, CompositePrimitive):
    def __init__(self, color=None, alpha=1, name=''):
        self.color = color
        self.alpha = alpha

        dc.PhysicalObject.__init__(self, name=name)

    def volmdlr_primitives(self):
        return [self]

    def babylon_param(self):
        babylon_param = {'alpha': self.alpha,
                         'name': self.name,
                         }
        if self.color is None:
            babylon_param['color'] = [0.8, 0.8, 0.8]
        else:
            babylon_param['color'] = list(self.color)

        return babylon_param

    def triangulation(self):
        raise NotImplementedError(
            'triangulation method should be implemented on class {}'.format(
                self.__class__.__name__))

    def babylon_meshes(self):
        mesh = self.triangulation()
        if mesh is None:
            return []
        positions, indices = mesh.to_babylon()

        babylon_mesh = {'positions': positions,
                        'indices': indices
                        }
        babylon_mesh.update(self.babylon_param())
        return [babylon_mesh]

    def babylon_points(self):
        points = []
        if hasattr(self, 'primitives'):
            points = [[self.primitives[0].start.x,
                       self.primitives[0].start.y,
                       self.primitives[0].start.z]]
            for primitive in self.primitives:
                if hasattr(primitive, 'curve'):
                    points.extend(primitive.curve.evalpts)
                else:
                    points.append([primitive.end.x, primitive.end.y,
                                   primitive.end.z])
        elif hasattr(self, 'curve'):
            points = self.curve.evalpts
        elif hasattr(self, 'polygon_points'):
            points = self.polygon_points(50)
        return points

    def babylon_lines(self, points=None):
        if points is None:
            points = self.babylon_points()
        babylon_lines = {'points': points}
        babylon_lines.update(self.babylon_param())
        return [babylon_lines]

    def babylon_curves(self):
        points = self.babylon_points()
        babylon_curves = self.babylon_lines(points)[0]
        return babylon_curves


class CompositePrimitive3D(Primitive3D):
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = ['basis_primitives']
    _non_data_eq_attributes = ['name', 'basis_primitives']
    _non_data_hash_attributes = []
    """
    A collection of simple primitives3D
    """

    def __init__(self, primitives: List[Primitive3D], name: str = ''):
        self.primitives = primitives

        Primitive3D.__init__(self, name=name)
        self._utd_primitives_to_index = False

    # def primitive_to_index(self, primitive):
    #     if not self._utd_primitives_to_index:
    #         self._primitives_to_index = {p: ip for ip, p in enumerate(self.primitives)}
    #         self._utd_primitives_to_index = True
    #     return self._primitives_to_index[primitive]

    # def update_basis_primitives(self):
    #     # TODO: This is a copy/paste from CompositePrimitive2D, in the future make a Common abstract class
    #     basis_primitives = []
    #     for primitive in self.primitives:
    #         if hasattr(primitive, 'basis_primitives'):
    #             basis_primitives.extend(primitive.primitives)
    #         else:
    #             basis_primitives.append(primitive)

    #     self.basis_primitives = basis_primitives

    # # def to_2d(self, plane_origin, x, y):
    # #     if name is None:
    # #         name = '2D of {}'.format(self.name)
    # #     primitives2d = [p.to_2d(plane_origin, x, y) for p in self.primitives]
    # #     return CompositePrimitive2D(primitives2d, name)

    def plot(self, ax=None, color='k', alpha=1, edge_details=False):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        for primitive in self.primitives:
            primitive.plot(ax=ax, color=color, alpha=alpha)
        return ax


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
        self.name = name

        self.center = volmdlr.Point3D(0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax))

    def __hash__(self):
        return sum(hash(point) for point in self.points)

    def __add__(self, other_bbox):
        return BoundingBox(min(self.xmin, other_bbox.xmin),
                           max(self.xmax, other_bbox.xmax),
                           min(self.ymin, other_bbox.ymin),
                           max(self.ymax, other_bbox.ymax),
                           min(self.zmin, other_bbox.zmin),
                           max(self.zmax, other_bbox.zmax))

    def to_dict(self, use_pointers: bool = True, memo=None, path: str = '#'):
        return {'object_class': 'volmdlr.edges.BoundingBox',
                'name': self.name,
                'xmin': self.xmin,
                'xmax': self.xmax,
                'ymin': self.ymin,
                'ymax': self.ymax,
                'zmin': self.zmin,
                'zmax': self.zmax,
                }

    @property
    def points(self):
        return [volmdlr.Point3D(self.xmin, self.ymin, self.zmin),
                volmdlr.Point3D(self.xmax, self.ymin, self.zmin),
                volmdlr.Point3D(self.xmax, self.ymax, self.zmin),
                volmdlr.Point3D(self.xmin, self.ymax, self.zmin),
                volmdlr.Point3D(self.xmin, self.ymin, self.zmax),
                volmdlr.Point3D(self.xmax, self.ymin, self.zmax),
                volmdlr.Point3D(self.xmax, self.ymax, self.zmax),
                volmdlr.Point3D(self.xmin, self.ymax, self.zmax)]

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
        # plt.show()
        return ax

    @classmethod
    def from_points(cls, points):
        # if len(points) == 0:
        #     return (0, 0, 0, 0, 0, 0)
        xmin = min(pt.x for pt in points)
        xmax = max(pt.x for pt in points)
        ymin = min(pt.y for pt in points)
        ymax = max(pt.y for pt in points)
        zmin = min(pt.z for pt in points)
        zmax = max(pt.z for pt in points)
        return cls(xmin, xmax, ymin, ymax, zmin, zmax)

    def to_frame(self):
        x = volmdlr.Vector3D((self.xmax - self.xmin), 0, 0)
        y = volmdlr.Vector3D(0, (self.ymax - self.ymin), 0)
        z = volmdlr.Vector3D(0, 0, (self.zmax - self.zmin))
        return volmdlr.Frame3D(self.center, x, y, z)

    def volume(self):
        return (self.xmax - self.xmin) * (self.ymax - self.ymin) * (
                    self.zmax - self.zmin)

    def bbox_intersection(self, bbox2):
        return self.xmin < bbox2.xmax and self.xmax > bbox2.xmin \
                and self.ymin < bbox2.ymax and self.ymax > bbox2.ymin \
                and self.zmin < bbox2.zmax and self.zmax > bbox2.zmin

    def is_inside_bbox(self, bbox2):
        return (self.xmin >= bbox2.xmin - 1e-6) and (self.xmax <= bbox2.xmax + 1e-6)\
                and (self.ymin >= bbox2.ymin - 1e-6) and (self.ymax <= bbox2.ymax + 1e-6) \
                and (self.zmin >= bbox2.zmin - 1e-6) and (self.zmax <= bbox2.zmax + 1e-6)

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
        dx = max(permute_bbox2.xmin - permute_bbox1.xmax, 0)

        if permute_bbox2.ymin < permute_bbox1.ymin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dy = max(permute_bbox2.ymin - permute_bbox1.ymax, 0)

        if permute_bbox2.zmin < permute_bbox1.zmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dz = max(permute_bbox2.zmin - permute_bbox1.zmax, 0)

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


class VolumeModel(dc.PhysicalObject):
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = ['shells', 'bounding_box']
    _non_data_eq_attributes = ['name', 'shells', 'bounding_box', 'contours',
                               'faces']
    _non_data_hash_attributes = ['name', 'shells', 'bounding_box', 'contours',
                                 'faces']
    _dessia_methods = ['to_stl_model']
    """
    :param groups: A list of two element tuple. The first element is a string naming the group and the second element is a list of primitives of the group
    """

    def __init__(self, primitives: List[Primitive3D], name: str = ''):
        self.primitives = primitives
        self.name = name
        self.shells = []
        # if self.primitives:
        #     self.shells = self._extract_shells()
        # if self.shells:

        self.bounding_box = self._bounding_box()
        # else:
        #     self.bounding_box = BoundingBox(-1, 1, -1, 1, -1, 1)
        dc.DessiaObject.__init__(self, name=name)

    def __hash__(self):
        return sum(hash(point) for point in self.primitives)

    # def _extract_shells(self):
    #     shells = []
    #     for primitive in self.primitives:
    #         if isinstance(primitive, Shell3D):
    #             shells.append(primitive)
    #     return shells

    def __eq__(self, other):
        if self.__class__.__name__ != other.__class__.__name__:
            return False
        equ = True
        if len(self.primitives) != len(other.primitives):
            return False
        for p1, p2 in zip(self.primitives, other.primitives):
            equ = equ and p1 == p2
        return equ

    def volume(self) -> float:
        """
        Return the sum of volumes of the primitives
        """
        volume = 0
        for primitive in self.primitives:
            volume += primitive.volume()
        return volume

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        VolumeModel rotation
        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated VolumeModel
        """
        new_primitives = [
            primitive.rotation(center, axis, angle) for
            primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        VolumeModel rotation. Object is updated inplace
        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        for primitive in self.primitives:
            primitive.rotation_inplace(center, axis, angle)
        self.bounding_box = self._bounding_box()

    def translation(self, offset: volmdlr.Vector3D):
        """
        VolumeModel translation
        :param offset: translation vector
        :return: A new translated VolumeModel
        """
        new_primitives = [primitive.translation(offset) for
                          primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        VolumeModel translation. Object is updated inplace
        :param offset: translation vector
        """
        for primitives in self.primitives:
            primitives.translation_inplace(offset)
        self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new VolumeModel
        side = 'old' or 'new'
        """
        new_primitives = [primitive.frame_mapping(frame, side)
                          for primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace
        side = 'old' or 'new'
        """
        for primitives in self.primitives:
            primitives.frame_mapping_inplace(frame, side)
        self.bounding_box = self._bounding_box()

    def copy(self, deep=True, memo=None):
        """
        Specific copy
        """
        new_primitives = [primitive.copy(deep=deep, memo=memo) for primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def _bounding_box(self) -> BoundingBox:
        """
        Computes the bounding box of the model
        """
        bboxes = []
        points = []
        for primitive in self.primitives:
            if hasattr(primitive, 'bounding_box'):
                bboxes.append(primitive.bounding_box)
            else:
                if primitive.__class__.__name__ == 'volmdlr.Point3D':
                    points.append(primitive)
        if bboxes:
            xmin = min(box.xmin for box in bboxes)
            xmax = max(box.xmax for box in bboxes)
            ymin = min(box.ymin for box in bboxes)
            ymax = max(box.ymax for box in bboxes)
            zmin = min(box.zmin for box in bboxes)
            zmax = max(box.zmax for box in bboxes)
        elif points:
            xmin = min(p[0] for p in points)
            xmax = max(p[0] for p in points)
            ymin = min(p[1] for p in points)
            ymax = max(p[1] for p in points)
            zmin = min(p[2] for p in points)
            zmax = max(p[2] for p in points)
        else:
            # raise ValueError('Bounding box cant be determined')
            return BoundingBox(-1, 1, -1, 1, 1 - 1, 1)
        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def plot2d(self, ax=None, color=None):
        fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')

        for i, shell in enumerate(self.shells):
            bbox = shell.bbox()
            bbox.plot(ax, color[i])

        return ax

    def plot(self, equal_aspect=True):
        """
        Matplotlib plot of model.
        To use for debug.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', adjustable='box')
        for primitive in self.primitives:
            primitive.plot(ax)
        if equal_aspect:
            ax.set_aspect('equal')

        ax.margins(0.1)
        return ax

    def freecad_script(self, fcstd_filepath,
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
                    # if isinstance(primitive, BSplineCurve3D) \
                    #         or isinstance(primitive, BSplineSurface3D) \
                    #         or isinstance(primitive, Circle3D) \
                    #         or isinstance(primitive, LineSegment3D) \
                    #         or isinstance(primitive, Ellipse3D):
                    #     #                            print(primitive)
                    #     #                            s += 'S = Part.Shape([primitive{}])\n'.format(ip)
                    #     #                            s += 'shapeobj.Shape = S\n'
                    #     s += 'shapeobj.Shape = primitive{}.toShape()\n'.format(
                    #         ip)
                    # else:
                    s += "shapeobj.Shape = primitive{}\n".format(ip)
                    s += 'part.addObject(shapeobj)\n\n'
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

    def freecad_export(self, fcstd_filepath,
                       python_path='python3',
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
        s = self.freecad_script(fcstd_filepath,
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

    # def babylon_script(self, use_cdn=True, debug=False):
    #     # env = Environment(loader=PackageLoader('volmdlr', 'templates'),
    #     #                   autoescape=select_autoescape(['html', 'xml']))
    #     #
    #     # template = env.get_template('babylon.html')
    #
    #     bbox = self._bounding_box()
    #     center = bbox.center
    #     max_length = max([bbox.xmax - bbox.xmin,
    #                       bbox.ymax - bbox.ymin,
    #                       bbox.zmax - bbox.zmin])
    #
    #     primitives_strings = []
    #     for primitive in self.primitives:
    #         if hasattr(primitive, 'babylon_script'):
    #             primitives_strings.append(primitive.babylon_script())
    #
    #     return template.render(name=self.name,
    #                            center=tuple(center),
    #                            length=2 * max_length,
    #                            primitives_strings=primitives_strings,
    #                            use_cdn=use_cdn,
    #                            debug=debug)
    #
    # def babylonjs_from_script(self, page_name=None, use_cdn=True, debug=False):
    #     script = self.babylon_script(use_cdn=use_cdn, debug=debug)
    #
    #     if page_name is None:
    #         with tempfile.NamedTemporaryFile(suffix=".html",
    #                                          delete=False) as file:
    #             file.write(bytes(script, 'utf8'))
    #         page_name = file.name
    #     else:
    #         page_name += '.html'
    #         with open(page_name, 'w')  as file:
    #             file.write(script)
    #
    #     webbrowser.open('file://' + os.path.realpath(page_name))

    def babylon_data(self):
        meshes = []
        lines = []
        for primitive in self.primitives:
            if hasattr(primitive, 'babylon_meshes'):
                meshes.extend(primitive.babylon_meshes())
            if hasattr(primitive, 'babylon_curves'):
                lines.append(primitive.babylon_curves())
        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])

        babylon_data = {'meshes': meshes,
                        'lines': lines,
                        'max_length': max_length,
                        'center': list(center)}

        return babylon_data

    @classmethod
    def babylonjs_script(cls, babylon_data, use_cdn=True,
                         debug=False):
        if use_cdn:
            script = volmdlr.templates.babylon_unpacker_cdn_header  # .substitute(name=page_name)
        else:
            script = volmdlr.templates.babylon_unpacker_embedded_header  # .substitute(name=page_name)

        script += volmdlr.templates.babylon_unpacker_body_template.substitute(
            babylon_data=babylon_data)
        return script

    def babylonjs(self, page_name=None, use_cdn=True, debug=False):
        babylon_data = self.babylon_data()
        script = self.babylonjs_script(babylon_data, use_cdn=use_cdn,
                                       debug=debug)
        if page_name is None:
            with tempfile.NamedTemporaryFile(suffix=".html",
                                             delete=False) as file:
                file.write(bytes(script, 'utf8'))
            page_name = file.name
        else:
            if not page_name.endswith('.html'):
                page_name += '.html'
            with open(page_name, 'w') as file:
                file.write(script)

        webbrowser.open('file://' + os.path.realpath(page_name))

        return page_name

    def save_babylonjs_to_file(self, filename: str = None,
                               use_cdn=True, debug=False):
        babylon_data = self.babylon_data()
        script = self.babylonjs_script(babylon_data, use_cdn=use_cdn,
                                       debug=debug)
        if filename is None:
            with tempfile.NamedTemporaryFile(suffix=".html",
                                             delete=False) as file:
                file.write(bytes(script, 'utf8'))
                return file.name

        if not filename.endswith('.html'):
            filename += '.html'

            with open(filename, 'w') as file:
                file.write(script)
            return filename

    def to_stl_model(self):
        mesh = self.primitives[0].triangulation()
        for primitive in self.primitives[1:]:
            mesh.merge_mesh(primitive.triangulation())
        stl = mesh.to_stl()
        return stl

    def to_stl(self, filepath: str):
        if not filepath.endswith('.stl'):
            filepath += '.stl'
        with open(filepath, 'wb') as file:
            self.to_stl_stream(file)

    def to_stl_stream(self, stream: dcf.BinaryFile):
        stl = self.to_stl_model()
        stl.save_to_stream(stream)
        return stream

    def to_step(self, filepath: str):
        if not (filepath.endswith('.step') or filepath.endswith('.stp')):
            filepath += '.step'
        with open(filepath, 'w') as file:
            self.to_step_stream(file)

    def to_step_stream(self, stream: dcf.StringFile):

        step_content = STEP_HEADER.format(name=self.name,
                                          filename='',
                                          timestamp=datetime.now().isoformat(),
                                          version=volmdlr.__version__)
        current_id = 8

        for primitive in self.primitives:
            primitive_content, primitive_id = primitive.to_step(current_id)

            step_content += primitive_content

            product_definition_context_id = primitive_id + 1
            step_content += "#{} = PRODUCT_DEFINITION_CONTEXT('part definition',#2,'design');\n"\
                .format(product_definition_context_id)

            product_context_id = product_definition_context_id + 1
            step_content += "#{} = PRODUCT_CONTEXT('',#2,'mechanical');\n".format(product_context_id)
            product_id = product_context_id + 1
            step_content += "#{} = PRODUCT('{}','{}','',(#{}));\n".format(product_id,
                                                                          primitive.name,
                                                                          primitive.name,
                                                                          product_context_id)
            product_definition_formation_id = product_id + 1
            step_content += "#{} = PRODUCT_DEFINITION_FORMATION('','',#{});\n".format(
                product_definition_formation_id, product_id)
            product_definition_id = product_definition_formation_id + 1
            step_content += "#{} = PRODUCT_DEFINITION('design','',#{},#{});\n".format(product_definition_id,
                                                                                      product_definition_formation_id,
                                                                                      product_definition_context_id)
            product_definition_shape_id = product_definition_id + 1
            step_content += "#{} = PRODUCT_DEFINITION_SHAPE('','',#{});\n".format(
                product_definition_shape_id, product_definition_id)
            shape_definition_repr_id = product_definition_shape_id + 1
            step_content += "#{} = SHAPE_DEFINITION_REPRESENTATION(#{},#{});\n".format(shape_definition_repr_id,
                                                                                       product_definition_shape_id,
                                                                                       primitive_id
                                                                                       )
            product_related_category = shape_definition_repr_id + 1
            step_content += "#{} = PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#{}));\n".format(
                product_related_category,
                product_id
            )
            draughting_id = product_related_category + 1
            step_content += "#{} = DRAUGHTING_PRE_DEFINED_CURVE_FONT('continuous');\n".format(
                draughting_id)
            color_id = draughting_id + 1
            primitive_color = (1, 1, 1)
            if hasattr(primitive, 'color'):
                primitive_color = primitive.color
            step_content += "#{} = COLOUR_RGB('',{}, {}, {});\n".format(
                color_id,
                round(float(primitive_color[0]), 4),
                round(float(primitive_color[1]), 4),
                round(float(primitive_color[2]), 4)
            )

            curve_style_id = color_id + 1
            step_content += "#{} = CURVE_STYLE('',#{},POSITIVE_LENGTH_MEASURE(0.1),#{});\n".format(
                    curve_style_id, draughting_id, color_id)

            fill_area_color_id = curve_style_id + 1
            step_content += "#{} = FILL_AREA_STYLE_COLOUR('',#{});\n".format(
                    fill_area_color_id, color_id)

            fill_area_id = fill_area_color_id + 1
            step_content += "#{} = FILL_AREA_STYLE('',#{});\n".format(
                    fill_area_id, fill_area_color_id)

            suface_fill_area_id = fill_area_id + 1
            step_content += "#{} = SURFACE_STYLE_FILL_AREA(#{});\n".format(
                    suface_fill_area_id, fill_area_id)

            suface_side_style_id = suface_fill_area_id + 1
            step_content += "#{} = SURFACE_SIDE_STYLE('',(#{}));\n".format(
                    suface_side_style_id, suface_fill_area_id)

            suface_style_usage_id = suface_side_style_id + 1
            step_content += "#{} = SURFACE_STYLE_USAGE(.BOTH.,#{});\n".format(
                    suface_style_usage_id, suface_side_style_id)

            presentation_style_id = suface_style_usage_id + 1

            step_content += "#{} = PRESENTATION_STYLE_ASSIGNMENT((#{},#{}));\n".format(
                    presentation_style_id, suface_style_usage_id, curve_style_id)

            styled_item_id = presentation_style_id + 1
            step_content += "#{} = STYLED_ITEM('color',(#{}),#{});\n".format(
                    styled_item_id, presentation_style_id, primitive_id)

            current_id = styled_item_id + 1

        step_content += STEP_FOOTER

        stream.write(step_content)

    def volmdlr_volume_model(self):
        return [self]

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
                primitive.frame_mapping(frame, side='old'))
        return VolumeModel(primitives)

    # def babylon_script(self, use_cdn=True, debug=False):
    #
    #     env = Environment(loader=PackageLoader('volmdlr', 'templates'),
    #                       autoescape=select_autoescape(['html', 'xml']))
    #
    #     template = env.get_template('babylon.html')
    #
    #     bbox = self._bounding_box()
    #     center = bbox.center
    #     max_length = max([bbox.xmax - bbox.xmin,
    #                       bbox.ymax - bbox.ymin,
    #                       bbox.zmax - bbox.zmin])
    #
    #     primitives_strings = []
    #     for primitive in self.primitives:
    #         if hasattr(primitive, 'babylon_script'):
    #             primitives_strings.append(primitive.babylon_script())
    #
    #     positions = []
    #     orientations = []
    #     for step in self.step_frames:
    #         step_positions = []
    #         step_orientations = []
    #
    #         for frame in step:
    #             step_positions.append(list(frame.origin))
    #             step_orientations.append([list(frame.u),
    #                                       list(frame.v),
    #                                       list(frame.w)])
    #
    #         positions.append(step_positions)
    #         orientations.append(step_orientations)
    #
    #     return template.render(name=self.name,
    #                            center=tuple(center),
    #                            length=2 * max_length,
    #                            primitives_strings=primitives_strings,
    #                            positions=positions,
    #                            orientations=orientations,
    #                            use_cdn=use_cdn,
    #                            debug=debug)

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
