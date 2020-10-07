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


def step_split_arguments(function_arg):
    """
    Split the arguments of a function that doesn't start with '(' but end with ')'
    ex: IN: '#123,#124,#125)'
       OUT: ['#123', '#124', '#125']
    """
    if len(function_arg) > 0 and function_arg[-1] != ')':
        function_arg += ')'
    arguments = []
    argument = ""
    parenthesis = 1
    for char in function_arg:
        if char == "(":
            parenthesis += 1

        if char != "," or parenthesis > 1:
            argument += char
        else:
            arguments.append(argument)
            argument = ""

        if char == ")":
            parenthesis -= 1
            if parenthesis == 0:
                arguments.append(argument[:-1])
                argument = ""
                break
    return arguments


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

class Edge(dc.DessiaObject):
    def __init__(self, start, end, name=''):
        self.start = start
        self.end = end
        dc.DessiaObject.__init__(self, name=name)

    @classmethod
    def from_step(cls, arguments, object_dict):
        if object_dict[arguments[3]].__class__ is Line3D:
            return LineSegment3D(object_dict[arguments[1]],
                                 object_dict[arguments[2]], arguments[0][1:-1])

        elif object_dict[arguments[3]].__class__ is Circle3D:
            # We supposed that STEP file is reading on trigo way
            center = object_dict[arguments[3]].center
            normal = object_dict[arguments[3]].normal
            normal.Normalize()
            radius = object_dict[arguments[3]].radius
            p1 = object_dict[arguments[1]]
            p2 = object_dict[arguments[2]]
            other_vec = object_dict[arguments[3]].other_vec
            if other_vec is None:
                other_vec = p1 - center
            other_vec.Normalize()
            frame = Frame3D(center, other_vec, normal.Cross(other_vec), normal)
            if p1 == p2:
                angle = math.pi
            else:
                # p1_new, p2_new = frame.NewCoordinates(p1), frame.NewCoordinates(p2)
                # #Angle for p1
                # u1, u2 = p1_new.vector[0]/radius, p1_new.vector[1]/radius
                # theta1 = sin_cos_angle(u1, u2)
                # #Angle for p2
                # u3, u4 = p2_new.vector[0]/radius, p2_new.vector[1]/radius
                # theta2 = sin_cos_angle(u3, u4)
                theta1, theta2 = posangle_arc(p1, p2, radius, frame)
                if theta1 > theta2:  # sens trigo
                    angle = math.pi + (theta1 + theta2) / 2
                else:
                    angle = (theta1 + theta2) / 2
            p_3 = Point3D(
                (radius * math.cos(angle), radius * math.sin(angle), 0))
            p3 = frame.OldCoordinates(p_3)
            if p1 == p3 or p2 == p3:
                p_3 = Point3D((radius * math.cos(0), radius * math.sin(0), 0))
                p3 = frame.OldCoordinates(p_3)
            arc = Arc3D(p1, p3, p2, normal, arguments[0][1:-1], other_vec)
            if math.isclose(arc.radius, 0, abs_tol=1e-9):
                if p1 == p2:
                    p_3 = Point3D(
                        (radius * math.cos(0), radius * math.sin(0), 0))
                    p3 = frame.OldCoordinates(p_3)
                    arc = Arc3D(p1, p3, p2, normal, arguments[0][1:-1],
                                other_vec)
            return arc

        elif object_dict[arguments[3]].__class__ is Ellipse3D:
            majorax = object_dict[arguments[3]].major_axis
            minorax = object_dict[arguments[3]].minor_axis
            center = object_dict[arguments[3]].center
            normal = object_dict[arguments[3]].normal
            normal.Normalize()
            majordir = object_dict[arguments[3]].major_dir
            majordir.Normalize()
            minordir = normal.Cross(majordir)
            minordir.Normalize()
            frame = Frame3D(center, majordir, minordir, normal)
            p1 = object_dict[
                arguments[1]]  # on part du principe que p1 suivant majordir
            p2 = object_dict[arguments[2]]
            if p1 == p2:
                angle = 5 * math.pi / 4
                xtra = Point3D((majorax * math.cos(math.pi / 2),
                                minorax * math.sin(math.pi / 2), 0))
                extra = frame.OldCoordinates(xtra)
            else:
                extra = None
                ## Positionnement des points dans leur frame
                p1_new, p2_new = frame.NewCoordinates(
                    p1), frame.NewCoordinates(p2)
                # Angle pour le p1
                u1, u2 = p1_new.vector[0] / majorax, p1_new.vector[1] / minorax
                theta1 = sin_cos_angle(u1, u2)
                # Angle pour le p2
                u3, u4 = p2_new.vector[0] / majorax, p2_new.vector[1] / minorax
                theta2 = sin_cos_angle(u3, u4)

                if theta1 > theta2:  # sens trigo
                    angle = math.pi + (theta1 + theta2) / 2
                else:
                    angle = (theta1 + theta2) / 2

            p_3 = Point3D(
                (majorax * math.cos(angle), minorax * math.sin(angle), 0))
            p3 = frame.OldCoordinates(p_3)

            arcellipse = ArcEllipse3D(p1, p3, p2, center, majordir, normal,
                                      arguments[0][1:-1], extra)

            return arcellipse

        elif object_dict[arguments[3]].__class__ is BSplineCurve3D:
            # print(object_dict[arguments[1]], object_dict[arguments[2]])
            # BSplineCurve3D à couper à gauche et à droite avec les points ci dessus ?
            return object_dict[arguments[3]]

        else:
            print(object_dict[arguments[3]])
            raise NotImplementedError

class Line(dc.DessiaObject):
    """
    Abstract class
    """
    def __init__(self, point1, point2, name=''):
        self.point1 = point1
        self.point2 = point2

    def unit_direction_vector(self):
        u = self.direction_vector()
        u.Normalize()
        return u

    def direction_vector(self):
        return self.point2 - self.point1

    def normal_vector(self):
        return self.unit_direction_vector().normal_vector()

    def point_projection(self, point):

        u = self.point2 - self.point1
        norm_u = u.Norm()
        t = (point - self.point1).Dot(u) / norm_u ** 2
        projection = self.point1 + t * u

        return projection, t * norm_u

    def split(self, split_point):
        return [self.__class__(self.point1, split_point),
                self.__class__(split_point, self.point2)]

class LineSegment(Edge):
    """
    Abstract class
    """

    def unit_direction_vector(self):
        u = self.direction_vector()
        u.Normalize()
        return u

    def direction_vector(self):
        return self.end - self.start

    def normal_vector(self):
        return self.unit_direction_vector().normal_vector()

    def point_projection(self, point):
        p1, p2 = self.points
        u = p2 - p1
        norm_u = u.Norm()
        t = (point - p1).Dot(u) / norm_u ** 2
        projection = p1 + t * u

        return projection, t * norm_u

    def split(self, split_point):
        return [self.__class__(self.start, split_point),
                self.__class__(split_point, self.end)]

class Line2D(Line):
    """
    Define an infinite line given by two points.
    """

    def __init__(self, point1, point2, *, name=''):
        Line.__init__(self, point1, point2, name=name)


    def To3D(self, plane_origin, x1, x2):
        p3D = [p.To3D(plane_origin, x1, x2) for p in self.points]
        return Line2D(*p3D, self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Line2D(
                *[p.Rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Line2D(
                *[p.Translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset, copy=False)

    def MPLPlot(self, ax=None, color='k', dashed=True):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')

        p1, p2 = self.points

        if version.parse(_mpl_version) >= version.parse('3.3.2'):
            if dashed:
                ax.axline(p1.vector, p2.vector, dashes=[30, 5, 10, 5])
            else:
                ax.axline(p1.vector, p2.vector)
        else:
            u = p2 - p1
            p3 = p1 - 3 * u
            p4 = p2 + 4 * u
            if dashed:
                ax.plot([p3[0], p4[0]], [p3[1], p4[1]], color=color,
                        dashes=[30, 5, 10, 5])
            else:
                ax.plot([p3[0], p4[0]], [p3[1], p4[1]], color=color)

        return ax

    def plot_data(self, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        p1, p2 = self.points
        u = p2 - p1
        p3 = p1 - 3 * u
        p4 = p2 + 4 * u
        return {'type': 'line',
                'data': [p3[0], p3[1],
                         p4[0], p4[1]],
                'color': color,
                'marker': marker,
                'size': stroke_width,
                'dash': dash,
                'opacity': opacity,
                'arrow': arrow
                }

    def CreateTangentCircle(self, point, other_line):
        """
        Computes the two circles that are tangent to 2 lines and intersect
        a point located on one of the two lines.
        """

        # point will be called I(x_I, y_I)
        # self will be (AB)
        # line will be (CD)

        if math.isclose(self.point_distance(point), 0, abs_tol=1e-10):
            I = Vector2D((point[0], point[1]))
            A = Vector2D((self.points[0][0], self.points[0][1]))
            B = Vector2D((self.points[1][0], self.points[1][1]))
            C = Vector2D((other_line.points[0][0], other_line.points[0][1]))
            D = Vector2D((other_line.points[1][0], other_line.points[1][1]))

        elif math.isclose(other_line.point_distance(point), 0, abs_tol=1e-10):
            I = Vector2D((point[0], point[1]))
            C = Vector2D((self.points[0][0], self.points[0][1]))
            D = Vector2D((self.points[1][0], self.points[1][1]))
            A = Vector2D((other_line.points[0][0], other_line.points[0][1]))
            B = Vector2D((other_line.points[1][0], other_line.points[1][1]))
        else:
            raise AttributeError("The point isn't on any of the two lines")

        # CHANGEMENT DE REPAIRE
        new_u = Vector2D((B - A))
        new_u.Normalize()
        new_v = new_u.NormalVector(unit=True)
        new_basis = Frame2D(I, new_u, new_v)

        new_A = new_basis.NewCoordinates(A)
        new_B = new_basis.NewCoordinates(B)
        new_C = new_basis.NewCoordinates(C)
        new_D = new_basis.NewCoordinates(D)

        if new_C[1] == 0 and new_D[1] == 0:
            # Segments are on the same line: no solution
            return None, None

        elif math.isclose(self.DirectionVector(unit=True).Dot(
                other_line.NormalVector(unit=True)), 0, abs_tol=1e-06):
            # Parallel segments: one solution

            segments_distance = abs(new_C[1] - new_A[1])
            r = segments_distance / 2
            new_circle_center = Point2D((0, npy.sign(new_C[1] - new_A[1]) * r))
            circle_center = new_basis.OldCoordinates(new_circle_center)
            circle = Circle2D(circle_center, r)

            return circle, None

        elif math.isclose(self.DirectionVector(unit=True).Dot(
                other_line.DirectionVector(unit=True)), 0, abs_tol=1e-06):
            # Perpendicular segments: 2 solution
            line_AB = Line2D(Point2D(new_A), Point2D(new_B))
            line_CD = Line2D(Point2D(new_C), Point2D(new_D))
            new_pt_K = Point2D.LinesIntersection(line_AB, line_CD)

            r = abs(new_pt_K[0])
            new_circle_center1 = Point2D((0, r))
            new_circle_center2 = Point2D((0, -r))
            circle_center1 = new_basis.OldCoordinates(new_circle_center1)
            circle_center2 = new_basis.OldCoordinates(new_circle_center2)
            circle1 = Circle2D(circle_center1, r)
            circle2 = Circle2D(circle_center2, r)

            return circle1, circle2

        # =============================================================================
        # LES SEGMENTS SONT QUELCONQUES
        #   => 2 SOLUTIONS
        # =============================================================================
        else:

            line_AB = Line2D(Point2D(new_A), Point2D(new_B))
            line_CD = Line2D(Point2D(new_C), Point2D(new_D))
            new_pt_K = Point2D.LinesIntersection(line_AB, line_CD)
            pt_K = Point2D(new_basis.OldCoordinates(new_pt_K))

            if pt_K == I:
                return None, None

            # CHANGEMENT DE REPERE:
            new_u2 = Vector2D(pt_K - I)
            new_u2.Normalize()
            new_v2 = new_u2.NormalVector(unit=True)
            new_basis2 = Frame2D(I, new_u2, new_v2)

            new_A = new_basis2.NewCoordinates(A)
            new_B = new_basis2.NewCoordinates(B)
            new_C = new_basis2.NewCoordinates(C)
            new_D = new_basis2.NewCoordinates(D)
            new_pt_K = new_basis2.NewCoordinates(pt_K)

            teta1 = math.atan2(new_C[1], new_C[0] - new_pt_K[0])
            teta2 = math.atan2(new_D[1], new_D[0] - new_pt_K[0])

            if teta1 < 0:
                teta1 += math.pi
            if teta2 < 0:
                teta2 += math.pi

            if not math.isclose(teta1, teta2, abs_tol=1e-08):
                if math.isclose(teta1, math.pi, abs_tol=1e-08) or math.isclose(
                        teta1, 0., abs_tol=1e-08):
                    teta = teta2
                elif math.isclose(teta2, math.pi,
                                  abs_tol=1e-08) or math.isclose(teta2, 0.,
                                                                 abs_tol=1e-08):
                    teta = teta1
            else:
                teta = teta1

            r1 = new_pt_K[0] * math.sin(teta) / (1 + math.cos(teta))
            r2 = new_pt_K[0] * math.sin(teta) / (1 - math.cos(teta))

            new_circle_center1 = Point2D((0, -r1))
            new_circle_center2 = Point2D((0, r2))

            circle_center1 = new_basis2.OldCoordinates(new_circle_center1)
            circle_center2 = new_basis2.OldCoordinates(new_circle_center2)

            if new_basis.NewCoordinates(circle_center1)[1] > 0:
                circle1 = Circle2D(circle_center1, r1)
                circle2 = Circle2D(circle_center2, r2)
            else:
                circle1 = Circle2D(circle_center2, r2)
                circle2 = Circle2D(circle_center1, r1)

            return circle1, circle2


class BSplineCurve2D(Primitive2D):
    def __init__(self, degree, control_points, knot_multiplicities, knots,
                 weights=None, periodic=False, name=''):
        Primitive2D.__init__(self, name=name)
        self.control_points = control_points
        self.degree = degree
        knots = standardize_knot_vector(knots)
        self.knots = knots
        self.knot_multiplicities = knot_multiplicities
        self.weights = weights
        self.periodic = periodic
        self.name = name

        curve = BSpline.Curve()
        curve.degree = degree
        if weights is None:
            P = [(control_points[i][0], control_points[i][1]) for i in
                 range(len(control_points))]
            curve.ctrlpts = P
        else:
            Pw = [(control_points[i][0] * weights[i],
                   control_points[i][1] * weights[i], weights[i]) for i in
                  range(len(control_points))]
            curve.ctrlptsw = Pw
        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot] * knot_multiplicities[i])
        curve.knotvector = knot_vector
        curve.delta = 0.1
        curve_points = curve.evalpts

        self.curve = curve
        self.points = [Point2D((p[0], p[1])) for p in curve_points]

    def Length(self):
        # Approximately
        length = 0
        for k in range(0, len(self.points) - 1):
            length += (self.points[k] - self.points[k + 1]).Norm()
        return length

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        # copy paste from wire3D
        length = 0.
        primitives = []
        for k in range(0, len(self.points) - 1):
            primitives.append(
                LineSegment2D(self.points[k], self.points[k + 1]))
        for primitive in primitives:
            primitive_length = primitive.Length()
            if length + primitive_length >= curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(
                    curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError

    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        x = [p.vector[0] for p in self.points]
        y = [p.vector[1] for p in self.points]
        ax.plot(x, y, 'o-k')
        return fig, ax

    def To3D(self, plane_origin, x1, x2):
        control_points3D = [p.To3D(plane_origin, x1, x2) for p in
                            self.control_points]
        return BSplineCurve3D(self.degree, control_points3D,
                              self.knot_multiplicities, self.knots,
                              self.weights, self.periodic, self.name)

    def tessellation_points(self):
        return self.points


class LineSegment2D(LineSegment):
    """
    Define a line segment limited by two points
    """

    def __init__(self, start, end, *, name=''):
        Edge.__init__(self, start, end, name=name)

    def Length(self):
        return self.end.point_distance(self.start)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.start + self.unit_direction_vector() * curvilinear_abscissa

    def point_distance(self, point, return_other_point=False):
        """
        Computes the distance of a point to segment of line
        """
        if self.point1 == self.point2:
            if return_other_point:
                return 0, Point2D(point)
            return 0
        distance, point = LineSegment2DPointDistance(
            [p.vector for p in self.points], point.vector)
        if return_other_point:
            return distance, Point2D(point)
        return distance

    def point_projection(self, point):
        """
        If the projection falls outside the LineSegment2D, returns None.
        """
        point, curv_abs = Line2D.point_projection(self, point)
        if curv_abs < 0 or curv_abs > self.Length():
            return None, curv_abs
        return point, curv_abs

    def line_intersections(self, line):
        point = Point2D.LinesIntersection(self, line)
        if point is not None:
            point_projection1, _ = self.point_projection(point)
            if point_projection1 is None:
                return []

            if line.__class__.__name__ == 'LineSegment2D':
                point_projection2, _ = line.point_projection(point)
                if point_projection2 is None:
                    return []

            return [point_projection1]
        else:
            return []

    def MPLPlot(self, ax=None, color='k', arrow=False, width=None,
                plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            # ax.set_aspect('equal')
        # else:
        #     fig = ax.figure

        p1, p2 = self.start, self.end
        if arrow:
            if plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color,
                        style='o-')
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color)

            length = ((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2) ** 0.5
            if width is None:
                width = length / 1000.
                head_length = length / 20.
                head_width = head_length / 2.
            else:
                head_width = 2 * width
                head_length = head_width
            ax.arrow(p1[0], p1[1],
                     (p2[0] - p1[0]) / length * (length - head_length),
                     (p2[1] - p1[1]) / length * (length - head_length),
                     head_width=head_width, fc='b', linewidth=0,
                     head_length=head_length, width=width, alpha=0.3)
        else:
            if width is None:
                width = 1
            if plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color,
                        marker='o', linewidth=width)
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color,
                        linewidth=width)
        return ax

    def To3D(self, plane_origin, x1, x2):
        start = self.start.To3D(plane_origin, x1, x2)
        end = self.end.To3D(plane_origin, x1, x2)
        return LineSegment3D(start, end, self.name)

    def reverse(self):
        return LineSegment2D(self.end.copy(), self.points[0].copy())

    def to_line(self):
        return Line2D(*self.points)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return LineSegment2D(
                *[p.Rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return LineSegment2D(
                *[p.Translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset, copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return LineSegment2D(
                    *[frame.OldCoordinates(p) for p in self.points])
            else:
                self.points = [frame.OldCoordinates(p) for p in self.points]
        if side == 'new':
            if copy:
                return LineSegment2D(
                    *[frame.NewCoordinates(p) for p in self.points])
            else:
                self.points = [frame.NewCoordinates(p) for p in self.points]

    def plot_data(self, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        return {'type': 'line',
                'data': [self.points[0].vector[0], self.points[0].vector[1],
                         self.end.vector[0], self.end.vector[1]],
                'color': color,
                'marker': marker,
                'size': stroke_width,
                'dash': dash,
                'opacity': opacity,
                'arrow': arrow
                }

    def CreateTangentCircle(self, point, other_line):
        circle1, circle2 = Line2D.CreateTangentCircle(other_line, point, self)
        if circle1 is not None:
            point_J1, curv_abs1 = Line2D.point_projection(self, circle1.center)
            if curv_abs1 < 0. or curv_abs1 > self.Length():
                circle1 = None
        if circle2 is not None:
            point_J2, curv_abs2 = Line2D.point_projection(self, circle2.center)
            if curv_abs2 < 0. or curv_abs2 > self.Length():
                circle2 = None
        return circle1, circle2

    def tessellation_points(self):
        return [self.start, self.end]

    def polygon_points(self, min_x_density=None, min_y_density=None):
        n = 0# Number of points to insert between start and end
        if min_x_density:
            dx = abs(self.start[0]-self.end[0])
            n = max(n, math.floor(dx*min_x_density))
        if min_y_density:
            dy = abs(self.start[1]-self.end[1])
            n = max(n, math.floor(dy*min_y_density))

        if n:
            l = self.Length()
            return [self.PointAtCurvilinearAbscissa(i*l/(n+1)) for i in range(n+2)]
        else:
            return [self.start, self.end]

class Arc2D(Primitive2D):
    """
    angle: the angle measure always >= 0
    """

    def __init__(self, start, interior, end, name=''):
        Primitive2D.__init__(self, name)
        self.interior = interior
        self.start = start
        self.end = end
        xi, yi = interior.vector
        xe, ye = end.vector
        xs, ys = start.vector
        try:
            A = Matrix22(2 * (xs - xi), 2 * (ys - yi),
                         2 * (xs - xe), 2 * (ys - ye))
            b = - Vector2D((xi ** 2 + yi ** 2 - xs ** 2 - ys ** 2,
                            xe ** 2 + ye ** 2 - xs ** 2 - ys ** 2))
            inv_A = A.inverse()
            x = inv_A.vector_multiplication(b)
            self.center = Point2D(x.vector)
        except ValueError:
            A = npy.array([[2 * (xs - xi), 2 * (ys - yi)],
                           [2 * (xs - xe), 2 * (ys - ye)]])
            b = - npy.array([xi ** 2 + yi ** 2 - xs ** 2 - ys ** 2,
                             xe ** 2 + ye ** 2 - xs ** 2 - ys ** 2])
            self.center = Point2D(solve(A, b))

        r1 = self.start - self.center
        r2 = self.end - self.center
        ri = self.interior - self.center

        self.radius = r1.Norm()
        angle1 = math.atan2(r1.vector[1], r1.vector[0])
        anglei = math.atan2(ri.vector[1], ri.vector[0])
        angle2 = math.atan2(r2.vector[1], r2.vector[0])

        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei + two_pi) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + two_pi

        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2 + two_pi) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + two_pi

        if clockwise_path > trigowise_path:
            self.is_trigo = True
            self.angle1 = angle1
            self.angle2 = angle2
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle1 = angle2
            self.angle2 = angle1
            self.angle = clockwise_path

    def _get_points(self):
        return [self.start, self.interior, self.end]

    points = property(_get_points)


    # def tessellation_points(self, resolution_for_circle=40):
    #     number_points_tesselation = math.ceil(resolution_for_circle*abs(self.angle)/2/math.pi)
    #     if number_points_tesselation == 1:
    #         number_points_tesselation += 1

    #     points = []
    #     if not self.is_trigo:
    #         delta_angle = -abs(self.angle1-self.angle2)/(number_points_tesselation-1)
    #         delta_angle = -self.angle/(number_points_tesselation-1)
    #     else:
    #         delta_angle =  abs(self.angle1-self.angle2)/(number_points_tesselation-1)
    #         delta_angle = self.angle/(number_points_tesselation-1)
    #     points.append(self.start)
    #     for i in range(number_points_tesselation-2):
    #         point_to_add = points[-1].Rotation(self.center, delta_angle)
    #         points.append(point_to_add)
    #     points.append(self.end)
    #     return points

    def tessellation_points(self, resolution_for_circle=40):
        number_points_tesselation = math.ceil(
            resolution_for_circle * abs(self.angle) / 2 / math.pi)
        number_points_tesselation = max(number_points_tesselation, 5)
        l = self.Length()
        return [self.PointAtCurvilinearAbscissa(
            i / (number_points_tesselation - 1) * l) for i in
                range(number_points_tesselation)]

    def point_belongs(self, point):
        """
        Computes if the point belongs to the pizza slice drawn by the arc and its center
        """
        circle = Circle2D(self.center, self.radius)
        if not circle.point_belongs(point):
            return False
        vector_start = self.start - self.center
        vector_point = point - self.center
        vector_end = self.end - self.center
        if self.is_trigo:
            vector_start, vector_end = vector_end, vector_start
        arc_angle = clockwise_angle(vector_start, vector_end)
        point_angle = clockwise_angle(vector_start, vector_point)
        if point_angle <= arc_angle:
            return True

    def point_distance(self, point):
        vector_start = self.start - self.center
        vector_point = point - self.center
        vector_end = self.end - self.center
        if self.is_trigo:
            vector_start, vector_end = vector_end, vector_start
        arc_angle = clockwise_angle(vector_start, vector_end)
        point_angle = clockwise_angle(vector_start, vector_point)
        if point_angle <= arc_angle:
            return abs(
                LineSegment2D(point, self.center).Length() - self.radius)
        else:
            return min(LineSegment2D(point, self.start).Length(),
                       LineSegment2D(point, self.end).Length())

    def line_intersections(self, line):
        circle = Circle2D(self.center, self.radius)
        circle_intersection_points = circle.line_intersections(line)

        if circle_intersection_points is None:
            return None

        intersection_points = []
        for pt in circle_intersection_points:
            if self.point_belongs(pt):
                intersection_points.append(pt)
        return intersection_points

    def Length(self):
        return self.radius * abs(self.angle)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        if self.is_trigo:
            return self.start.Rotation(self.center,
                                       curvilinear_abscissa / self.radius)
            # return self.start.Rotation(self.center, curvilinear_abscissa*self.angle)
        else:
            return self.start.Rotation(self.center,
                                       -curvilinear_abscissa / self.radius)
            # return self.start.Rotation(self.center, -curvilinear_abscissa*self.angle)

    def MiddlePoint(self):
        l = self.Length()
        return self.PointAtCurvilinearAbscissa(0.5 * l)

    def Area(self):
        if self.angle2 < self.angle1:
            angle = self.angle2 + two_pi - self.angle1
        else:
            angle = self.angle2 - self.angle1
        return self.radius ** 2 * angle / 2

    def CenterOfMass(self):
        #        u=self.middle.vector-self.center.vector
        u = self.MiddlePoint() - self.center
        u.Normalize()
        alpha = abs(self.angle)
        return self.center + 4 / (3 * alpha) * self.radius * math.sin(
            alpha * 0.5) * u

    def MPLPlot(self, ax=None, color='k', plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')

        if plot_points:
            for p in [self.center, self.start, self.interior, self.end]:
                p.MPLPlot(ax=ax)

        pc = self.center.vector
        ax.add_patch(Arc(pc, 2 * self.radius, 2 * self.radius, angle=0,
                         theta1=self.angle1 * 0.5 / math.pi * 360,
                         theta2=self.angle2 * 0.5 / math.pi * 360,
                         color=color))

        return ax

    def To3D(self, plane_origin, x, y):
        ps = self.start.To3D(plane_origin, x, y)
        pi = self.interior.To3D(plane_origin, x, y)
        pe = self.end.To3D(plane_origin, x, y)

        return Arc3D(ps, pi, pe, name=self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Arc2D(*[p.Rotation(center, angle, copy=True) for p in
                           [self.start, self.interior, self.end]])
        else:
            self.__init__(*[p.Rotation(center, angle, copy=True) for p in
                            [self.start, self.interior, self.end]])


    def Translation(self, offset, copy=True):
        if copy:
            return Arc2D(*[p.Translation(offset, copy=True) for p in
                           [self.start, self.interior, self.end]])
        else:
            self.__init__(*[p.Translation(offset, copy=True) for p in
                            [self.start, self.interior, self.end]])


    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            return Arc2D(*[p.frame_mapping(frame, side, copy=True) for p in
                           [self.start, self.interior, self.end]])
        else:
            self.__init__(*[p.frame_mapping(frame, side, copy=True) for p in
                            [self.start, self.interior, self.end]])


    def SecondMomentArea(self, point):
        """
        Second moment area of part of disk
        """
        if self.angle2 < self.angle1:
            angle2 = self.angle2 + two_pi

        else:
            angle2 = self.angle2
        angle1 = self.angle1

        Ix = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                    math.sin(2 * angle1) - math.sin(2 * angle2)))
        Iy = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                    math.sin(2 * angle2) - math.sin(2 * angle1)))
        Ixy = self.radius ** 4 / 8 * (
                    math.cos(angle1) ** 2 - math.cos(angle2) ** 2)
        Ic = npy.array([[Ix, Ixy], [Ixy, Iy]])
        return geometry.Huygens2D(Ic, self.Area(), self.center, point)

    def Discretise(self, num=10):
        list_node = []
        if (self.angle1 < 0) and (self.angle2 > 0):
            delta_angle = -self.angle1 + self.angle2
        elif (self.angle1 > 0) and (self.angle2 < 0):
            delta_angle = (2 * npy.pi + self.angle2) - self.angle1
        else:
            delta_angle = self.angle2 - self.angle1
        for angle in npy.arange(self.angle1, self.angle1 + delta_angle,
                                delta_angle / (num * 1.)):
            list_node.append(Point2D(self.center + self.radius * Vector2D(
                (npy.cos(angle), npy.sin(angle)))))
        list_node.append(Point2D(self.center + self.radius * Vector2D((npy.cos(
            self.angle1 + delta_angle), npy.sin(self.angle1 + delta_angle)))))
        if list_node[0] == self.start:
            return list_node
        else:
            return list_node[::-1]

    def plot_data(self, marker=None, color='black', stroke_width=1, dash=False,
                  opacity=1):
        list_node = self.Discretise()
        data = []
        for nd in list_node:
            data.append({'x': nd.vector[0], 'y': nd.vector[1]})
        return {'type': 'arc',
                'cx': self.center.vector[0],
                'cy': self.center.vector[1],
                'data': data,
                'r': self.radius,
                'color': color,
                'opacity': opacity,
                'size': stroke_width,
                'dash': None,
                'marker': marker,
                'angle1': self.angle1,
                'angle2': self.angle2, }

    def copy(self):
        return Arc2D(self.start.copy(),
                     self.interior.copy(),
                     self.end.copy())

    def split(self, split_point: Point2D):
        raise NotImplementedError
        return [Arc2D(self.start, self.split_point)]

    def polygon_points(self, points_per_radian=10, min_x_density=None,
                       min_y_density=None):

        number_points = math.ceil(self.angle*points_per_radian)
        densities = []
        for d in [min_x_density, min_y_density]:
            if d:
                densities.append(d)
        if densities:
            number_points = max(number_points,
                                min(densities)*self.angle*self.radius)
        l = self.Length()
        return [self.PointAtCurvilinearAbscissa(i*l/number_points)\
                for i in range(number_points+1)]


class ArcEllipse2D(Primitive2D):
    """

    """

    def __init__(self, start, interior, end, center, major_dir, name='',
                 extra=None):
        self.start = start
        self.interior = interior
        self.end = end
        self.center = center
        self.extra = extra
        self.major_dir = major_dir
        self.minor_dir = self.major_dir.deterministic_unit_normal_vector()

        frame = Frame2D(self.center, self.major_dir, self.minor_dir)
        start_new, end_new = frame.NewCoordinates(
            self.start), frame.NewCoordinates(self.end)
        interior_new, center_new = frame.NewCoordinates(
            self.interior), frame.NewCoordinates(self.center)

        #### from : https://math.stackexchange.com/questions/339126/how-to-draw-an-ellipse-if-a-center-and-3-arbitrary-points-on-it-are-given
        def theta_A_B(s, i, e,
                      c):  # theta=angle d'inclinaison ellipse par rapport à horizontal(sens horaire),A=demi grd axe, B=demi petit axe
            xs, ys, xi, yi, xe, ye = s[0] - c[0], s[1] - c[1], i[0] - c[0], i[
                1] - c[1], e[0] - c[0], e[1] - c[1]
            A = npy.array(([xs ** 2, ys ** 2, 2 * xs * ys],
                           [xi ** 2, yi ** 2, 2 * xi * yi],
                           [xe ** 2, ye ** 2, 2 * xe * ye]))
            invA = npy.linalg.inv(A)
            One = npy.array(([1],
                             [1],
                             [1]))
            C = npy.dot(invA, One)  # matrice colonne de taille 3
            theta = 0
            c1 = C[0] + C[1]
            c2 = (C[1] - C[0]) / math.cos(2 * theta)
            gdaxe = math.sqrt((2 / (c1 - c2)))
            ptax = math.sqrt((2 / (c1 + c2)))
            return theta, gdaxe, ptax

        if start == end:
            extra_new = frame.NewCoordinates(self.extra)
            theta, A, B = theta_A_B(start_new, extra_new, interior_new,
                                    center_new)
        else:
            theta, A, B = theta_A_B(start_new, interior_new, end_new,
                                    center_new)

        self.Gradius = A
        self.Sradius = B
        self.theta = theta

        # Angle pour start
        u1, u2 = start_new.vector[0] / self.Gradius, start_new.vector[
            1] / self.Sradius
        angle1 = sin_cos_angle(u1, u2)
        # Angle pour end
        u3, u4 = end_new.vector[0] / self.Gradius, end_new.vector[
            1] / self.Sradius
        angle2 = sin_cos_angle(u3, u4)
        # Angle pour interior
        u5, u6 = interior_new.vector[0] / self.Gradius, interior_new.vector[
            1] / self.Sradius
        anglei = sin_cos_angle(u5, u6)

        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei + two_pi) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + two_pi

        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2 + two_pi) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + two_pi

        if clockwise_path > trigowise_path:
            self.is_trigo = True
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle = clockwise_path

        if self.start == self.end or self.angle == 0:
            self.angle = two_pi

        if self.is_trigo:  # sens trigo
            self.offset_angle = angle1
        else:
            self.offset_angle = angle2

        Primitive2D.__init__(self, name=name)

    def _get_points(self):
        return self.tessellation_points()

    points = property(_get_points)

    def tessellation_points(self, resolution_for_ellipse=40):
        number_points_tesselation = math.ceil(
            resolution_for_ellipse * abs(0.5 * self.angle / math.pi))

        frame2d = Frame2D(self.center, self.major_dir, self.minor_dir)

        tessellation_points_2D = [(Point2D((self.Gradius * math.cos(
            self.offset_angle + self.angle * i / (number_points_tesselation)),
                                            self.Sradius * math.sin(
                                                self.offset_angle + self.angle * i / (
                                                    number_points_tesselation)))))
                                  for i in
                                  range(number_points_tesselation + 1)]

        global_points = []
        for pt in tessellation_points_2D:
            global_points.append(frame2d.OldCoordinates(pt))

        return global_points

    def To3D(self, plane_origin, x, y):
        ps = self.start.To3D(plane_origin, x, y)
        pi = self.interior.To3D(plane_origin, x, y)
        pe = self.end.To3D(plane_origin, x, y)
        pc = self.center.To3D(plane_origin, x, y)
        if self.extra is None:
            pextra = None
        else:
            pextra = self.extra.To3D(plane_origin, x, y)
        if ps == pe:
            p3 = pextra
        else:
            p3 = pe
        plane = volmdlr.surfaces3d.Plane3D.from_3_points(ps, pi, p3)
        n = plane.normal
        major_dir = self.major_dir.To3D(plane_origin, x, y)
        major_dir.Normalize()

        return ArcEllipse3D(ps, pi, pe, pc, major_dir, normal=n,
                            name=self.name, extra=pextra)

    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = ax.figure

        self.interior.MPLPlot(ax=ax, color='m')
        self.start.MPLPlot(ax=ax, color='r')
        self.end.MPLPlot(ax=ax, color='b')
        self.center.MPLPlot(ax=ax, color='y')

        x = []
        y = []
        for px, py in self.tessellation_points():
            x.append(px)
            y.append(py)

        plt.plot(x, y, 'k')
        return fig, ax


class DisplayMesh(dc.DessiaObject):
    def __init__(self, points, triangles, edges=None, name=''):

        self.points = points
        self.triangles = triangles
        if edges is None:
            edges = []
        self.edges = edges
        self.name = name

    def __add__(self, other_mesh):
        new_points = self.points[:]
        new_point_index = {p: i for i, p in enumerate(self.points)}
        ip = len(new_points)
        for point in other_mesh.points:
            if not point in new_point_index:
                new_point_index[point] = ip
                ip += 1
                new_points.append(point)

        new_triangles = self.triangles.copy()
        for i1, i2, i3 in other_mesh.triangles:
            p1 = other_mesh.points[i1]
            p2 = other_mesh.points[i2]
            p3 = other_mesh.points[i3]
            new_triangles.append((new_point_index[p1],
                                  new_point_index[p2],
                                  new_point_index[p3]))

        return self.__class__(new_points, new_triangles)

    def plot(self, ax=None, numbering=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')

        for p in self.points:
            p.MPLPlot(ax=ax)

        for i1, i2, i3 in self.triangles:
            self._linesegment_class(self.points[i1], self.points[i2]).MPLPlot(
                ax=ax)
            self._linesegment_class(self.points[i2], self.points[i3]).MPLPlot(
                ax=ax)
            self._linesegment_class(self.points[i1], self.points[i3]).MPLPlot(
                ax=ax)

        for i, (i1, i2) in enumerate(self.edges):
            self._linesegment_class(self.points[i1], self.points[i2]).MPLPlot(
                ax=ax)
            if numbering:
                ax.text(*0.5*(self.points[i1]+self.points[i2]), 'edge {}'.format(i+1),
                        ha='center', va='center')

        return ax


class DisplayMesh2D(DisplayMesh):
    _linesegment_class = LineSegment2D
    _point_class = Point2D

    def __init__(self, points: List[Point2D],
                 triangles: List[Tuple[int, int, int]],
                 edges: List[Tuple[int, int]]=None,
                 name: str=''):
        DisplayMesh.__init__(self, points, triangles, edges, name=name)


class Wire2D(CompositePrimitive2D):
    """
    A collection of simple primitives, following each other making a wire
    """

    def __init__(self, primitives, name=''):
        CompositePrimitive2D.__init__(self, primitives, name)

    # TODO: method to check if it is a wire

    def Length(self):
        length = 0.
        for primitive in self.primitives:
            length += primitive.Length()
        return length

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa: float):
        length = 0.
        for primitive in self.primitives:
            primitive_length = primitive.Length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(
                    curvilinear_abscissa - length)
            length += primitive_length
        return ValueError

    def plot_data(self, name: str = '', fill=None, color='black',
                  stroke_width: float = 1, opacity: float = 1):
        plot_data = {}
        plot_data['name'] = name
        plot_data['type'] = 'wire'
        plot_data['plot_data'] = []
        for item in self.primitives:
            plot_data['plot_data'].append(item.plot_data(color=color,
                                                         stroke_width=stroke_width,
                                                         opacity=opacity))
        return plot_data

    def line_intersections(self, line: Line2D):
        """
        Returns a list of intersection in ther form of a tuple (point, primitive)
        of the wire primitives intersecting with the line
        """
        intersection_points = []
        for primitive in self.primitives:
            for p in primitive.line_intersections(line):
                intersection_points.append((p, primitive))
        return intersection_points

    def tesselation_points(self):
        points = []
        for p in self.primitives:
            points.extend(p.tessellation_points())
        return points

class Contour2D(Wire2D):
    """
    A collection of 2D primitives forming a closed wire2D
    TODO : CenterOfMass and SecondMomentArea should be changed accordingly to
    Area considering the triangle drawn by the arcs
    """
    _non_serializable_attributes = ['internal_arcs', 'external_arcs',
                                    'polygon', 'straight_line_contour_polygon']

    def __init__(self, primitives, name=''):
        Wire2D.__init__(self, primitives, name)
        self._utd_analysis = False
        self.tessel_points = self.clean_points()

    def _primitives_analysis(self):
        """
        An internal arc is an arc that has his interior point inside the polygon
        """
        arcs = []
        internal_arcs = []
        external_arcs = []
        points_polygon = []
        points_straight_line_contour = []
        for primitive in self.primitives:
            if primitive.__class__.__name__ == 'LineSegment2D':
                points_polygon.append(primitive.start)
                points_straight_line_contour.append(primitive.start)
                points_straight_line_contour.append(primitive.end)
            elif primitive.__class__.__name__ == 'Arc2D':
                points_polygon.append(primitive.start)
                points_polygon.append(primitive.center)

                # points_polygon.append(primitive.end)
                arcs.append(primitive)
            elif primitive.__class__.__name__ == 'Circle2D':
                print(self.primitives)
                raise ValueError(
                    'Circle2D primitives should not be inserted in a contour, as a circle is already a contour. Use directcly the circle')
                # return None
            elif primitive.__class__.__name__ == 'OpenedRoundedLineSegments2D':
                for prim in primitive.primitives:
                    if prim.__class__.__name__ == 'LineSegment2D':
                        points_polygon.extend(prim.points)
                        points_straight_line_contour.extend(prim.points)
                    elif prim.__class__.__name__ == 'Arc2D':
                        #                points_polygon.append(primitive.center)
                        points_polygon.append(prim.start)
                        points_polygon.append(prim.end)
                        arcs.append(prim)
            else:
                raise NotImplementedError(
                    'primitive of type {} is not handled'.format(primitive))

        # points_polygon = list(set(points_polygon))
        polygon = Polygon2D(points_polygon)
        points_straight_line_contour = list(set(points_straight_line_contour))
        straight_line_contour_polygon = Polygon2D(points_straight_line_contour)

        for arc in arcs:
            if polygon.PointBelongs(arc.interior):
                internal_arcs.append(arc)
            else:
                external_arcs.append(arc)

        return internal_arcs, external_arcs, polygon, straight_line_contour_polygon

    def _get_internal_arcs(self):
        if not self._utd_analysis:
            (self._internal_arcs, self._external_arcs,
             self._polygon,
             self._straight_line_contour_polygon) = self._primitives_analysis()
            self._utd_analysis = True
        return self._internal_arcs

    internal_arcs = property(_get_internal_arcs)

    def _get_external_arcs(self):
        if not self._utd_analysis:
            (self._internal_arcs, self._external_arcs,
             self._polygon,
             self._straight_line_contour_polygon) = self._primitives_analysis()
            self._utd_analysis = True
        return self._external_arcs

    external_arcs = property(_get_external_arcs)

    def _get_polygon(self):
        if not self._utd_analysis:
            (self._internal_arcs, self._external_arcs,
             self._polygon,
             self._straight_line_contour_polygon) = self._primitives_analysis()
            self._utd_analysis = True
        return self._polygon

    polygon = property(_get_polygon)

    def _get_straight_line_contour_polygon(self):
        if not self._utd_analysis:
            (self._internal_arcs, self._external_arcs,
             self._polygon,
             self._straight_line_contour_polygon) = self._primitives_analysis()
            self._utd_analysis = True
        return self._straight_line_contour_polygon

    straight_line_contour_polygon = property(
        _get_straight_line_contour_polygon)

    def point_belongs(self, point):
        for arc in self.internal_arcs:
            if arc.point_belongs(point):
                return False
        if self.polygon.PointBelongs(point):
            return True
        for arc in self.external_arcs:
            if arc.point_belongs(point):
                return True
        return False

    def point_distance(self, point):
        min_distance = self.primitives[0].point_distance(point)
        for primitive in self.primitives[1:]:
            distance = primitive.point_distance(point)
            if distance < min_distance:
                min_distance = distance
        return min_distance

    def bounding_points(self):
        points = self.straight_line_contour_polygon.points[:]
        for arc in self.internal_arcs + self.external_arcs:
            points.extend(arc.tessellation_points())
        xmin = min([p[0] for p in points])
        xmax = max([p[0] for p in points])
        ymin = min([p[1] for p in points])
        ymax = max([p[1] for p in points])
        return (Point2D((xmin, ymin)), Point2D((xmax, ymax)))

    def to_2d(self, plane3d, name=None):
        primitives3D = [p.to_2d(plane3d) for p in self.primitives]
        return Contour3D(edges=primitives3D, name=name)

    def To3D(self, plane_origin, x, y, name=None):
        if name is None:
            name = '3D of {}'.format(self.name)
        primitives3D = [p.To3D(plane_origin, x, y) for p in self.primitives]
        return Contour3D(primitives=primitives3D, name=name)

    def Area(self):
        if len(self.primitives) == 1:
            return self.primitives[0].Area()

        A = self.polygon.Area()

        for arc in self.internal_arcs:
            triangle = Polygon2D([arc.start, arc.center, arc.end])
            A = A - arc.Area() + triangle.Area()
        for arc in self.external_arcs:
            triangle = Polygon2D([arc.start, arc.center, arc.end])
            A = A + arc.Area() - triangle.Area()

        return A

    def CenterOfMass(self):
        if len(self.primitives) == 1:
            return self.primitives[0].CenterOfMass()

        area = self.polygon.Area()
        if area > 0.:
            c = area * self.polygon.CenterOfMass()
        else:
            c = O2D

        for arc in self.internal_arcs:
            arc_area = arc.Area()
            c -= arc_area * arc.CenterOfMass()
            area -= arc_area
        for arc in self.external_arcs:
            arc_area = arc.Area()
            c += arc_area * arc.CenterOfMass()
            area += arc_area
        if area != 0:
            return c / area
        else:
            return False

    def SecondMomentArea(self, point):
        if len(self.primitives) == 1:
            return self.primitives[0].SecondMomentArea(point)

        A = self.polygon.SecondMomentArea(point)
        for arc in self.internal_arcs:
            A -= arc.SecondMomentArea(point)
        for arc in self.external_arcs:
            A += arc.SecondMomentArea(point)

        return A

    def plot_data(self, name='', fill=None, marker=None, color='black',
                  stroke_width=1, dash=False, opacity=1):

        plot_data = {}
        plot_data['fill'] = fill
        plot_data['name'] = name
        plot_data['type'] = 'contour'
        plot_data['plot_data'] = []
        for item in self.primitives:
            plot_data['plot_data'].append(item.plot_data(color=color,
                                                         stroke_width=stroke_width,
                                                         opacity=opacity))
        return plot_data

    def copy(self):
        primitives_copy = []
        for primitive in self.primitives:
            primitives_copy.append(primitive.copy())
        return Contour2D(primitives_copy)

    def average_center_point(self):
        nb = len(self.tessel_points)
        x = npy.sum([p[0] for p in self.tessel_points]) / nb
        y = npy.sum([p[1] for p in self.tessel_points]) / nb
        return Point2D((x, y))

    def clean_points(self):
        """
        This method is copy from Contour3D, if changes are done there or here,
        please change both method
        Be aware about primitives = 2D, edges = 3D
        """
        if hasattr(self.primitives[0], 'endpoints'):
            points = self.primitives[0].endpoints[:]
        else:
            points = self.primitives[0].tessellation_points()
        for primitive in self.primitives[1:]:
            if hasattr(primitive, 'endpoints'):
                points_to_add = primitive.endpoints[:]
            else:
                points_to_add = primitive.tessellation_points()
            if points[0] == points[
                -1]:  # Dans le cas où le (dernier) edge relie deux fois le même point
                points.extend(points_to_add[::-1])

            elif points_to_add[0] == points[-1]:
                points.extend(points_to_add[1:])
            elif points_to_add[-1] == points[-1]:
                points.extend(points_to_add[-2::-1])
            elif points_to_add[0] == points[0]:
                points = points[::-1]
                points.extend(points_to_add[1:])
            elif points_to_add[-1] == points[0]:
                points = points[::-1]
                points.extend(points_to_add[-2::-1])
            else:
                d1, d2 = (points_to_add[0] - points[0]).Norm(), (
                            points_to_add[0] - points[-1]).Norm()
                d3, d4 = (points_to_add[-1] - points[0]).Norm(), (
                            points_to_add[-1] - points[-1]).Norm()
                if math.isclose(d2, 0, abs_tol=1e-3):
                    points.extend(points_to_add[1:])
                elif math.isclose(d4, 0, abs_tol=1e-3):
                    points.extend(points_to_add[-2::-1])
                elif math.isclose(d1, 0, abs_tol=1e-3):
                    points = points[::-1]
                    points.extend(points_to_add[1:])
                elif math.isclose(d3, 0, abs_tol=1e-3):
                    points = points[::-1]
                    points.extend(points_to_add[-2::-1])

        if len(points) > 1:
            if points[0] == points[-1]:
                points.pop()
        return points

    def bounding_rectangle(self):
        # bounding rectangle
        tp = self.tesselation_points()
        xmin = tp[0][0]
        xmax = tp[0][0]
        ymin = tp[0][1]
        ymax = tp[0][1]
        for point in tp[1:]:
            xmin = min(point[0], xmin)
            xmax = max(point[0], xmax)
            ymin = min(point[1], ymin)
            ymax = max(point[1], ymax)
        return xmin, xmax, ymin, ymax


    def random_point_inside(self):
        xmin, xmax, ymin, ymax = self.bounding_rectangle()
        for i in range(1000):
            p = Point2D.random(xmin, xmax, ymin, ymax)
            if self.point_belongs(p):
                return p
    # def line_intersections(self, line:Line2D) -> List[Tuple[Point2D, Primitive2D]]:
    #     """
    #     Returns a list of points and lines of intersection with the contour
    #     """
    #     intersection_points = Wire2D.line_intersections(self, line)
    #     if not intersection_points:
    #         return []
    #     elif len(intersection_points) == 2:
    #         return [LineSegment2D(*intersection_points)]
    #     else:
    #         raise NotImplementedError('Non convex contour not supported yet')

    def cut_by_line(self, line):
        intersections = self.line_intersections(line)
        if not intersections:
            return [self]

        if len(intersections) < 2:
            return [self]
        elif len(intersections) == 2:
            if intersections[0][0].__class__.__name__ == 'Point2D' and \
                    intersections[1][0].__class__.__name__ == 'Point2D':
                ip1, ip2 = sorted([self.primitives.index(intersections[0][1]),
                                   self.primitives.index(intersections[1][1])])

                sp11, sp12 = intersections[0][1].split(intersections[0][0])
                sp21, sp22 = intersections[1][1].split(intersections[1][0])

                primitives1 = self.primitives[:ip1]
                primitives1.append(sp11)
                primitives1.append(LineSegment2D(intersections[0][0],
                                                 intersections[1][0]))
                primitives1.append(sp22)
                primitives1.extend(self.primitives[ip2 + 1:])

                primitives2 = self.primitives[ip1 + 1:ip2]
                primitives2.append(sp21)
                primitives2.append(LineSegment2D(intersections[1][0],
                                                 intersections[0][0]))
                primitives2.append(sp12)

                return Contour2D(primitives1), Contour2D(primitives2)

            else:
                print(intersections)
                raise NotImplementedError(
                    'Non convex contour not supported yet')

        raise NotImplementedError(
            '{} intersections not supported yet'.format(len(intersections)))

    def simple_triangulation(self):
        lpp = len(self.polygon.points)
        if lpp == 3:
            return self.polygon.points, [(0, 1, 2)]
        elif lpp == 4:
            return self.polygon.points, [(0, 1, 2), (0, 2, 3)]

        # Use delaunay triangulation
        tri = Delaunay([p.vector for p in self.polygon.points])
        indices = tri.simplices
        return self.polygon.points, tri.simplices

    def split_regularly(self, n):
        """
        Split in n slices
        """
        xmin, xmax, ymin, ymax = self.bounding_rectangle()
        cutted_contours = []
        iteration_contours = [self]
        for i in range(n - 1):
            # print(i)
            xi = xmin + (i + 1) * (xmax - xmin) / n
            # print(xi)
            cut_line = Line2D(Point2D((xi, 0)), Point2D((xi, 1)))

            iteration_contours2 = []
            for c in iteration_contours:
                sc = c.cut_by_line(cut_line)
                lsc = len(sc)
                # print('lsc', lsc)
                if lsc == 1:
                    cutted_contours.append(c)
                else:
                    iteration_contours2.extend(sc)

            iteration_contours = iteration_contours2[:]
        cutted_contours.extend(iteration_contours)
        return cutted_contours

    def triangulation(self):
        return self.grid_triangulation(number_points_x=20,
                                       number_points_y=20)

    def polygonization(self, min_x_density=None, min_y_density=None):
        # points = self.primitives[0].polygon_points()
        points = []
        for primitive in self.primitives:
            points.extend(primitive.polygon_points(min_x_density=min_x_density,
                                                   min_y_density=min_y_density)[1:])

        return Polygon2D(points)

    def grid_triangulation(self, x_density: float = None,
                           y_density: float = None,
                           min_points_x: int = 20,
                           min_points_y: int = 20,
                           number_points_x: int = None,
                           number_points_y: int = None):
        """
        Use a n by m grid to triangulize the contour
        """
        xmin, xmax, ymin, ymax = self.bounding_rectangle()
        dx = xmax - xmin
        dy = ymax - ymin
        if number_points_x is None:
            n = max(math.ceil(x_density * dx), min_points_x)
        else:
            n = number_points_x
        if number_points_y is None:
            m = max(math.ceil(y_density * dy), min_points_y)
        else:
            m = number_points_y

        x = [xmin + i * dx / n for i in range(n + 1)]
        y = [ymin + i * dy / m for i in range(m + 1)]

        point_is_inside = {}
        point_index = {}
        ip = 0
        points = []
        triangles = []
        for xi in x:
            for yi in y:
                p = Point2D((xi, yi))
                if self.point_belongs(p):
                    point_index[p] = ip
                    points.append(p)
                    ip += 1

        for i in range(n):
            for j in range(m):
                p1 = Point2D((x[i], y[j]))
                p2 = Point2D((x[i + 1], y[j]))
                p3 = Point2D((x[i + 1], y[j + 1]))
                p4 = Point2D((x[i], y[j + 1]))
                points_in = []
                for p in [p1, p2, p3, p4]:
                    if p in point_index:
                        points_in.append(p)
                if len(points_in) == 4:
                    triangles.append(
                        [point_index[p1], point_index[p2], point_index[p3]])
                    triangles.append(
                        [point_index[p1], point_index[p3], point_index[p4]])

                elif len(points_in) == 3:
                    triangles.append([point_index[p] for p in points_in])

        return DisplayMesh2D(points, triangles)


class Circle2D(Contour2D):
    _non_serializable_attributes = ['internal_arcs', 'external_arcs',
                                    'polygon', 'straight_line_contour_polygon',
                                    'primitives', 'basis_primitives']

    def __init__(self, center: Point2D, radius: float, name: str = ''):
        self.center = center
        self.radius = radius
        self.angle = two_pi
        self.utd_geo_points = False

        self.points = self.tessellation_points()

        Contour2D.__init__(self, [self], name=name)  # !!! this is dangerous

    def __hash__(self):
        return int(round(1e6 * (self.center.vector[0] + self.center.vector[
            1] + self.radius)))

    def __eq__(self, other_circle):
        return math.isclose(self.center.vector[0],
                            other_circle.center.vector[0], abs_tol=1e-06) \
               and math.isclose(self.center.vector[1],
                                other_circle.center.vector[1], abs_tol=1e-06) \
               and math.isclose(self.radius, other_circle.radius,
                                abs_tol=1e-06)

    def tessellation_points(self, resolution=40):
        return [self.center + self.radius * math.cos(teta) * Vector2D(
            (1, 0)) + self.radius * math.sin(teta) * Vector2D((0, 1)) \
                for teta in npy.linspace(0, two_pi, resolution + 1)][:-1]

    def point_belongs(self, point, tolerance=1e-9):
        return point.point_distance(self.center) <= self.radius + tolerance

    def line_intersections(self, line):
        V = Vector2D((line.points[1] - line.points[0]).vector)
        Q = Vector2D(self.center.vector)
        P1 = Vector2D(line.points[0].vector)

        a = V.Dot(V)
        b = 2 * V.Dot(P1 - Q)
        c = P1.Dot(P1) + Q.Dot(Q) - 2 * P1.Dot(Q) - self.radius ** 2

        disc = b ** 2 - 4 * a * c
        if disc < 0:
            return []

        if math.isclose(disc, 0, abs_tol=1e-8):
            t = -b / (2 * a)
            if line.__class__ is Line2D:
                return [Point2D((P1 + t * V).vector)]
            else:
                if 0 <= t <= 1:
                    return [Point2D((P1 + t * V).vector)]
                else:
                    return []

        sqrt_disc = math.sqrt(disc)
        t1 = (-b + sqrt_disc) / (2 * a)
        t2 = (-b - sqrt_disc) / (2 * a)
        if line.__class__ is Line2D:
            return [Point2D((P1 + t1 * V).vector),
                    Point2D((P1 + t2 * V).vector)]
        else:
            if not (0 <= t1 <= 1 or 0 <= t2 <= 1):
                return []
            elif 0 <= t1 <= 1 and not 0 <= t2 <= 1:
                return [Point2D((P1 + t1 * V).vector)]
            elif not 0 <= t1 <= 1 and 0 <= t2 <= 1:
                return [Point2D((P1 + t2 * V).vector)]
            else:
                [Point2D((P1 + t1 * V).vector), Point2D((P1 + t2 * V).vector)]

    def Length(self):
        return two_pi * self.radius

    def MPLPlot(self, ax=None, linestyle='-', color='k', linewidth=1):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        # else:
        #     fig = ax.figure

        pc = self.center.vector
        if self.radius > 0:
            ax.add_patch(Arc(pc,
                             2 * self.radius,
                             2 * self.radius,
                             angle=0,
                             theta1=0,
                             theta2=360,
                             color=color,
                             linestyle=linestyle,
                             linewidth=linewidth))
        return ax

    def To3D(self, plane_origin, x, y):
        normal = x.Cross(y)
        center3d = self.center.To3D(plane_origin, x, y)
        return Circle3D(Frame3D(center3d, x, y, normal),
                        self.radius, self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Circle2D(self.center.Rotation(center, angle, copy=True),
                            self.radius)
        else:
            self.center.Rotation(center, angle, copy=False)
            self.utd_geo_points = False

    def Translation(self, offset, copy=True):
        if copy:
            return Circle2D(self.center.Translation(offset, copy=True),
                            self.radius)
        else:
            self.center.Translation(offset, copy=False)
            self.utd_geo_points = False

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return Circle2D(frame.OldCoordinates(self.center), self.radius)
            else:
                self.center = frame.OldCoordinates(self.center)
        if side == 'new':
            if copy:
                return Circle2D(frame.NewCoordinates(self.center), self.radius)
            else:
                self.points = frame.NewCoordinates(self.center)

    def Area(self):
        return math.pi * self.radius ** 2

    def SecondMomentArea(self, point):
        """
        Second moment area of part of disk
        """
        I = math.pi * self.radius ** 4 / 4
        Ic = npy.array([[I, 0], [0, I]])
        return geometry.Huygens2D(Ic, self.Area(), self.center, point)

    def CenterOfMass(self):
        return self.center

    def point_symmetric(self, point):
        center = 2 * point - self.center
        return Circle2D(center, self.radius)

    def plot_data(self, marker=None, color='black', stroke_width=1, opacity=1,
                  fill=None):
        return {'type': 'circle',
                'cx': self.center.vector[0],
                'cy': self.center.vector[1],
                'r': self.radius,
                'color': color,
                'opacity': opacity,
                'size': stroke_width,
                'dash': None,
                'fill': fill}

    def copy(self):
        return Circle2D(self.center.copy(), self.radius)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        start = self.center + self.radius * X3D
        return start.Rotation(self.center,
                              curvilinear_abscissa / self.radius)

    def triangulation(self, n=35):
        l = self.Length()
        points = [self.PointAtCurvilinearAbscissa(l * i / n) for i in range(n)]
        points.append(self.center)
        triangles = [(i, i + 1, n) for i in range(n - 1)] + [(n - 1, 0, n)]
        return DisplayMesh(points, triangles)

    def polygon_points(self, points_per_radian=10, min_x_density=None,
                       min_y_density=None):
        return Arc2D.polygon_points(self, points_per_radian=points_per_radian,
                                    min_x_density=min_x_density,
                                    min_y_density=min_y_density)

class Polygon2D(Contour2D):

    def __init__(self, points, name=''):
        self.points = points
        self.line_segments = self._LineSegments()

        Contour2D.__init__(self, self.line_segments, name)

    def copy(self):
        points = [p.copy() for p in self.points]
        return Polygon2D(points, self.name)

    def __hash__(self):
        return sum([hash(p) for p in self.points])

    def __eq__(self, other_):
        equal = True
        for point, other_point in zip(self.points, other_.points):
            equal = (equal and point == other_point)
        return equal

    def Area(self):

        x = [point.vector[0] for point in self.points]
        y = [point.vector[1] for point in self.points]

        return 0.5 * npy.abs(
            npy.dot(x, npy.roll(y, 1)) - npy.dot(y, npy.roll(x, 1)))

    def CenterOfMass(self):

        x = [point.vector[0] for point in self.points]
        y = [point.vector[1] for point in self.points]

        xi_xi1 = x + npy.roll(x, -1)
        yi_yi1 = y + npy.roll(y, -1)
        xi_yi1 = npy.multiply(x, npy.roll(y, -1))
        xi1_yi = npy.multiply(npy.roll(x, -1), y)

        a = 0.5 * npy.sum(xi_yi1 - xi1_yi)  # signed area!
        #        a=self.Area()
        if not math.isclose(a, 0, abs_tol=1e-08):
            cx = npy.sum(npy.multiply(xi_xi1, (xi_yi1 - xi1_yi))) / 6. / a
            cy = npy.sum(npy.multiply(yi_yi1, (xi_yi1 - xi1_yi))) / 6. / a
            return Point2D((cx, cy))

        else:
            raise NotImplementedError

    def PointBelongs(self, point):
        """
        Ray casting algorithm copied from internet...
        """
        return PolygonPointBelongs(point.vector,
                                   [p.vector for p in self.points])

    def SecondMomentArea(self, point):
        Ix, Iy, Ixy = 0, 0, 0
        for pi, pj in zip(self.points, self.points[1:] + [self.points[0]]):
            xi, yi = (pi - point).vector
            xj, yj = (pj - point).vector
            Ix += (yi ** 2 + yi * yj + yj ** 2) * (xi * yj - xj * yi)
            Iy += (xi ** 2 + xi * xj + xj ** 2) * (xi * yj - xj * yi)
            Ixy += (xi * yj + 2 * xi * yi + 2 * xj * yj + xj * yi) * (
                        xi * yj - xj * yi)
        if Ix < 0:
            Ix = - Ix
            Iy = - Iy
            Ixy = - Ixy
        return npy.array([[Ix / 12., Ixy / 24.], [Ixy / 24., Iy / 12.]])

    def _LineSegments(self):
        lines = []
        for p1, p2 in zip(self.points, self.points[1:] + [self.points[0]]):
            lines.append(LineSegment2D(p1, p2))
        return lines

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Polygon2D(
                [p.Rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Polygon2D(
                [p.Translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset, copy=False)

    def PointBorderDistance(self, point, return_other_point=False):
        """
        Compute the distance to the border distance of polygon
        Output is always positive, even if the point belongs to the polygon
        """
        d_min, other_point_min = self.line_segments[0].point_distance(point,
                                                                      return_other_point=True)
        for line in self.line_segments[1:]:
            d, other_point = line.point_distance(point,
                                                 return_other_point=True)
            if d < d_min:
                d_min = d
                other_point_min = other_point
        if return_other_point:
            return d_min, other_point_min
        return d_min

    def SelfIntersect(self):
        epsilon = 0
        # BENTLEY-OTTMANN ALGORITHM
        # Sort the points along ascending x for the Sweep Line method
        sorted_index = sorted(range(len(self.points)), key=lambda p: (
        self.points[p][0], self.points[p][1]))
        nb = len(sorted_index)
        segments = []
        deleted = []

        while len(
                sorted_index) != 0:  # While all the points haven't been swept
            # Stock the segments between 2 consecutive edges
            # Ex: for the ABCDE polygon, if Sweep Line is on C, the segments
            #   will be (C,B) and (C,D)
            if sorted_index[0] - 1 < 0:
                segments.append((sorted_index[0], nb - 1))
            else:
                segments.append((sorted_index[0], sorted_index[0] - 1))
            if sorted_index[0] >= len(self.points) - 1:
                segments.append((sorted_index[0], 0))
            else:
                segments.append((sorted_index[0], sorted_index[0] + 1))

            # Once two edges linked by a segment have been swept, delete the
            # segment from the list
            to_del = []
            for index in deleted:
                if abs(index - sorted_index[0]) == 1 or abs(
                        index - sorted_index[0]) == nb - 1:
                    to_del.append((index, sorted_index[0]))
                    to_del.append((sorted_index[0], index))

            # Keep track of which edges have been swept
            deleted.append(sorted_index[0])
            sorted_index.pop(0)

            # Delete the segments that have just been swept
            index_to_del = []
            for i, segment in enumerate(segments):
                for seg_to_del in to_del:
                    if segment == seg_to_del:
                        index_to_del.append(i)
            for index in index_to_del[::-1]:
                segments.pop(index)

            # Checks if two segments are intersecting each other, returns True
            # if yes, otherwise the algorithm continues at WHILE
            for segment1 in segments:
                for segment2 in segments:
                    if segment1[0] != segment2[0] and segment1[1] != segment2[
                        1] and segment1[0] != segment2[1] and segment1[1] != \
                            segment2[0]:

                        line1 = LineSegment2D(
                            Point2D(self.points[segment1[0]]),
                            Point2D(self.points[segment1[1]]))
                        line2 = LineSegment2D(
                            Point2D(self.points[segment2[0]]),
                            Point2D(self.points[segment2[1]]))

                        p, a, b = Point2D.LinesIntersection(line1, line2, True)

                        if p is not None:
                            if a >= 0 + epsilon and a <= 1 - epsilon and b >= 0 + epsilon and b <= 1 - epsilon:
                                return True, line1, line2

        return False, None, None


    def plot_data(self, marker=None, color='black', stroke_width=1, opacity=1):
        data = []
        for nd in self.points:
            data.append({'x': nd.vector[0], 'y': nd.vector[1]})
        return {'type': 'wire',
                'data': data,
                'color': color,
                'size': stroke_width,
                'dash': None,
                'marker': marker,
                'opacity': opacity}

    @classmethod
    def points_convex_hull(cls, points):
        ymax, pos_ymax = max_pos([pt.vector[1] for pt in points])
        point_start = points[pos_ymax]
        hull, thetac = [point_start], 0  # thetac is the current theta

        barycenter = points[0]
        for pt in points[1:]:
            barycenter += pt
        barycenter = barycenter / (len(points))
        # second point of hull
        theta = []
        remaining_points = points
        del remaining_points[pos_ymax]

        vec1 = point_start - barycenter
        for pt in remaining_points:
            vec2 = pt - point_start
            theta_i = -clockwise_angle(vec1, vec2)
            theta.append(theta_i)

        min_theta, posmin_theta = min_pos(theta)
        thetac += min_theta
        next_point = remaining_points[posmin_theta]
        hull.append(next_point)
        del remaining_points[posmin_theta]
        # Adding first point to close the loop at the end
        remaining_points.append(hull[0])

        while next_point != point_start:
            vec1 = next_point - barycenter
            theta = []
            for pt in remaining_points:
                vec2 = pt - next_point
                theta_i = -clockwise_angle(vec1, vec2)
                theta.append(theta_i)

            min_theta, posmin_theta = min_pos(theta)
            thetac += min_theta
            next_point = remaining_points[posmin_theta]
            hull.append(next_point)
            del remaining_points[posmin_theta]

        hull.pop()

        return cls(hull)

    def MPLPlot(self, ax=None, color='k',
                plot_points=False, point_numbering=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')

        for ls in self.line_segments:
            ls.MPLPlot(ax=ax ,color=color)

        if plot_points or point_numbering:
            for point in self.points:
                point.MPLPlot(ax=ax, color=color)

        if point_numbering:
            for ip, point in enumerate(self.points):
                ax.text(*point, 'point {}'.format(ip+1),
                        ha='center', va='top')

        ax.margins(0.1)
        plt.show()

        return ax

class Surface2D(Primitive2D):
    """
    A surface bounded by an outer contour
    """
    def __init__(self, outer_contour: Contour2D,
                 inner_contours: List[Contour2D],
                 name:str='name'):
        self.outer_contour = outer_contour
        self.inner_contours = inner_contours

        Primitive2D.__init__(self, name=name)

    def triangulation(self, min_x_density=None, min_y_density=None):
        outer_polygon = self.outer_contour.polygonization(min_x_density=15, min_y_density=12)
        # ax2 = outer_polygon.MPLPlot(color='r', point_numbering=True)
        points = outer_polygon.points
        vertices = [p.vector for p in points]
        n = len(outer_polygon.points)
        segments = [(i, i+1) for i in range(n-1)]
        segments.append((n-1, 0))
        point_index = {p:i for i,p in enumerate(points)}
        holes = []

        for inner_contour in self.inner_contours:
            inner_polygon = inner_contour.polygonization()
            # inner_polygon.MPLPlot(ax=ax2)
            for point in inner_polygon.points:
                if not point in point_index:
                    points.append(point)
                    vertices.append(point.vector)
                    point_index[point] = n
                    n += 1
            for point1, point2 in zip(inner_polygon.points[:-1],
                                      inner_polygon.points[1:]):
                segments.append((point_index[point1],
                                 point_index[point2]))
            segments.append((point_index[inner_polygon.points[-1]],
                             point_index[inner_polygon.points[0]]))
            holes.append(inner_contour.random_point_inside().vector)


        tri = {'vertices': npy.array(vertices).reshape((-1, 2)),
               'segments': npy.array(segments).reshape((-1, 2)),
               }
        if holes:
            tri['holes'] = npy.array(holes).reshape((-1, 2))

        t = triangle.triangulate(tri, 'p')
        triangles = t['triangles'].tolist()

        return DisplayMesh2D(points, triangles=triangles, edges=None)

    def plot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        self.outer_contour.MPLPlot(ax=ax)
        for inner_contour in self.inner_contours:
            inner_contour.MPLPlot(ax=ax)

        ax.set_aspect('equal')
        ax.margins(0.1)
        return ax


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









class BSplineExtrusion(Primitive3D):

    def __init__(self, obj, vectorextru, name=''):
        self.obj = obj
        vectorextru.Normalize()
        self.vectorextru = vectorextru
        if obj.__class__ is Ellipse3D:
            self.points = obj.tessel_points
            # self.surface = obj.points
        else:
            self.points = obj.points

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]
        if object_dict[arguments[1]].__class__ is Ellipse3D:
            ell = object_dict[arguments[1]]
            vectextru = -object_dict[arguments[2]]
            return cls(ell, vectextru, name)

        elif object_dict[arguments[1]].__class__ is BSplineCurve3D:
            bsplinecurve = object_dict[arguments[1]]
            vectextru = object_dict[arguments[2]]
            return cls(bsplinecurve, vectextru, name)
        else:
            # surface = BSpline.Surface()
            # surface.degree_u = degree_u
            # surface.degree_v = degree_v
            # if weights is None:
            #     P = [(control_points[i][0], control_points[i][1], control_points[i][2]) for i in range(len(control_points))]
            #     surface.set_ctrlpts(P, nb_u, nb_v)
            # else:
            #     Pw = [(control_points[i][0]*weights[i], control_points[i][1]*weights[i], control_points[i][2]*weights[i], weights[i]) for i in range(len(control_points))]
            #     surface.set_ctrlpts(Pw, nb_u, nb_v)
            # knot_vector_u = []
            # for i, u_knot in enumerate(u_knots):
            #     knot_vector_u.extend([u_knot]*u_multiplicities[i])
            # knot_vector_v = []
            # for i, v_knot in enumerate(v_knots):
            #     knot_vector_v.extend([v_knot]*v_multiplicities[i])
            # surface.knotvector_u = knot_vector_u
            # surface.knotvector_v = knot_vector_v
            # surface.delta = 0.05
            # surface_points = surface.evalpts

            # self.surface = surface
            # self.points = [Point3D((p[0], p[1], p[2])) for p in surface_points]
            raise NotImplementedError  ## a adapter pour les bpsline

    # @classmethod
    # def from_step(cls, arguments, object_dict):
    #     name = arguments[0][1:-1]

    # degree_u = int(arguments[1])
    # degree_v = int(arguments[2])
    # points_sets = arguments[3][1:-1].split("),")
    # points_sets = [elem+")" for elem in points_sets[:-1]]+[points_sets[-1]]
    # control_points = []
    # for points_set in points_sets:
    #     points = [object_dict[int(i[1:])] for i in points_set[1:-1].split(",")]
    #     nb_v = len(points)
    #     control_points.extend(points)
    # nb_u = int(len(control_points) / nb_v)
    # surface_form = arguments[4]
    # if arguments[5] == '.F.':
    #     u_closed = False
    # elif arguments[5] == '.T.':
    #     u_closed = True
    # else:
    #     raise ValueError
    # if arguments[6] == '.F.':
    #     v_closed = False
    # elif arguments[6] == '.T.':
    #     v_closed = True
    # else:
    #     raise ValueError
    # self_intersect = arguments[7]
    # u_multiplicities = [int(i) for i in arguments[8][1:-1].split(",")]
    # v_multiplicities = [int(i) for i in arguments[9][1:-1].split(",")]
    # u_knots = [float(i) for i in arguments[10][1:-1].split(",")]
    # v_knots = [float(i) for i in arguments[11][1:-1].split(",")]
    # knot_spec = arguments[12]

    # if 13 in range(len(arguments)):
    #     weight_data = [float(i) for i in arguments[13][1:-1].replace("(", "").replace(")", "").split(",")]
    # else:
    #     weight_data = None

    # return cls(degree_u, degree_v, control_points, nb_u, nb_v, u_multiplicities, v_multiplicities, u_knots, v_knots, weight_data, name)


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
        # basis_primitives = []
        # for primitive in primitives:
        #     if hasattr(primitive, 'basis_primitives'):
        #         basis_primitives.extend(primitive.primitives)
        #     else:
        #         basis_primitives.append(primitive)

        Primitive3D.__init__(self, name=name)

    # def __eq__(self, other_):
    #     equal = True
    #     for primitive, other_primitive in zip(self.primitives, other_.primitives):
    #         equal = (equal and primitive == other_primitive)
    #     return equal

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

        # ax.set_aspect('equal')

        return ax


class Wire3D(CompositePrimitive3D):
    """
    A collection of simple primitives, following each other making a wire
    """

    def __init__(self, primitives, name=''):
        CompositePrimitive3D.__init__(self, primitives, name)

    def Length(self):
        length = 0.
        for primitive in self.primitives:
            length += primitive.Length()
        return length

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        length = 0.
        for primitive in self.primitives:
            primitive_length = primitive.Length()
            if length + primitive_length >= curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(
                    curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError

    # TODO: method to check if it is a wire
    def FreeCADExport(self, ip):
        name = 'primitive' + str(ip)

        s = 'E = []\n'
        for ip, primitive in enumerate(self.primitives):
            s += primitive.FreeCADExport('L{}'.format(ip))
            s += 'E.append(Part.Edge(L{}))\n'.format(ip)
        s += '{} = Part.Wire(E[:])\n'.format(name)

        return s

    def frame_mapping(self, frame, side, copy=True):
        new_wire = []
        if side == 'new':
            if copy:
                for primitive in self.primitives:
                    new_wire.append(primitive.frame_mapping(frame, side, copy))
                return Wire3D(new_wire)
            else:
                for primitive in self.primitives:
                    primitive.frame_mapping(frame, side, copy=False)

        if side == 'old':
            if copy:
                for primitive in self.primitives:
                    new_wire.append(primitive.frame_mapping(frame, side, copy))
                return Wire3D(new_wire)
            else:
                for primitive in self.primitives:
                    primitive.frame_mapping(frame, side, copy=False)

    def minimum_distance(self, wire2):
        distance = []
        for element in self.primitives:
            for element2 in wire2.primitives:
                distance.append(element.minimum_distance(element2))

        return min(distance)

    def copy(self):
        primitives_copy = []
        for primitive in self.primitives:
            primitives_copy.append(primitive.copy())
        return Wire3D(primitives_copy)


class LineSegment3D(LineSegment):
    """
    Define a line segment limited by two points
    """

    def __init__(self, start, end, name=''):
        LineSegment.__init__(self, start, end, name)
        self.bounding_box = self._bounding_box()

    def __hash__(self):
        return 2 + hash(self.start) + hash(self.end)

    def __eq__(self, other_linesegment3d):
        if other_linesegment3d.__class__ != self.__class__:
            return False
        return

    def _bounding_box(self):
        points = [self.start, self.end]

        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])

        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)


    def Length(self):
        return self.end.point_distance(self.start)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.start + curvilinear_abscissa * (
                    self.end - self.start) / self.Length()

    def frenet(self, curvilinear_abscissa):
        return self.unit_direction_vector(), None

    def middle_point(self):
        l = self.Length()
        return self.PointAtCurvilinearAbscissa(0.5 * l)

    def PlaneProjection2D(self, center, x, y):
        return LineSegment2D(self.start.PlaneProjection2D(center, x, y),
                             self.point2.PlaneProjection2D(center, x, y))

    def Intersection(self, segment2):
        x1 = self.start.vector[0]
        y1 = self.start.vector[1]
        z1 = self.start.vector[2]
        x2 = self.point2.vector[0]
        y2 = self.point2.vector[1]
        z2 = self.point2.vector[2]
        x3 = segment2.start.vector[0]
        y3 = segment2.start.vector[1]
        z3 = segment2.start.vector[2]
        x4 = segment2.end_point.vector[0]
        y4 = segment2.end_point.vector[1]
        z4 = segment2.end_point.vector[2]

        if x3 == 0 and x4 == 0 and y4 - y3 == 0:
            x5, y5, z5 = x3, y3, z3
            x6, y6, z6 = x4, y4, z4
            x3, y3, z3 = x1, y1, z1
            x4, y4, z4 = x2, y2, z2
            x1, y1, z1 = x5, y5, z5
            x2, y2, z2 = x6, y6, z6

        elif y3 == 0 and y4 == 0 and x4 - x3 == 0:
            x5, y5, z5 = x3, y3, z3
            x6, y6, z6 = x4, y4, z4
            x3, y3, z3 = x1, y1, z1
            x4, y4, z4 = x2, y2, z2
            x1, y1, z1 = x5, y5, z5
            x2, y2, z2 = x6, y6, z6

        res, list_t1 = [], []

        # 2 unknown 3eq with t1 et t2 unknown
        if (x2 - x1 + y1 - y2) != 0 and (y4 - y3) != 0:
            t1 = (x3 - x1 + (x4 - x3) * (y1 - y3) / (y4 - y3)) / (
                        x2 - x1 + y1 - y2)
            t2 = (y1 - y3 + (y2 - y1) * t1) / (y4 - y3)
            res1 = z1 + (z2 - z1) * t1
            res2 = z3 + (z4 - z3) * t2
            list_t1.append(t1)
            res.append([res1, res2])

        if (z2 - z1 + y1 - y2) != 0 and (y4 - y3) != 0:
            t1 = (z3 - z1 + (z4 - z3) * (y1 - y3) / (y4 - y3)) / (
                        z2 - z1 + y1 - y2)
            t2 = (y1 - y3 + (y2 - y1) * t1) / (y4 - y3)
            res1 = x1 + (x2 - x1) * t1
            res2 = x3 + (x4 - x3) * t2
            list_t1.append(t1)
            res.append([res1, res2])

        if (z2 - z1 + x1 - x2) != 0 and (x4 - x3) != 0:
            t1 = (z3 - z1 + (z4 - z3) * (x1 - x3) / (x4 - x3)) / (
                        z2 - z1 + x1 - x2)
            t2 = (x1 - x3 + (x2 - x1) * t1) / (x4 - x3)
            res1 = y1 + (y2 - y1) * t1
            res2 = y3 + (y4 - y3) * t2
            list_t1.append(t1)
            res.append([res1, res2])

        if len(res) == 0:
            return None

        for pair, t1 in zip(res, list_t1):
            res1, res2 = pair[0], pair[1]
            if math.isclose(res1, res2,
                            abs_tol=1e-7):  # if there is an intersection point
                if t1 >= 0 or t1 <= 1:
                    return Point3D([x1 + (x2 - x1) * t1, y1 + (y2 - y1) * t1,
                                    z1 + (z2 - z1) * t1])

        return None

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            return LineSegment3D(
                *[p.Rotation(center, axis, angle, copy=True) for p in
                  self.points])
        else:
            Edge.Rotation(self, center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def __contains__(self, point):
        point1, point2 = self.start, self.end
        axis = Vector3D(point2 - point1)
        test = point.Rotation(point1, axis, math.pi)
        if test == point:
            return True
        else:
            return False

    def Translation(self, offset, copy=True):
        if copy:
            return LineSegment3D(
                *[p.Translation(offset, copy=True) for p in self.points])
        else:
            Edge.Translation(self, offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return LineSegment3D(
                    *[frame.OldCoordinates(p) for p in self.points])
            else:
                Edge.frame_mapping(self, frame, side, copy=False)
                self.bounding_box = self._bounding_box()
        if side == 'new':
            if copy:
                return LineSegment3D(
                    *[frame.NewCoordinates(p) for p in self.points])
            else:
                Edge.frame_mapping(self, frame, side, copy=False)
                self.bounding_box = self._bounding_box()

    def copy(self):
        return LineSegment3D(*[p.copy() for p in self.points])

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        points = [self.start, self.end]
        x = [p.vector[0] for p in points]
        y = [p.vector[1] for p in points]
        z = [p.vector[2] for p in points]
        ax.plot(x, y, z, 'o-k')
        return ax

    def MPLPlot2D(self, x_3D, y_3D, ax=None, color='k', width=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        edge2D = self.PlaneProjection2D(O3D, x_3D, y_3D)
        edge2D.MPLPlot(ax=ax, color=color, width=width)
        return ax

    def plot_data(self, x_3D, y_3D, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        edge2D = self.PlaneProjection2D(O3D, x_3D, y_3D)
        return edge2D.plot_data(marker, color, stroke_width,
                                dash, opacity, arrow)

    def FreeCADExport(self, name, ndigits=6):
        name = 'primitive' + str(name)
        x1, y1, z1 = round(1000 * self.start, ndigits).vector
        x2, y2, z2 = round(1000 * self.end, ndigits).vector
        return '{} = Part.LineSegment(fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(
            name, x1, y1, z1, x2, y2, z2)

    def to_line(self):
        return Line3D(*self.points)

    def babylon_script(self, color=(1, 1, 1), name='line', type_='line',
                       parent=None):
        if type_ == 'line' or type_ == 'dashed':
            s = 'var myPoints = [];\n'
            s += 'var point1 = new BABYLON.Vector3({},{},{});\n'.format(
                *self.start)
            s += 'myPoints.push(point1);\n'
            s += 'var point2 = new BABYLON.Vector3({},{},{});\n'.format(
                *self.end)
            s += 'myPoints.push(point2);\n'
            if type_ == 'line':
                s += 'var {} = BABYLON.MeshBuilder.CreateLines("lines", {{points: myPoints}}, scene);\n'.format(
                    name)
            elif type_ == 'dashed':
                s += 'var {} = BABYLON.MeshBuilder.CreateDashedLines("lines", {{points: myPoints, dashNb:20}}, scene);'.format(
                    name)
            s += '{}.color = new BABYLON.Color3{};\n'.format(name,
                                                             tuple(color))
        elif type_ == 'tube':
            radius = 0.03 * self.start.point_distance(self.end)
            s = 'var points = [new BABYLON.Vector3({},{},{}), new BABYLON.Vector3({},{},{})];\n'.format(
                *self.start, *self.end)
            s += 'var {} = BABYLON.MeshBuilder.CreateTube("frame_U", {{path: points, radius: {}}}, {});'.format(
                name, radius, parent)
        #            s += 'line.material = red_material;\n'

        else:
            raise NotImplementedError

        if parent is not None:
            s += '{}.parent = {};\n'.format(name, parent)

        return s

    def To2D(self, plane_origin, x1, x2):
        p2D = [p.To2D(plane_origin, x1, x2) for p in (self.start, self.end)]
        return LineSegment2D(*p2D, name=self.name)

    def reverse(self):
        return LineSegment3D(self.end.copy(), self.start.copy())

    def MinimumDistancePoints(self, other_line):
        """
        Returns the points on this line and the other line that are the closest
        of lines
        """
        u = self.end - self.start
        v = other_line.end_point - other_line.start
        w = self.start - other_line.start
        a = u.Dot(u)
        b = u.Dot(v)
        c = v.Dot(v)
        d = u.Dot(w)
        e = v.Dot(w)
        if (a * c - b ** 2) != 0:
            s = (b * e - c * d) / (a * c - b ** 2)
            t = (a * e - b * d) / (a * c - b ** 2)
            p1 = self.start + s * u
            p2 = other_line.start + t * v
            return p1, p2, s, t
        else:
            return None, None, -1, -1

    def Matrix_distance(self, other_line):
        u = self.direction_vector()
        v = other_line.direction_vector()
        w = other_line.start - self.start

        a = u.Dot(u)
        b = -u.Dot(v)
        d = v.Dot(v)

        e = w.Dot(u)
        f = -w.Dot(v)

        A = npy.array([[a, b],
                       [b, d]])
        B = npy.array([e, f])

        res = scp.optimize.lsq_linear(A, B, bounds=(0, 1))
        p1 = self.PointAtCurvilinearAbscissa(res.x[0] * self.Length())
        p2 = other_line.PointAtCurvilinearAbscissa(
            res.x[1] * other_line.Length())
        return p1, p2

    def parallele_distance(self, other_linesegment):
        ptA, ptB, ptC = self.start, self.end, other_linesegment.points[0]
        u = Vector3D((ptA - ptB).vector)
        u.Normalize()
        plane1 = volmdlr.faces.Plane3D.from_3_points(ptA, ptB, ptC)
        v = u.Cross(plane1.normal)  # distance vector
        # ptA = k*u + c*v + ptC
        res = (ptA - ptC).vector
        x, y, z = res[0], res[1], res[2]
        u1, u2, u3 = u.vector[0], u.vector[1], u.vector[2]
        v1, v2, v3 = v.vector[0], v.vector[1], v.vector[2]

        if (u1 * v2 - v1 * u2) != 0 and u1 != 0:
            c = (y * u1 - x * u2) / (u1 * v2 - v1 * u2)
            k = (x - c * v1) / u1
            if math.isclose(k * u3 + c * v3, z, abs_tol=1e-7):
                return k
        elif (u1 * v3 - v1 * u3) != 0 and u1 != 0:
            c = (z * u1 - x * u3) / (u1 * v3 - v1 * u3)
            k = (x - c * v1) / u1
            if math.isclose(k * u2 + c * v2, y, abs_tol=1e-7):
                return k
        elif (v1 * u2 - v2 * u1) != 0 and u2 != 0:
            c = (u2 * x - y * u1) / (v1 * u2 - v2 * u1)
            k = (y - c * v2) / u2
            if math.isclose(k * u3 + c * v3, z, abs_tol=1e-7):
                return k
        elif (v3 * u2 - v2 * u3) != 0 and u2 != 0:
            c = (u2 * z - y * u3) / (v3 * u2 - v2 * u3)
            k = (y - c * v2) / u2
            if math.isclose(k * u1 + c * v1, x, abs_tol=1e-7):
                return k
        elif (u1 * v3 - v1 * u3) != 0 and u3 != 0:
            c = (z * u1 - x * u3) / (u1 * v3 - v1 * u3)
            k = (z - c * v3) / u3
            if math.isclose(k * u2 + c * v2, y, abs_tol=1e-7):
                return k
        elif (u2 * v3 - v2 * u3) != 0 and u3 != 0:
            c = (z * u2 - y * u3) / (u2 * v3 - v2 * u3)
            k = (z - c * v3) / u3
            if math.isclose(k * u1 + c * v1, x, abs_tol=1e-7):
                return k
        else:
            return NotImplementedError

    def minimum_distance(self, element, return_points=False):
        if element.__class__ is Arc3D or element.__class__ is Circle3D:
            pt1, pt2 = element.minimum_distance_points_line(self)
            if return_points:
                return pt1.point_distance(pt2), pt1, pt2
            else:
                return pt1.point_distance(pt2)

        elif element.__class__ is LineSegment3D:
            p1, p2 = self.Matrix_distance(element)
            if return_points:
                return p1.point_distance(p2), p1, p2
            else:
                return p1.point_distance(p2)

        else:
            return NotImplementedError

    def extrusion(self, extrusion_vector):
        u = self.unit_direction_vector()
        v = extrusion_vector.copy()
        v.Normalize()
        w = u.Cross(v)
        l1 = self.Length()
        l2 = extrusion_vector.Norm()
        # outer_contour = Polygon2D([O2D, Point2D((l1, 0.)),
        #                            Point2D((l1, l2)), Point2D((0., l2))])
        plane = volmdlr.surfaces3d.Plane3D(Frame3D(self.start, u, v, w))
        return plane.rectangular_cut(0, l1, 0, l2)

    def revolution(self, axis_point, axis, angle):
        axis_line3d = Line3D(axis_point, axis_point + axis)

        p1_proj, _ = axis_line3d.point_projection(self.start)
        p2_proj, _ = axis_line3d.point_projection(self.end)
        d1 = self.start.point_distance(p1_proj)
        d2 = self.end.point_distance(p2_proj)
        u1 = (self.start - p1_proj)  # Unit vector from p1_proj to p1
        u1.Normalize()
        if u1.is_colinear_to(self.direction_vector()):
            # Planar face
            surface = volmdlr.faces.Plane3D(Frame3D(p1_proj, u1, axis.Cross(u1), axis))
            if angle == two_pi:
                # Only 2 circles as countours
                r, R = sorted([d1, d2])
                v1 = axis.Cross(u1)
                outer_contour3d = Circle3D(Frame3D(p1_proj, u1, v1, axis),
                                           R)
                if not math.isclose(r, 0, abs_tol=1e-9):
                    inner_contours3d = [Circle3D(Frame3D(p1_proj, u1, v1, axis)
                                                 , r)]
                else:
                    inner_contours3d = []
            else:
                # Two arcs and lines
                arc1_i = self.end.Rotation(axis=axis,
                                                 axis_point=axis_point,
                                                 angle=0.5 * angle)
                arc1_e = self.end.Rotation(axis=axis,
                                                 axis_point=axis_point,
                                                 angle=angle)
                arc2_i = self.start.Rotation(axis=axis,
                                                 axis_point=axis_point,
                                                 angle=0.5 * angle)
                arc2_s = self.start.Rotation(axis=axis,
                                                 axis_point=axis_point,
                                                 angle=angle)

                arc1 = Arc3D(self.end, arc1_i, arc1_e)
                arc2 = Arc3D(arc2_s, arc1_i, self.start)
                line2 = LineSegment3D(arc1_e, arc2_s)
                outer_contour3d = Contour3D([self, arc1, line2, arc2])
                inner_contours3d = []
            return volmdlr.faces.PlaneFace3D.from_contours3d(
                outer_contour3d=outer_contour3d,
                inner_contours3d=inner_contours3d)

        elif d1 != d2:
            # Conical
            v = axis.Cross(u1)
            w = axis.Cross(v)
            u1 = self.direction_vector()
            semi_angle = math.asin(u1.Cross(axis).Norm())
            surface = volmdlr.surfaces.ConicalSurface3D(Frame3D(p1_proj, axis, v, w),
                                       semi_angle)
            return surface.rectangular_cut(0, self.Length(), 0, angle)
        else:
            v = axis.Cross(u1)
            surface = volmdlr.surfaces.CylindricalSurface3D(Frame3D(p1_proj, u1, v, axis), d1)
            return surface.rectangular_cut(0, angle,
                                           0, (self.end-self.start).Dot(axis))


class DisplayMesh3D(DisplayMesh):
    _linesegment_class = LineSegment3D
    _point_class = Point3D

    def __init__(self, points: List[Point3D],
                 triangles: List[Tuple[int, int, int]], name=''):
        DisplayMesh.__init__(self, points, triangles)

    def to_babylon(self):
        """
        return mesh in babylon format: https://doc.babylonjs.com/how_to/custom
        """
        positions = []
        for p in self.points:
            positions.extend([k for k in round(p, 6)])

        flatten_indices = []
        for i in self.triangles:
            flatten_indices.extend(i)
        return positions, flatten_indices


class Contour3D(Wire3D):
    _non_serializable_attributes = ['points']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['points', 'name']
    _generic_eq = True
    """
    A collection of 3D primitives forming a closed wire3D
    """

    def __init__(self, primitives, name=''):
        """

        """

        Wire3D.__init__(self, primitives=primitives, name=name)

    def __hash__(self):
        return sum([hash(e) for e in self.primitives]) + sum(
            [hash(p) for p in self.tessel_points])

    def __eq__(self, other_):
        equal = True
        for edge, other_edge in zip(self.primitives, other_.edges):
            equal = (equal and edge == other_edge)
        # for point, other_point in zip(self.points, other_.points):
        #     equal = (equal and point == other_point)
        #     print('contour', equal, point.vector, other_point.vector)
        return equal

    @classmethod
    def from_step(cls, arguments, object_dict):
        edges = []
        for edge in arguments[1]:
            # print(arguments[1])
            edges.append(object_dict[int(edge[1:])])

        return cls(edges, point_inside_contour=None, name=arguments[0][1:-1])

    def clean_points(self):
        """
        TODO : verifier si le dernier point est toujours le meme que le premier point
        lors d'un import step par exemple
        """
        # print('!' , self.primitives)
        if hasattr(self.primitives[0], 'points'):
            points = self.primitives[0].points[:]
        else:
            points = self.primitives[0].tessellation_points()
        for edge in self.primitives[1:]:
            if hasattr(edge, 'points'):
                points_to_add = edge.points[:]
            else:
                points_to_add = edge.tessellation_points()
            if points[0] == points[
                -1]:  # Dans le cas où le (dernier) edge relie deux fois le même point
                points.extend(points_to_add[::-1])

            elif points_to_add[0] == points[-1]:
                points.extend(points_to_add[1:])
            elif points_to_add[-1] == points[-1]:
                points.extend(points_to_add[-2::-1])
            elif points_to_add[0] == points[0]:
                points = points[::-1]
                points.extend(points_to_add[1:])
            elif points_to_add[-1] == points[0]:
                points = points[::-1]
                points.extend(points_to_add[-2::-1])
            else:
                d1, d2 = (points_to_add[0] - points[0]).Norm(), (
                            points_to_add[0] - points[-1]).Norm()
                d3, d4 = (points_to_add[-1] - points[0]).Norm(), (
                            points_to_add[-1] - points[-1]).Norm()
                if math.isclose(d2, 0, abs_tol=1e-3):
                    points.extend(points_to_add[1:])
                elif math.isclose(d4, 0, abs_tol=1e-3):
                    points.extend(points_to_add[-2::-1])
                elif math.isclose(d1, 0, abs_tol=1e-3):
                    points = points[::-1]
                    points.extend(points_to_add[1:])
                elif math.isclose(d3, 0, abs_tol=1e-3):
                    points = points[::-1]
                    points.extend(points_to_add[-2::-1])

        if len(points) > 1:
            if points[0] == points[-1]:
                points.pop()

        return points

    def average_center_point(self):
        nb = len(self.tessel_points)
        x = npy.sum([p[0] for p in self.points]) / nb
        y = npy.sum([p[1] for p in self.points]) / nb
        z = npy.sum([p[2] for p in self.points]) / nb
        return Point3D((x, y, z))

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_edges = [edge.Rotation(center, axis, angle, copy=True) for edge
                         in self.primitives]
            # new_points = [p.Rotation(center, axis, copy=True) for p in self.points]
            return Contour3D(new_edges, None, self.name)
        else:
            for edge in self.primitives:
                edge.Rotation(center, axis, angle, copy=False)
            for point in self.tessel_points:
                point.Rotation(center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            new_edges = [edge.Translation(offset, copy=True) for edge in
                         self.primitives]
            # new_points = [p.Translation(offset, copy=True) for p in self.points]
            return Contour3D(new_edges, None, self.name)
        else:
            for edge in self.primitives:
                edge.Translation(offset, copy=False)
            for point in self.tessel_points:
                point.Translation(offset, copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_edges = [edge.frame_mapping(frame, side, copy=True) for edge in
                         self.primitives]
            # new_points = [p.frame_mapping(frame, side, copy=True) for p in self.points]
            return Contour3D(new_edges, None, self.name)
        else:
            for edge in self.primitives:
                edge.frame_mapping(frame, side, copy=False)
            for point in self.tessel_points:
                point.frame_mapping(frame, side, copy=False)

    def copy(self):
        new_edges = [edge.copy() for edge in self.primitives]
        if self.point_inside_contour is not None:
            new_point_inside_contour = self.point_inside_contour.copy()
        else:
            new_point_inside_contour = None
        return Contour3D(new_edges, new_point_inside_contour, self.name)

    def Length(self):
        # TODO: this is duplicated code from Wire3D!
        length = 0.
        for edge in self.primitives:
            length += edge.Length()
        return length

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        # TODO: this is duplicated code from Wire3D!
        length = 0.
        for primitive in self.primitives:
            primitive_length = primitive.Length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(
                    curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        for edge in self.primitives:
            edge.MPLPlot(ax=ax)

        return ax

    def to_2d(self, plane3d, name=None):
        primitives2d = [plane3d.point3d_to_2d(p) for p in self.primitives]
        return Contour2D(primitives=primitives2d, name=name)

    def _bounding_box(self):
        """
        Flawed method, to be enforced by overloading
        """
        n = 50
        l = self.Length()
        points = [self.PointAtCurvilinearAbscissa(i/n*l)\
                  for i in range(n)]
        return BoundingBox.from_points(points)




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


class Measure2D(LineSegment2D):
    def __init__(self, point1, point2, label='', unit='mm', type_='distance'):
        """
        :param unit: 'mm', 'm' or None. If None, the distance won't be in the label

        """
        # TODO: offset parameter
        LineSegment2D.__init__(self, point1, point2)
        self.label = label
        self.unit = unit
        self.type_ = type_

    def MPLPlot(self, ax, ndigits=6):
        x1, y1 = self.points[0]
        x2, y2 = self.points[1]
        xm, ym = 0.5 * (self.points[0] + self.points[1])
        distance = self.points[1].point_distance(self.points[0])

        if self.label != '':
            label = '{}: '.format(self.label)
        else:
            label = ''
        if self.unit == 'mm':
            label += '{} mm'.format(round(distance * 1000, ndigits))
        else:
            label += '{} m'.format(round(distance, ndigits))

        if self.type_ == 'distance':
            arrow = FancyArrowPatch((x1, y1), (x2, y2),
                                    arrowstyle='<|-|>,head_length=10,head_width=5',
                                    shrinkA=0, shrinkB=0,
                                    color='k')
        elif self.type_ == 'radius':
            arrow = FancyArrowPatch((x1, y1), (x2, y2),
                                    arrowstyle='-|>,head_length=10,head_width=5',
                                    shrinkA=0, shrinkB=0,
                                    color='k')

        ax.add_patch(arrow)
        if x2 - x1 == 0.:
            theta = 90.
        else:
            theta = math.degrees(math.atan((y2 - y1) / (x2 - x1)))
        ax.text(xm, ym, label, va='bottom', ha='center', rotation=theta)





# class Group:
#     def __init__(self, primitives, name):
#         self.primitives = primitives
#         self.name = name





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



