#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing wires & contours.
"""

import itertools
import math
import sys
import warnings
# import random
from collections import deque
from functools import cached_property
from statistics import mean
from typing import List

import matplotlib.pyplot as plt
import networkx as nx
import numpy as npy
import plot_data.core as plot_data
import scipy.integrate as scipy_integrate
from mpl_toolkits.mplot3d import Axes3D
from scipy.spatial import ConvexHull, Delaunay
from triangle import triangulate

import volmdlr
import volmdlr.core
import volmdlr.display as vmd
import volmdlr.edges
import volmdlr.geometry
import volmdlr.utils.common_operations as vm_common_operations
import volmdlr.utils.intersections as vm_utils_intersections
from volmdlr.core_compiled import polygon_point_belongs
from volmdlr.core import EdgeStyle


def argmax(list_of_numbers):
    """
    Returns the max value and the argmax.

    """
    pos_max, max_value = 0, list_of_numbers[0]
    for pos, value in enumerate(list_of_numbers):
        if pos == 0:
            continue
        if value > max_value:
            max_value = value
            pos_max = pos
    return max_value, pos_max


def argmin(list_of_numbers):
    """
    Returns the minimum value from a list of numbers and its index.

    """
    pos_min, min_value = 0, list_of_numbers[0]
    for pos, value in enumerate(list_of_numbers):
        if pos == 0:
            continue
        if value < min_value:
            min_value = value
            pos_min = pos
    return min_value, pos_min


def bounding_rectangle_adjacent_contours(contours: List):
    """
    Compute the bounding box of a list of adjacent contours 2d.

    :param contours: A list of adjacent contours
    :type contours: List[:class:`volmdlr.wires.Contour2D`]
    :return: The bounding box
    :rtype: :class:`volmdlr.core.BoundingRectangle`
    """
    x_min, x_max, y_min, y_max = contours[0].bounding_rectangle.bounds()

    for i in range(1, len(contours)):
        xmin_contour, xmax_contour, ymin_contour, ymax_contour = contours[i].bounding_rectangle.bounds()
        x_min = min(x_min, xmin_contour)
        x_max = max(x_max, xmax_contour)
        y_min = min(y_min, ymin_contour)
        y_max = max(y_max, ymax_contour)

    return volmdlr.core.BoundingRectangle(x_min, x_max, y_min, y_max)


class WireMixin:
    """
    Abstract class for Wire, storing methods and attributes used by many classes in this module.

    """
    _non_data_hash_attributes = ['basis_primitives']
    _non_serializable_attributes = ['primitive_to_index',
                                    'basis_primitives']

    def _data_hash(self):
        return sum(hash(e) for e in self.primitives) + len(self.primitives)

    def length(self):
        length = 0.
        for primitive in self.primitives:
            length += primitive.length()
        return length

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """

        :param angle_resolution: distance between two discretized points.
        """
        length = self.length()
        if number_points:
            n = number_points - 1
        elif angle_resolution:
            n = int(length / angle_resolution) + 1

        return [self.point_at_abscissa(i / n * length) for i in
                range(n + 1)]

    def point_at_abscissa(self, curvilinear_abscissa: float):
        length = 0.
        for primitive in self.primitives:
            primitive_length = primitive.length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.point_at_abscissa(
                    curvilinear_abscissa - length)
            length += primitive_length
        # In case we did not find yet, ask last primitive its end
        if math.isclose(curvilinear_abscissa, length, abs_tol=1e-6):
            return self.primitives[-1].end  # point_at_abscissa(primitive_length)
        raise ValueError('abscissa out of contour length')

    def split_with_two_points(self, point1, point2):
        """
        Split a wire or contour in two points.

        :param point1: spliting point1.
        :param point2: spliting point2.
        :return: List of primitives in between these two points, and another list with the remaining primitives.
        """
        abscissa1 = self.abscissa(point1)
        abscissa2 = self.abscissa(point2)
        if abscissa1 > abscissa2:
            point1, point2 = point2, point1
            abscissa1, abscissa2 = abscissa2, abscissa1
        current_abscissa = 0
        primitives1 = []
        primitives2 = []
        for primitive in self.primitives:
            if abscissa1 < current_abscissa and current_abscissa + primitive.length() < abscissa2:
                primitives1.append(primitive)
            elif current_abscissa > abscissa2 or current_abscissa + primitive.length() < abscissa1:
                primitives2.append(primitive)
            elif current_abscissa <= abscissa1 <= current_abscissa + primitive.length() and \
                    current_abscissa <= abscissa2 <= current_abscissa + primitive.length():
                split_primitives1 = primitive.split(point1)
                if split_primitives1[0]:
                    primitives2.append(split_primitives1[0])
                split_primitives2 = primitive.split(point2)
                if split_primitives2[1]:
                    primitives2.append(split_primitives2[1])
                primitives1.append(primitive.split_between_two_points(point1, point2))
            elif current_abscissa <= abscissa1 <= current_abscissa + primitive.length():
                split_primitives = primitive.split(point1)
                if split_primitives[1]:
                    primitives1.append(split_primitives[1])
                if split_primitives[0]:
                    primitives2.append(split_primitives[0])
            elif current_abscissa <= abscissa2 <= current_abscissa + primitive.length():
                split_primitives = primitive.split(point2)
                if split_primitives[0]:
                    primitives1.append(split_primitives[0])
                if split_primitives[1]:
                    primitives2.append(split_primitives[1])
            else:
                raise NotImplementedError
            current_abscissa += primitive.length()
        return primitives1, primitives2

    def abscissa(self, point, tol=1e-6):
        """
        Compute the curvilinear abscissa of a point on a wire.

        """
        if self.point_over_wire(point, tol):
            length = 0
            for primitive in self.primitives:
                if primitive.point_belongs(point, tol):
                    length += primitive.abscissa(point)
                    break
                length += primitive.length()
            return length

        raise ValueError('Point is not on wire')

    def sort_points_along_wire(self, points):
        """ Sort given points along the wire with respect to the abscissa. """
        return sorted(points, key=self.abscissa)

    def is_ordered(self, tol=1e-6):
        """ Check if the wire's primitives are ordered or not. """

        for primitive_1, primitive_2 in zip(self.primitives, self.primitives[1:]):
            if primitive_1.end.point_distance(primitive_2.start) > tol:
                return False
        return True

    def order_wire(self, tol=1e-6):
        """ Order wire's primitives. """

        if self.is_ordered(tol=tol):
            return self.__class__(self.primitives[:])

        new_primitives = [self.primitives[0]]
        primitives = self.primitives[1:]
        length_primitives = len(primitives) + 1

        while len(new_primitives) < length_primitives:
            for primitive in primitives:
                if new_primitives[0].start.point_distance(primitive.start) < tol:
                    new_primitives.insert(0, primitive.reverse())
                    primitives.remove(primitive)
                elif new_primitives[-1].end.point_distance(primitive.start) < tol:
                    new_primitives.append(primitive)
                    primitives.remove(primitive)
                elif new_primitives[0].start.point_distance(primitive.end) < tol:
                    new_primitives.insert(0, primitive)
                    primitives.remove(primitive)
                elif new_primitives[-1].end.point_distance(primitive.end) < tol:
                    new_primitives.append(primitive.reverse())
                    primitives.remove(primitive)

        return self.__class__(new_primitives)

    @classmethod
    def from_wires(cls, wires):
        """
        Define a wire from successive wires.

        """

        primitives = []
        for wire in wires:
            primitives.extend(wire.primitives)

        wire = cls(primitives)

        if not wire.is_ordered():
            return wire.order_wire()
        return wire

    def inverted_primitives(self):
        """
        Invert wire's primitives.

        """

        new_primitives = []
        for prim in self.primitives[::-1]:
            new_primitives.append(prim.reverse())
        return new_primitives

    def is_followed_by(self, wire_2, tol=1e-6):
        """
        Check if the wire is followed by wire_2.

        """
        return self.primitives[-1].end.point_distance(wire_2.primitives[0].start) < tol

    def point_over_wire(self, point, abs_tol=1e-6):
        """
        Verifies if point is over wire.

        :param point: point to be verified.
        :param abs_tol: tolerance to be considered.
        :return: True or False
        """
        for primitive in self.primitives:
            if primitive.point_belongs(point, abs_tol):
                return True
        return False

    def primitive_over_wire(self, primitive, tol: float = 1e-6):
        """
        Verifies if point is over wire.

        :param primitive: point to be verified.
        :param tol: tolerance to be considered.
        :return: True or False
        """
        points = primitive.discretization_points(number_points=10)
        if all(self.point_over_wire(point, tol) for point in points):
            return True
        return False

    @classmethod
    def from_points(cls, points):
        """
        Create a contour from points with line_segments.

        """
        linesegment_name = 'LineSegment' + points[0].__class__.__name__[-2:]
        edges = []
        for i in range(0, len(points) - 1):
            edges.append(getattr(volmdlr.edges, linesegment_name)(points[i], points[i + 1]))
        contour = cls(edges)
        return contour

    @classmethod
    def from_edge(cls, edge, number_segments: int):
        """
        Creates a Wire object from an edge.

        :param edge: edge used to create Wire.
        :param number_segments: number of segment for the wire to have.
        :return: Wire object.
        """
        points = edge.discretization_points(number_points=number_segments + 1)
        class_name_ = 'Wire' + edge.__class__.__name__[-2:]
        class_ = getattr(sys.modules[__name__], class_name_)
        return class_.from_points(points)


class EdgeCollection3D(WireMixin):
    """
    A collection of simple edges 3D.
    """
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = ['basis_primitives']
    _non_data_eq_attributes = ['name', 'basis_primitives']
    _non_data_hash_attributes = []

    def __init__(self, primitives: List[volmdlr.edges.Edge], color=None, alpha=1, name: str = ''):
        self.primitives = primitives
        self.color = color
        self.alpha = alpha
        self._bbox = None
        self.name = name

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """ Plot edges with matplolib, not tested. """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        for primitive in self.primitives:
            primitive.plot(ax=ax, edge_style=edge_style)
        return ax

    def _bounding_box(self):
        """ Flawed method, to be enforced by overloading. """
        return volmdlr.core.BoundingBox.from_points(self.points())

    @property
    def bounding_box(self):
        """ Get big bounding box of all edges. """
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    def points(self):
        """ Get list of all points. """
        points = []
        for prim in self.primitives:
            points += [prim.start, prim.end]
        return points

    def babylon_param(self):
        """ Get dict for babylonjs object settings. """
        babylon_param = {'alpha': self.alpha,
                         'name': self.name,
                         'color': [0, 0, 0.6]}
        if self.color is None:
            babylon_param['edges_color'] = [0, 0, 0.6]
        else:
            babylon_param['edges_color'] = list(self.color)
        return babylon_param

    def babylon_points(self):
        """ Get list of points coordinates. """
        return [[point.x, point.y, point.z] for point in self.points()]

    def to_babylon(self):
        """ Generate a mesh from all edges for performance when drawing. """
        positions = []
        for prim in self.primitives:
            positions += [prim.start.x, prim.start.y, prim.start.z,
                          prim.end.x, prim.end.y, prim.end.z,
                          prim.end.x, prim.end.y, prim.end.z]

        indices = list(range(len(positions)))
        return positions, indices

    def babylon_meshes(self):
        """ Set the mesh for babylonjs. """
        positions, indices = self.to_babylon()
        babylon_mesh = {'positions': positions,
                        'indices': indices}
        babylon_mesh.update(self.babylon_param())
        return [babylon_mesh]


class Wire2D(volmdlr.core.CompositePrimitive2D, WireMixin):
    """
    A collection of simple primitives, following each other making a wire.

    """

    def __init__(self, primitives: List[volmdlr.core.Primitive2D],
                 name: str = ''):
        self._bounding_rectangle = None
        self._length = None
        volmdlr.core.CompositePrimitive2D.__init__(self, primitives, name)

    def __hash__(self):
        return hash(('wire2d', tuple(self.primitives)))

    def length(self):
        if not self._length:
            self._length = WireMixin.length(self)
        return self._length

    def to_3d(self, plane_origin, x, y):
        """
        Transforms a Wire2D into an Wire3D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Wire3D.
        """
        primitives3d = []
        for edge in self.primitives:
            primitives3d.append(edge.to_3d(plane_origin, x, y))
        return Wire3D(primitives3d)

    def extract_with_points(self, point1: volmdlr.Point2D, point2: volmdlr.Point2D, inside: bool = True):
        """
        Extract primitives between two given points.

        :param point1: extraction point 1.
        :param point2:extraction point2.
        :param inside: If True it'll Extract primitives from smaller point abscissa value
        to greater point abscissa value. If False, it'll return the contour primitives going from
        the greater point abscissa value to the smaller one.
        """
        inside_primitives, outside_primitives = self.split_with_two_points(point1, point2)
        if inside:
            return inside_primitives
        return outside_primitives

        # TODO: method to check if it is a wire

    def infinite_intersections(self, infinite_primitives):
        """
        Returns a list that contains the intersections between a succession of infinite primitives.

        There must be a method implemented to intersect the two infinite primitives.
        """
        offset_intersections = []

        for primitive_1, primitive_2 in zip(infinite_primitives,
                                            infinite_primitives[1:] + [infinite_primitives[0]]):

            i = infinite_primitives.index(primitive_1)
            # k = infinite_primitives.index(primitive_2)

            primitive_name = primitive_1.__class__.__name__.lower().replace('2d', '')
            intersection_method_name = f'{primitive_name}_intersections'
            next_primitive_name = primitive_2.__class__.__name__.lower().replace('2d', '')
            next_intersection_method_name = f'{next_primitive_name}_intersections'

            if hasattr(primitive_1, next_intersection_method_name):
                intersections = getattr(primitive_1, next_intersection_method_name)(
                    primitive_2)
                end = self.primitives[i].end
                if not intersections:
                    continue

                if len(intersections) == 1:
                    offset_intersections.append(intersections[0])

                else:
                    end = self.primitives[i].end
                    if intersections[0].point_distance(end) > intersections[1].point_distance(end):
                        intersections.reverse()
                    offset_intersections.append(intersections[0])

            elif hasattr(primitive_2, intersection_method_name):
                intersections = getattr(primitive_2, intersection_method_name)(primitive_1)
                if not intersections:
                    continue
                if len(intersections) == 1:
                    offset_intersections.append(intersections[0])
                else:
                    end = self.primitives[i].end
                    if intersections[0].point_distance(end) > intersections[
                            1].point_distance(end):
                        intersections.reverse()
                    offset_intersections.append(intersections[0])

            else:
                raise NotImplementedError(
                    f'No intersection method between {primitive_1.__class__.__name__} and'
                    f'{primitive_2.__class__.__name__}. Define {next_intersection_method_name} on '
                    f'{primitive_1.__class__.__name__} or {intersection_method_name} on'
                    f'{primitive_2.__class__.__name__}')

        return offset_intersections

    def offset(self, offset):
        """
        Generates an offset of a Wire2D.

        """
        offset_primitives = []
        infinite_primitives = []
        offset_intersections = []
        # ax = self.plot()
        for primitive in self.primitives:
            infinite_primitive = primitive.infinite_primitive(offset)
            if infinite_primitive is not None:
                infinite_primitives.append(infinite_primitive)
        offset_intersections += self.infinite_intersections(infinite_primitives)
        for i, (point1, point2) in enumerate(zip(offset_intersections,
                                                 offset_intersections[1:] + [offset_intersections[0]])):
            if i + 1 == len(offset_intersections):
                cutted_primitive = infinite_primitives[0].cut_between_two_points(point1, point2)
            else:
                cutted_primitive = infinite_primitives[i + 1].cut_between_two_points(point1, point2)
            offset_primitives.append(cutted_primitive)

        return self.__class__(offset_primitives)

    def plot_data(self, *args, **kwargs):
        data = []
        for item in self.primitives:
            data.append(item.plot_data())
        return data

    def line_intersections(self, line: 'volmdlr.edges.Line2D'):
        """
        Returns a list of intersection of the wire primitives intersecting with the line.

        :returns: a tuple (point, primitive)
        """
        intersection_points = []
        for primitive in self.primitives:
            for point in primitive.line_intersections(line):
                intersection_points.append((point, primitive))
        return intersection_points

    def linesegment_intersections(self,
                                  linesegment: 'volmdlr.edges.LineSegment2D'):
        """
        Returns a list of intersection of the wire primitives intersecting with the line segment.

        :returns: a tuple (point, primitive)
        """
        intersection_points = []
        for primitive in self.primitives:
            inters = primitive.linesegment_intersections(linesegment)
            for point in inters:
                intersection_points.append((point, primitive))
        return intersection_points

    def is_start_end_crossings_valid(self, line, intersections, primitive):
        """
        Returns if the crossings are valid.

        :param line: crossing line
        :param intersections: intersections results
         for primitive line intersections
        :param primitive: intersecting primitive
        :return: None if intersection not a start or
        end point of a contours primitives, or a volmdlr.Point2D if it is.
        """
        primitive_index = self.primitives.index(primitive)
        point1, point2 = None, None
        if intersections[0].is_close(primitive.start):
            point1 = primitive.point_at_abscissa(primitive.length() * 0.01)
            point2 = self.primitives[primitive_index - 1].point_at_abscissa(
                self.primitives[primitive_index - 1].length() * .99
            )

            # point2 = primitive.start + \
            #          self.primitives[primitive_index - 1].unit_direction_vector(0.5)
        elif intersections[0].is_close(primitive.end) and primitive != self.primitives[-1]:
            point1 = primitive.point_at_abscissa(primitive.length() * 0.99)
            point2 = self.primitives[primitive_index + 1].point_at_abscissa(
                self.primitives[primitive_index + 1].length() * .01)

            # point2 = primitive.end + \
            #          self.primitives[primitive_index + 1].unit_direction_vector(0.5)
        if point1 is not None and point2 is not None:
            return line.is_between_points(point1, point2)
        return False

    @staticmethod
    def is_crossing_start_end_point(intersections, primitive):
        """
        Returns True if the crossings provided are start or end of the Wire 2D.

        :param intersections: intersection results
         for primitive line intersections
        :param primitive: intersecting primitive
        :return: False if intersection not a start or
        end point of a contours primitives, or True if it is.
        """
        if intersections[0].is_close(primitive.start) or intersections[0].is_close(primitive.end):
            return True
        return False

    def line_crossings(self, line: volmdlr.edges.Line2D):
        """
        Calculates valid crossing intersections of a wire and an infinite line.

        :param line: line crossing the wire
        :type line: volmdlr.edges.Line2D
        returns a list of Tuples (point, primitive)
        of the wire primitives intersecting with the line
        """
        intersection_points = []
        intersection_points_primitives = []
        for primitive in self.primitives:
            intersections = primitive.line_intersections(line)
            for intersection in intersections:
                if not volmdlr.core.point_in_list(intersection, intersection_points):
                    if not self.is_crossing_start_end_point(intersections, primitive):
                        intersection_points.append(intersection)
                        intersection_points_primitives.append((intersection, primitive))
                    elif self.is_start_end_crossings_valid(line, intersections, primitive):
                        intersection_points.append(intersection)
                        intersection_points_primitives.append((intersection, primitive))
        return intersection_points_primitives

    def wire_intersections(self, wire):
        """
        Compute intersections between two wire 2d.

        :param wire : volmdlr.wires.Wire2D.
        :return: intersections : List[(volmdlr.Point2D, volmdlr.Primitive2D)]
        """
        intersections, intersections_points = [], []
        for primitive in wire.primitives:
            method_name = f'{primitive.__class__.__name__.lower()[0:-2]}_intersections'

            if hasattr(self, method_name):
                a_points = getattr(self, method_name)(primitive)
                # a_points = self.linesegment_intersections(primitive)
                if a_points:
                    for point1, point2 in a_points:
                        if not volmdlr.core.point_in_list(point1, intersections_points):
                            intersections.append([point1, point2])
                            intersections_points.append(point1)
            else:
                raise NotImplementedError(
                    f'Class {self.__class__.__name__} does not implement {method_name}')

        return intersections

    @classmethod
    def from_points(cls, points: List[volmdlr.Point2D]):
        """
        Define a wire based on points 2d with line_segments 2d.

        :param points: points to define wire 2d.
        """
        edges = []
        for i in range(0, len(points) - 1):
            edges.append(volmdlr.edges.LineSegment2D(points[i], points[i + 1]))

        return cls(edges)

    def linesegment_crossings(self,
                              linesegment: 'volmdlr.edges.LineSegment2D'):
        """
        Gets the wire primitives intersecting with the line.

        Returns a list of crossings in the form of a tuple (point, primitive).
        """
        results = self.line_crossings(linesegment.to_line())
        crossings_points = []
        for result in results:
            if linesegment.point_belongs(result[0]):
                crossings_points.append(result)
        return crossings_points

    def wire_crossings(self, wire):
        """
        Compute crossings between two wire 2d.

        :param wire: volmdlr.wires.Wire2D
        :type crossings: List[(volmdlr.Point2D, volmdlr.Primitive2D)]
        """
        crossings, crossings_points = [], []
        for primitive in wire.primitives:
            method_name = f'{primitive.__class__.__name__.lower()[0:-2]}_crossings'

            if hasattr(self, method_name):
                a_points = getattr(self, method_name)(primitive)
                # a_points = self.linesegment_crossings(primitive)
                if a_points:
                    for a in a_points:
                        if not volmdlr.core.point_in_list(a[0], crossings_points):
                            crossings.append([a[0], a[1]])
                            crossings_points.append(a[0])
            else:
                raise NotImplementedError(
                    f'Class {self.__class__.__name__} does not implement {method_name}')

        return crossings

    def to_wire_with_linesegments(self):
        """
        Convert a wire with different primitives to a wire with just line segments.
        """

        wires = []
        for primitive in self.primitives:
            if not isinstance(primitive, volmdlr.edges.LineSegment2D):
                wires.append(primitive.to_wire(10))
            else:
                wires.append(Wire2D([primitive]))

        return Wire2D.from_wires(wires)

    def invert(self):
        return Wire2D(self.inverted_primitives())

    def extend(self, point):
        """
        Extend a wire by adding a line segment connecting the given point to the nearest wire's extremities.
        """

        distances = [self.primitives[0].start.point_distance(point), self.primitives[-1].end.point_distance(point)]
        if distances.index(min(distances)) == 0:
            primitives = [volmdlr.edges.LineSegment2D(point, self.primitives[0].start)]
            primitives.extend(self.primitives)
        else:
            primitives = self.primitives
            primitives.append(volmdlr.edges.LineSegment2D(self.primitives[-1].end, point))

        return Wire2D(primitives)

    def point_distance(self, point):
        """
        Copied from Contour2D.

        """

        min_distance = self.primitives[0].point_distance(point)
        for primitive in self.primitives[1:]:
            distance = primitive.point_distance(point)
            if distance < min_distance:
                min_distance = distance
        return min_distance

    def nearest_primitive_to(self, point):
        """
        Search for the nearest primitive for a point.

        """

        primitives = self.primitives
        primitives_sorted = sorted(primitives, key=lambda primitive: primitive.point_distance(point))

        return primitives_sorted[0]

    def axial_symmetry(self, line):
        """
        Finds out the symmetric wire 2d according to a line.

        """

        primitives_symmetry = []
        for primitive in self.primitives:
            try:
                primitives_symmetry.append(primitive.axial_symmetry(line))
            except NotImplementedError:
                print(f'Class {self.__class__.__name__} does not implement symmetry method')

        return self.__class__(primitives=primitives_symmetry)

    def symmetry(self, line):
        """
        TODO: code this.
        """
        raise NotImplementedError('Not coded yet')

    def is_symmetric(self, wire2d, line):
        """
        Checks if the two wires 2d are symmetric or not according to line.

        """

        c_symmetry_0 = self.symmetry(line)
        c_symmetry_1 = wire2d.symmetry(line)

        if wire2d.is_superposing(c_symmetry_0) and self.is_superposing(c_symmetry_1):
            return True
        return False

    def bsplinecurve_crossings(self,
                               bsplinecurve: 'volmdlr.edges.BSplineCurve2D'):
        """
        Gets the wire primitives crossings with the bsplinecurve.

        Returns a list of crossings in the form of a tuple (point, primitive).
        """

        linesegments = bsplinecurve.to_wire(25).primitives
        crossings_points = []
        for linesegment in linesegments:
            crossings_linesegment = self.linesegment_crossings(linesegment)
            if crossings_linesegment:
                crossings_points.extend(crossings_linesegment)
        return crossings_points

    def bsplinecurve_intersections(self,
                                   bsplinecurve: 'volmdlr.edges.BSplineCurve2D'):
        """
        Gets the wire primitives intersections with the bsplinecurve.

        Returns a list of intersections in the form of a tuple (point, primitive).
        """

        linesegments = bsplinecurve.to_wire(25).primitives
        intersections_points = []
        for linesegment in linesegments:
            intersections_linesegments = self.linesegment_intersections(linesegment)
            if intersections_linesegments:
                intersections_points.extend(intersections_linesegments)
        return intersections_points

    @property
    def bounding_rectangle(self):
        if not self._bounding_rectangle:
            self._bounding_rectangle = self.get_bouding_rectangle()
        return self._bounding_rectangle

    def get_bouding_rectangle(self):
        x_min, x_max, y_min, y_max = self.primitives[0].bounding_rectangle.bounds()
        for edge in self.primitives[1:]:
            xmin_edge, xmax_edge, ymin_edge, ymax_edge = \
                edge.bounding_rectangle.bounds()
            x_min = min(x_min, xmin_edge)
            x_max = max(x_max, xmax_edge)
            y_min = min(y_min, ymin_edge)
            y_max = max(y_max, ymax_edge)
        return volmdlr.core.BoundingRectangle(x_min, x_max, y_min, y_max)

    def middle_point(self):
        return self.point_at_abscissa(self.length() / 2)


class Wire3D(volmdlr.core.CompositePrimitive3D, WireMixin):
    """
    A collection of simple primitives, following each other making a wire.

    """

    def __init__(self, primitives: List[volmdlr.core.Primitive3D],
                 name: str = ''):
        self._bbox = None
        volmdlr.core.CompositePrimitive3D.__init__(self, primitives=primitives, name=name)

    def _bounding_box(self):
        """
        Flawed method, to be enforced by overloading.

        """
        n = 20
        points = []
        for prim in self.primitives:
            points_ = prim.discretization_points(number_points=n)
            for point in points_:
                if point not in points:
                    points.append(point)
        return volmdlr.core.BoundingBox.from_points(points)

    @property
    def bounding_box(self):
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Wire3D.

        :param side: 'old' or 'new'
        """
        new_wire = []
        for primitive in self.primitives:
            new_wire.append(primitive.frame_mapping(frame, side))
        return Wire3D(new_wire)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated in-place.

        :param side: 'old' or 'new'
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for primitive in self.primitives:
            primitive.frame_mapping_inplace(frame, side)

    def minimum_distance(self, wire2):
        distance = []
        for element in self.primitives:
            for element2 in wire2.primitives:
                distance.append(element.minimum_distance(element2))

        return min(distance)

    def point_distance(self, point):
        distance, distance_point = math.inf, None
        for prim in self.primitives:
            prim_distance, prim_point = prim.point_distance(point)
            if prim_distance < distance:
                distance = prim_distance
                distance_point = prim_point
        return distance, distance_point

    def extrusion(self, extrusion_vector):
        faces = []
        for primitive in self.primitives:
            faces.extend(primitive.extrusion(extrusion_vector))
        return faces

    def to_bspline(self, discretization_parameter, degree):
        """
        Convert a wire 3d to a bspline curve 3d.

        """

        discretized_points = self.discretization_points(discretization_parameter)
        bspline_curve = volmdlr.edges.BSplineCurve3D.from_points_interpolation(discretized_points, degree)

        return bspline_curve

    def triangulation(self):
        return None

    def get_primitives_2d(self, plane_origin, x, y):
        """
        Pass primitives to 2d.

        :param plane_origin: plane origin.
        :param x: vector u.
        :param y: vector v.
        :return: list of 2d primitives.
        """
        z = x.cross(y)
        plane3d = volmdlr.faces.Plane3D(volmdlr.Frame3D(plane_origin, x, y, z))
        primitives2d = []
        for primitive in self.primitives:
            primitive2d = plane3d.point3d_to_2d(primitive)
            if primitive2d is not None:
                primitives2d.append(primitive2d)
        return primitives2d

    def to_2d(self, plane_origin, x, y):
        primitives2d = self.get_primitives_2d(plane_origin, x, y)
        return Wire2D(primitives=primitives2d)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Wire3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated Wire3D.
        """
        new_edges = [edge.rotation(center, axis, angle) for edge
                     in self.primitives]
        return Wire3D(new_edges, self.name)


# TODO: define an edge as an opened polygon and allow to compute area from this reference

class ContourMixin(WireMixin):
    """
    Abstract class for Contour, storing methods and attributes used by Contour2D and Contour3D.

    """

    def is_ordered(self, tol=1e-6):
        """
        Verifies if a contour is ordered (primitives following each other).

        :param tol: tolerance to be considered.
        :return: True if ordered, False if not.
        """
        for prim1, prim2 in zip(self.primitives, self.primitives[1:] + [self.primitives[0]]):
            if not prim1.end.is_close(prim2.start, tol):
                return False
        return True

    def ordering_contour(self, tol=1e-6):
        """
        Returns the points of the contour ordered.

        """
        list_point_pairs = [(prim.start, prim.end) for prim in self.primitives]
        length_list_points = len(list_point_pairs)
        points = [list_point_pairs[0]]
        primitives = self.primitives[:]
        new_primitives = [primitives[0]]
        primitives.remove(primitives[0])
        list_point_pairs.remove(
            (list_point_pairs[0][0], list_point_pairs[0][1]))
        finished = False
        counter = 0
        counter1 = 0

        while not finished:
            for i, (point1, point2) in enumerate(list_point_pairs):
                if point1.point_distance(point2) < tol:
                    list_point_pairs.remove((point1, point2))
                    primitives.remove(primitives[i])
                elif point1.point_distance(points[-1][-1]) < tol:
                    points.append((point1, point2))
                    new_primitives.append(primitives[i])
                    primitives.remove(primitives[i])
                    list_point_pairs.remove((point1, point2))
                elif point2.point_distance(points[-1][-1]) < tol:
                    points.append((point2, point1))
                    new_primitives.append(primitives[i].reverse())
                    primitives.remove(primitives[i])
                    list_point_pairs.remove((point1, point2))
                elif point1.point_distance(points[0][0]) < tol:
                    points = [(point2, point1)] + points
                    new_primitives = [primitives[i].reverse()] + new_primitives
                    primitives.remove(primitives[i])
                    list_point_pairs.remove((point1, point2))
                elif point2.point_distance(points[0][0]) < tol:
                    points = [(point1, point2)] + points
                    new_primitives = [primitives[i]] + new_primitives
                    primitives.remove(primitives[i])
                    list_point_pairs.remove((point1, point2))
            if len(list_point_pairs) == 0:
                finished = True
            counter1 += 1
            if counter1 >= 100 * length_list_points:
                self.plot()
                raise NotImplementedError
            if len(list_point_pairs) == 1:
                counter += 1
                if counter > 3:
                    # for point_pair in list_point_pairs:
                    if list_point_pairs[0][0] in points or list_point_pairs[0][::-1] in points:
                        finished = True
                        continue
                    warnings.warn('There may exist a problem with this'
                                  ' contour, it seems it cannot be reordered.'
                                  ' Please, verify its points')
                    finished = True
                    # ax = self.plot()
                    # for point_pair in list_point_pairs:
                    #     point_pair[0].plot(ax=ax, color='r')
                    #     point_pair[1].plot(ax=ax, color='r')
                    # raise NotImplementedError
        return new_primitives

    def order_contour(self):
        """
        Verifies if the contours'primitives are ordered (one after the other). If not, it will order it.

        """
        if self.is_ordered() or len(self.primitives) < 2:
            return self
        new_primitives = self.ordering_contour()
        self.primitives = new_primitives

        return self

    @staticmethod
    def touching_edges_pairs(edges):  # TO DO: move this to edges?
        touching_primitives = []
        for i, primitive1 in enumerate(edges):
            for j, primitive2 in enumerate(edges):
                if j > i:
                    if primitive1.unit_direction_vector(abscissa=0).is_colinear_to(
                            primitive2.unit_direction_vector(abscissa=0)):
                        continue
                    if not primitive2.end.is_close(primitive1.start) and \
                            not primitive1.start.is_close(primitive2.start) and \
                            not primitive2.end.is_close(primitive1.end) and \
                            not primitive1.end.is_close(primitive2.start):
                        if primitive1.point_belongs(primitive2.start) or primitive1.point_belongs(primitive2.end):
                            touching_primitives.append([primitive2, primitive1])
                        elif primitive2.point_belongs(primitive1.start) or primitive2.point_belongs(primitive1.end):
                            touching_primitives.append([primitive1, primitive2])
        return touching_primitives

    @staticmethod
    def contours_primitives_touching_primitives(touching_primitives):
        contours_primitives_lists = []
        for prim1, prim2 in touching_primitives:
            if prim2.point_belongs(prim1.start):
                intersection = prim1.start
            elif prim2.point_belongs(prim1.end):
                intersection = prim1.end
            prim2_split = prim2.split(intersection)
            for prim in prim2_split:
                if not prim:
                    continue
                if prim1.start == prim.start or prim1.end == prim.end:
                    prim = prim.reverse()
                if [prim1, prim] not in contours_primitives_lists:
                    contours_primitives_lists.append([prim1, prim])
        return contours_primitives_lists

    @staticmethod
    def connected_to_splited_primitives(edge, contours_list):
        """
        Verifies if edge is connected to one of the primitives inside contours.

        :param edge: edge for verification.
        :param contours_list: contours lists.
        :return: update contours_primitives_lists and a boolean to indicate if the edge should be removed or not.
        """
        remove = False
        for i, contour in enumerate(contours_list):
            if not contour.primitive_over_contour(edge):
                if volmdlr.core.point_in_list(contour.primitives[0].start, [edge.end, edge.start]):
                    contours_list[i].primitives = [edge.copy(deep=True)] + contour.primitives
                    remove = True
                elif volmdlr.core.point_in_list(contour.primitives[-1].end, [edge.end, edge.start]):
                    contours_list[i].primitives = contour.primitives + [edge.copy(deep=True)]
                    remove = True
        return contours_list, remove

    @staticmethod
    def is_edge_connected(contour_primitives, edge, tol):
        """
        Verifies if edge is connected to one of the primitives inside contour_primitives.

        :param contour_primitives: list of primitives to create a contour.
        :param edge: edge for verification.
        :param tol: tolerance use in verification.
        :return: returns the edge if true, and None if not connected.
        """
        edge_connected = None
        points = [point for prim in contour_primitives for point in prim]
        if (edge.start in points or edge.end in points) and edge not in contour_primitives:
            edge_connected = edge
            return edge_connected

        for point in points:
            if point.is_close(edge.start, tol=tol) and \
                    edge not in contour_primitives:
                edge.start = point
                edge_connected = edge
                return edge_connected
            if point.is_close(edge.end, tol=tol) and \
                    edge not in contour_primitives:
                edge.end = point
                edge_connected = edge
                return edge_connected
        return edge_connected

    @staticmethod
    def find_connected_edges(edges, contours_list, contour_primitives, tol):
        for line in edges:
            if contours_list:
                contours_list, remove = ContourMixin.connected_to_splited_primitives(line, contours_list)
                if remove:
                    edges.remove(line)
                    break
            if not contour_primitives:
                contour_primitives.append(line)
                edges.remove(line)
                break
            edge_connected = ContourMixin.is_edge_connected(contour_primitives, line, tol)
            if edge_connected is not None:
                contour_primitives.append(edge_connected)
                edges.remove(edge_connected)
                break
        return edges, contour_primitives, contours_list

    @staticmethod
    def get_edges_bifurcations(contour_primitives, edges, finished_loop):
        graph = nx.Graph()
        for prim in contour_primitives[:]:
            graph.add_edge(prim.start, prim.end)
        for node in graph.nodes:
            degree = graph.degree(node)
            if degree <= 2:
                continue
            for i, neihgbor in enumerate(graph.neighbors(node)):
                if graph.degree(neihgbor) == 1:
                    i_edge = volmdlr.edges.LineSegment2D(node, neihgbor)
                    if i_edge in contour_primitives:
                        contour_primitives.remove(i_edge)
                        edges.append(volmdlr.edges.LineSegment2D(node, neihgbor))
                        finished_loop = False
                        if i + 1 == degree - 2:
                            break
        return contour_primitives, edges, finished_loop

    @classmethod
    def contours_from_edges(cls, edges, tol=1e-7):
        if not edges:
            return []
        if len(edges) == 1:
            return [cls(edges)]
        # touching_primitives = cls.touching_edges_pairs(edges)
        # for prims in touching_primitives:
        #     if prims[0] in edges:
        #         edges.remove(prims[0])
        #     if prims[1] in edges:
        #         edges.remove(prims[1])
        # contours_primitives_lists = cls.contours_primitives_touching_primitives(touching_primitives)
        # contours_list = [cls(primitives) for primitives in contours_primitives_lists]
        # if not edges:
        #     return contours_list
        contours_primitives_lists = []
        contours_list = []
        list_contours = []
        finished = False
        contour_primitives = []

        while not finished:
            len1 = len(edges)
            edges, contour_primitives, contours_list = cls.find_connected_edges(
                edges, contours_list, contour_primitives, tol)
            if not edges:
                finished = True
            valid = False

            if edges and len(edges) == len1 and len(contour_primitives) != 0:
                valid = True
            elif len(edges) == 0 and len(contour_primitives) != 0:
                valid = True
                finished = True
            if valid:
                contour_primitives, edges, finished = cls.get_edges_bifurcations(contour_primitives,
                                                                                 edges, finished)
                if len(contour_primitives[:]) != 0:
                    contour_n = cls(contour_primitives[:])
                    contour_n.order_contour()
                    list_contours.append(contour_n)
                contour_primitives = []
        list_contours = list_contours + [cls(primitives).order_contour()
                                         for primitives
                                         in contours_primitives_lists]
        valid_contours = [list_contours[0]]
        list_contours.remove(list_contours[0])
        for contour in list_contours:
            for contour2 in valid_contours:
                if contour.is_superposing(contour2):
                    break
            else:
                valid_contours.append(contour)
        return valid_contours

    def discretized_primitives(self, number_points: float):
        """
        Discretize each contour's primitive and return a list of discretized primitives.

        """
        edges = []
        for primitive in self.primitives:
            auto_nb_pts = min(number_points, max(2, int(primitive.length() / 1e-6)))
            points = primitive.discretization_points(number_points=auto_nb_pts)
            for point1, point2 in zip(points[:-1], points[1:]):
                edges.append(volmdlr.edges.LineSegment2D(point1, point2))
        return edges

    def shares_primitives(self, contour):
        """
        Checks if two contour share primitives.

        """
        for prim1 in self.primitives:
            if contour.primitive_over_contour(prim1):
                return True
        return False

    def is_superposing(self, contour2):
        """
        Check if the contours are superposing (one on the other without necessarily having an absolute equality).

        """

        for primitive_2 in contour2.primitives:
            if not self.primitive_over_contour(primitive_2):
                return False
        return True

    def is_overlapping(self, contour2, intersecting_points=None):
        """
        Check if the contours are overlapping (a part of one is on the other).

        """

        if not intersecting_points:
            intersecting_points = self.contour_intersections(contour2)

        if len(intersecting_points) < 2:
            return False

        vec1_2 = volmdlr.edges.LineSegment2D(intersecting_points[0],
                                             intersecting_points[1])
        middle_point = vec1_2.middle_point()
        normal = vec1_2.normal_vector()
        point1 = middle_point + normal * 0.00001
        point2 = middle_point - normal * 0.00001
        if (self.point_belongs(point1) and contour2.point_belongs(point1)) or \
                (not self.point_belongs(point1) and not contour2.point_belongs(point1)) or \
                (self.point_belongs(point1) and self.point_belongs(point2)) or \
                (contour2.point_belongs(point1) and contour2.point_belongs(point2)):
            return True
        return False

    def is_sharing_primitives_with(self, contour):
        """
        Check if two contour are sharing primitives.

        """

        for prim1 in self.primitives:
            for prim2 in contour.primitives:
                shared_section = prim1.get_shared_section(prim2)
                if shared_section:
                    return True
        return False

    def shared_primitives_extremities(self, contour):
        """
        #todo: is this description correct?.

        Extract shared primitives extremities between two adjacent contours.

        """

        if self.is_superposing(contour):
            warnings.warn('The contours are superposing')
            return []

        list_p, edges1 = [], set()
        for edge_1, edge_2 in itertools.product(self.primitives, contour.primitives):
            edges = [edge_1, edge_2, edge_1]
            for edge1, edge2 in zip(edges, edges[1:]):
                for point in [edge2.start, edge2.end]:
                    if edge1.point_belongs(point, 1e-6):
                        if not list_p:
                            list_p.append(point)
                        if list_p and point.point_distance(point.nearest_point(list_p)) > 1e-4:
                            list_p.append(point)
                        try:
                            self.primitive_to_index(edge1)
                            edges1.add(edge1)
                        except KeyError:
                            edges1.add(edge2)

        if len(list_p) < 2:
            warnings.warn('The contours are not adjacent')
            return []

        if len(list_p) == 2:
            return list_p

        contours = self.__class__.contours_from_edges(list(edges1))
        points = []
        for contour_i in contours:
            points_ = contour_i.extremities_points(list_p)
            for point in points_:
                if not volmdlr.core.point_in_list(point, points):
                    points.append(point)

        return points

    def shared_primitives_with(self, contour):
        """
        Extract shared primitives between two adjacent contours.

        """
        shared_primitives_1 = []
        shared_primitives_2 = []

        for prim1 in self.primitives:
            for prim2 in contour.primitives:
                shared_section_1 = prim1.get_shared_section(prim2)
                shared_section_2 = prim2.get_shared_section(prim1)
                if shared_section_1:
                    shared_primitives_1.extend(shared_section_1)
                if shared_section_2:
                    shared_primitives_2.extend(shared_section_2)
        return shared_primitives_1, shared_primitives_2

    def delete_shared_contour_section(self, contour):
        """
        Delete shared primitives between two adjacent contours.

        :param contour: other contour.
        :return: list of new primitives, without those shared by both contours.
        """
        new_primitives_contour1 = self.primitives[:]
        new_primitives_contour2 = contour.primitives[:]
        for prim1 in self.primitives:
            for prim2 in contour.primitives:
                shared_section = prim1.get_shared_section(prim2)
                if shared_section:
                    prim1_delete_shared_section = prim1.delete_shared_section(shared_section[0])
                    prim2_delete_shared_section = prim2.delete_shared_section(shared_section[0])
                    if prim1 in new_primitives_contour1:
                        new_primitives_contour1.remove(prim1)
                    if prim2 in new_primitives_contour2:
                        new_primitives_contour2.remove(prim2)
                    new_primitives_contour1.extend(prim1_delete_shared_section)
                    new_primitives_contour2.extend(prim2_delete_shared_section)

        return new_primitives_contour1 + new_primitives_contour2

    def merge_primitives_with(self, contour):
        """
        Extract not shared primitives between two adjacent contours, to be merged.

        :param contour:
        :return:
        """
        merge_primitives = self.delete_shared_contour_section(contour)
        return merge_primitives

    def edges_order_with_adjacent_contour(self, contour):
        """
        Check if the shared edges between two adjacent contours are traversed with two \
        different directions along each contour.

        """

        contour1 = self
        contour2 = contour

        # shared_tuple = contour1.shared_edges_between2contours(contour2)
        shared_tuple = contour1.shared_primitives_with(contour2)
        # [shared_primitives_1, shared_primitives_2] = contour1.shared_primitives_with(contour2)

        # p1_start = contour1.primitives[shared_tuple[0][0]].start
        # p2_start = contour2.primitives[shared_tuple[0][1]].start
        # p2_end = contour2.primitives[shared_tuple[0][1]].end

        p1_start = shared_tuple[0][0].start
        p2_start = shared_tuple[1][-1].start
        p2_end = shared_tuple[1][-1].end

        if (p1_start.point_distance(p2_start)) < \
                (p1_start.point_distance(p2_end)):
            return False
        return True

    def extremities_points(self, list_p):
        """
        Return extremities points of a list of points on a contour.

        """
        # TODO: rewrite this awfull code!
        points = []
        primitives = self.primitives
        for prim in primitives:
            pts = []
            for point in list_p:  # due to errors
                if prim.point_belongs(point):
                    pts.append(point)
            if len(pts) == 1:
                points.append(pts[0])
                break
            if len(pts) > 1:
                points.append(prim.start.nearest_point(pts))
                break

        for i in range(len(primitives) - 1, -1, -1):
            pts = []
            for point in list_p:  # due to errors
                if primitives[i].point_belongs(point):
                    pts.append(point)
            if len(pts) == 1:
                if not volmdlr.core.point_in_list(pts[0], points):
                    points.append(pts[0])
                    break
            elif len(pts) > 1:
                point = primitives[i].end.nearest_point(pts)
                if not volmdlr.core.point_in_list(point, points):
                    points.append(point)
                    break
        return points

    def primitive_over_contour(self, primitive, tol: float = 1e-6):
        return self.primitive_over_wire(primitive, tol)

    def point_over_contour(self, point, abs_tol=1e-6):
        return self.point_over_wire(point, abs_tol)

    def get_geo_lines(self, tag: int, primitives_tags: List[int]):
        """
        Gets the lines that define a Contour in a .geo file.

        :param tag: The contour index
        :type tag: int
        :param primitives_tags: The contour's primitives index
        :type primitives_tags: List[int]

        :return: A line
        :rtype: str
        """

        return 'Line Loop(' + str(tag) + ') = {' + str(primitives_tags)[1:-1] + '};'

    def get_geo_points(self):
        points = set()
        for primitive in self.primitives:
            points.update(primitive.get_geo_points())
        return points

    def to_polygon(self, angle_resolution, discretize_line: bool = False):
        """
        Transform the contour_mixin to a polygon, COPY/PASTE from Contour2D.

        :param angle_resolution: Number of points per radians.
        :type angle_resolution: float
        :param discretize_line: Boolean indicating whether the line segments should be discretized or not.
        :type discretize_line: bool
        :return: The discretized version of the contour.
        :rtype: ClosedPolygon2D
        """

        polygon_points = []

        for primitive in self.primitives:
            if isinstance(primitive, volmdlr.edges.LineSegment2D) and not discretize_line:
                polygon_points.append(primitive.start)
            else:
                polygon_points.extend(primitive.discretization_points(angle_resolution=angle_resolution)[:-1])

        if isinstance(self, Contour2D):
            return ClosedPolygon2D(polygon_points)
        return ClosedPolygon3D(polygon_points)

    @classmethod
    def extract_contours(cls, contour, point1, point2, inside=False):

        new_primitives = contour.extract_with_points(point1, point2, inside)
        contours = [cls(new_primitives)]
        return contours


class Contour2D(ContourMixin, Wire2D):
    """
    A collection of 2D primitives forming a closed wire2D.

    TODO : center_of_mass and second_moment_area should be changed accordingly
    to area considering the triangle drawn by the arcs
    """
    _non_data_hash_attributes = ['_internal_arcs', '_external_arcs',
                                 '_polygon', '_straight_line_contour_polygon',
                                 'primitive_to_index',
                                 'basis_primitives', '_utd_analysis']
    _non_serializable_attributes = ['_internal_arcs', '_external_arcs',
                                    '_polygon',
                                    '_straight_line_contour_polygon',
                                    'primitive_to_index',
                                    'basis_primitives', '_utd_analysis']

    def __init__(self, primitives: List[volmdlr.edges.Edge],
                 name: str = ''):
        Wire2D.__init__(self, primitives, name)
        self._edge_polygon = None
        self._polygon_100_points = None
        self._area = None

    def __hash__(self):
        return hash(tuple(self.primitives))

    def __eq__(self, other_):
        if id(self) == id(other_):
            return True
        if other_.__class__.__name__ != self.__class__.__name__:
            return False
        if len(self.primitives) != len(other_.primitives) or self.length() != other_.length():
            return False
        equal = 0
        for prim1 in self.primitives:
            reverse1 = prim1.reverse()
            found = False
            for prim2 in other_.primitives:
                reverse2 = prim2.reverse()
                if (prim1 == prim2 or reverse1 == prim2
                        or reverse2 == prim1 or reverse1 == reverse2):
                    equal += 1
                    found = True
            if not found:
                return False
        if equal == len(self.primitives):
            return True
        return False

    @property
    def edge_polygon(self):
        if self._edge_polygon is None:
            self._edge_polygon = self._get_edge_polygon()
        return self._edge_polygon

    def _get_edge_polygon(self):
        points = []
        for edge in self.primitives:
            if points:
                if not edge.start.is_close(points[-1]):
                    points.append(edge.start)
            else:
                points.append(edge.start)
        closedpolygon = ClosedPolygon2D(points)
        return closedpolygon

    def to_3d(self, plane_origin, x, y):
        """
        Transforms a Contour2D into an Contour3D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Contour3D.
        """
        p3d = []
        for edge in self.primitives:
            p3d.append(edge.to_3d(plane_origin, x, y))

        return Contour3D(p3d)

    def point_belongs(self, point, include_edge_points: bool = False):
        # TODO: This is incomplete!!!
        x_min, x_max, y_min, y_max = self.bounding_rectangle
        if point.x < x_min or point.x > x_max or point.y < y_min or point.y > y_max:
            return False
        # if self.edge_polygon.point_belongs(point):
        #     return True
        # for edge in self.primitives:
        #     if hasattr(edge, 'straight_line_point_belongs'):
        #         if edge.straight_line_point_belongs(point):
        #             return True
        #     warnings.warn(f'{edge.__class__.__name__} does not implement straight_line_point_belongs yet')
        if not self._polygon_100_points:
            self._polygon_100_points = self.to_polygon(100)
        if self._polygon_100_points.point_belongs(point):
            return True
        return False

    def bounding_points(self):
        """Bounding points (x_min, y_min) (x_max, y_max)."""
        points = self.edge_polygon.points[:]
        for primitive in self.primitives:
            if hasattr(primitive, 'discretization_points'):
                points.extend(primitive.discretization_points(number_points=10))
        x_min = min(point[0] for point in points)
        x_max = max(point[0] for point in points)
        y_min = min(point[1] for point in points)
        y_max = max(point[1] for point in points)
        return volmdlr.Point2D(x_min, y_min), volmdlr.Point2D(x_max, y_max)

    def area(self):
        if not self._area:
            area = self.edge_polygon.area()
            classes = {prim.__class__ for prim in self.primitives}
            verify_classes = classes.issubset({volmdlr.edges.LineSegment2D, volmdlr.edges.Arc2D})
            if verify_classes:
                if self.edge_polygon.is_trigo:
                    trigo = 1
                else:
                    trigo = -1
                for edge in self.primitives:
                    area += trigo * edge.straight_line_area()
                    self._area = abs(area)
            else:
                polygon = self.to_polygon(angle_resolution=50)
                self._area = polygon.triangulation().area()
        return self._area

    def center_of_mass(self):
        """
        Calculates the center of mass of the Contour2D.

        :return: Contour's center of mass.
        """
        center = self.edge_polygon.area() * self.edge_polygon.center_of_mass()
        # ax = self.plot()
        # self.edge_polygon.center_of_mass().plot(ax=ax, color='b')
        if self.edge_polygon.is_trigo:
            trigo = 1
        else:
            trigo = -1
        for edge in self.primitives:
            # edge.straight_line_center_of_mass().plot(ax=ax, color='g')
            center += trigo * edge.straight_line_area() \
                      * edge.straight_line_center_of_mass()

        return center / self.area()

    def second_moment_area(self, point):

        second_moment_area_x, second_moment_area_y, second_moment_area_xy = self.edge_polygon.second_moment_area(point)
        for edge in self.primitives:
            second_moment_area_x_e, second_moment_area_y_e, second_moment_area_xy_e =\
                edge.straight_line_second_moment_area(point)
            if self.edge_polygon.is_trigo:
                second_moment_area_x += second_moment_area_x_e
                second_moment_area_y += second_moment_area_y_e
                second_moment_area_xy += second_moment_area_xy_e
            else:
                second_moment_area_x -= second_moment_area_x_e
                second_moment_area_y -= second_moment_area_y_e
                second_moment_area_xy -= second_moment_area_xy_e

        return second_moment_area_x, second_moment_area_y, second_moment_area_xy

    def plot_data(self, edge_style: plot_data.EdgeStyle = None,
                  surface_style: plot_data.SurfaceStyle = None):
        plot_data_primitives = [item.plot_data() for item in self.primitives]
        return plot_data.Contour2D(plot_data_primitives=plot_data_primitives,
                                   edge_style=edge_style,
                                   surface_style=surface_style,
                                   name=self.name)

    def is_inside(self, contour2):
        """
        Verifies if a contour is inside another contour perimiter, including the edges.

        :returns: True or False
        """
        if contour2.area() > self.area():
            return False
        points_contour2 = []
        for prim in contour2.primitives:
            if not volmdlr.core.point_in_list(prim.start, points_contour2):
                points_contour2.append(prim.start)
            if not volmdlr.core.point_in_list(prim.end, points_contour2):
                points_contour2.append(prim.end)
            points_contour2.extend(prim.discretization_points(number_points=10))
        for point in points_contour2:
            if not self.point_belongs(point) and not self.point_over_contour(point, abs_tol=1e-7):
                return False
        return True

    def inverted_primitives(self):
        new_primitives = []
        for prim in self.primitives[::-1]:
            new_primitives.append(prim.reverse())
        return new_primitives

    def invert(self):
        return Contour2D(self.inverted_primitives())

    def random_point_inside(self, include_edge_points: bool = False):
        """
        Finds a random point inside the polygon.

        :param include_edge_points: Choose True if you want to consider a point on the polygon inside.
        :type include_edge_points: bool
        :return: A random point inside the polygon
        :rtype: `volmdlr.Point2D`
        """
        x_min, x_max, y_min, y_max = self.bounding_rectangle.bounds()
        for _ in range(2000):
            point = volmdlr.Point2D.random(x_min, x_max, y_min, y_max)
            if self.point_belongs(point, include_edge_points):
                return point
        raise ValueError('Could not find a point inside')

    def cut_by_linesegments(self, lines: List[volmdlr.edges.LineSegment2D]):
        cut_lines = []
        for cut_ls in lines:
            cut_lines.append(cut_ls.to_line())

        contour_to_cut = [self]
        for line in cut_lines:
            new_contour_to_cut = []
            for contour in contour_to_cut:
                cutted_contour = contour.cut_by_line(line)
                new_contour_to_cut.extend(cutted_contour)
            contour_to_cut.extend(new_contour_to_cut)

        point1 = Contour2D(lines).center_of_mass()
        dist_min = math.inf
        c_opti = None
        for contour in contour_to_cut:
            if contour.area() > 1e-10:
                point0 = contour.center_of_mass()
                if point0.point_distance(point1) < dist_min:
                    c_opti = contour
                    dist_min = point0.point_distance(point1)
        return c_opti

    def repair_cut_contour(self, n, intersections, line):
        """
        Repair contour.

        Choose:
        n=0 for Side 1: opposite side of beginning of contour
        n=1 for Side 2: start of contour to first intersect (i=0) and
         i odd to i+1 even
        """
        if n not in [0, 1]:
            raise ValueError

        n_inter = len(intersections)
        contours = []
        # primitives_split = [primitive.split(point)
        #                     for point, primitive in intersections]
        x = [(ip, line.abscissa(point))
             for ip, (point, _) in enumerate(intersections)]
        # intersection_to_primitives_index = {
        #     i: self.primitives.index(primitive)
        #     for i, (_, primitive) in enumerate(intersections)}
        sorted_inter_index = [x[0] for x in sorted(x, key=lambda p: p[1])]
        sorted_inter_index_dict = {i: ii for ii, i in
                                   enumerate(sorted_inter_index)}
        sorted_inter_index_dict[n_inter] = sorted_inter_index_dict[0]
        if n == 1:
            intersections.append(intersections[0])

        remaining_transitions = list(range(n_inter // 2))
        # enclosing_transitions = {}
        while len(remaining_transitions) > 0:
            nb_max_enclosed_transitions = -1
            enclosed_transitions = {}
            for it in remaining_transitions:
                i1 = sorted_inter_index_dict[2 * it + n]
                i2 = sorted_inter_index_dict[2 * it + 1 + n]
                net = abs(i2 - i1) - 1
                if net > nb_max_enclosed_transitions:
                    nb_max_enclosed_transitions = net
                    best_transition = it
                    if i1 < i2:
                        enclosed_transitions[it] = [(i + abs(n - 1)) // 2 for i
                                                    in sorted_inter_index[
                                                       i2 - 1:i1:-2]]
                    else:
                        enclosed_transitions[it] = [(i + abs(n - 1)) // 2 for i
                                                    in sorted_inter_index[
                                                       i2 + 1:i1:2]]

            remaining_transitions.remove(best_transition)
            point_start, primitive1 = intersections[2 * best_transition + n]
            point2, primitive2 = intersections[2 * best_transition + 1 + n]
            primitives = self.extract_primitives(point_start, primitive1,
                                                 point2, primitive2,
                                                 inside=not n)
            last_point = point2
            for transition in enclosed_transitions[best_transition]:
                point1, primitive1 = intersections[2 * transition + n]
                point2, primitive2 = intersections[2 * transition + 1 + n]
                primitives.append(
                    volmdlr.edges.LineSegment2D(last_point, point1))
                primitives.extend(
                    self.extract_primitives(point1, primitive1, point2,
                                            primitive2, inside=not n))
                last_point = point2
                if transition in remaining_transitions:
                    remaining_transitions.remove(transition)

            primitives.append(
                volmdlr.edges.LineSegment2D(last_point, point_start))

            contour = Contour2D(primitives)
            contour.order_contour()
            contours.append(contour)
        return contours

    def cut_by_line(self, line: volmdlr.edges.Line2D) -> List['Contour2D']:
        """
        :param line: The line used to cut the contour.

        :return: A list of resulting contours
        """
        intersections = self.line_crossings(line)
        if not intersections or len(intersections) < 2:
            return [self]
        points_intersections = [point for point, prim in intersections]
        sorted_points = line.sort_points_along_line(points_intersections)
        list_contours = []
        contour_to_cut = self

        for point1, point2 in zip(sorted_points[:-1], sorted_points[1:]):
            closing_line = volmdlr.edges.LineSegment2D(point1, point2)
            if not contour_to_cut.point_belongs(closing_line.middle_point()):
                continue
            closing_contour = Contour2D([closing_line])
            contour1, contour2 = contour_to_cut.get_divided_contours(point1, point2, closing_contour, True)
            if sorted_points.index(point1) + 2 <= len(sorted_points) - 1:
                if contour1.point_over_contour(sorted_points[sorted_points.index(point1) + 2]):
                    contour_to_cut = contour1
                    list_contours.append(contour2)
                elif contour2.point_over_contour(sorted_points[sorted_points.index(point1) + 2]):
                    contour_to_cut = contour2
                    list_contours.append(contour1)
            else:
                list_contours.extend([contour1, contour2])

        return list_contours

    def split_by_line(self, line: volmdlr.edges.Line2D) -> List['Contour2D']:
        intersections = self.line_crossings(line)
        intersections = [point for point, prim in intersections]
        if not intersections:
            return [self]
        if len(intersections) < 2:
            extracted_outerpoints_contour1 = \
                volmdlr.wires.Contour2D.extract_contours(self, self.primitives[0].start, intersections[0], True)[0]
            extracted_innerpoints_contour1 = \
                volmdlr.wires.Contour2D.extract_contours(self, intersections[0], self.primitives[-1].end, True)[0]
            return extracted_outerpoints_contour1, extracted_innerpoints_contour1
        if len(intersections) == 2:
            extracted_outerpoints_contour1 = \
                volmdlr.wires.Contour2D.extract_contours(self, intersections[0], intersections[1], True)[0]
            extracted_innerpoints_contour1 = \
                volmdlr.wires.Contour2D.extract_contours(self, intersections[0], intersections[1], False)[0]
            return extracted_innerpoints_contour1, extracted_outerpoints_contour1
        raise NotImplementedError

    def split_regularly(self, n):
        """
        Split in n slices.

        """
        x_min, x_max, _, _ = self.bounding_rectangle.bounds()
        cutted_contours = []
        iteration_contours = [self]
        for i in range(n - 1):
            xi = x_min + (i + 1) * (x_max - x_min) / n
            cut_line = volmdlr.edges.Line2D(volmdlr.Point2D(xi, 0),
                                            volmdlr.Point2D(xi, 1))

            iteration_contours2 = []
            for c in iteration_contours:
                sc = c.cut_by_line(cut_line)
                lsc = len(sc)
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

    def grid_triangulation(self, x_density: float = None,
                           y_density: float = None,
                           min_points_x: int = 20,
                           min_points_y: int = 20,
                           number_points_x: int = None,
                           number_points_y: int = None):
        """
        Compute a triangulation using an n-by-m grid to triangulate the contour.
        """
        bounding_rectangle = self.bounding_rectangle
        # xmin, xmax, ymin, ymax = self.bounding_rectangle
        dx = bounding_rectangle[1] - bounding_rectangle[0]  # xmax - xmin
        dy = bounding_rectangle[3] - bounding_rectangle[2]  # ymax - ymin
        if number_points_x is None:
            n = max(math.ceil(x_density * dx), min_points_x)
        if number_points_y is None:
            m = max(math.ceil(y_density * dy), min_points_y)
        x = [bounding_rectangle[0] + i * dx / n for i in range(n + 1)]
        y = [bounding_rectangle[2] + i * dy / m for i in range(m + 1)]

        point_index = {}
        number_points = 0
        points = []
        triangles = []
        for xi in x:
            for yi in y:
                point = volmdlr.Point2D(xi, yi)
                if self.point_belongs(point):
                    point_index[point] = number_points
                    points.append(point)
                    number_points += 1

        for i in range(n):
            for j in range(m):
                point1 = volmdlr.Point2D(x[i], y[j])
                point2 = volmdlr.Point2D(x[i + 1], y[j])
                point3 = volmdlr.Point2D(x[i + 1], y[j + 1])
                point4 = volmdlr.Point2D(x[i], y[j + 1])
                points_in = []
                for point in [point1, point2, point3, point4]:
                    if point in point_index:
                        points_in.append(point)
                if len(points_in) == 4:
                    triangles.append(
                        [point_index[point1], point_index[point2], point_index[point3]])
                    triangles.append(
                        [point_index[point1], point_index[point3], point_index[point4]])

                elif len(points_in) == 3:
                    triangles.append([point_index[point] for point in points_in])

        return vmd.DisplayMesh2D(points, triangles)

    def contour_intersections(self, contour2d):
        intersecting_points = []
        for primitive1 in self.primitives:
            for primitive2 in contour2d.primitives:
                line_intersection = primitive1.linesegment_intersections(primitive2)
                if line_intersection:
                    if not volmdlr.core.point_in_list(line_intersection[0], intersecting_points):
                        intersecting_points.extend(line_intersection)
                else:
                    touching_points = primitive1.touching_points(primitive2)
                    for point in touching_points:
                        if not volmdlr.core.point_in_list(point, intersecting_points):
                            intersecting_points.append(point)
            if len(intersecting_points) == 2:
                break
        return intersecting_points

    def get_divided_contours(self, cutting_point1: volmdlr.Point2D,
                             cutting_point2: volmdlr.Point2D,
                             closing_contour,
                             inside: bool):
        extracted_outerpoints_contour1 = \
            volmdlr.wires.Contour2D.extract_contours(self,
                                                     cutting_point1,
                                                     cutting_point2,
                                                     inside)[0]
        extracted_innerpoints_contour1 = \
            volmdlr.wires.Contour2D.extract_contours(self,
                                                     cutting_point1,
                                                     cutting_point2,
                                                     not inside)[0]
        primitives1 = extracted_outerpoints_contour1.primitives + closing_contour.primitives
        primitives2 = extracted_innerpoints_contour1.primitives + closing_contour.primitives
        if extracted_outerpoints_contour1.primitives[0].start.is_close(closing_contour.primitives[0].start):
            cutting_contour_new = closing_contour.invert()
            primitives1 = cutting_contour_new.primitives + \
                extracted_outerpoints_contour1.primitives
        elif extracted_outerpoints_contour1.primitives[0].start.is_close(closing_contour.primitives[-1].end):
            primitives1 = closing_contour.primitives + \
                          extracted_outerpoints_contour1.primitives

        if extracted_innerpoints_contour1.primitives[0].start.is_close(closing_contour.primitives[0].start):
            cutting_contour_new = \
                closing_contour.invert()
            primitives2 = cutting_contour_new.primitives + \
                extracted_innerpoints_contour1.primitives
        elif extracted_innerpoints_contour1.primitives[0].start.is_close(closing_contour.primitives[-1].end):
            primitives2 = closing_contour.primitives + \
                          extracted_innerpoints_contour1.primitives
        contour1 = Contour2D(primitives1)
        contour1.order_contour()
        contour2 = Contour2D(primitives2)
        contour2.order_contour()
        return contour1, contour2

    def divide(self, contours, inside):
        new_base_contours = [self]
        finished = False
        counter = 0
        list_contour = contours[:]
        list_cutting_contours = contours[:]
        list_valid_contours = []
        while not finished:
            if not contours:
                break
            cutting_contour = contours[0]
            for base_contour in new_base_contours:
                cutting_points = []
                point1, point2 = [cutting_contour.primitives[0].start,
                                  cutting_contour.primitives[-1].end]
                if not any(base_contour.point_belongs(prim.middle_point()) for prim in cutting_contour.primitives):
                    continue
                if base_contour.point_over_contour(point1) and base_contour.point_over_contour(point2):
                    cutting_points = [point1, point2]
                elif len(new_base_contours) == 1:
                    contours.remove(cutting_contour)
                    continue
                if cutting_points:
                    contour1, contour2 = base_contour.get_divided_contours(
                        cutting_points[0], cutting_points[1], cutting_contour, inside)

                    new_base_contours_ = []
                    for cntr in [contour1, contour2]:
                        all_divided_contour = True
                        for cut_contour in list_cutting_contours:
                            points_at_abs = [prim.middle_point() for prim in cut_contour.primitives]
                            for point_at_abs in points_at_abs:
                                if cntr.point_belongs(point_at_abs) and \
                                        (not cntr.point_over_contour(point_at_abs) and
                                         True not in [cntr.primitive_over_contour(prim)
                                                      for prim in cut_contour.primitives]):
                                    all_divided_contour = False
                                    break
                            else:
                                continue
                            break
                        if all_divided_contour and not math.isclose(cntr.area(), 0.0, abs_tol=1e-6):
                            list_valid_contours.append(cntr)
                        else:
                            new_base_contours_.append(cntr)
                    contours.remove(cutting_contour)
                    break
            else:
                continue
            new_base_contours.remove(base_contour)
            new_base_contours.extend(new_base_contours_)
            if len(contours) == 1 and not new_base_contours:
                finished = True
                continue
            counter += 1
            if counter >= 100 * len(list_contour):
                # if base_contour.is_inside(contours[0]):
                #     contours.remove(cutting_contour)
                #     continue
                # list_valid_contours.append(base_contour)
                # finished = True
                contours = contours[::-1]
                if counter > 100 * len(list_contour) + len(contours):
                    # print('new_base_contours:', len(new_base_contours))
                    # print('len(contours):', len(contours))
                    # ax = contours[0].plot()
                    # base_contour.plot(ax=ax, color='b')
                    warnings.warn('There probably exists an open contour (two wires that could not be connected)')
                    finished = True

        return list_valid_contours

    def discretized_contour(self, n: float):
        """
        Discretize each contour's primitive and return a new contour with teses discretized primitives.
        """
        contour = Contour2D((self.discretized_primitives(n)))

        return contour.order_contour()

    @classmethod
    def from_bounding_rectangle(cls, x_min, x_max, y_min, y_max):
        """
        Create a contour 2d with bounding_box parameters, using line segments 2d.

        """

        edge0 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(x_min, y_min), volmdlr.Point2D(x_max, y_min))
        edge1 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(x_max, y_min), volmdlr.Point2D(x_max, y_max))
        edge2 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(x_max, y_max), volmdlr.Point2D(x_min, y_max))
        edge3 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(x_min, y_max), volmdlr.Point2D(x_min, y_min))

        edges = [edge0, edge1, edge2, edge3]

        return Contour2D(edges)

    def cut_by_bspline_curve(self, bspline_curve2d: volmdlr.edges.BSplineCurve2D):
        """
        Cut a contour 2d with bspline_curve 2d to define two different contours.

        """
        # TODO: BsplineCurve is descretized and defined with a wire. To be improved!

        contours = self.cut_by_wire(Wire2D.from_edge(bspline_curve2d, 20))

        return contours

    def clean_primitives(self):
        """
        Delete primitives with start=end, and return a new contour.
        """

        new_primitives = []
        for prim in self.primitives:
            if prim.start != prim.end:
                new_primitives.append(prim)

        return Contour2D(new_primitives)

    def merge_with(self, contour2d):
        """
        Merge two adjacent contours, and returns one outer contour and inner contours (if there are any).

        :param contour2d: contour to merge with.
        :return: merged contours.
        """
        is_sharing_primitive = self.is_sharing_primitives_with(contour2d)
        if self.is_inside(contour2d) and not is_sharing_primitive:
            return [self]
        if contour2d.is_inside(self) and not is_sharing_primitive:
            return [contour2d]

        merged_primitives = self.delete_shared_contour_section(contour2d)
        contours = Contour2D.contours_from_edges(merged_primitives)
        contours = sorted(contours, key=lambda contour: contour.area(),
                          reverse=True)
        return contours

    def union(self, contour2: 'Contour2D'):
        """
        Union two contours, if they adjacent, or overlap somehow.

        """
        if self.is_inside(contour2):
            return [self]
        if contour2.is_inside(self):
            return [contour2]
        contours_intersections = self.contour_intersections(contour2)
        if not self.is_sharing_primitives_with(contour2) and contours_intersections:
            resulting_primitives = []
            primitives1_inside = self.extract_with_points(contours_intersections[0], contours_intersections[1], True)
            primitives1_outside = self.extract_with_points(contours_intersections[0], contours_intersections[1], False)
            primitives2_inside = contour2.extract_with_points(contours_intersections[0],
                                                              contours_intersections[1], True)
            primitives2_outside = contour2.extract_with_points(contours_intersections[0],
                                                               contours_intersections[1], False)
            if contour2.point_belongs(primitives1_inside[0].middle_point()):
                resulting_primitives.extend(primitives1_outside)
            else:
                resulting_primitives.extend(primitives1_inside)
            if self.point_belongs(primitives2_inside[0].middle_point()):
                resulting_primitives.extend(primitives2_outside)
            else:
                resulting_primitives.extend(primitives2_inside)
            return [Contour2D(resulting_primitives).order_contour()]

        return self.merge_with(contour2)

    def cut_by_wire(self, wire: Wire2D):
        """
        Cut a contour2d with a wire2d and return a list of contours 2d.

        :param wire: volmdlr.wires.Wire2D
        :rtype: list[volmdlr.wires.Contour2D]

        :param wire: volmdlr.wires.Wire2D.
        :return: contours2d : list[volmdlr.wires.Contour2D].
        """

        intersections = self.wire_crossings(wire)  # crossings OR intersections (?)
        if not intersections or len(intersections) < 2:
            return [self]
        points_intersections = []
        for intersection, _ in intersections:
            if intersection not in points_intersections:
                points_intersections.append(intersection)
        if len(points_intersections) % 2 != 0:
            raise NotImplementedError(
                f'{len(points_intersections)} intersections not supported yet')

        # points_intersections = [point for point, prim in intersections]

        sorted_points = wire.sort_points_along_wire(points_intersections)
        list_contours = []
        contour_to_cut = self
        cutting_points_counter = 0
        while cutting_points_counter != len(sorted_points):

            point1 = sorted_points[cutting_points_counter]
            point2 = sorted_points[cutting_points_counter + 1]

            closing_wire = wire.extract_with_points(point1, point2, True)

            for point in points_intersections:
                if point not in [point1, point2] and Wire2D(closing_wire).point_over_wire(point):
                    closing_wire = wire.extract_with_points(point1, point2, False)
                    break

            closing_wire_prim = [closing_w for closing_w in closing_wire if closing_w]
            closing_contour = Contour2D(closing_wire_prim)
            contour1, contour2 = contour_to_cut.get_divided_contours(point1,
                                                                     point2,
                                                                     closing_contour,
                                                                     True)

            if sorted_points.index(point1) + 2 < len(sorted_points) - 1:
                if contour1.point_over_contour(
                        sorted_points[sorted_points.index(point1) + 2]):
                    contour_to_cut = contour1
                    list_contours.append(contour2)
                elif contour2.point_over_contour(
                        sorted_points[sorted_points.index(point1) + 2]):
                    contour_to_cut = contour2
                    list_contours.append(contour1)
            else:
                list_contours.extend([contour1, contour2])
            cutting_points_counter += 2

        return list_contours

    def get_furthest_point_to_point2(self, point2):
        """
        Search the furthest point from self to point2. It only considers the start or end or primitives.

        :param point2: other point.
        :return: the furthest point.
        """
        furthest_point = self.primitives[0].start
        furthest_distance = point2.point_distance(self.primitives[0].start)
        for prim in self.primitives:
            distance = point2.point_distance(prim.end)
            if distance > furthest_distance:
                furthest_distance = distance
                furthest_point = prim.end
        return furthest_point

    def closest_point_to_point2(self, point2):
        """
        Search the closest point from self to point2. It only considers the start or end or primitives.

        :param point2: other point.
        :return: the closest point to point2.
        """
        closest_point = self.primitives[0].start
        closest_distance = point2.point_distance(self.primitives[0].start)
        for prim in self.primitives:
            distance = point2.point_distance(prim.end)
            if distance < closest_distance:
                closest_distance = distance
                closest_point = prim.end
        return closest_point


class ClosedPolygonMixin:
    """
    Abstract class for ClosedPolygon, storing methods used by ClosedPolygon2D and ClosedPolygon3D.

    """

    def get_lengths(self):
        """
        Gets line segment lengths.

        """
        list_ = []
        for line_segment in self.line_segments:
            list_.append(line_segment.length())
        return list_

    def length(self):
        """
        Polygon length.

        :return: polygon length.
        """
        return sum(self.get_lengths())

    def min_length(self):
        """
        Gets the minimal length for a line segment in the polygon.

        """
        return min(self.get_lengths())

    def max_length(self):
        """
        Gets the minimal length for a line segment in the polygon.

        """
        return max(self.get_lengths())

    def edge_statistics(self):
        distances = []
        for i, point in enumerate(self.points):
            if i != 0:
                distances.append(point.point_distance(self.points[i - 1]))
        mean_distance = mean(distances)
        std = npy.std(distances)
        return mean_distance, std

    def simplify_polygon(self, min_distance: float = 0.01,
                         max_distance: float = 0.05, angle: float = 15):
        points = [self.points[0]]
        previous_point = None
        for point in self.points[1:]:
            distance = point.point_distance(points[-1])
            if distance > min_distance:
                if distance > max_distance:
                    number_segmnts = round(distance / max_distance) + 2
                    for n in range(number_segmnts):
                        new_point = points[-1] + (point - points[-1]) * (
                                n + 1) / number_segmnts
                        if new_point.point_distance(points[-1]) > max_distance:
                            points.append(new_point)
                else:
                    if not volmdlr.core.point_in_list(point, points):
                        points.append(point)
            if len(points) > 1:
                vector1 = points[-1] - points[-2]
                vector2 = point - points[-2]
                cos = vector1.dot(vector2) / (vector1.norm() * vector2.norm())
                cos = math.degrees(math.acos(round(cos, 6)))
                if abs(cos) > angle:
                    if not volmdlr.core.point_in_list(previous_point, points):
                        points.append(previous_point)
                    if not volmdlr.core.point_in_list(point, points):
                        points.append(point)
            if len(points) > 2:
                vector1 = points[-2] - points[-3]
                vector2 = points[-1] - points[-3]
                cos = vector1.dot(vector2) / (vector1.norm() * vector2.norm())
                cos = math.degrees(math.acos(round(cos, 6)))
                if points[-3].point_distance(points[-2]) < min_distance and cos < angle:
                    points = points[:-2] + [points[-1]]
            previous_point = point
        if points[0].point_distance(points[-1]) < min_distance:
            points.remove(points[-1])

        if math.isclose(volmdlr.wires.ClosedPolygon2D(points).area(), 0.0, abs_tol=1e-6):
            return self

        return self.__class__(points)

    @property
    def line_segments(self):
        if not self._line_segments:
            self._line_segments = self.get_line_segments()
        return self._line_segments

    def get_line_segments(self):
        raise NotImplementedError(
            f"get_line_segments method must be overloaded by {self.__class__.__name__}")


class ClosedPolygon2D(Contour2D, ClosedPolygonMixin):
    """
    A collection of points, connected by line segments, following each other.

    """
    _non_serializable_attributes = ['line_segments', 'primitives',
                                    'basis_primitives']

    def __init__(self, points: List[volmdlr.Point2D], name: str = ''):
        self.points = points
        self._line_segments = None

        Contour2D.__init__(self, self.line_segments, name)

    def copy(self, *args, **kwargs):
        points = [point.copy() for point in self.points]
        return ClosedPolygon2D(points, self.name)

    def __hash__(self):
        return sum(hash(point) for point in self.points)

    def __eq__(self, other_):
        if not isinstance(other_, self.__class__):
            return False
        equal = True
        for point, other_point in zip(self.points, other_.points):
            equal = (equal and point == other_point)
        return equal

    def area(self):
        # TODO: performance: cache number of points
        if len(self.points) < 3:
            return 0.

        x = [point.x for point in self.points]
        y = [point.y for point in self.points]

        x1 = [x[-1]] + x[0:-1]
        y1 = [y[-1]] + y[0:-1]
        return 0.5 * abs(sum(i * j for i, j in zip(x, y1))
                         - sum(i * j for i, j in zip(y, x1)))
        # return 0.5 * npy.abs(
        #     npy.dot(x, npy.roll(y, 1)) - npy.dot(y, npy.roll(x, 1)))

    def center_of_mass(self):
        lp = len(self.points)
        if lp == 0:
            return volmdlr.O2D
        if lp == 1:
            return self.points[0]
        if lp == 2:
            return 0.5 * (self.points[0] + self.points[1])

        x = [point.x for point in self.points]
        y = [point.y for point in self.points]

        xi_xi1 = x + npy.roll(x, -1)
        yi_yi1 = y + npy.roll(y, -1)
        xi_yi1 = npy.multiply(x, npy.roll(y, -1))
        xi1_yi = npy.multiply(npy.roll(x, -1), y)

        signed_area = 0.5 * npy.sum(xi_yi1 - xi1_yi)  # signed area!
        if not math.isclose(signed_area, 0, abs_tol=1e-08):
            cx = npy.sum(npy.multiply(xi_xi1, (xi_yi1 - xi1_yi))) / 6. / signed_area
            cy = npy.sum(npy.multiply(yi_yi1, (xi_yi1 - xi1_yi))) / 6. / signed_area
            return volmdlr.Point2D(cx, cy)

        self.plot()
        raise NotImplementedError

    def barycenter(self):
        """
        Calculates the geometric center of the polygon, which is the average position of all the points in it.

        :rtype: volmdlr.Point2D
        """
        barycenter1_2d = self.points[0]
        for point in self.points[1:]:
            barycenter1_2d += point
        return barycenter1_2d / len(self.points)

    def point_belongs(self, point, include_edge_points: bool = False):
        """
        Ray casting algorithm copied from internet.
        """
        return polygon_point_belongs((point.x, point.y), [(point_.x, point_.y) for point_ in self.points],
                                     include_edge_points=include_edge_points)

    def second_moment_area(self, point):
        second_moment_area_x, second_moment_area_y, second_moment_area_xy = 0., 0., 0.
        for point_i, point_j in zip(self.points, self.points[1:] + [self.points[0]]):
            xi, yi = point_i - point
            xj, yj = point_j - point
            second_moment_area_x += (yi ** 2 + yi * yj + yj ** 2) * (xi * yj - xj * yi)
            second_moment_area_y += (xi ** 2 + xi * xj + xj ** 2) * (xi * yj - xj * yi)
            second_moment_area_xy += (xi * yj + 2 * xi * yi + 2 * xj * yj + xj * yi) * (
                    xi * yj - xj * yi)
        if second_moment_area_x < 0:
            second_moment_area_x = - second_moment_area_x
            second_moment_area_y = - second_moment_area_y
            second_moment_area_xy = - second_moment_area_xy
        return second_moment_area_x / 12., second_moment_area_y / 12., second_moment_area_xy / 24.

    def get_line_segments(self):
        lines = []
        if len(self.points) > 1:
            for point1, point2 in zip(self.points,
                                      list(self.points[1:]) + [self.points[0]]):
                if not point1.is_close(point2):
                    lines.append(volmdlr.edges.LineSegment2D(point1, point2))
        return lines

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        ClosedPolygon2D rotation.

        :param center: rotation center
        :param angle: angle rotation
        :return: a new rotated ClosedPolygon2D
        """
        return ClosedPolygon2D(
            [point.rotation(center, angle) for point in self.points])

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        Line2D rotation, Object is updated in-place.

        :param center: rotation center
        :param angle: rotation angle
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for point in self.points:
            point.rotation_inplace(center, angle)

    def translation(self, offset: volmdlr.Vector2D):
        """
        ClosedPolygon2D translation.

        :param offset: translation vector
        :return: A new translated ClosedPolygon2D
        """
        return ClosedPolygon2D(
            [point.translation(offset) for point in self.points])

    def translation_inplace(self, offset: volmdlr.Vector2D):
        """
        ClosedPolygon2D translation. Object is updated in-place.

        :param offset: translation vector
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for point in self.points:
            point.translation_inplace(offset)

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        return self.__class__([point.frame_mapping(frame, side) for point in self.points])

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for point in self.points:
            point.frame_mapping_inplace(frame, side)

    def polygon_distance(self,
                         polygon: 'ClosedPolygon2D'):
        point_zero = self.points[0]
        distance = []
        for point in polygon.points:
            distance.append(point_zero.point_distance(point))
        index = distance.index(min(distance))
        return distance[index]

    @cached_property
    def is_trigo(self):
        if len(self.points) < 3:
            return True

        angle = 0.
        for ls1, ls2 in zip(self.line_segments,
                            self.line_segments[1:] + [self.line_segments[0]]):
            u = ls2.unit_direction_vector()
            x = u.dot(ls1.unit_direction_vector())
            y = u.dot(ls1.normal_vector())
            angle += math.atan2(y, x)
        return angle > 0

    def delaunay_triangulation(self):
        points = self.points
        new_points = []
        delaunay_triangles = []
        # ax=plt.subplot()
        for point in points:
            new_points.append([point[0], point[1]])

        delaunay = npy.array(new_points)

        tri = Delaunay(delaunay)

        for simplice in delaunay[tri.simplices]:
            triangle = Triangle2D(volmdlr.Point2D(simplice[0]),
                                  volmdlr.Point2D(simplice[1]),
                                  volmdlr.Point2D(simplice[2]))
            delaunay_triangles.append(triangle)

        return delaunay_triangles

    def offset(self, offset):
        x_min, x_max, y_min, y_max = self.bounding_rectangle.bounds()

        max_offset_len = min(x_max - x_min, y_max - y_min) / 2
        if offset <= -max_offset_len:
            print('Inadapted offset, '
                  'polygon might turn over. Offset must be greater than',
                  -max_offset_len)
            raise ValueError('inadapted offset')
        nb_points = len(self.points)
        vectors = []
        for i in range(nb_points - 1):
            v1 = self.points[i + 1] - self.points[i]
            v2 = self.points[i] - self.points[i + 1]
            v1.normalize()
            v2.normalize()
            vectors.append(v1)
            vectors.append(v2)

        v1 = self.points[0] - self.points[-1]
        v2 = self.points[-1] - self.points[0]
        v1.normalize()
        v2.normalize()
        vectors.append(v1)
        vectors.append(v2)

        offset_vectors = []
        offset_points = []

        for i in range(nb_points):

            # check = False
            vector_i = vectors[2 * i - 1] + vectors[2 * i]
            if vector_i == volmdlr.Vector2D(0, 0):
                vector_i = vectors[2 * i]
                vector_i = vector_i.normal_vector()
                offset_vectors.append(vector_i)
            else:
                vector_i.normalize()
                if vector_i.dot(vectors[2 * i - 1].normal_vector()) > 0:
                    vector_i = - vector_i
                    # check = True
                offset_vectors.append(vector_i)

            normal_vector1 = - vectors[2 * i - 1].normal_vector()
            normal_vector2 = vectors[2 * i].normal_vector()
            normal_vector1.normalize()
            normal_vector2.normalize()
            alpha = math.acos(normal_vector1.dot(normal_vector2))

            offset_point = self.points[i] + offset / math.cos(alpha / 2) * \
                (-offset_vectors[i])

            # ax=self.plot()
            # offset_point.plot(ax=ax, color='g')

            # if self.point_belongs(offset_point):
            #     offset_point = self.points[i] + offset / math.cos(alpha / 2) * \
            #                    (-offset_vectors[i])

            offset_points.append(offset_point)

            # self.points[i].plot(ax=ax, color='b')
            # offset_point.plot(ax=ax, color='r')

        return self.__class__(offset_points)

    def point_border_distance(self, point, return_other_point=False):
        """
        Compute the distance to the border distance of polygon.

        Output is always positive, even if the point belongs to the polygon.
        """
        d_min, other_point_min = self.line_segments[0].point_distance(
            point, return_other_point=True)
        for line in self.line_segments[1:]:
            dist_, other_point = line.point_distance(
                point, return_other_point=True)
            if dist_ < d_min:
                d_min = dist_
                other_point_min = other_point
        if return_other_point:
            return d_min, other_point_min
        return d_min

    def self_intersects(self):
        """
        Determines if a polygon self intersects using the Bentley-Ottmann algorithm.

        :return: True if the polygon self intersects, False otherwis. If True, returns two
            intersecting line segments as LineSegment2D objects. If False, returns two None values;
        :rtype: Tuple[bool, Union[volmdlr.edges.LineSegment2D, None], Union[volmdlr.edges.LineSegment2D, None]]
        """
        epsilon = 0
        segments = self._get_segments()

        for segment1 in segments:
            for segment2 in segments:
                if segment1 == segment2:
                    continue
                if self._segments_intersect(segment1, segment2, epsilon):
                    return True, segment1, segment2

        return False, None, None

    def _get_segments(self):
        """
        Helper function for self_intersects that generates segments for the Bentley-Ottmann algorithm.

        :return: A list of tuples representing the segments between consecutive edges.
        :rtype: List[Tuple[int, int]]
        """
        # Sort the points along ascending x for the Sweep Line method
        sorted_index = sorted(range(len(self.points)), key=lambda p: (self.points[p][0], self.points[p][1]))
        nb = len(sorted_index)
        segments = []

        for i, index in enumerate(sorted_index):
            # Stock the segments between 2 consecutive edges
            # Ex: for the ABCDE polygon, if Sweep Line is on C, the segments
            #   will be (C,B) and (C,D)
            if index - 1 < 0:
                segments.append((index, nb - 1))
            else:
                segments.append((index, sorted_index[i - 1]))
            if index >= len(self.points) - 1:
                segments.append((index, 0))
            else:
                segments.append((index, sorted_index[i + 1]))

        return segments

    def _segments_intersect(self, segment1, segment2, epsilon):
        """
        Helper function for self_intersects that determines if any segments in a list intersect.

        :param segment1: A tuple representing the index of the start and end point of the segments.
        :type segment1: Tuple[int, int]
        :param segment2: A tuple representing the index of the start and end point of the segments.
        :type segment2: Tuple[int, int]
        :param epsilon: A small positive value for numerical stability.
        :type epsilon: float
        :return: True if any segments intersect, False otherwise.
        :rtype: bool
        """
        line1 = volmdlr.edges.LineSegment2D(self.points[segment1[0]], self.points[segment1[1]])
        line2 = volmdlr.edges.LineSegment2D(self.points[segment2[0]], self.points[segment2[1]])
        point, param_a, param_b = volmdlr.Point2D.line_intersection(line1, line2, True)
        if point is not None and 0 + epsilon <= param_a <= 1 - epsilon and 0 + epsilon <= param_b <= 1 - epsilon:
            return True
        return False

    @classmethod
    def points_convex_hull(cls, points):
        if len(points) < 3:
            return None

        points_hull = [point.copy() for point in points]

        _, pos_ymax = argmax([point.y for point in points_hull])
        point_start = points_hull[pos_ymax]
        hull = [point_start]

        barycenter = points_hull[0]
        for point in points_hull[1:]:
            barycenter += point
        barycenter = barycenter / (len(points_hull))
        # second point of hull
        theta = []
        remaining_points = points_hull
        del remaining_points[pos_ymax]

        vec1 = point_start - barycenter
        for point in remaining_points:
            vec2 = point - point_start
            theta_i = -volmdlr.geometry.clockwise_angle(vec1, vec2)
            theta.append(theta_i)

        min_theta, posmin_theta = argmin(theta)
        next_point = remaining_points[posmin_theta]
        hull.append(next_point)
        del remaining_points[posmin_theta]
        # Adding first point to close the loop at the end
        remaining_points.append(hull[0])

        initial_vector = vec1.copy()
        total_angle = 0
        while not next_point.is_close(point_start):
            vec1 = next_point - hull[-2]
            theta = []
            for point in remaining_points:
                vec2 = point - next_point
                theta_i = -volmdlr.geometry.clockwise_angle(vec1, vec2)
                theta.append(theta_i)

            min_theta, posmin_theta = argmin(theta)
            if math.isclose(min_theta, -2 * math.pi, abs_tol=1e-6) \
                    or math.isclose(min_theta, 0, abs_tol=1e-6):
                if remaining_points[posmin_theta] == point_start:
                    break

            else:
                next_point = remaining_points[posmin_theta]

                vec_next_point = next_point - barycenter
                total_angle += (2 * math.pi - volmdlr.geometry.clockwise_angle(initial_vector, vec_next_point))

                if total_angle > 2 * math.pi:
                    break
                initial_vector = vec_next_point

                hull.append(next_point)

            del remaining_points[posmin_theta]

        hull.pop()

        return cls(hull)

    @classmethod
    def concave_hull(cls, points, concavity, scale_factor):
        """
        Calculates the concave hull from a cloud of points.

        i.e., it Unites all points under the smallest possible area.

        :param points: list of points corresponding to the cloud of points
        :type points: class: 'volmdlr.Point2D'
        :param concavity: Sets how sharp the concave angles can be. It goes from -1 (not concave at all. in fact,
                          the hull will be left convex) up to +1 (very sharp angles can occur. Setting concavity to
                          +1 might result in 0º angles!) concavity is defined as the cosine of the concave angles.
        :type concavity: float
        :param scale_factor: Sets how big is the area where concavities are going to be searched.
                             The bigger, the more sharp the angles can be. Setting it to a very high value might
                             affect the performance of the program.
                             This value should be relative to how close to each other the points to be connected are.
        :type scale_factor: float

        """

        def get_nearby_points(line, points, scale_factor):
            points_hull = [point.copy() for point in points]

            # print('i enter here')
            nearby_points = []
            line_midpoint = 0.5 * (line.start + line.end)
            # print(line_midpoint)
            tries = 0
            n = 5
            bounding_box = [line_midpoint.x - line.length() / 2,
                            line_midpoint.x + line.length() / 2,
                            line_midpoint.y - line.length() / 2,
                            line_midpoint.y + line.length() / 2]
            boundary = [int(bounding / scale_factor) for bounding in
                        bounding_box]
            while tries < n and len(nearby_points) == 0:
                for point in points_hull:
                    if not ((
                                    point.x == line.start.x and point.y == line.start.y) or (
                                    point.x == line.end.x and point.y == line.end.y)):
                        point_x_rel_pos = int(point.x / scale_factor)
                        point_y_rel_pos = int(point.y / scale_factor)
                        if boundary[0] <= point_x_rel_pos <= boundary[1] \
                                and point_y_rel_pos >= boundary[2] \
                                and point_y_rel_pos <= boundary[3]:
                            nearby_points.append(point)

                scale_factor *= 4 / 3
                tries += 1

            return nearby_points

        def line_colides_with_hull(line, concave_hull):
            for hull_line in concave_hull:
                if not line.start.is_close(hull_line.start) and not line.start.is_close(hull_line.end) and \
                        not line.end.is_close(hull_line.start) and not line.end.is_close(hull_line.end):
                    if line.line_intersections(hull_line.to_line()):
                        return True
            return False

        def get_divided_line(line, nearby_points, hull_concave_edges, concavity):
            divided_line = []
            ok_middle_points = []
            list_cossines = []
            for middle_point in nearby_points:
                vect1 = line.start - middle_point
                vect2 = line.end - middle_point
                if middle_point.is_close(line.start) or middle_point.is_close(line.end):
                    continue
                cos = round(vect1.dot(vect2) / (vect1.norm() * vect2.norm()),
                            4)
                if cos < concavity:
                    new_line_a = volmdlr.edges.LineSegment2D(start=line.start, end=middle_point)
                    new_line_b = volmdlr.edges.LineSegment2D(start=middle_point, end=line.end)
                    if not (line_colides_with_hull(line=new_line_a,
                                                   concave_hull=hull_concave_edges) and line_colides_with_hull(
                            line=new_line_b, concave_hull=hull_concave_edges)):
                        ok_middle_points.append(middle_point)
                        list_cossines.append(cos)
            if len(ok_middle_points) > 0:
                #  We want the middle-point to be the one with the widest angle (smallest cossine)
                min_cossine_index = list_cossines.index(min(list_cossines))
                divided_line.append(volmdlr.edges.LineSegment2D(line.start,
                                                                ok_middle_points[
                                                                    min_cossine_index]))
                divided_line.append(volmdlr.edges.LineSegment2D(
                    ok_middle_points[min_cossine_index], line.end))
            return divided_line

        hull_convex_edges = cls.points_convex_hull(points).line_segments
        hull_convex_edges.sort(key=lambda x: x.length(), reverse=True)
        hull_concave_edges = []
        hull_concave_edges.extend(hull_convex_edges)
        hull_points = list({point for line in hull_concave_edges for point in [line[0], line[1]]})
        unused_points = []
        for point in points:
            if not volmdlr.core.point_in_list(point, hull_points):
                unused_points.append(point)

        a_line_was_divided_in_the_iteration = True
        line = None
        divided_line = None
        while a_line_was_divided_in_the_iteration:
            a_line_was_divided_in_the_iteration = False
            for line in hull_concave_edges:
                nearby_points = get_nearby_points(line, unused_points,
                                                  scale_factor)
                divided_line = get_divided_line(line, nearby_points,
                                                hull_concave_edges, concavity)
                if len(divided_line) > 0:
                    a_line_was_divided_in_the_iteration = True
                    unused_points.remove(divided_line[0].end)
                    break
            else:
                continue
            hull_concave_edges.remove(line)
            hull_concave_edges.extend(divided_line)

            hull_concave_edges.sort(key=lambda x: x.length(), reverse=True)

        polygon_points = [(line.start, line.end) for line in hull_concave_edges]

        points = [polygon_points[0][0], polygon_points[0][1]]
        polygon_points.remove((polygon_points[0][0], polygon_points[0][1]))
        while True:
            if not polygon_points:
                break
            point1, point2 = None, None
            for point1, point2 in polygon_points:
                if point1 == points[-1] and point2 not in points:
                    points.append(point2)
                    break
                if point2 == points[-1] and point1 not in points:
                    points.append(point1)
                    break
            polygon_points.remove((point1, point2))

        return cls(points)  # , nearby_points

    @classmethod
    def convex_hull_points(cls, points):
        """
        Uses the scipy method ConvexHull to calculate the convex hull from a cloud of points.

        """

        points_hull = [point.copy() for point in points]

        numpy_points = npy.array([(point.x, point.y) for point in points_hull])
        hull = ConvexHull(numpy_points)
        polygon_points = []
        for simplex in hull.simplices:
            polygon_points.append((points_hull[simplex[0]], points_hull[simplex[1]]))

        points_hull = [polygon_points[0][0], polygon_points[0][1]]
        polygon_points.remove((polygon_points[0][0], polygon_points[0][1]))

        while True:
            if not polygon_points:
                break
            point1, point2 = None, None
            for point1, point2 in polygon_points:
                if point1 == points_hull[-1]:
                    points_hull.append(point2)
                    break
                if point2 == points_hull[-1]:
                    points_hull.append(point1)
                    break
            polygon_points.remove((point1, point2))

        points_hull.pop(-1)

        # the first point is the one with the lowest x value
        i_min = 0
        min_x = points_hull[0].x
        for i, point in enumerate(points_hull):
            if point.x < min_x:
                min_x = point.x
                i_min = i

        points_hull = points_hull[i_min:] + points_hull[:i_min]

        # we make sure that the points are ordered in the trigonometric direction
        if points_hull[0].y < points_hull[1].y:
            points_hull.reverse()

        return cls(points_hull)

    def to_3d(self, plane_origin, x, y):
        """
        Transforms a ClosedPolygon2D into an ClosedPolygon3D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: ClosedPolygon3D.
        """
        points3d = [point.to_3d(plane_origin, x, y) for point in self.points]
        return ClosedPolygon3D(points3d)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle(), point_numbering=False,
             fill=False, fill_color='w'):
        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect('equal')

        if fill:
            ax.fill([point[0] for point in self.points], [point[1] for point in self.points],
                    facecolor=fill_color)
        for line_segment in self.line_segments:
            line_segment.plot(ax=ax, edge_style=edge_style)

        if edge_style.plot_points or point_numbering:
            for point in self.points:
                point.plot(ax=ax, color=edge_style.color, alpha=edge_style.alpha)

        if point_numbering:
            for ip, point in enumerate(self.points):
                ax.text(*point, f'point {ip + 1}', ha='center', va='top')

        if edge_style.equal_aspect:
            ax.set_aspect('equal')
        else:
            ax.set_aspect('auto')

        ax.margins(0.1)
        plt.show()

        return ax

    def triangulation(self, tri_opt: str = 'pd'):
        """
        Perform triangulation on the polygon.

        To detail documentation, please refer to https://rufat.be/triangle/API.html

        :param tri_opt: (Optional) Triangulation preferences.
        :type tri_opt: str
        :return: A 2D mesh.
        :rtype: :class:`vmd.DisplayMesh2D`
        """
        # Converting points to nodes for performance
        nodes = [vmd.Node2D.from_point(point) for point in self.points]
        vertices = [(point.x, point.y) for point in nodes]
        n = len(nodes)
        segments = [(i, i + 1) for i in range(n - 1)]
        segments.append((n - 1, 0))

        tri = {'vertices': npy.array(vertices).reshape((-1, 2)),
               'segments': npy.array(segments).reshape((-1, 2)),
               }
        t = triangulate(tri, tri_opt)
        triangles = t['triangles'].tolist()
        np = t['vertices'].shape[0]
        points = [vmd.Node2D(*t['vertices'][i, :]) for i in range(np)]
        return vmd.DisplayMesh2D(points, triangles=triangles)

    def grid_triangulation_points(self, number_points_x: int = 25, number_points_y: int = 25):
        """
        Use an n by m grid to triangulize the contour.

        :param number_points_x: Number of discretization points in x direction.
        :type number_points_x: int
        :param number_points_y: Number of discretization points in y direction.
        :type number_points_y: int
        :return: Discretization data.
        :rtype: list
        """
        x_min, x_max, y_min, y_max = self.bounding_rectangle.bounds()

        n = number_points_x + 2
        m = number_points_y + 2

        x = npy.linspace(x_min, x_max, num=n)
        y = npy.linspace(y_min, y_max, num=m)

        grid_point_index = {}

        polygon_points = {vmd.Node2D.from_point(point) for point in self.points}
        points = []
        for i, xi in enumerate(x):
            for j, yi in enumerate(y):
                point = vmd.Node2D(xi, yi)
                if self.point_belongs(point, include_edge_points=True) and point not in polygon_points:
                    grid_point_index[(i, j)] = point
                    points.append(point)

        return points, x, y, grid_point_index

    def ear_clipping_triangulation(self):
        """
        Computes the triangulation of the polygon using ear clipping algorithm.

        Note: triangles have been inverted for a better rendering in babylonjs
        """
        # Converting to nodes for performance
        nodes = [vmd.Node2D.from_point(point) for point in self.points]

        initial_point_to_index = {point: i for i, point in enumerate(nodes)}
        triangles = []

        remaining_points = nodes[:]

        number_remaining_points = len(remaining_points)
        while number_remaining_points > 3:
            current_polygon = ClosedPolygon2D(remaining_points)

            found_ear = False
            for point1, point2, point3 in zip(remaining_points,
                                              remaining_points[1:] + remaining_points[0:1],
                                              remaining_points[2:] + remaining_points[0:2]):
                if not point1.is_close(point3):
                    line_segment = volmdlr.edges.LineSegment2D(point1, point3)

                # Checking if intersections does not contain the vertices
                # of line_segment
                intersect = False
                intersections = current_polygon.linesegment_intersections(line_segment)
                if intersections:
                    for inter in intersections:
                        if not volmdlr.core.point_in_list(inter[0], [line_segment.start, line_segment.end]):
                            intersect = True
                            break

                if not intersect:
                    if current_polygon.point_belongs(line_segment.middle_point()):

                        triangles.append((initial_point_to_index[point1],
                                          initial_point_to_index[point3],
                                          initial_point_to_index[point2]))
                        remaining_points.remove(point2)
                        number_remaining_points -= 1
                        found_ear = True

                        # Rolling the remaining list
                        if number_remaining_points > 4:
                            deq = deque(remaining_points)
                            # random.randint(1, number_remaining_points-1))
                            deq.rotate(int(0.3 * number_remaining_points))
                            remaining_points = list(deq)

                        break

            # Searching for a flat ear
            if not found_ear:
                remaining_polygon = ClosedPolygon2D(remaining_points)
                if remaining_polygon.area() > 0.:

                    found_flat_ear = False
                    for point1, point2, point3 in zip(remaining_points,
                                                      remaining_points[1:] + remaining_points[0:1],
                                                      remaining_points[2:] + remaining_points[0:2]):
                        triangle = Triangle2D(point1, point2, point3)
                        if math.isclose(triangle.area(), 0, abs_tol=1e-8):
                            remaining_points.remove(point2)
                            found_flat_ear = True
                            break

                    if not found_flat_ear:
                        print('Warning : There are no ear in the polygon, it seems malformed: skipping triangulation')
                        return vmd.DisplayMesh2D(nodes, triangles)
                else:
                    return vmd.DisplayMesh2D(nodes, triangles)

        if len(remaining_points) == 3:
            point1, point2, point3 = remaining_points
            triangles.append((initial_point_to_index[point1],
                              initial_point_to_index[point3],
                              initial_point_to_index[point2]))

        return vmd.DisplayMesh2D(nodes, triangles)

    def simplify(self, min_distance: float = 0.01, max_distance: float = 0.05):
        return ClosedPolygon2D(self.simplify_polygon(min_distance=min_distance,
                                                     max_distance=max_distance).points)

    def line_intersecting_closing_point(self, crossing_point):
        """
        Finds closing point for the sewing method using intersection of lines drawn from the barycenter.

        returns the closing point.
        """
        vec_dir = crossing_point.copy()
        vec_dir.normalize()

        line = volmdlr.edges.LineSegment2D(volmdlr.O2D,
                                           crossing_point + vec_dir * 5)
        # line.plot(ax=ax2d, color='b')

        point_intersections = {}
        for line_segment in self.line_segments:
            point_intersection = line_segment.linesegment_intersections(
                line)
            if point_intersection:
                point_intersections[line_segment] = point_intersection[
                    0]
            else:
                if line.point_belongs(line_segment.start):
                    point_intersections[line_segment] = line_segment.start
                if line.point_belongs(line_segment.end):
                    point_intersections[line_segment] = line_segment.end
        point_distance = list(point_intersections.values())[
            0].point_distance(crossing_point)
        point_intersection = list(point_intersections.values())[0]
        line_segment = list(point_intersections.keys())[0]
        for line, point in list(point_intersections.items())[1:]:
            dist = crossing_point.point_distance(point)
            if dist < point_distance:
                point_distance = dist
                point_intersection = point
                line_segment = line

        # point_intersection.plot(ax=ax2d)

        if point_intersection.point_distance(
                    line_segment.start) < point_intersection.point_distance(
                line_segment.end):
            closing_point = line_segment.start
        else:
            closing_point = line_segment.end

        return closing_point

    def point_in_polygon(self):
        """
        In case the barycenter of the polygon is outside, this method finds another point inside the polygon.

        """
        barycenter = self.barycenter()
        if self.point_belongs(barycenter):
            return barycenter
        intersetions1 = {}
        linex_pos = volmdlr.edges.LineSegment2D(volmdlr.O2D, volmdlr.X2D * 5)
        linex_neg = volmdlr.edges.LineSegment2D(volmdlr.O2D, -volmdlr.X2D * 5)
        liney_pos = volmdlr.edges.LineSegment2D(volmdlr.O2D, volmdlr.Y2D * 5)
        liney_neg = volmdlr.edges.LineSegment2D(volmdlr.O2D, -volmdlr.Y2D * 5)
        for line in [linex_pos, linex_neg, liney_pos, liney_neg]:
            intersections = []
            for line_segment in self.line_segments:
                point_intersection = line_segment.linesegment_intersections(
                    line)
                intersections.extend(point_intersection)
                if not point_intersection:
                    if line.point_belongs(line_segment.start):
                        intersections.append(line_segment.start)
                    if line.point_belongs(line_segment.end):
                        intersections.append(line_segment.end)
            intersetions1[line] = intersections[:]
        for i, value in enumerate(intersetions1.values()):
            if not value:
                if i % 2 == 0:
                    if len(list(intersetions1.values())[i + 1]) == 2:
                        translation1 = (list(intersetions1.values())[i + 1][0] +
                                        list(intersetions1.values())[
                                            i + 1][1]) * 0.5
                        break
                if i % 2 != 0:
                    if len(list(intersetions1.values())[i - 1]) == 2:
                        translation1 = (list(intersetions1.values())[i - 1][0]
                                        + list(intersetions1.values())[i - 1][1]) * 0.5
                        break

        return translation1

    def repositioned_polygon(self, x, y):
        linex = volmdlr.edges.LineSegment2D(-x.to_2d(volmdlr.O2D, x, y),
                                            x.to_2d(volmdlr.O2D, x, y))
        way_back = volmdlr.O3D
        barycenter = self.barycenter()
        if not self.point_belongs(barycenter):
            barycenter1_2d = self.point_in_polygon()
            self.translation_inplace(-barycenter1_2d)
            way_back = barycenter1_2d.to_3d(volmdlr.O3D, x, y)
        else:
            inters = self.linesegment_intersections(linex)
            distance = inters[0][0].point_distance(inters[-1][0])
            if distance / 2 > 3 * min(
                    self.point_distance(inters[0][0]),
                    self.point_distance(inters[-1][0])):
                mid_point = (inters[0][0] + inters[-1][0]) * 0.5
                self.translation(-mid_point)
                way_back = mid_point.to_3d(volmdlr.O3D, x, y)

        return self, way_back

    def get_possible_sewing_closing_points(self, polygon2, polygon_primitive,
                                           line_segment1: None, line_segment2: None):
        """
        Searches all possibles closing points available for the given primitive.

        """
        middle_point = polygon_primitive.middle_point()
        if line_segment1 is None and line_segment2 is None:
            normal_vector = polygon_primitive.unit_normal_vector()
            line_segment1 = volmdlr.edges.LineSegment2D(middle_point,
                                                        middle_point - normal_vector)
            line_segment2 = volmdlr.edges.LineSegment2D(middle_point,
                                                        middle_point + normal_vector)

        line_intersections = {line_segment1: [], line_segment2: []}
        for line_segment in [line_segment1, line_segment2
                             ]:
            inter_points = []
            for prim in polygon2.line_segments + self.line_segments[
                                                 :self.line_segments.index(
                                                     polygon_primitive)] + self.line_segments[
                                                                           self.line_segments.index(
                                                                               polygon_primitive) + 1:]:
                inters = prim.linesegment_intersections(line_segment)
                if inters:
                    line_intersections[line_segment].append((inters[0], prim))
                    inter_points.append(inters[0])
                elif line_segment.point_belongs(prim.start, 1e-7):
                    if not volmdlr.core.point_in_list(prim.start, inter_points):
                        line_intersections[line_segment].append((prim.start, prim))
                        inter_points.append(prim.start)
                elif line_segment.point_belongs(prim.end, 1e-7):
                    if not volmdlr.core.point_in_list(prim.end, inter_points):
                        line_intersections[line_segment].append((prim.end, prim))
                        inter_points.append(prim.end)
                elif prim.point_belongs(middle_point, 1e-7):
                    line_intersections[line_segment].append((prim.middle_point(), prim))
                    inter_points.append(prim.middle_point())
        return line_intersections

    def select_farthest_sewing_closing_point(self,
                                             line_segment: volmdlr.edges.LineSegment2D,
                                             polygon_primitive,
                                             possible_closing_points):
        """
        Searches the closest sewing closing point available.

        """
        closing_point = volmdlr.O2D
        middle_point = polygon_primitive.middle_point()
        distance = 0
        for intr_list in possible_closing_points:
            if intr_list[1] not in self.line_segments:
                dist = intr_list[0].point_distance(line_segment.start)
                if dist > distance:
                    distance = dist
                    closing_point = (intr_list[1].start if
                                     intr_list[0].point_distance(
                                         intr_list[1].start) <
                                     intr_list[0].point_distance(
                                         intr_list[1].end) else
                                     intr_list[1].end)

            elif intr_list[0].is_close(middle_point) and \
                    polygon_primitive.length() == intr_list[1].length():
                closing_point = intr_list[1].start
                distance = 0

        return closing_point

    def select_closest_sewing_closing_point(self,
                                            line_segment: volmdlr.edges.LineSegment2D,
                                            polygon_primitive,
                                            possible_closing_points):
        """
        Searches the closest sewing closing point available.

        """
        closing_point = volmdlr.O2D
        middle_point = polygon_primitive.middle_point()
        distance = math.inf
        for intr_list in possible_closing_points:
            if intr_list[1] not in self.line_segments:
                dist = intr_list[0].point_distance(line_segment.start)
                if dist < distance:
                    distance = dist
                    closing_point = (intr_list[1].start if
                                     intr_list[0].point_distance(
                                         intr_list[1].start) <
                                     intr_list[0].point_distance(
                                         intr_list[1].end) else
                                     intr_list[1].end)

            elif intr_list[0].is_close(middle_point) and \
                    polygon_primitive.length() == intr_list[1].length():
                closing_point = intr_list[1].start
                distance = 0

        return closing_point

    def search_farthest(self, interseting_point, possible_closing_points):
        """
        Chooses the closest of the farthest available.

        While Sewing two Polygons, and searching a face\'s closing point, this method verifies it
        :return: True if to search the farthest of False if not
        """
        distance = math.inf
        target_prim = None
        for intersection_point, prim in possible_closing_points:
            dist = interseting_point.point_distance(intersection_point)
            if dist < distance:
                distance = dist
                target_prim = prim
        if target_prim in self.line_segments:
            return True
        return False

    def get_closing_point(self, polygon2_2d, primitive, ax=None):
        """Gets sewing closing points for given primitive points."""
        closing_point = volmdlr.O2D
        middle_point = primitive.middle_point()

        normal_vector = primitive.unit_normal_vector()
        line_segment1 = volmdlr.edges.LineSegment2D(middle_point,
                                                    middle_point - normal_vector)
        line_segment2 = volmdlr.edges.LineSegment2D(middle_point,
                                                    middle_point + normal_vector)

        possible_sewing_closing_points_in_linesegment = \
            self.get_possible_sewing_closing_points(polygon2_2d, primitive,
                                                    line_segment1,
                                                    line_segment2)
        if possible_sewing_closing_points_in_linesegment[line_segment1] and \
                not possible_sewing_closing_points_in_linesegment[line_segment2]:
            closing_point = self.select_closest_sewing_closing_point(
                line_segment1, primitive,
                possible_sewing_closing_points_in_linesegment[line_segment1])
            if ax is not None:
                closing_point.plot(ax=ax, color='g')
        if possible_sewing_closing_points_in_linesegment[line_segment2] and \
                not possible_sewing_closing_points_in_linesegment[
                    line_segment1]:
            closing_point = self.select_closest_sewing_closing_point(
                line_segment2, primitive,
                possible_sewing_closing_points_in_linesegment[line_segment2])

        else:
            if len(possible_sewing_closing_points_in_linesegment[line_segment1]) == 1:
                closing_point = self.select_closest_sewing_closing_point(
                    line_segment1, primitive,
                    possible_sewing_closing_points_in_linesegment[
                        line_segment1])
                if closing_point.is_close(volmdlr.O2D):
                    closing_point = self.select_farthest_sewing_closing_point(
                        line_segment2, primitive,
                        possible_sewing_closing_points_in_linesegment[
                            line_segment2])
                if ax is not None:
                    closing_point.plot(ax=ax, color='c')
            elif len(possible_sewing_closing_points_in_linesegment[line_segment2]) == 1:
                closing_point = self.select_closest_sewing_closing_point(
                    line_segment2, primitive,
                    possible_sewing_closing_points_in_linesegment[
                        line_segment2])
                if closing_point.is_close(volmdlr.O2D):
                    closing_point = self.select_farthest_sewing_closing_point(
                        line_segment1, primitive,
                        possible_sewing_closing_points_in_linesegment[
                            line_segment1])
            else:
                if possible_sewing_closing_points_in_linesegment[line_segment1]:
                    if self.search_farthest(
                            middle_point,
                            possible_sewing_closing_points_in_linesegment[
                                line_segment2]):
                        closing_point = \
                            self.select_farthest_sewing_closing_point(
                                line_segment1, primitive,
                                possible_sewing_closing_points_in_linesegment[
                                    line_segment1])
                    else:
                        closing_point = \
                            self.select_closest_sewing_closing_point(
                                line_segment1, primitive,
                                possible_sewing_closing_points_in_linesegment[
                                    line_segment1])

                elif possible_sewing_closing_points_in_linesegment[
                        line_segment2]:
                    closing_point = self.select_closest_sewing_closing_point(
                        line_segment2, primitive,
                        possible_sewing_closing_points_in_linesegment[
                            line_segment2])
        if ax is not None:
            middle_point.plot(ax=ax, color='r')
            line_segment1.plot(ax=ax, edge_style=EdgeStyle(color='y'))
            line_segment2.plot(ax=ax, edge_style=EdgeStyle(color='b'))
            closing_point.plot(ax=ax)
            raise NotImplementedError('There should not be a plot inside this method')

        return closing_point

    def get_valid_sewing_polygon_primitive(self, polygon2_2d):
        """Get valid primitive to start sewing two polygons."""
        for primitive1 in self.line_segments:
            middle_point = primitive1.middle_point()
            normal_vector = primitive1.unit_normal_vector()
            line_segment1 = volmdlr.edges.LineSegment2D(middle_point,
                                                        middle_point - normal_vector)
            line_segment2 = volmdlr.edges.LineSegment2D(middle_point,
                                                        middle_point + normal_vector)
            possible_closing_points = self.get_possible_sewing_closing_points(
                polygon2_2d, primitive1, line_segment1, line_segment2)
            if len(possible_closing_points[line_segment1]) == 1 and \
                    possible_closing_points[line_segment1][0][1] in polygon2_2d.line_segments:
                closing_point = (possible_closing_points[
                                     line_segment1][0][1].start if
                                 possible_closing_points[
                                     line_segment1][0][0].point_distance(
                                     possible_closing_points[
                                         line_segment1][0][1].start) <
                                 possible_closing_points[
                                     line_segment1][0][0].point_distance(
                                     possible_closing_points[
                                         line_segment1][0][1].end) else
                                 possible_closing_points[
                                     line_segment1][0][1].end)

                if polygon2_2d.points.index(closing_point) >= len(polygon2_2d.points) * 2 / 4:
                    return primitive1

            if len(possible_closing_points[line_segment2]) == 1 and \
                    possible_closing_points[line_segment2][0][1] in polygon2_2d.line_segments:
                closing_point = (possible_closing_points[
                                     line_segment2][0][1].start if
                                 possible_closing_points[
                                     line_segment2][0][0].point_distance(
                                     possible_closing_points[
                                         line_segment2][0][1].start) <
                                 possible_closing_points[
                                     line_segment2][0][0].point_distance(
                                     possible_closing_points[
                                         line_segment2][0][1].end) else
                                 possible_closing_points[
                                     line_segment2][0][1].end)

                if polygon2_2d.points.index(closing_point) >= len(polygon2_2d.points) * 2 / 4:
                    return primitive1

        for primitive1 in self.line_segments:
            closing_point = self.get_closing_point(polygon2_2d,
                                                   primitive1)
            if not closing_point.is_close(volmdlr.O2D):
                return primitive1

        raise NotImplementedError('make sure the two polygons '
                                  'you are trying to sew are valid ones')

    def is_convex(self):
        """
        Verifies if a polygon is convex or Not.

        """
        for prim1, prim2 in zip(self.line_segments, self.line_segments[1:] + [self.line_segments[0]]):
            vector1 = prim1.direction_vector()
            vector2 = prim2.direction_vector()
            angle = volmdlr.geometry.clockwise_angle(vector1, vector2)
            if self.is_trigo:
                if angle < math.pi and angle != 0:
                    return False
            elif angle > math.pi and angle != 2 * math.pi:
                return False
        return True

    def axial_symmetry(self, line):
        """
        Finds out the symmetric closed_polygon2d according to a line.

        """

        axial_points = [point.axial_symmetry(line) for point in self.points]

        return self.__class__(points=axial_points)


class Triangle(ClosedPolygonMixin):
    """
    Defines a triangle from 3 points.

    It is a Super Class for Triangle2D and Triangle3D,
    storing their main attribute and methods.


    """

    def __init__(self, point1, point2,
                 point3, name: str = ''):
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
        self.name = name
        self._line_segments = None


class Triangle2D(ClosedPolygon2D):
    """
    Defines a triangle 2D.

    :param point1: triangle point 1.
    :param point2: triangle point 2.
    :param point3: triangle point 3.
    """

    def __init__(self, point1: volmdlr.Point2D, point2: volmdlr.Point2D,
                 point3: volmdlr.Point2D, name: str = ''):
        # TODO: This seems buggy. Is it still used?
        # self.point1 = point1
        # self.point2 = point2
        # self.point3 = point3
        # self.name = name

        ClosedPolygon2D.__init__(self, points=[point1, point2, point3], name=name)
        #
        # Triangle.__init__(self, point1,
        #                   point2,
        #                   point3,
        #                   name)

    def area(self):
        u = self.point2 - self.point1
        v = self.point3 - self.point1
        return abs(u.cross(v)) / 2

    def incircle_radius(self):
        param_a = self.point1.point_distance(self.point2)
        param_b = self.point1.point_distance(self.point3)
        param_c = self.point2.point_distance(self.point3)
        return 2 * self.area() / (param_a + param_b + param_c)

    def circumcircle_radius(self):
        param_a = self.point1.point_distance(self.point2)
        param_b = self.point1.point_distance(self.point3)
        param_c = self.point2.point_distance(self.point3)
        return param_a * param_b * param_c / (self.area() * 4.0)

    def ratio_circumr_length(self):
        return self.circumcircle_radius() / self.length()

    def ratio_incircler_length(self):
        return self.incircle_radius() / self.length()

    def aspect_ratio(self):
        param_a = self.point1.point_distance(self.point2)
        param_b = self.point1.point_distance(self.point3)
        param_c = self.point2.point_distance(self.point3)
        param_s = 0.5 * (param_a + param_b + param_c)
        try:
            return (0.125 * param_a * param_b * param_c / (param_s -
                                                           param_a) / (param_s - param_b) / (param_s - param_c))
        except ZeroDivisionError:
            return 1000000.

    def axial_symmetry(self, line):
        """
        Finds out the symmetric triangle 2d according to a line.

        """

        [point1, point2, point3] = [point.axial_symmetry(line)
                                    for point in [self.point1,
                                                  self.point2,
                                                  self.point3]]

        return self.__class__(point1, point2, point3)


class Circle2D(Contour2D):
    """
    Defines a Circle in two dimensions, with a center and a radius.

    """
    _non_serializable_attributes = ['internal_arcs', 'external_arcs',
                                    'polygon', 'straight_line_contour_polygon',
                                    'primitives', 'basis_primitives']

    def __init__(self, center: volmdlr.Point2D, radius: float, name: str = ''):
        self.center = center
        self.radius = radius
        self.angle = volmdlr.TWO_PI
        self.primitives = self._primitives()

        # self.points = self.tessellation_points()

        Contour2D.__init__(self, self.primitives, name=name)  # !!! this is dangerous

    def __hash__(self):
        return int(round(1e6 * (self.center.x + self.center.y + self.radius)))

    def __eq__(self, other_circle):
        if self.__class__.__name__ != other_circle.__class__.__name__:
            return False

        return math.isclose(self.center.x,
                            other_circle.center.x, abs_tol=1e-06) \
            and math.isclose(self.center.y,
                             other_circle.center.y, abs_tol=1e-06) \
            and math.isclose(self.radius, other_circle.radius,
                             abs_tol=1e-06)

    def _primitives(self):
        points = [
            volmdlr.Point2D(self.center.x + self.radius, self.center.y),
            volmdlr.Point2D(self.center.x, self.center.y - self.radius),
            volmdlr.Point2D(self.center.x - self.radius, self.center.y),
            volmdlr.Point2D(self.center.x, self.center.y + self.radius)]

        return [volmdlr.edges.Arc2D(points[0], points[1], points[2]),
                volmdlr.edges.Arc2D(points[2], points[3], points[0])]

    @classmethod
    def from_arc(cls, arc: volmdlr.edges.Arc2D):
        return cls(arc.center, arc.radius, arc.name + ' to circle')

    def point_belongs(self, point, include_edge_points: bool = False):
        """
        Verifies if a point is inside the Circle 2D.

        :param point: A 2D point to check if it is inside the Circle 2D.
        :type point: `volmdlr.Point2D`
        :param include_edge_points: A Boolean indicating whether points on the edge of the Circle 2D
            should be considered inside the circle.
        :type include_edge_points: bool
        :return: True if point inside the circle or false otherwise.
        :rtype: bool
        """

        if include_edge_points:
            return point.point_distance(self.center) <= self.radius
        return point.point_distance(self.center) < self.radius

    def point_distance(self, point):
        """
        Calculates the distance of given point to the circle.

        :param point: point to calculate distance.
        :return: the distance from the point to the circle 2D.
        """
        return point.point_distance(self.center) - self.radius

    def get_bounding_rectangle(self):

        x_min = self.center.x - self.radius
        x_max = self.center.x + self.radius
        y_min = self.center.y - self.radius
        y_max = self.center.y + self.radius
        return volmdlr.core.BoundingRectangle(x_min, x_max, y_min, y_max)

    def line_intersections(self, line: volmdlr.edges.Line2D, tol=1e-9):
        """
        Calculates the intersections between a circle 2D and Line 2D.

        :param line: line to calculate intersections
        :param tol: tolerance to consider in calculations.
        :return: circle and line intersections.
        """
        full_arc_2d = volmdlr.edges.FullArc2D(
            center=self.center, start_end=self.point_at_abscissa(0),
            name=self.name)
        return full_arc_2d.line_intersections(line, tol)

    def linesegment_intersections(self, linesegment: volmdlr.edges.LineSegment2D, tol=1e-9):
        """
        Calculates the intersections between a circle 2D and line segment 2D.

        :param linesegment: line segment to calculate intersections
        :param tol: tolerance to consider in calculations.
        :return: circle and line segment intersections.
        """
        full_arc_2d = volmdlr.edges.FullArc2D(
            center=self.center, start_end=self.point_at_abscissa(0),
            name=self.name)
        return full_arc_2d.linesegment_intersections(linesegment, tol)

    def cut_by_line(self, line: volmdlr.edges.Line2D):
        intersection_points = self.line_intersections(line)
        if not intersection_points:
            return [self]
        if len(intersection_points) == 1:
            raise NotImplementedError
        if len(intersection_points) == 2:
            linesegment = volmdlr.edges.LineSegment2D(intersection_points[0],
                                                      intersection_points[1])
            arc1, arc2 = self.split(intersection_points[0],
                                    intersection_points[1])
            contour1 = Contour2D([arc1, linesegment.copy()])
            contour2 = Contour2D([arc2, linesegment.copy()])
            return [contour1, contour2]
        raise ValueError

    def circle_intersections(self, circle: 'Circle2D'):
        x0, y0 = self.center
        x1, y1 = circle.center
        # r0 = self.radius
        # r1 = circle.radius

        distance = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)

        # non-intersecting
        if distance > self.radius + circle.radius:
            return []
        # One circle within other
        if distance < abs(self.radius - circle.radius):
            return []
        # coincident circles
        if distance == 0 and self.radius == circle.radius:
            return []
        a_param = (self.radius ** 2 - circle.radius ** 2 + distance ** 2) / (2 * distance)
        h_param = math.sqrt(self.radius ** 2 - a_param ** 2)
        x2 = x0 + a_param * (x1 - x0) / distance
        y2 = y0 + a_param * (y1 - y0) / distance
        x3 = x2 + h_param * (y1 - y0) / distance
        y3 = y2 - h_param * (x1 - x0) / distance

        x4 = x2 - h_param * (y1 - y0) / distance
        y4 = y2 + h_param * (x1 - x0) / distance

        return [volmdlr.Point2D(x3, y3), volmdlr.Point2D(x4, y4)]

    def arc_intersections(self, arc2d: volmdlr.edges.Arc2D):
        circle = Circle2D(arc2d.center, arc2d.radius)
        intersections = []

        for inter in self.circle_intersections(circle):
            try:
                arc2d.abscissa(inter)  # I guess it is a test?
                intersections.append(inter)
            except ValueError:
                pass
        return intersections

    def length(self):
        """
        Calculates the length of the Circle 2D.

        :return: the circle's length.
        """

        return volmdlr.TWO_PI * self.radius

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        return vm_common_operations.plot_circle(self, ax, edge_style)

    def to_3d(self, plane_origin, x, y):
        """
        Transforms a Circle2D into an Circle3D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Circle3D.
        """
        normal = x.cross(y)
        center3d = self.center.to_3d(plane_origin, x, y)
        return Circle3D(volmdlr.Frame3D(center3d, x, y, normal),
                        self.radius, self.name)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Circle2D rotation.

        :param center: rotation center.
        :param angle: angle rotation.
        :return: a new rotated Circle2D.
        """
        return Circle2D(self.center.rotation(center, angle), self.radius)

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        Circle2D rotation. Object is updated in-place.

        :param center: rotation center
        :param angle: rotation angle
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        self.center.rotation_inplace(center, angle)

    def translation(self, offset: volmdlr.Vector2D):
        """
        Circle2D translation.

        :param offset: translation vector
        :return: A new translated Circle2D
        """
        return Circle2D(self.center.translation(offset), self.radius)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Circle2D translation. Object is updated in-place.

        :param offset: translation vector
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        self.center.translation_inplace(offset)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Circle2D.

        side = 'old' or 'new'
        """
        if side == 'old':
            return Circle2D(frame.local_to_global_coordinates(self.center),
                            self.radius)
        if side == 'new':
            return Circle2D(frame.global_to_local_coordinates(self.center),
                            self.radius)
        raise ValueError('Side should be \'new\' \'old\'')

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated in-place.

        side = 'old' or 'new'
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        if side == 'old':
            self.center = frame.local_to_global_coordinates(self.center)
        elif side == 'new':
            self.center = frame.global_to_local_coordinates(self.center)
        else:
            raise ValueError('Side should be \'new\' \'old\'')

    def area(self):
        return math.pi * self.radius ** 2

    def second_moment_area(self, point):
        """Second moment area of part of disk."""
        sma = math.pi * self.radius ** 4 / 4
        return volmdlr.geometry.huygens2d(sma, sma, 0, self.area(), self.center, point)

    def center_of_mass(self):
        return self.center

    def point_symmetric(self, point):
        center = 2 * point - self.center
        return Circle2D(center, self.radius)

    def plot_data(self, edge_style: plot_data.EdgeStyle = None,
                  surface_style: plot_data.SurfaceStyle = None):
        return plot_data.Circle2D(cx=self.center.x, cy=self.center.y,
                                  r=self.radius,
                                  edge_style=edge_style,
                                  surface_style=surface_style)

    def copy(self, *args, **kwargs):
        return Circle2D(self.center.copy(), self.radius)

    def point_at_abscissa(self, curvilinear_abscissa):
        start = self.center + self.radius * volmdlr.X3D
        return start.rotation(self.center,
                              curvilinear_abscissa / self.radius)

    def split_by_line(self, line: volmdlr.edges.Line2D):
        """
        Split the Circle with a line into two Arc2D.
        """
        split_points = self.line_intersections(line)
        return self.split(split_points[0], split_points[1])

    def split(self, split_start, split_end):
        x1, y1 = split_start - self.center
        x2, y2 = split_end - self.center

        angle1 = math.atan2(y1, x1)
        angle2 = math.atan2(y2, x2)
        angle_i1 = 0.5 * (angle2 - angle1)
        angle_i2 = angle_i1 + math.pi
        interior_point1 = split_start.rotation(self.center, angle_i1)
        interior_point2 = split_start.rotation(self.center, angle_i2)

        return [volmdlr.edges.Arc2D(split_start, interior_point1,
                                    split_end),
                volmdlr.edges.Arc2D(split_start, interior_point2,
                                    split_end)]

    def axial_symmetry(self, line):
        """
        Finds out the symmetric circle 2d according to a line.
        """
        return self.__class__(center=self.center.axial_symmetry(line),
                              radius=self.radius)

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 40):
        """
        Discretize a Contour to have "n" points.

        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Useful
            to mesh an arc
        :return: a list of sampled points
        """
        if not number_points and angle_resolution:
            number_points = math.ceil(math.pi * angle_resolution) + 2
        step = self.length() / (number_points - 1)
        return [self.point_at_abscissa(i * step) for i in range(number_points)]

    def polygon_points(self, discretization_resolution: int):
        warnings.warn('polygon_points is deprecated,\
        please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(angle_resolution=discretization_resolution)

    def get_geo_points(self):
        return [volmdlr.Point3D(self.radius, self.center.y, 0),
                volmdlr.Point3D(self.center.x, self.center.y, 0),
                volmdlr.Point3D(-self.radius, self.center.y, 0)]

    def bsplinecurve_intersections(self, bsplinecurve: volmdlr.edges.BSplineCurve2D, abs_tol: float = 1e-7):
        """
        Calculates the intersections between a circle 2d and a BSpline Curve 2D.

        :param bsplinecurve: bsplinecurve to search for intersections.
        :param abs_tol: tolerance to be considered while validating an intersection.
        :return: a list with all intersections between circle and bsplinecurve.
        """
        circle_bounding_rectangle = self.bounding_rectangle
        bspline_discretized_points = bsplinecurve.discretization_points(number_points=10)
        param_intersections = []
        for point1, point2 in zip(bspline_discretized_points[:-1], bspline_discretized_points[1:]):
            line_seg = volmdlr.edges.LineSegment2D(point1, point2)
            abscissa1 = bsplinecurve.abscissa(point1)
            abscissa2 = bsplinecurve.abscissa(point2)
            if line_seg.bounding_rectangle.b_rectangle_intersection(circle_bounding_rectangle):
                intersection = self.linesegment_intersections(line_seg)
                if intersection:
                    param_intersections.append((abscissa1, abscissa2))
        intersections = []
        while True:
            if not param_intersections:
                break
            for abscissa1, abscissa2 in param_intersections:
                discretized_points_between_1_2 = [bsplinecurve.point_at_abscissa(abscissa) for abscissa
                                                  in npy.linspace(abscissa1, abscissa2, num=10)]
                break_flag = False
                for point1, point2 in zip(discretized_points_between_1_2[:-1], discretized_points_between_1_2[1:]):
                    line_seg = volmdlr.edges.LineSegment2D(point1, point2)
                    if line_seg.bounding_rectangle.b_rectangle_intersection(circle_bounding_rectangle):
                        intersection = self.linesegment_intersections(line_seg, 1e-12)
                        if not intersection:
                            continue
                        if bsplinecurve.point_distance(intersection[0]) > abs_tol:
                            param_intersections.insert(0, (bsplinecurve.abscissa(point1),
                                                           bsplinecurve.abscissa(point2)))
                        else:
                            intersections.append(intersection[0])
                        break_flag = True
                        break
                if break_flag:
                    break
            else:
                continue
            param_intersections.remove((abscissa1, abscissa2))
        return intersections


class Ellipse2D(Contour2D):
    """
    Defines an Ellipse in two-dimensions.

    Ellipse2D defined by a major axis (A), minor axis (B), a center and a vector
    representing the direction of the major axis.

    :param major_axis: ellipse's major axis (A)
    :type major_axis: float
    :param minor_axis: ellipse's minor axis (B)
    :type minor_axis: float
    :param center: ellipse's center
    :type center: volmdlr.Point3D
    :param major_dir: direction vector for major axis
    :type major_dir: volmdlr.Vector3D

    :Example:
    >>> ellipse2d = Ellipse2D(4, 2, volmdlr.O2D, volmdlr.Vector2D(1, 1))
    """

    def __init__(self, major_axis, minor_axis, center, major_dir, name=''):
        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.center = center
        self.major_dir = major_dir
        self.major_dir.normalize()
        self.minor_dir = - self.major_dir.normal_vector()
        self.theta = volmdlr.geometry.clockwise_angle(self.major_dir, volmdlr.X2D)
        if self.theta == math.pi * 2:
            self.theta = 0.0
        Contour2D.__init__(self, [self], name=name)

    def __hash__(self):
        return int(round(1e6 * (self.center.x + self.center.y + self.major_axis + self.minor_axis)))

    def area(self):
        """
        Calculates the ellipse's area.

        :return: ellipse's area, float.
        """
        return math.pi * self.major_axis * self.minor_axis

    def length(self):
        """
        Calculates the ellipse's length.

        :return: ellipse's length.
        """
        mid_point = self.center - self.major_axis * self.major_dir
        if self.theta != 0.0:
            mid_point = self.center - volmdlr.Point2D(self.major_axis, 0)
            mid_point = mid_point.rotation(self.center, self.theta)
        length = 2 * self.abscissa(mid_point)
        return length

    def to_3d(self, plane_origin, x, y):
        """
        Transforms a Ellipse2D into an Ellipse3D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Ellipse3D.
        """
        raise NotImplementedError

    def point_over_ellipse(self, point, abs_tol=1e-6):
        """
        Verifies if a point is on the ellipse.

        :param point: point to be verified.
         :param abs_tol: tolerance.
        :return: True or False.
        """
        return math.isclose(
            ((point.x - self.center.x) * math.cos(self.theta) +
             (point.y - self.center.y) * math.sin(self.theta)) ** 2 / self.major_axis ** 2 +
            ((point.x - self.center.x) * math.sin(self.theta) -
             (point.y - self.center.y) * math.cos(self.theta)) ** 2 / self.minor_axis ** 2, 1, abs_tol=abs_tol)

    def point_over_contour(self, point, abs_tol=1e-6):
        """
        Verifies if a point is on the ellipse.

        :param point: point to be verified.
        :param abs_tol: tolerance.
        :return: True or False.
        """
        return self.point_over_ellipse(point, abs_tol)

    def line_intersections(self, line: 'volmdlr.edges.Line2D'):
        """
        Calculates the intersections between a line and an ellipse.

        :param line: line to calculate intersections
        :return: list of points intersections, if there are any
        """
        intersections = vm_utils_intersections.ellipse2d_line_intersections(self, line)
        return intersections

    def linesegment_intersections(self, linesegment: 'volmdlr.edges.LineSegment2D'):
        """
        Calculates the intersections between a line segment and an ellipse.

        :param linesegment: line segment to calculate intersections.
        :return: list of points intersections, if there are any.
        """
        line_intersections = self.line_intersections(linesegment.to_line())
        intersections = []
        for intersection in line_intersections:
            if linesegment.point_belongs(intersection):
                intersections.append(intersection)
        return intersections

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """
        Calculates the discretized points for the ellipse.

        :param number_points: number of point to have in the discretized points.
        :param angle_resolution: the angle resolution to be used to discretize points.
        :return: discretized points.
        """
        if number_points:
            angle_resolution = number_points
        discretization_points = [self.center + volmdlr.Point2D(self.major_axis * math.cos(theta),
                                                               self.minor_axis * math.sin(theta))
                                 for theta in npy.linspace(0, volmdlr.TWO_PI, angle_resolution + 1)]
        discretization_points = [point.rotation(self.center, self.theta) for point in discretization_points]
        return discretization_points

    def abscissa(self, point: volmdlr.Point2D):
        """
        Calculates the abscissa for a given point.

        :param point: point to calculate the abscissa.
        :return: the corresponding abscissa, 0 < abscissa < ellipse's length.
        """
        if self.point_over_ellipse(point):
            angle_abscissa = self.point_angle_with_major_dir(point)

            def arc_length(theta):
                return math.sqrt((self.major_axis ** 2) * math.sin(theta) ** 2 +
                                 (self.minor_axis ** 2) * math.cos(theta) ** 2)

            res, _ = scipy_integrate.quad(arc_length, 0, angle_abscissa)
            return res
        raise ValueError(f'point {point} does not belong to ellipse')

    def point_angle_with_major_dir(self, point2d):
        """
        Given a point in the ellipse, calculates it angle with the major direction vector.

        """
        center2d_point2d = point2d - self.center
        angle_abscissa = volmdlr.geometry.clockwise_angle(center2d_point2d, self.major_dir)
        return angle_abscissa

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Matplotlib plot for an ellipse.

        """
        if ax is None:
            _, ax = plt.subplots()
        x = []
        y = []
        for point_x, point_y in self.discretization_points(number_points=50):
            x.append(point_x)
            y.append(point_y)
        plt.plot(x, y, color=edge_style.color, alpha=edge_style.alpha)
        if edge_style.equal_aspect:
            ax.set_aspect('equal')
        return ax

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Rotation of ellipse around a center and an angle.

        :param center: center of the rotation.
        :param angle: angle to rotated of.
        :return: a rotationed new ellipse.
        """
        rotationed_center = self.center.rotation(center, angle)
        point_major_dir = self.center + self.major_dir * self.major_axis
        rotationed_major_dir_point = point_major_dir.rotation(center, angle)
        major_dir = rotationed_major_dir_point - rotationed_center
        return Ellipse2D(self.major_axis, self.minor_axis, rotationed_center,
                         major_dir)

    def translation(self, offset: volmdlr.Vector2D):
        """
        Translation of ellipse from an offset vector.

        :param offset: corresponding translation vector.
        :return: translated new ellipse 2d.
        """
        return Ellipse2D(self.major_axis, self.minor_axis, self.center.translation(offset), self.major_dir)

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and return a new Ellipse2D.

        side = 'old' or 'new'.
        """
        if side == 'old':
            return Ellipse2D(self.major_axis, self.minor_axis, frame.local_to_global_coordinates(self.center),
                             self.major_dir)
        if side == 'new':
            point_major_dir = self.center + self.major_dir * self.major_axis
            major_dir = frame.global_to_local_coordinates(point_major_dir) - self.center
            return Ellipse2D(self.major_axis, self.minor_axis, frame.global_to_local_coordinates(self.center),
                             major_dir)
        raise ValueError('Side should be \'new\' \'old\'')


class Contour3D(ContourMixin, Wire3D):
    """
    A collection of 3D primitives forming a closed wire3D.

    """
    _non_serializable_attributes = ['points']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['points', 'name']
    _generic_eq = True

    def __init__(self, primitives: List[volmdlr.core.Primitive3D],
                 name: str = ''):
        """
        Defines a contour3D from a collection of edges following each other stored in primitives list.
        """

        Wire3D.__init__(self, primitives=primitives, name=name)
        self._edge_polygon = None
        self._utd_bounding_box = False

    def __hash__(self):
        return hash(tuple(self.primitives))

    def __eq__(self, other_):
        if other_.__class__.__name__ != self.__class__.__name__:
            return False
        if len(self.primitives) != len(other_.primitives):
            return False
        equal = 0
        for prim1 in self.primitives:
            reverse1 = prim1.reverse()
            found = False
            for prim2 in other_.primitives:
                reverse2 = prim2.reverse()
                if (prim1 == prim2 or reverse1 == prim2
                        or reverse2 == prim1 or reverse1 == reverse2):
                    equal += 1
                    found = True
            if not found:
                return False
        if equal == len(self.primitives):
            return True
        return False

    @property
    def edge_polygon(self):
        if self._edge_polygon is None:
            self._edge_polygon = self._get_edge_polygon()
        return self._edge_polygon

    def _get_edge_polygon(self):
        points = []
        for edge in self.primitives:
            if points:
                if not edge.start.is_close(points[-1]):
                    points.append(edge.start)
            else:
                points.append(edge.start)
        return ClosedPolygon3D(points)

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a Contour3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding Contour3D object.
        :rtype: :class:`volmdlr.wires.Contour3D`
        """
        step_id = kwargs.get("step_id", "#UNKNOW_ID")
        step_name = kwargs.get("name", "EDGE_LOOP")
        name = arguments[0][1:-1]
        raw_edges = []
        for edge_id in arguments[1]:
            edge = object_dict[int(edge_id[1:])]
            if edge:
                raw_edges.append(edge)

        if step_name == "POLY_LOOP":
            return cls.from_points(raw_edges)
        if (len(raw_edges)) == 1:
            if isinstance(raw_edges[0], cls):
                # Case of a circle, ellipse...
                return raw_edges[0]
            return cls(raw_edges, name=name)

        # if any(edge is None for edge in raw_edges):
        #     raw_edges = [edge for edge in raw_edges if edge is not None]
            # warnings.warn(f"Could not instantiate #{step_id} = {step_name}({arguments})"
            #               f" because some of the edges are NoneType."
            #               "See Contour3D.from_step method")
            # return None
        # Making things right for first 2 primitives
        distances = [raw_edges[0].end.point_distance(raw_edges[1].start),
                     raw_edges[0].start.point_distance(raw_edges[1].start),
                     raw_edges[0].end.point_distance(raw_edges[1].end),
                     raw_edges[0].start.point_distance(raw_edges[1].end)]
        index = distances.index(min(distances))
        if min(distances) > 6e-4:
            # Green color : well-placed and well-read
            ax = raw_edges[0].plot(edge_style=EdgeStyle(color='g'))
            ax.set_title(f"Step ID: #{step_id}")

            # Red color : can't be connected to green edge
            raw_edges[1].plot(ax=ax, edge_style=EdgeStyle(color='r'))
            # Black color : to be placed
            for re in raw_edges[2:]:
                re.plot(ax=ax)

            warnings.warn(
                f"Could not instantiate #{step_id} = {step_name}({arguments})"
                "because the first 2 edges of contour not following each other.\n"
                f'Number of edges: {len(raw_edges)}.\n'
                f'delta_x = {abs(raw_edges[0].start.x - raw_edges[1].end.x)}, '
                f' {abs(raw_edges[0].end.x - raw_edges[1].end.x)}.\n'
                f'delta_y = {abs(raw_edges[0].start.y - raw_edges[1].end.y)} ,'
                f' {abs(raw_edges[0].end.y - raw_edges[1].end.y)}.\n'
                f'delta_z = {abs(raw_edges[0].start.z - raw_edges[1].end.z)}, '
                f' {abs(raw_edges[0].end.z - raw_edges[1].end.z)}.\n'
                f'distance = {min(distances)}')
            return None

        if index == 0:
            edges = [raw_edges[0], raw_edges[1]]
        elif index == 1:
            edges = [raw_edges[0].reverse(), raw_edges[1]]
        elif index == 2:
            edges = [raw_edges[0], raw_edges[1].reverse()]
        elif index == 3:
            edges = [raw_edges[0].reverse(), raw_edges[1].reverse()]
        else:
            raise NotImplementedError

        # Connecting the next edges
        last_edge = edges[-1]
        for i, raw_edge in enumerate(raw_edges[2:]):
            distances = [raw_edge.start.point_distance(last_edge.end),
                         raw_edge.end.point_distance(last_edge.end)]
            index = distances.index(min(distances))
            if min(distances) > 6e-4:
                # Green color : well-placed and well-read
                ax = last_edge.plot(edge_style=EdgeStyle(color='g'))
                ax.set_title(f"Step ID: #{step_id}")

                for re in raw_edges[:2 + i]:
                    re.plot(ax=ax, edge_style=EdgeStyle(color='g'))
                    re.start.plot(ax=ax, color='g')
                    re.end.plot(ax=ax, color='g')
                last_edge.end.plot(ax=ax, color='g')
                # Red color : can't be connected to red dot
                raw_edge.plot(ax=ax, edge_style=EdgeStyle(color='g'))
                # Black color : to be placed
                for re in raw_edges[2 + i + 1:]:
                    re.plot(ax=ax)
                    re.start.plot(ax=ax)
                    re.end.plot(ax=ax)

                warnings.warn(
                    f"Could not instantiate #{step_id} = {step_name}({arguments})"
                    "because some Edges of contour are not following each other.\n"
                    f'Number of edges: {len(raw_edges)}.\n'
                    f'delta_x = {abs(raw_edge.start.x - last_edge.end.x)}, '
                    f' {abs(raw_edge.end.x - last_edge.end.x)}.\n'
                    f'delta_y = {abs(raw_edge.start.y - last_edge.end.y)}, '
                    f' {abs(raw_edge.end.y - last_edge.end.y)}.\n'
                    f'delta_z = {abs(raw_edge.start.z - last_edge.end.z)}, '
                    f' {abs(raw_edge.end.z - last_edge.end.z)}.\n'
                    f'distance = {min(distances)}')
                return None
            if index == 0:
                last_edge = raw_edge
            elif index == 1:
                last_edge = raw_edge.reverse()

            edges.append(last_edge)
        return cls(edges, name=name)

    def to_step(self, current_id, surface_id=None, surface3d=None):
        """
        Create a Circle3D step object.

        """
        content = ''
        edge_ids = []
        for primitive in self.primitives:
            if primitive.__class__.__name__ == 'BSplineCurve3D':
                method_name = f'{primitive.__class__.__name__.lower()}_to_2d'
                curve2d = getattr(surface3d, method_name)(primitive)[0]
                if curve2d.__class__.__name__ == 'LineSegment3D':
                    curve2d = curve2d.to_bspline_curve()
                primitive_content, primitive_ids = primitive.to_step(
                    current_id, surface_id=surface_id, curve2d=curve2d)
            else:
                primitive_content, primitive_ids = primitive.to_step(current_id, surface_id=surface_id)

            content += primitive_content
            current_id = primitive_ids[-1] + 1
            for primitive_id in primitive_ids:
                content += f"#{current_id} = ORIENTED_EDGE('{primitive.name}',*,*,#{primitive_id},.T.);\n"
                edge_ids.append(current_id)

                current_id += 1

        content += f"#{current_id} = EDGE_LOOP('{self.name}',({volmdlr.core.step_ids_to_str(edge_ids)}));\n"
        return content, current_id

    def average_center_point(self):
        nb = len(self.edge_polygon.points)
        x = sum(point[0] for point in self.edge_polygon.points) / nb
        y = sum(point[1] for point in self.edge_polygon.points) / nb
        z = sum(point[2] for point in self.edge_polygon.points) / nb

        return volmdlr.Point3D(x, y, z)

    def to_2d(self, plane_origin, x, y):
        primitives2d = self.get_primitives_2d(plane_origin, x, y)
        return Contour2D(primitives=primitives2d)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Contour3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated Contour3D.
        """
        new_edges = [edge.rotation(center, axis, angle) for edge
                     in self.primitives]
        return Contour3D(new_edges, self.name)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        Contour3D rotation. Object is updated in-place.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: rotation angle.
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for edge in self.primitives:
            edge.rotation_inplace(center, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Contour3D translation.

        :param offset: translation vector.
        :return: A new translated Contour3D.
        """
        new_edges = [edge.translation(offset) for edge in
                     self.primitives]
        return Contour3D(new_edges, self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Contour3D translation. Object is updated in-place.

        :param offset: translation vector.
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for edge in self.primitives:
            edge.translation_inplace(offset)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Contour3D.

        side = 'old' or 'new'.
        """
        new_edges = [edge.frame_mapping(frame, side) for edge in
                     self.primitives]
        return Contour3D(new_edges, self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated in-place.

        :param side: 'old' or 'new'
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for edge in self.primitives:
            edge.frame_mapping_inplace(frame, side)

    def copy(self, deep=True, memo=None):
        """
        Copies the Contour3D.
        """
        new_edges = [edge.copy(deep=deep, memo=memo) for edge in self.primitives]
        # if self.point_inside_contour is not None:
        #     new_point_inside_contour = self.point_inside_contour.copy()
        # else:
        #     new_point_inside_contour = None
        return Contour3D(new_edges, self.name)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            # ax = Axes3D(plt.figure())
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        for edge in self.primitives:
            edge.plot(ax=ax, edge_style=edge_style)

        return ax

    def _bounding_box(self):
        """
        Computes the bounding box of the contour3D.

        """
        return volmdlr.core.BoundingBox.from_bounding_boxes([prim.bounding_box for prim in self.primitives])

    @property
    def bounding_box(self):
        if not self._utd_bounding_box:
            self._bbox = self._bounding_box()
            self._utd_bounding_box = True
        return self._bbox

    def line_intersections(self, line: volmdlr.edges.Line3D):
        """
        Calculates intersections between a contour 3d and Line 3d.

        :param line: Line 3D to verify intersections.
        :return: list with the contour intersections with line
        """
        intersections = []
        for primitive in self.primitives:
            prim_line_intersections = primitive.line_intersections(line)
            if prim_line_intersections:
                for inters in prim_line_intersections:
                    if inters not in intersections:
                        intersections.append(inters)
        return intersections

    def linesegment_intersections(self, linesegment: volmdlr.edges.LineSegment3D):
        """
        Calculates intersections between a contour 3d and line segment 3D.

        :param linesegment: line segment 3D to verify intersections.
        :return: list with the contour intersections with line
        """
        intersections = []
        for primitive in self.primitives:
            prim_line_intersections = primitive.linesegment_intersections(linesegment)
            if prim_line_intersections:
                for inters in prim_line_intersections:
                    if inters not in intersections:
                        intersections.append(inters)
        return intersections

    def contour_intersection(self, contour3d):
        """
        Calculates intersections between two Contour3D.

        :param contour3d: second contour
        :return: list of points
        """
        dict_intersecting_points = {}
        for primitive in self.primitives:
            for primitive2 in contour3d.primitives:
                intersecting_point = primitive.linesegment_intersection(
                    primitive2)
                if intersecting_point is not None:
                    dict_intersecting_points[primitive2] = intersecting_point
        if dict_intersecting_points:
            return dict_intersecting_points
        return None

    @classmethod
    def from_points(cls, points: List[volmdlr.Point3D]):
        """
        Create a contour 3d from points with line_segments3D.

        """

        if len(points) < 3:
            raise ValueError('contour is defined at least with three points')
        edges = []
        for i in range(0, len(points) - 1):
            edges.append(volmdlr.edges.LineSegment3D(points[i], points[i + 1]))

        edges.append(volmdlr.edges.LineSegment3D(points[-1], points[0]))
        contour = cls(edges)

        return contour

    def clean_primitives(self):
        """
        Delete primitives with start=end, and return a new contour.

        """

        new_primitives = []
        for primitive in self.primitives:
            if not primitive.start.is_close(primitive.end):
                new_primitives.append(primitive)

        return Contour3D(new_primitives)

    def merge_with(self, contour3d):
        """
        Merge two adjacent contours, and returns one outer contour and inner contours (if there are any).

        """

        merged_primitives = self.delete_shared_contour_section(contour3d)
        contours = Contour3D.contours_from_edges(merged_primitives, tol=1e-6)

        return contours


class Circle3D(Contour3D):
    """
    Defines a Circle in three dimensions, with a center and a radius.

    frame.u, frame.v define the plane, frame.w the normal
    """
    _non_serializable_attributes = ['point', 'edges', 'point_inside_contour']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, frame: volmdlr.Frame3D, radius: float,
                 name: str = ''):
        self.radius = radius
        self.frame = frame
        self.angle = volmdlr.TWO_PI
        self._primitives = None
        self.primitives = self.get_primitives()
        Contour3D.__init__(self, self.primitives, name=name)

    @property
    def center(self):
        return self.frame.origin

    @property
    def normal(self):
        return self.frame.w

    def __hash__(self):
        return hash(self.frame.origin)

    def __eq__(self, other_circle):
        return self.frame.origin.is_close(other_circle.frame.origin) \
            and self.frame.w.is_colinear(other_circle.frame.w) \
            and math.isclose(self.radius,
                             other_circle.radius, abs_tol=1e-06)

    def get_primitives(self):
        """
        Calculates primitives to compose Circle: 2 Arc3D.

        :return: list containing two Arc3D
        """
        if not self._primitives:
            points = [self.center + self.frame.u * self.radius,
                      self.center - self.frame.v * self.radius,
                      self.center - self.frame.u * self.radius,
                      self.center + self.frame.v * self.radius]
            self._primitives = [volmdlr.edges.Arc3D(points[0], points[1], points[2]),
                                volmdlr.edges.Arc3D(points[2], points[3], points[0])]

        return self._primitives

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """
        Discretize a Circle to have "n" points.

        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Useful
            to mesh an arc
        :return: a list of sampled points
        """
        if number_points:
            angle_resolution = number_points
        discretization_points_3d = [self.center + self.radius * math.cos(teta) * self.frame.u +
                                    self.radius * math.sin(teta) * self.frame.v for teta in
                                    npy.linspace(0, volmdlr.TWO_PI, angle_resolution + 1)][:-1]
        return discretization_points_3d

    def abscissa(self, point: volmdlr.Point3D):
        """
        Calculates the abscissa a given point.

        :param point: point to calculate abscissa.
        :return: abscissa
        """
        x, y, _ = self.frame.global_to_local_coordinates(point)
        u1 = x / self.radius
        u2 = y / self.radius
        theta = volmdlr.geometry.sin_cos_angle(u1, u2)

        return self.radius * abs(theta)

    def length(self):
        return volmdlr.TWO_PI * self.radius

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Circle3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Circle3D
        """
        return Circle3D(self.frame.rotation(center, axis, angle),
                        self.radius, self.name)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Circle3D rotation. Object is updated in-place.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        self.frame.rotation_inplace(center, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Circle3D translation.

        :param offset: translation vector
        :return: A new translated Circle3D
        """
        return Circle3D(self.frame.translation(offset), self.radius, self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Circle3D translation. Object is updated in-place.

        :param offset: translation vector
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        self.frame.translation_inplace(offset)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = None

        x = []
        y = []
        z = []
        for point_x, point_y, point_z in self.discretization_points():
            x.append(point_x)
            y.append(point_y)
            z.append(point_z)
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])
        ax.plot(x, y, z, color=edge_style.color, alpha=edge_style.alpha)
        return ax

    def point_at_abscissa(self, curvilinear_abscissa):
        """ Start point is at intersection of frame.u axis. """
        start = self.frame.origin + self.radius * self.frame.u
        return start.rotation(self.frame.origin, self.frame.w,
                              curvilinear_abscissa / self.radius)

    def linesegment_intersections(self, linesegment: volmdlr.edges.LineSegment3D):
        """
        Calculates the intersections between the Circle3D and a line segment 3D.

        :param linesegment: line segment 3D to verify intersections
        :return: list of points intersecting Circle
        """
        intersections = []
        circle3d_line_intersections = vm_utils_intersections.circle_3d_line_intersections(self, linesegment.to_line())
        for intersection in circle3d_line_intersections:
            if linesegment.point_belongs(intersection):
                intersections.append(intersection)
        return intersections

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a Circle3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding Circle3D object.
        :rtype: :class:`volmdlr.wires.Circle3D`
        """
        length_conversion_factor = kwargs.get("length_conversion_factor", 1)

        center = object_dict[arguments[1]].origin
        radius = float(arguments[2]) * length_conversion_factor
        if object_dict[arguments[1]].u is not None:
            normal = object_dict[arguments[1]].u
            other_vec = object_dict[arguments[1]].v
            if other_vec is not None:
                other_vec.normalize()
        else:
            normal = object_dict[arguments[1]].v  # ou w
            other_vec = None
        normal.normalize()
        return cls.from_center_normal(center, normal, radius, arguments[0][1:-1])

    def to_step(self, current_id, surface_id=None, surface3d=None):
        circle_frame = volmdlr.Frame3D(self.center, self.frame.w, self.frame.u,
                                       self.frame.v)
        content, frame_id = circle_frame.to_step(current_id)
        curve_id = frame_id + 1
        content += f"#{curve_id} = CIRCLE('{self.name}',#{frame_id},{round(self.radius * 1000, 3)});\n"

        if surface_id:
            content += f"#{curve_id + 1} = SURFACE_CURVE('',#{curve_id},(#{surface_id}),.PCURVE_S1.);\n"
            curve_id += 1

        point1 = self.frame.origin + self.frame.u * self.radius
        point3 = self.frame.origin - self.frame.u * self.radius

        p1_content, p1_id = point1.to_step(curve_id + 1, vertex=True)
        p3_content, p3_id = point3.to_step(p1_id + 1, vertex=True)
        content += p1_content + p3_content

        arc1_id = p3_id + 1
        content += f"#{arc1_id} = EDGE_CURVE('{self.name}',#{p1_id},#{p3_id},#{curve_id},.T.);\n"
        oriented_edge1_id = arc1_id + 1
        content += f"#{oriented_edge1_id} = ORIENTED_EDGE('',*,*,#{arc1_id},.T.);\n"

        arc2_id = oriented_edge1_id + 1
        content += f"#{arc2_id} = EDGE_CURVE('{self.name}',#{p3_id},#{p1_id},#{curve_id},.T.);\n"
        oriented_edge2_id = arc2_id + 1
        content += f"#{oriented_edge2_id} = ORIENTED_EDGE('',*,*,#{arc2_id},.T.);\n"

        current_id = oriented_edge2_id + 1
        content += f"#{current_id} = EDGE_LOOP('{self.name}',(#{oriented_edge1_id},#{oriented_edge2_id}));\n"

        return content, current_id

    def _bounding_box(self):
        """
        Computes the bounding box.

        """
        points = [self.frame.origin + self.radius * v
                  for v in [self.frame.u, -self.frame.u,
                            self.frame.v, -self.frame.v]]
        return volmdlr.core.BoundingBox.from_points(points)

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a Circle3D into an Circle2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Circle2D.
        """
        z = x.cross(y)
        plane3d = volmdlr.faces.Plane3D(volmdlr.Frame3D(plane_origin, x, y, z))
        return Circle2D(plane3d.point3d_to_2d(self.center), self.radius)

    @classmethod
    def from_center_normal(cls, center: volmdlr.Point3D,
                           normal: volmdlr.Vector3D,
                           radius: float,
                           name: str = ''):
        u = normal.deterministic_unit_normal_vector()
        v = normal.cross(u)
        return cls(volmdlr.Frame3D(center, u, v, normal), radius, name)

    @classmethod
    def from_3_points(cls, point1, point2, point3):
        vector_u1 = point2 - point1
        vector_u2 = point2 - point3
        try:
            vector_u1.normalize()
            vector_u2.normalize()
        except ZeroDivisionError:
            raise ValueError(
                'the 3 points must be distincts')

        normal = vector_u2.cross(vector_u1)
        normal.normalize()

        if vector_u1.is_close(vector_u2):
            vector_u2 = normal.cross(vector_u1)
            vector_u2.normalize()

        vector_v1 = normal.cross(vector_u1)  # v1 is normal, equal u2
        vector_v2 = normal.cross(vector_u2)  # equal -u1

        point11 = 0.5 * (point1 + point2)  # Mid-point of segment s,m
        point21 = 0.5 * (point2 + point3)  # Mid-point of segment s,m

        line1 = volmdlr.edges.Line3D(point11, point11 + vector_v1)
        line2 = volmdlr.edges.Line3D(point21, point21 + vector_v2)

        try:
            center, _ = line1.minimum_distance_points(line2)
        except ZeroDivisionError:
            raise ValueError(
                'Start, end and interior points  of an arc must be distincts')

        radius = (center - point1).norm()
        return cls(frame=volmdlr.Frame3D(center, vector_u1, normal.cross(vector_u1), normal),
                   radius=radius)

    def extrusion(self, extrusion_vector):
        """
        Returns the cylindrical face generated by extrusion of the circle.
        """
        if self.normal.is_colinear_to(extrusion_vector):
            u = self.normal.deterministic_unit_normal_vector()
            v = self.normal.cross(u)
            w = extrusion_vector.copy()
            w.normalize()
            cylinder = volmdlr.faces.CylindricalSurface3D(
                volmdlr.Frame3D(self.center, u, v, w), self.radius)
            return [cylinder.rectangular_cut(0, volmdlr.TWO_PI,
                                             0, extrusion_vector.norm())]
        raise NotImplementedError(
            f'Extrusion along vector not colinar to normal for circle not '
            f'handled yet: dot={self.normal.dot(extrusion_vector)}')

    def revolution(self, axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                   angle: float):
        """
        Return the Toroidal face generated by the revolution of the circle.
        """
        line3d = volmdlr.edges.Line3D(axis_point, axis_point + axis)
        tore_center, _ = line3d.point_projection(self.center)
        u = self.center - tore_center
        u.normalize()
        v = axis.cross(u)
        if not math.isclose(self.normal.dot(u), 0., abs_tol=1e-9):
            raise NotImplementedError(
                'Outside of plane revolution not supported')

        tore_radius = tore_center.point_distance(self.center)
        surface = volmdlr.faces.ToroidalSurface3D(
            volmdlr.Frame3D(tore_center, u, v, axis),
            tore_radius, self.radius)
        return [surface.rectangular_cut(0, angle, 0, volmdlr.TWO_PI)]

    def point_belongs(self, point: volmdlr.Point3D, abs_tol: float = 1e-6):
        """
        Returns if given point belongs to the Circle3D.
        """
        distance = point.point_distance(self.center)
        vec = volmdlr.Vector3D(*point - self.center)
        dot = self.normal.dot(vec)
        if math.isclose(distance, self.radius, abs_tol=abs_tol) \
                and math.isclose(dot, 0, abs_tol=abs_tol):
            return True
        return False

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D):
        if not self.point_belongs(point1, 1e-4) or not self.point_belongs(point2, 1e-4):
            ax = self.plot()
            point1.plot(ax=ax, color='r')
            point2.plot(ax=ax, color='b')
            raise ValueError('Point not on circle for trim method')

        if point1.is_close(point2):
            return volmdlr.edges.FullArc3D(self.frame.origin, point1, self.frame.w)

        interior = volmdlr.geometry.clockwise_interior_from_circle3d(
            point1, point2, self)
        return volmdlr.edges.Arc3D(point1, interior, point2)


class Ellipse3D(Contour3D):
    """
    Defines a 3D ellipse.

    :param major_axis: Largest radius of the ellipse
    :type major_axis: float
    :param minor_axis: The Smallest radius of the ellipse
    :type minor_axis: float
    :param center: Ellipse's center
    :type center: Point3D
    :param normal: Ellipse's normal
    :type normal: Vector3D
    :param major_dir: Direction of the largest radius/major_axis
    :type major_dir: Vector3D
    """

    def __init__(self, major_axis: float, minor_axis: float,
                 center: volmdlr.Point3D, normal: volmdlr.Vector3D,
                 major_dir: volmdlr.Vector3D, name: str = ''):

        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.center = center
        normal.normalize()
        self.normal = normal
        major_dir.normalize()
        self.major_dir = major_dir
        self.minor_dir = normal.cross(normal)
        self._frame = None
        Contour3D.__init__(self, [self], name=name)

    @property
    def frame(self):
        """
        Gets the Ellipse's Frame3D.

        :return: Frame3D.
        """
        if not self._frame:
            self._frame = volmdlr.Frame3D(self.center, self.major_dir, self.normal.cross(self.major_dir), self.normal)
        return self._frame

    def point_belongs(self, point):
        """
        Verifies if a given point lies on the Ellipse3D.

        :param point: point to be verified.
        :return: True is point lies on the Ellipse, False otherwise
        """
        new_point = self.frame.global_to_local_coordinates(point)
        return math.isclose(new_point.x ** 2 / self.major_axis ** 2 +
                            new_point.y ** 2 / self.minor_axis ** 2, 1.0, abs_tol=1e-6)

    def length(self):
        """
        Calculates the length of the ellipse.

        Ramanujan's approximation for the perimeter of the ellipse.
        P = π (a + b) [ 1 + (3h) / (10 + √(4 - 3h) ) ], where h = (a - b)**2/(a + b)**2
        :return:
        """
        perimeter_formular_h = (self.major_axis - self.minor_axis) ** 2 / (self.major_axis + self.minor_axis) ** 2
        return math.pi * (self.major_axis + self.minor_axis) * \
            (1 + (3 * perimeter_formular_h / (10 + math.sqrt(4 - 3 * perimeter_formular_h))))

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """
        Discretize a Contour to have "n" points.

        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned.
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Useful
            to mesh an arc.
        :return: a list of sampled points.
        """
        if number_points:
            angle_resolution = number_points
        discretization_points_3d = [
                                      self.center + self.major_axis * math.cos(
                                          teta) * self.major_dir
                                      + self.minor_axis * math.sin(
                                          teta) * self.major_dir.cross(
                                          self.normal) for teta in
                                      npy.linspace(0, volmdlr.TWO_PI,
                                                   angle_resolution + 1)][:-1]
        return discretization_points_3d

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a Ellipse3D into an EllipseD, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Ellipse2D.
        """
        center = self.center.to_2d(plane_origin, x, y)
        major_dir_d2 = self.major_dir.to_2d(plane_origin, x, y)
        return Ellipse2D(self.major_axis, self.minor_axis, center, major_dir_d2)

    def abscissa(self, point: volmdlr.Point3D):
        """
        Calculates the abscissa a given point.

        :param point: point to calculate abscissa.
        :return: abscissa
        """
        vector_2 = self.normal.cross(self.major_dir)
        ellipse_2d = self.to_2d(self.center, self.major_dir, vector_2)
        point2d = point.to_2d(self.center, self.major_dir, vector_2)
        return ellipse_2d.abscissa(point2d)

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D):
        if point1.is_close(point2):
            return volmdlr.edges.FullArcEllipse3D(point1, self.major_axis, self.minor_axis, self.center, self.normal,
                                                  self.major_dir, self.name)

        p1_new, p2_new = self.frame.global_to_local_coordinates(point1), self.frame.global_to_local_coordinates(point2)

        theta1 = volmdlr.geometry.sin_cos_angle(p1_new.x / self.major_axis, p1_new.y / self.minor_axis)

        theta2 = volmdlr.geometry.sin_cos_angle(p2_new.x / self.major_axis, p2_new.y / self.minor_axis)

        if theta1 > theta2:  # sens trigo
            angle = math.pi + (theta1 + theta2) / 2
        else:
            angle = (theta1 + theta2) / 2

        point3 = self.frame.local_to_global_coordinates(volmdlr.Point3D(self.major_axis * math.cos(angle),
                                                                        self.minor_axis * math.sin(angle), 0))
        extra = None
        if math.isclose(angle % math.pi, 0.0, abs_tol=1e-6):
            extra = self.frame.local_to_global_coordinates(volmdlr.Point3D(self.major_axis * math.cos(0.125 * angle),
                                                                           self.minor_axis * math.sin(0.125 * angle),
                                                                           0))
        return volmdlr.edges.ArcEllipse3D(point1, point3, point2, self.center,
                                          self.major_dir, extra=extra)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Ellipse3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated Ellipse3D.
        """
        new_center = self.center.rotation(center, axis, angle)
        new_normal = self.normal.rotation(center, axis, angle)
        new_major_dir = self.major_dir.rotation(center, axis, angle)
        return Ellipse3D(self.major_axis, self.minor_axis, new_center,
                         new_normal, new_major_dir, self.name)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Ellipse3D rotation. Object is updated in-place.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        self.center.rotation_inplace(center, axis, angle)
        self.normal.rotation_inplace(center, axis, angle)
        self.major_dir.rotation_inplace(center, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Ellipse 3D translation.

        :param offset: translation vector.
        :return: A new translated Ellipse 3D.
        """
        new_center = self.center.translation(offset)
        # new_normal = self.normal.translation(offset)
        new_normal = self.normal
        new_major_dir = self.major_dir.translation(offset)
        return Ellipse3D(self.major_axis, self.minor_axis, new_center,
                         new_normal, new_major_dir, self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Ellipse3D translation. Object is updated in-place.

        :param offset: translation vector
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        self.center.translation_inplace(offset)
        self.normal.translation_inplace(offset)
        self.major_dir.translation_inplace(offset)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        x = []
        y = []
        z = []
        for point_x, point_y, point_z in self.discretization_points():
            x.append(point_x)
            y.append(point_y)
            z.append(point_z)
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])
        ax.plot(x, y, z, edge_style.color)
        return ax

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a Ellipse3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding Ellipse3D object.
        :rtype: :class:`volmdlr.wires.Ellipse3D`
        """
        length_conversion_factor = kwargs.get("length_conversion_factor", 1)

        center = object_dict[arguments[1]].origin
        normal = object_dict[arguments[1]].u  # ancien w
        major_dir = object_dict[arguments[1]].v  # ancien u
        major_axis = float(arguments[2]) * length_conversion_factor
        minor_axis = float(arguments[3]) * length_conversion_factor
        return cls(major_axis, minor_axis, center, normal, major_dir,
                   arguments[0][1:-1])


class ClosedPolygon3D(Contour3D, ClosedPolygonMixin):
    """
    A collection of points, connected by line segments, following each other.

    """
    _non_serializable_attributes = ['line_segments', 'primitives']
    _non_eq_attributes = ['line_segments', 'primitives']

    def __init__(self, points: List[volmdlr.Point3D], name: str = ''):
        self.points = points
        self._line_segments = None

        Contour3D.__init__(self, self.line_segments, name)

    def get_line_segments(self):
        lines = []
        if len(self.points) > 1:
            for point1, point2 in zip(self.points,
                                      list(self.points[1:]) + [self.points[0]]):
                if not point1.is_close(point2):
                    lines.append(volmdlr.edges.LineSegment3D(point1, point2))
        return lines

    def copy(self, *args, **kwargs):
        points = [point.copy() for point in self.points]
        return ClosedPolygon3D(points, self.name)

    def __hash__(self):
        return sum(hash(point) for point in self.points)

    def __eq__(self, other_):
        if not isinstance(other_, self.__class__):
            return False
        equal = True
        for point, other_point in zip(self.points, other_.points):
            equal = (equal and point.is_close(other_point))
        return equal

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        for line_segment in self.line_segments:
            ax = line_segment.plot(ax=ax, edge_style=edge_style)
        return ax

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        ClosedPolygon3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated ClosedPolygon3D.
        """
        return ClosedPolygon3D(
            [point.rotation(center, axis, angle) for point in
             self.points])

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        ClosedPolygon3D rotation. Object is updated in-place.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: rotation angle.
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for point in self.points:
            point.rotation_inplace(center, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        ClosedPolygon3D translation.

        :param offset: translation vector.
        :return: A new translated ClosedPolygon3D.
        """
        new_points = [point.translation(offset) for point in
                      self.points]
        return ClosedPolygon3D(new_points, self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        ClosedPolygon3D translation. Object is updated in-place.

        :param offset: translation vector.
        """
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for point in self.points:
            point.translation_inplace(offset)

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a ClosedPolygon3D into an ClosedPolygon2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: ClosedPolygon2D.
        """
        points2d = [point.to_2d(plane_origin, x, y) for point in self.points]
        return ClosedPolygon2D(points2d)

    def sewing_with(self, other_poly3d, x, y, resolution=20):
        self_center, other_center = self.average_center_point(), \
            other_poly3d.average_center_point()

        self_poly2d, other_poly2d = self.to_2d(self_center, x, y), \
            other_poly3d.to_2d(other_center, x, y)
        self_center2d, other_center2d = self_poly2d.center_of_mass(), \
            other_poly2d.center_of_mass()
        self_poly2d.translation_inplace(-self_center2d)
        other_poly2d.translation_inplace(-other_center2d)

        bbox_self2d, bbox_other2d = self_poly2d.bounding_rectangle.bounds(), \
            other_poly2d.bounding_rectangle.bounds()
        position = [abs(value) for value in bbox_self2d] \
            + [abs(value) for value in bbox_other2d]
        max_scale = 2 * max(position)

        lines = [volmdlr.edges.LineSegment2D(volmdlr.O2D, max_scale * (
                volmdlr.X2D * math.sin(n * 2 * math.pi / resolution) +
                volmdlr.Y2D * math.cos(n * 2 * math.pi / resolution))
                                             ) for n in range(resolution)]

        self_new_points, other_new_points = [], []
        for line in lines:
            for self_line in self_poly2d.line_segments:
                intersect = line.linesegment_intersections(self_line)
                if intersect:
                    self_new_points.extend(intersect)
                    break

            for other_line in other_poly2d.line_segments:
                intersect = line.linesegment_intersections(other_line)
                if intersect:
                    other_new_points.extend(intersect)
                    break

        new_self_poly2d, new_other_poly2d = ClosedPolygon2D(
            self_new_points), ClosedPolygon2D(other_new_points)
        new_self_poly2d.translation_inplace(self_center2d)
        new_other_poly2d.translation_inplace(other_center2d)

        new_poly1, new_poly2 = new_self_poly2d.to_3d(self_center, x, y), \
            new_other_poly2d.to_3d(other_center, x, y)

        triangles = []
        for point1, point2, other_point in zip(new_poly1.points,
                                               new_poly1.points[
                                                   1:] + new_poly1.points[:1],
                                               new_poly2.points):
            triangles.append([point1, point2, other_point])

        for point1, point2, other_point in zip(
                new_poly2.points, new_poly2.points[1:] + new_poly2.points[:1],
                new_poly1.points[1:] + new_poly1.points[:1]):
            triangles.append([other_point, point2, point1])

        return triangles

    def simplify(self, min_distance: float = 0.01, max_distance: float = 0.05):
        """
        Simplifies polygon 3d.

        :param min_distance: minimal allowed distance.
        :param max_distance: maximal allowed distance.
        :return: Simplified closed polygon 3d.
        """
        return ClosedPolygon3D(self.simplify_polygon(
            min_distance=min_distance, max_distance=max_distance).points)

    def convex_sewing(self, polygon2, x, y):
        """
        Sew to Convex Polygon.

        :param polygon2: other polygon to sew with.
        :param x: u vector for plane projection.
        :param y: v vector for plane projection.
        """
        center1, center2 = self.average_center_point(), polygon2.average_center_point()
        center1_, center2_ = volmdlr.Point3D(center1.x, center1.y, 0), volmdlr.Point3D(center2.x, center2.y, 0)
        new_polygon1, new_polygon2 = self.translation(-center1_), polygon2.translation(-center2_)
        new_center1, new_center2 = new_polygon1.average_center_point(), new_polygon2.average_center_point()

        new_polygon1_2d, new_polygon2_2d = \
            new_polygon1.to_2d(new_center1, x, y), new_polygon2.to_2d(new_center2, x, y)

        dict_closing_pairs = {}
        triangles = []
        list_closing_point_indexes = []
        new_polygon1_2d_points = new_polygon1_2d.points + [
            new_polygon1_2d.points[0]]
        for i, point_polygon1 in enumerate(
                new_polygon1.points + [new_polygon1.points[0]]):
            if i != 0:
                mean_point2d = 0.5 * (
                        new_polygon1_2d_points[i] + new_polygon1_2d_points[
                            i - 1])
                closing_point = new_polygon2_2d.line_intersecting_closing_point(
                    mean_point2d)
                closing_point_index = new_polygon2_2d.points.index(
                    closing_point)

                if i == 1:
                    previous_closing_point_index = closing_point_index
                if closing_point_index != previous_closing_point_index:
                    if closing_point_index in list_closing_point_indexes:
                        closing_point_index = previous_closing_point_index
                    else:
                        dict_closing_pairs[self.points[i - 1]] = (previous_closing_point_index, closing_point_index)

                if point_polygon1.is_close(new_polygon1.points[0]):
                    if list(dict_closing_pairs.values())[-1][-1] != list(dict_closing_pairs.values())[0][0]:
                        dict_closing_pairs[self.points[0]] = (list(dict_closing_pairs.values())[-1][-1],
                                                              list(dict_closing_pairs.values())[0][0])

                real_closing_point = polygon2.points[closing_point_index]

                face_points = [self.points[new_polygon1.points.index(
                    point_polygon1)], self.points[i - 1],
                               real_closing_point]
                triangles.append(face_points)

                list_closing_point_indexes.append(closing_point_index)
                previous_closing_point_index = closing_point_index
        triangles += polygon2.close_sewing(dict_closing_pairs)

        return triangles

    def get_valid_concave_sewing_polygon(self, polygon1_2d, polygon2_2d):
        polygon1_2d_valid__primitive = \
            polygon1_2d.get_valid_sewing_polygon_primitive(polygon2_2d)
        if polygon1_2d_valid__primitive == polygon1_2d.line_segments[0]:
            return self
        new_polygon_primitives = \
            self.line_segments[polygon1_2d.line_segments.index(polygon1_2d_valid__primitive):] + \
            self.line_segments[:polygon1_2d.line_segments.index(polygon1_2d_valid__primitive)]
        polygon1_3d_points = []
        for prim in new_polygon_primitives:
            if not volmdlr.core.point_in_list(prim.start, polygon1_3d_points):
                polygon1_3d_points.append(prim.start)
            if not volmdlr.core.point_in_list(prim.end, polygon1_3d_points):
                polygon1_3d_points.append(prim.end)
        return ClosedPolygon3D(polygon1_3d_points)

    def close_sewing(self, dict_closing_pairs):
        triangles_points = []
        for i, point_polygon2 in enumerate(
                self.points + [self.points[0]]):
            for j, index in enumerate(list(dict_closing_pairs.values())):
                if i != 0:
                    if i - 1 >= index[0] and i <= index[1]:
                        face_points = [self.points[i - 1],
                                       point_polygon2,
                                       list(dict_closing_pairs.keys())[j]]
                        triangles_points.append(face_points)
                    elif index[0] > index[1]:
                        if (i - 1 <= index[0] and i <= index[1]) or (
                                (i - 1 >= index[0]) and i >= index[1]):
                            face_points = [self.points[i - 1],
                                           point_polygon2,
                                           list(dict_closing_pairs.keys())[j]]
                            triangles_points.append(face_points)
        return triangles_points

    def check_sewing(self, polygon2, sewing_faces):
        if not len(self.line_segments) + len(polygon2.line_segments) == len(sewing_faces):
            return False
        return True

    def redefine_sewing_triangles_points(self, triangles_points,
                                         passed_by_zero_index,
                                         closing_point_index,
                                         previous_closing_point_index):
        for n, triangle_points in enumerate(triangles_points[::-1]):
            if (not passed_by_zero_index and
                self.points.index(
                    triangle_points[2]) > closing_point_index) or \
                    (passed_by_zero_index and
                     0 <= self.points.index(triangle_points[
                                                2]) <= previous_closing_point_index and
                     self.points.index(
                         triangle_points[2]) > closing_point_index):
                new_face_points = [triangles_points[-(n + 1)][0],
                                   triangles_points[-(n + 1)][1],
                                   self.points[
                                       closing_point_index]]
                triangles_points[-(n + 1)] = new_face_points

        return triangles_points

    @staticmethod
    def clean_sewing_closing_pairs_dictionary(dict_closing_pairs,
                                              closing_point_index,
                                              passed_by_zero_index):
        """
        Cleans the dictionary containing the sewing closing pairs information.

        In case it needs to be recalculated due to changing closing points.
        """
        dict_closing_pairs_values = list(dict_closing_pairs.values())
        dict_closing_pairs_keys = list(dict_closing_pairs.keys())
        previous_closing_point_index = dict_closing_pairs_values[-1][1]
        last_dict_value = previous_closing_point_index
        for i, key in enumerate(dict_closing_pairs_keys[::-1]):
            if (not passed_by_zero_index and
                last_dict_value > closing_point_index) or \
                    (passed_by_zero_index and
                     0 <= last_dict_value <= previous_closing_point_index and
                     last_dict_value > closing_point_index):
                lower_bounddary_closing_point = key
                del dict_closing_pairs[key]
                if not dict_closing_pairs:
                    break
                last_dict_value = dict_closing_pairs_values[-i - 2][1]

        return dict_closing_pairs, lower_bounddary_closing_point

    @staticmethod
    def is_sewing_forward(closing_point_index, list_closing_point_indexes) -> bool:
        if closing_point_index < list_closing_point_indexes[-1]:
            return False
        return True

    @staticmethod
    def sewing_closing_points_to_remove(closing_point_index, list_closing_point_indexes, passed_by_zero_index):
        list_remove_closing_points = []
        for idx in list_closing_point_indexes[::-1]:
            if not passed_by_zero_index:
                if idx > closing_point_index:
                    list_remove_closing_points.append(idx)
                else:
                    break
            else:
                if 0 < idx <= list_closing_point_indexes[-1] and \
                        idx > closing_point_index:
                    list_remove_closing_points.append(idx)
                else:
                    break
        return list_remove_closing_points

    @staticmethod
    def sewing_closing_point_past_point0(closing_point_index, list_closing_point_indexes,
                                         passed_by_zero_index, ratio_denominator):
        last_to_new_point_index_ratio = (list_closing_point_indexes[-1] -
                                         closing_point_index) / ratio_denominator
        if passed_by_zero_index:
            ratio = (list_closing_point_indexes[0] - closing_point_index) / ratio_denominator
            if math.isclose(ratio, 1, abs_tol=0.3):
                closing_point_index = list_closing_point_indexes[0]
            else:
                closing_point_index = list_closing_point_indexes[-1]
        else:
            if closing_point_index > list_closing_point_indexes[0]:
                ratio1 = (closing_point_index -
                          list_closing_point_indexes[0]) / ratio_denominator
                if math.isclose(ratio1, 0, abs_tol=0.3) and \
                        math.isclose(last_to_new_point_index_ratio, 1, abs_tol=0.3):
                    passed_by_zero_index = True
                    closing_point_index = list_closing_point_indexes[0]
                else:
                    closing_point_index = list_closing_point_indexes[-1]
            else:
                if closing_point_index < ratio_denominator / 4:
                    passed_by_zero_index = True
                elif ratio_denominator - list_closing_point_indexes[-1] >= 6:
                    closing_point_index = list_closing_point_indexes[-1] + 5
                else:
                    closing_point_index = list_closing_point_indexes[-1]
        return closing_point_index, passed_by_zero_index

    @staticmethod
    def validate_concave_closing_point(closing_point_index,
                                       list_closing_point_indexes,
                                       passed_by_zero_index,
                                       ratio_denominator, polygons_points_ratio):
        last_index = list_closing_point_indexes[-1]

        if closing_point_index == last_index:
            return closing_point_index, [], passed_by_zero_index

        list_remove_closing_points = []
        ratio = (last_index - closing_point_index) / ratio_denominator

        if not ClosedPolygon3D.is_sewing_forward(closing_point_index, list_closing_point_indexes):
            if closing_point_index > last_index - 10 and closing_point_index != last_index - 1:
                if closing_point_index - 1 in list_closing_point_indexes and \
                        closing_point_index + 1 in list_closing_point_indexes:
                    closing_point_index = last_index
                    return closing_point_index, list_remove_closing_points, passed_by_zero_index

                list_remove_closing_points = ClosedPolygon3D.sewing_closing_points_to_remove(
                    closing_point_index, list_closing_point_indexes, passed_by_zero_index)

            elif closing_point_index in list_closing_point_indexes:
                closing_point_index = last_index
            elif math.isclose(ratio, 0, abs_tol=0.3):
                closing_point_index = last_index
            else:
                closing_point_index, passed_by_zero_index = ClosedPolygon3D.sewing_closing_point_past_point0(
                    closing_point_index, list_closing_point_indexes, passed_by_zero_index, ratio_denominator)

        elif closing_point_index in list_closing_point_indexes:
            closing_point_index = last_index
        elif len(list_closing_point_indexes) > 2 and list_closing_point_indexes[0] < closing_point_index < last_index:
            closing_point_index = last_index
        elif passed_by_zero_index and closing_point_index > list_closing_point_indexes[0]:
            closing_point_index = last_index
        elif list_closing_point_indexes[0] == 0 and math.isclose(ratio, -1, abs_tol=0.3):
            closing_point_index = last_index
        elif math.isclose(ratio, -1, abs_tol=0.3):
            closing_point_index = last_index
        elif closing_point_index - last_index > 5 and list_closing_point_indexes[
                -1] + 4 <= ratio_denominator - 1 and polygons_points_ratio > 0.95:
            closing_point_index = last_index + 4

        return closing_point_index, list_remove_closing_points, passed_by_zero_index

    def concave_sewing(self, polygon2, x, y):
        polygon1_2d = self.to_2d(volmdlr.O2D, x, y)
        polygon2_2d = polygon2.to_2d(volmdlr.O2D, x, y)
        polygon1_3d = self
        polygon2_3d = polygon2
        if polygon2_2d.area() < polygon1_2d.area():
            polygon1_2d, polygon2_2d = polygon2_2d, polygon1_2d
            polygon1_3d = polygon2
            polygon2_3d = self
        polygon1_3d = polygon1_3d.get_valid_concave_sewing_polygon(
            polygon1_2d, polygon2_2d)
        polygon1_2d = polygon1_3d.to_2d(volmdlr.O2D, x, y)

        # ax=polygon1_2d.plot()
        # polygon2_2d.plot(ax=ax, color='r')

        dict_closing_pairs = {}
        triangles_points = []
        list_closing_point_indexes = []
        passed_by_zero_index = False
        ratio_denom = len(polygon2_2d.points)
        polygons_points_ratio = len(polygon1_2d.points) / ratio_denom
        previous_closing_point_index = None
        for i, primitive1 in enumerate(polygon1_2d.line_segments):
            list_remove_closing_points = []
            closing_point = polygon1_2d.get_closing_point(polygon2_2d,
                                                          primitive1)
            if closing_point.is_close(volmdlr.O2D):
                if previous_closing_point_index is not None:
                    closing_point_index = previous_closing_point_index
                else:
                    raise NotImplementedError(
                        'None of the normal lines intersect polygon2, '
                        'certify projection plane given is correct')
            else:
                closing_point_index = polygon2_2d.points.index(closing_point)

            if i == 0:
                previous_closing_point_index = closing_point_index
            else:
                closing_point_index, list_remove_closing_points, \
                    passed_by_zero_index = self.validate_concave_closing_point(
                        closing_point_index, list_closing_point_indexes,
                        passed_by_zero_index, ratio_denom, polygons_points_ratio)

            if list_remove_closing_points:
                new_list_closing_point_indexes = list(
                    dict.fromkeys(list_closing_point_indexes))
                new_list_remove_closing_indexes = list(
                    dict.fromkeys(list_remove_closing_points))
                if len(list_remove_closing_points) == len(triangles_points):
                    triangles_points = \
                        polygon2_3d.redefine_sewing_triangles_points(
                            triangles_points, passed_by_zero_index,
                            closing_point_index, previous_closing_point_index)
                    if dict_closing_pairs:
                        dict_closing_pairs, lower_bounddary_closing_point = \
                            self.clean_sewing_closing_pairs_dictionary(
                                dict_closing_pairs,
                                closing_point_index,
                                passed_by_zero_index)

                        if len(new_list_remove_closing_indexes) < \
                                len(new_list_closing_point_indexes):
                            dict_closing_pairs[
                                lower_bounddary_closing_point] = (
                                new_list_closing_point_indexes[
                                    -(len(new_list_remove_closing_indexes) + 1)],
                                closing_point_index)
                    for pt_index in list_remove_closing_points:
                        list_closing_point_indexes.remove(pt_index)
                    list_closing_point_indexes.append(closing_point_index)

                elif (not passed_by_zero_index and
                      closing_point_index > polygon2_3d.points.index(
                            triangles_points[-len(list_remove_closing_points) - 1][2])) or \
                        (passed_by_zero_index and closing_point_index >= 0):
                    triangles_points = \
                        polygon2_3d.redefine_sewing_triangles_points(
                            triangles_points, passed_by_zero_index,
                            closing_point_index, previous_closing_point_index)
                    dict_closing_pairs, lower_bounddary_closing_point = \
                        self.clean_sewing_closing_pairs_dictionary(
                            dict_closing_pairs, closing_point_index, passed_by_zero_index)

                    if not list(dict_closing_pairs.keys()) or dict_closing_pairs[
                        list(dict_closing_pairs.keys())[-1]][1] != \
                            closing_point_index:
                        dict_closing_pairs[lower_bounddary_closing_point] = \
                            (new_list_closing_point_indexes[
                                 -(len(new_list_remove_closing_indexes) + 1)],
                             closing_point_index)

                    for pt_index in list_remove_closing_points:
                        list_closing_point_indexes.remove(pt_index)
                    list_closing_point_indexes.append(closing_point_index)
                else:
                    closing_point_index = previous_closing_point_index

            elif closing_point_index != previous_closing_point_index:
                dict_closing_pairs[polygon1_3d.line_segments[i].start] = \
                    (previous_closing_point_index, closing_point_index)
            face_points = [polygon1_3d.line_segments[i].start,
                           polygon1_3d.line_segments[i].end,
                           polygon2_3d.points[closing_point_index]]
            triangles_points.append(face_points)
            list_closing_point_indexes.append(closing_point_index)
            previous_closing_point_index = closing_point_index
            if primitive1 == polygon1_2d.line_segments[-1]:
                if list_closing_point_indexes[-1] != list_closing_point_indexes[0]:
                    ratio = (list_closing_point_indexes[-1] -
                             list_closing_point_indexes[0]) / len(
                        polygon2_2d.points)
                    if math.isclose(ratio, -1,
                                    abs_tol=0.2) and passed_by_zero_index:
                        dict_closing_pairs[
                            polygon1_3d.points[0]] = (
                            list_closing_point_indexes[-2],
                            list_closing_point_indexes[0])
                        new_face_points = [triangles_points[-1][0],
                                           triangles_points[-1][1],
                                           polygon2_3d.points[
                                               list_closing_point_indexes[-2]]]
                        triangles_points.remove(triangles_points[-1])
                        triangles_points.append(new_face_points)
                    else:
                        dict_closing_pairs[polygon1_3d.points[0]] = (
                            list(dict_closing_pairs.values())[-1][-1],
                            list(dict_closing_pairs.values())[0][0])

        triangles_points += polygon2_3d.close_sewing(dict_closing_pairs)

        return triangles_points

    def sewing(self, polygon2, x, y):
        polygon1_2d = self.to_2d(volmdlr.O2D, x, y)
        polygon2_2d = polygon2.to_2d(volmdlr.O2D, x, y)
        if polygon1_2d.is_convex() and polygon2_2d.is_convex():
            return self.convex_sewing(polygon2, x, y)
        return self.concave_sewing(polygon2, x, y)


class Triangle3D(Triangle):
    """
    Defines a triangle 3D.

    :param point1: triangle point 1.
    :param point2: triangle point 2.
    :param point3: triangle point3.
    """

    def __init__(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D,
                 point3: volmdlr.Point3D, name: str = ''):
        Triangle.__init__(self, point1,
                          point2,
                          point3,
                          name)
