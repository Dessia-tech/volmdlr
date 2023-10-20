#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing wires & contours.
"""

import itertools
import math
import sys
import warnings
from collections import deque
from functools import cached_property
from statistics import mean
from typing import List

import matplotlib.pyplot as plt
import networkx as nx
import numpy as npy
import plot_data.core as plot_data
from scipy.spatial.qhull import ConvexHull, Delaunay
from triangle import triangulate

from dessia_common.core import PhysicalObject
import volmdlr
import volmdlr.core
import volmdlr.display as vmd
import volmdlr.geometry
from volmdlr import curves, edges
from volmdlr.core_compiled import polygon_point_belongs, points_in_polygon
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


def reorder_contour3d_edges_from_step(raw_edges, step_data):
    """Helper function to order edges from a 3D contour coming from a step file."""
    step_id, step_name, arguments = step_data
    reversed_distances = [edge1.start.point_distance(edge2.end)
                          for edge1, edge2 in zip(raw_edges[::-1][1:], raw_edges[::-1][:-1])]
    if all((dist < 1e-6) for dist in reversed_distances):
        return raw_edges[::-1]

    # Making things right for first 2 primitives
    distances = [raw_edges[0].end.point_distance(raw_edges[1].start),
                 raw_edges[0].start.point_distance(raw_edges[1].start),
                 raw_edges[0].end.point_distance(raw_edges[1].end),
                 raw_edges[0].start.point_distance(raw_edges[1].end)]
    index = distances.index(min(distances))
    if min(distances) > 1e-3:
        # # Green color : well-placed and well-read
        # ax = raw_edges[0].plot(edge_style=EdgeStyle(color='g'))
        # ax.set_title(f"Step ID: #{step_id}")
        #
        # # Red color : can't be connected to green edge
        # raw_edges[1].plot(ax=ax, edge_style=EdgeStyle(color='r'))
        # # Black color : to be placed
        # for re in raw_edges[2:]:
        #     re.plot(ax=ax)

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
        new_edges = [raw_edges[0], raw_edges[1]]
    elif index == 1:
        new_edges = [raw_edges[0].reverse(), raw_edges[1]]
    elif index == 2:
        new_edges = [raw_edges[0], raw_edges[1].reverse()]
    elif index == 3:
        new_edges = [raw_edges[0].reverse(), raw_edges[1].reverse()]
    else:
        raise NotImplementedError

    # Connecting the next edges
    last_edge = new_edges[-1]
    for raw_edge in raw_edges[2:]:
        if raw_edge.direction_independent_is_close(last_edge):
            continue
        distances = [raw_edge.start.point_distance(last_edge.end),
                     raw_edge.end.point_distance(last_edge.end)]
        index = distances.index(min(distances))
        if min(distances) > 1e-3:
            # Green color : well-placed and well-read
            # ax = last_edge.plot(edge_style=EdgeStyle(color='g'))
            # ax.set_title(f"Step ID: #{step_id}")
            #
            # for re in raw_edges[:2 + i]:
            #     re.plot(ax=ax, edge_style=EdgeStyle(color='g'))
            #     re.start.plot(ax=ax, color='g')
            #     re.end.plot(ax=ax, color='g')
            # last_edge.end.plot(ax=ax, color='g')
            # # Red color : can't be connected to red dot
            # raw_edge.plot(ax=ax, edge_style=EdgeStyle(color='g'))
            # # Black color : to be placed
            # for re in raw_edges[2 + i + 1:]:
            #     re.plot(ax=ax)
            #     re.start.plot(ax=ax)
            #     re.end.plot(ax=ax)

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

        new_edges.append(last_edge)
    return new_edges


class WireMixin:
    """
    Abstract class for Wire, storing methods and attributes used by many classes in this module.

    """
    _non_data_hash_attributes = ['basis_primitives']
    _non_serializable_attributes = ['primitive_to_index',
                                    'basis_primitives']

    def _data_hash(self):
        return sum(hash(e) for e in self.primitives) + len(self.primitives)

    def to_dict(self, *args, **kwargs):
        """Avoids storing points in memo that makes serialization slow."""
        return PhysicalObject.to_dict(self, use_pointers=False)

    def length(self):
        """Returns the wire's length."""
        if not self._length:
            length = 0.
            for primitive in self.primitives:
                length += primitive.length()
            self._length = length
        return self._length

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

        :param point1: splitting point1.
        :param point2: splitting point2.
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

    def abscissa(self, point, tol: float = 1e-6):
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

    def ordering_primitives(self, tol: float = 1e-6):
        """
        Ordering wire / contour primitives.

        :param tol: tolerance.
        :return:
        """
        primitives = [prim for prim in self.primitives[:] if not math.isclose(prim.length(), 0.0, abs_tol=tol)]
        new_primitives = [primitives[0]]
        primitives.remove(primitives[0])
        while True:
            if not primitives:
                break
            for primitive in primitives:
                if primitive.is_point_edge_extremity(new_primitives[-1].end, tol):
                    if new_primitives[-1].end.is_close(primitive.start, tol):
                        new_primitives.append(primitive)
                    else:
                        new_primitives.append(primitive.reverse())
                    primitives.remove(primitive)
                    break
                if primitive.is_point_edge_extremity(new_primitives[0].start, tol):
                    if new_primitives[0].start.is_close(primitive.end, tol):
                        new_primitives.insert(0, primitive)
                    else:
                        new_primitives.insert(0, primitive.reverse())
                    primitives.remove(primitive)
                    break
            else:
                # print(self, primitive)
                # ax = self.plot()
                # primitive.plot(edge_style=EdgeStyle(color='r'), ax=ax)
                raise NotImplementedError('There may exist a problem with this'
                                          ' contour, it seems it cannot be reordered.'
                                          ' Please, verify its points')
        return new_primitives

    def order_wire(self, tol=1e-6):
        """ Order wire's primitives. """

        if self.is_ordered(tol=tol):
            return self
        return self.__class__(self.ordering_primitives(tol), name=self.name)

    @classmethod
    def from_wires(cls, wires, name: str = ''):
        """
        Define a wire from successive wires.

        """

        primitives = []
        for wire in wires:
            primitives.extend(wire.primitives)

        wire = cls(primitives, name=name)

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

    def invert(self):
        """Gets the wire in the inverted direction."""
        return self.__class__(self.inverted_primitives())

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
        Verifies if primitive is over wire.

        :param primitive: primitive to be verified.
        :param tol: tolerance to be considered.
        :return: True or False
        """
        points = primitive.discretization_points(number_points=10)
        if all(self.point_over_wire(point, tol) for point in points):
            return True
        return False

    def is_primitive_section_over_wire(self, primitive, tol: float = 1e-6):
        """
        Verifies if primitive's section is over wire.

        :param primitive: primitive to be verified.
        :param tol: tolerance to be considered.
        :return: True or False
        """
        for edge in self.primitives:
            shared_section = edge.get_shared_section(primitive, tol)
            if shared_section:
                return True
        return False

    @classmethod
    def from_points(cls, points, name: str = ''):
        """
        Create a contour from points with line_segments.

        """
        linesegment_name = 'LineSegment' + points[0].__class__.__name__[-2:]
        primitives = []
        for i in range(0, len(points) - 1):
            primitives.append(getattr(edges, linesegment_name)(points[i], points[i + 1]))
        contour = cls(primitives, name=name)
        return contour

    @classmethod
    def from_edge(cls, edge, number_segments: int, name: str = ''):
        """
        Creates a Wire object from an edge.

        :param edge: edge used to create Wire.
        :param number_segments: number of segment for the wire to have.
        :param name: object's name.
        :return: Wire object.
        """
        points = edge.discretization_points(number_points=number_segments + 1)
        class_name_ = 'Wire' + edge.__class__.__name__[-2:]
        class_ = getattr(sys.modules[__name__], class_name_)
        return class_.from_points(points, name=name)

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

    @classmethod
    def extract(cls, contour, point1, point2, inside=False, name: str = ''):
        """Extracts a wire from another contour/wire, given two points."""
        new_primitives = contour.extract_with_points(point1, point2, inside)

        if cls.__name__[:-2] in ['Contour', 'Wire']:
            wires = [cls(new_primitives, name=name)]
        else:
            wire_class_ = getattr(sys.modules[__name__], 'Wire' + cls.__name__[-2:])
            wires = [wire_class_(new_primitives, name=name)]
        return wires

    def split_with_sorted_points(self, sorted_points):
        """
        Split contour in various sections using a list of sorted points along the contour.

        :param sorted_points: sorted list of points.
        :return: list of Contour sections.
        """
        self_start_equal_to_end = True
        if not self.primitives[0].start.is_close(self.primitives[-1].end):
            self_start_equal_to_end = False
            if not volmdlr.core.point_in_list(self.primitives[0].start, sorted_points):
                sorted_points = [self.primitives[0].start] + sorted_points
            if not volmdlr.core.point_in_list(self.primitives[-1].end, sorted_points):
                sorted_points.append(self.primitives[-1].end)
        if not self_start_equal_to_end:
            if len(sorted_points) == 2 and sorted_points[0].is_close(self.primitives[0].start) and \
                    sorted_points[1].is_close(self.primitives[-1].end):
                return [self]
        split_wires = []
        len_sorted_points = len(sorted_points)
        for i, (point1, point2) in enumerate(
                zip(sorted_points, sorted_points[1:] + [sorted_points[0]])):
            if i == len_sorted_points - 1:
                if self_start_equal_to_end:
                    split_wires.extend([wire.order_wire() for wire in
                                        self.__class__.extract(self, point1, point2, False)])
            else:
                split_wires.extend([wire.order_wire() for wire in
                                    self.__class__.extract(self, point1, point2, True)])
        return split_wires

    @classmethod
    def wires_from_edges(cls, list_edges, tol=1e-6, name: str = ''):
        """
        Defines a list of wires from edges, by ordering successive edges.

        :param list_edges: A list of edges
        :type list_edges: List[edges.Edge]
        :param tol: A tolerance, defaults to 1e-6
        :type tol: float, optional

        :return: A list of wires
        :param name: object's name.
        :rtype: List[wires.WireMixin]
        """

        if not list_edges:
            return []
        if len(list_edges) == 1:
            return [cls(list_edges, name=name)]

        new_primitives, i = [], -1
        index_primitive = 0
        while list_edges:
            i += 1
            new_primitives.append([list_edges[0]])
            list_edges.remove(list_edges[0])

            to_continue = True

            while to_continue:
                broke = False
                for index_primitive, primitive in enumerate(list_edges):

                    if primitive.is_point_edge_extremity(new_primitives[i][-1].end, tol):
                        if new_primitives[i][-1].end.is_close(primitive.start, tol):
                            new_primitives[i].append(primitive)
                        else:
                            new_primitives[i].append(primitive.reverse())
                        list_edges.remove(primitive)
                        broke = True
                        break

                    if primitive.is_point_edge_extremity(new_primitives[i][0].start, tol):
                        if new_primitives[i][0].start.is_close(primitive.end, tol):
                            new_primitives[i].insert(0, primitive)
                        else:
                            new_primitives[i].insert(0, primitive.reverse())
                        list_edges.remove(primitive)
                        broke = True
                        break

                if ((not broke) and (len(list_edges) == index_primitive + 1)) or len(list_edges) == 0:
                    to_continue = False

        wires = [cls(primitives_wire) for primitives_wire in new_primitives]

        return wires

    def to_wire_with_linesegments(self, number_segments: int):
        """
        Convert a wire with different primitives to a wire with just line segments by discretizing primitives.

        :param number_segments: number of segment for each primitive to be converted.
        :type number_segments: int
        """

        primitives = []
        class_name_ = 'Wire' + self.primitives[0].__class__.__name__[-2:]
        class_ = getattr(sys.modules[__name__], class_name_)

        for primitive in self.primitives:
            if primitive.__class__.__name__[0:-2] != 'LineSegment':
                primitives.extend(class_.from_edge(
                    edge=primitive, number_segments=number_segments).primitives)
            else:
                primitives.append(primitive)

        return class_(primitives)

    def get_connected_wire(self, list_wires):
        """
        Searches a wire in list_contour connected to self.

        :param list_wires: list of wires.
        :return:
        """
        connecting_contour_end = self.primitives[-1].end
        connecting_contour_start = self.primitives[0].start
        connected_contour = None
        for contour in list_wires:
            if connecting_contour_end.is_close(contour.primitives[0].start) or\
                    connecting_contour_end.is_close(contour.primitives[-1].end):
                connected_contour = contour
                break
            if connecting_contour_start.is_close(contour.primitives[0].start) or\
                    connecting_contour_start.is_close(contour.primitives[-1].end):
                connected_contour = contour
                break
        return connected_contour

    def is_sharing_primitives_with(self, contour, abs_tol: float = 1e-6):
        """
        Check if two contour are sharing primitives.

        """
        for prim1 in self.primitives:
            for prim2 in contour.primitives:
                shared_section = prim1.get_shared_section(prim2, abs_tol)
                if shared_section:
                    return True
        return False

    def middle_point(self):
        return self.point_at_abscissa(self.length() / 2)

    def is_superposing(self, contour2):
        """
        Check if the contours are superposing (one on the other without necessarily having an absolute equality).

        """

        for primitive_2 in contour2.primitives:
            if not self.primitive_over_wire(primitive_2):
                return False
        return True

    @classmethod
    def from_circle(cls, circle, name: str = ''):
        """
        Creates a Contour from a circle.

        :param circle: Circle.
        :param name: object's name.
        :return:
        """
        point = circle.point_at_abscissa(0.0)
        return cls([circle.trim(point, point)], name=name)

    def plot(self, ax=None, edge_style=EdgeStyle()):
        """Wire plot using Matplotlib."""
        if ax is None:
            ax = self._get_plot_ax()

        for element in self.primitives:
            element.plot(ax=ax, edge_style=edge_style)

        ax.margins(0.1)
        plt.show()

        return ax

    def edge_intersections(self, edge, abs_tol: float = 1e-6):
        """
        Compute intersections between a wire (2D or 3D) and an edge (2D or 3D).

        :param edge: edge to compute intersections.
        """
        edge_intersections = []
        for primitive in self.primitives:
            intersections = primitive.intersections(edge, abs_tol)
            for intersection in intersections:
                if not volmdlr.core.point_in_list(intersection, edge_intersections):
                    edge_intersections.append(intersection)
        return edge_intersections

    def wire_intersections(self, wire, abs_tol: float = 1e-6):
        """
        Compute intersections between two wire 2d.

        :param wire: volmdlr.wires.Wire2D
        """
        intersections_points = []
        for primitive in wire.primitives:
            edge_intersections = self.edge_intersections(primitive, abs_tol)
            for crossing in edge_intersections:
                if not volmdlr.core.point_in_list(crossing, intersections_points):
                    intersections_points.append(crossing)
        return intersections_points


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
        """ Plot edges with Matplotlib, not tested. """
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

    def babylon_meshes(self, **kwargs):
        """ Set the mesh for babylonjs. """
        positions, indices = self.to_babylon()
        babylon_mesh = {'positions': positions,
                        'indices': indices}
        babylon_mesh.update(self.babylon_param())
        return [babylon_mesh]


class Wire2D(WireMixin, PhysicalObject):
    """
    A collection of simple primitives, following each other making a wire.

    """

    def __init__(self, primitives, name: str = ''):
        self._bounding_rectangle = None
        self._length = None
        self.primitives = primitives
        PhysicalObject.__init__(self, name=name)

    def __hash__(self):
        return hash(('wire2d', tuple(self.primitives)))

    def area(self):
        """ Gets the area for a Wire2D."""
        return 0.0

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

    def line_intersections(self, line: 'curves.Line2D'):
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

    def is_start_end_crossings_valid(self, crossing_primitive, intersection, primitive):
        """
        Returns if the crossings are valid.

        :param crossing_primitive: crossing primitive.
        :param intersection: intersection result.
         for primitive line intersections
        :param primitive: intersecting primitive
        :return: None if intersection not a start or
        end point of a contours primitives, or a volmdlr.Point2D if it is.
        """
        primitive_index = self.primitives.index(primitive)
        point1, point2 = None, None
        if intersection.is_close(primitive.start):
            point1 = primitive.point_at_abscissa(primitive.length() * 0.01)
            point2 = self.primitives[primitive_index - 1].point_at_abscissa(
                self.primitives[primitive_index - 1].length() * .99
            )
        elif intersection.is_close(primitive.end) and primitive != self.primitives[-1]:
            point1 = primitive.point_at_abscissa(primitive.length() * 0.99)
            point2 = self.primitives[primitive_index + 1].point_at_abscissa(
                self.primitives[primitive_index + 1].length() * .01)
        if point1 is not None and point2 is not None:
            if not point1.is_close(point2):
                lineseg = volmdlr.edges.LineSegment2D(point1, point2)
                inter = crossing_primitive.linesegment_intersections(lineseg)
                if inter:
                    return True
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

    def line_crossings(self, line: curves.Line2D):
        """
        Calculates valid crossing intersections of a wire and an infinite line.

        :param line: line crossing the wire
        :type line: curves.Line2D
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
                    elif self.is_start_end_crossings_valid(line, intersection, primitive):
                        intersection_points.append(intersection)
                        intersection_points_primitives.append((intersection, primitive))
        return intersection_points_primitives

    def linesegment_crossings(self, linesegment: 'volmdlr.edges.LineSegment2D'):
        """
        Gets the wire primitives intersecting with the line.

        Returns a list of crossings in the form of a tuple (point, primitive).
        """
        results = self.line_crossings(linesegment.line)
        crossings_points = []
        for result in results:
            if linesegment.point_belongs(result[0]):
                crossings_points.append(result)
        return crossings_points

    def validate_edge_crossings(self, crossings):
        """
        Validates the crossings points from an edge and a wire.

        :param crossings: list of crossing points found.
        :return: list of valid crossing points.
        """
        crossings_ = []
        first_primitive = self.primitives[0]
        last_primitive = self.primitives[-1]
        for point in crossings:
            if not first_primitive.start.is_close(point) and not last_primitive.end.is_close(point):
                crossings_.append(point)
        return crossings_

    def edge_crossings(self, edge):
        """
        Gets the crossings between an edge and a Wire.

        :param edge: edge to search for crossings.
        :return: list of points containing all crossing points.
        """
        edge_crossings = []
        start_equal_to_end = self.primitives[0].start.is_close(self.primitives[-1].end)
        for primitive in self.primitives:
            crossings = primitive.intersections(edge)
            if not start_equal_to_end:
                crossings = self.validate_edge_crossings(crossings)
            for crossing in crossings:
                if not edge.is_point_edge_extremity(crossing) and\
                        not volmdlr.core.point_in_list(crossing, edge_crossings):
                    edge_crossings.append(crossing)
        return edge_crossings

    def validate_wire_crossing(self, crossing, current_wire_primitive, next_wire_primitive):
        """
        Validate the crossing point for the operation wire crossings.

        :param crossing: crossing point.
        :param current_wire_primitive: current wire primitive intersecting wire.
        :param next_wire_primitive: next wire primitive intersecting wire.
        :return:
        """
        self_primitives_to_test = [prim for prim in self.primitives if prim.is_point_edge_extremity(crossing)]
        if len(self_primitives_to_test) < 2:
            self_primitive = [prim for prim in self.primitives if prim.point_belongs(crossing)][0]
            crossing_abscissa = self_primitive.abscissa(crossing)
            vector_crossing = self_primitive.direction_vector(crossing_abscissa)
            current_vector = current_wire_primitive.direction_vector(current_wire_primitive.length())
            next_vector = next_wire_primitive.direction_vector(next_wire_primitive.length())
            if math.isclose(abs(vector_crossing.dot(current_vector)), 1, abs_tol=1e-6) or\
                    math.isclose(abs(vector_crossing.dot(next_vector)), 1, abs_tol=1e-6):
                return False
            return True
        if len(self_primitives_to_test) > 2:
            raise NotImplementedError
        if self_primitives_to_test[0] == self.primitives[0] and self_primitives_to_test[1] == self.primitives[-1]:
            point1 = self_primitives_to_test[0].point_at_abscissa(self_primitives_to_test[0].length() * 0.01)
            point2 = self_primitives_to_test[1].point_at_abscissa(self_primitives_to_test[1].length() * .99)
            point3 = current_wire_primitive.point_at_abscissa(current_wire_primitive.length() * .99)
            point4 = next_wire_primitive.point_at_abscissa(next_wire_primitive.length() * 0.01)
            linesegment1 = volmdlr.edges.LineSegment2D(point1, point2)
            linesegment2 = volmdlr.edges.LineSegment2D(point3, point4)
        else:
            point1 = self_primitives_to_test[0].point_at_abscissa(self_primitives_to_test[0].length() * .99)
            point2 = self_primitives_to_test[1].point_at_abscissa(self_primitives_to_test[1].length() * 0.01)
            point3 = current_wire_primitive.point_at_abscissa(current_wire_primitive.length() * .99)
            point4 = next_wire_primitive.point_at_abscissa(next_wire_primitive.length() * 0.01)
            linesegment1 = volmdlr.edges.LineSegment2D(point1, point2)
            linesegment2 = volmdlr.edges.LineSegment2D(point3, point4)
        inter = linesegment1.linesegment_intersections(linesegment2)
        if inter:
            return True
        return False

    def wire_crossings(self, wire):
        """
        Compute crossings between two wire 2d.

        :param wire: volmdlr.wires.Wire2D
        :return: crossing points: List[(volmdlr.Point2D)]
        """
        self_start_equal_to_end = self.primitives[0].start.is_close(self.primitives[-1].end)
        wire_start_equal_to_end = wire.primitives[0].start.is_close(wire.primitives[-1].end)
        wire_primitives = wire.primitives
        if wire_start_equal_to_end:
            wire_primitives = wire.primitives + [wire.primitives[0]]
        crossings_points = []
        len_wire_primitives = len(wire_primitives)
        invalid_crossings = []
        for i_prim, primitive in enumerate(wire_primitives):
            edge_intersections = self.edge_intersections(primitive)
            if not self_start_equal_to_end:
                edge_intersections = self.validate_edge_crossings(edge_intersections)
            if not wire_start_equal_to_end:
                edge_intersections = wire.validate_edge_crossings(edge_intersections)
            for crossing in edge_intersections:
                if i_prim != len_wire_primitives - 1:
                    if not self.validate_wire_crossing(crossing, primitive, wire_primitives[i_prim + 1]):
                        continue
                    if not volmdlr.core.point_in_list(crossing, crossings_points) and\
                            not volmdlr.core.point_in_list(crossing, invalid_crossings):
                        crossings_points.append(crossing)
        return crossings_points

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

        linesegments = Wire2D.from_edge(bsplinecurve, 25).primitives
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

        linesegments = Wire2D.from_edge(bsplinecurve, 25).primitives
        intersections_points = []
        for linesegment in linesegments:
            intersections_linesegments = self.linesegment_intersections(linesegment)
            if intersections_linesegments:
                intersections_points.extend(intersections_linesegments)
        return intersections_points

    @property
    def bounding_rectangle(self):
        """
        Returns the bounding rectangle of the wire.

        This property returns the bounding rectangle of the wire. If the bounding rectangle has not been calculated
        yet, it is computed using the `get_bouding_rectangle` method and stored in the `_bounding_rectangle` attribute.
        Subsequent calls to this property will return the pre-calculated bounding rectangle.

        :return: The bounding rectangle of the wire.
        :rtype: volmdlr.core.BoundingRectangle.
        """
        if not self._bounding_rectangle:
            self._bounding_rectangle = self.get_bouding_rectangle()
        return self._bounding_rectangle

    def get_bouding_rectangle(self):
        """
        Calculates the bounding rectangle of the wire.

        This method calculates the bounding rectangle of the wire. It initializes the minimum and maximum values
        for the x and y coordinates using the bounds of the first primitive. Then, it iterates over the remaining
        primitives and updates the minimum and maximum values based on their bounds. The resulting bounding rectangle
        is returned as a `volmdlr.core.BoundingRectangle` object.

        :return: The bounding rectangle of the wire.
        :rtype: volmdlr.core.BoundingRectangle.
        """

        x_min, x_max, y_min, y_max = self.primitives[0].bounding_rectangle.bounds()
        for edge in self.primitives[1:]:
            xmin_edge, xmax_edge, ymin_edge, ymax_edge = \
                edge.bounding_rectangle.bounds()
            x_min = min(x_min, xmin_edge)
            x_max = max(x_max, xmax_edge)
            y_min = min(y_min, ymin_edge)
            y_max = max(y_max, ymax_edge)
        return volmdlr.core.BoundingRectangle(x_min, x_max, y_min, y_max)

    def is_inside(self, other_contour):
        """
        Verifies if given contour is inside self contour perimeter, including its edges.

        :param other_contour: other contour.
        :returns: True or False
        """
        return False

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Rotates the wire 2D.

        :param center: rotation center.
        :param angle: angle rotation.
        :return: a new rotated Wire 2D.
        """
        return self.__class__([point.rotation(center, angle)
                               for point in self.primitives])

    def translation(self, offset: volmdlr.Vector2D):
        """
        Translates the Wire 2D.

        :param offset: translation vector
        :return: A new translated Wire 2D.
        """
        return self.__class__([primitive.translation(offset)
                               for primitive in self.primitives])

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and return a new Wire 2D.

        side = 'old' or 'new'
        """
        return self.__class__([primitive.frame_mapping(frame, side)
                               for primitive in self.primitives])

    def plot(self, ax=None, edge_style=EdgeStyle()):
        """Wire 2D plot using Matplotlib."""
        if ax is None:
            _, ax = plt.subplots()

        if edge_style.equal_aspect:
            ax.set_aspect('equal')

        for element in self.primitives:
            element.plot(ax=ax, edge_style=edge_style)

        ax.margins(0.1)
        b_rectangle = self.bounding_rectangle
        xlim, ylim = (b_rectangle[0] - 0.1, b_rectangle[1] + 0.1), (b_rectangle[2] - 0.1, b_rectangle[3] + 0.1)
        ax.set(xlim=xlim, ylim=ylim)
        plt.show()
        return ax

    def _get_plot_ax(self):
        _, ax = plt.subplots()
        return ax


class Wire3D(WireMixin, PhysicalObject):
    """
    A collection of simple primitives, following each other making a wire.

    """

    def __init__(self, primitives: List[volmdlr.core.Primitive3D], color=None, alpha: float = 1.0,
                 name: str = ''):
        self._bbox = None
        self._length = None
        self.primitives = primitives
        self.color = color
        self.alpha = alpha
        PhysicalObject.__init__(self, name=name)

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
        """
        Gets the distance from a point a Wire 3D object.

        :param point: other point.
        :return: the distance to wire and corresponding point.
        """
        distance, distance_point = math.inf, None
        for prim in self.primitives:
            prim_distance, prim_point = prim.point_distance(point)
            if prim_distance < distance:
                distance = prim_distance
                distance_point = prim_point
        return distance, distance_point

    def extrusion(self, extrusion_vector):
        """
        Extrudes a Wire 3D in a given direction.

        :param extrusion_vector: extrusion vector used.
        :return: A list of extruded faces.
        """
        faces = []
        for primitive in self.primitives:
            faces.extend(primitive.extrusion(extrusion_vector))
        return faces

    def to_bspline(self, discretization_parameter, degree):
        """
        Convert a wire 3d to a bspline curve 3d.

        """

        discretized_points = self.discretization_points(number_points=discretization_parameter)
        bspline_curve = volmdlr.edges.BSplineCurve3D.from_points_interpolation(discretized_points, degree, self.name)

        return bspline_curve

    def triangulation(self):
        """Triangulation method for a Wire3D."""
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
        plane3d = volmdlr.surfaces.Plane3D(volmdlr.Frame3D(plane_origin, x, y, z))
        primitives2d = []
        for primitive in self.primitives:
            primitive2d = plane3d.point3d_to_2d(primitive)
            if primitive2d:
                primitives2d.append(primitive2d)
        return primitives2d

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a Wire 3D into a Wire 2D, given a plane origin and an x and y vector.

        """
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

    def _get_plot_ax(self):
        _, ax = plt.subplots(subplot_kw={"projection": "3d"})
        return ax

    def babylon_points(self):
        """
        Returns a list of discretization points from the 3D primitive.
        """
        points = []
        if hasattr(self, 'primitives') and hasattr(self.primitives[0], "discretization_points"):
            for primitive in self.primitives:
                points.extend([*point] for point in primitive.discretization_points())
        elif hasattr(self, "discretization_points"):
            points.extend([*point] for point in self.discretization_points())
        return points

    def babylon_lines(self, points=None):
        if points is None:
            points = self.babylon_points()
        babylon_lines = {'points': points,
                         'alpha': self.alpha,
                         'name': self.name,
                         'color': list(self.color) if self.color is not None else [0.8, 0.8, 0.8]
                         }
        return [babylon_lines]

    def babylon_curves(self):
        points = self.babylon_points()
        if points:
            babylon_curves = self.babylon_lines(points)[0]
            return babylon_curves
        return None


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
        if len(self.primitives) == 2 and self.primitives[0].direction_independent_is_close(self.primitives[1]):
            return False
        for prim1, prim2 in zip(self.primitives, self.primitives[1:] + [self.primitives[0]]):
            if not prim1.end.is_close(prim2.start, tol):
                return False
        return True

    def order_contour(self, tol: float = 1e-6):
        """
        Verifies if the contours' primitives are ordered (one after the other). If not, it will order it.

        """
        if self.is_ordered(tol=tol) or len(self.primitives) < 2:
            return self
        new_primitives = self.ordering_primitives(tol)
        self.primitives = new_primitives

        return self

    @staticmethod
    def touching_edges_pairs(list_edges):  # TO DO: move this to edges?
        touching_primitives = []
        for i, primitive1 in enumerate(list_edges):
            for j, primitive2 in enumerate(list_edges):
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
    def find_connected_edges(list_edges, contours_list, contour_primitives, tol):
        for line in list_edges:
            if contours_list:
                contours_list, remove = ContourMixin.connected_to_splited_primitives(line, contours_list)
                if remove:
                    list_edges.remove(line)
                    break
            if not contour_primitives:
                contour_primitives.append(line)
                list_edges.remove(line)
                break
            edge_connected = ContourMixin.is_edge_connected(contour_primitives, line, tol)
            if edge_connected is not None:
                contour_primitives.append(edge_connected)
                list_edges.remove(edge_connected)
                break
        return list_edges, contour_primitives, contours_list

    @staticmethod
    def get_edges_bifurcations(contour_primitives, list_edges, finished_loop):
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
                        list_edges.append(volmdlr.edges.LineSegment2D(node, neihgbor))
                        finished_loop = False
                        if i + 1 == degree - 2:
                            break
        return contour_primitives, list_edges, finished_loop

    @classmethod
    def contours_from_edges(cls, list_edges, tol=1e-6, name: str = 'r'):
        if not list_edges:
            return []
        if len(list_edges) == 1:
            return [cls(list_edges, name=name)]
        list_contours = []
        points = [list_edges[0].start, list_edges[0].end]
        contour_primitives = [list_edges.pop(0)]
        while True:
            for i, edge in enumerate(list_edges):
                if (edge.is_point_edge_extremity(contour_primitives[-1].end, tol) and
                        not edge.direction_independent_is_close(contour_primitives[-1])):
                    if contour_primitives[-1].end.is_close(edge.start, tol):
                        contour_primitives.append(edge)
                    else:
                        contour_primitives.append(edge.reverse())
                    list_edges.pop(i)
                    validating_points = points[:]
                    validating_point = contour_primitives[-1].end
                    points.append(contour_primitives[-1].end)
                    break
                if (edge.is_point_edge_extremity(contour_primitives[0].start, tol) and
                        not edge.direction_independent_is_close(contour_primitives[0])):
                    if contour_primitives[0].start.is_close(edge.end, tol):
                        contour_primitives.insert(0, edge)
                    else:
                        contour_primitives.insert(0, edge.reverse())
                    validating_points = points[:]
                    validating_point = contour_primitives[0].start
                    points.insert(0, contour_primitives[0].start)
                    list_edges.pop(i)
                    break
            else:
                list_contours.append(cls(contour_primitives))
                if not list_edges:
                    break
                points = [list_edges[0].start, list_edges[0].end]
                contour_primitives = [list_edges.pop(0)]
                continue
            if volmdlr.core.point_in_list(validating_point, validating_points):
                if not validating_point.is_close(validating_points[0]):
                    spliting_primitives_index = volmdlr.core.get_point_index_in_list(
                        validating_point, validating_points)
                    if validating_point == points[0]:
                        new_contour = cls(contour_primitives[:spliting_primitives_index + 1])
                        contour_primitives = contour_primitives[spliting_primitives_index + 1:]
                        points = points[spliting_primitives_index + 1:]
                    else:
                        new_contour = cls(contour_primitives[spliting_primitives_index:])
                        contour_primitives = contour_primitives[:spliting_primitives_index]
                        points = points[:spliting_primitives_index + 1]
                    list_contours.append(new_contour)
                else:
                    list_contours.append(cls(contour_primitives))
                    if list_edges:
                        points = [list_edges[0].start, list_edges[0].end]
                        contour_primitives = [list_edges.pop(0)]
                    else:
                        break
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
        list_edges = []
        for primitive in self.primitives:
            auto_nb_pts = min(number_points, max(2, int(primitive.length() / 1e-6)))
            points = primitive.discretization_points(number_points=auto_nb_pts)
            for point1, point2 in zip(points[:-1], points[1:]):
                list_edges.append(edges.LineSegment2D(point1, point2))
        return list_edges

    def shares_primitives(self, contour):
        """
        Checks if two contour share primitives.

        """
        for prim1 in self.primitives:
            if contour.primitive_over_contour(prim1):
                return True
        return False

    def is_overlapping(self, contour2, intersecting_points=None):
        """
        Check if the contours are overlapping (a part of one is on the other).

        """

        if not intersecting_points:
            intersecting_points = self.intersection_points(contour2)

        if len(intersecting_points) < 2:
            return False

        vec1_2 = volmdlr.edges.LineSegment2D(intersecting_points[0],
                                             intersecting_points[1])
        middle_point = vec1_2.middle_point()
        normal = vec1_2.normal_vector()
        point1 = middle_point + normal * 0.00001
        point2 = middle_point - normal * 0.00001
        if (self.point_belongs(point1) and contour2.point_belongs(point1)) or \
                (not self.point_belongs(point1) and not contour2.point_belongs(point1)):
            return True
        if (self.point_belongs(point1) and self.point_belongs(point2)) or \
                (contour2.point_belongs(point1) and contour2.point_belongs(point2)):
            return True
        return False

    def is_adjacent(self, contour):
        """
        Check if two contour are adjacent.

        So: are sharing primitives but not superposing or none is inside the other.
        """

        if (self.is_inside(contour) or contour.is_inside(self)
                or self.is_overlapping(contour) or self.is_superposing(contour)):
            return False
        if self.is_sharing_primitives_with(contour):
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
            list_edges = [edge_1, edge_2, edge_1]
            for edge1, edge2 in zip(list_edges, list_edges[1:]):
                for point in [edge2.start, edge2.end]:
                    if edge1.point_belongs(point, 1e-6):
                        if not list_p:
                            list_p.append(point)
                        if list_p and point.point_distance(point.nearest_point(list_p)) > 1e-4:
                            list_p.append(point)
                        try:
                            # self.primitive_to_index(edge1)
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

    def delete_shared_contour_section(self, contour, abs_tol: float = 1e-6):
        """
        Delete shared primitives between two adjacent contours.

        :param contour: other contour.
        :param abs_tol: tolerance.
        :return: list of new primitives, without those shared by both contours.
        """
        new_primitives_contour1 = self.primitives[:]
        new_primitives_contour2 = contour.primitives[:]
        while True:
            for prim1 in new_primitives_contour1[:]:
                for prim2 in new_primitives_contour2[:]:
                    shared_section = prim1.get_shared_section(prim2, abs_tol)
                    if shared_section:
                        prim1_delete_shared_section = prim1.delete_shared_section(shared_section[0], abs_tol)
                        prim2_delete_shared_section = prim2.delete_shared_section(shared_section[0], abs_tol)
                        if prim1 in new_primitives_contour1:
                            new_primitives_contour1.remove(prim1)
                        if prim2 in new_primitives_contour2:
                            new_primitives_contour2.remove(prim2)
                        new_primitives_contour1.extend(prim1_delete_shared_section)
                        new_primitives_contour2.extend(prim2_delete_shared_section)
                        break
                else:
                    continue
                break
            else:
                break

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
        # TODO: rewrite this awful code!
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
        """
        Verifies if the entire primitive is over a contour.
        """
        return self.primitive_over_wire(primitive, tol)

    def primitive_section_over_contour(self, primitive, abs_tol: float = 1e-6):
        """
        Verifies if at least a small section of a primitive is over a contour, not necessarily the entire primitive.

        """
        for prim in self.primitives:
            shared_section = prim.get_shared_section(primitive, abs_tol)
            if shared_section:
                return True
        return False

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

    def to_polygon(self, angle_resolution, discretize_line: bool = False, discretize_line_direction: str = "xy"):
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
            if isinstance(primitive, volmdlr.edges.LineSegment2D):
                if not discretize_line:
                    polygon_points.append(primitive.start)
                else:
                    is_horizontal = math.isclose(primitive.start.y, primitive.end.y, abs_tol=1e-6)
                    is_vertical = math.isclose(primitive.start.x, primitive.end.x, abs_tol=1e-6)
                    should_discretize = discretize_line_direction == "xy" or \
                        (discretize_line_direction == "x" and is_horizontal) or \
                        (discretize_line_direction == "y" and is_vertical)
                    if should_discretize:
                        polygon_points.extend(primitive.discretization_points(angle_resolution=angle_resolution)[:-1])
                    else:
                        polygon_points.append(primitive.start)

            else:
                polygon_points.extend(primitive.discretization_points(angle_resolution=angle_resolution)[:-1])

        if isinstance(self, Contour2D):
            return ClosedPolygon2D(polygon_points)
        return ClosedPolygon3D(polygon_points)

    def invert(self):
        """Invert the Contour."""
        return self.__class__(self.inverted_primitives())

    @classmethod
    def from_points(cls, points, name: str = ''):
        """
        Create a contour from points with line_segments.
        """

        if len(points) < 3:
            raise ValueError('contour is defined at least with three points')

        linesegment_name = 'LineSegment' + points[0].__class__.__name__[-2:]
        list_edges = []
        for i in range(0, len(points) - 1):
            if points[i].is_close(points[i + 1]):
                continue
            list_edges.append(getattr(edges, linesegment_name)(points[i], points[i + 1]))
        if not points[-1].is_close(points[0]):
            list_edges.append(getattr(edges, linesegment_name)(points[-1], points[0]))

        contour = cls(list_edges, name=name)
        return contour

    def reorder_contour_at_point(self, point):
        """
        Create a new contour from self, but starting at given point.

        :param point: other point.
        :return: new contour
        """
        new_primitives_order = []
        for i, primitive in enumerate(self.primitives):
            if primitive.start.is_close(point, 1e-6):
                if i == 0:
                    return self
                new_primitives_order = self.primitives[i:] + self.primitives[:i]
                break
        new_contour = self.__class__(new_primitives_order)
        return new_contour

    def are_extremity_points_touching(self, wire):
        """
        Verifies if the extremities points of wire are touching contour.

        :param wire: other wire.
        :return: True if other contour is touching
        """
        return self.point_over_contour(wire.primitives[0].start) and self.point_over_contour(wire.primitives[-1].end)

    def is_contour_closed(self):
        return self.primitives[0].start.is_close(self.primitives[-1].end)


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

    def copy(self, deep=True, memo=None):
        """
        A specified copy of a Contour2D.
        """
        return self.__class__(primitives=[p.copy(deep, memo) for p in self.primitives],
                              name=self.name)

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
        """
        Returns the edge polygon of a contour.

        An edge polygon is the polygon generated by start and end points of each primitive of the contour.
        """
        if self._edge_polygon is None:
            self._edge_polygon = self._get_edge_polygon()
        return self._edge_polygon

    def _get_edge_polygon(self):
        """Helper function to get the edge polygon."""
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

    def point_belongs(self, point, include_edge_points: bool = False, tol: float = 1e-6):
        """
        Verifies if point belongs is within the contour.

        :param point: point to be verified.
        :param include_edge_points: consider bounds of contour or not.
        :param tol: tolerance to be considered.
        :return: True if point belongs, false otherwise.
        """
        # TODO: This is incomplete!!!
        x_min, x_max, y_min, y_max = self.bounding_rectangle
        if point.x < x_min - tol or point.x > x_max + tol or point.y < y_min - tol or point.y > y_max + tol:
            return False
        if include_edge_points:
            for primitive in self.primitives:
                if primitive.point_belongs(point, 1e-6):
                    return True
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
        """Returns the area of the contour."""
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
        """Returns the second moment of are of the contour."""
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

    def plot_data(self, edge_style: plot_data.EdgeStyle = None, surface_style: plot_data.SurfaceStyle = None):
        """Plot 2D representation of the contour."""
        plot_data_primitives = [item.plot_data() for item in self.primitives]
        return plot_data.Contour2D(plot_data_primitives=plot_data_primitives,
                                   edge_style=edge_style,
                                   surface_style=surface_style,
                                   name=self.name)

    def is_edge_inside(self, edge):
        """
        Verifies if given edge is inside self contour perimeter, including its edges.

        :param edge: other edge to verify if inside contour.
        :returns: True or False.
        """
        for point in edge.discretization_points(number_points=5):
            if not self.point_belongs(point, include_edge_points=True):
                return False
        return True

    def is_inside(self, other_contour):
        """
        Verifies if given contour is inside self contour perimeter, including its edges.

        :param other_contour: other contour.
        :returns: True or False
        """
        if other_contour.area() > self.area() and not math.isclose(other_contour.area(), self.area(), rel_tol=0.01):
            return False
        for edge in other_contour.primitives:
            if not self.is_edge_inside(edge):
                return False
        return True

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
            for i_transitions in remaining_transitions:
                i1 = sorted_inter_index_dict[2 * i_transitions + n]
                i2 = sorted_inter_index_dict[2 * i_transitions + 1 + n]
                net = abs(i2 - i1) - 1
                if net > nb_max_enclosed_transitions:
                    nb_max_enclosed_transitions = net
                    best_transition = i_transitions
                    if i1 < i2:
                        enclosed_transitions[i_transitions] = [(i + abs(n - 1)) // 2 for i
                                                               in sorted_inter_index[
                                                       i2 - 1:i1:-2]]
                    else:
                        enclosed_transitions[i_transitions] = [(i + abs(n - 1)) // 2 for i
                                                               in sorted_inter_index[
                                                       i2 + 1:i1:2]]

            remaining_transitions.remove(best_transition)
            point_start, _ = intersections[2 * best_transition + n]
            point2, _ = intersections[2 * best_transition + 1 + n]
            primitives = self.extract_with_points(point_start, point2, inside=not n)
            last_point = point2
            for transition in enclosed_transitions[best_transition]:
                point1, _ = intersections[2 * transition + n]
                point2, _ = intersections[2 * transition + 1 + n]
                primitives.append(
                    volmdlr.edges.LineSegment2D(last_point, point1))
                primitives.extend(
                    self.extract_with_points(point1, point2, inside=not n))
                last_point = point2
                if transition in remaining_transitions:
                    remaining_transitions.remove(transition)

            primitives.append(
                volmdlr.edges.LineSegment2D(last_point, point_start))

            contour = Contour2D(primitives)
            contour.order_contour()
            contours.append(contour)
        return contours

    def cut_by_line(self, line: curves.Line2D) -> List['Contour2D']:
        """
        :param line: The line used to cut the contour.

        :return: A list of resulting contours
        """
        intersections = self.line_crossings(line)
        if not intersections or len(intersections) < 2:
            return [self]
        points_intersections = [point for point, prim in intersections]
        sorted_points = line.sort_points_along_curve(points_intersections)
        list_contours = []
        contour_to_cut = self

        for point1, point2 in zip(sorted_points[:-1], sorted_points[1:]):
            closing_line = volmdlr.edges.LineSegment2D(point1, point2)
            if not contour_to_cut.point_belongs(closing_line.middle_point()):
                continue
            closing_contour = Contour2D([closing_line])
            contour1, contour2 = contour_to_cut.get_divided_contours(point1, point2, closing_contour)
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

    def split_by_line(self, line: curves.Line2D) -> List['Contour2D']:
        """Split the contour with the given line."""
        intersections = self.line_crossings(line)
        intersections = [point for point, prim in intersections]
        if not intersections:
            return [self]
        if len(intersections) < 2:
            extracted_outerpoints_contour1 = \
                Contour2D.extract(self, self.primitives[0].start, intersections[0], True)[0]
            extracted_innerpoints_contour1 = \
                Contour2D.extract(self, intersections[0], self.primitives[-1].end, True)[0]
            return extracted_outerpoints_contour1, extracted_innerpoints_contour1
        if len(intersections) == 2:
            extracted_outerpoints_contour1 = \
                Contour2D.extract(self, intersections[0], intersections[1], True)[0]
            extracted_innerpoints_contour1 = \
                Contour2D.extract(self, intersections[0], intersections[1], False)[0]
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
            cut_line = curves.Line2D(volmdlr.Point2D(xi, 0),
                                     volmdlr.Point2D(xi, 1))

            iteration_contours2 = []
            for contour in iteration_contours:
                split_contours = contour.cut_by_line(cut_line)
                if len(split_contours) == 1:
                    cutted_contours.append(contour)
                else:
                    iteration_contours2.extend(split_contours)

            iteration_contours = iteration_contours2[:]
        cutted_contours.extend(iteration_contours)
        return cutted_contours

    def triangulation(self):
        """Returns the triangulation of the contour 2d."""
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
            number_points_x = max(math.ceil(x_density * dx), min_points_x)
        if number_points_y is None:
            number_points_y = max(math.ceil(y_density * dy), min_points_y)
        x = [bounding_rectangle[0] + i * dx / number_points_x for i in range(number_points_x + 1)]
        y = [bounding_rectangle[2] + i * dy / number_points_y for i in range(number_points_y + 1)]

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

        for i in range(number_points_x):
            for j in range(number_points_y):
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

    def intersection_points(self, contour2d):
        """Returns the intersections points with other specified contour."""
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

    def get_divided_contours(self, cutting_point1: volmdlr.Point2D, cutting_point2: volmdlr.Point2D,
                             closing_contour):
        """Get divided contours."""
        extracted_innerpoints_contour1_prims, extracted_outerpoints_contour1_prims = self.split_with_two_points(
            cutting_point1, cutting_point2)
        extracted_outerpoints_contour1 = Contour2D(extracted_outerpoints_contour1_prims)
        extracted_innerpoints_contour1 = Contour2D(extracted_innerpoints_contour1_prims)
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

    def divide(self, contours):
        """Divide contour with other contours."""
        new_base_contours = [self]
        finished = False
        list_cutting_contours = contours[:]
        list_valid_contours = []
        while not finished:
            if not new_base_contours:
                break
            list_cutting_contours_modified = False
            for i, base_contour in enumerate(new_base_contours):
                for j, cutting_contour in enumerate(list_cutting_contours):
                    if base_contour.is_superposing(cutting_contour):
                        list_cutting_contours.pop(j)
                        list_cutting_contours_modified = True
                        break
                    contour_crossings = cutting_contour.wire_crossings(base_contour)
                    if contour_crossings:
                        sorted_points = cutting_contour.sort_points_along_wire(contour_crossings)
                        split_wires = cutting_contour.split_with_sorted_points(sorted_points)
                        list_cutting_contours.pop(j)
                        list_cutting_contours.extend(split_wires)
                        list_cutting_contours_modified = True
                        break
                    point1, point2 = [cutting_contour.primitives[0].start,
                                      cutting_contour.primitives[-1].end]
                    cutting_points = []
                    if base_contour.point_belongs(cutting_contour.middle_point()) and\
                            base_contour.point_over_contour(point1) and base_contour.point_over_contour(point2):
                        cutting_points = [point1, point2]
                    if cutting_points:
                        contour1, contour2 = base_contour.get_divided_contours(
                            cutting_points[0], cutting_points[1], cutting_contour)
                        new_base_contours.pop(i)
                        new_base_contours.extend([contour1, contour2])
                        break
                else:
                    list_valid_contours.append(base_contour)
                    new_base_contours.pop(i)
                    break
                if list_cutting_contours_modified:
                    break
                break

        return list_valid_contours

    def discretized_contour(self, n: float):
        """
        Discretize each contour's primitive and return a new contour with these discretized primitives.
        """
        contour = Contour2D((self.discretized_primitives(n)))

        return contour.order_contour()

    @classmethod
    def from_bounding_rectangle(cls, x_min, x_max, y_min, y_max, name: str = ''):
        """
        Create a contour 2d with bounding_box parameters, using line segments 2d.

        """

        edge0 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(x_min, y_min), volmdlr.Point2D(x_max, y_min))
        edge1 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(x_max, y_min), volmdlr.Point2D(x_max, y_max))
        edge2 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(x_max, y_max), volmdlr.Point2D(x_min, y_max))
        edge3 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(x_min, y_max), volmdlr.Point2D(x_min, y_min))

        return Contour2D([edge0, edge1, edge2, edge3], name=name)

    def cut_by_bspline_curve(self, bspline_curve2d: volmdlr.edges.BSplineCurve2D):
        """
        Cut a contour 2d with bspline_curve 2d to define two different contours.

        """
        # TODO: BsplineCurve is discretized and defined with a wire. To be improved!

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

    def merge_with(self, contour2d, abs_tol: float = 1e-6):
        """
        Merge two adjacent contours, and returns one outer contour and inner contours (if there are any).

        :param contour2d: contour to merge with.
        :param abs_tol: tolerance.
        :return: merged contours.
        """
        is_sharing_primitive = self.is_sharing_primitives_with(contour2d)
        if not is_sharing_primitive:
            if self.is_inside(contour2d):
                return [self]
            if contour2d.is_inside(self):
                return [contour2d]
            return [self, contour2d]

        merged_primitives = self.delete_shared_contour_section(contour2d, abs_tol)
        contours = Contour2D.contours_from_edges(merged_primitives, abs_tol)
        contours = sorted(contours, key=lambda contour: contour.area(),
                          reverse=True)
        return contours

    def union(self, contour2: 'Contour2D'):
        """
        Union two contours, if they are adjacent, or overlap somehow.

        """
        if self.is_inside(contour2):
            return [self]
        if contour2.is_inside(self):
            return [contour2]
        contours_intersections = self.intersection_points(contour2)
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
        merged_contours = self.merge_with(contour2)[::-1]
        merged_contours = sorted(merged_contours, key=lambda contour: contour.area(),
                                 reverse=True)
        return merged_contours

    def cut_by_wire(self, wire: Wire2D):
        """
        Cut a contour 2d with a wire 2d and return a list of contours 2d.

        :param wire: volmdlr.wires.Wire2D
        :rtype: list[volmdlr.wires.Contour2D]

        :param wire: volmdlr.wires.Wire2D.
        :return: contours2d : list[volmdlr.wires.Contour2D].
        """

        points_intersections = self.wire_intersections(wire)
        if len(points_intersections) < 2:
            return [self]
        if len(points_intersections) % 2 != 0:
            raise NotImplementedError(
                f'{len(points_intersections)} intersections not supported yet')
        sorted_points = wire.sort_points_along_wire(points_intersections)
        split_wires = wire.split_with_sorted_points(sorted_points)
        valid_cutting_wires = []
        for split_wire in split_wires:
            if self.is_superposing(split_wire) or not self.is_inside(split_wire):
                continue
            valid_cutting_wires.append(split_wire)
        divided_contours = self.divide(valid_cutting_wires)
        return divided_contours

    def intersection_contour_with(self, other_contour, abs_tol=1e-6):
        """
        Gets the contour(s) resulting from the intersections of two other contours.

        :param other_contour: other contour.
        :param abs_tol: tolerance.
        :return: list of resulting intersection contours.
        """
        contour_crossings = self.wire_crossings(other_contour)
        sorted_points_contour1 = sorted(contour_crossings, key=self.abscissa)
        sorted_points_contour2 = sorted(contour_crossings, key=other_contour.abscissa)
        split_wires1 = self.split_with_sorted_points(sorted_points_contour1)
        split_wires2 = other_contour.split_with_sorted_points(sorted_points_contour2)
        intersection_contour_primitives = []
        for section in split_wires1:
            if other_contour.is_inside(section):
                intersection_contour_primitives.extend(section.primitives)
        for section in split_wires2:
            if self.is_inside(section):
                intersection_contour_primitives.extend(section.primitives)
        return self.contours_from_edges(intersection_contour_primitives, abs_tol)

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
        """Returns polygon statistics."""
        distances = []
        for i, point in enumerate(self.points):
            if i != 0:
                distances.append(point.point_distance(self.points[i - 1]))
        mean_distance = mean(distances)
        std = npy.std(distances)
        return mean_distance, std

    def simplify_polygon(self, min_distance: float = 0.01, max_distance: float = 0.05, angle: float = 15):
        """Simplify polygon."""
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

    def invert(self):
        """Invert the polygon."""
        return self.__class__(self.points[::-1])

    @property
    def line_segments(self):
        """Polygon line segments."""
        if not self._line_segments:
            self._line_segments = self.get_line_segments()
        return self._line_segments

    def get_line_segments(self):
        """Abstract method."""
        raise NotImplementedError(
            f"get_line_segments method must be overloaded by {self.__class__.__name__}")


class ClosedPolygon2D(ClosedPolygonMixin, Contour2D):
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
        """Returns a copy of the object."""
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
        """Returns the area of the polygon."""
        # TODO: performance: cache number of points
        if len(self.points) < 3:
            return 0.

        x = [point.x for point in self.points]
        y = [point.y for point in self.points]

        x1 = [x[-1]] + x[0:-1]
        y1 = [y[-1]] + y[0:-1]
        return 0.5 * abs(sum(i * j for i, j in zip(x, y1))
                         - sum(i * j for i, j in zip(y, x1)))

    def center_of_mass(self):
        """Returns polygon's center of mass."""
        lngth_points = len(self.points)
        if lngth_points == 0:
            return volmdlr.O2D
        if lngth_points == 1:
            return self.points[0]
        if lngth_points == 2:
            return 0.5 * (self.points[0] + self.points[1])

        x = [point.x for point in self.points]
        y = [point.y for point in self.points]

        xi_xi1 = x + npy.roll(x, -1)
        yi_yi1 = y + npy.roll(y, -1)
        xi_yi1 = npy.multiply(x, npy.roll(y, -1))
        xi1_yi = npy.multiply(npy.roll(x, -1), y)

        signed_area = 0.5 * npy.sum(xi_yi1 - xi1_yi)  # signed area!
        if not math.isclose(signed_area, 0, abs_tol=1e-09):
            center_x = npy.sum(npy.multiply(xi_xi1, (xi_yi1 - xi1_yi))) / 6. / signed_area
            center_y = npy.sum(npy.multiply(yi_yi1, (xi_yi1 - xi1_yi))) / 6. / signed_area
            return volmdlr.Point2D(center_x, center_y)

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

    def point_belongs(self, point, include_edge_points: bool = False, tol: float = 1e-6):
        """
        Ray casting algorithm copied from internet.
        """
        return polygon_point_belongs(npy.array([(point_.x, point_.y) for point_ in self.points], dtype=npy.float64),
                                     npy.array([point.x, point.y], dtype=npy.float64),
                                     include_edge_points=include_edge_points, tol=tol)

    def points_in_polygon(self, points, include_edge_points: bool = False, tol: float = 1e-6):
        """
        Check if a list of points is inside the polygon using parallel computing.

        :param points: List of points in the form [(x1, y1), (x2, y2), ...]
        :type points: list or numpy.ndarray
        :param include_edge_points: Flag to include edge points as inside the polygon
        :type include_edge_points: bool, optional
        :param tol: Tolerance for numerical comparisons, defaults to 1e-6
        :type tol: float, optional
        :return: List of boolean values indicating whether each point is inside the polygon
        :rtype: numpy.ndarray
        """
        if isinstance(points, list):
            points = npy.array([(point[0], point[1]) for point in points], dtype=npy.float64)
        polygon = npy.array([(point_[0], point_[1]) for point_ in self.points], dtype=npy.float64)
        return points_in_polygon(polygon, points, include_edge_points=include_edge_points, tol=tol)

    def second_moment_area(self, point):
        """Returns the second moment of area of the polygon."""
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
        """Get polygon lines."""
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
        """Apply transformation to the object."""
        return self.__class__([point.frame_mapping(frame, side) for point in self.points])

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        warnings.warn("'in-place' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for point in self.points:
            point.frame_mapping_inplace(frame, side)

    def polygon_distance(self, polygon: 'ClosedPolygon2D') -> float:
        """Returns the minimum distance to other given polygon."""
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
            vector1 = self.points[i + 1] - self.points[i]
            vector2 = self.points[i] - self.points[i + 1]
            vector1 = vector1.unit_vector()
            vector2 = vector2.unit_vector()
            vectors.append(vector1)
            vectors.append(vector2)

        vector1 = self.points[0] - self.points[-1]
        vector2 = self.points[-1] - self.points[0]
        vector1 = vector1.unit_vector()
        vector2 = vector2.unit_vector()
        vectors.append(vector1)
        vectors.append(vector2)

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
                vector_i = vector_i.unit_vector()
                if vector_i.dot(vectors[2 * i - 1].normal_vector()) > 0:
                    vector_i = - vector_i
                    # check = True
                offset_vectors.append(vector_i)

            normal_vector1 = - vectors[2 * i - 1].normal_vector()
            normal_vector2 = vectors[2 * i].normal_vector()
            normal_vector1 = normal_vector1.unit_vector()
            normal_vector2 = normal_vector2.unit_vector()
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

        :return: True if the polygon self intersects, False otherwise. If True, returns two
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
        number = len(sorted_index)
        segments = []

        for i, index in enumerate(sorted_index):
            # Stock the segments between 2 consecutive edges
            # Ex: for the ABCDE polygon, if Sweep Line is on C, the segments
            #   will be (C,B) and (C,D)
            if index - 1 < 0:
                segments.append((index, number - 1))
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
    def points_convex_hull(cls, points, name: str = ''):
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

        return cls(hull, name=name)

    @classmethod
    def concave_hull(cls, points, concavity, scale_factor, name: str = ''):
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
        :type scale_factor: float.
        :param name: object's name.

        """

        def get_nearby_points(line, points, scale_factor):
            points_hull = [point.copy() for point in points]

            nearby_points = []
            line_midpoint = 0.5 * (line.start + line.end)
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
                        if boundary[1] >= point_x_rel_pos >= boundary[0] <= point_y_rel_pos <= boundary[3]:
                            nearby_points.append(point)

                scale_factor *= 4 / 3
                tries += 1

            return nearby_points

        def line_colides_with_hull(line, concave_hull):
            for hull_line in concave_hull:
                if not line.start.is_close(hull_line.start) and not line.start.is_close(hull_line.end) and \
                        not line.end.is_close(hull_line.start) and not line.end.is_close(hull_line.end):
                    if line.line_intersections(hull_line.line):
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
                #  We want the middle-point to be the one with the widest angle (smallest cosine)
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

        return cls(points, name=name)  # , nearby_points

    @classmethod
    def convex_hull_points(cls, points, name: str = ''):
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
                if point1.is_close(points_hull[-1]):
                    points_hull.append(point2)
                    break
                if point2.is_close(points_hull[-1]):
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

        return cls(points_hull, name=name)

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
            for index_point, point in enumerate(self.points):
                ax.text(*point, f'point {index_point + 1}', ha='center', va='top')

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
        if len(tri['vertices']) < 3:
            return None
        triangulate_result = triangulate(tri, tri_opt)
        triangles = triangulate_result['triangles'].tolist()
        number_points = triangulate_result['vertices'].shape[0]
        points = [vmd.Node2D(*triangulate_result['vertices'][i, :]) for i in range(number_points)]
        return vmd.DisplayMesh2D(points, triangles=triangles)

    def grid_triangulation_points(self, number_points_x: int = 25, number_points_y: int = 25,
                                  include_edge_points: bool = True):
        """
        Use an n by m grid to triangulate the contour.

        :param number_points_x: Number of discretization points in x direction.
        :type number_points_x: int
        :param number_points_y: Number of discretization points in y direction.
        :type number_points_y: int
        :param include_edge_points: Flag to include edge points as inside the polygon
        :type include_edge_points: bool, optional
        :return: Discretization data.
        :rtype: list
        """
        x_min, x_max, y_min, y_max = self.bounding_rectangle.bounds()

        x = npy.linspace(x_min, x_max, num=number_points_x + 2)
        y = npy.linspace(y_min, y_max, num=number_points_y + 2)

        grid_point_index = {}

        polygon_points = {vmd.Node2D.from_point(point) for point in self.points}

        # Generate all points in the grid
        grid_points = npy.array([(xi, yi) for xi in x for yi in y], dtype=npy.float64)

        # Use self.points_in_polygon to check if each point is inside the polygon
        points_in_polygon_ = self.points_in_polygon(grid_points, include_edge_points=include_edge_points)

        # Find the indices where points_in_polygon is True (i.e., points inside the polygon)
        indices = npy.where(points_in_polygon_)[0]

        points = []

        for i in indices:
            point = vmd.Node2D(*grid_points[i])
            if point not in polygon_points:
                grid_point_index[(i // (number_points_y + 2), i % (number_points_y + 2))] = point
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
        """Simplify polygon."""
        return ClosedPolygon2D(self.simplify_polygon(min_distance=min_distance,
                                                     max_distance=max_distance).points)

    def line_intersecting_closing_point(self, crossing_point):
        """
        Finds closing point for the sewing method using intersection of lines drawn from the barycenter.

        returns the closing point.
        """
        vec_dir = crossing_point.copy()
        vec_dir = vec_dir.unit_vector()

        line = volmdlr.edges.LineSegment2D(volmdlr.O2D,
                                           crossing_point + vec_dir * 5)
        # line.plot(ax=ax2d, color='b')

        point_intersections = {}
        for line_segment in self.line_segments:
            point_intersection = line_segment.linesegment_intersections(
                line)
            if point_intersection:
                point_intersections[line_segment] = point_intersection[0]
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
            new_polygon = self.translation(-barycenter1_2d)
            way_back = barycenter1_2d.to_3d(volmdlr.O3D, x, y)
        else:
            inters = self.linesegment_intersections(linex)
            distance = inters[0][0].point_distance(inters[-1][0])
            if distance / 2 > 3 * min(
                    self.point_distance(inters[0][0]),
                    self.point_distance(inters[-1][0])):
                mid_point = (inters[0][0] + inters[-1][0]) * 0.5
                new_polygon = self.translation(-mid_point)
                way_back = mid_point.to_3d(volmdlr.O3D, x, y)

        return new_polygon, way_back

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
                                                        (middle_point - normal_vector).to_point())
            line_segment2 = volmdlr.edges.LineSegment2D(middle_point,
                                                        (middle_point + normal_vector).to_point())
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


class Triangle2D(Triangle, ClosedPolygon2D):
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

        Triangle.__init__(self, point1, point2, point3, name)

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


class Contour3D(ContourMixin, Wire3D):
    """
    A collection of 3D primitives forming a closed wire3D.

    """
    _non_serializable_attributes = ['points']
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['points', 'name']
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
        step_id = kwargs.get("step_id")
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
        contour = cls(raw_edges, name=name)
        if contour.is_ordered():
            return contour
        list_contours = cls.contours_from_edges(raw_edges.copy())
        for contour_reordered in list_contours:
            if contour_reordered.is_ordered():
                return contour_reordered
        list_edges = reorder_contour3d_edges_from_step(raw_edges, [step_id, step_name, arguments])
        if list_edges:
            return cls(list_edges, name=name)
        return None

    def to_step(self, current_id, surface_id=None, surface3d=None):
        """
        Converts the object to a STEP representation.

        :param current_id: The ID of the last written primitive.
        :type current_id: int
        :return: The STEP representation of the object and the last ID.
        :rtype: tuple[str, list[int]]
        """
        content = ''
        edge_ids = []
        for primitive in self.primitives:
            primitive_content, primitive_id = primitive.to_step(current_id, surface_id=surface_id)

            content += primitive_content
            current_id = primitive_id + 1

            content += f"#{current_id} = ORIENTED_EDGE('{primitive.name}',*,*,#{primitive_id},.T.);\n"
            edge_ids.append(current_id)

            current_id += 1

        content += f"#{current_id} = EDGE_LOOP('{self.name}',({volmdlr.core.step_ids_to_str(edge_ids)}));\n"
        return content, current_id

    def average_center_point(self):
        """
        Returns the average center point.
        """
        number_points = len(self.edge_polygon.points)
        x = sum(point[0] for point in self.edge_polygon.points) / number_points
        y = sum(point[1] for point in self.edge_polygon.points) / number_points
        z = sum(point[2] for point in self.edge_polygon.points) / number_points

        return volmdlr.Point3D(x, y, z)

    def to_2d(self, plane_origin: volmdlr.Point3D, x: volmdlr.Vector3D, y: volmdlr.Vector3D) -> Contour2D:
        """
        Converts 3D contour into a 2D contour.

        :param plane_origin: 3D point representing the origin of the coordinates' system.
        :param x: 3D vector representing the x direction of the coordinates' system.
        :param y: 3D vector representing the y direction of the coordinates' system.
        :return: Equivalent 2D contour.
        """
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

    def line_intersections(self, line: curves.Line3D):
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

    def clean_primitives(self):
        """
        Delete primitives with start=end, and return a new contour.

        """

        new_primitives = []
        for primitive in self.primitives:
            if not primitive.start.is_close(primitive.end):
                new_primitives.append(primitive)

        return Contour3D(new_primitives)

    def merge_with(self, contour3d, abs_tol: float = 1e-6):
        """
        Merge two adjacent contours, and returns one outer contour and inner contours (if there are any).

        """

        merged_primitives = self.delete_shared_contour_section(contour3d, abs_tol)
        contours = Contour3D.contours_from_edges(merged_primitives, tol=abs_tol)

        return contours


class ClosedPolygon3D(Contour3D, ClosedPolygonMixin):
    """
    A collection of points, connected by line segments, following each other.

    """
    _non_serializable_attributes = ['line_segments', 'primitives']
    _non_data_eq_attributes = ['line_segments', 'primitives']

    def __init__(self, points: List[volmdlr.Point3D], name: str = ''):
        self.points = points
        self._line_segments = None

        Contour3D.__init__(self, self.line_segments, name)

    def get_line_segments(self):
        """Get polygon lines."""
        lines = []
        if len(self.points) > 1:
            for point1, point2 in zip(self.points,
                                      list(self.points[1:]) + [self.points[0]]):
                if not point1.is_close(point2):
                    lines.append(volmdlr.edges.LineSegment3D(point1, point2))
        return lines

    def copy(self, *args, **kwargs):
        """Returns a copy of the object."""
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

    @staticmethod
    def fix_sewing_normals(triangles, reference_linesegment):
        """ Fixes sewing triangle normal so it faces always outwards."""
        first_triangles_points = triangles[0]
        frame = volmdlr.Frame3D.from_3_points(*first_triangles_points)
        normal = frame.w
        middle_point = (first_triangles_points[0] + first_triangles_points[1] + first_triangles_points[2]) / 3
        point1 = middle_point + 0.05 * normal
        point2 = middle_point - 0.05 * normal
        if reference_linesegment.line.point_distance(point1) < reference_linesegment.line.point_distance(point2):
            return [points[::-1] for points in triangles]
        return triangles

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
        reference_linesegment = edges.LineSegment3D(center1, center2)
        triangles = self.fix_sewing_normals(triangles, reference_linesegment)
        closing_triangles = polygon2.close_sewing(dict_closing_pairs)
        closing_triangles = self.fix_sewing_normals(closing_triangles, reference_linesegment)
        triangles += closing_triangles

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
                        triangles_points.append(face_points[::-1])
                    elif index[0] > index[1]:
                        if (i - 1 <= index[0] and i <= index[1]) or (
                                (i - 1 >= index[0]) and i >= index[1]):
                            face_points = [self.points[i - 1],
                                           point_polygon2,
                                           list(dict_closing_pairs.keys())[j]]
                            triangles_points.append(face_points[::-1])
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
        """
        Chooses sewing closing point when point index passes through zero index again.
        """
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
        """Validates concave closing point."""
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
        polygon1_2d = self.to_2d(volmdlr.O3D, x, y)
        polygon2_2d = polygon2.to_2d(volmdlr.O3D, x, y)
        polygon1_3d = self
        polygon2_3d = polygon2
        need_fix_normal = False
        if polygon2_2d.area() < polygon1_2d.area():
            polygon1_2d, polygon2_2d = polygon2_2d, polygon1_2d
            polygon1_3d = polygon2
            polygon2_3d = self
            need_fix_normal = True
        polygon1_3d = polygon1_3d.get_valid_concave_sewing_polygon(
            polygon1_2d, polygon2_2d)
        polygon1_2d = polygon1_3d.to_2d(volmdlr.O3D, x, y)

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
                    triangles_points = polygon2_3d.redefine_sewing_triangles_points(triangles_points,
                                                                                    passed_by_zero_index,
                                                                                    closing_point_index,
                                                                                    previous_closing_point_index)
                    if dict_closing_pairs:
                        dict_closing_pairs, lower_bounddary_closing_point = \
                            self.clean_sewing_closing_pairs_dictionary(dict_closing_pairs,
                                                                       closing_point_index,
                                                                       passed_by_zero_index)

                        if len(new_list_remove_closing_indexes) < len(new_list_closing_point_indexes):
                            dict_closing_pairs[lower_bounddary_closing_point] = (
                                new_list_closing_point_indexes[-(len(new_list_remove_closing_indexes) + 1)],
                                closing_point_index)
                    for pt_index in list_remove_closing_points:
                        list_closing_point_indexes.remove(pt_index)
                    list_closing_point_indexes.append(closing_point_index)

                elif (not passed_by_zero_index and closing_point_index > polygon2_3d.points.index(
                            triangles_points[-len(list_remove_closing_points) - 1][2])) or \
                        (passed_by_zero_index and closing_point_index >= 0):
                    triangles_points = polygon2_3d.redefine_sewing_triangles_points(triangles_points,
                                                                                    passed_by_zero_index,
                                                                                    closing_point_index,
                                                                                    previous_closing_point_index)
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
        if need_fix_normal:
            center1, center2 = self.average_center_point(), polygon2.average_center_point()
            reference_segment = edges.LineSegment3D(center1, center2)
            triangles_points = self.fix_sewing_normals(triangles_points, reference_segment)
        return triangles_points

    def sewing(self, polygon2, x, y):
        """
        Sew two polygon3D together.

        :param x: The vector representing first direction to project polygons in
        :param y: The vector representing second direction to project polygons in
        """
        polygon1_2d = self.to_2d(volmdlr.O3D, x, y)
        polygon2_2d = polygon2.to_2d(volmdlr.O3D, x, y)
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
