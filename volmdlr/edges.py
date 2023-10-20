#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Edges related classes.
"""
import math
import sys
import warnings
from itertools import product
from typing import List, Union

import dessia_common.core as dc
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as npy
import plot_data.core as plot_data
import plot_data.colors
import scipy.integrate as scipy_integrate
from scipy.optimize import least_squares
from geomdl import NURBS, BSpline

from volmdlr.nurbs.operations import split_curve
from volmdlr.nurbs.core import evaluate_curve, derivatives_curve
from volmdlr.nurbs import fitting
import volmdlr.nurbs.helpers as nurbs_helpers

import volmdlr.core
import volmdlr.core_compiled
import volmdlr.geometry
from volmdlr import curves as volmdlr_curves
import volmdlr.utils.common_operations as vm_common_operations
import volmdlr.utils.intersections as vm_utils_intersections
from volmdlr.core import EdgeStyle


class Edge(dc.DessiaObject):
    """
    Defines a simple edge Object.
    """

    def __init__(self, start, end, name=''):
        self.start = start
        self.end = end
        self._length = None
        self._direction_vector_memo = None
        self._unit_direction_vector_memo = None
        self._reverse = None
        self._middle_point = None
        # Disabling super init call for performance
        # dc.DessiaObject.__init__(self, name=name)
        self.name = name

    def __getitem__(self, key):
        if key == 0:
            return self.start
        if key == 1:
            return self.end
        raise IndexError

    @property
    def periodic(self):
        """Returns True if edge is periodic."""
        return False

    def is_close(self, other_edge, tol: float = 1e-6):
        """
        Verify if two edges are equal, considering a certain tolerance.

        """
        raise NotImplementedError(f'is_close method not implemented by {self.__class__.__name__}')

    def get_reverse(self):
        """
        Gets the same edge, but in the opposite direction.

        """
        raise NotImplementedError(f'get_reverse method not implemented by {self.__class__.__name__}')

    def split(self, split_point, tol: float = 1e-6):
        """
        Gets the same edge, but in the opposite direction.

        """
        raise NotImplementedError(f'split method not implemented by {self.__class__.__name__}')

    def reverse(self):
        """Gets the edge in the reverse direction."""
        if self._reverse is None:
            self._reverse = self.get_reverse()
        return self._reverse

    def direction_independent_is_close(self, other_edge, tol: float = 1e-6):
        """
        Verifies if two line segments are the same, not considering its direction.

        """
        if not isinstance(self, other_edge.__class__):
            return False
        if self.is_close(other_edge, tol):
            return True
        return self.reverse().is_close(other_edge, tol)

    def length(self):
        """
        Calculates the edge's length.
        """
        raise NotImplementedError(f'length method not implemented by {self.__class__.__name__}')

    def point_at_abscissa(self, abscissa):
        """
        Calculates the point at given abscissa.

        """
        raise NotImplementedError(f'point_at_abscissa method not implemented by {self.__class__.__name__}')

    def middle_point(self):
        """
        Gets the middle point for an edge.

        :return:
        """
        if not self._middle_point:
            half_length = self.length() / 2
            self._middle_point = self.point_at_abscissa(abscissa=half_length)
        return self._middle_point

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = None):
        """
        Discretize an Edge to have "n" points.

        :param number_points: the number of points (including start and end
            points) if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have
            a controlled angular distance. Useful to mesh an arc
        :return: a list of sampled points
        """
        if angle_resolution:
            number_points = int(angle_resolution * (self.length() / math.pi))
        if number_points is None or number_points <= 1:
            number_points = 2
        step = self.length() / (number_points - 1)
        return [self.point_at_abscissa(i * step) for i in range(number_points)]

    def polygon_points(self, discretization_resolution: int):
        """
        Deprecated method of discretization_points.
        """
        warnings.warn('polygon_points is deprecated,\
        please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(number_points=discretization_resolution)

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to an Edge type object.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated
        :type object_dict: dict
        :return: The corresponding Edge object
        :rtype: :class:`volmdlr.edges.Edge`
        """
        obj = object_dict[arguments[3]]
        point1 = object_dict[arguments[1]]
        point2 = object_dict[arguments[2]]
        same_sense = bool(arguments[4] == ".T.")
        if obj.__class__.__name__ == 'LineSegment3D':
            if not point1.is_close(point2):
                return LineSegment3D(point1, point2, name=arguments[0][1:-1])
            return None
        if obj.__class__.__name__ == 'Line3D':
            if not same_sense:
                obj = obj.reverse()
            if not point1.is_close(point2):
                return LineSegment3D(point1, point2, obj, arguments[0][1:-1])
            return None

        if hasattr(obj, 'trim'):
            trimmed_edge = obj.trim(point1, point2, same_sense)
            if trimmed_edge:
                trimmed_edge.name = arguments[0][1:-1]
            return trimmed_edge

        raise NotImplementedError(f'Unsupported #{arguments[3]}: {object_dict[arguments[3]]}')

    def normal_vector(self, abscissa):
        """
        Calculates the normal vector the edge at given abscissa.

        :return: the normal vector
        """
        raise NotImplementedError('the normal_vector method must be'
                                  'overloaded by subclassing class')

    def unit_normal_vector(self, abscissa: float = 0.0):
        """
        Calculates the unit normal vector the edge at given abscissa.

        :param abscissa: edge abscissa
        :return: unit normal vector
        """
        vector = self.normal_vector(abscissa).copy(deep=True)
        vector = vector.unit_vector()
        return vector

    def direction_vector(self, abscissa):
        """
        Calculates the direction vector the edge at given abscissa.

        :param abscissa: edge abscissa
        :return: direction vector
        """
        raise NotImplementedError('the direction_vector method must be'
                                  'overloaded by subclassing class')

    def unit_direction_vector(self, abscissa: float = 0.0):
        """
        Calculates the unit direction vector the edge at given abscissa.

        :param abscissa: edge abscissa
        :return: unit direction vector
        """
        if not self._unit_direction_vector_memo:
            self._unit_direction_vector_memo = {}
        if abscissa not in self._unit_direction_vector_memo:
            vector = self.direction_vector(abscissa).copy(deep=True)
            vector = vector.unit_vector()
            self._unit_direction_vector_memo[abscissa] = vector
        return self._unit_direction_vector_memo[abscissa]

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point: Point to be verified
        :return: Return True if the point belongs to this surface,
            or False otherwise
        """
        raise NotImplementedError(f'the straight_line_point_belongs method must be'
                                  f' overloaded by {self.__class__.__name__}')

    def point_belongs(self, point, abs_tol: float = 1e-6):
        """
        Checks if a point belongs to the edge.

        :param point: The point to be checked
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :param abs_tol: The precision in terms of distance.
            Default value is 1e-6
        :type abs_tol: float, optional
        :return: `True` if the point belongs to the edge, `False` otherwise
        :rtype: bool
        """
        raise NotImplementedError(f'the point_belongs method must be'
                                  f' overloaded by {self.__class__.__name__}')

    def touching_points(self, edge2):
        """
        Verifies if two edges are touching each other.

        In case these two edges are touching each other, return these touching points.

        :param edge2: edge2 to verify touching points.
        :return: list of touching points.
        """
        point1, point2 = edge2.start, edge2.end
        point3, point4 = self.start, self.end
        touching_points = []
        for primitive, points in zip([self, edge2], [[point1, point2], [point3, point4]]):
            for point in points:
                if point not in touching_points and primitive.point_belongs(point):
                    touching_points.append(point)
        return touching_points

    def _get_intersection_sections(self, edge2):
        """
        Identify the sections where there may exist intersection between any two edges.

        :param edge2: other edge.
        :return: list containing the sections pairs to further search for intersections.
        """
        def edge3d_section_validator(line_seg1, line_seg2):
            return line_seg1.bounding_box.bbox_intersection(line_seg2.bounding_box)

        def edge2d_section_validator(line_seg1, line_seg2):
            return line_seg1.linesegment_intersections(line_seg2)

        # min_dist, pt1, pt2 = self.minimum_distance(edge2, True)
        lineseg_class_ = getattr(sys.modules[__name__], 'LineSegment' + self.__class__.__name__[-2:])
        section_validor_ = edge2d_section_validator
        if lineseg_class_ == LineSegment3D:
            section_validor_ = edge3d_section_validator
        bspline_discretized_points1 = []
        for point in self.discretization_points(number_points=30):
            if not volmdlr.core.point_in_list(point, bspline_discretized_points1):
                bspline_discretized_points1.append(point)
        line_segments1 = [lineseg_class_(point1, point2) for point1, point2 in
                          zip(bspline_discretized_points1[:-1], bspline_discretized_points1[1:])]
        edge_discretized_points2 = []
        for point in edge2.discretization_points(number_points=30):
            if not volmdlr.core.point_in_list(point, edge_discretized_points2):
                edge_discretized_points2.append(point)
        line_segments2 = [lineseg_class_(point1, point2) for point1, point2 in
                          zip(edge_discretized_points2[:-1], edge_discretized_points2[1:])]
        intersection_section_pairs = []
        for lineseg1, lineseg2 in product(line_segments1, line_segments2):
            valid_section = section_validor_(lineseg1, lineseg2)
            if valid_section:
                intersection_section_pairs.append((self.split_between_two_points(lineseg1.start, lineseg1.end),
                                                   edge2.split_between_two_points(lineseg2.start, lineseg2.end)))
        return intersection_section_pairs

    def _generic_edge_intersections(self, edge2, abs_tol: float = 1e-6):
        """
        General method to calculate the intersection of any two edges.

        :param edge2: other edge
        :param abs_tol: tolerance.
        :return: intersections between the two edges.
        """
        intersections = []
        for edge_extremity in [self.start, self.end]:
            if edge2.point_belongs(edge_extremity):
                intersections.append(edge_extremity)
        for edge_extremity in [edge2.start, edge2.end]:
            if self.point_belongs(edge_extremity):
                intersections.append(edge_extremity)
        intersection_section_pairs = self._get_intersection_sections(edge2)
        for bspline, edge2_ in intersection_section_pairs:
            min_dist, point_min_dist_1, _ = bspline.minimum_distance(edge2_, True)
            if not math.isclose(min_dist, 0.0, abs_tol=1e-6):
                continue
            intersections_points = vm_utils_intersections.get_bsplinecurve_intersections(
                edge2_, bspline, abs_tol=abs_tol)
            if not intersections_points:
                intersections.append(point_min_dist_1)
            for intersection in intersections_points:
                if not volmdlr.core.point_in_list(intersection, intersections):
                    intersections.append(intersection)
            # intersections.extend(intersections_points)
        return intersections

    def intersections(self, edge2: 'Edge', abs_tol: float = 1e-6):
        """
        Gets the intersections between two edges.

        :param edge2: other edge.
        :param abs_tol: tolerance.
        :return: list of intersection points.
        """
        method_name = f'{edge2.__class__.__name__.lower()[:-2]}_intersections'
        if hasattr(self, method_name):
            intersections = getattr(self, method_name)(edge2, abs_tol)
            return intersections
        method_name = f'{self.__class__.__name__.lower()[:-2]}_intersections'
        if hasattr(edge2, method_name):
            intersections = getattr(edge2, method_name)(self, abs_tol)
            return intersections
        return self._generic_edge_intersections(edge2, abs_tol)

    def validate_crossings(self, edge, intersection):
        """Validates the intersections as crossings: edge not touching the other at one end, or in a tangent point."""
        if not volmdlr.core.point_in_list(intersection, [self.start, self.end, edge.start, edge.end]):
            tangent1 = self.unit_direction_vector(self.abscissa(intersection))
            tangent2 = edge.unit_direction_vector(edge.abscissa(intersection))
            if math.isclose(abs(tangent1.dot(tangent2)), 1, abs_tol=1e-6):
                return None
        else:
            return None
        return intersection

    def crossings(self, edge):
        """
        Gets the crossings between two edges.

        """
        valid_crossings = []
        intersections = self.intersections(edge)
        for intersection in intersections:
            crossing = self.validate_crossings(edge, intersection)
            if crossing:
                valid_crossings.append(crossing)
        return valid_crossings

    def abscissa(self, point, tol: float = 1e-6):
        """
        Computes the abscissa of an Edge.

        :param point: The point located on the edge.
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`].
        :param tol: The precision in terms of distance. Default value is 1e-4.
        :type tol: float, optional.
        :return: The abscissa of the point.
        :rtype: float
        """
        raise NotImplementedError(f'the abscissa method must be overloaded by {self.__class__.__name__}')

    def local_discretization(self, point1, point2, number_points: int = 10):
        """
        Gets n discretization points between two given points of the edge.

        :param point1: point 1 on edge.
        :param point2: point 2 on edge.
        :param number_points: number of points to discretize locally.
        :return: list of locally discretized points.
        """
        abscissa1 = self.abscissa(point1)
        abscissa2 = self.abscissa(point2)
        return vm_common_operations.get_abscissa_discretization(self, abscissa1, abscissa2, number_points, False)

    def split_between_two_points(self, point1, point2):
        """
        Split edge between two points.

        :param point1: point 1.
        :param point2: point 2.
        :return: edge split.
        """
        split1 = self.split(point1)
        if split1[0] and split1[0].point_belongs(point2, abs_tol=1e-6):
            split2 = split1[0].split(point2)
        else:
            split2 = split1[1].split(point2)
        new_split_edge = None
        for split_edge in split2:
            if split_edge and split_edge.point_belongs(point1, 1e-4) and split_edge.point_belongs(point2, 1e-4):
                new_split_edge = split_edge
                break
        return new_split_edge

    def point_distance_to_edge(self, point):
        """
        Calculates the distance from a given point to an edge.

        :param point: point.
        :return: distance to edge.
        """
        return vm_common_operations.get_point_distance_to_edge(self, point, self.start, self.end)

    @property
    def simplify(self):
        """Search another simplified edge that can represent the edge."""
        return self

    def is_point_edge_extremity(self, other_point, abs_tol: float = 1e-6):
        """
        Verifies if a point is the start or the end of the edge.

        :param other_point: other point to verify if it is any end of the edge.
        :param abs_tol: tolerance.
        :return: True of False.
        """
        if self.start.is_close(other_point, abs_tol):
            return True
        if self.end.is_close(other_point, abs_tol):
            return True
        return False

    def _generic_minimum_distance(self, element, return_points=False):
        """
        Gets the minimum distance between two elements.

        This is a generalized method in a case an analytical method has not yet been defined.

        :param element: other element.
        :param return_points: Weather to return the corresponding points or not.
        :return: distance to edge.
        """
        n = max(1, int(self.length() / element.length()))
        best_distance = math.inf
        distance_points = None
        distance = best_distance

        point1_edge1_ = self.start
        point2_edge1_ = self.end

        point1_edge2_ = element.start
        point2_edge2_ = element.end
        # min_dist_point1 = None
        # min_dist_point2 = None
        linesegment_class_ = getattr(sys.modules[__name__], 'LineSegment' + self.__class__.__name__[-2:])
        while True:
            edge1_discretized_points_between_1_2 = self.local_discretization(point1_edge1_, point2_edge1_,
                                                                             number_points=10 * n)
            edge2_discretized_points_between_1_2 = element.local_discretization(point1_edge2_, point2_edge2_)
            if not edge1_discretized_points_between_1_2:
                break
            distance = edge2_discretized_points_between_1_2[0].point_distance(edge1_discretized_points_between_1_2[0])
            distance_points = [edge2_discretized_points_between_1_2[0], edge1_discretized_points_between_1_2[0]]
            for point1_edge1, point2_edge1 in zip(edge1_discretized_points_between_1_2[:-1],
                                                  edge1_discretized_points_between_1_2[1:]):
                lineseg1 = linesegment_class_(point1_edge1, point2_edge1)
                for point1_edge2, point2_edge2 in zip(edge2_discretized_points_between_1_2[:-1],
                                                      edge2_discretized_points_between_1_2[1:]):
                    lineseg2 = linesegment_class_(point1_edge2, point2_edge2)
                    dist, min_dist_point1_, min_dist_point2_ = lineseg1.minimum_distance(lineseg2, True)
                    if dist < distance:
                        point1_edge1_, point2_edge1_ = point1_edge1, point2_edge1
                        point1_edge2_, point2_edge2_ = point1_edge2, point2_edge2
                        distance = dist
                        distance_points = [min_dist_point1_, min_dist_point2_]
            if math.isclose(distance, best_distance, abs_tol=1e-6):
                break
            best_distance = distance
            # best_distance_points = distance_points
            n = 1
        if return_points:
            return distance, distance_points[0], distance_points[1]
        return distance

    def minimum_distance(self, element, return_points=False):
        """
        Evaluates the minimal distance between the edge and another specified primitive.

        :param element: Another primitive object to compute the distance to.
        :param return_points: (optional) If True, return the closest points on both primitives.
        :type return_points: bool

        :return: The minimum distance between the edge and the specified primitive.
            tuple, optional: A tuple containing the closest points if return_points is True.
        """
        method_name_ = 'distance_' + element.__class__.__name__.lower()[:-2]
        if hasattr(self, method_name_):
            return getattr(self, method_name_)(element, return_points)
        method_name_ = 'distance_to_' + self.__class__.__name__.lower()[:-2]
        if hasattr(element, method_name_):
            return getattr(element, method_name_)(self, return_points)
        return self._generic_minimum_distance(element, return_points)

    def abscissa_discretization(self, abscissa1, abscissa2, max_number_points: int = 10,
                                return_abscissas: bool = True):
        """
        Gets n discretization points between two given points of the edge.

        :param abscissa1: Initial abscissa.
        :param abscissa2: Final abscissa.
        :param max_number_points: Expected number of points to discretize locally.
        :param return_abscissas: By default, returns also a list of abscissas corresponding to the
            discretization points
        :return: list of locally discretized point and a list containing the abscissas' values.
        """
        return vm_common_operations.get_abscissa_discretization(self, abscissa1, abscissa2,
                                                                max_number_points, return_abscissas)

    def sort_points_along_curve(self, points: List[Union[volmdlr.Point2D, volmdlr.Point3D]]):
        """
        Sort point along a curve.

        :param points: list of points to be sorted.
        :return: sorted points.
        """
        return sorted(points, key=self.abscissa)


class LineSegment(Edge):
    """
    Abstract class.

    """

    def __init__(self, start: Union[volmdlr.Point2D, volmdlr.Point3D], end: Union[volmdlr.Point2D, volmdlr.Point3D],
                 line: [volmdlr_curves.Line2D, volmdlr_curves.Line3D] = None, name: str = ''):
        self.line = line
        Edge.__init__(self, start, end, name)

    def length(self):
        """Gets the length of a Line Segment."""
        if not self._length:
            self._length = self.end.point_distance(self.start)
        return self._length

    def abscissa(self, point, tol=1e-6):
        """
        Calculates the abscissa parameter of a Line Segment, at a point.

        :param point: point to verify abscissa.
        :param tol: tolerance.
        :return: abscissa parameter.
        """
        if point.point_distance(self.start) < tol:
            return 0
        if point.point_distance(self.end) < tol:
            return self.length()

        vector = self.end - self.start
        length = vector.norm()
        t_param = (point - self.start).dot(vector) / length
        if t_param < -1e-9 or t_param > length + 1e-9:
            raise ValueError(f'Point is not on linesegment: abscissa={t_param}')
        return t_param

    def direction_vector(self, abscissa=0.):
        """
        Returns a direction vector at a given abscissa, it is not normalized.

        :param abscissa: defines where in the line segment
            direction vector is to be calculated.
        :return: The direction vector of the LineSegment.
        """
        if not self._direction_vector_memo:
            self._direction_vector_memo = {}
        if abscissa not in self._direction_vector_memo:
            self._direction_vector_memo[abscissa] = self.end - self.start
        return self._direction_vector_memo[abscissa]

    def normal_vector(self, abscissa=0.):
        """
        Returns a normal vector at a given abscissa, it is not normalized.

        :param abscissa: defines where in the line_segment
        normal vector is to be calculated.
        :return: The normal vector of the LineSegment.
        """
        return self.direction_vector(abscissa).normal_vector()

    def point_projection(self, point):
        """
        Calculates the projection of a point on a Line Segment.

        :param point: point to be verified.
        :return: point projection.
        """
        point1, point2 = self.start, self.end
        vector = point2 - point1
        norm_u = vector.norm()
        t_param = (point - point1).dot(vector) / norm_u ** 2
        projection = point1 + t_param * vector

        return projection, t_param * norm_u

    def split(self, split_point, tol: float = 1e-6):
        """
        Split a Line Segment at a given point into two Line Segments.

        :param split_point: splitting point.
        :param tol: tolerance.
        :return: list with the two split line segments.
        """
        if split_point.is_close(self.start, tol):
            return [None, self.copy()]
        if split_point.is_close(self.end, tol):
            return [self.copy(), None]
        return [self.__class__(self.start, split_point),
                self.__class__(split_point, self.end)]

    def middle_point(self):
        """
        Calculates the middle point of a Line Segment.

        :return:
        """
        if not self._middle_point:
            self._middle_point = 0.5 * (self.start + self.end)
        return self._middle_point

    def point_at_abscissa(self, abscissa):
        """
        Calculates a point in the LineSegment at a given abscissa.

        :param abscissa: abscissa where in the curve the point should be calculated.
        :return: Corresponding point.
        """
        return self.start + self.unit_direction_vector() * abscissa

    def get_geo_lines(self, tag: int, start_point_tag: int, end_point_tag: int):
        """
        Gets the lines that define a LineSegment in a .geo file.

        :param tag: The linesegment index
        :type tag: int
        :param start_point_tag: The linesegment' start point index
        :type start_point_tag: int
        :param end_point_tag: The linesegment' end point index
        :type end_point_tag: int

        :return: A line
        :rtype: str
        """

        return 'Line(' + str(tag) + ') = {' + str(start_point_tag) + ', ' + str(end_point_tag) + '};'

    def get_geo_points(self):
        """Returns geo points."""
        return [self.start, self.end]

    def get_shared_section(self, other_linesegment, abs_tol: float = 1e-6):
        """
        Gets the shared section between two line segments.

        :param other_linesegment: other line segment to verify for shared section.
        :param abs_tol: tolerance.
        :return: shared line segment section.
        """
        if self.__class__ != other_linesegment.__class__:
            if self.__class__ == other_linesegment.simplify.__class__:
                return self.get_shared_section(other_linesegment.simplify)
            return []
        if not self.direction_vector().is_colinear_to(other_linesegment.direction_vector()) or \
                (not any(self.point_belongs(point, abs_tol)
                         for point in [other_linesegment.start, other_linesegment.end]) and
                 not any(other_linesegment.point_belongs(point, abs_tol) for point in [self.start, self.end])):
            return []
        if all(self.point_belongs(point) for point in other_linesegment.discretization_points(number_points=5)):
            return [other_linesegment]
        if all(other_linesegment.point_belongs(point) for point in self.discretization_points(number_points=5)):
            return [self]
        new_linesegment_points = []
        for point in [self.start, self.end]:
            if other_linesegment.point_belongs(point, abs_tol=abs_tol) and\
                    not volmdlr.core.point_in_list(point, new_linesegment_points):
                new_linesegment_points.append(point)
        for point in [other_linesegment.start, other_linesegment.end]:
            if self.point_belongs(point, abs_tol=abs_tol) and\
                    not volmdlr.core.point_in_list(point, new_linesegment_points):
                new_linesegment_points.append(point)
        if len(new_linesegment_points) == 1:
            return []
        if len(new_linesegment_points) != 2:
            raise ValueError
        class_ = self.__class__
        return [class_(new_linesegment_points[0], new_linesegment_points[1])]

    def delete_shared_section(self, other_linesegment, abs_tol: float = 1e-6):
        """
        Deletes from self, the section shared with the other line segment.

        :param other_linesegment:
        :param abs_tol: tolerance.
        :return:
        """
        shared_section = self.get_shared_section(other_linesegment, abs_tol)
        if not shared_section:
            return [self]
        points = []
        for point in [self.start, self.end, shared_section[0].start, shared_section[0].end]:
            if not volmdlr.core.point_in_list(point, points):
                points.append(point)
        points = sorted(points, key=self.start.point_distance)
        new_line_segments = []
        class_ = self.__class__
        for point1, point2 in zip(points[:-1], points[1:]):
            lineseg = class_(point1, point2)
            if not lineseg.direction_independent_is_close(shared_section[0]):
                new_line_segments.append(lineseg)
        return new_line_segments

    def straight_line_point_belongs(self, point):
        """
        Closing straight line point belongs verification.

        Verifies if a point belongs to the surface created by closing the edge with a
        line between its start and end points.

        :param point: Point to be verified.
        :return: Return True if the point belongs to this surface, or False otherwise.
        """
        return self.point_belongs(point)

    def point_belongs(self, point: Union[volmdlr.Point2D, volmdlr.Point3D], abs_tol: float = 1e-6):
        """
        Checks if a point belongs to the line segment. It uses the point_distance.

        :param point: The point to be checked
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :param abs_tol: The precision in terms of distance.
            Default value is 1e-6
        :type abs_tol: float, optional
        :return: `True` if the point belongs to the B-spline curve, `False`
            otherwise
        :rtype: bool
        """
        point_distance = self.point_distance(point)
        if math.isclose(point_distance, 0, abs_tol=abs_tol):
            return True
        return False

    def point_distance(self, point):
        """
        Abstract method.
        """
        raise NotImplementedError('the point_distance method must be'
                                  'overloaded by subclassing class')

    def to_step(self, current_id, *args, **kwargs):
        """Exports to STEP format."""
        line = self.line
        content, line_id = line.to_step(current_id)

        current_id = line_id + 1
        start_content, start_id = self.start.to_step(current_id, vertex=True)
        current_id = start_id + 1
        end_content, end_id = self.end.to_step(current_id + 1, vertex=True)
        content += start_content + end_content
        current_id = end_id + 1
        content += f"#{current_id} = EDGE_CURVE('{self.name}',#{start_id},#{end_id},#{line_id},.T.);\n"
        return content, current_id

    def is_close(self, other_edge, tol: float = 1e-6):
        """
        Checks if two line segments are the same considering the Euclidean distance.

        :param other_edge: other line segment.
        :param tol: The tolerance under which the Euclidean distance is considered equal to 0, defaults to 1e-6.
        :type tol: float, optional.
        """

        if isinstance(other_edge, self.__class__):
            if (self.start.is_close(other_edge.start, tol)
                    and self.end.is_close(other_edge.end, tol)):
                return True
        return False


class BSplineCurve(Edge):
    """
    An abstract class for B-spline curves.

    The following rule must be
    respected : `number of knots = number of control points + degree + 1`.

    :param degree: The degree of the B-spline curve.
    :type degree: int
    :param control_points: A list of 2 or 3-dimensional points
    :type control_points: Union[List[:class:`volmdlr.Point2D`],
        List[:class:`volmdlr.Point3D`]]
    :param knot_multiplicities: The vector of multiplicities for each knot
    :type knot_multiplicities: List[int]
    :param knots: The knot vector composed of values between 0 and 1
    :type knots: List[float]
    :param weights: The weight vector applied to the knot vector. Default
        value is None
    :type weights: List[float], optional
    :param periodic: If `True` the B-spline curve is periodic. Default value
        is False
    :type periodic: bool, optional
    :param name: The name of the B-spline curve. Default value is ''
    :type name: str, optional
    """

    def __init__(self,
                 degree: int,
                 control_points: Union[List[volmdlr.Point2D], List[volmdlr.Point3D]],
                 knot_multiplicities: List[int],
                 knots: List[float],
                 weights: List[float] = None,
                 name: str = ''):
        self.ctrlpts = [[*point] for point in control_points]
        self.degree = degree
        self.knots = nurbs_helpers.standardize_knot_vector(knots)
        self.knot_multiplicities = knot_multiplicities
        self.weights = weights

        Edge.__init__(self, control_points[0], control_points[-1], name=name)
        self._simplified = None
        self._delta = 0.01
        self._length = None
        self._eval_points = None
        self.ctrlptsw = None
        self.rational = False
        if self.weights:
            self.rational = True
            ctrlptsw = []
            for point, w in zip(self.control_points, weights):
                temp = [float(c * w) for c in point]
                temp.append(float(w))
                ctrlptsw.append(temp)
            self.ctrlptsw = ctrlptsw

    def __hash__(self):
        """
        Return a hash value for the B-spline curve.
        """
        return hash((tuple(self.control_points), self.degree, tuple(self.knots)))

    def __eq__(self, other):
        """
        Return True if the other B-spline curve has the same control points, degree, and knot vector, False otherwise.
        """
        if isinstance(other, self.__class__) and self.rational == other.rational:
            common_check = (self.control_points == other.control_points
                            and self.degree == other.degree
                            and self.knots == other.knots)
            if self.weights is None:
                return common_check
            return common_check and self.weights == other.weights
        return False

    def _data_eq(self, other):
        """
        Defines dessia common object equality.
        """
        return self == other

    @property
    def control_points(self):
        """Return the control points of the bspline curve."""
        point_name = "Point" + self.__class__.__name__[-2:]
        return [getattr(volmdlr, point_name)(*point) for point in self.ctrlpts]

    @property
    def knotvector(self):
        """Return the knot vector."""
        knot_vector = []
        for knot, knot_mut in zip(self.knots, self.knot_multiplicities):
            knot_vector.extend([knot] * knot_mut)
        return knot_vector

    @property
    def periodic(self):
        """Return True if the BSpline is periodic."""
        control_points = self.control_points
        return control_points[0].is_close(control_points[-1])

    @property
    def points(self):
        """
        Evaluate the BSpline points based on the set delta value of the curve.
        """
        if self._eval_points is None:
            self.evaluate()
        point_name = f'Point{self.__class__.__name__[-2:]}'
        return [getattr(volmdlr, point_name)(*point) for point in self._eval_points]

    @property
    def data(self):
        """
        Returns a dictionary of the BSpline data.
        """
        datadict = {
            "degree": self.degree,
            "knotvector": self.knotvector,
            "size": len(self.ctrlpts),
            "sample_size": self.sample_size,
            "rational": bool(self.weights),
            "dimension": 3 if self.__class__.__name__[-2:] == "3D" else 2,
            "precision": 18
        }
        if self.weights:
            datadict["control_points"] = tuple(self.ctrlptsw)
        else:
            datadict["control_points"] = tuple(self.ctrlpts)
        return datadict

    @property
    def sample_size(self):
        """
        Sample size.

        Sample size defines the number of evaluated points to generate. It also sets the ``delta`` property.

        :getter: Gets sample size
        :setter: Sets sample size
        :type: int
        """
        s_size = math.floor((1.0 / self.delta) + 0.5)
        return int(s_size)

    @sample_size.setter
    def sample_size(self, value):
        if not isinstance(value, int):
            raise ValueError("Sample size must be an integer value")

        # To make it operate like linspace, we have to know the starting and ending points.
        start = self.knotvector[self.degree]
        stop = self.knotvector[-(self.degree + 1)]

        # Set delta value
        self.delta = (stop - start) / float(value)

    @property
    def delta(self):
        """
        Evaluation delta.

        Evaluation delta corresponds to the *step size* while ``evaluate`` function iterates on the knot vector to
        generate curve points. Decreasing step size results in generation of more curve points.
        Therefore, smaller the delta value, smoother the curve.

        :getter: Gets the delta value
        :setter: Sets the delta value
        :type: float
        """
        return self._delta

    @delta.setter
    def delta(self, value):
        # Delta value for surface evaluation should be between 0 and 1
        if float(value) <= 0 or float(value) >= 1:
            raise ValueError("Curve evaluation delta should be between 0.0 and 1.0")

        # Clean up the curve points list
        self._points = None
        self._eval_points = None

        # Set new delta value
        self._delta = float(value)

    @property
    def domain(self):
        """
        Domain.

        Domain is determined using the knot vector(s).

        :getter: Gets the domain
        """
        return self.knotvector[self.degree], self.knotvector[-(self.degree + 1)]

    def to_geomdl(self):
        """Converts the BSpline curve into a geomdl curve."""
        if self.weights is None:
            curve = BSpline.Curve()
            curve.degree = self.degree
            curve.ctrlpts = self.ctrlpts
        else:
            curve = NURBS.Curve()
            curve.degree = self.degree
            curve.ctrlpts = self.ctrlpts
            curve.weights = self.weights
        curve.knotvector = self.knotvector
        curve.delta = self.delta
        return curve

    def to_dict(self, *args, **kwargs):
        """Avoids storing points in memo that makes serialization slow."""
        dict_ = self.base_dict()
        dict_['degree'] = self.degree
        dict_['control_points'] = [point.to_dict() for point in self.control_points]
        dict_['knot_multiplicities'] = self.knot_multiplicities
        dict_['knots'] = self.knots
        dict_['weights'] = self.weights
        return dict_

    def evaluate(self, **kwargs):
        """
        Evaluates the curve.

        The evaluated points are stored in :py:attr:`evalpts` property.

        Keyword Arguments:

            * ``start``: start parameter
            * ``stop``: stop parameter

        The ``start`` and ``stop`` parameters allow evaluation of a curve segment in the range *[start, stop]*, i.e.
        the curve will also be evaluated at the ``stop`` parameter value.

        The following examples illustrate the usage of the keyword arguments.

        """

        # Find evaluation start and stop parameter values
        start = kwargs.get('start', self.knotvector[self.degree])
        stop = kwargs.get('stop', self.knotvector[-(self.degree + 1)])

        # # Check parameters
        # if self._kv_normalize:
        #     if not utilities.check_params([start, stop]):
        #         raise GeomdlException("Parameters should be between 0 and 1")

        # Clean up the curve points
        self._points = None

        # Evaluate and cache
        self._eval_points = npy.asarray(evaluate_curve(self.data, start=start, stop=stop), dtype=npy.float64)

    def evaluate_single(self, u):
        """
        Calculates a point in the BSplineCurve at a given parameter u.

        :param u: Curve parameter. Must be a value between 0 and 1.
        :type u: float
        :return: Corresponding point.
        :rtype: Union[volmdlr.Point2D, Union[volmdlr.Point3D]
        """
        point_name = 'Point' + self.__class__.__name__[-2:]
        return getattr(volmdlr, point_name)(*evaluate_curve(self.data, u, u)[0])

    def derivatives(self, u, order):
        """
        Evaluates n-th order curve derivatives at the given parameter value.

        The output of this method is list of n-th order derivatives. If ``order`` is ``0``, then it will only output
        the evaluated point. Similarly, if ``order`` is ``2``, then it will output the evaluated point, 1st derivative
        and the 2nd derivative.

        :Example:

        Assuming a curve self is defined on a parametric domain [0.0, 1.0].
        Let's take the curve derivative at the parametric position u = 0.35.

        >>> derivatives = self.derivatives(u=0.35, order=2)
        >>> derivatives[0]  # evaluated point, equal to crv.evaluate_single(0.35)
        >>> derivatives[1]  # 1st derivative at u = 0.35
        >>> derivatives[2]  # 2nd derivative at u = 0.35

        :param u: parameter value
        :type u: float
        :param order: derivative order
        :type order: int
        :return: a list containing up to {order}-th derivative of the curve
        :rtype: Union[List[`volmdlr.Vector2D`], List[`volmdlr.Vector3D`]]
        """
        vector_name = 'Vector' + self.__class__.__name__[-2:]
        datadict = {
            "degree": self.degree,
            "knotvector": self.knotvector,
            "size": len(self.ctrlpts),
            "sample_size": self.sample_size,
            "rational": bool(self.weights),
            "dimension": 3 if vector_name == "Vector3D" else 2,
        }
        if self.weights:
            datadict["control_points"] = tuple(self.ctrlptsw)
        else:
            datadict["control_points"] = tuple(self.ctrlpts)
        return [getattr(volmdlr, vector_name)(*point)
                for point in derivatives_curve(datadict, u, order)]

    def split(self, point: Union[volmdlr.Point2D, volmdlr.Point3D],
              tol: float = 1e-6):
        """
        Splits of B-spline curve in two pieces using a 2D or 3D point.

        :param point: The point where the B-spline curve is split
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :param tol: The precision in terms of distance. Default value is 1e-6
        :type tol: float, optional
        :return: A list containing the first and second split of the B-spline
            curve
        :rtype: List[:class:`volmdlr.edges.BSplineCurve`]
        """
        if point.is_close(self.start, tol):
            return [None, self.copy()]
        if point.is_close(self.end, tol):
            return [self.copy(), None]
        parameter = round(self.point_to_parameter(point), 12)
        return split_curve(self, parameter)

    def get_reverse(self):
        """
        Reverses the BSpline's direction by reversing its control points.

        :return: A reversed B-Spline curve.
        :rtype: :class:`volmdlr.edges.BSplineCurve`.
        """
        return self.__class__(
            degree=self.degree,
            control_points=self.control_points[::-1],
            knot_multiplicities=self.knot_multiplicities[::-1],
            knots=self.knots[::-1],
            weights=self.weights,
        )

    @property
    def simplify(self):
        """Search another simplified edge that can represent the bspline."""
        if self.length() < 1e-6:
            return self
        class_sufix = self.__class__.__name__[-2:]
        if self._simplified is None:
            points = self.points
            if self.periodic:
                fullarc_class_ = getattr(sys.modules[__name__], 'FullArc' + class_sufix)
                n = len(points)
                try_fullarc = fullarc_class_.from_3_points(points[0], points[int(0.5 * n)],
                                                           points[int(0.75 * n)])

                if all(try_fullarc.point_belongs(point, 1e-6) for point in points):
                    self._simplified = try_fullarc
                    return try_fullarc
            else:
                lineseg_class = getattr(sys.modules[__name__], 'LineSegment' + class_sufix)
                lineseg = lineseg_class(points[0], points[-1])
                if all(lineseg.point_belongs(pt) for pt in points):
                    self._simplified = lineseg
                    return lineseg
                interior = self.point_at_abscissa(0.5 * self.length())
                vector1 = interior - self.start
                vector2 = interior - self.end
                if vector1.is_colinear_to(vector2) or vector1.norm() == 0 or vector2.norm() == 0:
                    return self
                arc_class_ = getattr(sys.modules[__name__], 'Arc' + class_sufix)
                try_arc = arc_class_.from_3_points(self.start, interior, self.end)
                if all(try_arc.point_belongs(point, 1e-6) for point in points):
                    self._simplified = try_arc
                    return try_arc
            self._simplified = self
        return self._simplified

    @classmethod
    def from_geomdl_curve(cls, curve, name: str = ""):
        """
        # TODO: to be completed.

        :param curve:
        :type curve:
        :param name: curve name.
        :return: A reversed B-spline curve
        :rtype: :class:`volmdlr.edges.BSplineCurve`
        """
        point_dimension = f'Point{cls.__name__[-2::]}'

        knots = list(sorted(set(curve.knotvector)))
        knot_multiplicities = [curve.knotvector.count(k) for k in knots]
        return cls(degree=curve.degree,
                   control_points=[getattr(volmdlr, point_dimension)(*point)
                                   for point in curve.ctrlpts],
                   knots=knots,
                   knot_multiplicities=knot_multiplicities,
                   weights=curve.weights, name=name)

    def length(self):
        """
        Returns the length of the B-spline curve.

        :return: The length of the B-spline curve.
        :rtype: float
        """
        if not self._length:
            if self._eval_points is None and self.delta == 0.01:
                self.evaluate()
                points = self._eval_points
            elif self.delta == 0.01:
                points = self._eval_points
            else:
                datadict = self.data
                datadict["sample_size"] = 100
                start, stop = self.domain
                points = npy.asarray(evaluate_curve(datadict, start=start, stop=stop), dtype=npy.float64)

            differences = npy.diff(points, axis=0)

            squared_distances = npy.sum(differences ** 2, axis=1)

            self._length = float(npy.sum(npy.sqrt(squared_distances)))
            # self._length = length_curve(self.curve)
        return self._length

    def normal_vector(self, abscissa):
        """
        Calculates the normal vector to the BSpline curve at given abscissa.

        :return: the normal vector
        """
        return self.direction_vector(abscissa).deterministic_unit_normal_vector()

    def direction_vector(self, abscissa):
        """
        Calculates the direction vector on the BSpline curve at given abscissa.

        :param abscissa: edge abscissa
        :return: direction vector
        """
        u = self.abscissa_to_parameter(abscissa)
        derivatives = self.derivatives(u, 1)
        return derivatives[1]

    def point_to_parameter(self, point: Union[volmdlr.Point2D, volmdlr.Point3D]):
        """
        Search for the value of the normalized evaluation parameter u.

        :return: the given point when the BSplineCurve3D is evaluated at the u value.
        """
        abscissa = self.abscissa(point)
        u = max(min(abscissa / self.length(), 1.), 0.0)
        u_min, u_max = self.domain
        if u_min != 0 or u_max != 1.0:
            u = u * (u_max - u_min) + u_min
        return u

    def abscissa_to_parameter(self, abscissa: float):
        """
        Search for the value of the normalized evaluation parameter u.

        :return: the given point when the BSplineCurve3D is evaluated at the u value.
        """
        u = max(min(abscissa / self.length(), 1.), 0.0)
        u_min, u_max = self.domain
        if u_min != 0 or u_max != 1.0:
            u = u * (u_max - u_min) + u_min
        return u

    def abscissa(self, point: Union[volmdlr.Point2D, volmdlr.Point3D],
                 tol: float = 1e-7):
        """
        Computes the abscissa of a 2D or 3D point using the least square method.

        :param point: The point located on the B-spline curve.
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`].
        :param tol: The precision in terms of distance. Default value is 1e-6.
        :type tol: float, optional.
        :return: The abscissa of the point.
        :rtype: float
        """
        if point.is_close(self.start):
            return 0
        if point.is_close(self.end):
            return self.length()
        length = self.length()
        point_array = npy.array(list(point), dtype=npy.float64)
        distances = npy.linalg.norm(self._eval_points - point_array, axis=1)
        index = npy.argmin(distances)
        u_min, u_max = self.domain
        u0 = u_min + index * (u_max - u_min) / (self.sample_size - 1)
        u, convergence_sucess = self.point_invertion(u0, point)
        if u_min != 0 or u_max != 1.0:
            u = (u - u_min) / (u_max - u_min)
        abscissa = u * length
        if convergence_sucess:  # sometimes we don't achieve convergence with a given initial guess
            return float(abscissa)

        def evaluate_point_distance(u_param):
            return (point - self.evaluate_single(u_param)).norm()
        results = [(abscissa, evaluate_point_distance(u))]
        initial_condition_list = npy.linspace(u_min, u_max, 10).tolist()
        initial_condition_list.sort(key=evaluate_point_distance)
        for u0 in initial_condition_list[:2]:
            u, convergence_sucess = self.point_invertion(u0, point)
            if u_min != 0 or u_max != 1.0:
                u = (u - u_min) / (u_max - u_min)
            abscissa = u * length
            if convergence_sucess:  # sometimes we don't achieve convergence with a given initial guess
                return float(abscissa)
            dist = evaluate_point_distance(u)
            if dist < tol:
                return float(abscissa)
            results.append((abscissa, dist))
        result = min(results, key=lambda r: r[1])[0]
        return float(result)

    def _point_inversion_funcs(self, u, point):
        """
        Helper function to evaluate Newton-Rapshon terms.
        """
        curve_derivatives = self.derivatives(u, 2)
        distance_vector = curve_derivatives[0] - point
        func = curve_derivatives[1].dot(distance_vector)
        func_first_derivative = curve_derivatives[2].dot(distance_vector) + curve_derivatives[1].norm() ** 2
        return func, func_first_derivative, curve_derivatives, distance_vector

    def point_invertion(self, u0: float, point, maxiter: int = 50, tol1: float = 1e-7, tol2: float = 1e-8):
        """
        Finds the equivalent B-Spline curve parameter u to a given a point 3D or 2D using an initial guess u0.

        :param u0: An initial guess between 0 and 1.
        :type u0: float
        :param point: Point to evaluation.
        :type point: Union[volmdlr.Point2D, volmdlr.Point3D]
        :param maxiter: Maximum number of iterations.
        :type maxiter: int
        :param tol1: Distance tolerance to stop.
        :type tol1: float
        :param tol2: Zero cos tolerance to stop.
        :type tol2: float
        :return: u parameter and convergence check
        :rtype: int, bool
        """
        if maxiter == 0:
            return u0, False
        func, func_first_derivative, curve_derivatives, distance_vector = self._point_inversion_funcs(u0, point)
        if self._check_convergence(curve_derivatives, distance_vector, tol1=tol1, tol2=tol2):
            return u0, True
        new_u = u0 - func / (func_first_derivative + 1e-18)
        new_u = self._check_bounds(new_u)
        residual = (new_u - u0) * curve_derivatives[1]
        if residual.norm() <= tol1:
            return u0, False
        u0 = new_u
        return self.point_invertion(u0, point, maxiter=maxiter - 1)

    @staticmethod
    def _check_convergence(curve_derivatives, distance_vector, tol1: float = 1e-7, tol2: float = 1e-8):
        """
        Helper function to check convergence of point_invertion method.
        """
        distance = distance_vector.norm()
        if distance <= tol1:
            return True
        if curve_derivatives[1].norm() == 0.0:
            return False
        zero_cos = abs(curve_derivatives[1].dot(distance_vector)) / curve_derivatives[1].norm() * distance
        if distance <= tol1 and zero_cos <= tol2:
            return True
        return False

    def _check_bounds(self, u):
        """
        Helper function to check if evaluated parameters in point_invertion method are contained in the bspline domain.
        """
        a, b = self.domain
        if self.periodic:
            if u < a:
                u = b - (a - u)
            elif u > b:
                u = a + (u - b)
        if u < a:
            u = a

        elif u > b:
            u = b
        return u

    def translation(self, offset: Union[volmdlr.Vector2D, volmdlr.Vector3D]):
        """
        Translates the B-spline curve.

        :param offset: The translation vector
        :type offset: Union[:class:`volmdlr.Vector2D`,
            :class:`volmdlr.Vector3D`]
        :return: A new translated BSplineCurve
        :rtype: :class:`volmdlr.edges.BSplineCurve`
        """
        control_points = [point.translation(offset)
                          for point in self.control_points]
        return self.__class__(self.degree, control_points,
                              self.knot_multiplicities, self.knots,
                              self.weights)

    def point_belongs(self, point: Union[volmdlr.Point2D, volmdlr.Point3D], abs_tol: float = 1e-6):
        """
        Checks if a 2D or 3D point belongs to the B-spline curve or not. It uses the point_distance.

        :param point: The point to be checked.
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :param abs_tol: The precision in terms of distance.
            Default value is 1e-6
        :type abs_tol: float, optional.
        :return: `True` if the point belongs to the B-spline curve, `False`
            otherwise
        :rtype: bool
        """
        if self.point_distance(point) < abs_tol:
            return True
        return False

    def point_distance(self, point: Union[volmdlr.Point2D, volmdlr.Point3D]):
        """
        Calculates the distance from a given point to a BSplineCurve2D or 3D.

        :param point: The point to be checked.
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :return: distance.
        """

        return self.point_distance_to_edge(point)

    def merge_with(self, bspline_curve: 'BSplineCurve'):
        """
        Merges consecutive B-spline curves to define a new merged one.

        :param bspline_curve: Another B-spline curve
        :type bspline_curve: :class:`volmdlr.edges.BSplineCurve`
        :return: A merged B-spline curve
        :rtype: :class:`volmdlr.edges.BSplineCurve`
        """
        point_dimension = f'Wire{self.__class__.__name__[-2::]}'
        wire = getattr(volmdlr.wires, point_dimension)(bspline_curve)
        ordered_wire = wire.order_wire()

        points, n = [], 10
        for primitive in ordered_wire.primitives:
            points.extend(primitive.discretization_points(n))
        points.pop(n + 1)

        return self.__class__.from_points_interpolation(points, min(self.degree, bspline_curve.degree))

    @classmethod
    def from_bsplines(cls, bsplines: List['BSplineCurve'],
                      discretization_points: int = 10, name: str = ''):
        """
        Creates a B-spline curve from a list of B-spline curves.

        :param bsplines: A list of B-spline curve
        :type bsplines: List[:class:`volmdlr.edges.BSplineCurve`]
        :param discretization_points: The number of points for the
            discretization. Default value is 10
        :type discretization_points: int, optional.
        :param name: object's name.
        :return: A merged B-spline curve
        :rtype: :class:`volmdlr.edges.BSplineCurve`
        """
        point_dimension = f'Wire{cls.__name__[-2::]}'
        wire = getattr(volmdlr.wires, point_dimension)(bsplines)
        ordered_wire = wire.order_wire()

        points, degree = [], []
        for i, primitive in enumerate(ordered_wire.primitives):
            degree.append(primitive.degree)
            if i == 0:
                points.extend(primitive.discretization_points(number_points=discretization_points))
            else:
                points.extend(
                    primitive.discretization_points(number_points=discretization_points)[1::])

        return cls.from_points_interpolation(points, min(degree), name=name)

    @classmethod
    def from_points_approximation(cls, points: Union[List[volmdlr.Point2D], List[volmdlr.Point3D]],
                                  degree: int, name: str = "", **kwargs):
        """
        Creates a B-spline curve approximation using least squares method with fixed number of control points.

        It is recommended to specify the
        number of control points.
        Please refer to The NURBS Book (2nd Edition), pp.410-413 for details.

        :param points: The data points
        :type points: Union[List[:class:`volmdlr.Point2D`],
            List[:class:`volmdlr.Point3D`]]
        :param degree: The degree of the output parametric curve
        :type degree: int
        :param name: (optional) Curve name.
        :param kwargs: See below
        :return: A B-spline curve from points approximation
        :rtype: :class:`volmdlr.edges.BSplineCurve`
        :keyword centripetal: Activates centripetal parametrization method.
            Default value is False
        :keyword ctrlpts_size: Number of control points. Default value is
            len(points) - 1
        """
        point_name = 'Point' + points[0].__class__.__name__[-2:]
        control_points, knots, knot_multiplicities = fitting.approximate_curve(
            npy.asarray([npy.asarray([*point], dtype=npy.float64) for point in points], dtype=npy.float64),
            degree, **kwargs)
        control_points = [getattr(volmdlr, point_name)(*point) for point in control_points]
        return cls(degree, control_points, knot_multiplicities, knots, name=name)

    def tangent(self, position: float = 0.0, normalize: bool = True):
        """
        Evaluates the tangent vector of the B-spline curve at the input parameter value.

        :param position: Value of the parameter, between 0 and 1
        :type position: float
        :param normalize: By default return a normalized tangent vector.
        :return: The tangent vector
        :rtype: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        """
        # 1st derivative of the curve gives the tangent
        ders = self.derivatives(position, 1)
        tangent = ders[1].unit_vector() if normalize else ders[1]
        return tangent

    @classmethod
    def from_points_interpolation(cls, points: Union[List[volmdlr.Point2D], List[volmdlr.Point3D]],
                                  degree: int, centripetal: bool = True, name: str = " "):
        """
        Creates a B-spline curve interpolation through the data points.

        Please refer to Algorithm A9.1 on The NURBS Book (2nd Edition),
        pp.369-370 for details.

        :param points: The data points
        :type points: Union[List[:class:`volmdlr.Point2D`],
            List[:class:`volmdlr.Point3D`]]
        :param degree: The degree of the output parametric curve
        :type degree: int
        :param centripetal: Please refer to Algorithm A9.1 on The NURBS Book (2nd Edition),
        pp.369-370 for details.
        :type centripetal: bool
        :param name: curve name.
        :return: A B-spline curve from points interpolation
        :rtype: :class:`volmdlr.edges.BSplineCurve`
        """
        if degree >= len(points):
            raise ValueError("Number of points for interpolation must be at least degree + 1")
        set_points = set(points)
        if len(set_points) < len(points) - 1:
            warnings.warn("Not able to perform point interpolation."
                          "There are repeated points not in the edges of the point list.")
            return None
        point_name = 'Point' + points[0].__class__.__name__[-2:]
        ctrlpts, knots, knot_multiplicities = fitting.interpolate_curve(
            npy.asarray([npy.asarray([*point], dtype=npy.float64) for point in points], dtype=npy.float64),
            degree, centripetal=centripetal)
        ctrlpts = [getattr(volmdlr, point_name)(*point) for point in ctrlpts]
        return cls(degree, ctrlpts, knot_multiplicities, knots, name=name)

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = None):
        """
        Linear spaced discretization of the curve.

        :param number_points: The number of points to include in the discretization.
        :type number_points: int
        :param angle_resolution: The resolution of the angle to use when calculating the number of points.
        :type angle_resolution: int
        :return: A list of discretized points on the B-spline curve.
        :rtype: List[`volmdlr.Point2D] or List[`volmdlr.Point3D]
        """

        if angle_resolution:
            number_points = int(math.pi * angle_resolution)

        if self.sample_size == number_points or (not number_points and not angle_resolution):
            return self.points

        datadict = self.data
        datadict["sample_size"] = number_points
        start, stop = self.domain
        points_list = evaluate_curve(datadict, start, stop)
        point_name = 'Point' + self.__class__.__name__[-2:]
        return [getattr(volmdlr, point_name)(*point) for point in points_list]

    def get_geo_lines(self, tag: int, control_points_tags: List[int]):
        """
        Gets the lines that define a BsplineCurve in a .geo file.

        :param tag: The BsplineCurve index
        :type tag: int
        :param start_point_tag: The linesegment' start point index
        :type start_point_tag: int
        :param end_point_tag: The linesegment' end point index
        :type end_point_tag: int

        :return: A line
        :rtype: str
        """

        return 'BSpline(' + str(tag) + ') = {' + str(control_points_tags)[1:-1] + '};'

    def get_geo_points(self):
        """Gets the points that define a BsplineCurve in a .geo file."""
        return list(self.discretization_points())

    def line_intersections(self, line, tol: float = 1e-6):
        """
        Calculates the intersections of a BSplineCurve (2D or 3D) with a Line (2D or 3D).

        :param line: line to verify intersections
        :return: list of intersections
        """
        # if self.
        polygon_points = []
        for point in self.points:
            if not volmdlr.core.point_in_list(point, polygon_points):
                polygon_points.append(point)
        list_intersections = []
        initial_abscissa = 0
        linesegment_name = 'LineSegment' + self.__class__.__name__[-2:]
        for points in zip(polygon_points[:-1], polygon_points[1:]):
            linesegment = getattr(sys.modules[__name__], linesegment_name)(points[0], points[1])
            intersections = linesegment.line_intersections(line)

            if not intersections and linesegment.direction_vector().is_colinear_to(line.direction_vector()):
                if line.point_distance(linesegment.middle_point()) < (tol * 0.01):
                    list_intersections.append(linesegment.middle_point())
            if intersections and intersections[0] not in list_intersections:
                if self.point_belongs(intersections[0], tol):
                    list_intersections.append(intersections[0])
                    continue
                abs1 = self.abscissa(linesegment.start)
                abs2 = self.abscissa(linesegment.end)
                list_abscissas = list(new_abscissa for new_abscissa in npy.linspace(abs1, abs2, 1000))
                intersection = self.select_intersection_point(list_abscissas, intersections, line)
                list_intersections.append(intersection)
            initial_abscissa += linesegment.length()
        return list_intersections

    def select_intersection_point(self, list_abscissas, intersections, line, abs_tol: float = 1e-7):
        """
        Select closest point in curve to intersection point obtained with discretized linesegment.

        :param list_abscissas: list of abscissas to verify the closest point.
        :param intersections: intersection with discretized line.
        :param line: other line.
        :param abs_tol: tolerance allowed.
        :return:
        """
        distance = npy.inf
        intersection = None
        for i_abscissa in list_abscissas:
            point_in_curve = BSplineCurve.point_at_abscissa(self, i_abscissa)
            if line.point_distance(point_in_curve) <= abs_tol:
                return point_in_curve
            dist = point_in_curve.point_distance(intersections[0])
            if dist < distance:
                distance = dist
                intersection = point_in_curve
            else:
                break
        return intersection

    def get_linesegment_intersections(self, linesegment):
        """
        Calculates intersections between a BSplineCurve and a LineSegment.

        :param linesegment: linesegment to verify intersections.
        :return: list with the intersections points.
        """
        results = self.line_intersections(linesegment.line)
        intersections_points = []
        for result in results:
            if linesegment.point_belongs(result, 1e-5):
                intersections_points.append(result)
        return intersections_points

    def point_at_abscissa(self, abscissa):
        """
        Calculates a point in the BSplineCurve at a given abscissa.

        :param abscissa: abscissa where in the curve the point should be calculated.
        :return: Corresponding point.
        """
        u = self.abscissa_to_parameter(abscissa)
        point_name = 'Point' + self.__class__.__name__[-2:]
        return getattr(volmdlr, point_name)(*self.evaluate_single(u))

    def get_shared_section(self, other_bspline2, abs_tol: float = 1e-6):
        """
        Gets the shared section between two BSpline curves.

        :param other_bspline2: other arc to verify for shared section.
        :param abs_tol: tolerance.
        :return: shared arc section.
        """
        if self.__class__ != other_bspline2.__class__:
            if self.simplify.__class__ == other_bspline2.__class__:
                return self.simplify.get_shared_section(other_bspline2, abs_tol)
            return []
        if not self.is_shared_section_possible(other_bspline2, 1e-7):
            return []
        if not any(self.point_belongs(point, abs_tol=abs_tol)
                   for point in other_bspline2.discretization_points(number_points=10)):
            return []
        if all(self.point_belongs(point, abs_tol=abs_tol) for point in other_bspline2.points):
            return [other_bspline2]
        if all(other_bspline2.point_belongs(point, abs_tol=abs_tol) for point in self.points):
            return [self]
        if self.point_belongs(other_bspline2.start, abs_tol=abs_tol):
            bspline1_, bspline2_ = self.split(other_bspline2.start, tol=abs_tol)
        elif self.point_belongs(other_bspline2.end, abs_tol=abs_tol):
            bspline1_, bspline2_ = self.split(other_bspline2.end, tol=abs_tol)
        else:
            return []
            # raise NotImplementedError
        return self._get_shared_section_from_split(bspline1_, bspline2_, other_bspline2, abs_tol)

    def is_shared_section_possible(self, other_bspline2, tol):
        """
        Verifies if it there is any possibility of the two bsplines share a section.

        :param other_bspline2: other bspline.
        :param tol: tolerance used.
        :return: True or False.
        """
        raise NotImplementedError(f"is_shared_section_possible is not yet implemented by {self.__class__.__name__}")

    @staticmethod
    def _get_shared_section_from_split(bspline1_, bspline2_, other_bspline2, abs_tol):
        """
        Helper function to get_shared_section.
        """
        shared_bspline_section = []
        for bspline in [bspline1_, bspline2_]:
            if bspline and all(other_bspline2.point_belongs(point, abs_tol=abs_tol)
                               for point in bspline.discretization_points(number_points=10)):
                shared_bspline_section.append(bspline)
                break
        return shared_bspline_section

    def delete_shared_section(self, other_bspline2, abs_tol: float = 1e-6):
        """
        Deletes from self, the section shared with the other arc.

        :param other_bspline2:
        :param abs_tol: tolerance.
        :return:
        """
        shared_section = self.get_shared_section(other_bspline2, abs_tol)
        if not shared_section:
            return [self]
        if shared_section == self:
            return []
        split_bspline1 = self.split(shared_section[0].start)
        split_bspline2 = self.split(shared_section[0].end)
        new_arcs = []
        shared_section_middle_point = shared_section[0].point_at_abscissa(0.5 * shared_section[0].length())
        for arc in split_bspline1 + split_bspline2:
            if arc and not arc.point_belongs(shared_section_middle_point, abs_tol=abs_tol):
                new_arcs.append(arc)
        return new_arcs

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point: Point to be verified
        :return: Return True if the point belongs to this surface,
            or False otherwise
        """
        raise NotImplementedError(f'the straight_line_point_belongs method must be'
                                  f' overloaded by {self.__class__.__name__}')

    def point_projection(self, point):
        """
        Calculates the projection of a point on the B-Spline.

        :param point: point to be verified.
        :return: point projection.
        """
        return [self.point_at_abscissa(self.abscissa(point))]

    def local_discretization(self, point1, point2, number_points: int = 10, tol: float = 1e-6):
        """
        Gets n discretization points between two given points of the edge.

        :param point1: point 1 on edge.
        :param point2: point 2 on edge.
        :param number_points: number of points to discretize locally.
        :return: list of locally discretized points.
        """
        abscissa1 = self.abscissa(point1)
        abscissa2 = self.abscissa(point2)
        # special case periodical bsplinecurve
        add_point_at_end = False
        if self.periodic:
            if math.isclose(abscissa2, 0.0, abs_tol=tol) or abscissa1 >= abscissa2:
                abscissa2 += self.length()
            if point1.is_close(point2):
                add_point_at_end = True

        discretized_points_between_1_2 = []
        length = self.length()
        for abscissa in npy.linspace(abscissa1, abscissa2, num=number_points):
            if self.periodic and abscissa > length:
                abscissa -= length
            abscissa_point = self.point_at_abscissa(abscissa)
            if not volmdlr.core.point_in_list(abscissa_point, discretized_points_between_1_2, tol=tol):
                discretized_points_between_1_2.append(abscissa_point)
        if add_point_at_end:
            discretized_points_between_1_2 += [discretized_points_between_1_2[0]]
        return discretized_points_between_1_2

    def is_close(self, other_edge, tol: float = 1e-6):
        """
        Checks if two bsplines are the same considering the Euclidean distance.

        :param other_edge: other bspline.
        :param tol: The tolerance under which the Euclidean distance is considered equal to 0, defaults to 1e-6.
        :type tol: float, optional
        """
        if isinstance(other_edge, self.__class__):
            if self.start.is_close(other_edge.start) and self.end.is_close(other_edge.end):
                is_true = True
                for point in other_edge.discretization_points(number_points=10):
                    if not self.point_belongs(point):
                        is_true = False
                        break
                if is_true:
                    return True
        return False

    def frame_mapping(self, frame: Union[volmdlr.Frame3D, volmdlr.Frame2D], side: str):
        """
        Returns a new Revolution Surface positioned in the specified frame.

        :param frame: Frame of reference
        :type frame: `volmdlr.Frame3D`
        :param side: 'old' or 'new'
        """
        new_control_points = [control_point.frame_mapping(frame, side) for control_point in self.control_points]
        return self.__class__(self.degree, new_control_points, self.knot_multiplicities, self.knots, self.weights,
                              self.name)


class BSplineCurve2D(BSplineCurve):
    """
    A class for 2-dimensional B-spline curves.

    The following rule must be
    respected : `number of knots = number of control points + degree + 1`.

    :param degree: The degree of the 2-dimensional B-spline curve
    :type degree: int
    :param control_points: A list of 2-dimensional points
    :type control_points: List[:class:`volmdlr.Point2D`]
    :param knot_multiplicities: The vector of multiplicities for each knot
    :type knot_multiplicities: List[int]
    :param knots: The knot vector composed of values between 0 and 1
    :type knots: List[float]
    :param weights: The weight vector applied to the knot vector. Default
        value is None
    :type weights: List[float], optional
    :param periodic: If `True` the B-spline curve is periodic. Default value
        is False
    :type periodic: bool, optional
    :param name: The name of the B-spline curve. Default value is ''
    :type name: str, optional
    """

    def __init__(self,
                 degree: int,
                 control_points: List[volmdlr.Point2D],
                 knot_multiplicities: List[int],
                 knots: List[float],
                 weights: List[float] = None,
                 name: str = ''):
        self._bounding_rectangle = None

        BSplineCurve.__init__(self, degree,
                              control_points,
                              knot_multiplicities,
                              knots,
                              weights,
                              name)
        self._bounding_rectangle = None
        self._length = None

    @property
    def bounding_rectangle(self):
        """
        Computes the bounding rectangle of the 2-dimensional B-spline curve.

        :return: The bounding rectangle.
        :rtype: :class:`volmdlr.core.BoundingRectangle`
        """
        if not self._bounding_rectangle:
            self._bounding_rectangle = volmdlr.core.BoundingRectangle.from_points(self.points)
        return self._bounding_rectangle

    def straight_line_area(self):
        """
        Uses shoelace algorithm for evaluating the area.
        """
        points = self.discretization_points(number_points=100)
        x = [point.x for point in points]
        y = [point.y for point in points]
        x1 = [x[-1]] + x[0:-1]
        y1 = [y[-1]] + y[0:-1]
        return 0.5 * abs(sum(i * j for i, j in zip(x, y1))
                         - sum(i * j for i, j in zip(y, x1)))

    def straight_line_center_of_mass(self):
        """Straight line center of mass."""
        polygon_points = self.discretization_points(number_points=100)
        cog = volmdlr.O2D
        for point in polygon_points:
            cog += point
        cog = cog / len(polygon_points)
        return cog

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plot a B-Spline curve 2D."""
        if ax is None:
            _, ax = plt.subplots()

        points = self.points

        x_points = [point.x for point in points]
        y_points = [point.y for point in points]
        ax.plot(x_points, y_points, color=edge_style.color, alpha=edge_style.alpha)
        if edge_style.plot_points:
            for point in points:
                point.plot(ax, color=edge_style.color)
        return ax

    def to_3d(self, plane_origin, x1, x2):
        """Transforms a B-Spline Curve 2D in 3D."""
        control_points3d = [point.to_3d(plane_origin, x1, x2) for point in
                            self.control_points]
        return BSplineCurve3D(self.degree, control_points3d,
                              self.knot_multiplicities, self.knots,
                              self.weights)

    def to_step(self, current_id, *args, **kwargs):
        """Exports to STEP format."""
        points_ids = []
        content = ''
        point_id = current_id
        for point in self.control_points:
            point_content, point_id = point.to_step(point_id,
                                                    vertex=False)
            content += point_content
            points_ids.append(point_id)
            point_id += 1

        content += f"#{point_id} = B_SPLINE_CURVE_WITH_KNOTS('{self.name}',{self.degree}," \
                   f"({volmdlr.core.step_ids_to_str(points_ids)})," \
                   f".UNSPECIFIED.,.F.,.F.,{tuple(self.knot_multiplicities)},{tuple(self.knots)},.UNSPECIFIED.);\n"
        return content, point_id + 1

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        BSplineCurve2D rotation.

        :param center: rotation center
        :param angle: angle rotation
        :return: a new rotated Line2D
        """
        control_points = [point.rotation(center, angle)
                          for point in self.control_points]
        return BSplineCurve2D(self.degree, control_points,
                              self.knot_multiplicities, self.knots,
                              self.weights)

    def line_crossings(self, line2d: volmdlr_curves.Line2D):
        """Bspline Curve crossings with a line 2d."""
        polygon_points = self.discretization_points(number_points=50)
        crossings = []
        for p1, p2 in zip(polygon_points[:-1], polygon_points[1:]):
            linesegment = LineSegment2D(p1, p2)
            crossings.extend(linesegment.line_crossings(line2d))
        return crossings

    def get_reverse(self):
        """
        Reverse the BSpline's direction by reversing its start and end points.

        """

        return self.__class__(degree=self.degree,
                              control_points=self.control_points[::-1],
                              knot_multiplicities=self.knot_multiplicities[::-1],
                              knots=self.knots[::-1],
                              weights=self.weights)

    def nearest_point_to(self, point):
        """
        Find out the nearest point on the linesegment to point.

        """

        points = self.discretization_points(number_points=500)
        return point.nearest_point(points)

    def linesegment_intersections(self, linesegment2d, abs_tol: float = 1e-6):
        """
        Calculates intersections between a BSpline Curve 2D and a Line Segment 2D.

        :param linesegment2d: line segment to verify intersections.
        :param abs_tol: tolerance.
        :return: list with the intersections points.
        """
        if self.bounding_rectangle.distance_to_b_rectangle(linesegment2d.bounding_rectangle) > abs_tol:
            return []
        intersections_points = vm_utils_intersections.get_bsplinecurve_intersections(
            linesegment2d, self, abs_tol=abs_tol)
        return intersections_points

    def arc_intersections(self, arc, abs_tol=1e-6):
        """
        Calculates intersections between a BSpline Curve 2D and an arc 2D.

        :param arc: arc to verify intersections.
        :param abs_tol: tolerance.
        :return: list with the intersections points.
        """
        if self.bounding_rectangle.distance_to_b_rectangle(arc.bounding_rectangle) > abs_tol:
            return []
        return self._generic_edge_intersections(arc, abs_tol)

    def bsplinecurve_intersections(self, bspline, abs_tol=1e-6):
        """
        Calculates intersections between a two BSpline Curve 2D.

        :param bspline: bspline to verify intersections.
        :param abs_tol: tolerance.
        :return: list with the intersections points.
        """
        if self.bounding_rectangle.distance_to_b_rectangle(bspline.bounding_rectangle) > abs_tol:
            return []
        return self._generic_edge_intersections(bspline, abs_tol)

    def axial_symmetry(self, line):
        """
        Finds out the symmetric bsplinecurve2d according to a line.

        """

        points_symmetry = [point.axial_symmetry(line) for point in self.control_points]

        return self.__class__(degree=self.degree,
                              control_points=points_symmetry,
                              knot_multiplicities=self.knot_multiplicities[::-1],
                              knots=self.knots[::-1],
                              weights=self.weights)

    def offset(self, offset_length: float):
        """
        Offsets a BSplineCurve2D in one of its normal direction.

        :param offset_length: the length taken to offset the BSpline. if positive, the offset is in the normal
            direction of the curve. if negative, in the opposite direction of the normal.
        :return: returns an offset bsplinecurve2D, created with from_points_interpolation.
        """
        points = self.points
        unit_normal_vectors = [self.unit_normal_vector(
            self.abscissa(point)) for point in points]
        offseted_points = [point.translation(normal_vector * offset_length) for point, normal_vector
                           in zip(points, unit_normal_vectors)]
        offseted_bspline = BSplineCurve2D.from_points_interpolation(offseted_points, self.degree)
        return offseted_bspline

    def is_shared_section_possible(self, other_bspline2, tol):
        """
        Verifies if it there is any possibility of the two bsplines share a section.

        :param other_bspline2: other bspline.
        :param tol: tolerance used.
        :return: True or False.
        """
        if self.bounding_rectangle.distance_to_b_rectangle(other_bspline2.bounding_rectangle) > tol:
            return False
        return True

    def normal(self, position: float = 0.0):
        """Normal vector to BPlineCurve2D."""
        der = self.derivatives(position, 1)
        tangent = der[1].unit_vector()
        return tangent.rotation(der[0], 0.5 * math.pi)


class BezierCurve2D(BSplineCurve2D):
    """
    A class for 2-dimensional Bézier curves.

    :param degree: The degree of the Bézier curve.
    :type degree: int
    :param control_points: A list of 2-dimensional points
    :type control_points: List[:class:`volmdlr.Point2D`]
    :param name: The name of the B-spline curve. Default value is ''
    :type name: str, optional
    """

    def __init__(self, degree: int, control_points: List[volmdlr.Point2D],
                 name: str = ''):
        knotvector = nurbs_helpers.generate_knot_vector(degree,
                                                        len(control_points))
        knot_multiplicity = [1] * len(knotvector)

        BSplineCurve2D.__init__(self, degree, control_points,
                                knot_multiplicity, knotvector,
                                None, name)


class LineSegment2D(LineSegment):
    """
    Define a line segment limited by two points.

    """

    def __init__(self, start: volmdlr.Point2D, end: volmdlr.Point2D, *,
                 line: volmdlr_curves.Line2D = None, name: str = ''):
        if start.is_close(end, 1e-6):
            raise NotImplementedError('Start & end of linesegment2D are equal')
        self._bounding_rectangle = None
        self.line = line
        if not line:
            self.line = volmdlr_curves.Line2D(start, end)
        LineSegment.__init__(self, start, end, self.line, name=name)

    def copy(self, deep=True, memo=None):
        """
        A specified copy of a LineSegment2D.
        """
        return self.__class__(start=self.start.copy(deep, memo), end=self.end.copy(deep, memo), name=self.name)

    def __hash__(self):
        return hash(('linesegment2d', self.start, self.end, self.line))

    def _data_hash(self):
        return self.start._data_hash() + self.end._data_hash()

    def _data_eq(self, other_object):
        if self.__class__.__name__ != other_object.__class__.__name__:
            return False
        return self.start == other_object.start and self.end == other_object.end

    def __eq__(self, other_object):
        if self.__class__.__name__ != other_object.__class__.__name__:
            return False
        return self.start == other_object.start and self.end == other_object.end

    def to_dict(self, *args, **kwargs):
        """Stores all Line Segment 2D attributes in a dict object."""
        return {'object_class': 'volmdlr.edges.LineSegment2D',
                'name': self.name,
                'start': self.start.to_dict(),
                'end': self.end.to_dict()
                }

    @property
    def bounding_rectangle(self):
        """
        Evaluates the bounding rectangle of the Line segment.
        """
        if not self._bounding_rectangle:
            self._bounding_rectangle = volmdlr.core.BoundingRectangle(
                min(self.start.x, self.end.x), max(self.start.x, self.end.x),
                min(self.start.y, self.end.y), max(self.start.y, self.end.y))
        return self._bounding_rectangle

    def straight_line_area(self):
        """
        Calculates the area of the LineSegment2D, with line drawn from start to end.

        :return: straight_line_area.
        """
        return 0.

    def straight_line_second_moment_area(self, *args, **kwargs):
        """Straight line second moment area for a line segment."""
        return 0, 0, 0

    def straight_line_center_of_mass(self):
        """Straight line center of mass."""
        return 0.5 * (self.start + self.end)

    def point_distance(self, point, return_other_point=False):
        """
        Computes the distance of a point to segment of line.

        :param point: point to calculate distance.
        :param return_other_point: Boolean variable to return line segment's corresponding point or not.
        """
        distance, point = volmdlr.linesegment2d_point_distance((self.start.x, self.start.y),
                                                               (self.end.x, self.end.y), (point.x, point.y))
        if return_other_point:
            return distance, volmdlr.Point2D(*point)
        return distance

    def point_projection(self, point):
        """
        If the projection falls outside the LineSegment2D, returns None.
        """
        point, curv_abs = volmdlr_curves.Line2D.point_projection(self.line, point)
        if curv_abs < 0 or curv_abs > self.length():
            if abs(curv_abs) < 1e-6 or math.isclose(curv_abs, self.length(),
                                                    abs_tol=1e-6):
                return point, curv_abs
            return None, curv_abs
        return point, curv_abs

    def line_intersections(self, line: volmdlr_curves.Line2D):
        """Line Segment intersections with volmdlr_curves.Line2D."""
        if self.direction_vector().is_colinear_to(line.direction_vector()):
            return []
        point = volmdlr.Point2D.line_intersection(self, line)
        if point is not None:
            point_projection1, _ = self.point_projection(point)
            intersections = [point_projection1]
            if point_projection1 is None:
                intersections = []

            elif line.__class__.__name__ == 'LineSegment2D':
                point_projection2, _ = line.point_projection(point)
                if point_projection2 is None:
                    intersections = []

            return intersections
        if line.point_belongs(self.start):
            return [self.start]
        if line.point_belongs(self.end):
            return [self.end]
        return []

    def linesegment_intersections(self, linesegment2d: 'LineSegment2D', abs_tol=1e-6):
        """
        Touching line segments does not intersect.
        """
        if self.bounding_rectangle.distance_to_b_rectangle(linesegment2d.bounding_rectangle) > abs_tol:
            return []
        if self.direction_vector(0.0).is_colinear_to(linesegment2d.direction_vector(0.0), abs_tol=abs_tol):
            return []
        point = volmdlr.Point2D.line_intersection(self, linesegment2d)
        # TODO: May be these commented conditions should be used for linesegment_crossings
        if point:  # and (point != self.start) and (point != self.end):
            point_projection1, _ = self.point_projection(point)
            if point_projection1 is None:
                return []

            point_projection2, _ = linesegment2d.point_projection(point)
            if point_projection2 is None:
                return []

            return [point_projection1]
        return []

    def line_crossings(self, line: 'volmdlr.curves.Line2D'):
        """Line Segment crossings with line 2d."""
        if self.direction_vector().is_colinear_to(line.direction_vector()):
            return []
        line_intersection = self.line_intersections(line)
        if line_intersection and (line_intersection[0].is_close(self.end) or
                                  line_intersection[0].is_close(self.start)):
            return []
        return line_intersection

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Plots the Linesegment2D.
        """
        width = edge_style.width

        if ax is None:
            _, ax = plt.subplots()

        p1, p2 = self.start, self.end
        if edge_style.arrow:
            if edge_style.plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=edge_style.color,
                        alpha=edge_style.alpha, style='o-')
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=edge_style.color,
                        alpha=edge_style.alpha)

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
            if edge_style.plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=edge_style.color,
                        marker='o', linewidth=width, alpha=edge_style.alpha)
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=edge_style.color,
                        linewidth=width, alpha=edge_style.alpha)
        return ax

    def to_3d(self, plane_origin, x1, x2):
        """
        Transforms the Line segment 2D into a 3D line segment.

        :param plane_origin: The origin of plane to draw the Line segment 3D.
        :type plane_origin: volmdlr.Point3D
        :param x1: First direction of the plane
        :type x1: volmdlr.Vector3D
        :param x2: Second direction of the plane.
        :type x2: volmdlr.Vector3D
        :return: A 3D line segment.
        :rtype: LineSegment3D
        """
        start = self.start.to_3d(plane_origin, x1, x2)
        end = self.end.to_3d(plane_origin, x1, x2)
        return LineSegment3D(start, end, name=self.name)

    def get_reverse(self):
        """
        Invert the sense of the line segment.
        """
        return LineSegment2D(self.end.copy(), self.start.copy())

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        LineSegment2D rotation.

        :param center: rotation center
        :param angle: angle rotation
        :return: a new rotated LineSegment2D
        """
        return LineSegment2D(self.start.rotation(center, angle), self.end.rotation(center, angle))

    def translation(self, offset: volmdlr.Vector2D):
        """
        LineSegment2D translation.

        :param offset: translation vector.
        :return: A new translated LineSegment2D.
        """
        return LineSegment2D(self.start.translation(offset), self.end.translation(offset))

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes vector frame_mapping and return a new LineSegment2D.

        side = 'old' or 'new'.
        """
        if side == 'old':
            new_start = frame.local_to_global_coordinates(self.start)
            new_end = frame.local_to_global_coordinates(self.end)
        elif side == 'new':
            new_start = frame.global_to_local_coordinates(self.start)
            new_end = frame.global_to_local_coordinates(self.end)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        return LineSegment2D(new_start, new_end)

    def plot_data(self, edge_style: plot_data.EdgeStyle = None):
        """
        Plot data method for a LineSegment2D.

        :param edge_style: edge style.
        :return: plot_data.LineSegment2D object.
        """
        return plot_data.LineSegment2D([self.start.x, self.start.y],
                                       [self.end.x, self.end.y],
                                       edge_style=edge_style)

    def create_tangent_circle(self, point, other_line):
        """Create a circle tangent to a LineSegment."""
        circle1, circle2 = other_line.create_tangent_circle(point, self.line)
        if circle1 is not None:
            _, curv_abs1 = self.line.point_projection(circle1.center)
            if curv_abs1 < 0. or curv_abs1 > self.length():
                circle1 = None
        if circle2 is not None:
            _, curv_abs2 = self.line.point_projection(circle2.center)
            if curv_abs2 < 0. or curv_abs2 > self.length():
                circle2 = None
        return circle1, circle2

    def infinite_primitive(self, offset):
        """Get an infinite primitive."""
        n = -self.unit_normal_vector()
        offset_point_1 = self.start + offset * n
        offset_point_2 = self.end + offset * n

        return volmdlr_curves.Line2D(offset_point_1, offset_point_2)

    def nearest_point_to(self, point):
        """
        Find out the nearest point on the linesegment to point.

        """

        points = self.discretization_points(number_points=500)
        return point.nearest_point(points)

    def axial_symmetry(self, line):
        """
        Finds out the symmetric linesegment2d according to a line.
        """

        points_symmetry = [point.axial_symmetry(line) for point in [self.start, self.end]]

        return self.__class__(points_symmetry[0], points_symmetry[1])

    def closest_point_on_segment(self, point):
        """Gets the closest point on the line segment and another given point."""
        segment_vector = self.direction_vector()
        p_vector = point - self.start
        p_vector = p_vector.to_vector()
        t_param = p_vector.dot(segment_vector) / segment_vector.dot(segment_vector)
        t_param = max(0, min(t_param, 1))
        closest_point = volmdlr.Point2D(self.start.x + t_param * segment_vector[0],
                                        self.start.y + t_param * segment_vector[1])
        return closest_point

    def distance_linesegment(self, linesegment, return_points=False):
        """
        Calculates the minimum distance between two line segments.

        :param linesegment: other line segment.
        :param return_points: boolean weather to return the minimum distance corresponding points or not.
        :return: minimum distance / minimal distance with corresponding points.
        """
        intersections = self.linesegment_intersections(linesegment)
        if intersections:
            if return_points:
                return 0.0, intersections[0], intersections[0]
            return 0.0
        closest_point_on_self_to_start = self.closest_point_on_segment(linesegment.start)
        closest_point_on_self_to_end = self.closest_point_on_segment(linesegment.end)
        closest_point_on_lineseg_to_start = linesegment.closest_point_on_segment(self.start)
        closest_point_on_lineseg_to_end = linesegment.closest_point_on_segment(self.end)
        min_dist, min_dist_point1, min_dist_point2 = math.inf, None, None
        for point1, point2 in zip([closest_point_on_self_to_start, closest_point_on_self_to_end,
                                   closest_point_on_lineseg_to_start, closest_point_on_lineseg_to_end],
                                  [linesegment.start, linesegment.end, self.start, self.end]):
            dist = point1.point_distance(point2)
            if dist < min_dist:
                min_dist = dist
                min_dist_point1, min_dist_point2 = point1, point2
        if return_points:
            return min_dist, min_dist_point1, min_dist_point2
        return min_dist


class ArcMixin:
    """
    Abstract class representing an arc.

    :param circle: arc related circle curve.
    :type circle: Union['volmdlr.curves.Circle2D', 'volmdlr.curves.Circle2D'].
    # :param start: The starting point
    # :type start: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
    # :param end: The finish point
    # :type end: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
    # :param name: The name of the arc. Default value is an empty string
    # :type name: str, optional
    """

    def __init__(self, circle, start, end, is_trigo: bool = True, name: str = ''):
        # Edge.__init__(self, start=start, end=end, name=name)
        self.start = start
        self.end = end
        self.circle = circle
        self.center = circle.center
        self.is_trigo = is_trigo
        self._length = None
        self.name = name

    def length(self):
        """
        Calculates the length of the Arc, with its radius, and its arc angle.

        :return: the length of the Arc.
        """
        if not self._length:
            self._length = self.circle.radius * abs(self.angle)
        return self._length

    def point_at_abscissa(self, abscissa):
        """
        Calculates a point in the Arc at a given abscissa.

        :param abscissa: abscissa where in the curve the point should be calculated.
        :return: Corresponding point.
        """
        if self.is_trigo:
            return self.start.rotation(self.circle.center, abscissa / self.circle.radius)
        return self.start.rotation(self.circle.center, -abscissa / self.circle.radius)

    def normal_vector(self, abscissa: float):
        """
        Get the normal vector of the Arc2D.

        :param abscissa: defines where in the Arc2D the
        normal vector is to be calculated
        :return: The normal vector of the Arc2D
        """
        point = self.point_at_abscissa(abscissa)
        normal_vector = self.circle.center - point
        normal_vector = normal_vector.to_vector()
        return normal_vector

    def direction_vector(self, abscissa: float):
        """
        Get direction vector of the Arc2D.

        :param abscissa: defines where in the Arc2D the
        direction vector is to be calculated
        :return: The direction vector of the Arc2D
        """
        return -self.normal_vector(abscissa=abscissa).normal_vector()

    def point_distance(self, point):
        """Returns the minimal distance to a point."""
        if self.point_belongs(point):
            return 0
        if self.circle.center.is_close(point):
            return self.circle.radius
        class_sufix = self.__class__.__name__[-2:]
        linesegment_class = getattr(sys.modules[__name__], 'LineSegment' + class_sufix)
        linesegment = linesegment_class(self.circle.center, point)
        if linesegment.length() > self.circle.radius:
            if self.linesegment_intersections(linesegment):
                return linesegment.length() - self.circle.radius
            return min(self.start.point_distance(point), self.end.point_distance(point))
        vector_to_point = point - self.circle.center
        vector_to_point = vector_to_point.unit_vector()
        projected_point = self.circle.center + self.circle.radius * vector_to_point
        if self.point_belongs(projected_point):
            return self.circle.radius - linesegment.length()
        return min(self.start.point_distance(point), self.end.point_distance(point))

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """
        Discretize an Edge to have "n" points.

        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Useful
            to mesh an arc
        :return: a list of sampled points
        """
        if not number_points:
            if not angle_resolution:
                number_points = 2
            else:
                number_points = max(math.ceil(self.angle * angle_resolution) + 1, 2)

        step = self.length() / (number_points - 1)
        return [self.point_at_abscissa(i * step)
                for i in range(number_points)]

    def get_geo_lines(self, tag: int, start_point_tag: int, center_point_tag: int, end_point_tag: int):
        """
        Gets the lines that define an Arc in a .geo file.

        :param tag: The linesegment index
        :type tag: int
        :param start_point_tag: The linesegment' start point index
        :type start_point_tag: int
        :param center_point_tag: The linesegment' center point index
        :type center_point_tag: int
        :param end_point_tag: The line segment's end point index
        :type end_point_tag: int

        :return: A line
        :rtype: str
        """

        return 'Circle(' + str(tag) + ') = {' + str(start_point_tag) + ', ' + \
            str(center_point_tag) + ', ' + str(end_point_tag) + '};'

    def get_geo_points(self):
        """
        Gets the points that define an Arc to use them in a .geo file.

        :return: A list of characteristic arc points
        :rtype: List

        """
        return [self.start, self.circle.center, self.end]

    def get_reverse(self):
        """
        Gets the reverse version of an arc.

        :return: An arc
        """

        return self.__class__(self.circle, start=self.end, end=self.start, is_trigo=not self.is_trigo)

    def split(self, split_point, tol: float = 1e-6):
        """
        Splits arc at a given point.

        :param split_point: splitting point.
        :param tol: tolerance.
        :return: list of two Arc.
        """
        if split_point.is_close(self.start, tol):
            return [None, self.copy()]
        if split_point.is_close(self.end, tol):
            return [self.copy(), None]
        if self.__class__.__name__[-2:] == "2D":
            return [self.__class__(self.circle, self.start, split_point, self.is_trigo),
                    self.__class__(self.circle, split_point, self.end, self.is_trigo)]
        return [self.__class__(self.circle, self.start, split_point),
                self.__class__(self.circle, split_point, self.end)]

    def get_shared_section(self, other_arc2, abs_tol: float = 1e-6):
        """
        Gets the shared section between two arcs.

        :param other_arc2: other arc to verify for shared section.
        :param abs_tol: tolerance.
        :return: shared arc section.
        """
        if self.__class__ != other_arc2.__class__:
            if self.__class__ == other_arc2.simplify.__class__:
                return self.get_shared_section(other_arc2.simplify, abs_tol)
            return []
        if not self.circle.center.is_close(other_arc2.circle.center) or self.circle.radius != self.circle.radius or \
                not any(self.point_belongs(point) for point in [other_arc2.start,
                                                                other_arc2.middle_point(), other_arc2.end]):
            return []
        if all(self.point_belongs(point, abs_tol) for point in
               [other_arc2.start, other_arc2.middle_point(), other_arc2.end]):
            return [other_arc2]
        if all(other_arc2.point_belongs(point, abs_tol) for point in
               [self.start, self.point_at_abscissa(self.length() * .5), self.end]):
            return [self]
        if self.point_belongs(other_arc2.start, abs_tol):
            arc1_, arc2_ = self.split(other_arc2.start, abs_tol)
        elif self.point_belongs(other_arc2.end, abs_tol):
            arc1_, arc2_ = self.split(other_arc2.end, abs_tol)
        else:
            raise NotImplementedError
        shared_arc_section = []
        for arc in [arc1_, arc2_]:
            if arc and all(other_arc2.point_belongs(point, abs_tol)
                           for point in [arc.start, arc.middle_point(), arc.end]):
                shared_arc_section.append(arc)
                break
        return shared_arc_section

    def delete_shared_section(self, other_arc2, abs_tol: float = 1e-6):
        """
        Deletes from self, the section shared with the other arc.

        :param other_arc2:
        :param abs_tol: tolerance.
        :return:
        """
        shared_section = self.get_shared_section(other_arc2, abs_tol)
        if not shared_section:
            return [self]
        if shared_section == self:
            return []
        split_arcs1 = self.split(shared_section[0].start)
        split_arcs2 = self.split(shared_section[0].end)
        new_arcs = []
        for arc in split_arcs1 + split_arcs2:
            if arc and not arc.point_belongs(shared_section[0].middle_point(), abs_tol):
                new_arcs.append(arc)
        return new_arcs

    def is_close(self, other_edge, tol: float = 1e-6):
        """
        Checks if two arc are the same considering the Euclidean distance.

        :param other_edge: other arc.
        :param tol: The tolerance under which the Euclidean distance is considered equal to 0, defaults to 1e-6
        :type tol: float, optional
        """

        if isinstance(other_edge, self.__class__):
            if (self.start.is_close(other_edge.start, tol) and self.end.is_close(other_edge.end, tol)
                    and self.circle.center.is_close(other_edge.circle.center, tol)
                    and self.point_belongs(other_edge.middle_point(), tol)):
                return True
        return False


class FullArcMixin(ArcMixin):
    """
    Abstract class for representing a circle with a start and end points that are the same.
    """

    def __init__(self, circle: Union[volmdlr.curves.Circle2D, volmdlr.curves.Circle3D],
                 start_end: Union[volmdlr.Point2D, volmdlr.Point3D], name: str = ''):
        self.circle = circle
        self.start_end = start_end
        ArcMixin.__init__(self, circle=circle, start=start_end, end=start_end, name=name)  # !!! this is dangerous

    @property
    def angle(self):
        """Angle of Full Arc. """
        return volmdlr.TWO_PI

    @property
    def periodic(self):
        """Return True if an edge is periodic."""
        return True

    def split(self, split_point, tol: float = 1e-6):
        """
        Splits arc at a given point.

        :param split_point: splitting point.
        :param tol: tolerance.
        :return: list of two Arc.
        """
        if split_point.is_close(self.start, tol):
            return [None, self.copy()]
        if split_point.is_close(self.end, tol):
            return [self.copy(), None]
        class_ = getattr(sys.modules[__name__], 'Arc' + self.__class__.__name__[-2:])
        return [class_(self.circle, self.start, split_point, self.is_trigo),
                class_(self.circle, split_point, self.end, self.is_trigo)]

    @classmethod
    def from_curve(cls, circle, name: str = ''):
        """Creates A full arc, 2d or 3d, from circle."""
        return cls(circle, circle.center + circle.frame.u * circle.radius, name=name)


class Arc2D(ArcMixin, Edge):
    """
    Class to draw Arc2D.

    angle: the angle measure always >= 0
    """

    def __init__(self, circle: 'volmdlr.curves.Circle2D',
                 start: volmdlr.Point2D,
                 end: volmdlr.Point2D,
                 is_trigo: bool = True,
                 name: str = ''):
        # self._center = center
        self.circle = circle
        self.is_trigo = is_trigo
        self._angle = None
        self._bounding_rectangle = None
        ArcMixin.__init__(self, circle, start, end, is_trigo, name=name)
        Edge.__init__(self, start=start, end=end, name=name)
        start_to_center = start - self.circle.center
        end_to_center = end - self.circle.center
        angle1 = math.atan2(start_to_center.y, start_to_center.x)
        angle2 = math.atan2(end_to_center.y, end_to_center.x)
        if self.is_trigo:
            self.angle1 = angle1
            self.angle2 = angle2
            if self.angle2 == 0.0:
                self.angle2 = volmdlr.TWO_PI
        else:
            self.angle1 = angle2
            self.angle2 = angle1

    def __hash__(self):
        return hash(('arc2d', self.circle, self.start, self.end, self.is_trigo))

    def __eq__(self, other_arc):
        if self.__class__.__name__ != other_arc.__class__.__name__:
            return False
        return (self.circle == other_arc.circle and self.start == other_arc.start
                and self.end == other_arc.end and self.is_trigo == other_arc.is_trigo)

    @classmethod
    def from_3_points(cls, point1, point2, point3, name: str = ''):
        """
        Creates a circle 2d from 3 points.

        :return: circle 2d.
        """
        circle = volmdlr_curves.Circle2D.from_3_points(point1, point2, point3)
        arc = cls(circle, point1, point3)
        if not arc.point_belongs(point2):
            return cls(circle, point1, point3, False, name=name)
        return arc

    @property
    def angle(self):
        """
        Returns the angle in radians of the arc.
        """
        if not self._angle:
            self._angle = self.get_angle()
        return self._angle

    def get_angle(self):
        """
        Gets arc angle.

        """
        clockwise_arc = self.reverse() if self.is_trigo else self
        vector_start = clockwise_arc.start - clockwise_arc.circle.center
        vector_end = clockwise_arc.end - clockwise_arc.circle.center
        arc_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_end)
        return arc_angle

    def _get_points(self):
        return [self.start, self.end]

    points = property(_get_points)

    def point_belongs(self, point, abs_tol=1e-6):
        """
        Check if a Point2D belongs to the Arc2D.

        """
        distance_point_to_center = point.point_distance(self.circle.center)
        if not math.isclose(distance_point_to_center, self.circle.radius, rel_tol=0.005):
            return False
        if point.is_close(self.start, abs_tol) or point.is_close(self.end, abs_tol):
            return True
        clockwise_arc = self.reverse() if self.is_trigo else self
        vector_start = clockwise_arc.start - clockwise_arc.circle.center
        vector_end = clockwise_arc.end - clockwise_arc.circle.center
        vector_point = point - clockwise_arc.circle.center
        arc_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_end)
        point_start_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_point)
        point_end_angle = volmdlr.geometry.clockwise_angle(vector_point, vector_end)
        if math.isclose(arc_angle, point_start_angle + point_end_angle, rel_tol=0.01):
            return True
        return False

    def to_full_arc_2d(self):
        """
        Convert to a full arc2d.
        """
        return FullArc2D(circle=self.circle, start_end=self.point_at_abscissa(0), name=self.name)

    def line_intersections(self, line2d: volmdlr_curves.Line2D):
        """
        Calculates the intersection between a line and an Arc2D.

        :param line2d: Line2D to verify intersections.
        :return: a list with intersections points.
        """
        full_arc_2d = self.to_full_arc_2d()
        fa2d_intersection_points = full_arc_2d.line_intersections(line2d)
        intersection_points = []
        for point in fa2d_intersection_points:
            if self.point_belongs(point):
                intersection_points.append(point)
        return intersection_points

    def linesegment_intersections(self, linesegment2d: LineSegment2D, abs_tol=1e-6):
        """
        Calculates the intersection between a LineSegment2D and an Arc2D.

        :param linesegment2d: LineSegment2D to verify intersections.
        :param abs_tol: tolerance.
        :return: a list with intersections points.
        """
        if self.bounding_rectangle.distance_to_b_rectangle(linesegment2d.bounding_rectangle) > abs_tol:
            return []
        full_arc_2d = self.to_full_arc_2d()
        fa2d_intersection_points = full_arc_2d.linesegment_intersections(linesegment2d, abs_tol)
        intersection_points = []
        for point in fa2d_intersection_points:
            if self.point_belongs(point, abs_tol):
                intersection_points.append(point)
        return intersection_points

    def bsplinecurve_intersections(self, bspline, abs_tol: float = 1e-6):
        """
        Intersections between an arc 2d and bspline curve 2d.

        :param bspline: bspline curve 2d.
        :param abs_tol: tolerance.
        :return: list of intersection points.
        """
        intersections = bspline.arc_intersections(self, abs_tol)
        return intersections

    def arc_intersections(self, arc, abs_tol: float = 1e-6):
        """Intersections between two arc 2d."""
        circle_intersections = vm_utils_intersections.get_circle_intersections(self.circle, arc.circle)
        arc_intersections = [inter for inter in circle_intersections if self.point_belongs(inter, abs_tol)]
        return arc_intersections

    def arcellipse_intersections(self, arcellipse, abs_tol: float = 1e-6):
        """
        Intersections between an arc 2d and arc-ellipse 2d.

        :param arcellipse: arc-ellipse 2d.
        :param abs_tol: tolerance
        :return: list of intersection points.
        """
        if self.bounding_rectangle.distance_to_b_rectangle(arcellipse.bounding_rectangle) > abs_tol:
            return []
        intersections = vm_utils_intersections.get_bsplinecurve_intersections(arcellipse, self, abs_tol)
        return intersections

    def abscissa(self, point: volmdlr.Point2D, tol=1e-6):
        """
        Returns the abscissa of a given point 2d.

        """
        if not math.isclose(point.point_distance(self.circle.center), self.circle.radius, abs_tol=tol):
            raise ValueError('Point not in arc')
        if point.point_distance(self.start) < tol:
            return 0
        if point.point_distance(self.end) < tol:
            return self.length()
        clockwise_arc = self.reverse() if self.is_trigo else self
        vector_start = clockwise_arc.start - clockwise_arc.circle.center
        vector_end = clockwise_arc.end - clockwise_arc.circle.center
        vector_point = point - clockwise_arc.circle.center
        arc_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_end)
        point_start_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_point)
        point_end_angle = volmdlr.geometry.clockwise_angle(vector_point, vector_end)
        if math.isclose(arc_angle, point_start_angle + point_end_angle, abs_tol=tol):
            if self.is_trigo:
                return self.length() - self.circle.radius * point_start_angle
            return self.circle.radius * point_start_angle
        raise ValueError('Point not in arc')

    def area(self):
        """
        Calculates the area of the Arc2D.

        :return: the area of the Arc2D.
        """
        return self.circle.radius ** 2 * self.angle / 2

    def center_of_mass(self):
        """
        Calculates the center of mass of the Arc2D.

        :return: center of mass point.
        """
        u = self.middle_point() - self.circle.center
        u = u.unit_vector()
        return self.circle.center + 4 / (3 * self.angle) * self.circle.radius * math.sin(
            self.angle * 0.5) * u

    @property
    def bounding_rectangle(self):
        """Gets the bounding rectangle for an Arc 2D."""
        if not self._bounding_rectangle:
            discretization_points = self.discretization_points(number_points=20)
            x_values, y_values = [], []
            for point in discretization_points:
                x_values.append(point.x)
                y_values.append(point.y)
            self._bounding_rectangle = volmdlr.core.BoundingRectangle(min(x_values), max(x_values),
                                                                      min(y_values), max(y_values))
        return self._bounding_rectangle

    def straight_line_area(self):
        """
        Calculates the area of the arc 2D, with line drawn from start to end.

        :return: straight_line_area.
        """
        if self.angle >= math.pi:
            angle = volmdlr.TWO_PI - self.angle
            area = math.pi * self.circle.radius ** 2 - 0.5 * self.circle.radius ** 2 * (
                    angle - math.sin(angle))
        else:
            angle = self.angle
            area = 0.5 * self.circle.radius ** 2 * (angle - math.sin(angle))

        if self.is_trigo:
            return area
        return -area

    def straight_line_second_moment_area(self, point: volmdlr.Point2D):
        """Straight line second moment area for an Arc 2D."""
        if self.angle2 < self.angle1:
            angle2 = self.angle2 + volmdlr.TWO_PI

        else:
            angle2 = self.angle2
        angle1 = self.angle1

        # Full arc section
        moment_area_x1 = self.circle.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle1) - math.sin(2 * angle2)))
        moment_area_y1 = self.circle.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle2) - math.sin(2 * angle1)))
        moment_area_xy1 = self.circle.radius ** 4 / 8 * (
                math.cos(angle1) ** 2 - math.cos(angle2) ** 2)

        # Triangle
        moment_area_x2, moment_area_y2, moment_area_xy2 = self._triangle_moment_inertia()
        if moment_area_x2 < 0.:
            moment_area_x2, moment_area_y2, moment_area_xy2 = -moment_area_x2, -moment_area_y2, -moment_area_xy2
        if self.angle < math.pi:
            if self.is_trigo:
                moment_area_x = moment_area_x1 - moment_area_x2
                moment_area_y = moment_area_y1 - moment_area_y2
                moment_area_xy = moment_area_xy1 - moment_area_xy2
            else:
                moment_area_x = moment_area_x2 - moment_area_x1
                moment_area_y = moment_area_y2 - moment_area_y1
                moment_area_xy = moment_area_xy2 - moment_area_xy1
        else:
            if self.is_trigo:
                moment_area_x = moment_area_x1 + moment_area_x2
                moment_area_y = moment_area_y1 + moment_area_y2
                moment_area_xy = moment_area_xy1 + moment_area_xy2
            else:
                moment_area_x = -moment_area_x2 - moment_area_x1
                moment_area_y = -moment_area_y2 - moment_area_y1
                moment_area_xy = -moment_area_xy2 - moment_area_xy1

        return volmdlr.geometry.huygens2d(moment_area_x, moment_area_y, moment_area_xy,
                                          self.straight_line_area(),
                                          self.circle.center,
                                          point)

    def _full_arc_moment_inertia(self, angle1, angle2):
        moment_inertia_x1 = self.circle.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle1) - math.sin(2 * angle2)))
        moment_inertia_y1 = self.circle.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle2) - math.sin(2 * angle1)))
        moment_inertia_xy1 = self.circle.radius ** 4 / 8 * (
                math.cos(angle1) ** 2 - math.cos(angle2) ** 2)
        return moment_inertia_x1, moment_inertia_y1, moment_inertia_xy1

    def _triangle_moment_inertia(self):
        xi, yi = self.start - self.circle.center
        xj, yj = self.end - self.circle.center
        moment_inertia_x2 = (yi ** 2 + yi * yj + yj ** 2) * (xi * yj - xj * yi) / 12.
        moment_inertia_y2 = (xi ** 2 + xi * xj + xj ** 2) * (xi * yj - xj * yi) / 12.
        moment_inertia_xy2 = (xi * yj + 2 * xi * yi + 2 * xj * yj + xj * yi) * (
                xi * yj - xj * yi) / 24.
        return moment_inertia_x2, moment_inertia_y2, moment_inertia_xy2

    def straight_line_center_of_mass(self):
        """Straight line center of mass."""
        if self.angle == math.pi:
            return self.center_of_mass()

        u = self.middle_point() - self.circle.center
        u = u.unit_vector()
        if self.angle >= math.pi:
            u = -u
        bissec = volmdlr_curves.Line2D(self.circle.center, self.circle.center + u)
        string = volmdlr_curves.Line2D(self.start, self.end)
        point = volmdlr.Point2D.line_intersection(bissec, string)
        a = point.point_distance(self.start)
        height = point.point_distance(self.circle.center)
        triangle_area = height * a
        triangle_cog = self.circle.center + 2 / 3. * height * u
        if self.angle < math.pi:
            cog = (
                          self.center_of_mass() * self.area() - triangle_area * triangle_cog) / abs(
                self.straight_line_area())
        else:
            cog = (
                          self.center_of_mass() * self.area() + triangle_area * triangle_cog) / abs(
                self.straight_line_area())

        return cog

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point: Point to be verified.
        :return: Return True if the point belongs to this surface, or False otherwise.
        """
        if self.point_belongs(point):
            return True
        if self.start == self.end:
            if point.point_distance(self.circle.center) <= self.circle.radius:
                return True
        center_distance_point = self.circle.center.point_distance(point)
        straight_line = LineSegment2D(self.start, self.end)
        for edge in [self, straight_line]:
            line_passing_trough_point = volmdlr_curves.Line2D(self.circle.center, point)
            straight_line_intersections = edge.line_intersections(line_passing_trough_point)
            if straight_line_intersections:
                if self.circle.center.point_distance(straight_line_intersections[0]) > center_distance_point:
                    return True
        return False

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plot arc 2d with Matplotlib."""
        if ax is None:
            _, ax = plt.subplots()

        if edge_style.plot_points:
            for point in [self.circle.center, self.circle.start, self.circle.end]:
                point.plot(ax=ax, color=edge_style.color, alpha=edge_style.alpha)

        ax.add_patch(matplotlib.patches.Arc((self.circle.center.x, self.circle.center.y), 2 * self.circle.radius,
                                            2 * self.circle.radius, angle=0,
                                            theta1=self.angle1 * 0.5 / math.pi * 360,
                                            theta2=self.angle2 * 0.5 / math.pi * 360,
                                            color=edge_style.color,
                                            alpha=edge_style.alpha))
        x_min, x_max = self.circle.center[0] - self.circle.radius*1.2, self.circle.center[0] + self.circle.radius*1.2
        y_min, y_max = self.circle.center[1] - self.circle.radius*1.2, self.circle.center[1] + self.circle.radius*1.2
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        return ax

    def to_3d(self, plane_origin, x, y):
        """
        Transforms the arc 2D into a 3D arc.

        :param plane_origin: The origin of plane to draw the arc 3D.
        :type plane_origin: volmdlr.Point3D
        :param x: First direction of the plane
        :type x: volmdlr.Vector3D
        :param y: Second direction of the plane.
        :type y: volmdlr.Vector3D
        :return: A 3D arc.
        :type: Arc3D.
        """
        circle3d = self.circle.to_3d(plane_origin, x, y)
        point_start = self.start.to_3d(plane_origin, x, y)
        point_interior = self.middle_point().to_3d(plane_origin, x, y)
        point_end = self.end.to_3d(plane_origin, x, y)
        arc = Arc3D(circle3d, point_start, point_end, name=self.name)
        if not arc.point_belongs(point_interior):
            circle3d = volmdlr_curves.Circle3D(volmdlr.Frame3D(
                circle3d.center, circle3d.frame.u, -circle3d.frame.v, circle3d.frame.u.cross(-circle3d.frame.v)),
                circle3d.radius)
            arc = Arc3D(circle3d, point_start, point_end, name=self.name)
        return arc

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Arc2D rotation.

        :param center: rotation center
        :param angle: angle rotation.
        :return: a new rotated Arc2D.
        """
        return Arc2D(*[point.rotation(center, angle) if point else point for point in
                       [self.circle, self.start, self.end]])

    def translation(self, offset: volmdlr.Vector2D):
        """
        Arc2D translation.

        :param offset: translation vector.
        :return: A new translated Arc2D.
        """
        return Arc2D(*[point.translation(offset) if point else point for point in
                       [self.circle, self.start, self.end]])

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes vector frame_mapping and return a new Arc2D.

        side = 'old' or 'new'
        """
        return Arc2D(self.circle.frame_mapping(frame, side), self.start.frame_mapping(frame, side),
                     self.end.frame_mapping(frame, side))

    def second_moment_area(self, point):
        """
        Second moment area of part of disk.

        """
        if self.angle2 < self.angle1:
            angle2 = self.angle2 + volmdlr.TWO_PI

        else:
            angle2 = self.angle2
        angle1 = self.angle1
        moment_area_x = self.circle.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle1) - math.sin(2 * angle2)))
        moment_area_y = self.circle.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle2) - math.sin(2 * angle1)))
        moment_area_xy = self.circle.radius ** 4 / 8 * (
                math.cos(angle1) ** 2 - math.cos(angle2) ** 2)

        # Must be computed at center, so huygens related to center
        return volmdlr.geometry.huygens2d(moment_area_x, moment_area_y, moment_area_xy, self.area(),
                                          self.circle.center, point)

    def plot_data(self, edge_style: plot_data.EdgeStyle = None, anticlockwise: bool = None):
        """
        Plot data method for a Arc2D.

        :param edge_style: edge style.
        :return: plot_data.Arc2D object.
        """
        list_node = self.discretization_points(number_points=20)
        data = []
        for node in list_node:
            data.append({'x': node.x, 'y': node.y})
        return plot_data.Arc2D(cx=self.circle.center.x,
                               cy=self.circle.center.y,
                               r=self.circle.radius,
                               start_angle=self.angle1,
                               end_angle=self.angle2,
                               edge_style=edge_style,
                               data=data,
                               anticlockwise=anticlockwise,
                               name=self.name)

    def copy(self, *args, **kwargs):
        """
        Creates and returns a deep copy of the Arc2D object.

        :param *args: Variable-length argument list.
        :param **kwargs: Arbitrary keyword arguments.
        :return: A new Arc2D object that is a deep copy of the original.

        """
        return Arc2D(self.circle.copy(), self.start.copy(), self.end.copy(), self.is_trigo)

    def cut_between_two_points(self, point1, point2):
        """
        Cuts Arc between two points, and return a new arc between these two points.
        """
        if (point1.is_close(self.start) and point2.is_close(self.end)) or \
                (point2.is_close(self.start) and point1.is_close(self.end)):
            return self
        raise NotImplementedError

    def infinite_primitive(self, offset):
        """Create an offset curve from a distance of the original curve."""
        vector_start_center = self.start - self.circle.center
        vector_start_center = vector_start_center.unit_vector()
        vector_end_center = self.end - self.circle.center
        vector_end_center = vector_end_center.unit_vector()
        if self.is_trigo:
            radius = self.circle.radius + offset
            center = self.circle.center
        else:
            radius = self.circle.radius - offset
            if radius < 0:
                return None
            center = self.circle.center
        new_circle = volmdlr_curves.Circle2D(center, radius)
        start = center + radius * vector_start_center
        end = center + radius * vector_end_center
        return Arc2D(new_circle, start, end, self.is_trigo)

    def complementary(self):
        """Gets the complementary Arc 2D. """
        return Arc2D(self.circle, self.end, self.start, self.is_trigo)

    def axial_symmetry(self, line):
        """ Finds out the symmetric arc 2D according to a line. """
        points_symmetry = [point.axial_symmetry(line) for point in [self.start, self.end]]

        return self.__class__(self.circle, start=points_symmetry[0],
                              end=points_symmetry[1], is_trigo=self.is_trigo)


class FullArc2D(FullArcMixin, Arc2D):
    """ An edge that starts at start_end, ends at the same point after having described a circle. """

    def __init__(self, circle: 'volmdlr.curves.Circle2D', start_end: volmdlr.Point2D,
                 name: str = ''):
        # self.interior = start_end.rotation(center, math.pi)
        self._bounding_rectangle = None
        FullArcMixin.__init__(self, circle=circle, start_end=start_end, name=name)
        Arc2D.__init__(self, circle=circle, start=start_end, end=start_end, name=name)
        self.angle1 = 0.0
        self.angle2 = volmdlr.TWO_PI

    def to_dict(self, use_pointers: bool = False, memo=None, path: str = '#', id_method=True, id_memo=None):
        """Stores all Full Arc 2D attributes in a dict object."""
        dict_ = self.base_dict()
        dict_['circle'] = self.circle.to_dict(use_pointers=use_pointers, memo=memo,
                                              id_method=id_method, id_memo=id_memo, path=path + '/circle')
        dict_['angle'] = self.angle
        dict_['is_trigo'] = self.is_trigo
        dict_['start_end'] = self.start.to_dict(use_pointers=use_pointers, memo=memo,
                                                id_method=id_method, id_memo=id_memo, path=path + '/start_end')
        return dict_

    def copy(self, *args, **kwargs):
        """Creates a copy of a fullarc 2d."""
        return FullArc2D(self.circle.copy(), self.start.copy())

    @classmethod
    def dict_to_object(cls, dict_, *args, **kwargs):
        """
        Create a FullArc2D object from a dictionary representation.

        This class method takes a dictionary containing the necessary data for
        creating a FullArc2D object and returns an instance of the FullArc2D class.
        It expects the dictionary to have the following keys:

        :param cls: The FullArc2D class itself (automatically passed).
        :param dict_: A dictionary containing the required data for object creation.
        :param args: Additional positional arguments (if any).
        :param kwargs: Additional keyword arguments (if any).

        :return: FullArc2D: An instance of the FullArc2D class created from the provided dictionary.
        """
        circle = volmdlr_curves.Circle2D.dict_to_object(dict_['circle'])
        start_end = volmdlr.Point2D.dict_to_object(dict_['start_end'])

        return cls(circle, start_end, name=dict_['name'])

    def __hash__(self):
        return hash((self.__class__.__name__, self.circle, self.start_end))

    def __eq__(self, other_arc):
        if self.__class__.__name__ != other_arc.__class__.__name__:
            return False
        return (self.circle == other_arc.circle) \
            and (self.start_end == other_arc.start_end)

    @property
    def bounding_rectangle(self):
        """Gets the bounding rectangle for a full arc 2d."""
        if not self._bounding_rectangle:
            self._bounding_rectangle = volmdlr.core.BoundingRectangle(
                self.circle.center.x - self.circle.radius, self.circle.center.x + self.circle.radius,
                self.circle.center.y - self.circle.radius, self.circle.center.y + self.circle.radius)
        return self._bounding_rectangle

    def straight_line_area(self):
        """
        Calculates the area of the full arc, with line drawn from start to end.

        :return: straight_line_area.
        """
        area = self.area()
        return area

    def center_of_mass(self):
        """Gets the center of the full arc 2d."""
        return self.circle.center

    def straight_line_center_of_mass(self):
        """Straight line center of mass."""
        return self.center_of_mass()

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point: Point to be verified.
        :return: Return True if the point belongs to this surface, or False otherwise.
        """
        if point.point_distance(self.circle.center) <= self.circle.radius:
            return True
        return False

    def to_3d(self, plane_origin, x, y):
        """
        Transforms the full arc 2D into a 3D full arc.

        :param plane_origin: The origin of plane to draw the full arc 3D.
        :type plane_origin: volmdlr.Point3D
        :param x: First direction of the plane
        :type x: volmdlr.Vector3D
        :param y: Second direction of the plane.
        :type y: volmdlr.Vector3D
        :return: A 3D full arc.
        :type: Full Arc 3D.
        """
        circle = self.circle.to_3d(plane_origin, x, y)
        start = self.start.to_3d(plane_origin, x, y)
        return FullArc3D(circle, start)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """Rotation of a full arc 2D."""
        new_circle = self.circle.rotation(center, angle)
        new_start_end = self.start.rotation(center, angle)
        return FullArc2D(new_circle, new_start_end)

    def translation(self, offset: volmdlr.Vector2D):
        """Translation of a full arc 2D."""
        new_circle = self.circle.translation(offset)
        new_start_end = self.start.translation(offset)
        return FullArc2D(new_circle, new_start_end)

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Map the 2D full arc to a new frame or its original frame.

        :param frame: The target frame for the mapping.
        :type frame: :class:`volmdlr.Frame2D`
        :param side: Specify whether to map the arc to the new frame ('new')
            or its original frame ('old').
        :type side: str
        :return: The full arc in the specified frame.
        :rtype: :class:`volmdlr.edges.FullArc2D`
        """
        return FullArc2D(*[point.frame_mapping(frame, side) for point in
                           [self.circle, self.start]])

    def polygonization(self):
        """Creates a Polygon from a full arc 2d."""
        return volmdlr.wires.ClosedPolygon2D(self.discretization_points(angle_resolution=15))

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plots a fullarc using Matplotlib."""
        return vm_common_operations.plot_circle(self.circle, ax, edge_style)

    def cut_between_two_points(self, point1, point2):
        """
        Cuts a full arc between two points on the fullarc.

        This method calculates the angles between the circle's center and the two points
        in order to determine the starting and ending angles of the arc. It then creates
        an Arc2D object representing the cut arc. If the original arc and the cut arc have
        opposite rotation directions, the cut arc is flipped to match the original arc's
        direction.

        :param point1: The first point defining the cut arc.
        :param point2: The second point defining the cut arc.

        :return: The cut arc between the two points.
        :rtype: Arc2D.
        """
        x1, y1 = point1 - self.circle.center
        x2, y2 = point2 - self.circle.center
        angle1 = math.atan2(y1, x1)
        angle2 = math.atan2(y2, x2)
        if angle2 < angle1:
            angle2 += volmdlr.TWO_PI
        arc = Arc2D(self.circle, point1, point2, self.is_trigo)
        if self.is_trigo != arc.is_trigo:
            arc = arc.complementary()
        return arc

    def line_intersections(self, line2d: volmdlr_curves.Line2D, tol=1e-9):
        """Full Arc 2D intersections with a Line 2D."""
        return self.circle.line_intersections(line2d, tol)

    def linesegment_intersections(self, linesegment2d: LineSegment2D, abs_tol=1e-9):
        """Full arc 2D intersections with a line segment."""
        return self.circle.linesegment_intersections(linesegment2d, abs_tol)

    def get_reverse(self):
        """Reverse of full arc 2D."""
        return self

    def point_belongs(self, point: volmdlr.Point2D, abs_tol: float = 1e-6):
        """
        Returns if given point belongs to the FullArc2D.
        """
        distance = point.point_distance(self.circle.center)
        return math.isclose(distance, self.circle.radius, abs_tol=abs_tol)


class ArcEllipse2D(Edge):
    """
    An 2-dimensional elliptical arc.

    :param ellipse: An ellipse curve, as base for the arc ellipse.
    :type ellipse: volmdlr.curves.Ellipse2D.
    :param start: The starting point of the elliptical arc
    :type start: :class:`volmdlr.Point2D`
    :param end: The end point of the elliptical arc
    :type end: :class:`volmdlr.Point2D`
    :param name: The name of the elliptical arc. Default value is ''
    :type name: str, optional
    """

    def __init__(self, ellipse: volmdlr_curves.Ellipse2D, start: volmdlr.Point2D,
                 end: volmdlr.Point2D, name: str = ''):
        Edge.__init__(self, start, end, name)
        self.ellipse = ellipse
        self.angle_start, self.angle_end = self.get_start_end_angles()
        self.angle = self.angle_end - self.angle_start
        self.center = ellipse.center
        self._bounding_rectangle = None
        self._reverse = None

    def get_start_end_angles(self):
        """Gets start and end angles for start and end points, respectively."""
        local_start_point = self.ellipse.frame.global_to_local_coordinates(self.start)
        u1, u2 = local_start_point.x / self.ellipse.major_axis, local_start_point.y / self.ellipse.minor_axis
        start_angle = volmdlr.geometry.sin_cos_angle(u1, u2)
        local_end_point = self.ellipse.frame.global_to_local_coordinates(self.end)
        u1, u2 = local_end_point.x / self.ellipse.major_axis, local_end_point.y / self.ellipse.minor_axis
        end_angle = volmdlr.geometry.sin_cos_angle(u1, u2)
        if self.ellipse.is_trigo and end_angle == 0.0:
            end_angle = volmdlr.TWO_PI
        return start_angle, end_angle

    @classmethod
    def from_3_points_and_center(cls, start, interior, end, center, name: str = ''):
        """
        Creates an arc ellipse using 3 points and a center.

        :param start: start point.
        :param interior: interior point.
        :param end: end point.
        :param center: ellipse's point.
        :param name: object's name.
        :return: An arc-ellipse2D object.
        """
        vector_center_start = start - center
        vector_center_end = end - center
        if vector_center_start.norm() >= vector_center_end.norm():
            x1 = start.x - center.x
            y1 = start.y - center.y
            x2 = end.x - center.x
            y2 = end.y - center.y
        else:
            x2 = start.x - center.x
            y2 = start.y - center.y
            x1 = end.x - center.x
            y1 = end.y - center.y
        if vector_center_start.is_colinear_to(vector_center_end) or abs(x1) == abs(x2):
            x2 = interior.x - center.x
            y2 = interior.y - center.y
            if abs(x1) == abs(x2):
                raise ValueError(f"Interior point{interior} is not valid. Try specifying another interior point.")
        minor_axis = math.sqrt((x1 ** 2 * y2 ** 2 - x2 ** 2 * y1 ** 2) / (x1 ** 2 - x2 ** 2))
        if abs(y1) != minor_axis:
            major_axis = math.sqrt(x1 ** 2 / (1 - (y1 ** 2 / minor_axis ** 2)))
        elif abs(y2) != minor_axis:
            major_axis = math.sqrt(x2 ** 2 / (1 - (y2 ** 2 / minor_axis ** 2)))
        else:
            raise NotImplementedError
        arcellipse = cls(volmdlr_curves.Ellipse2D(
            major_axis, minor_axis, volmdlr.Frame2D(center, volmdlr.X2D, volmdlr.Y2D)), start, end, name=name)
        if not arcellipse.point_belongs(interior):
            arcellipse = cls(volmdlr_curves.Ellipse2D(
                major_axis, minor_axis, volmdlr.Frame2D(center, volmdlr.X2D, -volmdlr.Y2D)), start, end, name=name)

        return arcellipse

    def _get_points(self):
        return self.discretization_points(number_points=20)

    points = property(_get_points)

    def length(self):
        """
        Calculates the length of the arc-ellipse 2d.

        :return: arc ellipse 2d's length
        """
        if not self._length:
            self._length = self.abscissa(self.end)
        return self._length

    def point_belongs(self, point, abs_tol: float = 1e-6):
        """
        Verifies if a point belongs to the arc ellipse 2d.

        :param point: point to be verified
        :param abs_tol: tolerance applied during calculations
        :return: True if the point belongs, False otherwise
        """
        if self.start.is_close(point, abs_tol) or self.end.is_close(point, abs_tol):
            return True
        point_in_local_coords = self.ellipse.frame.global_to_local_coordinates(point)
        local_ellipse = self.frame_mapping(self.ellipse.frame, 'new')
        if not math.isclose((point_in_local_coords.x -
                             local_ellipse.ellipse.center.x) ** 2 / local_ellipse.ellipse.major_axis ** 2 +
                            (point_in_local_coords.y -
                             local_ellipse.ellipse.center.y) ** 2 / local_ellipse.ellipse.minor_axis ** 2,
                            1, abs_tol=abs_tol) and\
                not math.isclose((point_in_local_coords.x -
                                  local_ellipse.ellipse.center.x) ** 2 / local_ellipse.ellipse.minor_axis ** 2 +
                                 (point_in_local_coords.y -
                                  local_ellipse.ellipse.center.y) ** 2 / local_ellipse.ellipse.major_axis ** 2,
                                 1, abs_tol=abs_tol):
            return False
        clockwise_arcellipse = self.reverse() if self.ellipse.is_trigo else self
        vector_start = clockwise_arcellipse.start - clockwise_arcellipse.ellipse.center
        vector_end = clockwise_arcellipse.end - clockwise_arcellipse.ellipse.center
        vector_point = point - clockwise_arcellipse.ellipse.center
        arc_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_end)
        point_start_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_point)
        point_end_angle = volmdlr.geometry.clockwise_angle(vector_point, vector_end)
        if math.isclose(arc_angle, point_start_angle + point_end_angle, abs_tol=1e-5):
            return True
        return False

    def valid_abscissa_start_end_angle(self, angle_abscissa):
        """Get valid abscissa angle for start and end."""
        angle_start = self.angle_start
        angle_end = angle_abscissa
        if self.angle_start > angle_abscissa >= self.angle_end:
            if angle_abscissa >= 0.0:
                angle_abscissa += 2 * math.pi
                angle_end = angle_abscissa
            else:
                angle_start = angle_abscissa
                angle_end = self.angle_start
        elif self.angle_start > self.angle_end >= angle_abscissa:
            angle_start = self.angle_start - 2 * math.pi
        return angle_start, angle_end

    def point_at_abscissa(self, abscissa):
        """Get a point at given abscissa."""
        if math.isclose(abscissa, 0.0, abs_tol=1e-6):
            return self.start
        if math.isclose(abscissa, self.length(), abs_tol=1e-6):
            return self.end
        if not self.ellipse.is_trigo:
            arc_ellipse_trigo = self.reverse()
            new_abscissa = self.length() - abscissa
            return arc_ellipse_trigo.point_at_abscissa(new_abscissa)
        discretized_points = self.discretization_points(number_points=100)
        aproximation_abscissa = 0
        aproximation_point = None
        for point1, point2 in zip(discretized_points[:-1], discretized_points[1:]):
            dist1 = point1.point_distance(point2)
            if aproximation_abscissa + dist1 > abscissa:
                aproximation_point = point1
                break
            aproximation_abscissa += dist1
        initial_point = self.ellipse.frame.global_to_local_coordinates(aproximation_point)
        u1, u2 = initial_point.x / self.ellipse.major_axis, initial_point.y / self.ellipse.minor_axis
        initial_angle = volmdlr.geometry.sin_cos_angle(u1, u2)
        angle_start, initial_angle = self.valid_abscissa_start_end_angle(initial_angle)
        abscissa_angle = vm_common_operations.ellipse_abscissa_angle_integration(
            self.ellipse, abscissa, angle_start, initial_angle)
        x = self.ellipse.major_axis * math.cos(abscissa_angle)
        y = self.ellipse.minor_axis * math.sin(abscissa_angle)
        return self.ellipse.frame.local_to_global_coordinates(volmdlr.Point2D(x, y))

    def abscissa(self, point: volmdlr.Point2D, tol: float = 1e-6):
        """
        Calculates the abscissa of a given point.

        :param point: point for calculating abscissa
        :param tol: tolerance.
        :return: a float, between 0 and the arc ellipse 2d's length
        """
        if self.start.is_close(point, tol):
            return 0.0
        if self.end.is_close(point, tol):
            if self._length:
                return self._length
            if not self.ellipse.is_trigo:
                arc_ellipse_trigo = self.reverse()
                abscissa_end = arc_ellipse_trigo.abscissa(self.start)
                return abscissa_end
        if self.point_belongs(point, 1e-4):
            if not self.ellipse.is_trigo:
                arc_ellipse_trigo = self.reverse()
                abscissa_point = arc_ellipse_trigo.abscissa(point)
                return self.length() - abscissa_point
            new_point = self.ellipse.frame.global_to_local_coordinates(point)
            u1, u2 = new_point.x / self.ellipse.major_axis, new_point.y / self.ellipse.minor_axis
            angle_abscissa = volmdlr.geometry.sin_cos_angle(u1, u2)
            if angle_abscissa == 0.0 and point.is_close(self.end):
                angle_abscissa = 2 * math.pi
            angle_start, angle_end = self.valid_abscissa_start_end_angle(angle_abscissa)

            def ellipse_arc_length(theta):
                return math.sqrt((self.ellipse.major_axis ** 2) * math.sin(theta) ** 2 +
                                 (self.ellipse.minor_axis ** 2) * math.cos(theta) ** 2)

            res, _ = scipy_integrate.quad(ellipse_arc_length, angle_start, angle_end)
            return res
        raise ValueError(f'point {point} does not belong to ellipse')

    @property
    def bounding_rectangle(self):
        """
        Calculates the bounding rectangle for the arc ellipse 2d.

        :return: Bounding Rectangle object.
        """
        if not self._bounding_rectangle:
            discretization_points = self.discretization_points(number_points=20)
            x_values, y_values = [], []
            for point in discretization_points:
                x_values.append(point.x)
                y_values.append(point.y)
            self._bounding_rectangle = volmdlr.core.BoundingRectangle(min(x_values), max(x_values),
                                                                      min(y_values), max(y_values))
        return self._bounding_rectangle

    def straight_line_area(self):
        """
        Calculates the area of the elliptic arc, with line drawn from start to end.

        :return: straight_line_area.
        """
        if self.angle >= math.pi:
            angle = volmdlr.TWO_PI - self.angle
            area = math.pi * self.ellipse.major_axis * self.ellipse.minor_axis -\
                0.5 * self.ellipse.major_axis * self.ellipse.minor_axis * (angle - math.sin(angle))
        else:
            angle = self.angle
            area = 0.5 * self.ellipse.major_axis * self.ellipse.minor_axis * (angle - math.sin(angle))

        if self.ellipse.is_trigo:
            return area
        return -area

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = None):
        """
        Discretization of an Edge to have "n" points.

        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned.
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Useful
            to mesh an arc.
        :return: a list of sampled points.
        """

        if not number_points:
            if not angle_resolution:
                number_points = 2
            else:
                number_points = math.ceil(angle_resolution * abs(self.angle / math.pi)) + 2
        if self.angle_start > self.angle_end:
            angle_end = self.angle_end + volmdlr.TWO_PI
            angle_start = self.angle_start
        elif self.angle_start == self.angle_end:
            angle_start = self.angle_start
            angle_end = angle_start + 2 * math.pi
        else:
            angle_end = self.angle_end
            angle_start = self.angle_start
        discretization_points = [self.ellipse.frame.local_to_global_coordinates(
            volmdlr.Point2D(self.ellipse.major_axis * math.cos(angle), self.ellipse.minor_axis * math.sin(angle)))
            for angle in npy.linspace(angle_start, angle_end, number_points)]
        return discretization_points

    def to_3d(self, plane_origin, x, y):
        """
        Transforms the arc of ellipse 2D into a 3D arc of ellipse.

        :param plane_origin: The origin of plane to draw the arc of ellipse 3D.
        :type plane_origin: volmdlr.Point3D
        :param x: First direction of the plane
        :type x: volmdlr.Vector3D
        :param y: Second direction of the plane.
        :type y: volmdlr.Vector3D
        :return: A 3D arc of ellipse.
        :type: ArcEllipse3D.
        """
        interior2d = self.point_at_abscissa(self.length() * 0.5)
        ellipse3d = self.ellipse.to_3d(plane_origin, x, y)
        start3d = self.start.to_3d(plane_origin, x, y)
        end3d = self.end.to_3d(plane_origin, x, y)
        interior3d = interior2d.to_3d(plane_origin, x, y)
        arcellipse = ArcEllipse3D(ellipse3d, start3d, end3d)
        if not arcellipse.point_belongs(interior3d):
            raise NotImplementedError
        return ArcEllipse3D(ellipse3d, start3d, end3d)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Plot arc-ellipse 2d using Matplotlib.

        :param ax: Matplotlib plot if there exists any.
        :param edge_style: edge styles.
        :return: Matplotlib plot
        """
        if ax is None:
            _, ax = plt.subplots()

        self.start.plot(ax=ax, color='r')
        self.end.plot(ax=ax, color='b')
        self.ellipse.center.plot(ax=ax, color='y')

        return vm_common_operations.plot_from_discretization_points(ax, edge_style, self, number_points=100)

    def normal_vector(self, abscissa):
        """
        Calculates the normal vector to an ellipse at a given abscissa.

        :param abscissa: The abscissa value at which the normal vector is to be calculated.
        :type abscissa: float.
        :return: The normal vector to the ellipse at the given abscissa.
        :rtype: volmdlr.Vector2D.

        :raises: ValueError If the abscissa is out of range.
        """
        tangent_vector = self.direction_vector(abscissa)
        return tangent_vector.normal_vector()

    def direction_vector(self, abscissa):
        """
        Calculates the tangent vector to an ellipse at a given abscissa.

        :param abscissa: The abscissa value at which the tangent vector is to be calculated.
        :type abscissa: float.
        :return: The tangent vector to the ellipse at the given abscissa.
        :rtype: volmdlr.Vector2D.

        :raises: ValueError If the abscissa is out of range.
        """
        point_at_abscissa = self.point_at_abscissa(abscissa)

        # Convert the point to local coordinates within the ellipse's frame
        point_at_abscissa_at_local_coord = self.ellipse.frame.global_to_local_coordinates(point_at_abscissa)

        # Calculate the slope of the tangent line at the given abscissa
        dy_dx = -(self.ellipse.minor_axis ** 2 * point_at_abscissa_at_local_coord.x) / (
                self.ellipse.major_axis ** 2 * point_at_abscissa_at_local_coord.y)

        # Construct the second point on the tangent line still on ellipse's frame.
        tangent_second_point = point_at_abscissa_at_local_coord + 1 * volmdlr.Point2D(1, dy_dx)

        # Convert the second point back to global coordinates
        global_coord_second_point = self.ellipse.frame.local_to_global_coordinates(tangent_second_point)

        tangent_vector = global_coord_second_point - point_at_abscissa
        tangent_vector = tangent_vector.to_vector()

        return tangent_vector

    def get_reverse(self):
        """Gets the arc ellipse in the reverse direction."""
        ellipse = self.ellipse.__class__(self.ellipse.major_axis, self.ellipse.minor_axis,
                                         volmdlr.Frame2D(self.ellipse.center, self.ellipse.frame.u,
                                                         -self.ellipse.frame.v))
        return self.__class__(ellipse, self.end.copy(), self.start.copy(), self.name + '_reverse')

    def line_intersections(self, line2d: volmdlr_curves.Line2D, tol: float = 1e-6):
        """
        Intersections between an Arc Ellipse 2D and a Line 2D.

        :param line2d: Line 2D to verify intersections
        :param tol: maximum tolerance.
        :return: List with all intersections
        """
        ellipse2d_linesegment_intersections = vm_utils_intersections.ellipse2d_line_intersections(
            self.ellipse, line2d)
        linesegment_intersections = []
        for inter in ellipse2d_linesegment_intersections:
            if self.point_belongs(inter, tol):
                linesegment_intersections.append(inter)
        return linesegment_intersections

    def linesegment_intersections(self, linesegment2d: LineSegment2D, abs_tol=1e-6):
        """
        Intersections between an Arc Ellipse 2D and a Line Segment 2D.

        :param linesegment2d: LineSegment 2D to verify intersections.
        :param abs_tol: tolerance.
        :return: List with all intersections.
        """
        if self.bounding_rectangle.distance_to_b_rectangle(linesegment2d.bounding_rectangle) > abs_tol:
            return []
        intersections = self.line_intersections(linesegment2d.line)
        linesegment_intersections = []
        for inter in intersections:
            if linesegment2d.point_belongs(inter, abs_tol):
                linesegment_intersections.append(inter)
        return linesegment_intersections

    def bsplinecurve_intersections(self, bspline, abs_tol: float = 1e-6):
        """
        Intersections between an Arc Ellipse 2D and a bSpline 2D.

        :param bspline: bspline 2D to verify intersections.
        :param abs_tol: tolerance.
        :return: List with all intersections.
        """
        if self.bounding_rectangle.distance_to_b_rectangle(bspline.bounding_rectangle) > abs_tol:
            return []
        intersections = vm_utils_intersections.get_bsplinecurve_intersections(self, bspline, abs_tol)
        return intersections

    def rotation(self, center, angle: float):
        """
        Rotation of ellipse around a center and an angle.

        :param center: center of the rotation.
        :param angle: angle to rotated of.
        :return: a rotated new ellipse.
        """
        return ArcEllipse2D(self.ellipse.rotation(center, angle), self.start.rotation(center, angle),
                            self.end.rotation(center, angle))

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and return a new Arc Ellipse 2D.

        side = 'old' or 'new'
        """
        return ArcEllipse2D(self.ellipse.frame_mapping(frame, side),
                            self.start.frame_mapping(frame, side),
                            self.end.frame_mapping(frame, side))

    def translation(self, offset: volmdlr.Vector2D):
        """
        Translates the Arc ellipse given an offset vector.

        :param offset: offset vector
        :return: new translated arc ellipse 2d.
        """
        return ArcEllipse2D(self.ellipse.translation(offset),
                            self.start.translation(offset),
                            self.end.translation(offset))

    def point_distance(self, point):
        """
        Calculates the distance from a given point to an Arc Ellipse 2d.

        :param point: point 2d.
        :return: distance.
        """
        return self.point_distance_to_edge(point)

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point: Point to be verified
        :return: Return True if the point belongs to this surface,
            or False otherwise
        """
        raise NotImplementedError(f'the straight_line_point_belongs method must be'
                                  f' overloaded by {self.__class__.__name__}')

    def split(self, split_point, tol: float = 1e-6):
        """
        Splits arc-ellipse at a given point.

        :param split_point: splitting point.
        :param tol: tolerance.
        :return: list of two Arc-Ellipse.
        """
        if split_point.is_close(self.start, tol):
            return [None, self.copy()]
        if split_point.is_close(self.end, tol):
            return [self.copy(), None]
        return [self.__class__(self.ellipse, self.start, split_point),
                self.__class__(self.ellipse, split_point, self.end)]

    def is_close(self, other_edge, tol: float = 1e-6):
        """
        Checks if two arc-ellipse are the same considering the Euclidean distance.

        :param other_edge: other arc-ellipse.
        :param tol: The tolerance under which the Euclidean distance is considered equal to 0, defaults to 1e-6.
        :type tol: float, optional
        """

        if isinstance(other_edge, self.__class__):
            if (self.start.is_close(other_edge.start, tol) and self.end.is_close(other_edge.end, tol)
                    and self.ellipse.center.is_close(other_edge.ellipse.center, tol) and
                    self.point_belongs(other_edge.point_at_abscissa(other_edge.length() * 0.5), tol)):
                return True
        return False

    def complementary(self):
        """Gets the complementary arc of ellipse."""
        return self.__class__(self.ellipse, self.end, self.start, name=self.name + '_complementary')


class FullArcEllipse(Edge):
    """
    Abstract class to define an ellipse.
    """

    def __init__(self, ellipse: Union[volmdlr_curves.Ellipse2D, volmdlr_curves.Ellipse3D],
                 start_end: Union[volmdlr.Point2D, volmdlr.Point3D], name: str = ''):
        self.start_end = start_end
        self.ellipse = ellipse
        self.is_trigo = True
        self.angle_start = 0.0
        self.center = ellipse.center
        self.angle_end = volmdlr.TWO_PI
        Edge.__init__(self, start=start_end, end=start_end, name=name)

    @property
    def periodic(self):
        """Returns True if edge is periodic."""
        return True

    def length(self):
        """
        Calculates the length of the ellipse.

        Ramanujan's approximation for the perimeter of the ellipse.
        P = math.pi * (a + b) [ 1 + (3h) / (10 + √(4 - 3h) ) ], where h = (a - b)**2/(a + b)**2.

        :return: Perimeter of the ellipse
        :rtype: float
        """
        return self.ellipse.length()

    def point_belongs(self, point: Union[volmdlr.Point2D, volmdlr.Point3D], abs_tol: float = 1e-6):
        """
        Verifies if a given point lies on the ellipse.

        :param point: point to be verified.
        :param abs_tol: Absolute tolerance to consider the point on the ellipse.
        :return: True is point lies on the ellipse, False otherwise
        """
        new_point = self.ellipse.frame.global_to_local_coordinates(point)
        return math.isclose(round(new_point.x ** 2 / self.ellipse.major_axis ** 2 +
                            new_point.y ** 2 / self.ellipse.minor_axis ** 2, 2), 1.0, abs_tol=abs_tol)

    def get_reverse(self):
        """
        Defines a new FullArcEllipse, identical to self, but in the opposite direction.

        """
        ellipse = self.ellipse.reverse()
        return self.__class__(ellipse, self.start_end)

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point: Point to be verified
        :return: Return True if the point belongs to this surface,
            or False otherwise
        """
        raise NotImplementedError(f'the straight_line_point_belongs method must be'
                                  f' overloaded by {self.__class__.__name__}')

    def abscissa(self, point, tol: float = 1e-6):
        """
        Computes the abscissa of an Edge.

        :param point: The point located on the edge.
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`].
        :param tol: The precision in terms of distance. Default value is 1e-4.
        :type tol: float, optional.
        :return: The abscissa of the point.
        :rtype: float
        """
        raise NotImplementedError(f'the abscissa method must be overloaded by {self.__class__.__name__}')

    @classmethod
    def from_curve(cls, ellipse, name: str = ''):
        """Creates a fullarc ellipse from a ellipse curve."""
        return cls(ellipse, ellipse.center + ellipse.frame.u * ellipse.major_axis, name=name)


class FullArcEllipse2D(FullArcEllipse, ArcEllipse2D):
    """
    Defines a FullArcEllipse2D.
    """

    def __init__(self, ellipse: volmdlr_curves.Ellipse2D, start_end: volmdlr.Point2D, name: str = ''):
        FullArcEllipse.__init__(self, ellipse, start_end, name)
        ArcEllipse2D.__init__(self, ellipse, start_end, start_end, name)
        self.theta = volmdlr.geometry.clockwise_angle(self.ellipse.major_dir, volmdlr.X2D)
        if self.theta == math.pi * 2:
            self.theta = 0.0
        self._bounding_rectangle = None

    def to_3d(self, plane_origin, x, y):
        """
        Transforms the full arc of ellipse 2D into a 3D full arc of ellipse.

        :param plane_origin: The origin of plane to draw the full arc of ellipse 3D.
        :type plane_origin: volmdlr.Point3D
        :param x: First direction of the plane
        :type x: volmdlr.Vector3D
        :param y: Second direction of the plane.
        :type y: volmdlr.Vector3D
        :return: A 3D full arc of ellipse.
        :rtype: FullArcEllipse3D
        """
        point_start_end3d = self.start_end.to_3d(plane_origin, x, y)
        ellipse = self.ellipse.to_3d(plane_origin, x, y)
        return FullArcEllipse3D(ellipse, point_start_end3d, name=self.name + "_3D")

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and return a new FullArcEllipse2D.

        :param frame: Local coordinate system.
        :type frame: volmdlr.Frame2D
        :param side: 'old' will perform a transformation from local to global coordinates. 'new' will
            perform a transformation from global to local coordinates.
        :type side: str
        :return: A new transformed FulLArcEllipse2D.
        :rtype: FullArcEllipse2D
        """
        return FullArcEllipse2D(self.ellipse.frame_mapping(frame, side),
                                self.start_end.frame_mapping(frame, side))

    def translation(self, offset: volmdlr.Vector2D):
        """
        Full ArcEllipse 2D translation.

        :param offset: translation vector.
        :type offset: volmdlr.Vector2D
        :return: A new translated FullArcEllipse2D.
        :rtype: FullArcEllipse2D
        """
        return FullArcEllipse2D(self.ellipse.translation(offset), self.start_end.translation(offset), self.name)

    def abscissa(self, point: Union[volmdlr.Point2D, volmdlr.Point3D], tol: float = 1e-6):
        """
        Calculates the abscissa of a given point.

        :param point: point for calculating abscissa.
        :param tol: tolerance.
        :return: a float, between 0 and the ellipse's length.
        """
        return self.ellipse.abscissa(point, tol)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Matplotlib plot for an ellipse.

        """
        if ax is None:
            _, ax = plt.subplots()
        ax = vm_common_operations.plot_from_discretization_points(
            ax, edge_style=edge_style, element=self, number_points=50)
        if edge_style.equal_aspect:
            ax.set_aspect('equal')
        return ax


class LineSegment3D(LineSegment):
    """
    Define a line segment limited by two points.

    """

    def __init__(self, start: volmdlr.Point3D, end: volmdlr.Point3D, line: volmdlr_curves.Line3D = None,
                 name: str = ''):
        if start.is_close(end):
            raise NotImplementedError('Start and end of Linesegment3D are equal')
        self.line = line
        if not line:
            self.line = volmdlr_curves.Line3D(start, end)
        LineSegment.__init__(self, start=start, end=end, line=self.line, name=name)
        self._bbox = None

    @property
    def bounding_box(self):
        """Gets bounding box for Line Segment 3D. """
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        """Sets new value to Line Segment's Bounding box."""
        self._bbox = new_bounding_box

    def __hash__(self):
        return hash((self.__class__.__name__, self.start, self.end))

    def __eq__(self, other_linesegment3d):
        if other_linesegment3d.__class__ != self.__class__:
            return False
        return (self.start == other_linesegment3d.start
                and self.end == other_linesegment3d.end)

    def _bounding_box(self):
        """
        Calculates the bounding box for a line segment 3D.

        :return: Bounding box for line segment 3d.
        """

        xmin = min(self.start.x, self.end.x)
        xmax = max(self.start.x, self.end.x)
        ymin = min(self.start.y, self.end.y)
        ymax = max(self.start.y, self.end.y)
        zmin = min(self.start.z, self.end.z)
        zmax = max(self.start.z, self.end.z)

        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def to_dict(self, *args, **kwargs):
        """Stores all Line Segment 3D in a dict object."""
        return {'object_class': 'volmdlr.edges.LineSegment3D',
                'name': self.name,
                'start': self.start.to_dict(),
                'end': self.end.to_dict()
                }

    def normal_vector(self, abscissa=0.):
        """
        Returns the normal vector to the curve at the specified abscissa.
        """
        direction_vector = self.direction_vector()
        return direction_vector.deterministic_normal_vector()

    def unit_normal_vector(self, abscissa=0.):
        """Calculates the Line segment's unit normal vector."""
        return self.normal_vector().unit_vector()

    def point_distance(self, point):
        """Returns the minimal distance to a point."""
        distance, point = volmdlr.linesegment3d_point_distance((self.start.x, self.start.y, self.start.z),
                                                               (self.end.x, self.end.y, self.end.z),
                                                               (point.x, point.y, point.z))
        return distance

    def plane_projection2d(self, center, x, y):
        """
        Calculates the projection of a line segment 3d on to a plane.

        :param center: plane center.
        :param x: plane u direction.
        :param y: plane v direction.
        :return: line segment 3d.
        """
        start, end = self.start.plane_projection2d(center, x, y), self.end.plane_projection2d(center, x, y)
        if not start.is_close(end):
            return LineSegment2D(start, end)
        return None

    def line_intersections(self, line, tol: float = 1e-6):
        """
        Gets the intersection between a line segment 3d and line3D.

        :param line: other line.
        :param tol: maximum tolerance.
        :return: a list with the intersection points.
        """
        line_self = self.line
        if line_self.skew_to(line):
            return []
        intersection = line_self.intersection(line, tol=tol)
        if intersection and self.point_belongs(intersection):
            return [intersection]
        if line.point_belongs(self.start):
            return [self.start]
        if line.point_belongs(self.end):
            return [self.end]
        return []

    def linesegment_intersections(self, linesegment, abs_tol: float = 1e-6):
        """
        Gets the intersection between a line segment 3d and another line segment 3D.

        :param linesegment: other line segment.
        :param abs_tol: tolerance.
        :return: a list with the intersection points.
        """
        intersection = self.line.intersection(linesegment.line)
        if intersection and self.point_belongs(intersection, abs_tol=abs_tol) and\
                linesegment.point_belongs(intersection, abs_tol=abs_tol):
            return [intersection]
        return []

    def rotation(self, center: volmdlr.Point3D,
                 axis: volmdlr.Vector3D, angle: float):
        """
        LineSegment3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated LineSegment3D
        """
        start = self.start.rotation(center, axis, angle)
        end = self.end.rotation(center, axis, angle)
        return LineSegment3D(start, end)

    def __contains__(self, point):

        point1, point2 = self.start, self.end
        axis = point2 - point1
        test = point.rotation(point1, axis, math.pi)
        if test.is_close(point):
            return True

        return False

    def translation(self, offset: volmdlr.Vector3D):
        """
        LineSegment3D translation.

        :param offset: translation vector
        :return: A new translated LineSegment3D
        """
        return LineSegment3D(
            self.start.translation(offset), self.end.translation(offset))

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes LineSegment3D frame_mapping and return a new LineSegment3D.

        side = 'old' or 'new'
        """
        if side == 'old':
            return LineSegment3D(
                *[frame.local_to_global_coordinates(point) for point in [self.start, self.end]])
        if side == 'new':
            return LineSegment3D(
                *[frame.global_to_local_coordinates(point) for point in [self.start, self.end]])
        raise ValueError('Please Enter a valid side: old or new')

    def copy(self, *args, **kwargs):
        """Returns a copy of the line segment."""
        return LineSegment3D(self.start.copy(), self.end.copy())

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plots the Line segment 3d using matplotlib."""
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        points = [self.start, self.end]
        x = [point.x for point in points]
        y = [point.y for point in points]
        z = [point.z for point in points]
        if edge_style.edge_ends:
            ax.plot(x, y, z, color=edge_style.color, alpha=edge_style.alpha, marker='o')
        else:
            ax.plot(x, y, z, color=edge_style.color, alpha=edge_style.alpha)
        if edge_style.edge_direction:
            x, y, z = self.point_at_abscissa(0.5 * self.length())
            u, v, w = 0.05 * self.direction_vector()
            ax.quiver(x, y, z, u, v, w, length=self.length() / 100,
                      arrow_length_ratio=5, normalize=True,
                      pivot='tip', color=edge_style.color)
        return ax

    def plot2d(self, x_3d, y_3d, ax=None, color='k', width=None):
        """Creates a 2d plot of the Line segment 3d using matplotlib."""
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        edge2d = self.plane_projection2d(volmdlr.O3D, x_3d, y_3d)
        edge2d.plot(ax=ax, edge_style=EdgeStyle(color=color, width=width))
        return ax

    def plot_data(self, x_3d, y_3d, edge_style=plot_data.EdgeStyle(color_stroke=plot_data.colors.BLACK,
                                                                   line_width=1, dashline=None)):
        """Plot a Line Segment 3D object using dessia's plot_data library."""
        edge2d = self.plane_projection2d(volmdlr.O3D, x_3d, y_3d)
        return edge2d.plot_data(edge_style)

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a LineSegment3D into an LineSegment2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: LineSegment2D.
        """
        p2d = [point.to_2d(plane_origin, x, y) for point in (self.start, self.end)]
        if p2d[0].is_close(p2d[1]):
            return None
        return LineSegment2D(*p2d, name=self.name)

    def to_bspline_curve(self, resolution=10):
        """
        Convert a LineSegment3D to a BSplineCurve3D.
        """
        degree = 1
        points = [self.point_at_abscissa(abscissa / self.length())
                  for abscissa in range(resolution + 1)]
        bspline_curve = BSplineCurve3D.from_points_interpolation(points, degree)
        return bspline_curve

    def get_reverse(self):
        """
        Gets the reverse of the Line Segment.
        """
        return LineSegment3D(self.end.copy(), self.start.copy())

    def minimum_distance_points(self, other_line):
        """
        Returns the points on this line and on the other line that are the closest of lines.
        """
        u = self.end - self.start
        v = other_line.end - other_line.start
        w = self.start - other_line.start
        u_dot_u = u.dot(u)
        u_dot_v = u.dot(v)
        v_dot_v = v.dot(v)
        u_dot_w = u.dot(w)
        v_dot_w = v.dot(w)
        if (u_dot_u * v_dot_v - u_dot_v ** 2) != 0:
            s_param = (u_dot_v * v_dot_w - v_dot_v * u_dot_w) / (u_dot_u * v_dot_v - u_dot_v ** 2)
            t_param = (u_dot_u * v_dot_w - u_dot_v * u_dot_w) / (u_dot_u * v_dot_v - u_dot_v ** 2)
            point1 = self.start + s_param * u
            point2 = other_line.start + t_param * v
            return point1, point2
        return self.start, other_line.start

    def matrix_distance(self, other_line):
        """
        Gets the points corresponding to the distance between to lines using matrix distance.

        :param other_line: Other line.
        :return: Two points corresponding to the distance between to lines.
        """
        return volmdlr.core_compiled.LineSegment3DDistance([self.start, self.end], [other_line.start, other_line.end])

    def parallel_distance(self, other_linesegment):
        """Calculates the parallel distance between two Line Segments 3D."""
        pt_a, pt_b, pt_c = self.start, self.end, other_linesegment.start
        vector = volmdlr.Vector3D((pt_a - pt_b).vector)
        vector = vector.unit_vector()
        plane1 = volmdlr.surfaces.Plane3D.from_3_points(pt_a, pt_b, pt_c)
        v = vector.cross(plane1.frame.w)  # distance vector
        # pt_a = k*u + c*v + pt_c
        res = (pt_a - pt_c).vector
        x, y, z = res[0], res[1], res[2]
        u1, u2, u3 = vector.x, vector.y, vector.z
        v1, v2, v3 = v.x, v.y, v.z

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
        raise NotImplementedError

    def distance_linesegment(self, linesegment, return_points=False):
        """
        Calculates the minimum distance between two line segments in 3d.

        :param linesegment: other line segment.
        :param return_points: boolean weather to return the minimum distance corresponding points or not.
        :return: minimum distance / minimal distance with corresponding points.
        """
        p1, p2 = self.minimum_distance_points(linesegment)
        if not self.point_belongs(p1):
            p1 = self.start if self.start.point_distance(p1) < self.end.point_distance(p1) else self.end
        if not linesegment.point_belongs(p2):
            p2 = linesegment.start if linesegment.start.point_distance(p2) <\
                                      linesegment.end.point_distance(p2) else linesegment.end
        if return_points:
            return p1.point_distance(p2), p1, p2
        return p1.point_distance(p2)

    def distance_arc(self, arc3d, return_points=False):
        """
        Calculates the minimum distance between a line segment and an arc in 3d.

        :param arc3d: other line segment.
        :param return_points: boolean weather to return the minimum distance corresponding points or not.
        :return: minimum distance / minimal distance with corresponding points.
        """
        return arc3d.distance_linesegment(self, return_points)

    def extrusion(self, extrusion_vector):
        """
        Extrusion of a Line Segment 3D, in a specific extrusion direction.

        :param extrusion_vector: the extrusion vector used.
        :return: An extruded Plane Face 3D.
        """
        u = self.unit_direction_vector()
        v = extrusion_vector.copy()
        v = v.unit_vector()
        w = u.cross(v)
        length_1 = self.length()
        length_2 = extrusion_vector.norm()
        plane = volmdlr.surfaces.Plane3D(volmdlr.Frame3D(self.start, u, v, w))
        return [volmdlr.faces.PlaneFace3D.from_surface_rectangular_cut(plane, 0, length_1, 0, length_2)]

    def _conical_revolution(self, params):
        """Creates a conical revolution of a Line Segment 3D."""
        axis, u, dist1, dist2, angle, cone_origin = params
        v = axis.cross(u)
        direction_vector = self.direction_vector()
        direction_vector = direction_vector.unit_vector()

        semi_angle = math.atan2(direction_vector.dot(u), direction_vector.dot(axis))
        if semi_angle > 0.5 * math.pi:
            semi_angle = math.pi - semi_angle
            cone_frame = volmdlr.Frame3D(cone_origin, u, -v, -axis)
            angle2 = - angle
        else:
            angle2 = angle
            cone_frame = volmdlr.Frame3D(cone_origin, u, v, axis)

        surface = volmdlr.surfaces.ConicalSurface3D(cone_frame, semi_angle)
        return [volmdlr.faces.ConicalFace3D.from_surface_rectangular_cut(
            surface, 0, angle2, z1=dist1 / math.tan(semi_angle), z2=dist2 / math.tan(semi_angle))]

    def _cylindrical_revolution(self, params):
        """Creates a cylindrical revolution of a Line Segment 3D."""
        axis, u, p1_proj, dist1, _, angle = params
        v = axis.cross(u)
        surface = volmdlr.surfaces.CylindricalSurface3D(volmdlr.Frame3D(p1_proj, u, v, axis), dist1)
        return [volmdlr.faces.CylindricalFace3D.from_surface_rectangular_cut(
            surface, 0, angle, 0, (self.end - self.start).dot(axis))]

    def _plane_revolution(self, params):
        """
        Creates Plane Revolution of a Line Segment 3D.

        :param params: needed parameters.
        :return: List of plane revolution faces.
        """
        axis, angle, p1_proj, u, distance_1, distance_2, line_intersection = params
        v = axis.cross(u)
        surface = volmdlr.surfaces.Plane3D(
            volmdlr.Frame3D(p1_proj, u, v, axis))
        if self.point_belongs(line_intersection):  # Linesegment intersects revolution axis
            faces = []
            for i, radius in enumerate([distance_1, distance_2]):
                if math.isclose(radius, 0, abs_tol=1e-9):
                    continue
                if i == 0:
                    arc_point1 = volmdlr.O2D + volmdlr.X2D * radius
                else:
                    arc_point1 = volmdlr.O2D - volmdlr.X2D * radius
                arc_point2 = arc_point1.rotation(volmdlr.O2D, angle / 2)
                arc_point3 = arc_point1.rotation(volmdlr.O2D, angle)
                arc = Arc2D.from_3_points(arc_point1, arc_point2, arc_point3)
                outer_contour = volmdlr.wires.Contour2D([LineSegment2D(volmdlr.O2D, arc_point1), arc,
                                                         LineSegment2D(arc_point3, volmdlr.O2D)])
                face = volmdlr.faces.PlaneFace3D(surface, volmdlr.surfaces.Surface2D(outer_contour, []))
                faces.append(face)
            return faces
        smaller_r, bigger_r = sorted([distance_1, distance_2])
        inner_contours2d = []
        if angle == volmdlr.TWO_PI:
            # Only 2 circles as contours
            bigger_circle = volmdlr_curves.Circle2D(volmdlr.O2D, bigger_r)
            outer_contour2d = volmdlr.wires.Contour2D(
                bigger_circle.split_at_abscissa(bigger_circle.length() * 0.5))
            if not math.isclose(smaller_r, 0, abs_tol=1e-9):
                smaller_circle = volmdlr_curves.Circle2D(volmdlr.O2D, smaller_r)
                inner_contours2d = [volmdlr.wires.Contour2D(
                    smaller_circle.split_at_abscissa(smaller_circle.length() * 0.5))]
            return [volmdlr.faces.PlaneFace3D(surface,
                                              volmdlr.surfaces.Surface2D(outer_contour2d, inner_contours2d))]
        # Two arcs and lines
        arc1_s = volmdlr.Point2D(bigger_r, 0)
        arc1_i = arc1_s.rotation(center=volmdlr.O2D,
                                 angle=0.5 * angle)
        arc1_e = arc1_s.rotation(center=volmdlr.O2D, angle=angle)
        arc1 = Arc2D.from_3_points(arc1_s, arc1_i, arc1_e)

        arc2_e = volmdlr.Point2D(smaller_r, 0)
        arc2_i = arc2_e.rotation(center=volmdlr.O2D,
                                 angle=0.5 * angle)
        arc2_s = arc2_e.rotation(center=volmdlr.O2D, angle=angle)
        arc2 = Arc2D.from_3_points(arc2_s, arc2_i, arc2_e)

        line1 = LineSegment2D(arc1_e, arc2_s)
        line2 = LineSegment2D(arc2_e, arc1_s)

        outer_contour2d = volmdlr.wires.Contour2D([arc1, line1,
                                                   arc2, line2])

        return [volmdlr.faces.PlaneFace3D(surface, volmdlr.surfaces.Surface2D(outer_contour2d, inner_contours2d))]

    def revolution(self, axis_point, axis, angle):
        """
        Returns the face generated by the revolution of the line segments.
        """
        axis_line3d = volmdlr_curves.Line3D(axis_point, axis_point + axis)
        if axis_line3d.point_belongs(self.start) and axis_line3d.point_belongs(
                self.end):
            return []
        line_intersection = self.line.intersection(axis_line3d)
        p1_proj, _ = axis_line3d.point_projection(self.start)
        p2_proj, _ = axis_line3d.point_projection(self.end)
        distance_1 = self.start.point_distance(p1_proj)
        distance_2 = self.end.point_distance(p2_proj)
        if not math.isclose(distance_1, 0., abs_tol=1e-9):
            u = self.start - p1_proj  # Unit vector from p1_proj to p1
            u = u.unit_vector()
        elif not math.isclose(distance_2, 0., abs_tol=1e-9):
            u = self.end - p2_proj  # Unit vector from p1_proj to p1
            u = u.unit_vector()
        else:
            return []
        if u.is_colinear_to(self.direction_vector()):
            # Planar face
            return self._plane_revolution([axis, angle, p1_proj, u, distance_1, distance_2, line_intersection])
        if line_intersection and self.point_belongs(line_intersection):
            if not math.isclose(distance_1, 0., abs_tol=1e-9) and not math.isclose(distance_2, 0., abs_tol=1e-9):
                u1 = self.start - p1_proj  # Unit vector from p1_proj to p1
                u1 = u1.unit_vector()
                u2 = self.end - p2_proj  # Unit vector from p1_proj to p1
                u2 = u2.unit_vector()
                faces = self._conical_revolution([axis, u1, 0, distance_1, angle, line_intersection]) + \
                        self._conical_revolution([axis, u2, 0, distance_2, angle, line_intersection])
                return faces
        if not math.isclose(distance_1, distance_2, abs_tol=1e-9):
            # Conical
            return self._conical_revolution([axis, u, distance_1, distance_2, angle, line_intersection])
        # Cylindrical face
        return self._cylindrical_revolution([axis, u, p1_proj, distance_1, distance_2, angle])

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D):
        """Trims a Line Segment at two given points."""
        if not self.point_belongs(point1) or not self.point_belongs(point2):
            raise ValueError('Point not on curve')

        return LineSegment3D(point1, point2)

    def sweep(self, *args):
        """
        Line Segment 3D is used as path for sweeping given section through it.

        :return:
        """
        section_contour2d, section_contour3d = args
        if section_contour3d is None:
            start_tangent = self.unit_direction_vector(0.)
            normal = self.unit_normal_vector(0.)
            if normal is None:
                normal = start_tangent.deterministic_unit_normal_vector()
            tangent_normal_orthonormal = start_tangent.cross(normal)
            section_contour3d = section_contour2d.to_3d(self.start, normal, tangent_normal_orthonormal)
        new_faces = []
        for contour_primitive in section_contour3d.primitives:
            new_faces.extend(contour_primitive.extrusion(self.length()
                                                         * self.unit_direction_vector()))
        return new_faces


class BSplineCurve3D(BSplineCurve):
    """
    A class for 3-dimensional B-spline curves.

    The following rule must be respected : `number of knots = number of control points + degree + 1`

    :param degree: The degree of the 3-dimensional B-spline curve
    :type degree: int
    :param control_points: A list of 3-dimensional points
    :type control_points: List[:class:`volmdlr.Point3D`]
    :param knot_multiplicities: The vector of multiplicities for each knot
    :type knot_multiplicities: List[int]
    :param knots: The knot vector composed of values between 0 and 1
    :type knots: List[float]
    :param weights: The weight vector applied to the knot vector. Default
        value is None
    :type weights: List[float], optional
    :param periodic: If `True` the B-spline curve is periodic. Default value
        is False
    :type periodic: bool, optional
    :param name: The name of the B-spline curve. Default value is ''
    :type name: str, optional
    """

    def __init__(self,
                 degree: int,
                 control_points: List[volmdlr.Point3D],
                 knot_multiplicities: List[int],
                 knots: List[float],
                 weights: List[float] = None,
                 name: str = ''):

        BSplineCurve.__init__(self, degree,
                              control_points,
                              knot_multiplicities,
                              knots,
                              weights,
                              name)

        self._bbox = None

    @property
    def bounding_box(self):
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def _bounding_box(self):
        """Creates a bounding box from the bspline points."""
        return volmdlr.core.BoundingBox.from_points(self.discretization_points())

    def look_up_table(self, resolution: int = 20, start_parameter: float = 0,
                      end_parameter: float = 1):
        """
        Creates a table of equivalence between parameter t (evaluation of BSplineCurve) and the cumulative distance.

        :param resolution: The precision of the table. Auto-adjusted by the
            algorithm. Default value set to 20
        :type resolution: int, optional
        :param start_parameter: First parameter evaluated in the table.
            Default value set to 0
        :type start_parameter: float, optional
        :param end_parameter: Last parameter evaluated in the table.
            Default value set to 1
        :type start_parameter: float, optional
        :return: Yields a list of tuples containing the parameter and the
            cumulated distance along the BSplineCruve3D from the evaluation of
            start_parameter
        :rtype: Tuple[float, float]
        """
        resolution = max(10, min(resolution, int(self.length() / 1e-4)))
        delta_param = 1 / resolution * (end_parameter - start_parameter)
        distance = 0
        for i in range(resolution + 1):
            if i == 0:
                yield start_parameter, 0
            else:
                param1 = start_parameter + (i - 1) * delta_param
                param2 = start_parameter + i * delta_param
                point1 = self.evaluate_single(param1)
                point2 = self.evaluate_single(param2)
                distance += point1.point_distance(point2)
                yield param2, distance

    def normal(self, position: float = 0.0):
        """Returns the normal vector at a given parameter of the curve."""
        der = self.derivatives(position, 1)
        point1 = self.evaluate_single(0.0)
        points = [point1]
        count = 1
        u = 0.1
        while count < 3 and u <= 0.9:
            point = self.evaluate_single(u)
            if not volmdlr.core.point_in_list(point, points):
                points.append(point)
                count += 1
            u += 0.1
        if count < 3:
            raise NotImplementedError("BSplineCurve3D is a line segment")
        vec1 = points[1] - point1
        vec2 = points[2] - point1
        plane_normal = vec1.cross(vec2)
        normal = plane_normal.cross(der[1])
        return normal.unit_vector()

    def get_direction_vector(self, abscissa=0.0):
        """
        Calculates direction vector at given abscissa value (value between o and bspline length).

        """
        length = self.length()
        if abscissa >= length:
            abscissa2 = length
            abscissa = abscissa2 - 0.001 * length

        else:
            abscissa2 = min(abscissa + 0.001 * length, length)

        tangent = self.point_at_abscissa(abscissa2) - self.point_at_abscissa(
            abscissa)
        return tangent

    def direction_vector(self, abscissa=0.):
        """
        Gets direction vector at given abscissa value (value between o and bspline length).

        """
        if not self._direction_vector_memo:
            self._direction_vector_memo = {}
        if abscissa not in self._direction_vector_memo:
            self._direction_vector_memo[abscissa] = self.get_direction_vector(abscissa)
        return self._direction_vector_memo[abscissa]

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a BSplineCurve3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated
        :type object_dict: dict
        :return: The corresponding BSplineCurve3D.
        :rtype: :class:`volmdlr.edges.BSplineCurve3D`
        """
        name = arguments[0][1:-1]
        degree = int(arguments[1])
        points = [object_dict[int(i[1:])] for i in arguments[2]]
        lines = [LineSegment3D(pt1, pt2) for pt1, pt2 in zip(points[:-1], points[1:]) if not pt1.is_close(pt2)]
        if lines and not points[0].is_close(points[-1]):
            # quick fix. Real problem: Tolerance too low (1e-6 m = 0.001mm)
            dir_vector = lines[0].unit_direction_vector()
            if all(line.unit_direction_vector() == dir_vector for line in lines):
                return LineSegment3D(points[0], points[-1])

        knot_multiplicities = [int(i) for i in arguments[6][1:-1].split(",")]
        knots = [float(i) for i in arguments[7][1:-1].split(",")]
        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot] * knot_multiplicities[i])

        if 9 in range(len(arguments)):
            weight_data = [float(i) for i in arguments[9][1:-1].split(",")]
        else:
            weight_data = None

        return cls(degree, points, knot_multiplicities, knots, weight_data, name)

    def to_step(self, current_id, surface_id=None, curve2d=None):
        """Exports to STEP format."""
        points_ids = []
        content = ''
        point_id = current_id
        for point in self.control_points:
            point_content, point_id = point.to_step(point_id,
                                                    vertex=False)
            content += point_content
            points_ids.append(point_id)
            point_id += 1

        curve_id = point_id
        content += f"#{curve_id} = B_SPLINE_CURVE_WITH_KNOTS('{self.name}',{self.degree}," \
                   f"({volmdlr.core.step_ids_to_str(points_ids)})," \
                   f".UNSPECIFIED.,.F.,.F.,{tuple(self.knot_multiplicities)},{tuple(self.knots)}," \
                   f".UNSPECIFIED.);\n"

        if surface_id and curve2d:
            content += f"#{curve_id + 1} = SURFACE_CURVE('',#{curve_id},(#{curve_id + 2}),.PCURVE_S1.);\n"
            content += f"#{curve_id + 2} = PCURVE('',#{surface_id},#{curve_id + 3});\n"

            # 2D parametric curve
            curve2d_content, curve2d_id = curve2d.to_step(curve_id + 3)  # 5

            # content += f"#{curve_id + 3} = DEFINITIONAL_REPRESENTATION('',(#{curve2d_id - 1}),#{curve_id + 4});\n"
            # content += f"#{curve_id + 4} = ( GEOMETRIC_REPRESENTATION_CONTEXT(2)" \
            #            f"PARAMETRIC_REPRESENTATION_CONTEXT() REPRESENTATION_CONTEXT('2D SPACE','') );\n"

            content += curve2d_content
            current_id = curve2d_id
        else:
            current_id = curve_id + 1

        start_content, start_id = self.start.to_step(current_id, vertex=True)
        current_id = start_id + 1
        end_content, end_id = self.end.to_step(current_id + 1, vertex=True)
        content += start_content + end_content
        current_id = end_id + 1
        if surface_id:
            content += f"#{current_id} = EDGE_CURVE('{self.name}',#{start_id},#{end_id},#{curve_id},.T.);\n"
        else:
            content += f"#{current_id} = EDGE_CURVE('{self.name}',#{start_id},#{end_id},#{curve_id},.T.);\n"
        return content, current_id

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        BSplineCurve3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated BSplineCurve3D
        """
        new_control_points = [point.rotation(center, axis, angle) for point in
                              self.control_points]
        new_bsplinecurve3d = BSplineCurve3D(self.degree, new_control_points,
                                            self.knot_multiplicities,
                                            self.knots, self.weights, self.name)
        return new_bsplinecurve3d

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D, same_sense: bool = True):
        """
        Trims a bspline curve between two points.

        :param point1: point 1 used to trim.
        :param point2: point2 used to trim.
        :param same_sense: Used for periodical curves only. Indicates whether the curve direction agrees with (True)
            or is in the opposite direction (False) to the edge direction. By default, it's assumed True
        :return: New BSpline curve between these two points.
        """
        if self.periodic:
            return self.trim_with_interpolation(point1, point2, same_sense)
        bsplinecurve = self
        if not same_sense:
            bsplinecurve = self.reverse()
        parameter1 = bsplinecurve.point_to_parameter(point1)
        parameter2 = bsplinecurve.point_to_parameter(point2)

        if (point1.is_close(bsplinecurve.start) and point2.is_close(bsplinecurve.end)) \
                or (point1.is_close(bsplinecurve.end) and point2.is_close(bsplinecurve.start)):
            return bsplinecurve

        if point1.is_close(bsplinecurve.start) and not point2.is_close(bsplinecurve.end):
            return bsplinecurve.cut_after(parameter2)

        if point2.is_close(bsplinecurve.start) and not point1.is_close(bsplinecurve.end):
            bsplinecurve = bsplinecurve.cut_after(parameter1)
            return bsplinecurve

        if not point1.is_close(bsplinecurve.start) and point2.is_close(bsplinecurve.end):
            return bsplinecurve.cut_before(parameter1)

        if not point2.is_close(bsplinecurve.start) and point1.is_close(bsplinecurve.end):
            bsplinecurve = bsplinecurve.cut_before(parameter2)
            return bsplinecurve

        if parameter1 is None or parameter2 is None:
            raise ValueError('Point not on BSplineCurve for trim method')

        if parameter1 > parameter2:
            parameter1, parameter2 = parameter2, parameter1
            point1, point2 = point2, point1

        bsplinecurve = bsplinecurve.cut_before(parameter1)
        new_param2 = bsplinecurve.point_to_parameter(point2)
        trimmed_bspline_cruve = bsplinecurve.cut_after(new_param2)
        return trimmed_bspline_cruve

    def trim_with_interpolation(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D, same_sense: bool = True):
        """
        Creates a new BSplineCurve3D between point1 and point2 using interpolation method.
        """
        bspline_curve = self
        if not same_sense:
            bspline_curve = self.reverse()
        n = len(bspline_curve.control_points)
        local_discretization = bspline_curve.local_discretization(point1, point2, n, tol=1e-8)
        if len(local_discretization) <= bspline_curve.degree:
            return bspline_curve
        return bspline_curve.__class__.from_points_interpolation(local_discretization, bspline_curve.degree)

    def trim_between_evaluations(self, parameter1: float, parameter2: float):
        """
        Trims the Bspline between two abscissa evaluation parameters.

        :param parameter1: evaluation parameter 1, bigger than 0 and smaller than its length.
        :param parameter2: evaluation parameter 2, bigger than 0 and smaller than its length.
        """
        warnings.warn('Use BSplineCurve3D.trim instead of trim_between_evaluation')
        parameter1, parameter2 = min([parameter1, parameter2]), \
            max([parameter1, parameter2])

        if math.isclose(parameter1, 0, abs_tol=1e-7) \
                and math.isclose(parameter2, 1, abs_tol=1e-7):
            return self
        if math.isclose(parameter1, 0, abs_tol=1e-7):
            return self.cut_after(parameter2)
        if math.isclose(parameter2, 1, abs_tol=1e-7):
            return self.cut_before(parameter1)

        # Cut before
        bspline_curve = self.insert_knot(parameter1, num=self.degree)
        if bspline_curve.weights is not None:
            raise NotImplementedError

        # Cut after
        bspline_curve = bspline_curve.insert_knot(parameter2, num=self.degree)
        if bspline_curve.weights is not None:
            raise NotImplementedError

        new_ctrlpts = bspline_curve.control_points[bspline_curve.degree:
                                                   -bspline_curve.degree]
        new_multiplicities = bspline_curve.knot_multiplicities[1:-1]
        # new_multiplicities = bspline_curve.knot_multiplicities[2:-5]
        new_multiplicities[-1] += 1
        new_multiplicities[0] += 1
        new_knots = bspline_curve.knots[1:-1]
        # new_knots = bspline_curve.knots[2:-5]
        new_knots = nurbs_helpers.standardize_knot_vector(new_knots)

        return BSplineCurve3D(degree=bspline_curve.degree,
                              control_points=new_ctrlpts,
                              knot_multiplicities=new_multiplicities,
                              knots=new_knots,
                              weights=None,
                              name=bspline_curve.name)

    def cut_before(self, parameter: float):
        """
        Returns the right side of the split curve at a given parameter.

        :param parameter: parameter value that specifies where to split the curve.
        :type parameter: float
        """
        point3d = self.evaluate_single(parameter)
        if self.start.is_close(point3d):
            return self.copy()
        if self.end.is_close(point3d):
            return self.reverse()

        curves = volmdlr.nurbs.operations.split_curve(self, round(parameter, 7))
        return curves[1]

    def cut_after(self, parameter: float):
        """
        Returns the left side of the split curve at a given parameter.

        :param parameter: parameter value that specifies where to split the curve.
        :type parameter: float
        """
        # Is a value of parameter below 5e-6 a real need for precision ?
        point3d = self.evaluate_single(parameter)
        if self.start.is_close(point3d):
            return self.reverse()
        if self.end.is_close(point3d):
            return self.copy()
        curves = volmdlr.nurbs.operations.split_curve(self, round(parameter, 7))
        return curves[0]

    def insert_knot(self, knot: float, num: int = 1):
        """
        Returns a new BSplineCurve3D.

        """
        return volmdlr.nurbs.operations.insert_knot_curve(self, [knot], num=[num])

    # Copy paste du LineSegment3D
    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        points = self.points
        x = [point.x for point in points]
        y = [point.y for point in points]
        z = [point.z for point in points]
        ax.plot(x, y, z, color=edge_style.color, alpha=edge_style.alpha)
        if edge_style.edge_ends:
            ax.plot(x, y, z, 'o', color=edge_style.color, alpha=edge_style.alpha)
        return ax

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a BSplineCurve3D into an BSplineCurve2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: BSplineCurve2D.
        """
        control_points2d = [point.to_2d(plane_origin, x, y) for point in
                            self.control_points]
        return BSplineCurve2D(self.degree, control_points2d,
                              self.knot_multiplicities, self.knots,
                              self.weights, self.name)

    def curvature(self, u: float, point_in_curve: bool = False):
        """
        Returns the curvature of a curve and the point where it is located.
        """
        ders = self.derivatives(u, 3)  # 3 first derivative
        c1, c2 = ders[1], ders[2]
        denom = c1.cross(c2)
        if c1.is_close(volmdlr.O3D) or c2.is_close(volmdlr.O3D) or denom.norm() == 0.0:
            if point_in_curve:
                return 0., volmdlr.Point3D(*ders[0])
            return 0.
        r_c = ((c1.norm()) ** 3) / denom.norm()
        point = volmdlr.Point3D(*ders[0])
        if point_in_curve:
            return 1 / r_c, point
        return 1 / r_c

    def global_maximum_curvature(self, nb_eval: int = 21, point_in_curve: bool = False):
        """
        Returns the global maximum curvature of a curve and the point where it is located.
        """
        check = [i / (nb_eval - 1) for i in range(nb_eval)]
        curvatures = []
        for u in check:
            curvatures.append(self.curvature(u, point_in_curve))
        return curvatures

    def maximum_curvature(self, point_in_curve: bool = False):
        """
        Returns the maximum curvature of a curve and the point where it is located.
        """
        if point_in_curve:
            maximum_curvarture, point = max(self.global_maximum_curvature(nb_eval=21, point_in_curve=point_in_curve))
            return maximum_curvarture, point
        maximum_curvarture = max(self.global_maximum_curvature(nb_eval=21, point_in_curve=point_in_curve))
        return maximum_curvarture

    def minimum_radius(self, point_in_curve=False):
        """
        Returns the minimum curvature radius of a curve and the point where it is located.
        """
        if point_in_curve:
            maximum_curvarture, point = self.maximum_curvature(point_in_curve)
            return 1 / maximum_curvarture, point
        maximum_curvarture = self.maximum_curvature(point_in_curve)
        return 1 / maximum_curvarture

    # def global_minimum_curvature(self, nb_eval: int = 21):
    #     check = [i / (nb_eval - 1) for i in range(nb_eval)]
    #     radius = []
    #     for u in check:
    #         radius.append(self.minimum_curvature(u))
    #     return radius

    def triangulation(self):
        """Triangulation method for a BSplineCurve3D."""
        return None

    def linesegment_intersections(self, linesegment3d: LineSegment3D, abs_tol: float = 1e-6):
        """
        Calculates intersections between a BSplineCurve3D and a LineSegment3D.

        :param linesegment3d: linesegment to verify intersections.
        :param abs_tol: tolerance.
        :return: list with the intersections points.
        """
        if not self.bounding_box.bbox_intersection(linesegment3d.bounding_box, abs_tol):
            return []
        intersection_section_pairs = self._get_intersection_sections(linesegment3d)
        intersections = []
        for bspline, edge2_ in intersection_section_pairs:
            intersections_points = bspline.get_linesegment_intersections(edge2_)
            for inter in intersections_points:
                if not volmdlr.core.point_in_list(inter, intersections, abs_tol):
                    intersections.append(inter)
        return intersections

    def arc_intersections(self, arc, abs_tol=1e-6):
        """
        Calculates intersections between a BSpline Curve 3D and an arc 3D.

        :param arc: arc to verify intersections.
        :param abs_tol: tolerance.
        :return: list with the intersections points.
        """
        if self.bounding_box.distance_to_bbox(arc.bounding_box) > abs_tol:
            return []
        return self._generic_edge_intersections(arc, abs_tol)

    def is_shared_section_possible(self, other_bspline2, tol):
        """
        Verifies if it there is any possibility of the two bsplines share a section.

        :param other_bspline2: other bspline.
        :param tol: tolerance used.
        :return: True or False.
        """
        if self.bounding_box.distance_to_bbox(other_bspline2.bounding_box) > tol:
            return False
        return True

    def sweep(self, *args):
        """
        Bspline 3D is used as path for sweeping given section through it.

        :return:
        """
        new_faces = []
        tangents = []
        section_contour2d, _ = args
        points = self.points
        for k, _ in enumerate(points):
            position = k / (len(points) - 1)
            tangents.append(self.unit_direction_vector(position * self.length()))

        contours = []
        for point, tan in zip(points, tangents):
            normal = tan.deterministic_unit_normal_vector()
            v_vector = tan.cross(normal)
            section_contour3d = section_contour2d.to_3d(point, normal, v_vector)
            contours.append(section_contour3d)

        polys = [volmdlr.wires.ClosedPolygon3D(c.discretization_points(number_points=36)) for c in contours]

        size_v, size_u = len(polys[0].points), len(polys)
        degree_u, degree_v = 3, 3

        points_3d = []
        for poly in polys:
            points_3d.extend(poly.points)

        bspline_surface3d = volmdlr.surfaces.BSplineSurface3D.from_points_interpolation(points_3d, size_u,
                                                                                       size_v,degree_u, degree_v)

        outer_contour = volmdlr.wires.Contour2D([volmdlr.edges.LineSegment2D(volmdlr.O2D, volmdlr.X2D.to_point()),
                                                 volmdlr.edges.LineSegment2D(
                                                     volmdlr.X2D.to_point(), (volmdlr.X2D + volmdlr.Y2D).to_point()),
                                                 volmdlr.edges.LineSegment2D(
                                                     (volmdlr.X2D + volmdlr.Y2D).to_point(), volmdlr.Y2D.to_point()),
                                                 volmdlr.edges.LineSegment2D(volmdlr.Y2D.to_point(), volmdlr.O2D)])
        surf2d = volmdlr.surfaces.Surface2D(outer_contour, [])

        bsface3d = volmdlr.faces.BSplineFace3D(bspline_surface3d, surf2d)
        new_faces.append(bsface3d)
        return new_faces

    def revolution(self, axis_point, axis, angle):
        """
        Returns the face generated by the revolution of the BSpline Curve 3D.
        """
        surface = volmdlr.surfaces.RevolutionSurface3D(self, axis_point, axis)
        face = volmdlr.faces.RevolutionFace3D.from_surface_rectangular_cut(surface, 0, angle, 0, self.length())
        return face

    def split_between_two_points(self, point1, point2):
        """
        Split edge between two points.

        :param point1: point 1.
        :param point2: point 2.
        :return: edge split.
        """
        return self.trim(point1, point2)

class BezierCurve3D(BSplineCurve3D):
    """
    A class for 3-dimensional Bézier curves.

    :param degree: The degree of the Bézier curve
    :type degree: int
    :param control_points: A list of 3-dimensional points
    :type control_points: List[:class:`volmdlr.Point3D`]
    :param name: The name of the B-spline curve. Default value is ''
    :type name: str, optional
    """

    def __init__(self, degree: int, control_points: List[volmdlr.Point3D],
                 name: str = ''):
        knotvector = nurbs_helpers.generate_knot_vector(degree,
                                                        len(control_points))
        knot_multiplicity = [1] * len(knotvector)

        BSplineCurve3D.__init__(self, degree, control_points,
                                knot_multiplicity, knotvector,
                                None, name)


class Arc3D(ArcMixin, Edge):
    """
    An arc is defined by a starting point, an end point and an interior point.

    """

    def __init__(self, circle, start, end, name=''):
        ArcMixin.__init__(self, circle, start=start, end=end, name=name)
        Edge.__init__(self, start=start, end=end, name=name)
        self._angle = None
        self.frame = self.circle.frame
        self.radius = self.circle.radius
        self.angle_start, self.angle_end = self.get_start_end_angles()
        self._bbox = None

    def __hash__(self):
        return hash(('arc3d', self.circle, self.start, self.end, self.is_trigo))

    def __eq__(self, other_arc):
        if self.__class__.__name__ != other_arc.__class__.__name__:
            return False
        return (self.circle == other_arc.circle and self.start == other_arc.start
                and self.end == other_arc.end and self.is_trigo == other_arc.is_trigo)

    def to_dict(self, use_pointers: bool = False, memo=None, path: str = '#', id_method=True, id_memo=None):
        """Saves the object parameters into a dictionary."""
        dict_ = self.base_dict()
        dict_['circle'] = self.circle.to_dict(use_pointers=use_pointers, memo=memo,
                                              id_method=id_method, id_memo=id_memo, path=path + '/circle')
        dict_['start'] = self.start.to_dict(use_pointers=use_pointers, memo=memo,
                                            id_method=id_method, id_memo=id_memo, path=path + '/start')
        dict_['end'] = self.end.to_dict(use_pointers=use_pointers, memo=memo,
                                        id_method=id_method, id_memo=id_memo, path=path + '/end')
        return dict_

    def _arc_point_angle(self, point):
        """Helper function to calculate the angle of point on a trigonometric arc."""
        local_start_point = self.circle.frame.global_to_local_coordinates(point)
        u1, u2 = local_start_point.x / self.radius, local_start_point.y / self.radius
        point_angle = volmdlr.geometry.sin_cos_angle(u1, u2)
        return point_angle

    def get_arc_point_angle(self, point):
        """Returns the angle of point on a trigonometric arc."""
        point_theta = self._arc_point_angle(point)
        if self.angle_start > point_theta:
            point_theta += volmdlr.TWO_PI
        return point_theta

    def get_start_end_angles(self):
        """Returns the start and end angle of the arc."""
        start_angle = self._arc_point_angle(self.start)
        end_angle = self._arc_point_angle(self.end)
        if start_angle >= end_angle:
            end_angle += volmdlr.TWO_PI
        return start_angle, end_angle

    @property
    def bounding_box(self):
        """Bounding box for Arc 3D."""
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def get_bounding_box(self):
        """
        Calculates the bounding box of the Arc3D.

        :return: Bounding Box object.
        """
        # TODO: implement exact calculation

        points = self.discretization_points(angle_resolution=5)
        xmin = min(point.x for point in points)
        xmax = max(point.x for point in points)
        ymin = min(point.y for point in points)
        ymax = max(point.y for point in points)
        zmin = min(point.z for point in points)
        zmax = max(point.z for point in points)
        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    @classmethod
    def from_angle(cls, start: volmdlr.Point3D, angle: float,
                   axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D, name: str = ''):
        """Gives the arc3D from a start, an angle and an axis."""
        start_gen = start
        end_gen = start_gen.rotation(axis_point, axis, angle)
        line = volmdlr_curves.Line3D(axis_point, axis_point + axis)
        center, _ = line.point_projection(start)
        radius = center.point_distance(start)
        u = start - center
        v = axis.cross(u)
        circle = volmdlr.curves.Circle3D(volmdlr.Frame3D(center, u, v, axis), radius)
        if angle == volmdlr.TWO_PI:
            return circle
        return cls(circle, start_gen, end_gen, name=name)

    @classmethod
    def from_3_points(cls, point1, point2, point3, name: str = ''):
        """
        Creates an Arc 3d using three points.

        :param point1: start point.
        :param point2: interior point.
        :param point3: end point.
        :param name: object's name.
        :return: Arc 3D.
        """
        circle = volmdlr_curves.Circle3D.from_3_points(point1, point2, point3)
        arc = cls(circle, point1, point3, name=name)
        return arc

    @property
    def angle(self):
        """
        Arc angle property.

        :return: arc angle.
        """
        if not self._angle:
            self._angle = self.angle_end - self.angle_start
        return self._angle

    @property
    def points(self):
        return [self.start, self.end]

    def get_reverse(self):
        """
        Defines a new Arc3D, identical to self, but in the opposite direction.

        """
        circle3d = self.circle.reverse()
        return self.__class__(circle3d, self.end, self.start, self.name + '_reverse')

    def abscissa(self, point: volmdlr.Point3D, tol: float = 1e-6):
        """
        Calculates the abscissa given a point in the Arc3D.

        :param point: point to calculate the abscissa.
        :param tol: (Optional) Confusion distance to consider points equal. Default 1e-6.
        :return: corresponding abscissa.
        """
        if point.point_distance(self.start) <= tol:
            return 0
        if point.point_distance(self.end) <= tol:
            return self.length()
        point_theta = self.get_arc_point_angle(point)
        if not self.angle_start <= point_theta <= self.angle_end:
            raise ValueError(f"{point} not in Arc3D.")
        return self.radius * abs(point_theta - self.angle_start)

    def point_at_abscissa(self, abscissa):
        """
        Calculates a point in the Arc3D at a given abscissa.

        :param abscissa: abscissa where in the curve the point should be calculated.
        :return: Corresponding point.
        """
        if abscissa > self.length() + 1e-6:
            raise ValueError(f"{abscissa} abscissa is not on the curve. max length of arc is {self.length()}.")
        return self.start.rotation(self.circle.center, self.circle.normal, abscissa / self.radius)

    def direction_vector(self, abscissa):
        """
        Calculates a direction vector at a given abscissa of the Arc3D.

        :param abscissa: abscissa where in the curve the direction vector should be calculated.
        :return: Corresponding direction vector.
        """
        normal_vector = self.normal_vector(abscissa)
        tangent = normal_vector.cross(self.circle.normal)
        return tangent

    def rotation(self, center: volmdlr.Point3D,
                 axis: volmdlr.Vector3D, angle: float):
        """
        Arc3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Arc3D
        """
        circle = self.circle.rotation(center, axis, angle)
        new_start = self.start.rotation(center, axis, angle)
        new_end = self.end.rotation(center, axis, angle)
        return Arc3D(circle, new_start, new_end, name=self.name)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Arc3D translation.

        :param offset: translation vector.
        :return: A new translated Arc3D.
        """
        new_circle = self.circle.translation(offset)
        new_start = self.start.translation(offset)
        new_end = self.end.translation(offset)
        return Arc3D(new_circle, new_start, new_end, name=self.name)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and return a new Arc3D.

        side = 'old' or 'new'
        """
        new_circle = self.circle.frame_mapping(frame, side)
        new_start = self.start.frame_mapping(frame, side)
        new_end = self.end.frame_mapping(frame, side)
        return Arc3D(new_circle, new_start, new_end, name=self.name)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plot method for Arc 3D using Matplotlib."""
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        ax = vm_common_operations.plot_from_discretization_points(
            ax, edge_style=edge_style, element=self, number_points=25)
        if edge_style.edge_ends:
            self.start.plot(ax=ax, color='r')
            self.end.plot(ax=ax, color='b')

        if edge_style.edge_direction:
            x, y, z = self.point_at_abscissa(0.5 * self.length())
            u, v, w = 0.05 * self.unit_direction_vector(0.5 * self.length())
            ax.quiver(x, y, z, u, v, w, length=self.length() / 100,
                      arrow_length_ratio=5, normalize=True,
                      pivot='tip', color=edge_style.color)
        return ax

    def plot2d(self, center: volmdlr.Point3D = volmdlr.O3D,
               x3d: volmdlr.Vector3D = volmdlr.X3D, y3d: volmdlr.Vector3D = volmdlr.Y3D,
               ax=None, color='k'):
        """Plot data."""

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        # TODO: Enhance this plot
        length = self.length()
        x = []
        y = []
        for i in range(30):
            point = self.point_at_abscissa(i / 29. * length)
            xi, yi = point.plane_projection2d(center, x3d, y3d)
            x.append(xi)
            y.append(yi)
        ax.plot(x, y, color=color)

        return ax

    def copy(self, *args, **kwargs):
        return Arc3D(self.circle.copy(), self.start.copy(), self.end.copy())

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a Arc3D into an Arc2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Arc2D.
        """
        circle2d = self.circle.to_2d(plane_origin, x, y)
        point_start = self.start.to_2d(plane_origin, x, y)
        point_interior = self.middle_point().to_2d(plane_origin, x, y)
        point_end = self.end.to_2d(plane_origin, x, y)
        arc = Arc2D(circle2d, point_start, point_end, self.is_trigo, name=self.name)
        if not arc.point_belongs(point_interior):
            arc = Arc2D(circle2d, point_start, point_end, False, name=self.name)
        return arc

    def minimum_distance_points_arc(self, other_arc):
        """Calculates the minimum distance points between two arcs."""
        u1 = self.start - self.circle.center
        u1 = u1.unit_vector()
        u2 = self.circle.normal.cross(u1)

        w = other_arc.circle.center - self.circle.center

        u3 = other_arc.start - other_arc.circle.center
        u3 = u3.unit_vector()
        u4 = other_arc.circle.normal.cross(u3)

        radius1, radius2 = self.radius, other_arc.radius

        u1_u1, u1_u2, u1_u3, u1_u4 = u1.dot(u1), u1.dot(u2), u1.dot(u3), u1.dot(u4)
        u2_u2, u2_u3, u2_u4 = u2.dot(u2), u2.dot(u3), u2.dot(u4)
        u3_u3, u3_u4 = u3.dot(u3), u3.dot(u4)
        u4_u4 = u4.dot(u4)
        w_u1, w_u2, w_u3, w_u4, w_w = w.dot(u1), w.dot(u2), w.dot(u3), w.dot(u4), w.dot(w)

        def distance_squared(x):
            return (u1_u1 * ((math.cos(x[0])) ** 2) * radius1 ** 2 + u2_u2 * (
                    (math.sin(x[0])) ** 2) * radius1 ** 2
                    + w_w + u3_u3 * ((math.cos(x[1])) ** 2) * radius2 ** 2 + u4_u4 * (
                            (math.sin(x[1])) ** 2) * radius2 ** 2
                    + u1_u2 * math.sin(2 * x[0]) * radius1 ** 2 - 2 * radius1 * math.cos(
                        x[0]) * w_u1
                    - 2 * radius1 * radius2 * math.cos(x[0]) * math.cos(x[1]) * u1_u3
                    - 2 * radius1 * radius2 * math.cos(x[0]) * math.sin(
                        x[1]) * u1_u4 - 2 * radius1 * math.sin(x[0]) * w_u2
                    - 2 * radius1 * radius2 * math.sin(x[0]) * math.cos(x[1]) * u2_u3
                    - 2 * radius1 * radius2 * math.sin(x[0]) * math.sin(
                        x[1]) * u2_u4 + 2 * radius2 * math.cos(x[1]) * w_u3
                    + 2 * radius2 * math.sin(x[1]) * w_u4 + u3_u4 * math.sin(
                        2 * x[1]) * radius2 ** 2)

        x01 = npy.array([self.angle / 2, other_arc.angle / 2])

        res1 = least_squares(distance_squared, x01, bounds=[(0, 0), (self.angle, other_arc.angle)])

        point1 = self.point_at_abscissa(res1.x[0] * radius1)
        point2 = other_arc.point_at_abscissa(res1.x[1] * radius2)

        return point1, point2

    def distance_linesegment(self, linesegment3d, return_points=False):
        """
        Gets the minimum distance between an Arc 3D and Line Segment 3D.

        :param linesegment3d: other line segment 3d.
        :param return_points: boolean to decide weather to return the corresponding minimal distance points or not.
        :return: minimum distance / minimal distance with corresponding points.
        """
        point1, point2 = vm_common_operations.minimum_distance_points_circle3d_linesegment3d(self, linesegment3d)
        if return_points:
            return point1.point_distance(point2), point1, point2
        return point1.point_distance(point2)

    def distance_arc(self, arc3d, return_points=False):
        """
        Gets the minimum distance between two Arcs 3D.

        :param arc3d: other arc 3d.
        :param return_points: boolean to decide weather to return the corresponding minimal distance points or not.
        :return: minimum distance / minimal distance with corresponding points.
        """
        p1, p2 = self.minimum_distance_points_arc(arc3d)
        if return_points:
            return p1.point_distance(p2), p1, p2
        return p1.point_distance(p2)

    def extrusion(self, extrusion_vector):
        """Extrudes an arc 3d in the given extrusion vector direction."""
        if self.circle.normal.is_colinear_to(extrusion_vector):
            u = self.start - self.circle.center
            u = u.unit_vector()
            w = extrusion_vector.copy()
            w = w.unit_vector()
            v = w.cross(u)
            arc2d = self.to_2d(self.circle.center, u, v)
            angle1, angle2 = arc2d.angle1, arc2d.angle2
            if angle2 < angle1:
                angle2 += volmdlr.TWO_PI
            cylinder = volmdlr.surfaces.CylindricalSurface3D(
                volmdlr.Frame3D(self.circle.center, u, v, w),
                self.radius
            )
            return [volmdlr.faces.CylindricalFace3D.from_surface_rectangular_cut(
                cylinder, angle1, angle2, 0., extrusion_vector.norm())]
        raise NotImplementedError(f'Elliptic faces not handled: dot={self.circle.normal.dot(extrusion_vector)}')

    def revolution(self, axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                   angle: float):
        """
        Revolution of Arc 3D around an axis.

        :param axis_point: revolution axis point.
        :param axis: revolution axis.
        :param angle: revolution angle.
        """
        line3d = volmdlr_curves.Line3D(axis_point, axis_point + axis)
        tore_center, _ = line3d.point_projection(self.circle.center)

        # Sphere
        if math.isclose(tore_center.point_distance(self.circle.center), 0.,
                        abs_tol=1e-6):

            start_p, _ = line3d.point_projection(self.start)
            u = self.start - start_p

            if math.isclose(u.norm(), 0, abs_tol=1e-6):
                end_p, _ = line3d.point_projection(self.end)
                u = self.end - end_p
                if math.isclose(u.norm(), 0, abs_tol=1e-6):
                    interior_p, _ = line3d.point_projection(self.middle_point())
                    u = self.middle_point - interior_p

            u = u.unit_vector()
            v = axis.cross(u)
            arc2d = self.to_2d(self.circle.center, u, axis)

            surface = volmdlr.surfaces.SphericalSurface3D(
                volmdlr.Frame3D(self.circle.center, u, v, axis), self.radius)

            return [volmdlr.faces.SphericalFace3D.from_surface_rectangular_cut(surface, 0, angle,
                                                                               arc2d.angle1, arc2d.angle2)]

        # Toroidal
        u = self.circle.center - tore_center
        u = u.unit_vector()
        v = axis.cross(u)
        if not math.isclose(self.circle.normal.dot(u), 0., abs_tol=1e-6):
            raise NotImplementedError(
                'Outside of plane revolution not supported')

        radius = tore_center.point_distance(self.circle.center)
        # from volmdlr import surfaces, faces
        surface = volmdlr.surfaces.ToroidalSurface3D(
            volmdlr.Frame3D(tore_center, u, v, axis), radius,
            self.radius)
        arc2d = self.to_2d(tore_center, u, axis)
        return [volmdlr.faces.ToroidalFace3D.from_surface_rectangular_cut(
            surface, 0, angle, arc2d.angle1, arc2d.angle2)]

    def to_step(self, current_id, *args, **kwargs):
        """
        Converts the object to a STEP representation.

        :param current_id: The ID of the last written primitive.
        :type current_id: int
        :return: The STEP representation of the object and the last ID.
        :rtype: tuple[str, list[int]]
        """
        content, frame_id = self.circle.frame.to_step(current_id)
        curve_id = frame_id + 1
        content += f"#{curve_id} = CIRCLE('{self.name}', #{frame_id}, {self.radius * 1000});\n"

        current_id = curve_id
        start_content, start_id = self.start.to_step(current_id, vertex=True)
        end_content, end_id = self.end.to_step(start_id + 1, vertex=True)
        content += start_content + end_content
        current_id = end_id + 1
        content += f"#{current_id} = EDGE_CURVE('{self.name}',#{start_id},#{end_id},#{curve_id},.T.);\n"
        return content, current_id

    def point_belongs(self, point, abs_tol: float = 1e-6):
        """
        Check if a point 3d belongs to the arc_3d or not.

        :param point: point to be verified is on arc.
        :param abs_tol: tolerance allowed.
        :return: True if point is on Arc, False otherwise.
        """
        # point_local_coordinates = self.circle.frame.global_to_local_coordinates(point)
        if not math.isclose(point.point_distance(self.circle.center), self.radius, abs_tol=abs_tol):
            return False
        vector = point - self.circle.center
        if not math.isclose(vector.dot(self.circle.frame.w), 0.0, abs_tol=abs_tol):
            return False
        point_theta = self.get_arc_point_angle(point)
        if not self.angle_start <= point_theta <= self.angle_end:
            return False
        return True

    def triangulation(self):
        """
        Triangulation for an Arc3D.

        """
        return None

    def line_intersections(self, line3d: volmdlr_curves.Line3D, tol: float = 1e-6):
        """
        Calculates intersections between an Arc3D and a Line3D.

        :param line3d: line to verify intersections.
        :param tol: maximum tolerance.
        :return: list with intersections points between line and Arc3D.
        """
        if line3d.point_belongs(self.start):
            return [self.start]
        if line3d.point_belongs(self.end):
            return [self.end]
        circle3d_lineseg_inters = vm_utils_intersections.circle_3d_line_intersections(self.circle, line3d)
        linesegment_intersections = []
        for intersection in circle3d_lineseg_inters:
            if self.point_belongs(intersection, tol):
                linesegment_intersections.append(intersection)
        return linesegment_intersections

    def linesegment_intersections(self, linesegment3d: LineSegment3D, abs_tol: float = 1e-6):
        """
        Calculates intersections between an Arc3D and a LineSegment3D.

        :param linesegment3d: linesegment to verify intersections.
        :param abs_tol: tolerance to be considered while validating an intersection.
        :return: list with intersections points between linesegment and Arc3D.
        """
        linesegment_intersections = []
        intersections = self.line_intersections(linesegment3d.line)
        for intersection in intersections:
            if linesegment3d.point_belongs(intersection, abs_tol=abs_tol):
                linesegment_intersections.append(intersection)
        return linesegment_intersections

    def arc_intersections(self, other_arc, abs_tol: 1e-6):
        """
        Calculates intersections between two Arc3D.

        :param other_arc: Arc 3D to verify intersections.
        :param abs_tol: tolerance.
        :return: list with intersections points between the two Arc3D.
        """
        circle_intersections = self.circle.circle_intersections(other_arc.circle)
        intersections = []
        for intersection in circle_intersections:
            if self.point_belongs(intersection, abs_tol) and other_arc.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    def arcellipse_intersections(self, arcellipse3d, abs_tol: float = 1e-6):
        """
        Calculates intersections between two Arc3D.

        :param arcellipse3d: ArcEllipse 3D to verify intersections.
        :param abs_tol: Tolerance.
        :return: List with intersections points between ArcEllipse3D and Arc3D.
        """
        ellipse_intersections = self.circle.ellipse_intersections(arcellipse3d.ellipse, abs_tol)
        intersections = []
        for intersection in ellipse_intersections:
            if self.point_belongs(intersection, abs_tol) and arcellipse3d.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    def complementary(self):
        """Creates the corresponding complementary arc."""
        return Arc3D(self.circle, self.end, self.start)

    def sweep(self, *args):
        """
        Arc 3D is used as path for sweeping given section through it.

        :return:
        """
        new_faces = []
        section_contour2d, section_contour3d = args
        if section_contour3d is None:
            start_tangent = self.unit_direction_vector(0.)
            normal = self.unit_normal_vector(0.)
            tangent_normal_orthonormal = start_tangent.cross(normal)
            section_contour3d = section_contour2d.to_3d(self.start, normal, tangent_normal_orthonormal)
        for contour_primitive in section_contour3d.primitives:
            new_faces.extend(contour_primitive.revolution(
                self.circle.center, self.circle.normal, self.angle))
        return new_faces


class FullArc3D(FullArcMixin, Arc3D):
    """
    An edge that starts at start_end, ends at the same point after having described a circle.

    """

    def __init__(self, circle: volmdlr.curves.Circle3D, start_end: volmdlr.Point3D,
                 name: str = ''):
        self._utd_frame = None
        self._bbox = None
        FullArcMixin.__init__(self, circle=circle, start_end=start_end, name=name)
        Arc3D.__init__(self, circle=circle, start=start_end, end=start_end, name=name)

    def __hash__(self):
        return hash(('Fullarc3D', self.circle, self.start_end))

    def __eq__(self, other_arc):
        return (self.circle == other_arc.circle) \
            and (self.start == other_arc.start)

    def copy(self, *args, **kwargs):
        """Returns a new instance with the same parameters."""
        return FullArc3D(self.circle.copy(), self.end.copy())

    def to_dict(self, use_pointers: bool = False, memo=None, path: str = '#'):
        """Object serialization."""
        dict_ = self.base_dict()
        dict_['circle'] = self.circle.to_dict(use_pointers=use_pointers, memo=memo, path=path + '/circle')
        dict_['angle'] = self.angle
        dict_['is_trigo'] = self.is_trigo
        dict_['start_end'] = self.start.to_dict(use_pointers=use_pointers, memo=memo, path=path + '/start_end')
        dict_['name'] = self.name
        return dict_

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a FullArc3D into an FullArc2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: FullArc2D.
        """
        circle = self.circle.to_2d(plane_origin, x, y)
        start_end = self.start.to_2d(plane_origin, x, y)
        return FullArc2D(circle, start_end)

    def to_step(self, current_id, surface_id=None):
        """Exports to STEP format."""
        content, frame_id = self.circle.frame.to_step(current_id)
        # Not calling Circle3D.to_step because of circular imports
        u = self.start - self.circle.center
        u = u.unit_vector()
        curve_id = frame_id + 1
        # Not calling Circle3D.to_step because of circular imports
        content += f"#{curve_id} = CIRCLE('{self.name}',#{frame_id},{self.radius * 1000});\n"

        point1 = (self.circle.center + u * self.radius).to_point()

        p1_content, p1_id = point1.to_step(curve_id, vertex=True)
        content += p1_content

        edge_curve = p1_id + 1
        content += f"#{edge_curve} = EDGE_CURVE('{self.name}',#{p1_id},#{p1_id},#{curve_id},.T.);\n"

        return content, edge_curve

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle(), show_frame=False):
        """
        Plot fullarc3d using matplotlib.
        """
        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')
        if show_frame:
            self.circle.frame.plot(ax, ratio=self.radius)
        ax = vm_common_operations.plot_from_discretization_points(
            ax, edge_style=edge_style, element=self, number_points=25, close_plot=True)
        if edge_style.edge_ends:
            self.start.plot(ax=ax)
            self.end.plot(ax=ax)
        if edge_style.edge_direction:
            half_length = 0.5 * self.length()
            x, y, z = self.point_at_abscissa(half_length)
            tangent = self.unit_direction_vector(half_length)
            arrow_length = 0.15 * half_length
            ax.quiver(x, y, z, *arrow_length * tangent, pivot='tip')

        return ax

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Rotates the FullArc3D object around a specified axis by a given angle.

        :param center: The center point of rotation.
        :type center: (volmdlr.Point3D)
        :param axis: The axis of rotation.
        :type axis: (volmdlr.Vector3D)
        :param angle: The angle of rotation in radians.
        :type angle: (float)

        :return: A new FullArc3D object that is the result of the rotation.
        :rtype: FullArc3D:
        """
        new_start_end = self.start.rotation(center, axis, angle)
        new_circle = self.circle.rotation(center, axis, angle)
        return FullArc3D(new_circle, new_start_end, name=self.name)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Translates the FullArc3D object by a specified offset.

        :param offset: The translation offset vector.
        :type offset: (volmdlr.Vector3D).
        :return: A new FullArc3D object that is the result of the translation.
        :rtype: FullArc3D.
        """
        new_start_end = self.start.translation(offset)
        new_circle = self.circle.translation(offset)
        return FullArc3D(new_circle, new_start_end, name=self.name)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and return a new FullArc3D.

        side = 'old' or 'new'
        """
        new_circle = self.circle.frame_mapping(frame, side)
        new_start_end = self.start_end.frame_mapping(frame, side)
        return FullArc3D(new_circle, new_start_end, name=self.name)

    def linesegment_intersections(self, linesegment3d: LineSegment3D, abs_tol=1e-6):
        """
        Calculates the intersections between a full arc 3d and a line segment 3d.

        :param linesegment3d: linesegment 3d to verify intersections.
        :param abs_tol: tolerance.
        :return: list of points 3d, if there are any intersections, an empty list if otherwise.
        """
        distance_center_lineseg = linesegment3d.point_distance(self.circle.frame.origin)
        if distance_center_lineseg > self.radius:
            return []
        return self.circle.linesegment_intersections(linesegment3d)

    def get_reverse(self):
        """
        Defines a new FullArc3D, identical to self, but in the opposite direction.

        """
        circle = self.circle.reverse()
        return self.__class__(circle, self.start_end)

    def point_belongs(self, point: volmdlr.Point3D, abs_tol: float = 1e-6):
        """
        Returns if given point belongs to the FullArc3D.
        """
        distance = point.point_distance(self.circle.center)
        vec = volmdlr.Vector3D(*point - self.circle.center)
        dot = self.circle.normal.dot(vec)
        return math.isclose(distance, self.radius, abs_tol=abs_tol) \
            and math.isclose(dot, 0, abs_tol=abs_tol)

    @classmethod
    def from_3_points(cls, point1, point2, point3, name: str = ''):
        """Creates a full arc 3D from 3 given points."""
        fullarc = cls(volmdlr_curves.Circle3D.from_3_points(point1, point2, point3), point1, name=name)
        return fullarc

    def split(self, split_point, tol: float = 1e-6):
        """
        Splits the circle into two arcs at a given point.

        :param split_point: splitting point.
        :param tol: tolerance.
        :return: list of two arcs.
        """
        if split_point.is_close(self.start, tol) or split_point.is_close(self.end, tol):
            raise ValueError("Point should be different of start and end.")
        if not self.point_belongs(split_point, 1e-5):
            raise ValueError("Point not on the circle.")
        return [Arc3D(self.circle, self.start, split_point),
                Arc3D(self.circle, split_point, self.end)]

    @classmethod
    def from_center_normal(cls, center: volmdlr.Point3D, normal: volmdlr.Vector3D,
                           start_end: volmdlr.Point3D, name: str = ''):
        u_vector = normal.deterministic_unit_normal_vector()
        v_vector = normal.cross(u_vector)
        circle = volmdlr_curves.Circle3D(volmdlr.Frame3D(center, u_vector, v_vector, normal),
                                         center.point_distance(start_end))
        return cls(circle, start_end, name=name)

    @classmethod
    def from_curve(cls, circle, start_end=None, name: str = ''):
        """
        Initialize a full arc from a circle.
        """
        if start_end is None:
            start_end = circle.center + circle.frame.u * circle.radius
        return cls(circle, start_end, name=name)


class ArcEllipse3D(Edge):
    """
    An arc is defined by a starting point, an end point and an interior point.

    """

    def __init__(self, ellipse: volmdlr_curves.Ellipse3D, start: volmdlr.Point3D, end: volmdlr.Point3D, name=''):
        Edge.__init__(self, start=start, end=end, name=name)
        self.ellipse = ellipse
        self.angle_start, self.angle_end = self.get_start_end_angles()
        self.angle = self.angle_end - self.angle_start
        self.center = ellipse.center
        self._self_2d = None
        self._length = None
        self._bbox = None

    def get_start_end_angles(self):
        """
        Calculate the start and end angles of the ArcEllipse3D in radians.

        This method computes the start and end angles of the ArcEllipse3D, which represent
        the angles, in radians, between the major axis of the ellipse and the start and end points
        on the ellipse's boundary.

        :return: tuple of floats
            A tuple containing the start and end angles in radians.
        """
        local_start_point = self.ellipse.frame.global_to_local_coordinates(self.start)
        u1, u2 = local_start_point.x / self.ellipse.major_axis, local_start_point.y / self.ellipse.minor_axis
        start_angle = volmdlr.geometry.sin_cos_angle(u1, u2)
        local_end_point = self.ellipse.frame.global_to_local_coordinates(self.end)
        u1, u2 = local_end_point.x / self.ellipse.major_axis, local_end_point.y / self.ellipse.minor_axis
        end_angle = volmdlr.geometry.sin_cos_angle(u1, u2)
        if math.isclose(end_angle, 0.0, abs_tol=1e-6):
            end_angle = volmdlr.TWO_PI
        return start_angle, end_angle

    @property
    def self_2d(self):
        if not self._self_2d:
            self._self_2d = self.to_2d(self.ellipse.center, self.ellipse.frame.u, self.ellipse.frame.v)
        return self._self_2d

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """
        Discretization of a Contour to have "n" points.

        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Useful
            to mesh an arc
        :return: a list of sampled points
        """
        if not number_points:
            if not angle_resolution:
                number_points = 2
            else:
                number_points = math.ceil(angle_resolution * abs(0.5 * self.angle / math.pi)) + 1
        angle_end = self.angle_end
        angle_start = self.angle_start
        if self.angle_start == self.angle_end:
            angle_end = angle_start + 2 * math.pi
        else:
            if angle_end < angle_start:
                angle_end = self.angle_end + volmdlr.TWO_PI

        discretization_points = [self.ellipse.frame.local_to_global_coordinates(
            volmdlr.Point3D(self.ellipse.major_axis * math.cos(angle),
                            self.ellipse.minor_axis * math.sin(angle), 0))
            for angle in npy.linspace(angle_start, angle_end, number_points)]
        return discretization_points

    def to_2d(self, plane_origin, x, y):
        """
        Transforms an Arc Ellipse 3D into an Arc Ellipse 2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: ArcEllipse2D.
        """
        point_start2d = self.start.to_2d(plane_origin, x, y)
        point_end2d = self.end.to_2d(plane_origin, x, y)
        ellipse2d = self.ellipse.to_2d(plane_origin, x, y)
        return ArcEllipse2D(ellipse2d, point_start2d, point_end2d)

    def length(self):
        """Computes the length."""
        if not self._length:
            self._length = self.self_2d.length()
        return self._length

    def normal_vector(self, abscissa):
        """Returns the normal vector at a given abscissa."""
        return self.direction_vector(abscissa).deterministic_normal_vector()

    def direction_vector(self, abscissa):
        """
        Returns the tangent vector at a given abscissa along the ArcEllipse3D.

        This method calculates and returns the tangent vector at a specific abscissa
        along the ArcEllipse3D, which represents the direction of the curve at that point.

        :param abscissa: The parameter value (abscissa) along the curve.
        :type abscissa: float

        :return: Vector3D
            A Vector3D object representing the tangent vector at the given abscissa.

        :raises:
            - ValueError: If the abscissa is out of the valid range of the curve.
        """
        if abscissa > self.length():
            raise ValueError('The abscissa is out of the valid range of the curve.')
        direction_vector_2d = self.self_2d.direction_vector(abscissa)
        direction_vector_3d = direction_vector_2d.to_3d(
            self.ellipse.center, self.ellipse.frame.u, self.ellipse.frame.v)
        return direction_vector_3d

    def abscissa(self, point: volmdlr.Point3D, tol: float = 1e-6):
        """
        Calculates the abscissa a given point.

        :param point: point to calculate abscissa.
        :param tol: tolerance allowed.
        :return: abscissa
        """
        if point.point_distance(self.start) < tol:
            return 0
        point2d = point.to_2d(self.ellipse.center, self.ellipse.major_dir, self.ellipse.minor_dir)
        return self.self_2d.abscissa(point2d)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plot the arc ellipse."""
        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')

        ax.plot([self.start[0]], [self.start[1]], [self.start[2]], c='r')
        ax.plot([self.end[0]], [self.end[1]], [self.end[2]], c='b')
        ax = vm_common_operations.plot_from_discretization_points(
            ax, edge_style=edge_style, element=self, number_points=25)
        if edge_style.edge_ends:
            self.start.plot(ax, 'r')
            self.end.plot(ax, 'b')
        return ax

    def plot2d(self, x3d: volmdlr.Vector3D = volmdlr.X3D, y3d: volmdlr.Vector3D = volmdlr.Y3D,
               ax=None, color='k'):
        """
        Plot 2d for an arc ellipse 3d.

        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        # TODO: Enhance this plot
        length = self.length()
        x = []
        y = []
        number_points = 30
        for i in range(number_points):
            point = self.point_at_abscissa(i / (number_points - 1) * length)
            xi, yi = point.plane_projection2d(x3d, y3d)
            x.append(xi)
            y.append(yi)
        ax.plot(x, y, color=color)
        return ax

    def triangulation(self):
        """
        Triangulation for an Arc Ellipse 3D.

        """
        return None

    @property
    def bounding_box(self):
        """
        Getter Bounding Box for an arc ellipse 3d.

        :return: bounding box.
        """
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        """
        Bounding Box setter.

        :param new_bounding_box: new bounding box.
        """
        self._bbox = new_bounding_box

    def get_bounding_box(self):
        """
        Calculates the bounding box of the Arc3D.

        :return: Bounding Box object.
        """
        # TODO: implement exact calculation

        points = self.discretization_points(angle_resolution=10)
        xmin = min(point.x for point in points)
        xmax = max(point.x for point in points)
        ymin = min(point.y for point in points)
        ymax = max(point.y for point in points)
        zmin = min(point.z for point in points)
        zmax = max(point.z for point in points)
        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Arc-Ellipse3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated Arc-Ellipse3D.
        """
        new_start = self.start.rotation(center, axis, angle)
        new_end = self.end.rotation(center, axis, angle)
        new_ellipse3d = self.ellipse.rotation(center, axis, angle)
        return ArcEllipse3D(new_ellipse3d, new_start, new_end)

    def translation(self, offset: volmdlr.Vector3D):
        """
        ArcEllipse3D translation.

        :param offset: translation vector.
        :return: A new translated ArcEllipse3D.
        """
        new_start = self.start.translation(offset)
        new_end = self.end.translation(offset)
        new_ellipse3d = self.ellipse.translation(offset)
        return ArcEllipse3D(new_ellipse3d, new_start, new_end)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new ArcEllipse3D.

        :param frame: Local coordinate system.
        :type frame: volmdlr.Frame3D
        :param side: 'old' will perform a transformation from local to global coordinates. 'new' will
            perform a transformation from global to local coordinates.
        :type side: str
        :return: A new transformed ArcEllipse3D.
        :rtype: ArcEllipse3D
        """
        return ArcEllipse3D(self.ellipse.frame_mapping(frame, side), self.start.frame_mapping(frame, side),
                            self.end.frame_mapping(frame, side))

    def point_belongs(self, point, abs_tol: float = 1e-6):
        """
        Verifies if a given point lies on the arc of ellipse 3D.

        :param point: point to be verified.
        :param abs_tol: Absolute tolerance to consider the point on the curve.
        :return: True is point lies on the arc of ellipse, False otherwise
        """
        point2d = point.to_2d(self.ellipse.center, self.ellipse.major_dir, self.ellipse.minor_dir)
        return self.self_2d.point_belongs(point2d, abs_tol=abs_tol)

    def is_close(self, other_edge, tol: float = 1e-6):
        """
        Checks if two arc-ellipse are the same considering the Euclidean distance.

        :param other_edge: other arc-ellipse.
        :param tol: The tolerance under which the Euclidean distance is considered equal to 0, defaults to 1e-6.
        :type tol: float, optional
        """

        if isinstance(other_edge, self.__class__):
            if (self.start.is_close(other_edge.start, tol) and self.end.is_close(other_edge.end, tol)
                    and self.ellipse.center.is_close(other_edge.ellipse.center, tol)
                    and self.point_belongs(other_edge.point_at_abscissa(other_edge.length() * 0.5), tol)):
                return True
        return False

    def complementary(self):
        """Gets the complementary arc of ellipse."""
        return self.__class__(self.ellipse, self.end, self.start)

    def point_at_abscissa(self, abscissa):
        """
        Calculates the point at a given abscissa.

        :param abscissa: abscissa to calculate point.
        :return: volmdlr.Point3D
        """
        point2d = self.self_2d.point_at_abscissa(abscissa)
        return point2d.to_3d(self.ellipse.center, self.ellipse.major_dir, self.ellipse.minor_dir)

    def split(self, split_point, tol: float = 1e-6):
        """
        Splits arc-ellipse at a given point.

        :param split_point: splitting point.
        :param tol: tolerance.
        :return: list of two Arc-Ellipse.
        """
        if split_point.is_close(self.start, tol):
            return [None, self.copy()]
        if split_point.is_close(self.end, tol):
            return [self.copy(), None]
        return [self.__class__(self.ellipse, self.start, split_point),
                self.__class__(self.ellipse, split_point, self.end)]

    def get_reverse(self):
        """Gets the same ellipse but in the reverse direction."""
        new_frame = volmdlr.Frame3D(self.ellipse.frame.origin, self.ellipse.frame.u, -self.ellipse.frame.v,
                                    self.ellipse.frame.u.cross(-self.ellipse.frame.v))
        ellipse3d = volmdlr_curves.Ellipse3D(self.ellipse.major_axis, self.ellipse.minor_axis, new_frame)
        return self.__class__(ellipse3d, self.end, self.start, self.name + '_reverse')

    def linesegment_intersections(self, linesegment, abs_tol: float = 1e-6):
        """
        Gets the intersections between an Ellipse 3D and a Line Segment 3D.

        :param linesegment: The other linesegment.
        :param abs_tol: The absolute tolerance.
        :return: A list with the intersections points between the two edges.
        """
        ellipse_linesegment_intersections = self.ellipse.linesegment_intersections(linesegment, abs_tol)
        intersections = []
        for intersection in ellipse_linesegment_intersections:
            if self.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections

    def arcellipse_intersections(self, arcellipse3d, abs_tol: float = 1e-6):
        """
        Gets the intersections between an Ellipse 3D and a Line Segment 3D.

        :param arcellipse3d: The other linesegment.
        :param abs_tol: The absolute tolerance.
        :return: A list with the intersections points between the two edges.
        """
        ellipse_intersections = self.ellipse.ellipse_intersections(arcellipse3d.ellipse, abs_tol)
        intersections = []
        for intersection in ellipse_intersections:
            if self.point_belongs(intersection, abs_tol) and arcellipse3d.point_belongs(intersection, abs_tol):
                intersections.append(intersection)
        return intersections


class FullArcEllipse3D(FullArcEllipse, ArcEllipse3D):
    """
    Defines a FullArcEllipse3D.
    """

    def __init__(self, ellipse: volmdlr_curves.Ellipse3D, start_end: volmdlr.Point3D, name: str = ''):
        self.ellipse = ellipse
        self.normal = self.ellipse.normal
        center2d = self.ellipse.center.to_2d(self.ellipse.center,
                                             self.ellipse.major_dir, self.ellipse.minor_dir)
        point_major_dir = self.ellipse.center + self.ellipse.major_axis * self.ellipse.major_dir
        point_major_dir_2d = point_major_dir.to_2d(
            self.ellipse.center, self.ellipse.major_dir, self.ellipse.minor_dir)
        vector_major_dir_2d = (point_major_dir_2d - center2d).to_vector()
        self.theta = volmdlr.geometry.clockwise_angle(vector_major_dir_2d, volmdlr.X2D)
        if self.theta == math.pi * 2:
            self.theta = 0.0
        self._bbox = None

        FullArcEllipse.__init__(self, self.ellipse, start_end, name)
        ArcEllipse3D.__init__(self, self.ellipse, start_end, start_end)

    def to_dict(self, use_pointers: bool = False, memo=None, path: str = '#'):
        dict_ = self.base_dict()
        dict_["ellipse"] = self.ellipse.to_dict(use_pointers=use_pointers, memo=memo, path=path + '/ellipse')
        dict_['start_end'] = self.start_end.to_dict(use_pointers=use_pointers, memo=memo, path=path + '/start_end')
        return dict_

    @classmethod
    def dict_to_object(cls, dict_, *args, **kwargs):
        """
        Create a FullArcEllipse3D object from a dictionary representation.

        This class method takes a dictionary containing the necessary data for
        creating a FullArcEllipse3D object and returns an instance of the FullArcEllipse3D class.
        It expects the dictionary to have the following keys:

        :param cls: The FullArcEllipse3D class itself (automatically passed).
        :param dict_: A dictionary containing the required data for object creation.
        :param args: Additional positional arguments (if any).
        :param kwargs: Additional keyword arguments (if any).

        :return: FullArcEllipse3D: An instance of the FullArcEllipse3D class created from the provided dictionary.
        """
        ellipse = volmdlr_curves.Ellipse3D.dict_to_object(dict_['ellipse'])
        start_end = volmdlr.Point3D.dict_to_object(dict_['start_end'])

        return cls(ellipse, start_end, name=dict_['name'])

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a FullArcEllipse3D into an FullArcEllipse2D, given an plane origin and a u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: FullArcEllipse2D.
        """
        point_start_end2d = self.start_end.to_2d(plane_origin, x, y)
        ellipse2d = self.ellipse.to_2d(plane_origin, x, y)
        return FullArcEllipse2D(ellipse2d, point_start_end2d, name=self.name)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new FullArcEllipse3D.

        :param frame: Local coordinate system.
        :type frame: volmdlr.Frame3D
        :param side: 'old' will perform a transformation from local to global coordinates. 'new' will
            perform a transformation from global to local coordinates.
        :type side: str
        :return: A new transformed FulLArcEllipse3D.
        :rtype: FullArcEllipse3D
        """
        return FullArcEllipse3D(self.ellipse.frame_mapping(frame, side),
                                self.start_end.frame_mapping(frame, side), name=self.name)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Ellipse3D translation.

        :param offset: translation vector.
        :type offset: volmdlr.Vector3D
        :return: A new translated FullArcEllipse3D.
        :rtype: FullArcEllipse3D
        """
        return FullArcEllipse3D(self.ellipse.translation(offset), self.start_end.translation(offset), self.name)

    def abscissa(self, point: volmdlr.Point3D, tol: float = 1e-6):
        """
        Calculates the abscissa a given point.

        :param point: point to calculate abscissa.
        :param tol: tolerance allowed.
        :return: abscissa
        """
        point2d = point.to_2d(self.ellipse.center, self.ellipse.major_dir, self.ellipse.minor_dir)
        return self.self_2d.abscissa(point2d, tol=tol)

    def split(self, split_point, tol: float = 1e-6):
        """
        Splits the ellipse into two arc of ellipse at a given point.

        :param split_point: splitting point.
        :param tol: tolerance.
        :return: list of two Arc of ellipse.
        """
        if split_point.is_close(self.start, tol) or split_point.is_close(self.end, tol):
            return [self, None]
        if not self.point_belongs(split_point, 1e-5):
            raise ValueError("Point not on the ellipse.")
        return [ArcEllipse3D(self.ellipse, self.start_end, split_point),
                ArcEllipse3D(self.ellipse, split_point, self.start_end)]

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Ellipse plot."""
        return self.ellipse.plot(ax, edge_style)
