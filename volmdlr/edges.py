#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Edges related classes.
"""

import math
import sys
import warnings
from typing import Any, Dict, List, Union

import dessia_common.core as dc
import matplotlib.patches
import matplotlib.pyplot as plt
import numpy as npy
import plot_data.core as plot_data
import scipy.integrate as scipy_integrate
from scipy.optimize import least_squares, minimize, lsq_linear
from geomdl import NURBS, BSpline, fitting, operations, utilities
from geomdl.operations import length_curve, split_curve
from matplotlib import __version__ as _mpl_version
from mpl_toolkits.mplot3d import Axes3D
from packaging import version

import volmdlr.core
import volmdlr.core_compiled
import volmdlr.geometry
import volmdlr.utils.common_operations as vm_common_operations
import volmdlr.utils.intersections as vm_utils_intersections
from volmdlr import bspline_fitting
from volmdlr.core import EdgeStyle


def standardize_knot_vector(knot_vector):
    """
    Standardize a knot vector to range from 0 to 1.
    """
    first_knot = knot_vector[0]
    last_knot = knot_vector[-1]
    standard_u_knots = []
    if first_knot != 0 or last_knot != 1:
        x = 1 / (last_knot - first_knot)
        y = first_knot / (first_knot - last_knot)
        for u in knot_vector:
            standard_u_knots.append(u * x + y)
        return standard_u_knots
    return knot_vector


def insert_knots_and_mutiplicity(knots, knot_mutiplicities, knot_to_add, num):
    """
    Compute knot-elements and multiplicities based on the global knot vector.

    """
    new_knots = []
    new_knot_mutiplicities = []
    i = 0
    for i, knot in enumerate(knots):
        if knot > knot_to_add:
            new_knots.extend([knot_to_add])
            new_knot_mutiplicities.append(num)
            new_knots.extend(knots[i:])
            new_knot_mutiplicities.extend(knot_mutiplicities[i:])
            break
        new_knots.append(knot)
        new_knot_mutiplicities.append(knot_mutiplicities[i])
    return new_knots, new_knot_mutiplicities, i


class Edge(dc.DessiaObject):
    """
    Defines a simple edge Object.
    """

    def __init__(self, start, end, name=''):
        self.start = start
        self.end = end
        self._length = None
        self._direction_vector = None

        # Disabling super init call for performance
        # dc.DessiaObject.__init__(self, name=name)
        self.name = name

    def __getitem__(self, key):
        if key == 0:
            return self.start
        if key == 1:
            return self.end
        raise IndexError

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
        half_length = self.length() / 2
        middle_point = self.point_at_abscissa(abscissa=half_length)
        return middle_point

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = None):
        """
        Discretize an Edge to have "n" points.

        :param number_points: the number of points (including start and end
            points) if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have
            a controlled angular distance. Useful to mesh an arc
        :return: a list of sampled points
        """
        if number_points is None or number_points == 1:
            number_points = 2
        if angle_resolution:
            number_points = int(math.pi * angle_resolution)
        step = self.length() / (number_points - 1)
        return [self.point_at_abscissa(i * step) for i in range(number_points)]

    def polygon_points(self, discretization_resolution: int):
        """
        Deprecated method of discretization_points.
        """
        warnings.warn('polygon_points is deprecated,\
        please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(discretization_resolution)

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
        # obj can be an instance of wires or edges.
        obj = object_dict[arguments[3]]
        point1 = object_dict[arguments[1]]
        point2 = object_dict[arguments[2]]
        orientation = arguments[4]
        if orientation == '.F.':
            point1, point2 = point2, point1
        if obj.__class__.__name__ == 'LineSegment3D':
            return object_dict[arguments[3]]
        if obj.__class__.__name__ == 'Line3D':
            if not point1.is_close(point2):
                return LineSegment3D(point1, point2, arguments[0][1:-1])
            return None
        if hasattr(obj, 'trim'):
            if obj.__class__.__name__ == 'Circle3D':
                point1, point2 = point2, point1
            return obj.trim(point1, point2)

        raise NotImplementedError(f'Unsupported: {object_dict[arguments[3]]}')

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
        vector = self.normal_vector(abscissa)
        vector.normalize()
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
        vector = self.direction_vector(abscissa)
        vector.normalize()
        return vector

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point: Point to be verified
        :return: Return True if the point belongs to this surface,
            or False otherwise
        """
        raise NotImplementedError(f'the straight_line_point_belongs method must be'
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

    def local_discretization(self, point1, point2, number_points):
        """
        Gets n discretization points between two given points of the edge.

        :param point1: point 1 on edge.
        :param point2: point 2 on edge.
        :param number_points: number of points to discretize locally.
        :return: list of locally discretized points.
        """
        abscissa1 = self.abscissa(point1)
        abscissa2 = self.abscissa(point2)
        discretized_points_between_1_2 = [self.point_at_abscissa(abscissa) for abscissa
                                          in npy.linspace(abscissa1, abscissa2, num=number_points)]
        return discretized_points_between_1_2

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
            if split_edge.point_belongs(point1, 1e-4) and split_edge.point_belongs(point2, 1e-4):
                new_split_edge = split_edge
                break
        return new_split_edge


class Line(dc.DessiaObject):
    """
    Abstract class representing a line.

    :param point1: The first point defining the line
    :type point1: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
    :param point2: The second point defining the line
    :type point2: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
    :param name: Name of the line. Default value is an empty string
    :type name: str, optional
    """

    def __init__(self, point1, point2, name=''):
        self.point1 = point1
        self.point2 = point2
        self._direction_vector = None
        dc.DessiaObject.__init__(self, name=name)

    def __getitem__(self, key):
        """
        Get a point of the line by its index.
        """
        if key == 0:
            return self.point1
        if key == 1:
            return self.point2
        raise IndexError

    def unit_direction_vector(self, *args, **kwargs):
        """
        Get the unit direction vector of the line.

        :return: The unit direction vector of the line
        :rtype:  Union[:class:`volmdlr.Vector2D`, :class:`volmdlr.Vector3D`]
        """
        vector = self.direction_vector()
        vector.normalize()
        return vector

    def direction_vector(self, *args, **kwargs):
        """
        Get the direction vector of the line.

        :return: The direction vector of the line
        :rtype: Union[:class:`volmdlr.Vector2D`, :class:`volmdlr.Vector3D`]
        """
        if not self._direction_vector:
            self._direction_vector = self.point2 - self.point1
        return self._direction_vector

    def normal_vector(self, *args, **kwargs):
        """
        Get the normal vector of the line.

        :return: The normal vector of the line
        :rtype: Union[:class:`volmdlr.Vector2D`, :class:`volmdlr.Vector3D`]
        """
        return self.direction_vector().normal_vector()

    def unit_normal_vector(self, *args, **kwargs):
        """
        Get the unit normal vector of the line.

        :return: The unit normal vector of the line
        :rtype: Union[:class:`volmdlr.Vector2D`, :class:`volmdlr.Vector3D`]
        """
        return self.unit_direction_vector().normal_vector()

    def point_projection(self, point):
        """
        Calculate the projection of a point onto the line.

        :param point: The point to project
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :return: The projection of the point onto the line and the distance
            between the point and the projection
        :rtype: Tuple(Union[:class:`volmdlr.Point2D`,
            :class:`volmdlr.Point3D`], float)
        """
        vector = self.point2 - self.point1
        norm_u = vector.norm()
        t = (point - self.point1).dot(vector) / norm_u ** 2
        projection = self.point1 + t * vector
        projection = projection.to_point()
        return projection, t * norm_u

    def abscissa(self, point):
        """
        Calculate the abscissa of a point on the line.

        :param point: The point for which to calculate the abscissa
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :return: The abscissa of the point
        :rtype: float
        """
        vector = self.point2 - self.point1
        norm_u = vector.norm()
        t_param = (point - self.point1).dot(vector) / norm_u
        return t_param

    def point_at_abscissa(self, abscissa):
        """
        Returns the point that corresponds to the given abscissa.

        :param abscissa: The abscissa
        :type abscissa: float
        :return: The point that corresponds to the given abscissa.
        :rtype: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        """
        return self.point1 + (self.point2 - self.point1) * abscissa

    def sort_points_along_line(self, points):
        """
        Sort point along a line.

        :param points: list of points to be sorted.
        :return: sorted points.
        """
        return sorted(points, key=self.abscissa)

    def split(self, split_point):
        """
        Split a line into two lines.

        :param split_point: The point where to split the line
        :type split_point: Union[:class:`volmdlr.Point2D`,
            :class:`volmdlr.Point3D`]
        :return: A list containing two lines
        """
        return [self.__class__(self.point1, split_point),
                self.__class__(split_point, self.point2)]

    def is_between_points(self, point1: Union[volmdlr.Point2D, volmdlr.Point3D],
                          point2: Union[volmdlr.Point2D, volmdlr.Point3D]):
        """
        Verifies if a line is between two points.

        :param point1: The first point
        :type point1: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :param point2: The second point
        :type point2: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :return: True if the line is between the two points, False otherwise
        :rtype: bool
        """

        if point1.is_close(point2):
            return False

        line_segment = LineSegment2D(point1, point2)
        if line_segment.line_intersections(self):
            return True
        return False

    def to_step(self, current_id, *args, **kwargs):
        """Exports to STEP format."""
        p1_content, p1_id = self.point1.to_step(current_id)
        # p2_content, p2_id = self.point2.to_step(current_id+1)
        current_id = p1_id + 1
        u_content, u_id = self.unit_direction_vector().to_step(current_id)
        current_id = u_id + 1
        content = p1_content + u_content
        content += f"#{current_id} = LINE('{self.name}',#{p1_id},#{u_id});\n"
        return content, [current_id]


class LineSegment(Edge):
    """
    Abstract class.

    """

    def direction_independent_eq(self, linesegment2):
        """
        Verifies if two line segments are the same, not considering its direction.

        """
        if self.start.is_close(linesegment2.start) and self.end.is_close(linesegment2.end):
            return True
        return self.start.is_close(linesegment2.end) and self.end.is_close(linesegment2.start)

    def length(self):
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
        t = (point - self.start).dot(vector) / length
        if t < -1e-9 or t > length + 1e-9:
            raise ValueError(f'Point is not on linesegment: abscissa={t}')
        return t

    def direction_vector(self, abscissa=0.):
        """
        Returns a direction vector at a given abscissa, it is not normalized.

        :param abscissa: defines where in the line segment
            direction vector is to be calculated.
        :return: The direction vector of the LineSegment.
        """
        if not self._direction_vector:
            self._direction_vector = self.end - self.start
        return self._direction_vector

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

    def split(self, split_point):
        """
        Split a Line Segment at a given point into two Line Segments.

        :param split_point: splitting point.
        :return: list with the two split line segments.
        """
        if split_point.is_close(self.start):
            return [None, self.copy()]
        if split_point.is_close(self.end):
            return [self.copy(), None]
        return [self.__class__(self.start, split_point),
                self.__class__(split_point, self.end)]

    def middle_point(self):
        """
        Calculates the middle point of a Line Segment.

        :return:
        """
        return 0.5 * (self.start + self.end)

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
        return [self.start, self.end]

    def get_shared_section(self, other_linesegment):
        """
        Gets the shared section between two line segments.

        :param other_linesegment: other line segment to verify for shared section.
        :return: shared line segment section.
        """
        if self.__class__ != other_linesegment.__class__:
            return []
        if not self.direction_vector().is_colinear_to(other_linesegment.direction_vector()) or \
                (not any(self.point_belongs(point, 1e-6)
                         for point in [other_linesegment.start, other_linesegment.end]) and
                 not any(other_linesegment.point_belongs(point, 1e-6) for point in [self.start, self.end])):
            return []
        if all(self.point_belongs(point) for point in other_linesegment.discretization_points(number_points=5)):
            return [other_linesegment]
        if all(other_linesegment.point_belongs(point) for point in self.discretization_points(number_points=5)):
            return [self]
        new_linesegment_points = []
        for point in [self.start, self.end]:
            if other_linesegment.point_belongs(point, abs_tol=1e-6) and\
                    not volmdlr.core.point_in_list(point, new_linesegment_points):
                new_linesegment_points.append(point)
        for point in [other_linesegment.start, other_linesegment.end]:
            if self.point_belongs(point, abs_tol=1e-6) and\
                    not volmdlr.core.point_in_list(point, new_linesegment_points):
                new_linesegment_points.append(point)
        if len(new_linesegment_points) == 1:
            return []
        if len(new_linesegment_points) != 2:
            raise ValueError
        class_ = self.__class__
        return [class_(new_linesegment_points[0], new_linesegment_points[1])]

    def delete_shared_section(self, other_linesegment):
        """
        Deletes from self, the section shared with the other line segment.

        :param other_linesegment:
        :return:
        """
        shared_section = self.get_shared_section(other_linesegment)
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
            if not lineseg.direction_independent_eq(shared_section[0]):
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
        line = self.to_line()
        content, (line_id,) = line.to_step(current_id)
        current_id = line_id + 1
        start_content, start_id = self.start.to_step(current_id, vertex=True)
        current_id = start_id + 1
        end_content, end_id = self.end.to_step(current_id + 1, vertex=True)
        content += start_content + end_content
        current_id = end_id + 1
        content += f"#{current_id} = EDGE_CURVE('{self.name}',#{start_id},#{end_id},#{line_id},.T.);\n"
        return content, [current_id]


class BSplineCurve(Edge):
    """
    An abstract class for B-spline curves.

    The following rule must be
    respected : `number of knots = number of control points + degree + 1`.

    :param degree: The degree of the B-spline curve.
    :type degree: int
    :param control_points: A list of 2 or 3 dimensional points
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
    _non_serializable_attributes = ['curve']

    def __init__(self,
                 degree: int,
                 control_points: Union[List[volmdlr.Point2D], List[volmdlr.Point3D]],
                 knot_multiplicities: List[int],
                 knots: List[float],
                 weights: List[float] = None,
                 periodic: bool = False,
                 name: str = ''):
        self.control_points = control_points
        self.degree = degree
        knots = standardize_knot_vector(knots)
        self.knots = knots
        self.knot_multiplicities = knot_multiplicities
        self.weights = weights
        self.periodic = periodic

        points = [[*point] for point in control_points]
        if weights is None:
            curve = BSpline.Curve()
            curve.degree = degree
            curve.ctrlpts = points
        else:
            curve = NURBS.Curve()
            curve.degree = degree
            curve.ctrlpts = points
            curve.weights = weights

        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot] * knot_multiplicities[i])
        curve.knotvector = knot_vector
        curve.delta = 0.01
        curve_points = curve.evalpts
        self.curve = curve

        self._length = None
        self.points = [getattr(volmdlr,
                               f'Point{self.__class__.__name__[-2::]}')(*point)
                       for point in curve_points]

        Edge.__init__(self, self.points[0], self.points[-1], name=name)

    def to_dict(self, *args, **kwargs):
        """Avoids storing points in memo that makes serialization slow."""
        dict_ = self.base_dict()
        dict_['degree'] = self.degree
        dict_['control_points'] = [point.to_dict() for point in self.control_points]
        dict_['knot_multiplicities'] = self.knot_multiplicities
        dict_['knots'] = self.knots
        dict_['weights'] = self.weights
        dict_['periodic'] = self.periodic
        return dict_

    def __hash__(self):
        """
        Return a hash value for the B-spline curve.
        """
        return hash((tuple(self.control_points), self.degree, tuple(self.knots)))

    def __eq__(self, other):
        """
        Return True if the other B-spline curve has the same control points, degree, and knot vector, False otherwise.
        """
        if isinstance(other, self.__class__):
            return (self.control_points == other.control_points
                    and self.degree == other.degree
                    and self.knots == other.knots)
        return False

    def reverse(self):
        """
        Reverses the B-spline's direction by reversing its control points.

        :return: A reversed B-spline curve.
        :rtype: :class:`volmdlr.edges.BSplineCurve`.
        """
        return self.__class__(
            degree=self.degree,
            control_points=self.control_points[::-1],
            knot_multiplicities=self.knot_multiplicities[::-1],
            knots=self.knots[::-1],
            weights=self.weights,
            periodic=self.periodic)

    @classmethod
    def from_geomdl_curve(cls, curve, name: str = ""):
        """
        # TODO: to be completed.

        :param curve:
        :type curve:
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
                   knot_multiplicities=knot_multiplicities, name=name)

    def length(self):
        """
        Returns the length of the B-spline curve.

        :return: The length of the B-spline curve.
        :rtype: float
        """
        if not self._length:
            self._length = length_curve(self.curve)
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
        u = abscissa / self.length()
        derivatives = self.derivatives(u, 1)
        return derivatives[1]

    def middle_point(self):
        """
        Computes the 2D or 3D middle point of the B-spline curve.

        :return: The middle point
        :rtype: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        """
        return self.point_at_abscissa(self.length() * 0.5)

    def abscissa(self, point: Union[volmdlr.Point2D, volmdlr.Point3D],
                 tol: float = 1e-6):
        """
        Computes the abscissa of a 2D or 3D point using the least square method.

        :param point: The point located on the B-spline curve.
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`].
        :param tol: The precision in terms of distance. Default value is 1e-6.
        :type tol: float, optional.
        :return: The abscissa of the point.
        :rtype: float
        """
        length = self.length()
        initial_condition_list = [0, 0.25, 0.5, 0.75, 1]

        def evaluate_point_distance(u):
            return (point - self.evaluate_single(u)).norm()
        results = []
        initial_condition_list.sort(key=evaluate_point_distance)
        for u0 in initial_condition_list:
            u, convergence_sucess = self.point_invertion(u0, point)
            abscissa = u * length
            if convergence_sucess:  # sometimes we don't achieve convergence with a given initial guess
                return abscissa
            dist = evaluate_point_distance(u)
            if dist < tol:
                return abscissa
            results.append((abscissa, dist))
        result = min(results, key=lambda r: r[1])[0]
        return result

    def _point_inversion_funcs(self, u, point):
        """
        Helper function to evaluate Newton-Rapshon terms.
        """
        curve_derivatives = self.derivatives(u, 2)
        distance_vector = curve_derivatives[0] - point
        func = curve_derivatives[1].dot(distance_vector)
        func_first_derivative = curve_derivatives[2].dot(distance_vector) + curve_derivatives[1].norm() ** 2
        return func, func_first_derivative, curve_derivatives, distance_vector

    def point_invertion(self, u0: float, point, maxiter: int = 50, tol1: float = 1e-6, tol2: float = 1e-6):
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
        new_u = u0 - func / func_first_derivative
        new_u = self._check_bounds(new_u)
        residual = (new_u - u0) * curve_derivatives[1]
        if residual.norm() <= 1e-6:
            return u0, False
        u0 = new_u
        return self.point_invertion(u0, point, maxiter=maxiter - 1)

    @staticmethod
    def _check_convergence(curve_derivatives, distance_vector, tol1: float = 1e-6, tol2: float = 1e-6):
        """
        Helper function to check convergence of point_invertion method.
        """
        distance = distance_vector.norm()
        if distance <= tol1:
            return True
        zero_cos = abs(curve_derivatives[1].dot(distance_vector)) / curve_derivatives[1].norm() * distance
        if distance <= tol1 and zero_cos <= tol2:
            return True
        return False

    def _check_bounds(self, u):
        """
        Helper function to check if evaluated parameters in point_invertion method are contained in the bspline domain.
        """
        a, b = self.curve.domain
        if self.periodic:
            if u < a:
                u = b - (a - u)
            elif u > b:
                u = a + (u - b)
        else:
            if u < a:
                u = a

            elif u > b:
                u = b
        return u

    def split(self, point: Union[volmdlr.Point2D, volmdlr.Point3D],
              tol: float = 1e-5):
        """
        Splits of B-spline curve in two pieces using a 2D or 3D point.

        :param point: The point where the B-spline curve is split
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :param tol: The precision in terms of distance. Default value is 1e-4
        :type tol: float, optional
        :return: A list containing the first and second split of the B-spline
            curve
        :rtype: List[:class:`volmdlr.edges.BSplineCurve`]
        """
        if point.point_distance(self.start) < tol:
            return [None, self.copy()]
        if point.point_distance(self.end) < tol:
            return [self.copy(), None]
        adim_abscissa = self.abscissa(point) / self.length()
        curve1, curve2 = split_curve(self.curve, adim_abscissa)

        return [self.__class__.from_geomdl_curve(curve1),
                self.__class__.from_geomdl_curve(curve2)]

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
                              self.weights, self.periodic)

    def translation_inplace(self, offset: Union[volmdlr.Vector2D, volmdlr.Vector3D]):
        """
        Translates the B-spline curve and its parameters are modified inplace.

        :param offset: The translation vector
        :type offset: Union[:class:`volmdlr.Vector2D`,
            :class:`volmdlr.Vector3D`]
        :return: None
        :rtype: None
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in self.control_points:
            point.translation_inplace(offset)

    def point_belongs(self, point: Union[volmdlr.Point2D, volmdlr.Point3D],
                      abs_tol: float = 1e-6):
        """
        Checks if a 2D or 3D point belongs to the B-spline curve or not. It uses the least square method.

        :param point: The point to be checked
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :param abs_tol: The precision in terms of distance.
            Default value is 1e-4
        :type abs_tol: float, optional
        :return: `True` if the point belongs to the B-spline curve, `False`
            otherwise
        :rtype: bool
        """
        point_dimension = f'Point{self.__class__.__name__[-2::]}'

        def fun(x):
            return (point - getattr(volmdlr, point_dimension)(*self.curve.evaluate_single(x))).norm()

        x = npy.linspace(0, 1, 5)
        x_init = []
        for xi in x:
            x_init.append(xi)

        for x0 in x_init:
            z = least_squares(fun, x0=x0, bounds=([0, 1]))
            if z.fun < abs_tol:
                return True
        return False

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

        return self.__class__.from_points_interpolation(
            points, min(self.degree, bspline_curve.degree))

    @classmethod
    def from_bsplines(cls, bsplines: List['BSplineCurve'],
                      discretization_points: int = 10):
        """
        Creates a B-spline curve from a list of B-spline curves.

        :param bsplines: A list of B-spline curve
        :type bsplines: List[:class:`volmdlr.edges.BSplineCurve`]
        :param discretization_points: The number of points for the
            discretization. Default value is 10
        :type discretization_points: int, optional
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

        return cls.from_points_interpolation(points, min(degree))

    @classmethod
    def from_points_approximation(cls, points: Union[List[volmdlr.Point2D], List[volmdlr.Point3D]],
                                  degree: int, **kwargs):
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
        :param kwargs: See below
        :return: A B-spline curve from points approximation
        :rtype: :class:`volmdlr.edges.BSplineCurve`
        :keyword centripetal: Activates centripetal parametrization method.
            Default value is False
        :keyword ctrlpts_size: Number of control points. Default value is
            len(points) - 1
        """
        curve = fitting.approximate_curve([[*point] for point in points],
                                          degree, **kwargs)
        return cls.from_geomdl_curve(curve)

    def tangent(self, position: float = 0.0):
        """
        Evaluates the tangent vector of the B-spline curve at the input parameter value.

        :param position: Value of the parameter, between 0 and 1
        :type position: float
        :return: The tangent vector
        :rtype: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        """
        _, tangent = operations.tangent(self.curve, position, normalize=True)

        dimension = f'Vector{self.__class__.__name__[-2::]}'
        tangent = getattr(volmdlr, dimension)(*tangent)

        return tangent

    @classmethod
    def from_points_interpolation(cls, points: Union[List[volmdlr.Point2D], List[volmdlr.Point3D]],
                                  degree: int, periodic: bool = False, name: str = ""):
        """
        Creates a B-spline curve interpolation through the data points.

        Please refer to Algorithm A9.1 on The NURBS Book (2nd Edition),
        pp.369-370 for details.

        :param points: The data points
        :type points: Union[List[:class:`volmdlr.Point2D`],
            List[:class:`volmdlr.Point3D`]]
        :param degree: The degree of the output parametric curve
        :type degree: int
        :param periodic: `True` if the curve should be periodic. Default value
            is `False`
        :type periodic: bool, optional
        :return: A B-spline curve from points interpolation
        :rtype: :class:`volmdlr.edges.BSplineCurve`
        """
        curve = bspline_fitting.interpolate_curve([[*point] for point in points], degree, centripetal=True)

        bsplinecurve = cls.from_geomdl_curve(curve, name=name)
        if not periodic:
            return bsplinecurve
        bsplinecurve.periodic = True
        return bsplinecurve

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

        if len(self.points) == number_points or (not number_points and not angle_resolution):
            return self.points
        curve = self.curve
        curve.delta = 1 / number_points
        curve_points = curve.evalpts

        point_dimension = f'Point{self.__class__.__name__[-2::]}'
        return [getattr(volmdlr, point_dimension)(*point) for point in curve_points]

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

        return [getattr(volmdlr, f'Vector{self.__class__.__name__[-2::]}')(*point)
                for point in self.curve.derivatives(u, order)]

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

    def line_intersections(self, line):
        """
        Calculates the intersections of a BSplineCurve (2D or 3D) with a Line (2D or 3D).

        :param line: line to verify intersections
        :return: list of intersections
        """
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
                if line.point_distance(linesegment.middle_point()) < 1e-8:
                    list_intersections.append(linesegment.middle_point())
            if intersections and intersections[0] not in list_intersections:
                if self.point_belongs(intersections[0], 1e-6):
                    list_intersections.append(intersections[0])
                    continue
                abs1 = self.abscissa(linesegment.start)
                abs2 = self.abscissa(linesegment.end)
                list_abscissas = list(new_abscissa for new_abscissa in npy.linspace(abs1, abs2, 1000))
                intersection = self.select_intersection_point(list_abscissas, intersections)
                list_intersections.append(intersection)
            initial_abscissa += linesegment.length()
        return list_intersections

    def select_intersection_point(self, list_abscissas, intersections):
        """
        Select closest point in curve to intersection point obtained with discretized linesegment.

        :param list_abscissas: list of abscissas to verify the closest point.
        :param intersections: intersection with discretised line.
        :return:
        """
        distance = npy.inf
        intersection = None
        for i_abscissa in list_abscissas:
            point_in_curve = BSplineCurve.point_at_abscissa(self, i_abscissa)
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
        results = self.line_intersections(linesegment.to_line())
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
        length = self.length()
        adim_abs = max(min(abscissa / length, 1.), 0.)
        point_name = 'Point' + self.__class__.__name__[-2:]
        return getattr(volmdlr, point_name)(*self.curve.evaluate_single(adim_abs))

    def get_shared_section(self, other_bspline2):
        """
        Gets the shared section between two BSpline curves.

        :param other_bspline2: other arc to verify for shared section.
        :return: shared arc section.
        """
        if self.__class__ != other_bspline2.__class__:
            return []
        if self.__class__.__name__[-2:] == '3D':
            if self.bounding_box.distance_to_bbox(other_bspline2.bounding_box) > 1e-7:
                return []
        elif self.bounding_rectangle.distance_to_b_rectangle(other_bspline2.bounding_rectangle) > 1e-7:
            return []
        if not any(self.point_belongs(point, abs_tol=1e-6)
                   for point in other_bspline2.discretization_points(number_points=10)):
            return []
        if all(self.point_belongs(point, abs_tol=1e-6) for point in other_bspline2.points):
            return [other_bspline2]
        if all(other_bspline2.point_belongs(point, abs_tol=1e-6) for point in self.points):
            return [self]
        if self.point_belongs(other_bspline2.start, abs_tol=1e-6):
            bspline1_, bspline2_ = self.split(other_bspline2.start)
        elif self.point_belongs(other_bspline2.end, abs_tol=1e-6):
            bspline1_, bspline2_ = self.split(other_bspline2.end)
        else:
            raise NotImplementedError
        shared_bspline_section = []
        for bspline in [bspline1_, bspline2_]:
            if bspline and all(other_bspline2.point_belongs(point)
                               for point in bspline.discretization_points(number_points=10)):
                shared_bspline_section.append(bspline)
                break
        return shared_bspline_section

    def delete_shared_section(self, other_bspline2):
        """
        Deletes from self, the section shared with the other arc.

        :param other_bspline2:
        :return:
        """
        shared_section = self.get_shared_section(other_bspline2)
        if not shared_section:
            return [self]
        if shared_section == self:
            return []
        split_bspline1 = self.split(shared_section[0].start)
        split_bspline2 = self.split(shared_section[0].end)
        new_arcs = []
        shared_section_middle_point = shared_section[0].point_at_abscissa(0.5 * shared_section[0].length())
        for arc in split_bspline1 + split_bspline2:
            if arc and not arc.point_belongs(shared_section_middle_point, abs_tol=1e-6):
                new_arcs.append(arc)
        return new_arcs

    def evaluate_single(self, u):
        """
        Calculates a point in the BSplineCurve at a given parameter u.

        :param u: Curve parameter. Must be a value between 0 and 1.
        :type u: float
        :return: Corresponding point.
        :rtype: Union[volmdlr.Point2D, Union[volmdlr.Point3D]
        """
        point_name = 'Point' + self.__class__.__name__[-2:]
        return getattr(volmdlr, point_name)(*self.curve.evaluate_single(u))

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point: Point to be verified
        :return: Return True if the point belongs to this surface,
            or False otherwise
        """
        raise NotImplementedError(f'the straight_line_point_belongs method must be'
                                  f' overloaded by {self.__class__.__name__}')


class Line2D(Line):
    """
    Define an infinite line given by two points.

    """

    def __init__(self, point1: volmdlr.Point2D,
                 point2: volmdlr.Point2D, *, name=''):
        # self.points = [point1, point2]
        Line.__init__(self, point1, point2, name=name)

    def to_3d(self, plane_origin, x1, x2):
        """
        Convert the line to a 3D line.

        :param plane_origin: Origin of the plane in which the line is.
        :type plane_origin: :class:`volmdlr.Point3D`
        :param x1: First direction of the plane in which the line is.
        :type x1: :class:`volmdlr.Vector3D`
        :param x2: Second direction of the plane in which the line is.
        :type x2: :class:`volmdlr.Vector3D`
        :return: The 3D line.
        :rtype: :class:`volmdlr.edges.Line3D`
        """
        points_3d = [point.to_3d(plane_origin, x1, x2) for point in [self.point1, self.point2]]
        return Line3D(*points_3d, self.name)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Line2D rotation.

        :param center: rotation center.
        :param angle: angle rotation.
        :return: a new rotated Line2D.
        """
        return Line2D(*[point.rotation(center, angle)
                        for point in [self.point1, self.point2]])

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        Line2D rotation. Object is updated inplace.

        :param center: rotation center.
        :param angle: rotation angle.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in [self.point1, self.point2]:
            point.rotation_inplace(center, angle)

    def translation(self, offset: volmdlr.Vector2D):
        """
        Line2D translation.

        :param offset: translation vector.
        :return: A new translated Line2D.
        """
        return Line2D(*[point.translation(offset) for point in [self.point1, self.point2]])

    def translation_inplace(self, offset: volmdlr.Vector2D):
        """
        Line2D translation. Object is updated inplace.

        :param offset: translation vector.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in [self.point1, self.point2]:
            point.translation_inplace(offset)

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Map the line to a new coordinate frame.

        :param frame: The new coordinate frame.
        :type frame: :class:`volmdlr.Frame2D`
        :param side: The side to which the mapping is made. 'old' for the
            original coordinate frame, 'new' for the new one.
        :type side: str
        :return: The mapped line.
        :rtype: :class:`volmdlr.edges.Line2D`
        """
        return Line2D(*[point.frame_mapping(frame, side) for point in [self.point1, self.point2]])

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Plot the line.

        :param ax: Matplotlib axis on which to plot the line. If none,
            a new figure is created.
        :type ax: matplotlib.axes._subplots.AxesSubplot, optional
        :param edge_style: data class instance, containing all parameters needed to plot Line 2D.
        :return: The matplotlib axis.
        :rtype: matplotlib.axes._subplots.AxesSubplot
        """
        if ax is None:
            _, ax = plt.subplots()

        if version.parse(_mpl_version) >= version.parse('3.3.2'):
            if edge_style.dashed:
                ax.axline((self.point1.x, self.point1.y),
                          (self.point2.x, self.point2.y),
                          dashes=[30, 5, 10, 5],
                          color=edge_style.color)
            else:
                ax.axline((self.point1.x, self.point1.y),
                          (self.point2.x, self.point2.y),
                          color=edge_style.color)
        else:
            direction_vector = self.direction_vector()
            point3 = self.point1 - 3 * direction_vector
            point4 = self.point2 + 4 * direction_vector
            if edge_style.dashed:
                ax.plot([point3[0], point4[0]], [point3[1], point4[1]], color=edge_style.color,
                        dashes=[30, 5, 10, 5])
            else:
                ax.plot([point3[0], point4[0]], [point3[1], point4[1]], color=edge_style.color)

        return ax

    def plot_data(self, edge_style=None):
        """
        Get plot data for the line.

        :param edge_style: Plotting style for the line.
        :type edge_style: :class:`plot_data.EdgeStyle`, optional
        :return: Plot data for the line.
        :rtype: :class:`plot_data.Line2D`
        """
        return plot_data.Line2D([self.point1.x, self.point1.y],
                                [self.point2.x, self.point2.y],
                                edge_style=edge_style)

    def line_intersections(self, line):
        """
        Calculate the intersection between the two lines.

        :param line: The line to calculate intersections with.
        :type line: :class:`volmdlr.Line2D`
        :return: A list of at most one intersection point between
            the two lines.
        :rtype: List[:class:`volmdlr.Point2D`]
        """

        point = volmdlr.Point2D.line_intersection(self, line)
        if point is not None:
            point_projection1, _ = self.point_projection(point)
            if point_projection1 is None:
                return []

            if line.__class__.__name__ == 'Line2D':
                point_projection2, _ = line.point_projection(point)
                if point_projection2 is None:
                    return []

            return [point_projection1]
        return []

    @staticmethod
    def _compute_data_create_tangent_circle(line, point, other_line):
        """
        Static helper method to compute some data used in create_tangent_circle method.
        """
        if math.isclose(line.point_distance(point), 0, abs_tol=1e-10):
            vector_i = volmdlr.Vector2D(point.x, point.y)
            vector_a = volmdlr.Vector2D(line.point1.x, line.point1.y)
            vector_b = volmdlr.Vector2D(line.point2.x, line.point2.y)
            vector_c = volmdlr.Vector2D(other_line.point1.x, other_line.point1.y)
            vector_d = volmdlr.Vector2D(other_line.point2.x, other_line.point2.y)
        elif math.isclose(other_line.point_distance(point), 0, abs_tol=1e-10):
            vector_i = volmdlr.Vector2D(line.x, point.y)
            vector_c = volmdlr.Vector2D(line.point1.x, line.point1.y)
            vector_d = volmdlr.Vector2D(line.point2.x, line.point2.y)
            vector_a = volmdlr.Vector2D(other_line.point1.x, other_line.point1.y)
            vector_b = volmdlr.Vector2D(other_line.point2.x, other_line.point2.y)
        else:
            raise AttributeError("The point isn't on any of the two lines")
        return vector_i, vector_a, vector_b, vector_c, vector_d

    @staticmethod
    def _change_reference_frame(vector_i, vector_a, vector_b, vector_c, vector_d):
        new_u = volmdlr.Vector2D((vector_b - vector_a))
        new_u.normalize()
        new_v = new_u.unit_normal_vector()
        new_basis = volmdlr.Frame2D(vector_i, new_u, new_v)

        new_a = new_basis.global_to_local_coordinates(vector_a)
        new_b = new_basis.global_to_local_coordinates(vector_b)
        new_c = new_basis.global_to_local_coordinates(vector_c)
        new_d = new_basis.global_to_local_coordinates(vector_d)

        return new_basis, new_a, new_b, new_c, new_d

    @staticmethod
    def compute_tangent_circle_for_parallel_segments(new_basis, new_a, new_c):
        """
        Compute tangent circle betwen parallel segments.

        """
        segments_distance = abs(new_c[1] - new_a[1])
        radius = segments_distance / 2
        new_circle_center = volmdlr.Point2D((0, npy.sign(new_c[1] - new_a[1]) * radius))
        circle_center = new_basis.local_to_global_coordinates(new_circle_center)
        circle = volmdlr.wires.Circle2D(circle_center, radius)
        return circle, None

    @staticmethod
    def compute_tangent_circles_for_perpendicular_segments(new_basis, new_a, new_b, new_c, new_d):
        """
        Compute tangent circle betwen perpendicular segments.

        """
        line_ab = Line2D(volmdlr.Point2D(new_a), volmdlr.Point2D(new_b))
        line_cd = Line2D(volmdlr.Point2D(new_c), volmdlr.Point2D(new_d))
        new_pt_k = volmdlr.Point2D.line_intersection(line_ab, line_cd)

        radius = abs(new_pt_k[0])
        new_circle_center1 = volmdlr.Point2D((0, radius))
        new_circle_center2 = volmdlr.Point2D((0, -radius))
        circle_center1 = new_basis.local_to_global_coordinates(new_circle_center1)
        circle_center2 = new_basis.local_to_global_coordinates(new_circle_center2)
        circle1 = volmdlr.wires.Circle2D(circle_center1, radius)
        circle2 = volmdlr.wires.Circle2D(circle_center2, radius)

        return circle1, circle2

    def create_tangent_circle(self, point, other_line):
        """
        Computes the two circles that are tangent to 2 lines and intersect a point located on one of the two lines.
        """
        # point will be called I(x_I, y_I)
        # self will be (AB)
        # line will be (CD)
        vector_i, vector_a, vector_b, vector_c, vector_d = self._compute_data_create_tangent_circle(
            self, point, other_line)
        # Basis change
        new_basis, new_a, new_b, new_c, new_d = self._change_reference_frame(vector_i, vector_a, vector_b,
                                                                             vector_c, vector_d)

        if new_c[1] == 0 and new_d[1] == 0:
            # Segments are on the same line: no solution
            return None, None

        if math.isclose(self.unit_direction_vector().dot(
                other_line.unit_normal_vector()), 0, abs_tol=1e-06):
            # Parallel segments: one solution
            return self.compute_tangent_circle_for_parallel_segments(new_basis, new_a, new_c)

        if math.isclose(self.unit_direction_vector().dot(
                other_line.unit_direction_vector()), 0, abs_tol=1e-06):
            # Perpendicular segments: 2 solution
            return self.compute_tangent_circles_for_perpendicular_segments(new_basis, new_a, new_b, new_c, new_d)

        # =============================================================================
        # LES SEGMENTS SONT QUELCONQUES
        #   => 2 SOLUTIONS
        # =============================================================================

        line_ab = Line2D(volmdlr.Point2D(new_a), volmdlr.Point2D(new_b))
        line_cd = Line2D(volmdlr.Point2D(new_c), volmdlr.Point2D(new_d))
        new_pt_k = volmdlr.Point2D.line_intersection(line_ab, line_cd)
        pt_k = volmdlr.Point2D(new_basis.local_to_global_coordinates(new_pt_k))

        if pt_k.is_close(vector_i):
            return None, None

        # CHANGEMENT DE REPERE:
        new_u2 = volmdlr.Vector2D(pt_k - vector_i)
        new_u2.normalize()
        new_v2 = new_u2.normal_vector(unit=True)
        new_basis2 = volmdlr.Frame2D(vector_i, new_u2, new_v2)

        new_c = new_basis2.global_to_local_coordinates(vector_c)
        new_d = new_basis2.global_to_local_coordinates(vector_d)
        new_pt_k = new_basis2.global_to_local_coordinates(pt_k)

        teta1 = math.atan2(new_c[1], new_c[0] - new_pt_k[0])
        teta2 = math.atan2(new_d[1], new_d[0] - new_pt_k[0])

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

        radius1 = new_pt_k[0] * math.sin(teta) / (1 + math.cos(teta))
        radius2 = new_pt_k[0] * math.sin(teta) / (1 - math.cos(teta))

        new_circle_center1 = volmdlr.Point2D(0, -radius1)
        new_circle_center2 = volmdlr.Point2D(0, radius2)

        circle_center1 = new_basis2.local_to_global_coordinates(new_circle_center1)
        circle_center2 = new_basis2.local_to_global_coordinates(new_circle_center2)

        if new_basis.global_to_local_coordinates(circle_center1)[1] > 0:
            circle1 = volmdlr.wires.Circle2D(circle_center1, radius1)
            circle2 = volmdlr.wires.Circle2D(circle_center2, radius2)
        else:
            circle1 = volmdlr.wires.Circle2D(circle_center2, radius2)
            circle2 = volmdlr.wires.Circle2D(circle_center1, radius1)

        return circle1, circle2

    def cut_between_two_points(self, point1: volmdlr.Point2D,
                               point2: volmdlr.Point2D):
        """
        Cut the line between two points to create a linesegment.

        :param point1: The first point defining the linesegment
        :type point1: :class:`volmdlr.Point2D`
        :param point2: The second point defining the linesegment
        :type point2: :class:`volmdlr.Point2D`
        :return: The created linesegment
        :rtype: :class:`volmdlr.edges.LineSegment2D`
        """
        return LineSegment2D(point1, point2)

    def point_belongs(self, point2d, abs_tol: float = 1e-6):
        """
        Verifies if the point 2D belongs to the line.

        :param point2d: point to be verified.
        :param abs_tol: absolute tolerance to consider in calculus.
        :return: True if point belongs to line, False otherwise.
        """
        return math.isclose(self.point_distance(point2d), 0, abs_tol=abs_tol)

    def point_distance(self, point2d):
        """
        Calculate the shortest distance between a line and a point.

        :param point2d: Point to calculate distance.
        :type point2d: :class:`volmdlr.Point2D`.
        :return: Distance to point.
        :rtype: float.
        """
        vector_r = self.point1 - point2d
        vector_v = self.normal_vector()
        return abs(vector_v.dot(vector_r)) / vector_v.norm()


class BSplineCurve2D(BSplineCurve):
    """
    A class for 2 dimensional B-spline curves.

    The following rule must be
    respected : `number of knots = number of control points + degree + 1`.

    :param degree: The degree of the 2 dimensional B-spline curve
    :type degree: int
    :param control_points: A list of 2 dimensional points
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

    _non_serializable_attributes = ['curve']

    def __init__(self,
                 degree: int,
                 control_points: List[volmdlr.Point2D],
                 knot_multiplicities: List[int],
                 knots: List[float],
                 weights: List[float] = None,
                 periodic: bool = False,
                 name: str = ''):
        self._bounding_rectangle = None

        BSplineCurve.__init__(self, degree,
                              control_points,
                              knot_multiplicities,
                              knots,
                              weights,
                              periodic,
                              name)
        self._bounding_rectangle = None
        self._length = None

    @property
    def bounding_rectangle(self):
        """
        Computes the bounding rectangle of the 2 dimensional B-spline curve.

        :return: The bounding rectangle.
        :rtype: :class:`volmdlr.core.BoundingRectangle`
        """
        if not self._bounding_rectangle:
            self._bounding_rectangle = volmdlr.core.BoundingRectangle.from_points(self.points)
        return self._bounding_rectangle

    def tangent(self, position: float = 0.0):
        """
        Computes the tangent at a given parameter between 0 and 1.

        :param position: The parameter at which the tangent is computed.
        :type position: float
        :return: A 2 dimensional point representing the tangent
        :rtype: :class:`volmdlr.Point2D`
        """
        _, tangent = operations.tangent(self.curve, position,
                                        normalize=True)
        tangent = volmdlr.Point2D(tangent[0], tangent[1])
        return tangent

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
                              self.weights, self.periodic)

    def to_step(self, current_id, surface_id=None):
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
        return content, [point_id + 1]

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
                              self.weights, self.periodic)

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        BSplineCurve2D rotation. Object is updated inplace.

        :param center: rotation center
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in self.control_points:
            point.rotation_inplace(center, angle)

    def line_crossings(self, line2d: Line2D):
        polygon_points = self.discretization_points(number_points=50)
        crossings = []
        for p1, p2 in zip(polygon_points[:-1], polygon_points[1:]):
            linesegment = LineSegment2D(p1, p2)
            crossings.extend(linesegment.line_crossings(line2d))
        return crossings

    def reverse(self):
        """
        Reverse the Bspline's direction by reversing its start and end points.

        """

        return self.__class__(degree=self.degree,
                              control_points=self.control_points[::-1],
                              knot_multiplicities=self.knot_multiplicities[::-1],
                              knots=self.knots[::-1],
                              weights=self.weights,
                              periodic=self.periodic)

    def point_distance(self, point):
        """
        Calculates the distance from a given point to a BSplineCurve2D.

        :param point: point 2d.
        :return: distance.
        """
        best_distance = math.inf
        abscissa1 = 0
        abscissa2 = self.abscissa(self.end)
        distance = best_distance
        point1_ = None
        point2_ = None
        while True:
            discretized_points_between_1_2 = []
            for abscissa in npy.linspace(abscissa1, abscissa2, num=8):
                abscissa_point = self.point_at_abscissa(abscissa)
                if not volmdlr.core.point_in_list(abscissa_point, discretized_points_between_1_2):
                    discretized_points_between_1_2.append(abscissa_point)
            if not discretized_points_between_1_2:
                break
            distance = point.point_distance(discretized_points_between_1_2[0])
            for point1, point2 in zip(discretized_points_between_1_2[:-1], discretized_points_between_1_2[1:]):
                line = LineSegment2D(point1, point2)
                dist = line.point_distance(point)
                if dist < distance:
                    point1_ = point1
                    point2_ = point2
                    distance = dist
            if not point1_ or math.isclose(distance, best_distance, abs_tol=1e-6):
                break
            abscissa1 = self.abscissa(point1_)
            abscissa2 = self.abscissa(point2_)
            best_distance = distance
            if math.isclose(abscissa1, abscissa2, abs_tol=1e-6):
                break
        return distance

    def nearest_point_to(self, point):
        """
        Find out the nearest point on the linesegment to point.

        """

        points = self.discretization_points(number_points=500)
        return point.nearest_point(points)

    def linesegment_intersections(self, linesegment2d):
        """
        Calculates intersections between a BSplineCurve2D and a LineSegment2D.

        :param linesegment2d: linesegment to verify intersections.
        :return: list with the intersections points.
        """
        if not self.bounding_rectangle.b_rectangle_intersection(linesegment2d.bounding_rectangle):
            return []
        intersections_points = self.get_linesegment_intersections(linesegment2d)
        return intersections_points

    def axial_symmetry(self, line):
        """
        Finds out the symmetric bsplinecurve2d according to a line.

        """

        points_symmetry = [point.axial_symmetry(line) for point in self.control_points]

        return self.__class__(degree=self.degree,
                              control_points=points_symmetry,
                              knot_multiplicities=self.knot_multiplicities[::-1],
                              knots=self.knots[::-1],
                              weights=self.weights,
                              periodic=self.periodic)

    def offset(self, offset_length: float):
        """
        Offsets a BSplineCurve2D in one of its normal direction.

        :param offset_length: the length taken to offset the BSpline. if positive, the offset is in the normal
            direction of the curve. if negative, in the opposite direction of the normal.
        :return: returns an offset bsplinecurve2D, created with from_points_interpolation.
        """
        unit_normal_vectors = [self.unit_normal_vector(
            self.abscissa(point)) for point in self.points]
        offseted_points = [point.translation(normal_vector * offset_length) for point, normal_vector
                           in zip(self.points, unit_normal_vectors)]
        offseted_bspline = BSplineCurve2D.from_points_interpolation(offseted_points, self.degree,
                                                                    self.periodic)
        return offseted_bspline

    def point_belongs(self, point: volmdlr.Point2D, abs_tol: float = 1e-6):
        """
        Checks if a 2D point belongs to the B-spline curve 2D or not. It uses the point_distance.

        :param point: The point to be checked
        :type point: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
        :param abs_tol: The precision in terms of distance.
            Default value is 1e-7
        :type abs_tol: float, optional
        :return: `True` if the point belongs to the B-spline curve, `False`
            otherwise
        :rtype: bool
        """
        if self.point_distance(point) < abs_tol:
            return True
        return False


class BezierCurve2D(BSplineCurve2D):
    """
    A class for 2 dimensional Bezier curves.

    :param degree: The degree of the Bezier curve.
    :type degree: int
    :param control_points: A list of 2 dimensional points
    :type control_points: List[:class:`volmdlr.Point2D`]
    :param name: The name of the B-spline curve. Default value is ''
    :type name: str, optional
    """

    def __init__(self, degree: int, control_points: List[volmdlr.Point2D],
                 name: str = ''):
        knotvector = utilities.generate_knot_vector(degree,
                                                    len(control_points))
        knot_multiplicity = [1] * len(knotvector)

        BSplineCurve2D.__init__(self, degree, control_points,
                                knot_multiplicity, knotvector,
                                None, False, name)


class LineSegment2D(LineSegment):
    """
    Define a line segment limited by two points.

    """

    def __init__(self, start: volmdlr.Point2D, end: volmdlr.Point2D, *, name: str = ''):
        if start.is_close(end):
            raise NotImplementedError
        self._bounding_rectangle = None
        LineSegment.__init__(self, start, end, name=name)

    def __hash__(self):
        # return self._data_hash()
        # tolerance = 1e-6
        # factor = 1 / tolerance
        # return hash(math.floor(component * factor) / factor for point in [self.start, self.end]
        #             for component in [point.x, point.y])
        return hash(('linesegment2d', self.start, self.end))

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
        return {'object_class': 'volmdlr.edges.LineSegment2D',
                'name': self.name,
                'start': self.start.to_dict(),
                'end': self.end.to_dict()
                }

    # def middle_point(self):
    #     return 0.5 * (self.start + self.end)
    #
    # def point_at_abscissa(self, abscissa):
    #     return self.start + self.unit_direction_vector() * abscissa

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

    def straight_line_second_moment_area(self, point: volmdlr.Point2D):
        return 0, 0, 0

    def straight_line_center_of_mass(self):
        """Straight line center of mass."""
        return 0.5 * (self.start + self.end)

    def point_distance(self, point, return_other_point=False):
        """
        Computes the distance of a point to segment of line.

        :param point: point to calculate distance.
        :param return_other_points: Boolean variable to return linesegment's corresponding point or not.
        """
        distance, point = volmdlr.LineSegment2DPointDistance(
            [(self.start.x, self.start.y), (self.end.x, self.end.y)],
            (point.x, point.y))
        if return_other_point:
            return distance, volmdlr.Point2D(*point)
        return distance

    def point_projection(self, point):
        """
        If the projection falls outside the LineSegment2D, returns None.
        """
        point, curv_abs = Line2D.point_projection(Line2D(self.start, self.end),
                                                  point)
        # print('curv_abs :', curv_abs, 'length :', self.length())
        if curv_abs < 0 or curv_abs > self.length():
            if abs(curv_abs) < 1e-6 or math.isclose(curv_abs, self.length(),
                                                    abs_tol=1e-6):
                return point, curv_abs
            return None, curv_abs
        return point, curv_abs

    def line_intersections(self, line: Line2D):
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

    def linesegment_intersections(self, linesegment2d: 'LineSegment2D'):
        """
        Touching line segments does not intersect.
        """
        if not self.bounding_rectangle.b_rectangle_intersection(linesegment2d.bounding_rectangle):
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

    def line_crossings(self, line: 'Line2D'):
        if self.direction_vector().is_colinear_to(line.direction_vector()):
            return []
        line_intersection = self.line_intersections(line)
        if line_intersection and (line_intersection[0].is_close(self.end) or
                                  line_intersection[0].is_close(self.start)):
            return []
        return line_intersection

    def linesegment_crossings(self, linesegment: 'LineSegment2D'):
        """
        Gives the crossings with a linesegment.
        """
        if self.direction_vector().is_colinear_to(
                linesegment.direction_vector()):
            return []
        return self.linesegment_intersections(linesegment)

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

    def reverse(self):
        """
        Invert the sense of the line segment.
        """
        return LineSegment2D(self.end.copy(), self.start.copy())

    def to_line(self):
        return Line2D(self.start, self.end)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        LineSegment2D rotation.

        :param center: rotation center
        :param angle: angle rotation
        :return: a new rotated LineSegment2D
        """
        return LineSegment2D(self.start.rotation(center, angle),
                             self.end.rotation(center, angle))

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        LineSegment2D rotation. Object is updated inplace.

        :param center: rotation center.
        :param angle: rotation angle.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in [self.start, self.end]:
            point.rotation_inplace(center, angle)

    def translation(self, offset: volmdlr.Vector2D):
        """
        LineSegment2D translation.

        :param offset: translation vector.
        :return: A new translated LineSegment2D.
        """
        return LineSegment2D(self.start.translation(offset),
                             self.end.translation(offset))

    def translation_inplace(self, offset: volmdlr.Vector2D):
        """
        LineSegment2D translation. Object is updated inplace.

        :param offset: translation vector.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in [self.start, self.end]:
            point.translation_inplace(offset)

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

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes vector frame_mapping and the object is updated inplace.

        :param frame: frame to execute the frame mapping.
        :param side: 'old' or 'new'.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        if side == 'old':
            new_start = frame.local_to_global_coordinates(self.start)
            new_end = frame.local_to_global_coordinates(self.end)
        elif side == 'new':
            new_start = frame.global_to_local_coordinates(self.start)
            new_end = frame.global_to_local_coordinates(self.end)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        self.start = new_start
        self.end = new_end

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
        circle1, circle2 = Line2D.create_tangent_circle(other_line, point, self)
        if circle1 is not None:
            _, curv_abs1 = Line2D.point_projection(self, circle1.center)
            if curv_abs1 < 0. or curv_abs1 > self.length():
                circle1 = None
        if circle2 is not None:
            _, curv_abs2 = Line2D.point_projection(self, circle2.center)
            if curv_abs2 < 0. or curv_abs2 > self.length():
                circle2 = None
        return circle1, circle2

    def infinite_primitive(self, offset):
        n = -self.unit_normal_vector()
        offset_point_1 = self.start + offset * n
        offset_point_2 = self.end + offset * n

        return Line2D(offset_point_1, offset_point_2)

    def to_wire(self, n: int):
        """
        Convert a linesegment2d to a wire 2D defined with 'n' line_segments.

        """
        warnings.warn('To avoid Circular imports, a new method was created in Wire2D called from_edge.'
                      'You can use it instead of to_wire.')
        raise AttributeError

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


class Arc(Edge):
    """
    Abstract class representing an arc.

    :param start: The starting point
    :type start: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
    :param end: The finish point
    :type end: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
    :param interior: An interior point
    :type interior: Union[:class:`volmdlr.Point2D`, :class:`volmdlr.Point3D`]
    :param name: The name of the arc. Default value is an empty string
    :type name: str, optional
    """

    def __init__(self, start,
                 end,
                 interior,
                 name: str = ''):
        Edge.__init__(self, start=start, end=end, name=name)
        self.interior = interior
        self._utd_clockwise_and_trigowise_paths = False
        self._clockwise_and_trigowise_paths = None
        self._radius = None
        self._length = None

    @property
    def center(self):
        """
        Gets the arc's center.

        :return: The center of the arc.
        """
        raise NotImplementedError(
            'the property method center must be overloaded by subclassing'
            'class if not a given parameter')

    @property
    def angle(self):
        """
        Gets the angle of the arc.

        :return: The angle of the arc.
        """
        return NotImplementedError(
            'the property method angle must be overloaded by subclassing'
            'class if not a given parameter')

    @property
    def is_trigo(self):
        """
        Verifies if arc is trigonometric wise or clockwise.

        :return: True if trigonometric wise or False otherwise.
        """
        return NotImplementedError(
            'the property method is_trigo must be overloaded by subclassing'
            'class if not a given parameter')

    @property
    def radius(self):
        if not self._radius:
            self._radius = (self.start - self.center).norm()
        return self._radius

    def length(self):
        """
        Calculates the length of the Arc, with its radius, and its arc angle.

        :return: the length of the Arc.
        """
        if not self._length:
            self._length = self.radius * abs(self.angle)
        return self._length

    def point_at_abscissa(self, abscissa):
        """
        Calculates a point in the Arc at a given abscissa.

        :param abscissa: abscissa where in the curve the point should be calculated.
        :return: Corresponding point.
        """
        if self.is_trigo:
            return self.start.rotation(self.center, abscissa / self.radius)
        return self.start.rotation(self.center, -abscissa / self.radius)

    def normal_vector(self, abscissa: float):
        """
        Get the normal vector of the Arc2D.

        :param abscissa: defines where in the Arc2D the
        normal vector is to be calculated
        :return: The normal vector of the Arc2D
        """
        point = self.point_at_abscissa(abscissa)
        normal_vector = self.center - point
        return normal_vector

    def direction_vector(self, abscissa: float):
        """
        Get direction vector of the Arc2D.

        :param abscissa: defines where in the Arc2D the
        direction vector is to be calculated
        :return: The direction vector of the Arc2D
        """
        return -self.normal_vector(abscissa=abscissa).normal_vector()

    @staticmethod
    def get_clockwise_and_trigowise_paths(radius_1, radius_2, radius_i):
        """
        Get arc paths from radius.

        :param radius_1: radius from center to start point.
        :param radius_2: radius form center to end point.
        :param radius_i: radius from center to interior point.
        :return: the clockwise and trigowise paths.
        """
        angle1 = math.atan2(radius_1.y, radius_1.x)
        anglei = math.atan2(radius_i.y, radius_i.x)
        angle2 = math.atan2(radius_2.y, radius_2.x)

        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei + volmdlr.TWO_PI) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + volmdlr.TWO_PI

        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2 + volmdlr.TWO_PI) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + volmdlr.TWO_PI
        return clockwise_path, trigowise_path

    def middle_point(self):
        """
        Get point in the middle of Arc.
        """
        return self.point_at_abscissa(0.5 * self.length())

    def point_distance(self, point):
        """Returns the minimal distance to a point."""
        points = self.discretization_points(angle_resolution=100)
        return point.point_distance(point.nearest_point(points))

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = None):
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
                number_points = math.ceil(self.angle * angle_resolution) + 2

        step = self.length() / (number_points - 1)
        return [self.point_at_abscissa(i * step)
                for i in range(number_points)]

    def polygon_points(self, discretization_resolution: int):
        warnings.warn('polygon_points is deprecated,\
        please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(number_points=discretization_resolution)

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
        return [self.start, self.center, self.end]

    def reverse(self):
        """
        Gets the reverse version of an arc.

        :return: An arc
        :rtype: Arc
        """

        return self.__class__(start=self.end, interior=self.interior, end=self.start)

    def split(self, split_point):
        """
        Splits arc at a given point.

        :param split_point: splitting point.
        :return: list of two Arc.
        """
        if split_point.is_close(self.start, 1e-6):
            return [None, self.copy()]
        if split_point.is_close(self.end, 1e-6):
            return [self.copy(), None]
        abscissa = self.abscissa(split_point)
        return [self.__class__(self.start, self.point_at_abscissa(0.5 * abscissa), split_point),
                self.__class__(split_point, self.point_at_abscissa((self.abscissa(self.end) -
                                                                    abscissa) * 0.5 + abscissa), self.end)]

    def get_shared_section(self, other_arc2):
        """
        Gets the shared section between two arcs.

        :param arc2: other arc to verify for shared section.
        :return: shared arc section.
        """
        if self.__class__ != other_arc2.__class__:
            return []
        if not self.center.is_close(other_arc2.center) or self.radius != self.radius or \
                not any(self.point_belongs(point) for point in [other_arc2.start,
                                                                other_arc2.interior, other_arc2.end]):
            return []
        if all(self.point_belongs(point) for point in [other_arc2.start, other_arc2.interior, other_arc2.end]):
            return [other_arc2]
        if all(other_arc2.point_belongs(point) for point in [self.start, self.interior, self.end]):
            return [self]
        if self.point_belongs(other_arc2.start):
            arc1_, arc2_ = self.split(other_arc2.start)
        elif self.point_belongs(other_arc2.end):
            arc1_, arc2_ = self.split(other_arc2.end)
        else:
            raise NotImplementedError
        shared_arc_section = []
        for arc in [arc1_, arc2_]:
            if arc and all(other_arc2.point_belongs(point) for point in [arc.start, arc.interior, arc.end]):
                shared_arc_section.append(arc)
                break
        return shared_arc_section

    def delete_shared_section(self, other_arc2):
        """
        Deletes from self, the section shared with the other arc.

        :param other_arc2:
        :return:
        """
        shared_section = self.get_shared_section(other_arc2)
        if not shared_section:
            return [self]
        if shared_section == self:
            return []
        split_arcs1 = self.split(shared_section[0].start)
        split_arcs2 = self.split(shared_section[0].end)
        new_arcs = []
        for arc in split_arcs1 + split_arcs2:
            if arc and not arc.point_belongs(shared_section[0].interior):
                new_arcs.append(arc)
        return new_arcs


class Arc2D(Arc):
    """
    Class to draw Arc2D.

    angle: the angle measure always >= 0
    """

    def __init__(self,
                 start: volmdlr.Point2D,
                 interior: volmdlr.Point2D,
                 end: volmdlr.Point2D,
                 name: str = ''):
        self._center = None
        self._is_trigo = None
        self._angle = None
        self._bounding_rectangle = None
        Arc.__init__(self, start=start, end=end, interior=interior, name=name)
        start_to_center = start - self.center
        end_to_center = end - self.center
        angle1 = math.atan2(start_to_center.y, start_to_center.x)
        angle2 = math.atan2(end_to_center.y, end_to_center.x)
        if self.is_trigo:
            self.angle1 = angle1
            self.angle2 = angle2
        else:
            self.angle1 = angle2
            self.angle2 = angle1

    def __hash__(self):
        return hash(('arc2d', self.start, self.interior, self.end))

    def __eq__(self, other_arc):
        if self.__class__.__name__ != other_arc.__class__.__name__:
            return False
        return (self.center == other_arc.center
                and self.start == other_arc.start
                and self.end == other_arc.end
                and self.interior == other_arc.interior)

    @property
    def center(self):
        if not self._center:
            self._center = self.get_center()
        return self._center

    def get_center(self):
        """
        Calculates the center of the Arc.

        :return: asc's center.
        """
        x_interior, y_interior = self.interior.x, self.interior.y
        x_end, y_end = self.end.x, self.end.y
        x_start, y_start = self.start.x, self.start.y
        try:
            matrix_a = volmdlr.Matrix22(2 * (x_start - x_interior), 2 * (y_start - y_interior),
                                        2 * (x_start - x_end), 2 * (y_start - y_end))
            b_vector = - volmdlr.Vector2D(x_interior ** 2 + y_interior ** 2 - x_start ** 2 - y_start ** 2,
                                          x_end ** 2 + y_end ** 2 - x_start ** 2 - y_start ** 2)
            inv_matrix_a = matrix_a.inverse()
            x = inv_matrix_a.vector_multiplication(b_vector)
            center = volmdlr.Point2D(x.x, x.y)
        except ValueError:
            matrix_a = npy.array([[2 * (x_start - x_interior), 2 * (y_start - y_interior)],
                                  [2 * (x_start - x_end), 2 * (y_start - y_end)]])
            b_vector = - npy.array([x_interior ** 2 + y_interior ** 2 - x_start ** 2 - y_start ** 2,
                                    x_end ** 2 + y_end ** 2 - x_start ** 2 - y_start ** 2])
            center = volmdlr.Point2D(*npy.linalg.solve(matrix_a, b_vector))
        return center

    @property
    def is_trigo(self):
        """
        Gives if the edge goes in the trigo direction.
        """
        if not self._is_trigo:
            self._is_trigo = self.get_arc_direction()
        return self._is_trigo

    @property
    def clockwise_and_trigowise_paths(self):
        if not self._clockwise_and_trigowise_paths:
            radius_1 = self.start - self.center
            radius_2 = self.end - self.center
            radius_i = self.interior - self.center
            self._clockwise_and_trigowise_paths = \
                self.get_clockwise_and_trigowise_paths(radius_1,
                                                       radius_2,
                                                       radius_i)
            self._utd_clockwise_and_trigowise_paths = True
        return self._clockwise_and_trigowise_paths

    def get_arc_direction(self):
        """
        Gets arc direction: clockwise or trigonometric.

        :return: True if clockwise. False if counterclockwise.
        """
        clockwise_path, trigowise_path = \
            self.clockwise_and_trigowise_paths
        if clockwise_path > trigowise_path:
            return True
        return False

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
        clockwise_path, trigowise_path = \
            self.clockwise_and_trigowise_paths
        if self.is_trigo:
            return trigowise_path
        return clockwise_path

    def _get_points(self):
        return [self.start, self.interior, self.end]

    points = property(_get_points)

    def point_distance(self, point):
        """
        Returns the distance between a point and the edge.
        """
        vector_start = self.start - self.center
        vector_point = point - self.center
        vector_end = self.end - self.center
        if self.is_trigo:
            vector_start, vector_end = vector_end, vector_start
        arc_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_end)
        point_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_point)
        if point_angle <= arc_angle:
            return abs(
                LineSegment2D(point, self.center).length() - self.radius)
        return min(point.point_distance(self.start), point.point_distance(self.end))

    def point_belongs(self, point2d, abs_tol=1e-6):
        """
        Check if a Point2D belongs to the Arc2D.

        """
        distance_point_to_center = point2d.point_distance(self.center)
        if not math.isclose(distance_point_to_center, self.radius, abs_tol=abs_tol):
            return False
        if point2d.is_close(self.start) or point2d.is_close(self.end):
            return True
        clockwise_arc = self.reverse() if self.is_trigo else self
        vector_start = clockwise_arc.start - clockwise_arc.center
        vector_end = clockwise_arc.end - clockwise_arc.center
        vector_point = point2d - clockwise_arc.center
        arc_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_end)
        point_start_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_point)
        point_end_angle = volmdlr.geometry.clockwise_angle(vector_point, vector_end)
        if math.isclose(arc_angle, point_start_angle + point_end_angle, abs_tol=abs_tol):
            return True
        return False

    def to_full_arc_2d(self):
        """
        Convert to a full arc2d.
        """
        return FullArc2D(center=self.center,
                         start_end=self.point_at_abscissa(0),
                         name=self.name)

    def line_intersections(self, line2d: Line2D):
        """
        Calculates the intersection between a line and an Arc2D.

        :param line2d: Line2D to verify intersections.
        :return: a list with intersections points.
        """
        # circle = self.to_circle()
        # circle_intersection_points = circle.line_intersections(line2d)
        full_arc_2d = self.to_full_arc_2d()
        fa2d_intersection_points = full_arc_2d.line_intersections(line2d)
        intersection_points = []
        for pt in fa2d_intersection_points:
            if self.point_belongs(pt):
                intersection_points.append(pt)
        return intersection_points

    def linesegment_intersections(self, linesegment2d: LineSegment2D):
        """
        Calculates the intersection between a LineSegment2D and an Arc2D.

        :param line2d: LineSegment2D to verify intersections.
        :return: a list with intersections points.
        """
        if not self.bounding_rectangle.b_rectangle_intersection(linesegment2d.bounding_rectangle):
            return []
        full_arc_2d = self.to_full_arc_2d()
        fa2d_intersection_points = full_arc_2d.linesegment_intersections(
            linesegment2d)
        intersection_points = []
        for point in fa2d_intersection_points:
            if self.point_belongs(point):
                intersection_points.append(point)
        return intersection_points

    def abscissa(self, point: volmdlr.Point2D, tol=1e-6):
        """
        Returns the abscissa of a given point 2d.
        """
        if not math.isclose(point.point_distance(self.center), self.radius, abs_tol=tol):
            raise ValueError('Point not in arc')
        if point.point_distance(self.start) < tol:
            return 0
        if point.point_distance(self.end) < tol:
            return self.length()
        clockwise_arc = self.reverse() if self.is_trigo else self
        vector_start = clockwise_arc.start - clockwise_arc.center
        vector_end = clockwise_arc.end - clockwise_arc.center
        vector_point = point - clockwise_arc.center
        arc_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_end)
        point_start_angle = volmdlr.geometry.clockwise_angle(vector_start, vector_point)
        point_end_angle = volmdlr.geometry.clockwise_angle(vector_point, vector_end)
        if math.isclose(arc_angle, point_start_angle + point_end_angle, abs_tol=tol):
            if self.is_trigo:
                return self.length() - self.radius * point_start_angle
            return self.radius * point_start_angle
        raise ValueError('Point not in arc')

    def area(self):
        """
        Calculates the area of the Arc2D.

        :return: the area of the Arc2D.
        """
        return self.radius ** 2 * self.angle / 2

    def center_of_mass(self):
        """
        Calculates the center of mass of the Arc2D.

        :return: center of mass point.
        """
        #        u=self.middle.vector-self.center.vector
        u = self.middle_point() - self.center
        u.normalize()
        # alpha = abs(self.angle)
        return self.center + 4 / (3 * self.angle) * self.radius * math.sin(
            self.angle * 0.5) * u

    @property
    def bounding_rectangle(self):
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
            area = math.pi * self.radius ** 2 - 0.5 * self.radius ** 2 * (
                    angle - math.sin(angle))
        else:
            angle = self.angle
            area = 0.5 * self.radius ** 2 * (angle - math.sin(angle))

        if self.is_trigo:
            return area
        return -area

    def straight_line_second_moment_area(self, point: volmdlr.Point2D):

        if self.angle2 < self.angle1:
            angle2 = self.angle2 + volmdlr.TWO_PI

        else:
            angle2 = self.angle2
        angle1 = self.angle1

        # Full arc section
        moment_area_x1 = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle1) - math.sin(2 * angle2)))
        moment_area_y1 = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle2) - math.sin(2 * angle1)))
        moment_area_xy1 = self.radius ** 4 / 8 * (
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
                                          self.center,
                                          point)

    def _full_arc_moment_inertia(self, angle1, angle2):
        moment_inertia_x1 = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle1) - math.sin(2 * angle2)))
        moment_inertia_y1 = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle2) - math.sin(2 * angle1)))
        moment_inertia_xy1 = self.radius ** 4 / 8 * (
                math.cos(angle1) ** 2 - math.cos(angle2) ** 2)
        return moment_inertia_x1, moment_inertia_y1, moment_inertia_xy1

    def _triangle_moment_inertia(self):
        xi, yi = self.start - self.center
        xj, yj = self.end - self.center
        moment_inertia_x2 = (yi ** 2 + yi * yj + yj ** 2) * (xi * yj - xj * yi) / 12.
        moment_inertia_y2 = (xi ** 2 + xi * xj + xj ** 2) * (xi * yj - xj * yi) / 12.
        moment_inertia_xy2 = (xi * yj + 2 * xi * yi + 2 * xj * yj + xj * yi) * (
                xi * yj - xj * yi) / 24.
        return moment_inertia_x2, moment_inertia_y2, moment_inertia_xy2

    def straight_line_center_of_mass(self):
        """Straight line center of mass."""
        if self.angle == math.pi:
            return self.center_of_mass()

        u = self.middle_point() - self.center
        u.normalize()
        if self.angle >= math.pi:
            u = -u
        bissec = Line2D(self.center, self.center + u)
        string = Line2D(self.start, self.end)
        point = volmdlr.Point2D.line_intersection(bissec, string)
        a = point.point_distance(self.start)
        height = point.point_distance(self.center)
        triangle_area = height * a
        # alpha = abs(self.angle)
        triangle_cog = self.center + 2 / 3. * height * u
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

        :param point_2d: Point to be verified.
        :return: Return True if the point belongs to this surface, or False otherwise.
        """
        if self.point_belongs(point):
            return True
        if self.start == self.end:
            if point.point_distance(self.center) <= self.radius:
                return True
        center_distance_point = self.center.point_distance(point)
        straight_line = LineSegment2D(self.start, self.end)
        for edge in [self, straight_line]:
            line_passing_trough_point = Line2D(self.center, point)
            straight_line_intersections = edge.line_intersections(line_passing_trough_point)
            if straight_line_intersections:
                if self.center.point_distance(straight_line_intersections[0]) > center_distance_point:
                    return True
        return False

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            _, ax = plt.subplots()

        if edge_style.plot_points:
            for point in [self.center, self.start, self.interior, self.end]:
                point.plot(ax=ax, color=edge_style.color, alpha=edge_style.alpha)

        ax.add_patch(matplotlib.patches.Arc((self.center.x, self.center.y), 2 * self.radius,
                                            2 * self.radius, angle=0,
                                            theta1=self.angle1 * 0.5 / math.pi * 360,
                                            theta2=self.angle2 * 0.5 / math.pi * 360,
                                            color=edge_style.color,
                                            alpha=edge_style.alpha))
        return ax

    def to_3d(self, plane_origin, x, y):
        point_start = self.start.to_3d(plane_origin, x, y)
        point_interior = self.interior.to_3d(plane_origin, x, y)
        point_end = self.end.to_3d(plane_origin, x, y)

        return Arc3D(point_start, point_interior, point_end, name=self.name)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Arc2D rotation.

        :param center: rotation center
        :param angle: angle rotation.
        :return: a new rotated Arc2D.
        """
        return Arc2D(*[point.rotation(center, angle, ) for point in
                       [self.start, self.interior, self.end]])

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        Arc2D rotation. Object is updated inplace.

        :param center: rotation center.
        :param angle: rotation angle.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.start.rotation_inplace(center, angle)
        self.interior.rotation_inplace(center, angle)
        self.end.rotation_inplace(center, angle)
        self._angle = None
        self._is_trigo = None
        self._center = None
        self._clockwise_and_trigowise_paths = None

    def translation(self, offset: volmdlr.Vector2D):
        """
        Arc2D translation.

        :param offset: translation vector.
        :return: A new translated Arc2D.
        """
        return Arc2D(*[point.translation(offset) for point in
                       [self.start, self.interior, self.end]])

    def translation_inplace(self, offset: volmdlr.Vector2D):
        """
        Arc2D translation. Object is updated inplace.

        :param offset: translation vector.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.start.translation_inplace(offset)
        self.interior.translation_inplace(offset)
        self.end.translation_inplace(offset)
        self._angle = None
        self._is_trigo = None
        self._center = None
        self._clockwise_and_trigowise_paths = None

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes vector frame_mapping and return a new Arc2D.

        side = 'old' or 'new'
        """
        return Arc2D(*[point.frame_mapping(frame, side) for point in
                       [self.start, self.interior, self.end]])

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes vector frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.__init__(*[point.frame_mapping(frame, side) for point in
                        [self.start, self.interior, self.end]])

    def second_moment_area(self, point):
        """
        Second moment area of part of disk.

        """
        if self.angle2 < self.angle1:
            angle2 = self.angle2 + volmdlr.TWO_PI

        else:
            angle2 = self.angle2
        angle1 = self.angle1
        moment_area_x = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle1) - math.sin(2 * angle2)))
        moment_area_y = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                math.sin(2 * angle2) - math.sin(2 * angle1)))
        moment_area_xy = self.radius ** 4 / 8 * (
                math.cos(angle1) ** 2 - math.cos(angle2) ** 2)

        # Must be computed at center, so huygens related to center
        return volmdlr.geometry.huygens2d(moment_area_x, moment_area_y, moment_area_xy, self.area(),
                                          self.center, point)

    def plot_data(self, edge_style: plot_data.EdgeStyle = None,
                  anticlockwise: bool = None):
        """
        Plot data method for a Arc2D.

        :param edge_style: edge style.
        :return: plot_data.Arc2D object.
        """
        list_node = self.discretization_points(number_points=20)
        data = []
        for node in list_node:
            data.append({'x': node.x, 'y': node.y})
        return plot_data.Arc2D(cx=self.center.x,
                               cy=self.center.y,
                               r=self.radius,
                               start_angle=self.angle1,
                               end_angle=self.angle2,
                               edge_style=edge_style,
                               data=data,
                               anticlockwise=anticlockwise,
                               name=self.name)

    def copy(self, *args, **kwargs):
        return Arc2D(self.start.copy(),
                     self.interior.copy(),
                     self.end.copy())

    def cut_between_two_points(self, point1, point2):
        """
        Cuts Arc between two points, and return a new arc between these two points.
        """
        if (point1.is_close(self.start) and point2.is_close(self.end)) or \
                (point2.is_close(self.start) and point1.is_close(self.end)):
            return self
        raise NotImplementedError

    def infinite_primitive(self, offset):
        vector_start_center = self.start - self.center
        vector_start_center.normalize()
        vector_end_center = self.end - self.center
        vector_end_center.normalize()
        vector_interior_center = self.interior - self.center
        vector_interior_center.normalize()
        if self.is_trigo:
            radius = self.radius + offset
            center = self.center

        else:
            radius = self.radius - offset
            if radius < 0:
                return None
            center = self.center
            # mid_point = self.middle_point()
            # vec1 = self.center - mid_point
            # vec1.normalize()
            # vec1 = 2 * offset * math.sqrt(2) * vec1
            # center = self.center.translation(vec1)
        start = center + radius * vector_start_center
        end = center + radius * vector_end_center
        interior = center + radius * vector_interior_center
        return Arc2D(start, interior, end)

    def complementary(self):

        interior = self.middle_point().rotation(self.center, math.pi)
        return Arc2D(self.start, interior, self.end)

    def axial_symmetry(self, line):
        """ Finds out the symmetric arc 2D according to a line. """
        points_symmetry = [point.axial_symmetry(line) for point in [self.start, self.interior, self.end]]

        return self.__class__(start=points_symmetry[0],
                              interior=points_symmetry[1],
                              end=points_symmetry[2])


class FullArc2D(Arc2D):
    """ An edge that starts at start_end, ends at the same point after having described a circle. """

    def __init__(self, center: volmdlr.Point2D, start_end: volmdlr.Point2D,
                 name: str = ''):
        self.__center = center
        self.start_end = start_end
        interior = start_end.rotation(center, math.pi)
        Arc2D.__init__(self, start=start_end, interior=interior,
                       end=start_end, name=name)  # !!! this is dangerous

    @property
    def is_trigo(self):
        return True

    @property
    def center(self):
        return self.__center

    @property
    def angle(self):
        return volmdlr.TWO_PI

    def to_dict(self, use_pointers: bool = False, memo=None, path: str = '#'):
        dict_ = self.base_dict()
        dict_['center'] = self.center.to_dict(use_pointers=use_pointers, memo=memo, path=path + '/center')
        dict_['radius'] = self.radius
        dict_['angle'] = self.angle
        dict_['is_trigo'] = self.is_trigo
        dict_['start_end'] = self.start.to_dict(use_pointers=use_pointers, memo=memo, path=path + '/start_end')
        return dict_

    def copy(self, *args, **kwargs):
        return FullArc2D(self.center.copy(), self.start.copy())

    @classmethod
    def dict_to_object(cls, dict_, global_dict=None, pointers_memo: Dict[str, Any] = None, path: str = '#'):
        center = volmdlr.Point2D.dict_to_object(dict_['center'])
        start_end = volmdlr.Point2D.dict_to_object(dict_['start_end'])

        return cls(center, start_end, name=dict_['name'])

    # def __hash__(self):
    #     return hash(("fullarc", self.center, self.radius, self.start, self.end))
    #     # return hash(self.center) + 5*hash(self.start)

    def __hash__(self):
        return hash((self.__class__.__name__, self.center, self.radius, self.start_end))

    def __eq__(self, other_arc):
        if self.__class__.__name__ != other_arc.__class__.__name__:
            return False
        return (self.center == other_arc.center) \
            and (self.start_end == other_arc.start_end)

    @property
    def bounding_rectangle(self):
        if not self._bounding_rectangle:
            self._bounding_rectangle = volmdlr.core.BoundingRectangle(
                self.center.x - self.radius, self.center.x + self.radius,
                self.center.y - self.radius, self.center.y + self.radius)
        return self._bounding_rectangle

    def straight_line_area(self):
        """
        Calculates the area of the fullarc, with line drawn from start to end.

        :return: straight_line_area.
        """
        area = self.area()
        return area

    def center_of_mass(self):
        return self.center

    def straight_line_center_of_mass(self):
        """Straight line center of mass."""
        return self.center_of_mass()

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point2d: Point to be verified.
        :return: Return True if the point belongs to this surface, or False otherwise.
        """
        if point.point_distance(self.center) <= self.radius:
            return True
        return False

    def to_3d(self, plane_origin, x, y):
        center = self.center.to_3d(plane_origin, x, y)
        start = self.start.to_3d(plane_origin, x, y)
        z = x.cross(y)
        z.normalize()

        return FullArc3D(center, start, z)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        new_center = self._center.rotation(center, angle, True)
        new_start_end = self.start.rotation(center, angle, True)
        return FullArc2D(new_center, new_start_end)

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self._center.rotation(center, angle, False)
        self.start.rotation(center, angle, False)
        self.interior.rotation(center, angle, False)
        self.end.rotation(center, angle, False)

    def translation(self, offset: volmdlr.Vector2D):
        new_center = self._center.translation(offset)
        new_start_end = self.start.translation(offset)
        return FullArc2D(new_center, new_start_end)

    def translation_inplace(self, offset: volmdlr.Vector2D):
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self._center.translation_inplace(offset)
        self.start.translation_inplace(offset)
        self.end.translation_inplace(offset)
        self.interior.translation_inplace(offset)

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
                           [self._center, self.start]])

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in [self._center, self.start, self.end, self.interior]:
            point.frame_mapping_inplace(frame, side)

    def polygonization(self):
        return volmdlr.wires.ClosedPolygon2D(self.discretization_points(angle_resolution=15))

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        return vm_common_operations.plot_circle(self, ax, edge_style)

    def cut_between_two_points(self, point1, point2):

        x1, y1 = point1 - self.center
        x2, y2 = point2 - self.center

        angle1 = math.atan2(y1, x1)
        angle2 = math.atan2(y2, x2)
        if angle2 < angle1:
            angle2 += volmdlr.TWO_PI
        angle_i = 0.5 * (angle1 + angle2)
        interior = point1.rotation(self.center, angle_i)
        arc = Arc2D(point1, interior, point2)
        if self.is_trigo != arc.is_trigo:
            arc = arc.complementary()

        return arc

    def line_intersections(self, line2d: Line2D, tol=1e-9):
        try:
            if line2d.start.is_close(self.center):
                pt1 = line2d.end
                vec = line2d.start - line2d.end
            else:
                pt1 = line2d.start
                vec = line2d.end - line2d.start
        except AttributeError:
            if line2d.point1.is_close(self.center):
                pt1 = line2d.point2
                vec = line2d.point1 - line2d.point2
            else:
                pt1 = line2d.point1
                vec = line2d.point2 - line2d.point1
        vector1 = vec.dot(vec)
        vector2 = 2 * vec.dot(pt1 - self.center)
        vector3 = pt1.dot(pt1) + self.center.dot(self.center) \
            - 2 * pt1.dot(self.center) - self.radius ** 2

        disc = vector2 ** 2 - 4 * vector1 * vector3
        if math.isclose(disc, 0., abs_tol=tol):
            t_param = -vector2 / (2 * vector1)
            return [pt1 + t_param * vec]

        if disc > 0:
            sqrt_disc = math.sqrt(disc)
            t_param = (-vector2 + sqrt_disc) / (2 * vector1)
            s_param = (-vector2 - sqrt_disc) / (2 * vector1)
            return [pt1 + t_param * vec,
                    pt1 + s_param * vec]

        return []

    def linesegment_intersections(self, linesegment2d: LineSegment2D, tol=1e-9):
        if not self.bounding_rectangle.b_rectangle_intersection(linesegment2d.bounding_rectangle):
            return []
        try:
            if linesegment2d.start.is_close(self.center):
                pt1 = linesegment2d.end
                vec = linesegment2d.start - linesegment2d.end
            else:
                pt1 = linesegment2d.start
                vec = linesegment2d.end - linesegment2d.start
        except AttributeError:
            if linesegment2d.point1.is_close(self.center):
                pt1 = linesegment2d.point2
                vec = linesegment2d.point1 - linesegment2d.point2
            else:
                pt1 = linesegment2d.point1
                vec = linesegment2d.point2 - linesegment2d.point1
        vector1 = vec.dot(vec)
        vector2 = 2 * vec.dot(pt1 - self.center)
        vector3 = pt1.dot(pt1) + self.center.dot(self.center) \
            - 2 * pt1.dot(self.center) - self.radius ** 2

        disc = vector2 ** 2 - 4 * vector1 * vector3
        if math.isclose(disc, 0., abs_tol=tol):
            t_param = -vector2 / (2 * vector1)
            points = [pt1 + t_param * vec]
            if linesegment2d.point_belongs(points[0]):
                return points
            return []

        if disc > 0:
            sqrt_disc = math.sqrt(disc)
            t_param = (-vector2 + sqrt_disc) / (2 * vector1)
            s_param = (-vector2 - sqrt_disc) / (2 * vector1)
            points = [pt1 + t_param * vec, pt1 + s_param * vec]
            valid_points = [pt for pt in points if
                            linesegment2d.point_belongs(pt)]
            return valid_points

        return []

    def reverse(self):
        return self


class ArcEllipse2D(Edge):
    """
    An 2 dimensional elliptical arc.

    :param start: The starting point of the elliptical arc
    :type start: :class:`volmdlr.Point2D`
    :param interior: An interior point of the elliptical arc
    :type interior: :class:`volmdlr.Point2D`
    :param end: The end point of the elliptical arc
    :type end: :class:`volmdlr.Point2D`
    :param center: The center of the ellipse
    :type center: :class:`volmdlr.Point2D`
    :param major_dir: The major direction of the ellipse
    :type major_dir: :class:`volmdlr.Vector2D`
    :param name: The name of the elliptical arc. Default value is ''
    :type name: str, optional
    :param extra: An extra interior point if start is equal to end. Default
        value is None
    :type extra: :class:`volmdlr.Point2D`, optional
    """

    def __init__(self, start: volmdlr.Point2D, interior: volmdlr.Point2D,
                 end: volmdlr.Point2D, center: volmdlr.Point2D,
                 major_dir: volmdlr.Vector2D, extra: volmdlr.Point2D = None, name: str = '',
                 ):
        Edge.__init__(self, start, end, name)
        self.interior = interior
        self.center = center
        self.extra = extra
        self.major_dir = major_dir
        self.minor_dir = self.major_dir.deterministic_unit_normal_vector()
        frame = volmdlr.Frame2D(self.center, self.major_dir, self.minor_dir)
        self.frame = frame
        start_new = frame.global_to_local_coordinates(self.start)
        end_new = frame.global_to_local_coordinates(self.end)
        interior_new = frame.global_to_local_coordinates(self.interior)
        center_new = frame.global_to_local_coordinates(self.center)
        self._bounding_rectangle = None

        def theta_a_b(start_, iterior_, end_, center_):
            """
            Calculates the major, minor and theta.

            From : https://math.stackexchange.com/questions/339126/how-to-draw-an-ellipse-if-a- \
            center-and-3-arbitrary-points-on-it-are-given.
            theta= ellipse's inclination angle related to the horizontal
            (clockwise), A=semi major axis, B=semi minor axis.

            """
            x_start, y_start, x_interior, y_interior, x_end, y_end = start_[0] - center_[0], start_[1] - center_[1], \
                iterior_[0] - center_[0], iterior_[1] - center_[
                                                                         1], end_[0] - center_[0], end_[1] - center_[1]
            matrix_a = npy.array(([x_start ** 2, y_start ** 2, 2 * x_start * y_start],
                                  [x_interior ** 2, y_interior ** 2, 2 * x_interior * y_interior],
                                  [x_end ** 2, y_end ** 2, 2 * x_end * y_end]))
            inv_matrix_a = npy.linalg.inv(matrix_a)
            matriz_one = npy.array(([1],
                                    [1],
                                    [1]))
            vector_c = npy.dot(inv_matrix_a, matriz_one)
            theta = 0.5 * math.atan(2 * vector_c[2] / (vector_c[1] - vector_c[0]))
            c1 = vector_c[0] + vector_c[1]
            c2 = (vector_c[1] - vector_c[0]) / math.cos(2 * theta)
            gdaxe = math.sqrt((2 / (c1 - c2)))
            ptax = math.sqrt((2 / (c1 + c2)))
            return theta, gdaxe, ptax

        if start.is_close(end):
            extra_new = frame.global_to_local_coordinates(self.extra)
            theta, major_axis, minor_axis = theta_a_b(start_new, extra_new, interior_new,
                                                      center_new)
        else:
            theta, major_axis, minor_axis = theta_a_b(start_new, interior_new, end_new,
                                                      center_new)

        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.theta = theta

        # Angle pour start
        u1, u2 = start_new.x / self.major_axis, start_new.y / self.minor_axis
        angle1 = volmdlr.geometry.sin_cos_angle(u1, u2)
        self.angle_start = angle1
        # Angle pour end
        u3, u4 = end_new.x / self.major_axis, end_new.y / self.minor_axis
        angle2 = volmdlr.geometry.sin_cos_angle(u3, u4)
        self.angle_end = angle2
        # Angle pour interior
        u5, u6 = interior_new.x / self.major_axis, interior_new.y / self.minor_axis
        anglei = volmdlr.geometry.sin_cos_angle(u5, u6)
        self.angle_interior = anglei
        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei + volmdlr.TWO_PI) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + volmdlr.TWO_PI

        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2 + volmdlr.TWO_PI) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + volmdlr.TWO_PI

        if clockwise_path > trigowise_path:
            self.is_trigo = True
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle = clockwise_path

        if self.start.is_close(self.end) or self.angle == 0:
            self.angle = volmdlr.TWO_PI

        if self.is_trigo:  # sens trigo
            self.offset_angle = angle1
        else:
            self.offset_angle = angle2

    def _get_points(self):
        return self.discretization_points(number_points=20)

    points = property(_get_points)

    def length(self):
        """
        Calculates the length of the arcellipse 2d.

        :return: arc ellipse 2d's length
        """
        length = self.abscissa(self.end)
        return length

    def point_belongs(self, point, abs_tol: float = 1e-6):
        """
        Verifies if a point belongs to the arc ellipse 2d.

        :param point: point to be verified
        :param abs_tol: tolerance applied during calculations
        :return: True if the point belongs, False otherwise
        """
        if not math.isclose((point.x - self.center.x) ** 2 / self.major_axis ** 2 +
                            (point.y - self.center.y) ** 2 / self.minor_axis ** 2, 1, abs_tol=abs_tol) and not \
                math.isclose((point.x - self.center.x) ** 2 / self.minor_axis ** 2 +
                             (point.y - self.center.y) ** 2 / self.major_axis ** 2, 1, abs_tol=abs_tol):
            return False
        new_point = self.frame.global_to_local_coordinates(point)
        u1, u2 = new_point.x / self.major_axis, new_point.y / self.minor_axis
        angle_new_point = volmdlr.geometry.sin_cos_angle(u1, u2)
        if self.angle_start < self.angle_end and self.angle_end >= angle_new_point >= self.angle_start:
            return True
        if self.angle_start > self.angle_end and self.angle_end <= angle_new_point <= self.angle_start:
            return True
        return False

    def abscissa(self, point: volmdlr.Point2D, tol: float = 1e-6):
        """
        Calculates the abscissa of a given point.

        :param point: point for calculating abscissa
        :return: a float, between 0 and the arc ellipse 2d's length
        """
        if point.point_distance(self.start) < tol:
            return 0
        if self.point_belongs(point):
            angle_abscissa = volmdlr.geometry.clockwise_angle(point - self.center, self.major_dir)
            angle_start = self.angle_start
            angle_end = angle_abscissa
            if self.angle_start > angle_abscissa > self.angle_end:
                angle_start = angle_abscissa
                angle_end = self.angle_start

            def arc_length(theta):
                return math.sqrt((self.major_axis ** 2) * math.sin(theta) ** 2 +
                                 (self.minor_axis ** 2) * math.cos(theta) ** 2)

            res, _ = scipy_integrate.quad(arc_length, angle_start, angle_end)
            return res
        raise ValueError(f'point {point} does not belong to ellipse')

    @property
    def bounding_rectangle(self):
        """
        Calculates the bounding rectangle for the arc ellipse 2d.

        :return: Bouding Rectangle object.
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
            area = math.pi * self.major_axis * self.minor_axis - 0.5 * self.major_axis * self.minor_axis * (
                    angle - math.sin(angle))
        else:
            angle = self.angle
            area = 0.5 * self.major_axis * self.minor_axis * (angle - math.sin(angle))

        if self.is_trigo:
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
        is_trigo = True
        if self.angle_start > self.angle_end:
            if self.angle_start >= self.angle_interior >= self.angle_end:
                angle_start = self.angle_end
                angle_end = self.angle_start
                is_trigo = False
            else:
                angle_end = self.angle_end + volmdlr.TWO_PI
                angle_start = self.angle_start
        elif self.angle_start == self.angle_end:
            angle_start = 0
            angle_end = 2 * math.pi
        else:
            angle_end = self.angle_end
            angle_start = self.angle_start

        discretization_points = [self.frame.local_to_global_coordinates(
            volmdlr.Point2D(self.major_axis * math.cos(angle), self.minor_axis * math.sin(angle)))
            for angle in npy.linspace(angle_start, angle_end, number_points)]
        if not is_trigo:
            discretization_points = discretization_points[::-1]
        return discretization_points

    def polygon_points(self, discretization_resolution: int):
        warnings.warn('polygon_points is deprecated,\
                please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(angle_resolution=discretization_resolution)

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
        :rtype: ArcEllipse3D
        """
        point_start3d = self.start.to_3d(plane_origin, x, y)
        point_interior3d = self.interior.to_3d(plane_origin, x, y)
        point_end3d = self.end.to_3d(plane_origin, x, y)
        point_center3d = self.center.to_3d(plane_origin, x, y)

        a_max2d = self.center + self.major_dir * self.major_axis
        a_max3d = a_max2d.to_3d(plane_origin, x, y)
        new_major_dir = a_max3d - point_center3d
        new_major_dir.normalize()
        extra3d = self.extra
        if extra3d:
            extra3d = self.extra.to_3d(plane_origin, x, y)
        return ArcEllipse3D(point_start3d, point_interior3d, point_end3d,
                            point_center3d, new_major_dir, extra3d, name=self.name)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            _, ax = plt.subplots()

        self.interior.plot(ax=ax, color='m')
        self.start.plot(ax=ax, color='r')
        self.end.plot(ax=ax, color='b')
        self.center.plot(ax=ax, color='y')

        x = []
        y = []
        for x_component, y_component in self.discretization_points(number_points=100):
            x.append(x_component)
            y.append(y_component)

        plt.plot(x, y, color=edge_style.color, alpha=edge_style.alpha)
        return ax

    def normal_vector(self, abscissa):
        raise NotImplementedError

    def direction_vector(self, abscissa):
        raise NotImplementedError

    def reverse(self):
        return self.__class__(self.end.copy(), self.interior.copy(), self.start.copy(),
                              self.center.copy(), self.major_dir.copy(), self.name)

    def line_intersections(self, line2d: Line2D):
        """
        Intersections between an Arc Ellipse 2D and a Line 2D.

        :param line2d: Line 2D to verify intersections
        :return: List with all intersections
        """
        ellipse2d_linesegment_intersections = vm_utils_intersections.ellipse2d_line_intersections(self, line2d)
        linesegment_intersections = []
        for inter in ellipse2d_linesegment_intersections:
            if self.point_belongs(inter):
                linesegment_intersections.append(inter)
        return linesegment_intersections

    def linesegment_intersections(self, linesegment2d: LineSegment2D):
        """
        Intersections between an Arc Ellipse 2D and a Line Segment 2D.

        :param linesegment2d: LineSegment 2D to verify intersections.
        :return: List with all intersections.
        """
        if not self.bounding_rectangle.b_rectangle_intersection(linesegment2d.bounding_rectangle):
            return []
        intersections = self.line_intersections(linesegment2d.to_line())
        linesegment_intersections = []
        for inter in intersections:
            if linesegment2d.point_belongs(inter):
                linesegment_intersections.append(inter)
        return linesegment_intersections

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and return a new Arc Ellipse 2D.

        side = 'old' or 'new'
        """
        if side == 'old':
            return ArcEllipse2D(frame.local_to_global_coordinates(self.start),
                                frame.local_to_global_coordinates(self.interior),
                                frame.local_to_global_coordinates(self.end),
                                frame.local_to_global_coordinates(self.center),
                                self.major_dir)
        if side == 'new':
            point_major_dir = self.center + self.major_dir * self.major_axis
            major_dir = frame.global_to_local_coordinates(point_major_dir).to_vector()
            major_dir.normalize()
            return ArcEllipse2D(frame.global_to_local_coordinates(self.start),
                                frame.global_to_local_coordinates(self.interior),
                                frame.global_to_local_coordinates(self.end),
                                frame.global_to_local_coordinates(self.center),
                                major_dir)
        raise ValueError('Side should be \'new\' \'old\'')

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point: Point to be verified
        :return: Return True if the point belongs to this surface,
            or False otherwise
        """
        raise NotImplementedError(f'the straight_line_point_belongs method must be'
                                  f' overloaded by {self.__class__.__name__}')


class FullArcEllipse(Edge):
    """
    Abstract class to define an ellipse.
    """

    def __init__(self, start_end: Union[volmdlr.Point2D, volmdlr.Point3D], major_axis: float, minor_axis: float,
                 center: Union[volmdlr.Point2D, volmdlr.Point3D],
                 major_dir: Union[volmdlr.Vector2D, volmdlr.Vector3D], name: str = ''):
        self.start_end = start_end
        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.center = center
        self.major_dir = major_dir

        Edge.__init__(self, start=start_end, end=start_end, name=name)

    def length(self):
        """
        Calculates the length of the ellipse.

        Ramanujan's approximation for the perimeter of the ellipse.
        P = math.pi * (a + b) [ 1 + (3h) / (10 + √(4 - 3h) ) ], where h = (a - b)**2/(a + b)**2
        :return: Perimeter of the ellipse
        :rtype: float
        """
        perimeter_formular_h = (self.major_axis - self.minor_axis) ** 2 / (self.major_axis + self.minor_axis) ** 2
        return math.pi * (self.major_axis + self.minor_axis) * \
            (1 + (3 * perimeter_formular_h / (10 + math.sqrt(4 - 3 * perimeter_formular_h))))

    def point_belongs(self, point: Union[volmdlr.Point2D, volmdlr.Point3D], abs_tol: float = 1e-6):
        """
        Verifies if a given point lies on the ellipse.

        :param point: point to be verified.
        :param abs_tol: Absolute tolerance to consider the point on the ellipse.
        :return: True is point lies on the ellipse, False otherwise
        """
        new_point = self.frame.global_to_local_coordinates(point)
        return math.isclose(new_point.x ** 2 / self.major_axis ** 2 +
                            new_point.y ** 2 / self.minor_axis ** 2, 1.0, abs_tol=abs_tol)

    def point_at_abscissa(self, abscissa: float, resolution: int = 2500):
        """
        Calculates a point on the FullArcEllipse at a given abscissa.

        :param abscissa: abscissa where in the curve the point should be calculated.
        :return: Corresponding point.
        """
        # TODO: enhance this method to a more precise method
        points = self.discretization_points(number_points=resolution)
        approx_abscissa = 0
        last_point = None
        for p1, p2 in zip(points[:-1], points[1:]):
            if approx_abscissa <= abscissa:
                approx_abscissa += p1.point_distance(p2)
                last_point = p2
            else:
                break
        return last_point

    def reverse(self):
        """
        Defines a new FullArcEllipse, identical to self, but in the opposite direction.

        """
        return self

    def straight_line_point_belongs(self, point):
        """
        Verifies if a point belongs to the surface created by closing the edge.

        :param point: Point to be verified
        :return: Return True if the point belongs to this surface,
            or False otherwise
        """
        raise NotImplementedError(f'the straight_line_point_belongs method must be'
                                  f' overloaded by {self.__class__.__name__}')

    def normal_vector(self, abscissa):
        """
        Calculates the normal vector the edge at given abscissa.

        :return: the normal vector
        """
        raise NotImplementedError

    def direction_vector(self, abscissa):
        """
        Calculates the direction vector the edge at given abscissa.

        :param abscissa: edge abscissa
        :return: direction vector
        """
        raise NotImplementedError

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


class FullArcEllipse2D(FullArcEllipse, ArcEllipse2D):
    """
    Defines a FullArcEllipse2D.
    """

    def __init__(self, start_end: volmdlr.Point2D, major_axis: float, minor_axis: float,
                 center: volmdlr.Point2D, major_dir: volmdlr.Vector2D, name: str = ''):
        major_dir.normalize()
        self.minor_dir = major_dir.deterministic_unit_normal_vector()
        self.frame = volmdlr.Frame2D(center, major_dir, self.minor_dir)
        self.theta = volmdlr.geometry.clockwise_angle(major_dir, volmdlr.X2D)
        if self.theta == math.pi * 2:
            self.theta = 0.0
        self._bounding_rectangle = None

        FullArcEllipse.__init__(self, start_end, major_axis, minor_axis, center, major_dir, name)

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """
        Calculates the discretized points for the ellipse.

        :param number_points: number of point to have in the discretized points.
        :param angle_resolution: the angle resolution to be used to discretize points.
        :return: discretized points.
        """
        if not number_points:
            number_points = math.ceil(volmdlr.TWO_PI * angle_resolution) + 2

        discretization_points = [self.center + volmdlr.Point2D(self.major_axis * math.cos(theta),
                                                               self.minor_axis * math.sin(theta))
                                 for theta in npy.linspace(0, volmdlr.TWO_PI, number_points)]
        discretization_points = [point.rotation(self.center, self.theta) for point in discretization_points]
        return discretization_points

    def to_3d(self, plane_origin, x, y):
        point_start_end3d = self.start_end.to_3d(plane_origin, x, y)
        point_center3d = self.center.to_3d(plane_origin, x, y)

        a_max2d = self.center + self.major_dir * self.major_axis
        a_max3d = a_max2d.to_3d(plane_origin, x, y)
        new_major_dir = (a_max3d - point_center3d).to_vector()
        new_major_dir.normalize()
        normal = x.cross(y)
        return FullArcEllipse3D(point_start_end3d, self.major_axis, self.minor_axis,
                                point_center3d, normal, new_major_dir, name=self.name)

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
        if side == 'old':
            return FullArcEllipse2D(frame.local_to_global_coordinates(self.start_end),
                                    self.major_axis, self.minor_axis,
                                    frame.local_to_global_coordinates(self.center),
                                    self.major_dir, self.name)
        if side == 'new':
            point_major_dir = self.center + self.major_dir * self.major_axis
            major_dir = frame.global_to_local_coordinates(point_major_dir).to_vector()
            major_dir.normalize()
            return FullArcEllipse2D(frame.global_to_local_coordinates(self.start_end),
                                    self.major_axis, self.minor_axis,
                                    frame.global_to_local_coordinates(self.center),
                                    major_dir, self.name)
        raise ValueError('Side should be \'new\' \'old\'')

    def translation(self, offset: volmdlr.Vector2D):
        """
        FullArcEllipse2D translation.

        :param offset: translation vector.
        :type offset: volmdlr.Vector2D
        :return: A new translated FullArcEllipse2D.
        :rtype: FullArcEllipse2D
        """
        return FullArcEllipse2D(self.start_end.translation(offset), self.major_axis, self.minor_axis,
                                self.center.translation(offset), self.major_dir, self.name)

    def abscissa(self, point: Union[volmdlr.Point2D, volmdlr.Point3D], tol: float = 1e-6):
        """
        Calculates the abscissa of a given point.

        :param point: point for calculating abscissa.
        :param tol: tolerance.
        :return: a float, between 0 and the ellipse's length.
        """
        if self.point_belongs(point):
            angle_abscissa = volmdlr.geometry.clockwise_angle(point - self.center, self.major_dir)
            angle_start = 0.0

            if angle_abscissa == volmdlr.TWO_PI:
                return self.length()

            def arc_length(theta):
                return math.sqrt((self.major_axis ** 2) * math.sin(theta) ** 2 +
                                 (self.minor_axis ** 2) * math.cos(theta) ** 2)

            res, _ = scipy_integrate.quad(arc_length, angle_start, angle_abscissa)
            return res
        raise ValueError(f'point {point} does not belong to ellipse')

    def normal_vector(self, abscissa):
        """
        Calculates the normal vector the edge at given abscissa.

        :return: the normal vector
        """
        raise NotImplementedError

    def direction_vector(self, abscissa):
        """
        Calculates the direction vector the edge at given abscissa.

        :param abscissa: edge abscissa
        :return: direction vector
        """
        raise NotImplementedError


class Line3D(Line):
    """
    Define an infinite line passing through the 2 points.

    """
    _non_eq_attributes = ['name', 'basis_primitives', 'bounding_box']

    def __init__(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D,
                 name: str = ''):
        Line.__init__(self, point1, point2, name=name)
        # self.points = [point1, point2]
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
        xmin = min([self.point1[0], self.point2[0]])
        xmax = max([self.point1[0], self.point2[0]])
        ymin = min([self.point1[1], self.point2[1]])
        ymax = max([self.point1[1], self.point2[1]])
        zmin = min([self.point1[2], self.point2[2]])
        zmax = max([self.point1[2], self.point2[2]])

        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def point_belongs(self, point3d):
        if point3d.is_close(self.point1):
            return True
        return self.direction_vector().is_colinear_to(point3d - self.point1)

    def point_distance(self, point):
        """Returns the minimal distance to a point."""
        vector1 = point - self.point1
        vector1.to_vector()
        vector2 = self.point2 - self.point1
        vector2.to_vector()
        return vector1.cross(vector2).norm() / vector2.norm()

    def line_distance(self, line2):
        """
        Calculates the distance between two Line3D.

        :param line2: other Line3D.
        :return: The distance between the two lines.
        """
        direction_vector1 = self.direction_vector()
        direction_vector2 = line2.direction_vector()
        if direction_vector1.is_colinear_to(direction_vector2):
            return direction_vector1.cross(line2.point1 - self.point1).norm() / direction_vector1.norm()
        vector = line2.point1 - self.point1
        line_distance = abs(vector.dot(direction_vector1.cross(direction_vector2))) / direction_vector1.cross(
            direction_vector2).norm()
        return line_distance

    def skew_to(self, line):
        """
        Verifies if two Line3D are skew to each other, that is, they are not parallel and never intersect.

        :param line: other line.
        :return: True if they are skew, False otherwise.
        """
        if self.direction_vector().is_colinear_to(line.direction_vector()):
            return False
        if math.isclose(self.line_distance(line), 0, abs_tol=1e-6):
            return False
        return True

    def intersection(self, line2):
        """
        Calculates the intersection between to Line3D, if there is an intersection.

        :param line: other Line3D
        :return: None if there is no intersection between Lines. A volmdlr.Point3D if there existes an intersection
        """
        direction_vector1 = self.direction_vector()
        direction_vector2 = line2.direction_vector()
        distance_to_line = self.line_distance(line2)
        if direction_vector1.is_colinear_to(direction_vector2) or \
                not math.isclose(distance_to_line, 0, abs_tol=1e-6):
            return None
        if math.isclose(distance_to_line, 0, abs_tol=1e-6) and \
                math.isclose(direction_vector1.dot(direction_vector2), 0, abs_tol=1e-6):
            projected_point, _ = self.point_projection(line2.point1)
            return projected_point
        vector = self.point1 - line2.point1
        t_coefficient = (
                                vector.dot(direction_vector2) * direction_vector2.dot(direction_vector1) -
                                vector.dot(direction_vector1) * direction_vector2.dot(direction_vector2)) / (
                                direction_vector1.dot(direction_vector1) * direction_vector2.dot(direction_vector2) -
                                direction_vector1.dot(direction_vector2) * direction_vector2.dot(direction_vector1))
        # u_coefficient = (vector.dot(direction_vector2) + t_coefficient * direction_vector1.dot(
        # direction_vector2)) / direction_vector2.dot(direction_vector2)
        intersection = self.point1 + t_coefficient * direction_vector1
        return intersection

    def plot(self, ax=None, color='k', alpha=1, dashed=True):
        if ax is None:
            ax = Axes3D(plt.figure())

        # Line segment
        ax.plot([self.point1.x, self.point2.x], [self.point1.y, self.point2.y],
                [self.point1.z, self.point2.z], color=color, alpha=alpha)

        # Drawing 3 times length of segment on each side
        u = self.point2 - self.point1
        v1 = self.point1 - 3 * u
        x1, y1, z1 = v1.x, v1.y, v1.z
        v2 = self.point2 - 3 * u
        x2, y2, z2 = v2.x, v2.y, v2.z
        if dashed:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=color,
                    dashes=[30, 5, 10, 5])
        else:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=color)
        return ax

    def plane_projection2d(self, center, x, y):
        return Line2D(self.point1.plane_projection2d(center, x, y),
                      self.point2.plane_projection2d(center, x, y))

    def minimum_distance_points(self, other_line):
        """
        Returns the points on this line and the other line that are the closest of lines.
        """
        u = self.point2 - self.point1
        v = other_line.point2 - other_line.point1
        w = self.point1 - other_line.point1
        u_dot_u = u.dot(u)
        u_dot_v = u.dot(v)
        v_dot_v = v.dot(v)
        u_dot_w = u.dot(w)
        v_dot_w = v.dot(w)

        s_param = (u_dot_v * v_dot_w - v_dot_v * u_dot_w) / (u_dot_u * v_dot_v - u_dot_v ** 2)
        t_param = (u_dot_u * v_dot_w - u_dot_v * u_dot_w) / (u_dot_u * v_dot_v - u_dot_v ** 2)
        point1 = self.point1 + s_param * u
        point2 = other_line.point1 + t_param * v
        return point1, point2

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Line3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Line3D
        """

        return Line3D(*[point.rotation(center, axis, angle) for point in
                        [self.point1, self.point2]])

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Line3D rotation. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in [self.point1, self.point2]:
            point.rotation_inplace(center, axis, angle)
        self._bbox = None

    def translation(self, offset: volmdlr.Vector3D):
        """
        Line3D translation.

        :param offset: translation vector
        :return: A new translated Line3D
        """
        return Line3D(*[point.translation(offset) for point in
                        [self.point1, self.point2]])

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Line3D translation. Object is updated inplace.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in [self.point1, self.point2]:
            point.translation_inplace(offset)
        self._bbox = None

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and return a new Line3D.

        side = 'old' or 'new'
        """
        if side == 'old':
            new_start = frame.local_to_global_coordinates(self.point1)
            new_end = frame.local_to_global_coordinates(self.point2)
        elif side == 'new':
            new_start = frame.global_to_local_coordinates(self.point1)
            new_end = frame.global_to_local_coordinates(self.point2)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        return Line3D(new_start, new_end)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes Line3D frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        if side == 'old':
            new_start = frame.local_to_global_coordinates(self.point1)
            new_end = frame.local_to_global_coordinates(self.point2)
        elif side == 'new':
            new_start = frame.global_to_local_coordinates(self.point1)
            new_end = frame.global_to_local_coordinates(self.point2)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        self.point1 = new_start
        self.point2 = new_end
        self._bbox = None

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D):
        if not self.point_belongs(point1) or not self.point_belongs(point2):
            raise ValueError('Point not on curve')

        return LineSegment3D(point1, point2)

    def copy(self, *args, **kwargs):
        return Line3D(*[point.copy() for point in [self.point1, self.point2]])

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to an Line3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated
        :type object_dict: dict
        :return: The corresponding Line3D object
        :rtype: :class:`volmdlr.edges.Line3D`
        """
        point1 = object_dict[arguments[1]]
        direction = object_dict[arguments[2]]
        point2 = point1 + direction
        return cls(point1, point2, arguments[0][1:-1])

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a Line3D into an Line2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Line2D.
        """
        p2d = [point.to_2d(plane_origin, x, y) for point in (self.point1, self.point2)]
        if p2d[0] == p2d[1]:
            return None
        return Line2D(*p2d, name=self.name)


class LineSegment3D(LineSegment):
    """
    Define a line segment limited by two points.

    """

    def __init__(self, start: volmdlr.Point3D, end: volmdlr.Point3D,
                 name: str = ''):
        if start.is_close(end):
            raise NotImplementedError
        # self.points = [start, end]
        LineSegment.__init__(self, start=start, end=end, name=name)
        self._bbox = None

    @property
    def bounding_box(self):
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def __hash__(self):
        return hash((self.__class__.__name__, self.start, self.end))

    def __eq__(self, other_linesegment3d):
        if other_linesegment3d.__class__ != self.__class__:
            return False
        return (self.start == other_linesegment3d.start
                and self.end == other_linesegment3d.end)

    def _bounding_box(self):

        xmin = min(self.start.x, self.end.x)
        xmax = max(self.start.x, self.end.x)
        ymin = min(self.start.y, self.end.y)
        ymax = max(self.start.y, self.end.y)
        zmin = min(self.start.z, self.end.z)
        zmax = max(self.start.z, self.end.z)

        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def to_dict(self, *args, **kwargs):
        return {'object_class': 'volmdlr.edges.LineSegment3D',
                'name': self.name,
                'start': self.start.to_dict(),
                'end': self.end.to_dict()
                }

    def normal_vector(self, abscissa=0.):
        return None

    def unit_normal_vector(self, abscissa=0.):
        return None

    # def middle_point(self):
    #     return self.point_at_abscissa(0.5 * self.length())

    def point_distance(self, point):
        """Returns the minimal distance to a point."""
        distance, point = volmdlr.LineSegment3DPointDistance(
            [(self.start.x, self.start.y, self.start.z),
             (self.end.x, self.end.y, self.end.z)],
            (point.x, point.y, point.z))
        return distance

    def plane_projection2d(self, center, x, y):
        start, end = self.start.plane_projection2d(center, x, y), self.end.plane_projection2d(center, x, y)
        if not start.is_close(end):
            return LineSegment2D(start, end)
        return None

    def line_intersections(self, line):
        line_self = self.to_line()
        if line_self.skew_to(line):
            return []
        intersection = line_self.intersection(line)
        if intersection and self.point_belongs(intersection):
            return [intersection]
        return []

    def linesegment_intersections(self, linesegment):
        line1 = self.to_line()
        line2 = linesegment.to_line()
        intersection = line1.intersection(line2)
        if intersection and self.point_belongs(intersection) and linesegment.point_belongs(intersection):
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

    def rotation_inplace(self, center: volmdlr.Point3D,
                         axis: volmdlr.Vector3D, angle: float):
        """
        Line2D rotation. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in self.points:
            point.rotation_inplace(center, axis, angle)
        self._bbox = None

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

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        LineSegment3D translation. Object is updated inplace.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in self.points:
            point.translation_inplace(offset)
        self._bbox = None

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

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        if side == 'old':
            new_start = frame.local_to_global_coordinates(self.start)
            new_end = frame.local_to_global_coordinates(self.end)
        elif side == 'new':
            new_start = frame.global_to_local_coordinates(self.start)
            new_end = frame.global_to_local_coordinates(self.end)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        self.start = new_start
        self.end = new_end
        self._bbox = None

    def copy(self, *args, **kwargs):
        """Returns a copy of the line segment."""
        return LineSegment3D(self.start.copy(), self.end.copy())

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

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
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        edge2d = self.plane_projection2d(volmdlr.O3D, x_3d, y_3d)
        edge2d.plot(ax=ax, edge_style=EdgeStyle(color=color, width=width))
        return ax

    def plot_data(self, x_3d, y_3d, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        edge2d = self.plane_projection2d(volmdlr.O3D, x_3d, y_3d)
        return edge2d.plot_data(marker, color, stroke_width,
                                dash, opacity, arrow)

    def to_line(self):
        """
        Converts the line segment into a line object.
        """
        return Line3D(self.start, self.end)

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
        bspline_curve = BSplineCurve3D.from_points_interpolation(points,
                                                                 degree)
        return bspline_curve

    def reverse(self):
        return LineSegment3D(self.end.copy(), self.start.copy())

    def minimum_distance_points(self, other_line):
        """
        Returns the points on this line and the other line that are the closest of lines.
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
        u = self.direction_vector()
        v = other_line.direction_vector()
        w = other_line.start - self.start

        a11 = u.dot(u)
        a12 = -u.dot(v)
        a22 = v.dot(v)

        a_matrix = npy.array([[a11, a12],
                              [a12, a22]])
        b_matrix = npy.array([w.dot(u), -w.dot(v)])

        res = lsq_linear(a_matrix, b_matrix, bounds=(0, 1))
        point1 = self.point_at_abscissa(res.x[0] * self.length())
        point2 = other_line.point_at_abscissa(res.x[1] * other_line.length())
        return point1, point2

    def parallel_distance(self, other_linesegment):
        pt_a, pt_b, pt_c = self.start, self.end, other_linesegment.start
        vector = volmdlr.Vector3D((pt_a - pt_b).vector)
        vector.normalize()
        plane1 = volmdlr.faces.Plane3D.from_3_points(pt_a, pt_b, pt_c)
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

    def minimum_distance(self, element, return_points=False):
        if element.__class__ is Arc3D or element.__class__ is volmdlr.wires.Circle3D:
            pt1, pt2 = element.minimum_distance_points_line(self)
            if return_points:
                return pt1.point_distance(pt2), pt1, pt2
            return pt1.point_distance(pt2)

        if element.__class__ is LineSegment3D:
            p1, p2 = self.matrix_distance(element)
            if return_points:
                return p1.point_distance(p2), p1, p2
            return p1.point_distance(p2)

        if element.__class__ is BSplineCurve3D:
            points = element.points
            lines = []
            dist_min = math.inf
            for p1, p2 in zip(points[0:-1], points[1:]):
                lines.append(LineSegment3D(p1, p2))
            for line in lines:
                p1, p2 = self.matrix_distance(line)
                dist = p1.point_distance(p2)
                if dist < dist_min:
                    dist_min = dist
                    min_points = (p1, p2)
            if return_points:
                p1, p2 = min_points
                return dist_min, p1, p2
            return dist_min

        raise NotImplementedError

    def extrusion(self, extrusion_vector):
        u = self.unit_direction_vector()
        v = extrusion_vector.copy()
        v.normalize()
        w = u.cross(v)
        length_1 = self.length()
        length_2 = extrusion_vector.norm()
        # outer_contour = Polygon2D([O2D, Point2D((l1, 0.)),
        #                            Point2D((l1, l2)), Point2D((0., l2))])
        plane = volmdlr.faces.Plane3D(volmdlr.Frame3D(self.start, u, v, w))
        return [plane.rectangular_cut(0, length_1, 0, length_2)]

    def _revolution_conical(self, params):
        axis, u, p1_proj, dist1, dist2, angle = params
        v = axis.cross(u)
        direction_vector = self.direction_vector()
        direction_vector.normalize()

        semi_angle = math.atan2(direction_vector.dot(u), direction_vector.dot(axis))
        cone_origin = p1_proj - dist1 / math.tan(semi_angle) * axis
        if semi_angle > 0.5 * math.pi:
            semi_angle = math.pi - semi_angle

            cone_frame = volmdlr.Frame3D(cone_origin, u, -v, -axis)
            angle2 = - angle
        else:
            angle2 = angle
            cone_frame = volmdlr.Frame3D(cone_origin, u, v, axis)

        surface = volmdlr.faces.ConicalSurface3D(cone_frame, semi_angle)
        return [surface.rectangular_cut(0, angle2, z1=dist1 / math.tan(semi_angle), z2=dist2 / math.tan(semi_angle))]

    def _cylindrical_revolution(self, params):
        axis, u, p1_proj, dist1, dist2, angle = params
        v = axis.cross(u)
        surface = volmdlr.faces.CylindricalSurface3D(volmdlr.Frame3D(p1_proj, u, v, axis), dist1)
        return [surface.rectangular_cut(0, angle, 0, (self.end - self.start).dot(axis))]

    def revolution(self, axis_point, axis, angle):
        """
        Returns the face generated by the revolution of the line segments.
        """
        axis_line3d = Line3D(axis_point, axis_point + axis)
        if axis_line3d.point_belongs(self.start) and axis_line3d.point_belongs(
                self.end):
            return []

        p1_proj, _ = axis_line3d.point_projection(self.start)
        p2_proj, _ = axis_line3d.point_projection(self.end)
        distance_1 = self.start.point_distance(p1_proj)
        distance_2 = self.end.point_distance(p2_proj)
        if not math.isclose(distance_1, 0., abs_tol=1e-9):
            u = self.start - p1_proj  # Unit vector from p1_proj to p1
            u.normalize()
        elif not math.isclose(distance_2, 0., abs_tol=1e-9):
            u = self.end - p2_proj  # Unit vector from p1_proj to p1
            u.normalize()
        else:
            return []
        if u.is_colinear_to(self.direction_vector()):
            # Planar face
            v = axis.cross(u)
            surface = volmdlr.faces.Plane3D(
                volmdlr.Frame3D(p1_proj, u, v, axis))
            smaller_r, bigger_r = sorted([distance_1, distance_2])
            if angle == volmdlr.TWO_PI:
                # Only 2 circles as contours
                outer_contour2d = volmdlr.wires.Circle2D(volmdlr.O2D, bigger_r)
                if not math.isclose(smaller_r, 0, abs_tol=1e-9):
                    inner_contours2d = [volmdlr.wires.Circle2D(volmdlr.O2D, smaller_r)]
                else:
                    inner_contours2d = []
            else:
                inner_contours2d = []
                if math.isclose(smaller_r, 0, abs_tol=1e-9):
                    # One arc and 2 lines (pizza slice)
                    arc2_e = volmdlr.Point2D(bigger_r, 0)
                    arc2_i = arc2_e.rotation(center=volmdlr.O2D,
                                             angle=0.5 * angle)
                    arc2_s = arc2_e.rotation(center=volmdlr.O2D, angle=angle)
                    arc2 = Arc2D(arc2_s, arc2_i, arc2_e)
                    line1 = LineSegment2D(arc2_e, volmdlr.O2D)
                    line2 = LineSegment2D(volmdlr.O2D, arc2_s)
                    outer_contour2d = volmdlr.wires.Contour2D([arc2, line1,
                                                               line2])

                else:
                    # Two arcs and lines
                    arc1_s = volmdlr.Point2D(bigger_r, 0)
                    arc1_i = arc1_s.rotation(center=volmdlr.O2D,
                                             angle=0.5 * angle)
                    arc1_e = arc1_s.rotation(center=volmdlr.O2D, angle=angle)
                    arc1 = Arc2D(arc1_s, arc1_i, arc1_e)

                    arc2_e = volmdlr.Point2D(smaller_r, 0)
                    arc2_i = arc2_e.rotation(center=volmdlr.O2D,
                                             angle=0.5 * angle)
                    arc2_s = arc2_e.rotation(center=volmdlr.O2D, angle=angle)
                    arc2 = Arc2D(arc2_s, arc2_i, arc2_e)

                    line1 = LineSegment2D(arc1_e, arc2_s)
                    line2 = LineSegment2D(arc2_e, arc1_s)

                    outer_contour2d = volmdlr.wires.Contour2D([arc1, line1,
                                                               arc2, line2])

            return [volmdlr.faces.PlaneFace3D(surface,
                                              volmdlr.faces.Surface2D(
                                                  outer_contour2d,
                                                  inner_contours2d))]

        if not math.isclose(distance_1, distance_2, abs_tol=1e-9):
            # Conical
            return self._revolution_conical([axis, u, p1_proj, distance_1, distance_2, angle])

        # Cylindrical face
        return self._cylindrical_revolution([axis, u, p1_proj, distance_1, distance_2, angle])


class BSplineCurve3D(BSplineCurve):
    """
    A class for 3 dimensional B-spline curves.

    The following rule must be respected : `number of knots = number of control points + degree + 1`

    :param degree: The degree of the 3 dimensional B-spline curve
    :type degree: int
    :param control_points: A list of 3 dimensional points
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
    _non_serializable_attributes = ['curve']

    def __init__(self,
                 degree: int,
                 control_points: List[volmdlr.Point3D],
                 knot_multiplicities: List[int],
                 knots: List[float],
                 weights: List[float] = None,
                 periodic: bool = False,
                 name: str = ''):

        BSplineCurve.__init__(self, degree,
                              control_points,
                              knot_multiplicities,
                              knots,
                              weights,
                              periodic,
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
        bbox = self.curve.bbox
        return volmdlr.core.BoundingBox(bbox[0][0], bbox[1][0],
                                        bbox[0][1], bbox[1][1],
                                        bbox[0][2], bbox[1][2])

    def look_up_table(self, resolution: int = 20, start_parameter: float = 0,
                      end_parameter: float = 1):
        """
        Creates a table of equivalence between the parameter t (eval. of the BSplineCurve) and the cumulative distance.

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
                point1 = volmdlr.Point3D(*self.curve.evaluate_single(param1))
                point2 = volmdlr.Point3D(*self.curve.evaluate_single(param2))
                distance += point1.point_distance(point2)
                yield param2, distance

    def normal(self, position: float = 0.0):
        _, normal = operations.normal(self.curve, position, normalize=True)
        normal = volmdlr.Point3D(normal[0], normal[1], normal[2])
        return normal

    def direction_vector(self, abscissa=0.):
        length = self.length()
        if abscissa >= length:
            abscissa2 = length
            abscissa = abscissa2 - 0.001 * length

        else:
            abscissa2 = min(abscissa + 0.001 * length, length)

        tangent = self.point_at_abscissa(abscissa2) - self.point_at_abscissa(
            abscissa)
        return tangent

    def point3d_to_parameter(self, point: volmdlr.Point3D):
        """
        Search for the value of the normalized evaluation parameter t (between 0 and 1).

        :return: the given point when the BSplineCurve3D is evaluated at the t value.
        """

        def fun(param):
            p3d = volmdlr.Point3D(*self.curve.evaluate_single(param))
            return point.point_distance(p3d)

        res = minimize(fun=fun, x0=0.5, bounds=[(0, 1)], tol=1e-9)
        return res.x[0]

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
        # curve_form = arguments[3]
        if arguments[4] == '.F.':
            closed_curve = False
        elif arguments[4] == '.T.':
            closed_curve = True
        else:
            raise ValueError
        # self_intersect = arguments[5]
        knot_multiplicities = [int(i) for i in arguments[6][1:-1].split(",")]
        knots = [float(i) for i in arguments[7][1:-1].split(",")]
        # knot_spec = arguments[8]
        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot] * knot_multiplicities[i])

        if 9 in range(len(arguments)):
            weight_data = [float(i) for i in arguments[9][1:-1].split(",")]
        else:
            weight_data = None

        # FORCING CLOSED_CURVE = FALSE:
        # closed_curve = False
        return cls(degree, points, knot_multiplicities, knots, weight_data,
                   closed_curve, name)

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

        if surface_id:
            content += f"#{curve_id + 1} = SURFACE_CURVE('',#{curve_id},(#{curve_id + 2}),.PCURVE_S1.);\n"
            content += f"#{curve_id + 2} = PCURVE('',#{surface_id},#{curve_id + 3});\n"

            # 2D parametric curve
            curve2d_content, (curve2d_id,) = curve2d.to_step(curve_id + 3)  # 5

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
            content += f"#{current_id} = EDGE_CURVE('{self.name}',#{start_id},#{end_id},#{curve_id + 1},.T.);\n"
        else:
            content += f"#{current_id} = EDGE_CURVE('{self.name}',#{start_id},#{end_id},#{curve_id},.T.);\n"
        return content, [current_id]

    def point_distance(self, pt1):
        """Returns the minimal distance to a point."""
        distances = []
        for point in self.points:
            distances.append(pt1.point_distance(point))
        return min(distances)

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
                                            self.knots, self.weights,
                                            self.periodic, self.name)
        return new_bsplinecurve3d

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        BSplineCurve3D rotation. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_control_points = [point.rotation(center, axis, angle) for point in
                              self.control_points]
        new_bsplinecurve3d = BSplineCurve3D(self.degree, new_control_points,
                                            self.knot_multiplicities,
                                            self.knots, self.weights,
                                            self.periodic, self.name)
        self.control_points = new_control_points
        self.curve = new_bsplinecurve3d.curve
        self.points = new_bsplinecurve3d.points
        self._bbox = None

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D):
        if (point1.is_close(self.start) and point2.is_close(self.end)) \
                or (point1.is_close(self.end) and point2.is_close(self.start)):
            return self

        if point1.is_close(self.start) and not point2.is_close(self.end):
            return self.cut_after(self.point3d_to_parameter(point2))

        if point2.is_close(self.start) and not point1.is_close(self.end):
            return self.cut_after(self.point3d_to_parameter(point1))

        if not point1.is_close(self.start) and point2.is_close(self.end):
            return self.cut_before(self.point3d_to_parameter(point1))

        if not point2.is_close(self.start) and point1.is_close(self.end):
            return self.cut_before(self.point3d_to_parameter(point2))

        parameter1 = self.point3d_to_parameter(point1)
        parameter2 = self.point3d_to_parameter(point2)
        if parameter1 is None or parameter2 is None:
            raise ValueError('Point not on BSplineCurve for trim method')

        if parameter1 > parameter2:
            parameter1, parameter2 = parameter2, parameter1
            point1, point2 = point2, point1

        bspline_curve = self.cut_before(parameter1)
        new_param2 = bspline_curve.point3d_to_parameter(point2)
        trimmed_bspline_cruve = bspline_curve.cut_after(new_param2)
        return trimmed_bspline_cruve

    def trim_between_evaluations(self, parameter1: float, parameter2: float):
        print('Use BSplineCurve3D.trim instead of trim_between_evaluation')
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
        new_knots = standardize_knot_vector(new_knots)

        return BSplineCurve3D(degree=bspline_curve.degree,
                              control_points=new_ctrlpts,
                              knot_multiplicities=new_multiplicities,
                              knots=new_knots,
                              weights=None,
                              periodic=bspline_curve.periodic,
                              name=bspline_curve.name)

    def cut_before(self, parameter: float):
        """
        Returns the right side of the splitted curve at a given parameter.

        :param parameter: parameter value that specifies where to split the curve.
        :type parameter: float
        """
        # Is a value of parameter below 4e-3 a real need for precision ?
        if math.isclose(parameter, 0, abs_tol=4e-3):
            return self
        if math.isclose(parameter, 1, abs_tol=4e-3):
            return self.reverse()
        #     raise ValueError('Nothing will be left from the BSplineCurve3D')

        curves = operations.split_curve(self.curve, parameter)
        return self.from_geomdl_curve(curves[1])

    def cut_after(self, parameter: float):
        """
        Returns the left side of the splitted curve at a given parameter.

        :param parameter: parameter value that specifies where to split the curve.
        :type parameter: float
        """
        # Is a value of parameter below 4e-3 a real need for precision ?
        if math.isclose(parameter, 0, abs_tol=1e-6):
            #     # raise ValueError('Nothing will be left from the BSplineCurve3D')
            #     curves = operations.split_curve(operations.refine_knotvector(self.curve, [4]), parameter)
            #     return self.from_geomdl_curve(curves[0])
            return self.reverse()
        if math.isclose(parameter, 1, abs_tol=4e-3):
            return self
        curves = operations.split_curve(self.curve, parameter)
        return self.from_geomdl_curve(curves[0])

    def insert_knot(self, knot: float, num: int = 1):
        """
        Returns a new BSplineCurve3D.

        """
        curve_copy = self.curve.__deepcopy__({})
        modified_curve = operations.insert_knot(curve_copy, [knot], num=[num])
        return self.from_geomdl_curve(modified_curve)

    # Copy paste du LineSegment3D
    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        x = [point.x for point in self.points]
        y = [point.y for point in self.points]
        z = [point.z for point in self.points]
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
                              self.weights, self.periodic, self.name)

    def polygon_points(self, discretization_resolution: int):
        warnings.warn('polygon_points is deprecated,\
                please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(angle_resolution=discretization_resolution)

    def curvature(self, u: float, point_in_curve: bool = False):
        # u should be in the interval [0,1]
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
        # print(self.global_maximum_curvature(point_in_curve))
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

    def global_minimum_curvature(self, nb_eval: int = 21):
        check = [i / (nb_eval - 1) for i in range(nb_eval)]
        radius = []
        for u in check:
            radius.append(self.minimum_curvature(u))
        return radius

    def triangulation(self):
        return None

    def linesegment_intersections(self, linesegment3d: LineSegment3D):
        """
        Calculates intersections between a BSplineCurve3D and a LineSegment3D.

        :param linesegment3d: linesegment to verify intersections.
        :return: list with the intersections points.
        """
        if not self.bounding_box.bbox_intersection(linesegment3d.bounding_box):
            return []
        intersections_points = self.get_linesegment_intersections(linesegment3d)
        return intersections_points


class BezierCurve3D(BSplineCurve3D):
    """
    A class for 3 dimensional Bézier curves.

    :param degree: The degree of the Bézier curve
    :type degree: int
    :param control_points: A list of 3 dimensional points
    :type control_points: List[:class:`volmdlr.Point3D`]
    :param name: The name of the B-spline curve. Default value is ''
    :type name: str, optional
    """

    def __init__(self, degree: int, control_points: List[volmdlr.Point3D],
                 name: str = ''):
        knotvector = utilities.generate_knot_vector(degree,
                                                    len(control_points))
        knot_multiplicity = [1] * len(knotvector)

        BSplineCurve3D.__init__(self, degree, control_points,
                                knot_multiplicity, knotvector,
                                None, False, name)


class Arc3D(Arc):
    """
    An arc is defined by a starting point, an end point and an interior point.

    """

    def __init__(self, start, interior, end, name=''):
        self._utd_normal = False
        self._utd_center = False
        self._utd_frame = False
        self._utd_is_trigo = False
        self._utd_angle = False
        self._normal = None
        self._frame = None
        self._center = None
        self._is_trigo = None
        self._angle = None
        # self._utd_clockwise_and_trigowise_paths = False
        Arc.__init__(self, start=start, end=end, interior=interior, name=name)
        self._bbox = None

    def __hash__(self):
        return hash(('arc3d', self.interior, self.start, self.end))

    def __eq__(self, other_arc):
        if self.__class__.__name__ != other_arc.__class__.__name__:
            return False
        return (self.center == other_arc.center
                and self.start == other_arc.start
                and self.end == other_arc.end
                and self.interior == other_arc.interior)

    @property
    def bounding_box(self):
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

        points = self.discretization_points(angle_resolution=10)
        xmin = min(point.x for point in points)
        xmax = max(point.x for point in points)
        ymin = min(point.y for point in points)
        ymax = max(point.y for point in points)
        zmin = min(point.z for point in points)
        zmax = max(point.z for point in points)
        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    @classmethod
    def from_angle(cls, start: volmdlr.Point3D, angle: float,
                   axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D):
        """Gives the arc3D from a start, an angle and an axis."""
        start_gen = start
        int_gen = start_gen.rotation(axis_point, axis, angle / 2)
        end_gen = start_gen.rotation(axis_point, axis, angle)
        if angle == volmdlr.TWO_PI:
            line = Line3D(axis_point, axis_point + axis)
            center, _ = line.point_projection(start)
            radius = center.point_distance(start)
            u = start - center
            v = axis.cross(u)
            return volmdlr.wires.Circle3D(volmdlr.Frame3D(center, u, v, axis),
                                          radius)
        return cls(start_gen, int_gen, end_gen, axis)

    @property
    def normal(self):
        if not self._utd_normal:
            self._normal = self.get_normal()
            self._utd_normal = True
        return self._normal

    def get_normal(self):
        u1 = self.interior - self.start
        u2 = self.interior - self.end
        try:
            u1.normalize()
            u2.normalize()
        except ZeroDivisionError:
            raise ValueError(
                'Start, end and interior points of an arc must be distincts') from ZeroDivisionError

        normal = u2.cross(u1)
        normal.normalize()
        return normal

    @property
    def center(self):
        if not self._utd_center:
            self._center = self.get_center()
            self._utd_center = True
        return self._center

    def get_center(self):
        vector_u1 = self.interior - self.start
        vector_u2 = self.interior - self.end
        if vector_u1.is_close(vector_u2):
            vector_u2 = self.normal.cross(vector_u1)
            vector_u2.normalize()

        vector_v1 = self.normal.cross(vector_u1)  # v1 is normal, equal u2
        vector_v2 = self.normal.cross(vector_u2)  # equal -u1

        point11 = 0.5 * (self.start + self.interior)  # Mid-point of segment s,m
        point12 = point11 + vector_v1
        point21 = 0.5 * (self.end + self.interior)  # Mid-point of segment s,m
        point22 = point21 + vector_v2

        line_1 = Line3D(point11, point12)
        line_2 = Line3D(point21, point22)

        try:
            center, _ = line_1.minimum_distance_points(line_2)
        except ZeroDivisionError:
            raise ValueError(
                'Start, end and interior points  of an arc must be distincts') from ZeroDivisionError

        return center

    @property
    def frame(self):
        if not self._utd_frame:
            self._frame = self.get_frame()
            self._utd_frame = True
        return self._frame

    def get_frame(self):
        vec1 = self.start - self.center
        vec1.normalize()
        vec2 = self.normal.cross(vec1)
        frame = volmdlr.Frame3D(self.center, vec1, vec2, self.normal)
        return frame

    @property
    def is_trigo(self):
        if not self._utd_is_trigo:
            self._is_trigo = self.get_arc_direction()
            self._utd_is_trigo = True
        return self._is_trigo

    def get_arc_direction(self):
        """
        Verifies if arc is clockwise or counterclockwise.

        :return: True if clockwise, False if counterclockwise.
        """
        clockwise_path, trigowise_path = self.clockwise_and_trigowise_paths
        if clockwise_path > trigowise_path:
            return True
        return False

    @property
    def clockwise_and_trigowise_paths(self):
        """
        :return: clockwise path and trigonometric path property.
        """
        if not self._utd_clockwise_and_trigowise_paths:
            vec1 = self.start - self.center
            vec1.normalize()
            vec2 = self.normal.cross(vec1)
            radius_1 = self.start.to_2d(self.center, vec1, vec2)
            radius_2 = self.end.to_2d(self.center, vec1, vec2)
            radius_i = self.interior.to_2d(self.center, vec1, vec2)
            self._clockwise_and_trigowise_paths = \
                self.get_clockwise_and_trigowise_paths(radius_1,
                                                       radius_2,
                                                       radius_i)
            self._utd_clockwise_and_trigowise_paths = True
        return self._clockwise_and_trigowise_paths

    @property
    def angle(self):
        """
        Arc angle property.

        :return: arc angle.
        """
        if not self._utd_angle:
            self._angle = self.get_angle()
            self._utd_angle = True
        return self._angle

    def get_angle(self):
        """
        Gets the arc angle.

        :return: arc angle.
        """
        clockwise_path, trigowise_path = \
            self.clockwise_and_trigowise_paths
        if self.is_trigo:
            return trigowise_path
        return clockwise_path

    @property
    def points(self):
        return [self.start, self.interior, self.end]

    def reverse(self):
        """
        Defines a new Arc3D, identical to self, but in the opposite direction.

        """
        return self.__class__(self.end.copy(),
                              self.interior.copy(),
                              self.start.copy())

    def point_at_abscissa(self, abscissa):
        """
        Calculates a point in the Arc3D at a given abscissa.

        :param abscissa: abscissa where in the curve the point should be calculated.
        :return: Corresponding point.
        """
        return self.start.rotation(self.center, self.normal, abscissa / self.radius)

    def direction_vector(self, abscissa):
        """
        Calculates a direction vector at a given abscissa of the Arc3D.

        :param abscissa: abscissa where in the curve the direction vector should be calculated.
        :return: Corresponding direction vector.
        """
        normal_vector = self.normal_vector(abscissa)
        tangent = normal_vector.cross(self.normal)
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
        new_start = self.start.rotation(center, axis, angle)
        new_interior = self.interior.rotation(center, axis, angle)
        new_end = self.end.rotation(center, axis, angle)
        return Arc3D(new_start, new_interior, new_end, name=self.name)

    def rotation_inplace(self, center: volmdlr.Point3D,
                         axis: volmdlr.Vector3D, angle: float):
        """
        Arc3D rotation. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.center.rotation_inplace(center, axis, angle)
        self.start.rotation_inplace(center, axis, angle)
        self.interior.rotation_inplace(center, axis, angle)
        self.end.rotation_inplace(center, axis, angle)
        self._bbox = None

    def translation(self, offset: volmdlr.Vector3D):
        """
        Arc3D translation.

        :param offset: translation vector.
        :return: A new translated Arc3D.
        """
        new_start = self.start.translation(offset)
        new_interior = self.interior.translation(offset)
        new_end = self.end.translation(offset)
        return Arc3D(new_start, new_interior, new_end, name=self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Arc3D translation. Object is updated inplace.

        :param offset: translation vector.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.center.translation_inplace(offset)
        self.start.translation_inplace(offset)
        self.interior.translation_inplace(offset)
        self.end.translation_inplace(offset)
        self._bbox = None

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')
        # if plot_points:
        #     ax.plot([self.interior[0]], [self.interior[1]], [self.interior[2]],
        #             color='b')
        #     ax.plot([self.start[0]], [self.start[1]], [self.start[2]], c='r')
        #     ax.plot([self.end[0]], [self.end[1]], [self.end[2]], c='r')
        #     ax.plot([self.interior[0]], [self.interior[1]], [self.interior[2]],
        #             c='g')
        x = []
        y = []
        z = []
        for pointx, pointy, pointz in self.discretization_points(number_points=25):
            x.append(pointx)
            y.append(pointy)
            z.append(pointz)

        ax.plot(x, y, z, color=edge_style.color, alpha=edge_style.alpha)
        if edge_style.edge_ends:
            self.start.plot(ax=ax)
            self.end.plot(ax=ax)

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

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

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
        return Arc3D(self.start.copy(), self.interior.copy(), self.end.copy())

    def frame_mapping_parameters(self, frame: volmdlr.Frame3D, side: str):
        if side == 'old':
            new_start = frame.local_to_global_coordinates(self.start.copy())
            new_interior = frame.local_to_global_coordinates(self.interior.copy())
            new_end = frame.local_to_global_coordinates(self.end.copy())
        elif side == 'new':
            new_start = frame.global_to_local_coordinates(self.start.copy())
            new_interior = frame.global_to_local_coordinates(self.interior.copy())
            new_end = frame.global_to_local_coordinates(self.end.copy())
        else:
            raise ValueError('side value not valid, please specify'
                             'a correct value: \'old\' or \'new\'')
        return new_start, new_interior, new_end

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and return a new Arc3D.

        side = 'old' or 'new'
        """
        new_start, new_interior, new_end = \
            self.frame_mapping_parameters(frame, side)

        return Arc3D(new_start, new_interior, new_end, name=self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_start, new_interior, new_end = \
            self.frame_mapping_parameters(frame, side)
        self.start, self.interior, self.end = new_start, new_interior, new_end
        self._bbox = None

    def abscissa(self, point: volmdlr.Point3D, tol: float = 1e-6):
        """
        Calculates the abscissa given a point in the Arc3D.

        :param point: point to calculate the abscissa.
        :param tol: (Optional) Confusion distance to consider points equal. Default 1e-6.
        :return: corresponding abscissa.
        """
        if point.point_distance(self.start) < tol:
            return 0
        if point.point_distance(self.end) < tol:
            return self.length()
        x, y, _ = self.frame.global_to_local_coordinates(point)
        u1 = x / self.radius
        u2 = y / self.radius
        theta = volmdlr.geometry.sin_cos_angle(u1, u2)

        return self.radius * abs(theta)

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a Arc3D into an Arc2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: Arc2D.
        """
        point_start = self.start.to_2d(plane_origin, x, y)
        point_interior = self.interior.to_2d(plane_origin, x, y)
        point_end = self.end.to_2d(plane_origin, x, y)
        return Arc2D(point_start, point_interior, point_end, name=self.name)

    def minimum_distance_points_arc(self, other_arc):

        u1 = self.start - self.center
        u1.normalize()
        u2 = self.normal.cross(u1)

        w = other_arc.center - self.center

        u3 = other_arc.start - other_arc.center
        u3.normalize()
        u4 = other_arc.normal.cross(u3)

        r1, r2 = self.radius, other_arc.radius

        a, b, c, d = u1.dot(u1), u1.dot(u2), u1.dot(u3), u1.dot(u4)
        e, f, g = u2.dot(u2), u2.dot(u3), u2.dot(u4)
        h, i = u3.dot(u3), u3.dot(u4)
        j = u4.dot(u4)
        k, l, m, n, o = w.dot(u1), w.dot(u2), w.dot(u3), w.dot(u4), w.dot(w)

        def distance_squared(x):
            return (a * ((math.cos(x[0])) ** 2) * r1 ** 2 + e * (
                    (math.sin(x[0])) ** 2) * r1 ** 2
                    + o + h * ((math.cos(x[1])) ** 2) * r2 ** 2 + j * (
                            (math.sin(x[1])) ** 2) * r2 ** 2
                    + b * math.sin(2 * x[0]) * r1 ** 2 - 2 * r1 * math.cos(
                        x[0]) * k
                    - 2 * r1 * r2 * math.cos(x[0]) * math.cos(x[1]) * c
                    - 2 * r1 * r2 * math.cos(x[0]) * math.sin(
                        x[1]) * d - 2 * r1 * math.sin(x[0]) * l
                    - 2 * r1 * r2 * math.sin(x[0]) * math.cos(x[1]) * f
                    - 2 * r1 * r2 * math.sin(x[0]) * math.sin(
                        x[1]) * g + 2 * r2 * math.cos(x[1]) * m
                    + 2 * r2 * math.sin(x[1]) * n + i * math.sin(
                        2 * x[1]) * r2 ** 2)

        x01 = npy.array([self.angle / 2, other_arc.angle / 2])

        res1 = least_squares(distance_squared, x01, bounds=[(0, 0), (self.angle, other_arc.angle)])

        point1 = self.point_at_abscissa(res1.x[0] * r1)
        point2 = other_arc.point_at_abscissa(res1.x[1] * r2)

        return point1, point2

    def distance_squared(self, x, u, v, k, w):
        radius = self.radius
        return (u.dot(u) * x[0] ** 2 + w.dot(w) + v.dot(v) * (
                (math.sin(x[1])) ** 2) * radius ** 2 + k.dot(k) * ((math.cos(x[1])) ** 2) * radius ** 2
                - 2 * x[0] * w.dot(u) - 2 * x[0] * radius * math.sin(x[1]) * u.dot(v) - 2 * x[
                    0] * radius * math.cos(x[1]) * u.dot(k)
                + 2 * radius * math.sin(x[1]) * w.dot(v) + 2 * radius * math.cos(x[1]) * w.dot(k)
                + math.sin(2 * x[1]) * v.dot(k) * radius ** 2)

    def minimum_distance_points_line(self, other_line):
        u = other_line.direction_vector()
        k = self.start - self.center
        k.normalize()
        w = self.center - other_line.start
        v = self.normal.cross(k)

        radius = self.radius

        x01 = npy.array([0.5, self.angle / 2])
        x02 = npy.array([0.5, 0])
        x03 = npy.array([0.5, self.angle])

        res1 = least_squares(self.distance_squared, x01, bounds=[(0, 0), (1, self.angle)], args=(u, v, k, w))
        res2 = least_squares(self.distance_squared, x02, bounds=[(0, 0), (1, self.angle)], args=(u, v, k, w))
        res3 = least_squares(self.distance_squared, x03, bounds=[(0, 0), (1, self.angle)], args=(u, v, k, w))

        point1 = other_line.point_at_abscissa(res1.x[0] * other_line.length())
        point2 = self.point_at_abscissa(res1.x[1] * radius)

        res = [res2, res3]
        for couple in res:
            ptest1 = other_line.point_at_abscissa(
                couple.x[0] * other_line.length())
            ptest2 = self.point_at_abscissa(couple.x[1] * radius)
            dtest = ptest1.point_distance(ptest2)
            if dtest < v.dot(v):
                point1, point2 = ptest1, ptest2

        return point1, point2

    def minimum_distance(self, element, return_points=False):
        if element.__class__ is Arc3D or element.__class__.__name__ == 'Circle3D':
            p1, p2 = self.minimum_distance_points_arc(element)
            if return_points:
                return p1.point_distance(p2), p1, p2
            return p1.point_distance(p2)

        if element.__class__ is LineSegment3D:
            pt1, pt2 = self.minimum_distance_points_line(element)
            if return_points:
                return pt1.point_distance(pt2), pt1, pt2
            return pt1.point_distance(pt2)

        return NotImplementedError

    def extrusion(self, extrusion_vector):
        if self.normal.is_colinear_to(extrusion_vector):
            u = self.start - self.center
            u.normalize()
            w = extrusion_vector.copy()
            w.normalize()
            v = w.cross(u)
            arc2d = self.to_2d(self.center, u, v)
            angle1, angle2 = arc2d.angle1, arc2d.angle2
            if angle2 < angle1:
                angle2 += volmdlr.TWO_PI
            cylinder = volmdlr.faces.CylindricalSurface3D(
                volmdlr.Frame3D(self.center,
                                u,
                                v,
                                w),
                self.radius
            )
            return [cylinder.rectangular_cut(angle1, angle2, 0., extrusion_vector.norm())]
        raise NotImplementedError(f'Elliptic faces not handled: dot={self.normal.dot(extrusion_vector)}')

    def revolution(self, axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                   angle: float):
        line3d = Line3D(axis_point, axis_point + axis)
        tore_center, _ = line3d.point_projection(self.center)

        # Sphere
        if math.isclose(tore_center.point_distance(self.center), 0.,
                        abs_tol=1e-6):

            start_p, _ = line3d.point_projection(self.start)
            u = self.start - start_p

            if math.isclose(u.norm(), 0, abs_tol=1e-6):
                end_p, _ = line3d.point_projection(self.end)
                u = self.end - end_p
                if math.isclose(u.norm(), 0, abs_tol=1e-6):
                    interior_p, _ = line3d.point_projection(self.interior)
                    u = self.interior - interior_p

            u.normalize()
            v = axis.cross(u)
            arc2d = self.to_2d(self.center, u, axis)

            surface = volmdlr.faces.SphericalSurface3D(
                volmdlr.Frame3D(self.center, u, v, axis), self.radius)

            return [surface.rectangular_cut(0, angle,
                                            arc2d.angle1, arc2d.angle2)]

        # Toroidal
        u = self.center - tore_center
        u.normalize()
        v = axis.cross(u)
        if not math.isclose(self.normal.dot(u), 0., abs_tol=1e-6):
            raise NotImplementedError(
                'Outside of plane revolution not supported')

        radius = tore_center.point_distance(self.center)
        surface = volmdlr.faces.ToroidalSurface3D(
            volmdlr.Frame3D(tore_center, u, v, axis), radius,
            self.radius)
        arc2d = self.to_2d(tore_center, u, axis)
        return [surface.rectangular_cut(0, angle,
                                        arc2d.angle1, arc2d.angle2)]

    def to_step(self, current_id, surface_id=None):
        """Exports to STEP format."""
        if self.angle >= math.pi:
            length = self.length()
            arc1, arc2 = self.split(self.point_at_abscissa(0.33 * length))
            arc2, arc3 = arc2.split(self.point_at_abscissa(0.66 * length))
            content, arcs1_id = arc1.to_step_without_splitting(current_id)
            arc2_content, arcs2_id = arc2.to_step_without_splitting(
                arcs1_id[0] + 1)
            arc3_content, arcs3_id = arc3.to_step_without_splitting(
                arcs2_id[0] + 1)
            content += arc2_content + arc3_content
            return content, [arcs1_id[0], arcs2_id[0], arcs3_id[0]]
        return self.to_step_without_splitting(current_id)

    def to_step_without_splitting(self, current_id, surface_id=None):
        u = self.start - self.center
        u.normalize()
        v = self.normal.cross(u)
        frame = volmdlr.Frame3D(self.center, self.normal, u, v)

        content, frame_id = frame.to_step(current_id)
        curve_id = frame_id + 1
        content += f"#{curve_id} = CIRCLE('{self.name}', #{frame_id}, {self.radius * 1000});\n"

        if surface_id:
            content += f"#{curve_id + 1} = SURFACE_CURVE('',#{curve_id},(#{surface_id}),.PCURVE_S1.);\n"
            curve_id += 1

        current_id = curve_id + 1
        start_content, start_id = self.start.to_step(current_id, vertex=True)
        end_content, end_id = self.end.to_step(start_id + 1, vertex=True)
        content += start_content + end_content
        current_id = end_id + 1
        content += f"#{current_id} = EDGE_CURVE('{self.name}',#{start_id},#{end_id},#{curve_id},.T.);\n"
        return content, [current_id]

    def point_belongs(self, point3d, abs_tol: float = 1e-6):
        """
        Check if a point 3d belongs to the arc_3d or not.

        :param point3d: point to be verified is on arc
        :return: True if point is on Arc, False otherwise.
        """
        if not math.isclose(point3d.point_distance(self.center), self.radius, abs_tol=abs_tol):
            return False
        # vector1 = self.start - self.center
        # vector2 = self.interior - self.center
        vector = point3d - self.center
        if not math.isclose(vector.dot(self.frame.w), 0.0, abs_tol=abs_tol):
            return False
        point_abscissa = self.abscissa(point3d)
        abscissa_start = self.abscissa(self.start)
        abscissa_end = self.abscissa(self.end)
        if abscissa_start <= point_abscissa <= abscissa_end:
            return True
        return False

    def triangulation(self):
        """
        Triangulation for an Arc3D.

        """
        return None

    def middle_point(self):
        return self.point_at_abscissa(self.length() / 2)

    def line_intersections(self, line3d: Line3D):
        """
        Calculates intersections between an Arc3D and a Line3D.

        :param linesegment3d: linesegment to verify intersections.
        :return: list with intersections points between line and Arc3D.
        """
        circle3d_lineseg_inters = vm_utils_intersections.circle_3d_line_intersections(self, line3d)
        linesegment_intersections = []
        for intersection in circle3d_lineseg_inters:
            if self.point_belongs(intersection, 1e-6):
                linesegment_intersections.append(intersection)
        return linesegment_intersections

    def linesegment_intersections(self, linesegment3d: LineSegment3D):
        """
        Calculates intersections between an Arc3D and a LineSegment3D.

        :param linesegment3d: linesegment to verify intersections.
        :return: list with intersections points between linesegment and Arc3D.
        """
        linesegment_intersections = []
        intersections = self.line_intersections(linesegment3d.to_line())
        for intersection in intersections:
            if linesegment3d.point_belongs(intersection):
                linesegment_intersections.append(intersection)
        return linesegment_intersections


class FullArc3D(Arc3D):
    """
    An edge that starts at start_end, ends at the same point after having described a circle.

    """

    def __init__(self, center: volmdlr.Point3D, start_end: volmdlr.Point3D,
                 normal: volmdlr.Vector3D,
                 name: str = ''):
        self.__center = center
        self.__normal = normal
        self.start_end = start_end
        interior = start_end.rotation(center, normal, math.pi)
        Arc3D.__init__(self, start=start_end, end=start_end,
                       interior=interior, name=name)  # !!! this is dangerous

    def __hash__(self):
        return hash(self.center) + 5 * hash(self.start_end)

    def __eq__(self, other_arc):
        return (self.center == other_arc.center) \
            and (self.start == other_arc.start)

    @property
    def center(self):
        return self.__center

    @property
    def angle(self):
        return volmdlr.TWO_PI

    @property
    def normal(self):
        return self.__normal

    @property
    def is_trigo(self):
        return True

    def copy(self, *args, **kwargs):
        return FullArc3D(self._center.copy(), self.end.copy(), self._normal.copy())

    def to_dict(self, use_pointers: bool = False, memo=None, path: str = '#'):
        dict_ = self.base_dict()
        dict_['center'] = self.center.to_dict(use_pointers=use_pointers, memo=memo, path=path + '/center')
        dict_['radius'] = self.radius
        dict_['angle'] = self.angle
        dict_['is_trigo'] = self.is_trigo
        dict_['start_end'] = self.start.to_dict(use_pointers=use_pointers, memo=memo, path=path + '/start_end')
        dict_['normal'] = self.normal.to_dict(use_pointers=use_pointers, memo=memo, path=path + '/normal')
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
        center = self.center.to_2d(plane_origin, x, y)
        start_end = self.start.to_2d(plane_origin, x, y)
        return FullArc2D(center, start_end)

    def to_step(self, current_id, surface_id=None):
        """Exports to STEP format."""
        # Not calling Circle3D.to_step because of circular imports
        u = self.start - self.center
        u.normalize()
        v = self.normal.cross(u)
        frame = volmdlr.Frame3D(self.center, self.normal, u, v)
        content, frame_id = frame.to_step(current_id)
        curve_id = frame_id + 1
        # Not calling Circle3D.to_step because of circular imports
        content += f"#{curve_id} = CIRCLE('{self.name}',#{frame_id},{self.radius * 1000});\n"

        if surface_id:
            content += f"#{curve_id + 1} = SURFACE_CURVE('',#{curve_id},(#{surface_id}),.PCURVE_S1.);\n"
            curve_id += 1

        point1 = (self.center + u * self.radius).to_point()

        p1_content, p1_id = point1.to_step(curve_id + 1, vertex=True)
        content += p1_content

        edge_curve = p1_id + 1
        content += f"#{edge_curve} = EDGE_CURVE('{self.name}',#{p1_id},#{p1_id},#{curve_id},.T.);\n"
        curve_id += 1

        return content, [edge_curve]

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            ax = Axes3D(plt.figure())

        x = []
        y = []
        z = []
        for x_component, y_component, z_component in self.discretization_points(number_points=20):
            x.append(x_component)
            y.append(y_component)
            z.append(z_component)
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])
        ax.plot(x, y, z, color=edge_style.color, alpha=edge_style.alpha)

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
        new_start_end = self.start.rotation(center, axis, angle)
        new_center = self._center.rotation(center, axis, angle)
        new_normal = self._normal.rotation(center, axis, angle)
        return FullArc3D(new_center, new_start_end,
                         new_normal, name=self.name)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.start.rotation(center, axis, angle, False)
        self.end.rotation(center, axis, angle, False)
        self._center.rotation(center, axis, angle, False)
        self.interior.rotation(center, axis, angle, False)
        self._bbox = None

    def translation(self, offset: volmdlr.Vector3D):
        new_start_end = self.start.translation(offset, True)
        new_center = self._center.translation(offset, True)
        new_normal = self._normal.translation(offset, True)
        return FullArc3D(new_center, new_start_end,
                         new_normal, name=self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.start.translation(offset, False)
        self.end.translation(offset, False)
        self._center.translation(offset, False)
        self.interior.translation(offset, False)
        self._bbox = None

    def linesegment_intersections(self, linesegment3d: LineSegment3D):
        """
        Calculates the intersections between a full arc 3d and a line segment 3d.

        :param linesegment3d: linesegment 3d to verify intersections.
        :return: list of points 3d, if there are any intersections, an empty list if otherwise.
        """
        distance_center_lineseg = linesegment3d.point_distance(self.frame.origin)
        if distance_center_lineseg > self.radius:
            return []
        direction_vector = linesegment3d.direction_vector()
        if math.isclose(self.frame.w.dot(direction_vector), 0, abs_tol=1e-6) and \
                not math.isclose(linesegment3d.start.z - self.frame.origin.z, 0, abs_tol=1e-6):
            return []

        if linesegment3d.start.z == linesegment3d.end.z == self.frame.origin.z:
            quadratic_equation_a = 1 + (direction_vector.y ** 2 / direction_vector.x ** 2)
            quadratic_equation_b = (-2 * (direction_vector.y ** 2 / direction_vector.x ** 2) * linesegment3d.start.x +
                                    2 * (direction_vector.y / direction_vector.x) * linesegment3d.start.y)
            quadratic_equation_c = (linesegment3d.start.y - (direction_vector.y / direction_vector.x) *
                                    linesegment3d.start.x) ** 2 - self.radius ** 2
            delta = quadratic_equation_b ** 2 - 4 * quadratic_equation_a * quadratic_equation_c
            x1 = (- quadratic_equation_b + math.sqrt(delta)) / (2 * quadratic_equation_a)
            x2 = (- quadratic_equation_b - math.sqrt(delta)) / (2 * quadratic_equation_a)
            y1 = (direction_vector.y / direction_vector.x) * (x1 - linesegment3d.start.x) + linesegment3d.start.y
            y2 = (direction_vector.y / direction_vector.x) * (x2 - linesegment3d.start.x) + linesegment3d.start.y
            return [volmdlr.Point3D(x1, y1, self.frame.origin.z), volmdlr.Point3D(x2, y2, self.frame.origin.z)]
        constant = (self.frame.origin.z - linesegment3d.start.z) / direction_vector.z
        x_coordinate = constant * direction_vector.x + linesegment3d.start.x
        y_coordinate = constant * direction_vector.y + linesegment3d.start.y
        if math.isclose((x_coordinate - self.frame.origin.x) ** 2 + (y_coordinate - self.frame.origin.y) ** 2,
                        self.radius ** 2, abs_tol=1e-6):
            return [volmdlr.Point3D(x_coordinate, y_coordinate, self.frame.origin.z)]
        return []

    def reverse(self):
        """
        Defines a new FullArc3D, identical to self, but in the opposite direction.

        """
        return self


class ArcEllipse3D(Edge):
    """
    An arc is defined by a starting point, an end point and an interior point.

    """

    def __init__(self, start: volmdlr.Point3D, interior: volmdlr.Point3D, end: volmdlr.Point3D,
                 center: volmdlr.Point3D, major_dir: volmdlr.Vector3D, normal: volmdlr.Vector3D = None,
                 extra: volmdlr.Point3D = None, name=''):
        Edge.__init__(self, start=start, end=end, name=name)
        self.interior = interior
        self.center = center
        major_dir.normalize()
        self.major_dir = major_dir  # Vector for Gradius
        self.normal = normal
        self.extra = extra
        if not normal:
            u1 = self.interior - self.start
            u2 = self.interior - self.end
            u1.normalize()
            u2.normalize()

            if u1.is_close(u2):
                u2 = self.interior - self.extra
                u2.normalize()

            n = u2.cross(u1)
            n.normalize()
            self.normal = n

        self.minor_dir = self.normal.cross(self.major_dir)

        frame = volmdlr.Frame3D(self.center, self.major_dir, self.minor_dir, self.normal)
        self.frame = frame
        start_new, end_new = frame.global_to_local_coordinates(
            self.start), frame.global_to_local_coordinates(self.end)
        interior_new, center_new = frame.global_to_local_coordinates(
            self.interior), frame.global_to_local_coordinates(self.center)
        self._bbox = None

        # from :
        # https://math.stackexchange.com/questions/339126/how-to-draw-an-ellipse-if-a-center-and-3-arbitrary-points-on-it-are-given

        def theta_a_b(start_, iterior_, end_, center_):
            """
            center-and-3-arbitrary-points-on-it-are-given.

            theta= ellipse's inclination angle related to the horizontal
            (clockwise),a=semi major axis, B=semi minor axis.

            """
            x_start, y_start, x_interior, y_interior, x_end, y_end = start_[0] - center_[0], start_[1] - center_[1], \
                iterior_[0] - center_[0], iterior_[1] - center_[
                                                                         1], end_[0] - center_[0], end_[1] - center_[1]
            matrix_a = npy.array(([x_start ** 2, y_start ** 2, 2 * x_start * y_start],
                                  [x_interior ** 2, y_interior ** 2, 2 * x_interior * y_interior],
                                  [x_end ** 2, y_end ** 2, 2 * x_end * y_end]))
            inv_matrix_a = npy.linalg.inv(matrix_a)
            identity = npy.array(([1], [1], [1]))
            r1, r2, r3 = npy.dot(inv_matrix_a, identity)  # 3 item column matrix
            theta = 0.5 * math.atan(2 * r3 / (r2 - r1))
            c1 = r1 + r2
            c2 = (r2 - r1) / math.cos(2 * theta)
            major_axis = math.sqrt((2 / (c1 - c2)))
            minor_axis = math.sqrt((2 / (c1 + c2)))
            return theta, major_axis, minor_axis

        if start.is_close(end):
            extra_new = frame.global_to_local_coordinates(self.extra)
            theta, major_axis, minor_axis = theta_a_b(start_new, interior_new, extra_new, center_new)

        else:
            if not self.extra:
                theta, major_axis, minor_axis = theta_a_b(start_new, interior_new, end_new, center_new)
            else:
                extra_new = frame.global_to_local_coordinates(self.extra)
                theta, major_axis, minor_axis = theta_a_b(start_new, interior_new, extra_new, center_new)

        self.Gradius = major_axis
        self.Sradius = minor_axis
        self.theta = theta

        # Angle start
        start_u1, start_u2 = start_new.x / self.Gradius, start_new.y / self.Sradius
        # angle1 = volmdlr.geometry.sin_cos_angle(start_u1, start_u2)
        angle1 = math.atan2(start_u2, start_u1)
        self.angle_start = angle1
        # Angle end
        end_u3, end_u4 = end_new.x / self.Gradius, end_new.y / self.Sradius
        # angle2 = volmdlr.geometry.sin_cos_angle(end_u3, end_u4)
        angle2 = math.atan2(end_u4, end_u3)
        self.angle_end = angle2
        # Angle interior
        interior_u5, interior_u6 = interior_new.x / self.Gradius, interior_new.y / self.Sradius
        # anglei = volmdlr.geometry.sin_cos_angle(interior_u5, interior_u6)
        anglei = math.atan2(interior_u6, interior_u5)
        self.angle_interior = anglei
        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei + volmdlr.TWO_PI) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + volmdlr.TWO_PI

        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2 + volmdlr.TWO_PI) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + volmdlr.TWO_PI

        if clockwise_path > trigowise_path:
            self.is_trigo = True
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle = clockwise_path

        if self.start.is_close(self.end):
            self.angle = volmdlr.TWO_PI

        if self.is_trigo:
            self.offset_angle = angle1
        else:
            self.offset_angle = angle2

        volmdlr.core.CompositePrimitive3D.__init__(self,
                                                   primitives=self.discretization_points(number_points=20),
                                                   name=name)

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
        if angle_start > self.angle_interior > angle_end or angle_start < self.angle_interior < angle_end:
            angle_end = self.angle_end
            angle_start = self.angle_start
        elif self.angle_start == self.angle_end:
            angle_start = 0
            angle_end = 2 * math.pi
        else:
            if angle_end < angle_start:
                angle_end = self.angle_end + volmdlr.TWO_PI
            elif angle_start < angle_end:
                angle_end = self.angle_end - volmdlr.TWO_PI

        discretization_points = [self.frame.local_to_global_coordinates(
            volmdlr.Point3D(self.Gradius * math.cos(angle), self.Sradius * math.sin(angle), 0))
            for angle in npy.linspace(angle_start, angle_end, number_points)]
        return discretization_points

    def polygon_points(self, discretization_resolution: int):
        warnings.warn('polygon_points is deprecated,\
        please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(angle_resolution=discretization_resolution)

    def _get_points(self):
        return self.discretization_points(number_points=20)

    points = property(_get_points)

    def to_2d(self, plane_origin, x, y):
        """
        Transforms an Arc Ellipse 3D into an Arc Ellipse 2D, given a plane origin and an u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: ArcEllipse2D.
        """
        point_start2d = self.start.to_2d(plane_origin, x, y)
        point_interior2d = self.interior.to_2d(plane_origin, x, y)
        point_end2d = self.end.to_2d(plane_origin, x, y)
        center = self.center.to_2d(plane_origin, x, y)
        point_major_dir = self.center + self.Gradius * self.major_dir
        point_major_dir_2d = point_major_dir.to_2d(plane_origin, x, y)
        vector_major_dir_2d = point_major_dir_2d - center
        vector_major_dir_2d.normalize()
        extra = self.extra
        if extra:
            extra = self.extra.to_2d(plane_origin, x, y)
        return ArcEllipse2D(point_start2d, point_interior2d, point_end2d, center, vector_major_dir_2d, extra,
                            name=self.name)

    def length(self):
        """Computes the length."""
        return self.angle * math.sqrt(
            (self.Gradius ** 2 + self.Sradius ** 2) / 2)

    def normal_vector(self, abscissa):
        raise NotImplementedError

    def direction_vector(self, abscissa):
        raise NotImplementedError

    def abscissa(self, point: volmdlr.Point3D, tol: float = 1e-6):
        """
        Calculates the abscissa a given point.

        :param point: point to calculate abscissa.
        :return: abscissa
        """
        if point.point_distance(self.start) < tol:
            return 0
        vector_2 = self.normal.cross(self.major_dir)
        ellipse_2d = self.to_2d(self.center, self.major_dir, vector_2)
        point2d = point.to_2d(self.center, self.major_dir, vector_2)
        return ellipse_2d.abscissa(point2d)

    def reverse(self):
        """
        Reverse the Arc Ellipse 3D.

        :return:
        """
        normal = None
        extra = None
        if self.normal:
            normal = self.normal.copy()
        if self.extra:
            extra = self.extra.copy()
        return self.__class__(self.end.copy(),
                              self.interior.copy(),
                              self.start.copy(),
                              self.center.copy(),
                              self.major_dir.copy(),
                              normal,
                              extra,
                              self.name)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """Plot the arc ellipse."""
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        ax.plot([self.interior[0]], [self.interior[1]], [self.interior[2]],
                color='b')
        ax.plot([self.start[0]], [self.start[1]], [self.start[2]], c='r')
        ax.plot([self.end[0]], [self.end[1]], [self.end[2]], c='r')
        ax.plot([self.interior[0]], [self.interior[1]], [self.interior[2]],
                c='g')
        x = []
        y = []
        z = []
        for x_component, y_component, z_component in self.discretization_points(number_points=20):
            x.append(x_component)
            y.append(y_component)
            z.append(z_component)

        ax.plot(x, y, z, edge_style.color, alpha=edge_style.alpha)
        if edge_style.edge_ends:
            self.start.plot(ax)
            self.end.plot(ax)
        return ax

    def plot2d(self, x3d: volmdlr.Vector3D = volmdlr.X3D, y3d: volmdlr.Vector3D = volmdlr.Y3D,
               ax=None, color='k'):
        """
        Plot 2d for an arc ellipse 3d.

        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

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
        Triangulation for an ArcEllipse3D.

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

        :return: a volmdlr.core.BoundingBox object.
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

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new ArcEllipse3D.

        :param frame: Local coordinate system.
        :type frame: volmdlr.Frame3D
        :param side: 'old' will perform a transformation from local to global coordinates. 'new' will
            perform a tranformation from global to local coordinates.
        :type side: str
        :return: A new transformed ArcEllipse3D.
        :rtype: ArcEllipse3D
        """
        if side == 'old':
            return ArcEllipse3D(frame.local_to_global_coordinates(self.start),
                                frame.local_to_global_coordinates(self.interior),
                                frame.local_to_global_coordinates(self.end),
                                frame.local_to_global_coordinates(self.center),
                                self.major_dir)
        if side == 'new':
            point_major_dir = self.center + self.major_dir * self.major_axis
            major_dir = frame.global_to_local_coordinates(point_major_dir).to_vector()
            major_dir.normalize()
            return ArcEllipse3D(frame.global_to_local_coordinates(self.start),
                                frame.global_to_local_coordinates(self.interior),
                                frame.global_to_local_coordinates(self.end),
                                frame.global_to_local_coordinates(self.center),
                                major_dir)
        raise ValueError('Side should be \'new\' \'old\'')

    def point_belongs(self, point, abs_tol: float = 1e-6):
        """
        Verifies if a given point lies on the arc of ellipse 3D.

        :param point: point to be verified.
        :param abs_tol: Absolute tolerance to consider the point on the curve.
        :return: True is point lies on the arc of ellipse, False otherwise
        """
        vector_2 = self.normal.cross(self.major_dir)
        ellipse_2d = self.to_2d(self.center, self.major_dir, vector_2)
        point2d = point.to_2d(self.center, self.major_dir, vector_2)
        return ellipse_2d.point_belongs(point2d, abs_tol=abs_tol)


class FullArcEllipse3D(FullArcEllipse, ArcEllipse3D):
    """
    Defines a FullArcEllipse3D.
    """

    def __init__(self, start_end: volmdlr.Point3D, major_axis: float, minor_axis: float,
                 center: volmdlr.Point3D, normal: volmdlr.Vector3D, major_dir: volmdlr.Vector3D, name: str = ''):
        normal.normalize()
        self.normal = normal
        major_dir.normalize()
        self.minor_dir = normal.cross(major_dir)
        frame = volmdlr.Frame3D(center, major_dir, self.minor_dir, normal)
        self.frame = frame
        center2d = center.to_2d(center, major_dir, self.minor_dir)
        point_major_dir = center + major_axis * major_dir
        point_major_dir_2d = point_major_dir.to_2d(center, major_dir, self.minor_dir)
        vector_major_dir_2d = (point_major_dir_2d - center2d).to_vector()
        self.theta = volmdlr.geometry.clockwise_angle(vector_major_dir_2d, volmdlr.X2D)
        if self.theta == math.pi * 2:
            self.theta = 0.0
        self._bbox = None

        FullArcEllipse.__init__(self, start_end, major_axis, minor_axis, center, major_dir, name)

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = 20):
        """
        Discretize a Contour to have "n" points.

        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned.
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Useful
            to mesh an arc.
        :return: a list of sampled points.
        """
        if not number_points:
            number_points = math.ceil(volmdlr.TWO_PI * angle_resolution) + 2
        discretization_points_3d = [
                                       self.center + self.major_axis * math.cos(
                                           teta) * self.major_dir
                                       + self.minor_axis * math.sin(
                                           teta) * self.major_dir.cross(
                                           self.normal) for teta in
                                       npy.linspace(0, volmdlr.TWO_PI,
                                                    number_points)][:-1]
        return discretization_points_3d

    def to_2d(self, plane_origin, x, y):
        """
        Transforms a FullArcEllipse3D into an FullArcEllipse2D, given an plane origin and a u and v plane vector.

        :param plane_origin: plane origin.
        :param x: plane u vector.
        :param y: plane v vector.
        :return: FullArcEllipse2D.
        """
        point_start_end2d = self.start_end.to_2d(plane_origin, x, y)
        center2d = self.center.to_2d(plane_origin, x, y)
        point_major_dir = self.center + self.major_axis * self.major_dir
        point_major_dir_2d = point_major_dir.to_2d(plane_origin, x, y)
        vector_major_dir_2d = (point_major_dir_2d - center2d).to_vector()
        vector_major_dir_2d.normalize()
        return FullArcEllipse2D(point_start_end2d, self.major_axis, self.minor_axis, center2d,
                                vector_major_dir_2d, name=self.name)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new FullArcEllipse3D.

        :param frame: Local coordinate system.
        :type frame: volmdlr.Frame3D
        :param side: 'old' will perform a tranformation from local to global coordinates. 'new' will
            perform a tranformation from global to local coordinates.
        :type side: str
        :return: A new transformed FulLArcEllipse3D.
        :rtype: FullArcEllipse3D
        """
        if side == 'old':
            return FullArcEllipse3D(frame.local_to_global_coordinates(self.start_end),
                                    self.major_axis, self.minor_axis,
                                    frame.local_to_global_coordinates(self.center),
                                    frame.local_to_global_coordinates(self.normal), self.major_dir, self.name)
        if side == 'new':
            point_major_dir = self.center + self.major_dir * self.major_axis
            major_dir = frame.global_to_local_coordinates(point_major_dir).to_vector()
            major_dir.normalize()
            return FullArcEllipse3D(frame.global_to_local_coordinates(self.start_end),
                                    self.major_axis, self.minor_axis,
                                    frame.global_to_local_coordinates(self.center),
                                    frame.global_to_local_coordinates(self.normal), major_dir, self.name)
        raise ValueError('Side should be \'new\' \'old\'')

    def translation(self, offset: volmdlr.Vector3D):
        """
        FullArcEllipse3D translation.

        :param offset: translation vector.
        :type offset: volmdlr.Vector3D
        :return: A new translated FullArcEllipse3D.
        :rtype: FullArcEllipse3D
        """
        return FullArcEllipse3D(self.start_end.translation(offset), self.major_axis, self.minor_axis,
                                self.center.translation(offset), self.normal, self.major_dir, self.name)

    def abscissa(self, point: volmdlr.Point3D, tol: float = 1e-6):
        """
        Calculates the abscissa a given point.

        :param point: point to calculate abscissa.
        :return: abscissa
        """
        vector_2 = self.normal.cross(self.major_dir)
        ellipse_2d = self.to_2d(self.center, self.major_dir, vector_2)
        point2d = point.to_2d(self.center, self.major_dir, vector_2)
        return ellipse_2d.abscissa(point2d)

    def normal_vector(self, abscissa):
        """
        Calculates the normal vector the edge at given abscissa.

        :return: the normal vector
        """
        raise NotImplementedError

    def direction_vector(self, abscissa):
        """
        Calculates the direction vector the edge at given abscissa.

        :param abscissa: edge abscissa
        :return: direction vector
        """
        raise NotImplementedError
