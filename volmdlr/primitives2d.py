#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extended primitives 2D classes.
"""

import math
from typing import List, Dict

import matplotlib.patches

import volmdlr
from volmdlr import wires, edges, curves, PATH_ROOT
from volmdlr.core import EdgeStyle
from volmdlr.primitives import RoundedLineSegments


class RoundedLineSegments2D(RoundedLineSegments):
    """
    A class representing a series of rounded line segments in 2D.

    This class inherits from the `RoundedLineSegments` class,
    and provides methods to work with rounded line segments in 2D.

    :param points: The list of points defining the line segments.
    :type points: List[volmdlr.Point2D]
    :param radius: The dictionary mapping segment indices to their respective radii.
    :type radius:Dict[int, float]
    :param adapt_radius: Flag indicating whether to adapt the radius based on segment length.
    Defaults to False.
    :type adapt_radius: bool, optional
    :param name: The name of the rounded line segments. Defaults to ''.
    :type name: str, optional
    """
    line_class = volmdlr.edges.LineSegment2D
    arc_class = volmdlr.edges.Arc2D

    def __init__(self, points: List[volmdlr.Point2D], radius: Dict[int, float],
                 adapt_radius: bool = False, reference_path: str = PATH_ROOT, name: str = ''):
        RoundedLineSegments.__init__(self, points=points, radius=radius, adapt_radius=adapt_radius,
                                     reference_path=reference_path, name=name)

    def arc_features(self, point_index: int):
        """
        Returns the arc features for point at index.
        """
        # raise NotImplementedError
        radius = self.radius[point_index]
        pt1, pti, pt2 = self.get_points(point_index)
        # TODO: change to point_distance ------> done
        point_distance1 = (pt1 - pti).norm()
        point_distance2 = (pt2 - pti).norm()
        point_distance3 = (pt1 - pt2).norm()
        alpha = math.acos(
            -(point_distance3 ** 2 - point_distance1 ** 2 - point_distance2 ** 2) / (2 * point_distance1
                                                                                     * point_distance2)) / 2.
        point_distance = radius / math.tan(alpha)

        u1 = (pt1 - pti) / point_distance1
        u2 = (pt2 - pti) / point_distance2

        p3 = pti + u1 * point_distance
        p4 = pti + u2 * point_distance

        w = (u1 + u2).to_vector()
        if not w.is_close(volmdlr.Vector2D(0, 0)):
            w = w.unit_vector()

        v1 = u1.deterministic_unit_normal_vector()
        if v1.dot(w) < 0:
            v1 = -v1

        point_curvature = p3 + v1 * radius
        point_interior = point_curvature - radius * w

        return p3, point_interior, p4, point_distance, alpha

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        OpenedRoundedLineSegments2D rotation.

        :param center: rotation center
        :param angle: angle rotation
        :return: a new rotated OpenedRoundedLineSegments2D
        """
        return self.__class__([point.rotation(center, angle)
                               for point in self.points],
                              self.radius,
                              adapt_radius=self.adapt_radius,
                              name=self.name)

    def translation(self, offset: volmdlr.Vector2D):
        """
        OpenedRoundedLineSegments2D translation.

        :param offset: translation vector
        :return: A new translated OpenedRoundedLineSegments2D
        """
        return self.__class__(
            [point.translation(offset) for point in self.points],
            self.radius, adapt_radius=self.adapt_radius, name=self.name)

    def _helper_offset_points_and_radii(self, vectors, number_points, offset):
        """
        Helper method to get offset points and new radii for offset method.

        :param vectors: offset vectors.
        :param number_points: number of points.
        :param offset: offset.
        :return: list of offset points and radii.
        """
        offset_vectors = []
        new_radii = {}
        offset_points = []

        for i in range((not self.closed), number_points - (not self.closed)):
            check = False
            normal_i = vectors[2 * i - 1] + vectors[2 * i]
            if normal_i.is_close(volmdlr.Vector2D(0, 0)):
                normal_i = vectors[2 * i]
                normal_i = normal_i.normal_vector()
                offset_vectors.append(normal_i)
            else:
                normal_i = normal_i.unit_vector()
                if normal_i.dot(vectors[2 * i - 1].normal_vector()) > 0:
                    normal_i = - normal_i
                    check = True
                offset_vectors.append(normal_i)

            if i in self.radius:
                if (check and offset > 0) or (not check and offset < 0):
                    new_radius = self.radius[i] + abs(offset)
                else:
                    new_radius = self.radius[i] - abs(offset)
                if new_radius > 0:
                    new_radii[i] = new_radius
                else:
                    if self.adapt_radius:
                        new_radii[i] = 1e-6

            normal_vector1 = - vectors[2 * i - 1].normal_vector().unit_vector()
            normal_vector2 = vectors[2 * i].normal_vector().unit_vector()
            alpha = math.acos(normal_vector1.dot(normal_vector2))

            offset_point = self.points[i] + offset / math.cos(alpha / 2) * offset_vectors[i - (not self.closed)]
            offset_points.append(offset_point)
        return offset_points, new_radii

    def offset(self, offset):
        """
        Return a new rounded line segment with the specified offset.

        This method creates a new rounded line segment by offsetting the current one by a given distance.
        The offset can be both positive and negative, moving the line segments outward or inward.

        :param offset: The offset distance for the new rounded line segment.
        :type offset: float

        :return: A new RoundedLineSegments2D instance with the specified offset.
        :rtype: RoundedLineSegments2D
        """
        number_points = len(self.points)
        vectors = []
        for i in range(number_points - 1):
            v1 = (self.points[i + 1] - self.points[i]).unit_vector()
            v2 = (self.points[i] - self.points[i + 1]).unit_vector()
            vectors.extend([v1, v2])

        if self.closed:
            v1 = (self.points[0] - self.points[-1]).unit_vector()
            v2 = (self.points[-1] - self.points[0]).unit_vector()
            vectors.extend([v1, v2])
        offset_points, new_radii = self._helper_offset_points_and_radii(vectors, number_points, offset)

        if not self.closed:
            normal_1 = vectors[0].normal_vector()
            offset_points.insert(0, self.points[0] + offset * normal_1)

            n_last = vectors[-1].normal_vector()
            n_last = - n_last
            offset_points.append(self.points[-1] + offset * n_last)

        return self.__class__(offset_points, new_radii, adapt_radius=self.adapt_radius)

    def _offset_directive_vector_helper(self, line_indexes):
        """
        Computes the directive vectors between which the offset will be drawn.

        :param line_indexes: A list of consecutive line indexes.
        :return: directive vector 1 and 2.
        """
        dir_vec_1 = None
        dir_vec_2 = None

        if line_indexes[0] == 0 and not self.closed:
            pass
        else:
            dir_vec_1 = self.points[line_indexes[0]] - self.points[line_indexes[0] - 1]

        if line_indexes[-1] == len(self.points) - (2 - self.closed):
            if not self.closed:
                pass
            else:
                dir_vec_2 = volmdlr.Vector2D((self.points[0] - self.points[1]))
        elif self.closed and line_indexes[-1] == len(self.points) - 2:
            dir_vec_2 = volmdlr.Vector2D(
                (self.points[line_indexes[-1] + 1] - self.points[0]))
        else:
            dir_vec_2 = self.points[line_indexes[-1] + 1] - self.points[line_indexes[-1] + 2]

        if dir_vec_1 is None:
            dir_vec_1 = dir_vec_2
        if dir_vec_2 is None:
            dir_vec_2 = dir_vec_1

        dir_vec_1 = dir_vec_1.unit_vector()
        dir_vec_2 = dir_vec_2.unit_vector()
        return dir_vec_1, dir_vec_2

    def get_offset_normal_vectors(self, line_indexes):
        """
        Gets offset normal vectors.

        :param line_indexes:
        :return: A list of consecutive line indexes.
        """
        normal_vectors = []
        for index in line_indexes:
            if index == len(self.points) - 1:
                normal_vectors.append(volmdlr.Vector2D(
                    self.points[0] - self.points[index]).normalVector(
                    unit=True))
            else:
                normal_vectors.append(
                    (self.points[index + 1] - self.points[index]).unit_normal_vector())
        return normal_vectors

    def _helper_get_offset_vectors(self, line_indexes, directive_vector1, directive_vector2):
        """
        Helper method: get offset vectors.

        """
        intersection = volmdlr.Point2D.line_intersection(
            curves.Line2D(self.points[line_indexes[0]], self.points[line_indexes[0]] + directive_vector1),
            curves.Line2D(self.points[line_indexes[-1] + 1], self.points[line_indexes[-1] + 1] + directive_vector2))
        vec1 = intersection.point_distance(self.points[line_indexes[0]]) * directive_vector1
        vec2 = intersection.point_distance(self.points[line_indexes[-1] + 1]) * directive_vector2
        return vec1, vec2

    def get_offset_new_points(self, line_indexes, offset, distance_dir1, distance_dir2, directive_vector1,
                              directive_vector2, normal_vectors):
        """
        Get Offset new points.

        """
        if len(line_indexes) <= 1:
            return []
        vec1, vec2 = self._helper_get_offset_vectors(line_indexes, directive_vector1, directive_vector2)
        new_points = {line_indexes[0]: self.points[line_indexes[0]] + distance_dir1 * directive_vector1}

        for i, index in enumerate(line_indexes[1:]):
            coeff_vec_2 = volmdlr.Point2D.point_distance(
                self.points[line_indexes[0]], self.points[index]) / volmdlr.Point2D.point_distance(
                self.points[line_indexes[0]], self.points[line_indexes[-1] + 1])
            coeff_vec_1 = 1 - coeff_vec_2
            if directive_vector1.dot(normal_vectors[i + 1]) < 0:
                coeff_vec_1 = - coeff_vec_1
            if directive_vector2.dot(normal_vectors[i + 1]) < 0:
                coeff_vec_2 = - coeff_vec_2
            index_dir_vector = coeff_vec_1 * vec1 + coeff_vec_2 * vec2
            index_dot = index_dir_vector.dot(normal_vectors[i + 1])
            new_points[index] = self.points[index] + (offset / index_dot) * index_dir_vector
        if self.closed and line_indexes[-1] == len(self.points) - 1:
            new_points[0] = self.points[0] + distance_dir2 * directive_vector2
        else:
            new_points[line_indexes[-1] + 1] = self.points[line_indexes[-1] + 1] + distance_dir2 * directive_vector2
        return new_points

    def offset_lines(self, line_indexes, offset):
        """
        Line indexes is a list of consecutive line indexes.

        These line should all be aligned line_indexes = 0 being the 1st line.

        if self.close last line_index can be len(self.points)-1
        if not, last line_index can be len(self.points)-2
        """
        new_linesegment2d_points = []

        dir_vec_1, dir_vec_2 = self._offset_directive_vector_helper(line_indexes)

        normal_vectors = self.get_offset_normal_vectors(line_indexes)

        # =============================================================================
        # COMPUTES THE ANGLE BETWEEN THE NORMAL VECTOR OF THE SURFACE TO OFFSET AND
        # THE DIRECTIVE VECTOR IN ORDER TO SET THE NEW POINT AT THE RIGHT DISTANCE
        # =============================================================================
        dot1 = dir_vec_1.dot(normal_vectors[0])
        dot2 = dir_vec_2.dot(normal_vectors[-1])

        if math.isclose(dot1, 0, abs_tol=1e-9):
            # call function considering the line before, because the latter and
            # the first offset segment are parallel
            return self.offset_lines([line_indexes[0] - 1] + line_indexes, offset)
        if math.isclose(dot2, 0, abs_tol=1e-9):
            # call function considering the line after, because the latter and
            # the last offset segment are parallel
            return self.offset_lines(line_indexes + [line_indexes[-1] + 1], offset)

        distance_dir1 = offset / dot1
        distance_dir2 = offset / dot2

        new_points = self.get_offset_new_points(line_indexes, offset, distance_dir1, distance_dir2,
                                                dir_vec_1, dir_vec_2, normal_vectors)
        for i, point in enumerate(self.points):
            if i in new_points:
                new_linesegment2d_points.append(new_points[i])
            else:
                new_linesegment2d_points.append(point)

        rls_2d = self.__class__(new_linesegment2d_points, self.radius,
                                adapt_radius=self.adapt_radius)

        return rls_2d


class OpenedRoundedLineSegments2D(RoundedLineSegments2D, wires.Wire2D):
    """
    Opened Rounded LineSegment2D class.

    :param points: Points used to draw the wire.
    :type points: List of Point2D.
    :param radius: Radius used to connect different parts of the wire.
    :type radius: {position1(n): float which is the radius linked the n-1 and n+1 points, position2(n+1):...}.
    """

    def __init__(self, points: List[volmdlr.Point2D], radius: Dict[int, float], adapt_radius: bool = False,
                 reference_path: str = PATH_ROOT, name: str = ''):
        RoundedLineSegments2D.__init__(self, points, radius, adapt_radius=adapt_radius,
                                       reference_path=reference_path, name='')
        self.closed = False
        wires.Wire2D.__init__(self, self._primitives(), reference_path=reference_path, name=name)


class ClosedRoundedLineSegments2D(RoundedLineSegments2D, wires.Contour2D):
    """
    Defines a polygon with some rounded corners.

    :param points: Points used to draw the wire
    :type points: List of Point2D
    :param radius: Radius used to connect different parts of the wire
    :type radius: {position1(n): float which is the radius linked the n-1 and n+1 points, position2(n+1):...}
    """

    def __init__(self, points: List[volmdlr.Point2D], radius: Dict[int, float],
                 adapt_radius: bool = False, reference_path: str = PATH_ROOT, name: str = ''):
        RoundedLineSegments2D.__init__(self, points, radius, adapt_radius=adapt_radius,
                                       reference_path=reference_path, name='')
        self.closed = True
        # RoundedLineSegments.__init__(self, points, radius, closed=True, adapt_radius=adapt_radius, name=name)
        wires.Contour2D.__init__(self, self._primitives(), reference_path=reference_path, name=name)

    def copy(self, deep=True, memo=None):
        """Returns a copy of the object."""
        return self.__class__([point.copy(deep, memo) for point in self.points], self.radius.copy(),
                              self.adapt_radius, name='copy_' + self.name)


class Measure2D(edges.LineSegment2D):
    """
    Measure 2D class.

    :param unit: 'mm', 'm' or None. If None, the distance won't be in the label.
    """

    def __init__(self, point1, point2, label='', unit: str = 'mm', type_: str = 'distance'):
        """
        :param unit: 'mm', 'm' or None. If None, the distance won't be in the label.

        """
        # TODO: offset parameter
        edges.LineSegment2D.__init__(self, point1, point2)
        self.label = label
        self.unit = unit
        self.type_ = type_

    def plot(self, ax, edge_style: EdgeStyle()):
        """Plots the Measure2D."""
        ndigits = 6
        x1, y1 = self.start
        x2, y2 = self.end
        x_middle, y_middle = 0.5 * (self.start + self.end)
        distance = self.end.point_distance(self.start)

        if self.label != '':
            label = f'{self.label}: '
        else:
            label = ''
        if self.unit == 'mm':
            label += f'{round(distance * 1000, ndigits)} mm'
        else:
            label += f'{round(distance, ndigits)} m'

        if self.type_ == 'distance':
            arrow = matplotlib.patches.FancyArrowPatch((x1, y1), (x2, y2),
                                                       arrowstyle='<|-|>,head_length=10,head_width=5',
                                                       shrinkA=0, shrinkB=0,
                                                       color=edge_style.color)
        elif self.type_ == 'radius':
            arrow = matplotlib.patches.FancyArrowPatch((x1, y1), (x2, y2),
                                                       arrowstyle='-|>,head_length=10,head_width=5',
                                                       shrinkA=0, shrinkB=0,
                                                       color=edge_style.color)

        ax.add_patch(arrow)
        if x2 - x1 == 0.:
            theta = 90.
        else:
            theta = math.degrees(math.atan((y2 - y1) / (x2 - x1)))
        ax.text(x_middle, y_middle, label, va='bottom', ha='center', rotation=theta)
