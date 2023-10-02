#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extended primitives 2D classes.
"""

import math
from typing import List, Dict

import matplotlib.patches

import volmdlr
from volmdlr import wires, edges, curves
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
                 adapt_radius: bool = False, name: str = ''):
        RoundedLineSegments.__init__(self, points, radius,
                                     adapt_radius=adapt_radius,
                                     name='')

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

    def offset(self, offset):
        number_points = len(self.points)
        vectors = []
        for i in range(number_points - 1):
            v1 = self.points[i + 1] - self.points[i]
            v2 = self.points[i] - self.points[i + 1]
            v1 = v1.unit_vector()
            v2 = v2.unit_vector()
            vectors.append(v1)
            vectors.append(v2)

        if self.closed:
            v1 = self.points[0] - self.points[-1]
            v2 = self.points[-1] - self.points[0]
            v1 = v1.unit_vector()
            v2 = v2.unit_vector()
            vectors.append(v1)
            vectors.append(v2)

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

            normal_vector1 = - vectors[2 * i - 1].normal_vector()
            normal_vector2 = vectors[2 * i].normal_vector()
            normal_vector1 = normal_vector1.unit_vector()
            normal_vector2 = normal_vector2.unit_vector()
            alpha = math.acos(normal_vector1.dot(normal_vector2))

            offset_point = self.points[i] + offset / math.cos(alpha / 2) * \
                offset_vectors[i - (not self.closed)]
            offset_points.append(offset_point)

        if not self.closed:
            normal_1 = vectors[0].normal_vector()
            offset_vectors.insert(0, normal_1)
            offset_points.insert(0,
                                 self.points[0] + offset * offset_vectors[0])

            n_last = vectors[-1].normal_vector()
            n_last = - n_last
            offset_vectors.append(n_last)
            offset_points.append(self.points[-1] + offset * offset_vectors[-1])

        return self.__class__(offset_points, new_radii,
                              adapt_radius=self.adapt_radius)

    def offset_single_line(self, line_index, offset):
        """
        Offsets a single line.

        :param line_index: 0 being the 1st line
        """
        new_linesegment2d_points = []
        dont_add_last_point = False

        for i, point in enumerate(
                self.points[:-1] + (self.closed) * [self.points[-1]]):

            if i == line_index:
                # Not closed RLS2D and the offset line is the last one
                if i == len(self.points) - 2:
                    dir_vec_1 = volmdlr.Vector2D(point - self.points[i - 1])
                    dir_vec_1 = dir_vec_1.unit_vector()
                    dir_vec_2 = dir_vec_1
                    dont_add_last_point = True
                # The offset line is the first one
                elif i == 0:
                    dir_vec_2 = volmdlr.Vector2D(
                        self.points[i + 1] - self.points[i + 2])
                    dir_vec_2 = dir_vec_2.unit_vector()
                    if not self.closed:
                        dir_vec_1 = dir_vec_2
                    else:
                        dir_vec_1 = volmdlr.Vector2D(
                            point - self.points[i - 1])
                        dir_vec_1 = dir_vec_1.unit_vector()
                # Closed RLS2D and the offset line is the last one
                elif i == len(self.points) - 1:
                    dir_vec_1 = volmdlr.Vector2D(point - self.points[i - 1])
                    dir_vec_1 = dir_vec_1.unit_vector()
                    dir_vec_2 = volmdlr.Vector2D(
                        self.points[0] - self.points[1])
                    dir_vec_2 = dir_vec_2.unit_vector()
                    dont_add_last_point = True
                else:
                    dir_vec_1 = volmdlr.Vector2D(point - self.points[i - 1])
                    dir_vec_1 = dir_vec_1.unit_vector()
                    dir_vec_2 = volmdlr.Vector2D(
                        self.points[i + 1] - self.points[i + 2])
                    dir_vec_2 = dir_vec_2.unit_vector()

                if self.closed and line_index == len(self.points) - 1:
                    normal_vector = volmdlr.Vector2D(
                        self.points[0] - point).unit_normal_vector()
                else:
                    normal_vector = volmdlr.Vector2D(
                        self.points[i + 1] - point).unit_normal_vector()

                alpha1 = math.acos(dir_vec_1.dot(normal_vector))
                alpha2 = math.acos(dir_vec_2.dot(normal_vector))

                # If 3 segments are aligned and the middle one have to be offset
                if math.isclose(math.cos(alpha1), 0,
                                abs_tol=1e-9) or math.isclose(math.cos(alpha2),
                                                              0, abs_tol=1e-9):
                    return self
                #                    distance_dir1 = offset
                #                    distance_dir2 = offset

                distance_dir1 = offset / math.cos(alpha1)
                distance_dir2 = offset / math.cos(alpha2)

                new_point1 = point + distance_dir1 * dir_vec_1
                if self.closed and line_index == len(self.points) - 1:
                    new_point2 = self.points[0] + distance_dir2 * dir_vec_2
                else:
                    new_point2 = self.points[i + 1] + distance_dir2 * dir_vec_2

                new_linesegment2d_points.append(new_point1)
                new_linesegment2d_points.append(new_point2)

            elif i == line_index + 1:
                pass

            elif line_index == len(self.points) - 1 and i == 0:
                pass
            else:
                new_linesegment2d_points.append(point)

        if not dont_add_last_point and not self.closed:
            new_linesegment2d_points.append(self.points[-1])

        rls_2d = self.__class__(new_linesegment2d_points, self.radius,
                                adapt_radius=self.adapt_radius)

        return rls_2d

    def offset_lines(self, line_indexes, offset):
        """
        line_indexes is a list of consecutive line indexes.

        These line should all be aligned line_indexes = 0 being the 1st line.

        if self.close last line_index can be len(self.points)-1
        if not, last line_index can be len(self.points)-2
        """
        new_linesegment2d_points = []

        # =============================================================================
        # COMPUTES THE DIRECTIVE VECTORS BETWEEN WHICH THE OFFSET WILL BE DRAWN
        # =============================================================================
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

        # =============================================================================
        # COMPUTES THE ANGLE BETWEEN THE NORMAL VECTOR OF THE SURFACE TO OFFSET AND
        # THE DIRECTIVE VECTOR IN ORDER TO SET THE NEW POINT AT THE RIGHT DISTANCE
        # =============================================================================
        normal_vectors = []
        for index in line_indexes:
            if index == len(self.points) - 1:
                normal_vectors.append(volmdlr.Vector2D(
                    self.points[0] - self.points[index]).normalVector(
                    unit=True))
            else:
                normal_vectors.append(
                    (self.points[index + 1] - self.points[index]).unit_normal_vector())

        dot1 = dir_vec_1.dot(normal_vectors[0])
        dot2 = dir_vec_2.dot(normal_vectors[-1])

        if math.isclose(dot1, 0, abs_tol=1e-9):
            # call function considering the line before, because the latter and
            # the first offset segment are parallel
            return self.offset_lines([line_indexes[0] - 1] + line_indexes,
                                     offset)
        if math.isclose(dot2, 0, abs_tol=1e-9):
            # call function considering the line after, because the latter and
            # the last offset segment are parallel
            return self.offset_lines(line_indexes + [line_indexes[-1] + 1],
                                     offset)

        distance_dir1 = offset / dot1
        distance_dir2 = offset / dot2

        if len(line_indexes) > 1:
            intersection = volmdlr.Point2D.line_intersection(
                curves.Line2D(self.points[line_indexes[0]],
                              self.points[line_indexes[0]] + dir_vec_1),
                curves.Line2D(self.points[line_indexes[-1] + 1],
                              self.points[line_indexes[-1] + 1] + dir_vec_2))
            vec1 = intersection.point_distance(
                self.points[line_indexes[0]]) * dir_vec_1
            vec2 = intersection.point_distance(
                self.points[line_indexes[-1] + 1]) * dir_vec_2

        # =============================================================================
        # COMPUTES THE NEW POINTS AFTER THE OFFSET
        # =============================================================================
        new_points = {}

        new_points[line_indexes[0]] = self.points[line_indexes[
            0]] + distance_dir1 * dir_vec_1

        for i, index in enumerate(line_indexes[1:]):
            coeff_vec_2 = volmdlr.Point2D.point_distance(
                self.points[line_indexes[0]],
                self.points[index]) / volmdlr.Point2D.point_distance(
                self.points[line_indexes[0]],
                self.points[line_indexes[-1] + 1])
            coeff_vec_1 = 1 - coeff_vec_2
            if dir_vec_1.dot(normal_vectors[i + 1]) < 0:
                coeff_vec_1 = - coeff_vec_1
            if dir_vec_2.dot(normal_vectors[i + 1]) < 0:
                coeff_vec_2 = - coeff_vec_2
            index_dir_vector = coeff_vec_1 * vec1 + coeff_vec_2 * vec2
            index_dot = index_dir_vector.dot(normal_vectors[i + 1])
            index_distance_dir = offset / index_dot
            new_points[index] = self.points[
                                    index] + index_distance_dir * index_dir_vector

        if self.closed and line_indexes[-1] == len(self.points) - 1:
            new_points[0] = self.points[0] + distance_dir2 * dir_vec_2
        else:
            new_points[line_indexes[-1] + 1] = self.points[line_indexes[
                                                               -1] + 1] + distance_dir2 * dir_vec_2

        # =============================================================================
        # CREATE THE NEW POINTS' LIST
        # =============================================================================
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

    def __init__(self, points: List[volmdlr.Point2D], radius: Dict[int, float],
                 adapt_radius: bool = False, name: str = ''):
        RoundedLineSegments2D.__init__(self, points, radius, adapt_radius=adapt_radius, name='')
        self.closed = False
        wires.Wire2D.__init__(self, self._primitives(), name)


class ClosedRoundedLineSegments2D(RoundedLineSegments2D, wires.Contour2D):
    """
    Defines a polygon with some rounded corners.

    :param points: Points used to draw the wire
    :type points: List of Point2D
    :param radius: Radius used to connect different parts of the wire
    :type radius: {position1(n): float which is the radius linked the n-1 and n+1 points, position2(n+1):...}
    """

    def __init__(self, points: List[volmdlr.Point2D], radius: Dict[int, float],
                 adapt_radius: bool = False, name: str = ''):
        RoundedLineSegments2D.__init__(self, points, radius, adapt_radius=adapt_radius, name='')
        self.closed = True
        # RoundedLineSegments.__init__(self, points, radius, closed=True, adapt_radius=adapt_radius, name=name)
        wires.Contour2D.__init__(self, self._primitives(), name)

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
