#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from typing import List
import math
import numpy as npy
import matplotlib.patches
import volmdlr
# from volmdlr.core_compiled import polygon_point_belongs
from volmdlr.primitives import RoundedLineSegments
import volmdlr.edges
import volmdlr.wires
import matplotlib.pyplot as plt


class OpenedRoundedLineSegments2D(RoundedLineSegments, volmdlr.wires.Wire2D):
    closed = False

    def __init__(self, points, radius, adapt_radius=False, name=''):
        primitives = RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.edges.LineSegment2D,
                                                  volmdlr.edges.Arc2D,
                                                  closed=False,
                                                  adapt_radius=adapt_radius,
                                                  name='')

        volmdlr.wires.Wire2D.__init__(self, primitives, name)

    def polygon_points(self, angle_resolution=5):
        points = []
        for primitive in self.primitives:
            points.extend(primitive.polygon_points())
        return points

    def arc_features(self, ipoint):
        radius = self.radius[ipoint]
        if self.closed:
            if ipoint == 0:
                pt1 = self.points[-1]
            else:
                pt1 = self.points[ipoint - 1]
            pti = self.points[ipoint]
            if ipoint < self.npoints - 1:
                pt2 = self.points[ipoint + 1]
            else:
                pt2 = self.points[0]
        else:
            pt1 = self.points[ipoint - 1]
            pti = self.points[ipoint]
            pt2 = self.points[ipoint + 1]

        # TODO: change to point_distance
        dist1 = (pt1 - pti).norm()
        dist2 = (pt2 - pti).norm()
        dist3 = (pt1 - pt2).norm()
        alpha = math.acos(
            -(dist3 ** 2 - dist1 ** 2 - dist2 ** 2) / (2 * dist1 * dist2)) / 2.
        dist = radius / math.tan(alpha)

        u1 = (pt1 - pti) / dist1
        u2 = (pt2 - pti) / dist2

        p3 = pti + u1 * dist
        p4 = pti + u2 * dist

        w = (u1 + u2)
        if w != volmdlr.Vector2D(0, 0):
            w.normalize()

        v1 = u1.deterministic_unit_normal_vector()
        if v1.dot(w) < 0:
            v1 = -v1

        pc = p3 + v1 * radius
        pm = pc - radius * w

        return p3, pm, p4, dist, alpha

    def rotation(self, center, angle, copy=True):
        if copy:
            return self.__class__([p.rotation(center, angle, copy=True) \
                                   for p in self.points],
                                  self.radius,
                                  adapt_radius=self.adapt_radius,
                                  name=self.name)
        else:
            self.__init__(
                [p.rotation(center, angle, copy=True) for p in self.points],
                self.radius,
                adapt_radius=self.adapt_radius, name=self.name)

    def translation(self, offset, copy=True):
        if copy:
            return self.__class__(
                [p.translation(offset, copy=True) for p in self.points],
                self.radius, adapt_radius=self.adapt_radius, name=self.name)
        else:
            self.__init__(
                [p.translation(offset, copy=True) for p in self.points],
                self.radius, adapt_radius=self.adapt_radius, name=self.name)

    def offset(self, offset):
        nb = len(self.points)
        vectors = []
        for i in range(nb - 1):
            v1 = self.points[i + 1] - self.points[i]
            v2 = self.points[i] - self.points[i + 1]
            v1.normalize()
            v2.normalize()
            vectors.append(v1)
            vectors.append(v2)

        if self.closed:
            v1 = self.points[0] - self.points[-1]
            v2 = self.points[-1] - self.points[0]
            v1.normalize()
            v2.normalize()
            vectors.append(v1)
            vectors.append(v2)

        offset_vectors = []
        new_radii = {}
        offset_points = []

        for i in range((not self.closed), nb - (not self.closed)):

            check = False
            ni = vectors[2 * i - 1] + vectors[2 * i]
            if ni == volmdlr.Vector2D(0, 0):
                ni = vectors[2 * i]
                ni = ni.normalVector()
                offset_vectors.append(ni)
            else:
                ni.normalize()
                if ni.dot(vectors[2 * i - 1].normal_vector()) > 0:
                    ni = - ni
                    check = True
                offset_vectors.append(ni)

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
            normal_vector1.normalize()
            normal_vector2.normalize()
            alpha = math.acos(normal_vector1.dot(normal_vector2))

            offset_point = self.points[i] + offset / math.cos(alpha / 2) * \
                           offset_vectors[i - (not self.closed)]
            offset_points.append(offset_point)

        if not self.closed:
            n1 = vectors[0].normalVector(unit=True)
            offset_vectors.insert(0, n1)
            offset_points.insert(0,
                                 self.points[0] + offset * offset_vectors[0])

            n_last = vectors[-1].normalVector(unit=True)
            n_last = - n_last
            offset_vectors.append(n_last)
            offset_points.append(self.points[-1] + offset * offset_vectors[-1])

        return self.__class__(offset_points, new_radii,
                              adapt_radius=self.adapt_radius)

    def offset_single_line(self, line_index, offset):
        """
        line_index = 0 being the 1st line
        """
        new_linesegment2D_points = []
        dont_add_last_point = False

        for i, point in enumerate(
                self.points[:-1] + (self.closed) * [self.points[-1]]):

            if i == line_index:
                # Not closed RLS2D and the offset line is the last one
                if i == len(self.points) - 2:
                    dir_vec_1 = volmdlr.Vector2D(point - self.points[i - 1])
                    dir_vec_1.normalize()
                    dir_vec_2 = dir_vec_1
                    dont_add_last_point = True
                # The offset line is the first one
                elif i == 0:
                    dir_vec_2 = volmdlr.Vector2D(
                        self.points[i + 1] - self.points[i + 2])
                    dir_vec_2.normalize()
                    if not self.closed:
                        dir_vec_1 = dir_vec_2
                    else:
                        dir_vec_1 = volmdlr.Vector2D(
                            point - self.points[i - 1])
                        dir_vec_1.normalize()
                # Closed RLS2D and the offset line is the last one
                elif i == len(self.points) - 1:
                    dir_vec_1 = volmdlr.Vector2D(point - self.points[i - 1])
                    dir_vec_1.normalize()
                    dir_vec_2 = volmdlr.Vector2D(
                        self.points[0] - self.points[1])
                    dir_vec_2.normalize()
                    dont_add_last_point = True
                else:
                    dir_vec_1 = volmdlr.Vector2D(point - self.points[i - 1])
                    dir_vec_1.normalize()
                    dir_vec_2 = volmdlr.Vector2D(
                        self.points[i + 1] - self.points[i + 2])
                    dir_vec_2.normalize()

                if self.closed and line_index == len(self.points) - 1:
                    normal_vector = volmdlr.Vector2D(
                        self.points[0] - point).normalVector(unit=True)
                else:
                    normal_vector = volmdlr.Vector2D(
                        self.points[i + 1] - point).normalVector(unit=True)

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

                new_linesegment2D_points.append(new_point1)
                new_linesegment2D_points.append(new_point2)

            elif i == line_index + 1:
                pass

            elif line_index == len(self.points) - 1 and i == 0:
                pass
            else:
                new_linesegment2D_points.append(point)

        if not dont_add_last_point and not self.closed:
            new_linesegment2D_points.append(self.points[-1])

        rls2D = self.__class__(new_linesegment2D_points, self.radius,
                               self.closed, adapt_radius=self.adapt_radius)

        return rls2D

    def offset_lines(self, line_indexes, offset):
        """
        line_indexes is a list of consecutive line indexes
        These line should all be aligned
        line_indexes = 0 being the 1st line

        if self.close last line_index can be len(self.points)-1
        if not, last line_index can be len(self.points)-2
        """
        new_linesegment2D_points = []

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

        dir_vec_1.normalize()
        dir_vec_2.normalize()

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
                volmdlr.edges.Line2D(self.points[line_indexes[0]],
                               self.points[line_indexes[0]] + dir_vec_1),
                volmdlr.edges.Line2D(self.points[line_indexes[-1] + 1],
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

        for nb, index in enumerate(line_indexes[1:]):
            coeff_vec_2 = volmdlr.Point2D.point_distance(
                self.points[line_indexes[0]],
                self.points[index]) / volmdlr.Point2D.point_distance(
                self.points[line_indexes[0]],
                self.points[line_indexes[-1] + 1])
            coeff_vec_1 = 1 - coeff_vec_2
            if dir_vec_1.dot(normal_vectors[nb + 1]) < 0:
                coeff_vec_1 = - coeff_vec_1
            if dir_vec_2.dot(normal_vectors[nb + 1]) < 0:
                coeff_vec_2 = - coeff_vec_2
            index_dir_vector = coeff_vec_1 * vec1 + coeff_vec_2 * vec2
            index_dot = index_dir_vector.dot(normal_vectors[nb + 1])
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
        for i in range(len(self.points)):
            if i in new_points.keys():
                new_linesegment2D_points.append(new_points[i])
            else:
                new_linesegment2D_points.append(self.points[i])

        rls2D = self.__class__(new_linesegment2D_points, self.radius,
                               adapt_radius=self.adapt_radius)

        return rls2D


    
class ClosedRoundedLineSegments2D(OpenedRoundedLineSegments2D,
                                  volmdlr.wires.Contour2D):
    """
    :param points: Points used to draw the wire 
    :type points: List of Point2D
    :param radius: Radius used to connect different parts of the wire
    :type radius: {position1(n): float which is the radius linked the n-1 and n+1 points, position2(n+1):...}
    """
    closed = True
    def __init__(self, points, radius, adapt_radius=False, name=''):
        primitives = RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.edges.LineSegment2D,
                                                  volmdlr.edges.Arc2D,
                                                  closed=True,
                                                  adapt_radius=adapt_radius, name='')

        volmdlr.wires.Contour2D.__init__(self, primitives, name)

class Measure2D(volmdlr.edges.LineSegment2D):
    def __init__(self, point1, point2, label='', unit='mm', type_='distance'):
        """
        :param unit: 'mm', 'm' or None. If None, the distance won't be in the label

        """
        # TODO: offset parameter
        volmdlr.edges.LineSegment2D.__init__(self, point1, point2)
        self.label = label
        self.unit = unit
        self.type_ = type_

    def plot(self, ax, ndigits=6):
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
            arrow = matplotlib.patches.FancyArrowPatch((x1, y1), (x2, y2),
                                    arrowstyle='<|-|>,head_length=10,head_width=5',
                                    shrinkA=0, shrinkB=0,
                                    color='k')
        elif self.type_ == 'radius':
            arrow = matplotlib.patches.FancyArrowPatch((x1, y1), (x2, y2),
                                    arrowstyle='-|>,head_length=10,head_width=5',
                                    shrinkA=0, shrinkB=0,
                                    color='k')

        ax.add_patch(arrow)
        if x2 - x1 == 0.:
            theta = 90.
        else:
            theta = math.degrees(math.atan((y2 - y1) / (x2 - x1)))
        ax.text(xm, ym, label, va='bottom', ha='center', rotation=theta)


