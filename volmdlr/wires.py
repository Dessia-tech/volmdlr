#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing wires & contours
"""

import warnings
import itertools
import math
from statistics import mean
from typing import List
import numpy as npy
from scipy.spatial import Delaunay, ConvexHull
import matplotlib.pyplot as plt
import matplotlib.patches
from mpl_toolkits.mplot3d import Axes3D
import plot_data.core as plot_data

import volmdlr
from volmdlr.core_compiled import polygon_point_belongs
import volmdlr.core
import volmdlr.edges
import volmdlr.display as vmd


def bounding_rectangle_adjacent_contours(contours: List):
    '''
    compute the bounding_box of a list of adjacent contours2d

    Parameters
    ----------
    contours : List[volmdlr.wires.Contour2D]

    Returns
    -------
    xmin : float
    xmax : float
    ymin : float
    ymax : float

    '''
    xmin, xmax, ymin, ymax = contours[0].bounding_rectangle()

    for i in range(1, len(contours)):
        xmin_contour, xmax_contour, ymin_contour, ymax_contour = contours[i].bounding_rectangle()
        xmin = min(xmin, xmin_contour)
        xmax = max(xmax, xmax_contour)
        ymin = min(ymin, ymin_contour)
        ymax = max(ymax, ymax_contour)

    return xmin, xmax, ymin, ymax


class Wire:
    _non_data_hash_attributes = ['basis_primitives']
    _non_serializable_attributes = ['primitive_to_index',
                                    'basis_primitives']

    # def __init__(self):
    #     raise TypeError ('It cannot be instantiated directly, see Wire2D, Wire3D, Contour2D or Contour3D')

    def length(self):
        length = 0.
        for primitive in self.primitives:
            length += primitive.length()
        return length

    def discretization_points(self, resolution: float):
        '''
        resolution: distance between two discretized points
        '''

        length = self.length()
        n = int(length / resolution) + 1
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

        if curvilinear_abscissa < length + 1e-9:
            return self.primitives[-1].end
        raise ValueError(
            'abscissa over length: {}>{}'.format(curvilinear_abscissa, length))

    def extract_primitives(self, point1, primitive1, point2, primitive2,
                           inside: bool = True):
        """
        inside: extracted contour is between the two points if True and outside
        these points if False
        """
        primitives = []

        # TODO: Check if it is: self.primitive_to_index[primitive1] OR self.primitive_to_index(primitive1)
        ip1 = self.primitive_to_index(primitive1)
        ip2 = self.primitive_to_index(primitive2)

        if inside:
            if ip1 < ip2:
                pass
            elif ip1 == ip2:  # primitive1 == primitive2
                if point1.point_distance(
                            primitive1.start) < point2.point_distance(
                        primitive1.start):
                    pass
                else:
                    primitive1, primitive2 = primitive2, primitive1
                    point1, point2 = point2, point1

            else:
                primitive1, primitive2 = primitive2, primitive1
                point1, point2 = point2, point1
        else:
            if ip1 > ip2:
                pass
            elif ip1 == ip2:  # primitive1 == primitive2
                if point1.point_distance(
                            primitive1.start) > point2.point_distance(
                        primitive1.start):
                    pass
                else:
                    primitive1, primitive2 = primitive2, primitive1
                    point1, point2 = point2, point1
            else:
                primitive1, primitive2 = primitive2, primitive1
                point1, point2 = point2, point1

        if ip1 < ip2:
            primitives.append(primitive1.split(point1)[1])
            primitives.extend(self.primitives[ip1 + 1:ip2])
            primitives.append(primitive2.split(point2)[0])
        else:
            primitives.append(primitive2.split(point2)[1])
            primitives.extend(self.primitives[ip2 + 1:ip1])
            primitives.append(primitive2.split(point2)[0])

        return primitives

    def extract_without_primitives(self, point1, point2, inside: bool = True):
        """
        inside: extracted contour is between the two points if True and outside
        these points if False
        """
        split_primitives = []
        primitives = self.primitives
        for point in [point1, point2]:
            dist_min = math.inf
            for primitive in primitives:
                dist = primitive.point_distance(point)
                if dist < dist_min:
                    dist_min = dist
                    prim_opt = primitive
            split_primitives.append(prim_opt)
        return self.extract_primitives(point1, split_primitives[0], point2,
                                       split_primitives[1], inside)

    def point_belongs(self, point, abs_tol=1e-7):
        '''
        find out if a point is on the wire or not. If it belongs, we return the primitive's index
        '''

        for i, primitive in enumerate(self.primitives):
            belongs = primitive.point_belongs(point, abs_tol=abs_tol)
            if belongs:
                return i
        return False

    def abscissa(self, point):
        '''
        compute the curvilinear abscisse of a point on a wire
        '''

        index = self.point_belongs(point)

        length = self.primitives[index].abscissa(point)
        for primitive in self.primitives[0:index]:
            length += primitive.length()

        return length


class Wire2D(volmdlr.core.CompositePrimitive2D, Wire):
    """
    A collection of simple primitives, following each other making a wire
    """

    def __init__(self, primitives: List[volmdlr.core.Primitive2D],
                 name: str = ''):
        volmdlr.core.CompositePrimitive2D.__init__(self, primitives, name)

    def to_3d(self, plane_origin, x, y):
        primitives3d = []
        for edge in self.primitives:
            primitives3d.append(edge.to_3d(plane_origin, x, y))

        return Wire3D(primitives3d)

    def extract(self, point1, primitive1, point2, primitive2,
                inside: bool = True):
        """
        inside: extracted contour is between the two points if True and outside these points if False
        """
        return Wire2D(
            self.extract_primitives(point1, primitive1, point2, primitive2,
                                    inside))

    def extract_with_points(self, point1: volmdlr.Point2D,
                            point2: volmdlr.Point2D, inside: bool = True):
        """
        inside: extracted contour is between the two points if True and outside these points if False
        """
        return self.extract_without_primitives(point1, point2, inside)

        # TODO: method to check if it is a wire

    def infinite_intersections(self, infinite_primitives):
        """
        returns a list  that contains:
        the intersections between a succession of infinite primitives (line,
        circle). There must be a method implemented to intersect the two
        infinite primitives.
        """
        offset_intersections = []

        for primitive_1, primitive_2 in zip(infinite_primitives,
                                            infinite_primitives[1:]):

            i = infinite_primitives.index(primitive_1)
            # k = infinite_primitives.index(primitive_2)

            primitive_name = primitive_1.__class__.__name__.lower().replace(
                '2d', '')
            intersection_method_name = '{}_intersections'.format(
                primitive_name)
            next_primitive_name = primitive_2.__class__.__name__.lower(). \
                replace('2d', '')
            next_intersection_method_name = '{}_intersections'.format(
                next_primitive_name)

            if hasattr(primitive_1, next_intersection_method_name):
                intersections = getattr(primitive_1,
                                        next_intersection_method_name)(
                    primitive_2)
                end = self.primitives[i].end

                if len(intersections) == 1:
                    offset_intersections.append(intersections[0])

                else:
                    end = self.primitives[i].end
                    if intersections[0].point_distance(end) > intersections[
                            1].point_distance(end):
                        intersections.reverse()
                    offset_intersections.append(intersections[0])

            elif hasattr(primitive_2, intersection_method_name):
                intersections = getattr(primitive_2, intersection_method_name)(
                    primitive_1)
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
                    'No intersection method between {} and {}. Define {} on {} or {} on {}'.format(
                        primitive_1.__class__.__name__,
                        primitive_2.__class__.__name__,
                        next_intersection_method_name,
                        primitive_1.__class__.__name__,
                        intersection_method_name,
                        primitive_2.__class__.__name__
                    ))

        return offset_intersections

    def offset(self, offset):
        """"
        generates an offset of a Wire2D

        """
        offset_primitives = []
        infinite_primitives = []
        offset_intersections = []
        # ax = self.plot()
        for primitive in self.primitives:
            infinite_primitive = primitive.infinite_primitive(offset)
            infinite_primitives.append(infinite_primitive)
            # infinite_primitive.plot(ax=ax, color='grey')

        offset_intersections += self.infinite_intersections(
            infinite_primitives)

        # [p.plot(ax=ax, color='r') for p in offset_intersections]

        # offset_primitives.append(
        #     self.primitives[0].border_primitive(infinite_primitives[0],
        #                                         offset_intersections[0], 0))

        # offset_primitives.append(
        #     self.primitives[-1].border_primitive(infinite_primitives[-1],
        #                                          offset_intersections[-1],
        #                                          -1))

        for j in range(len(offset_intersections) - 1):
            p1 = offset_intersections[j]
            p2 = offset_intersections[j + 1]
            cutted_primitive = infinite_primitives[
                j + 1].cut_between_two_points(p1, p2)
            offset_primitives.append(cutted_primitive)

        return Wire2D(offset_primitives)

    def plot_data(self, name: str = '', fill=None, color='black',
                  stroke_width: float = 1, opacity: float = 1):
        data = []
        for item in self.primitives:
            data.append(item.plot_data())
        return data

    def line_intersections(self, line: 'volmdlr.edges.Line2D'):
        """
        Returns a list of intersection in ther form of a tuple (point,
        primitive) of the wire primitives intersecting with the line
        """
        intersection_points = []
        for primitive in self.primitives:
            for p in primitive.line_intersections(line):
                intersection_points.append((p, primitive))
        return intersection_points

    def linesegment_intersections(self,
                                  linesegment: 'volmdlr.edges.LineSegment2D'):
        """
        Returns a list of intersection in ther form of a tuple (point,
        primitive) of the wire primitives intersecting with the line
        """
        intersection_points = []
        for primitive in self.primitives:
            for p in primitive.linesegment_intersections(linesegment):
                intersection_points.append((p, primitive))
        return intersection_points

    def is_start_end_crossings_valid(self, line, intersections, primitive):
        """
        :param line: crossing line
        :param intersections: intersections results
         for primitive line intersections
        :param primitive: intersecting primitive
        :return: None if intersection not a start or
        end point of a contours primitives, or a volmdlr.Point2D if it is.
        """
        primitive_index = self.primitives.index(primitive)
        point1, point2 = None, None
        if intersections[0] == primitive.start:
            point1 = primitive.point_at_abscissa(primitive.length() * 0.01)
            point2 = self.primitives[primitive_index - 1].point_at_abscissa(
                self.primitives[primitive_index - 1].length() * .99
            )

            # point2 = primitive.start + \
            #          self.primitives[primitive_index - 1].unit_direction_vector(0.5)
        elif intersections[0] == primitive.end and \
                primitive != self.primitives[-1]:
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
        :param intersections: intersections results
         for primitive line intersections
        :param primitive: intersecting primitive
        :return: False if intersection not a start or
        end point of a contours primitives, or True if it is.
        """
        if intersections[0] == primitive.start or intersections[0] == primitive.end:
            return True
        return False

    def line_crossings(self, line: volmdlr.edges.Line2D):
        """
        Calculates valid crossing intersections of a wire and an infinit line
        :param line: line crossing the wire
        :type line: volmdlr.edges.Line2D
        returns a list of Tuples (point, primitive)
        of the wire primitives intersecting with the line
        """
        intersection_points = []
        intersection_points_primitives = []
        for primitive in self.primitives:
            intersections = primitive.line_intersections(line)
            if intersections and intersections[0] not in intersection_points:
                if not self.is_crossing_start_end_point(intersections,
                                                        primitive):
                    intersection_points.append(intersections[0])
                    intersection_points_primitives.append(
                        (intersections[0],
                         primitive))
                elif self.is_start_end_crossings_valid(line, intersections,
                                                       primitive):
                    intersection_points.append(intersections[0])
                    intersection_points_primitives.append(
                        (intersections[0],
                         primitive))
        return intersection_points_primitives

    def wire_intersections(self, wire):
        '''
        compute intersections between two wire2d.

        Parameters
        ----------
        wire : volmdlr.wires.Wire2D

        Returns
        -------
        intersections : List[(volmdlr.Point2D, volmdlr.Primitive2D)]

        '''

        intersections = []
        for primitive in wire.primitives:
            a = self.linesegment_intersections(primitive)
            if a:
                intersections.append([a[0][0], a[0][1]])

        return intersections

    def cut_by_wire(self, wire):
        '''
        cut a wire2d (contour2d) with a wire2d
        '''
        # TODO: To be checked: it is an adapted version of 'cut_by_line'

        intersections = self.wire_intersections(wire)

        n_inter = len(intersections)
        if not intersections or n_inter < 2:
            return [self]
        if n_inter % 2 == 0:

            contours = []
            # primitives_split = [primitive.split(point) \
            #                     for point, primitive in intersections]
            x = [(ip, wire.abscissa(point))
                 for ip, (point, _) in enumerate(intersections)]
            # intersection_to_primitives_index = {
            #     i: self.primitives.index(primitive) \
            #     for i, (_, primitive) in enumerate(intersections)}
            sorted_inter_index = [x[0] for x in sorted(x, key=lambda x: x[1])]
            sorted_inter_index_dict = {i: ii for ii, i in
                                       enumerate(sorted_inter_index)}
            sorted_inter_index_dict[n_inter] = sorted_inter_index_dict[0]

            # Side 1: opposite side of begining of contour

            remaining_transitions1 = list(range(n_inter // 2))
            # enclosing_transitions = {}
            while len(remaining_transitions1) > 0:
                nb_max_enclosed_transitions = -1
                enclosed_transitions = {}
                for it in remaining_transitions1:
                    i1 = sorted_inter_index_dict[2 * it]
                    i2 = sorted_inter_index_dict[2 * it + 1]
                    net = abs(i2 - i1) - 1
                    if net > nb_max_enclosed_transitions:
                        nb_max_enclosed_transitions = net
                        best_transition = it
                        if i1 < i2:
                            enclosed_transitions[it] = [(i + 1) // 2 for i in
                                                        sorted_inter_index[
                                                        i2 - 1:i1:-2]]
                        else:
                            enclosed_transitions[it] = [(i + 1) // 2 for i in
                                                        sorted_inter_index[
                                                        i2 + 1:i1:2]]

                remaining_transitions1.remove(best_transition)
                point_start, primitive1 = intersections[2 * best_transition]
                point2, primitive2 = intersections[2 * best_transition + 1]
                primitives = self.extract_primitives(point_start, primitive1,
                                                     point2, primitive2, inside=True)
                last_point = point2
                for transition in enclosed_transitions[best_transition]:
                    point1, primitive1 = intersections[2 * transition]
                    point2, primitive2 = intersections[2 * transition + 1]
                    primitives.append(
                        volmdlr.edges.LineSegment2D(last_point, point1))
                    primitives.extend(
                        self.extract_primitives(point1, primitive1, point2,
                                                primitive2, inside=True))
                    last_point = point2
                    remaining_transitions1.remove(transition)

                primitives.append(
                    volmdlr.edges.LineSegment2D(last_point, point_start))

                # points = (volmdlr.edges.LineSegment2D(last_point, point_start)).discretise(5)
                # line_segment=volmdlr.edges.LineSegment2D(last_point, point_start)
                # for p in points[1:-1]:
                #     r = line_segment.split(p)
                #     primitives.append(r[0])
                #     line_segment = r[1]
                #     if p == points[-2]:
                #         primitives.append(r[1])

                # primitives.extend(
                #         wire.extract_with_points(last_point, point_start ,inside=True))

                contour = Contour2D(primitives)
                # contour = contour.order_contour()
                contours.append(contour)

            # Side 2: start of contour to first intersect (i=0) and  i odd to i+1 even
            intersections.append(intersections[0])

            remaining_transitions2 = list(range(n_inter // 2))
            while len(remaining_transitions2) > 0:
                nb_max_enclosed_transitions = -1
                enclosed_transitions = {}
                for it in remaining_transitions2:
                    i1 = sorted_inter_index_dict[2 * it + 1]
                    i2 = sorted_inter_index_dict[2 * it + 2]
                    net = abs(i2 - i1) - 1
                    if net > nb_max_enclosed_transitions:
                        nb_max_enclosed_transitions = net
                        best_transition = it
                        if i1 < i2:
                            enclosed_transitions[it] = [i // 2 for i in
                                                        sorted_inter_index[
                                                            i2 - 1:i1:-2]]
                        else:
                            enclosed_transitions[it] = [i // 2 for i in
                                                        sorted_inter_index[
                                                            i2 + 1:i1:2]]

                remaining_transitions2.remove(best_transition)
                point_start, primitive1 = intersections[
                    2 * best_transition + 1]
                point2, primitive2 = intersections[2 * best_transition + 2]
                primitives = self.extract_primitives(point_start, primitive1,
                                                     point2, primitive2, inside=False)
                last_point = point2
                for transition in enclosed_transitions[best_transition]:
                    point1, primitive1 = intersections[2 * transition + 1]
                    point2, primitive2 = intersections[2 * transition + 2]
                    primitives.append(
                        volmdlr.edges.LineSegment2D(last_point, point1))
                    primitives.extend(
                        self.extract_primitives(point1, primitive1, point2,
                                                primitive2, inside=False))
                    last_point = point2
                    remaining_transitions2.remove(transition)

                primitives.append(
                    volmdlr.edges.LineSegment2D(last_point, point_start))

                # points = (volmdlr.edges.LineSegment2D(last_point, point_start)).discretise(5)
                # line_segment=volmdlr.edges.LineSegment2D(last_point, point_start)
                # for p in points[1:-1]:
                #     r = line_segment.split(p)
                #     primitives.append(r[0])
                #     line_segment = r[1]
                #     if p == points[-2]:
                #         primitives.append(r[1])

                # primitives.extend(
                #         wire.extract_with_points(last_point, point_start ,inside=True))

                contour = Contour2D(primitives)
                # contour = contour.order_contour()
                contours.append(contour)

            return contours

        # ax = self.plot(equal_aspect=False)
        # # line.plot(ax=ax, color='b')
        # for point, prim in intersections:
        #     point.plot(ax=ax, color='r')
        # ax = self.plot()
        # for p in intersections:
        #     p[0].plot(ax=ax, color='r')
        # ax.set_aspect('auto')
        raise NotImplementedError(
            '{} intersections not supported yet'.format(len(intersections)))

    @classmethod
    def from_points(cls, points: List[volmdlr.Point2D]):
        '''
        define a wire based on points2d with line_segments2d
        '''

        edges = []
        for i in range(0, len(points) - 1):
            edges.append(volmdlr.edges.LineSegment2D(points[i], points[i + 1]))

        return cls(edges)


class Wire3D(volmdlr.core.CompositePrimitive3D, Wire):
    """
    A collection of simple primitives, following each other making a wire
    """

    def __init__(self, primitives: List[volmdlr.core.Primitive3D],
                 name: str = ''):
        volmdlr.core.CompositePrimitive3D.__init__(self, primitives, name)

    def extract(self, point1, primitive1, point2, primitive2):
        return Wire3D(self.extract_primitives(self, point1, primitive1, point2,
                                              primitive2))

    def extract_with_points(self, point1: volmdlr.Point3D,
                            point2: volmdlr.Point3D, inside):
        return self.extract_without_primitives(point1, point2, inside)

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

    def extrusion(self, extrusion_vector):
        faces = []
        for primitive in self.primitives:
            faces.extend(primitive.extrusion(extrusion_vector))
        return faces

    def to_bspline(self, discretization_parameter, degree):
        '''
        convert a wire3d to a bspline curve3d
        '''

        discretized_points = self.discretization_points(discretization_parameter)
        bspline_curve = volmdlr.edges.BSplineCurve3D.from_points_interpolation(discretized_points, degree)

        return bspline_curve

    # def copy(self, deep=True, memo=None):
    #     primitives_copy = []
    #     for primitive in self.primitives:
    #         primitives_copy.append(primitive.copy())
    #     return Wire3D(primitives_copy)


# TODO: define an edge as an opened polygon and allow to compute area from this reference

class Contour(Wire):

    # def __init__(self):
    #     Wire.__init__(self)

    def extract_primitives(self, point1, primitive1, point2, primitive2, inside: bool = True):
        """
        inside: extracted contour is between the two points if True and outside these points if False
        """
        primitives = []
        ip1 = self.primitive_to_index(primitive1)
        ip2 = self.primitive_to_index(primitive2)
        if inside:
            if ip1 < ip2:
                pass
            elif ip1 == ip2:  # primitive1 == primitive2
                if point1.point_distance(
                            primitive1.start) < point2.point_distance(
                        primitive1.start):
                    pass
                else:
                    primitive1, primitive2 = primitive2, primitive1
                    point1, point2 = point2, point1

            else:
                primitive1, primitive2 = primitive2, primitive1
                point1, point2 = point2, point1
        else:
            if ip1 > ip2:
                pass
            elif ip1 == ip2:  # primitive1 == primitive2
                if point1.point_distance(
                            primitive1.start) > point2.point_distance(
                        primitive1.start):
                    pass
                else:
                    primitive1, primitive2 = primitive2, primitive1
                    point1, point2 = point2, point1
            else:
                primitive1, primitive2 = primitive2, primitive1
                point1, point2 = point2, point1

        ip1 = self.primitive_to_index(primitive1)
        ip2 = self.primitive_to_index(primitive2)

        if ip1 < ip2:
            if primitive1.start == point1:
                primitives.append(primitive1)
            elif primitive1.end == point1:
                pass
            else:
                primitives.append(primitive1.split(point1)[1])
            primitives.extend(self.primitives[ip1 + 1:ip2])
            if primitive2.start == point2:
                pass
            elif primitive2.end == point2:
                primitives.append(primitive2)
            else:
                primitives.append(primitive2.split(point2)[0])
        elif ip1 > ip2 or (ip1 == ip2 and point1.point_distance(
                primitive1.start) > point2.point_distance(primitive1.start)):
            if primitive1.start == point1:
                primitives.append(primitive1)
            elif primitive1.end == point1:
                pass
            else:
                primitives.append(primitive1.split(point1)[1])
            # primitives.append(primitive1.split(point1)[1])
            primitives.extend(self.primitives[ip1 + 1:])
            primitives.extend(self.primitives[:ip2])
            if primitive2.start == point2:
                pass
            elif primitive2.end == point2:
                primitives.append(primitive2)
            else:
                primitives.append(primitive2.split(point2)[0])
        elif (ip1 == ip2 and point1.point_distance(
                primitive1.start) < point2.point_distance(primitive1.start)):
            if primitive1.start != point1:
                primitive = primitive1.split(point1)[1]
                primitive = primitive.split(point2)[0]
                primitives.append(primitive)
            elif primitive1.end != point2:
                primitive = primitive1.split(point2)[0]
                primitives.append(primitive)
            else:
                primitives.append(primitive2)

        return primitives

    def is_ordered(self):
        for prim1, prim2 in zip(
                self.primitives, self.primitives[1:] + [self.primitives[0]]):
            if prim1.end != prim2.start:
                return False
        return True

    def ordering_contour(self, tol=1e-6):
        """
        returns the points of the contour ordered
        """
        list_point_pairs = [(prim.start, prim.end) for prim in self.primitives]
        length_list_points = len(list_point_pairs)
        points = [list_point_pairs[0]]
        list_point_pairs.remove(
            (list_point_pairs[0][0], list_point_pairs[0][1]))
        finished = False
        counter = 0
        counter1 = 0
        while not finished:
            for p1, p2 in list_point_pairs:
                if p1.point_distance(p2) < tol:
                    list_point_pairs.remove((p1, p2))
                elif p1.point_distance(points[-1][-1]) < tol:
                    points.append((p1, p2))
                    list_point_pairs.remove((p1, p2))
                elif p2.point_distance(points[-1][-1]) < tol:
                    points.append((p2, p1))
                    list_point_pairs.remove((p1, p2))
                elif p1.point_distance(points[0][0]) < tol:
                    points = [(p2, p1)] + points
                    list_point_pairs.remove((p1, p2))
                elif p2.point_distance(points[0][0]) < tol:
                    points = [(p1, p2)] + points
                    list_point_pairs.remove((p1, p2))
            if len(list_point_pairs) == 0:
                finished = True
            counter1 += 1
            if counter1 >= 100 * length_list_points:
                raise NotImplementedError
            if len(list_point_pairs) == 1:
                counter += 1
                if counter > 3:
                    warnings.warn('There may exist a problem with this'
                                  ' contour, it seems it cannot be reordered.'
                                  ' Please, verify its points')
                    raise NotImplementedError

        return points

    @classmethod
    def contours_from_edges(cls, edges, tol=5e-5):
        list_contours = []
        finished = False
        contour = []
        while not finished:
            len1 = len(edges)
            for line in edges:
                points = [p for prim in contour for p in prim]
                if not contour:
                    contour.append(line)
                    edges.remove(line)
                    break
                elif line.start in points or line.end in points:
                    contour.append(line)
                    edges.remove(line)
                    break
                else:
                    for point in points:

                        if point.is_close(line.start, tol=tol):
                            line.start = point
                            contour.append(line)
                            edges.remove(line)
                            break
                        elif point.is_close(line.end, tol=tol):
                            line.end = point
                            contour.append(line)
                            edges.remove(line)
                            break

            if len(edges) != 0 and len(edges) == len1 and len(contour) != 0:
                contour_n = cls(contour[:])
                contour_n.order_contour()
                list_contours.append(contour_n)
                contour = []
            elif len(edges) == 0 and len(contour) != 0:
                contour_n = cls(contour[:])
                contour_n.order_contour()
                list_contours.append(contour_n)
                finished = True
        return list_contours

    def discretized_primitives(self, n: float):
        """
        discretize each contour's primitive and return a list of discretized primitives
        """

        edges = []
        for primitive in self.primitives:
            points = primitive.discretise(n)
            for p1, p2 in zip(points[:-1], points[1:]):
                edges.append(volmdlr.edges.LineSegment2D(p1, p2))

        return edges

    def shares_primitives(self, contour):
        """checks if two contour share primitives"""
        for prim1 in self.primitives:
            if contour.primitive_over_contour(prim1):
                return True
        return False

    def is_sharing_primitives_with(self, contour, all_points=False):
        '''
        check if two contour are sharing primitives
        "all_points" is by default False. Turn it True if you need to get points and edges used to find out shared primitives
        '''

        list_p = []
        edges1 = set()
        # edges2 = set()

        for edge_1, edge_2 in itertools.product(self.primitives,
                                                contour.primitives):
            edges = [edge_1, edge_2, edge_1]
            for edge1, edge2 in zip(edges, edges[1:]):
                for point in [edge2.start, edge2.end]:
                    if edge1.point_belongs(point, 1e-6):
                        # list_p.append(point)
                        # instead of point not in list_p (due to errors)
                        if list_p == []:
                            list_p.append(point)
                        if list_p != [] and point.point_distance(point.nearest_point(list_p)) > 1e-4:
                            list_p.append(point)

                        # edges2.add(edge2)
                        try:
                            self.primitive_to_index(edge1)
                            edges1.add(edge1)
                        except KeyError:
                            edges1.add(edge2)

                    if len(list_p) == 2 and all_points is False:
                        return True
        if len(list_p) < 2:
            return False
        if len(list_p) >= 2 and all_points is True:
            return (edges1, list_p)
        else:
            return False

    def shared_primitives_extremities(self, contour):
        '''
        extract shared primitives extremities between two adjacent contours
        '''

        results = self.is_sharing_primitives_with(contour, all_points=True)
        if results is False:
            raise ValueError(
                    'The contours are not adjacent. They dont share primitives')
        elif len(results[1]) == 2:
            return results[1]
        else:
            edges1, list_p = results
            if isinstance(self, volmdlr.wires.Contour2D):
                contours1 = volmdlr.wires.Contour2D.contours_from_edges(edges1)
                # contours2 = volmdlr.wires.Contour2D.contours_from_edges(edges2)
            elif isinstance(self, volmdlr.wires.Contour3D):
                contours1 = volmdlr.wires.Contour3D.contours_from_edges(edges1)
                # contours2 = volmdlr.wires.Contour3D.contours_from_edges(edges2)
            points = []

            for contour_i in contours1:
                points.extend(contour_i.extremities_points(list_p))

            return points

    def shared_primitives_with(self, contour):
        '''
        extract shared primitives between two adjacent contours
        '''

        shared_primitives_1 = []
        shared_primitives_2 = []

        points = self.shared_primitives_extremities(contour)
        for i in range(0, len(points), 2):
            point1, point2 = points[i], points[i + 1]

            shared_primitives_prim = self.extract_without_primitives(point1,
                                                                     point2,
                                                                     False)
            if contour.point_over_contour(shared_primitives_prim[0].middle_point(), 1e-4) is False:
                shared_primitives_1.extend(self.extract_without_primitives(point1,
                                                                           point2,
                                                                           True))
            else:
                shared_primitives_1.extend(shared_primitives_prim)

            shared_primitives_prim = contour.extract_without_primitives(point1,
                                                                        point2,
                                                                        False)
            if self.point_over_contour(shared_primitives_prim[0].middle_point(), 1e-4) is False:
                shared_primitives_2.extend(contour.extract_without_primitives(point1,
                                                                              point2,
                                                                              True))
            else:
                shared_primitives_2.extend(shared_primitives_prim)

        return [shared_primitives_1, shared_primitives_2]

    def merge_primitives_with(self, contour):
        '''
        extract not shared primitives between two adjacent contours, to be merged
        '''

        points = self.shared_primitives_extremities(contour)
        merge_primitives = []

        for i in range(1, len(points) + 1, 2):
            if i == (len(points) - 1):
                point1, point2 = points[i], points[0]
            else:
                point1, point2 = points[i], points[i + 1]

            merge_primitives_prim = self.extract_without_primitives(point1,
                                                                    point2,
                                                                    False)
            if contour.point_over_contour(merge_primitives_prim[0].middle_point(), 1e-4) is True:
                merge_primitives_prim = self.extract_without_primitives(point1,
                                                                        point2,
                                                                        True)
                merge_primitives.extend(merge_primitives_prim)
            else:
                merge_primitives.extend(merge_primitives_prim)

            merge_primitives_prim = contour.extract_without_primitives(point1,
                                                                       point2,
                                                                       False)
            if self.point_over_contour(merge_primitives_prim[0].middle_point(), 1e-4) is True:
                merge_primitives_prim = contour.extract_without_primitives(point1,
                                                                           point2,
                                                                           True)
                merge_primitives.extend(merge_primitives_prim)
            else:
                merge_primitives.extend(merge_primitives_prim)

        return merge_primitives

    def edges_order_with_adjacent_contour(self, contour):
        """
        check if the shared edges between two adjacent contours are traversed with two different directions along each contour
        """

        contour1 = self
        contour2 = contour

        # shared_tuple = contour1.shared_edges_between2contours(contour2)
        shared_tuple = contour1.shared_primitives_with(contour2)
        # [shared_primitives_1, shared_primitives_2] = contour1.shared_primitives_with(contour2)

        p1_start = contour1.primitives[shared_tuple[0][0]].start
        p2_start = contour2.primitives[shared_tuple[0][1]].start
        p2_end = contour2.primitives[shared_tuple[0][1]].end

        if (p1_start.point_distance(p2_start)) < (p1_start.point_distance(p2_end)):
            return False
        else:
            return True

    def extremities_points(self, list_p):
        '''
        return extremitises points of a list of points on a contour
        '''
        points = []
        primitives = self.primitives
        for i in range(0, len(primitives)):
            pts = []
            for p in list_p:  # due to errors
                if primitives[i].point_belongs(p):
                    pts.append(p)
            if len(pts) == 1:
                points.append(pts[0])
                break
            elif len(pts) > 1:
                points.append(primitives[i].start.nearest_point(pts))
                break

        for i in range(len(primitives) - 1, -1, -1):
            pts = []
            for p in list_p:  # due to errors
                if primitives[i].point_belongs(p):
                    pts.append(p)
            if len(pts) == 1:
                if pts[0] not in points:
                    points.append(pts[0])
                    break
            elif len(pts) > 1:
                point = primitives[i].end.nearest_point(pts)
                if point not in points:
                    points.append(point)
                    break
        return points


class Contour2D(Contour, Wire2D):
    """
    A collection of 2D primitives forming a closed wire2D
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

    def __init__(self, primitives: List[volmdlr.core.Primitive2D],
                 name: str = ''):
        Wire2D.__init__(self, primitives, name)
        self._utd_edge_polygon = False

    # def __eq__(self, other_):
    #     if other_.__class__.__name__ != self.__class__.__name__:
    #         return False
    #     if len(self.primitives) != len(other_.primitives):
    #         return False
    #     equal = True
    #     for prim1, prim2 in zip(self.primitives, other_.primitives):
    #         equal = (equal and prim1 == prim2)
    #     return equal

    @property
    def edge_polygon(self):
        if not self._utd_edge_polygon:
            self._edge_polygon = self._get_edge_polygon()
            self._utd_edge_polygon = True
        return self._edge_polygon

    def _get_edge_polygon(self):
        points = []
        for edge in self.primitives:

            if points:
                if edge.start != points[-1]:
                    points.append(edge.start)
            else:
                points.append(edge.start)
        return ClosedPolygon2D(points)

    def to_3d(self, plane_origin, x, y):
        p3d = []
        for edge in self.primitives:
            p3d.append(edge.to_3d(plane_origin, x, y))

        return Contour3D(p3d)

    def point_belongs(self, point):
        if self.edge_polygon.point_belongs(point):
            return True
        # TODO: This is incomplete!!!
        return False

    def point_over_contour(self, point, abs_tol=1e-6):
        belongs = False
        for primitive in self.primitives:
            if primitive.point_belongs(point, abs_tol):
                belongs = True
        return belongs

    def primitive_over_contour(self, primitive, tol: float = 1e-6):
        for prim in self.primitives:
            if not hasattr(prim, 'unit_direction_vector') and \
                    hasattr(prim, 'tangent'):
                vector1 = prim.tangent(0.5)
            else:
                vector1 = prim.unit_direction_vector(0.5)

            if not hasattr(primitive, 'unit_direction_vector') and \
                    hasattr(primitive, 'tangent'):
                vector2 = primitive.tangent(0.5)
            else:
                vector2 = primitive.unit_direction_vector(0.5)
            if vector1.is_colinear_to(vector2):
                mid_point = primitive.middle_point()
                if self.point_over_contour(mid_point, tol):
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
        points = self.edge_polygon.points[:]
        for primitive in self.primitives:
            if hasattr(primitive, 'polygon_points'):
                points.extend(primitive.polygon_points())
        xmin = min([p[0] for p in points])
        xmax = max([p[0] for p in points])
        ymin = min([p[1] for p in points])
        ymax = max([p[1] for p in points])
        return (volmdlr.Point2D(xmin, ymin), volmdlr.Point2D(xmax, ymax))

    def area(self):
        area = self.edge_polygon.area()
        if self.edge_polygon.is_trigo():
            trigo = 1
        else:
            trigo = -1
        for edge in self.primitives:
            area += trigo * edge.straight_line_area()

        return area

    def center_of_mass(self):
        center = self.edge_polygon.area() * self.edge_polygon.center_of_mass()
        # ax = self.plot()
        # self.edge_polygon.center_of_mass().plot(ax=ax, color='b')
        if self.edge_polygon.is_trigo():
            trigo = 1
        else:
            trigo = -1
        for edge in self.primitives:
            # edge.straight_line_center_of_mass().plot(ax=ax, color='g')
            center += trigo * edge.straight_line_area() \
                      * edge.straight_line_center_of_mass()

        return center / self.area()

    def second_moment_area(self, point):

        Ix, Iy, Ixy = self.edge_polygon.second_moment_area(point)
        for edge in self.primitives:
            Ix_e, Iy_e, Ixy_e = edge.straight_line_second_moment_area(point)
            if self.edge_polygon.is_trigo():
                Ix += Ix_e
                Iy += Iy_e
                Ixy += Ixy_e
            else:
                Ix -= Ix_e
                Iy -= Iy_e
                Ixy -= Ixy_e

        return Ix, Iy, Ixy

    def plot_data(self, edge_style: plot_data.EdgeStyle = None,
                  surface_style: plot_data.SurfaceStyle = None):
        plot_data_primitives = [item.plot_data() for item in self.primitives]
        return plot_data.Contour2D(plot_data_primitives=plot_data_primitives,
                                   edge_style=edge_style,
                                   surface_style=surface_style,
                                   name=self.name)

    def is_inside(self, contour2):
        """
        verifies if a contour is inside another contour perimiter,
        including the edges
        return True or False
        """
        points_contour2 = []
        for prim in contour2.primitives:
            if prim.start not in points_contour2:
                points_contour2.append(prim.start)
            if prim.end not in points_contour2:
                points_contour2.append(prim.end)
            points_contour2.append(prim.middle_point())
        for point in points_contour2:
            if not self.point_belongs(point) and not self.point_over_contour(point, abs_tol=1e-7):
                return False
        return True

    def bounding_rectangle(self):
        xmin, xmax, ymin, ymax = self.primitives[0].bounding_rectangle()
        for edge in self.primitives[1:]:
            xmin_edge, xmax_edge, ymin_edge, ymax_edge = \
                edge.bounding_rectangle()
            xmin = min(xmin, xmin_edge)
            xmax = max(xmax, xmax_edge)
            ymin = min(ymin, ymin_edge)
            ymax = max(ymax, ymax_edge)
        return xmin, xmax, ymin, ymax

    def inverted_primitives(self):
        new_primitives = []
        for prim in self.primitives[::-1]:
            new_primitives.append(
                volmdlr.edges.LineSegment2D(prim.end, prim.start))
        return new_primitives

    def invert(self):
        return Contour2D(self.inverted_primitives())

    def invert_inplace(self):
        self.primitives = self.inverted_primitives()

    def random_point_inside(self):
        xmin, xmax, ymin, ymax = self.bounding_rectangle()
        for i in range(2000):
            p = volmdlr.Point2D.random(xmin, xmax, ymin, ymax)
            if self.point_belongs(p):
                return p

    def order_contour(self):
        if self.is_ordered() or len(self.primitives) < 2:
            return self

        initial_points = []
        for p in self.primitives:
            initial_points.append((p.start, p.end))

        new_primitives = []
        points = self.ordering_contour()
        for p1, p2 in points:
            try:
                index = initial_points.index((p1, p2))
            except ValueError:
                index = initial_points.index((p2, p1))

            if isinstance(self.primitives[index], volmdlr.edges.LineSegment2D):
                new_primitives.append(volmdlr.edges.LineSegment2D(p1, p2))
            elif isinstance(self.primitives[index], volmdlr.edges.Arc2D):
                new_primitives.append(volmdlr.edges.Arc2D(p1, self.primitives[index].interior, p2))
            elif isinstance(self.primitives[index], volmdlr.edges.BSplineCurve2D):
                if (p1, p2) == (self.primitives[index].start, self.primitives[index].end):
                    new_primitives.append(self.primitives[index])
                else:
                    new_primitives.append(self.primitives[index].reverse())

        self.primitives = new_primitives

        return self

    @classmethod
    def extract_contours(cls, contour, point1: volmdlr.Point3D,
                         point2: volmdlr.Point3D, inside=False):

        new_primitives = contour.extract_with_points(point1, point2, inside)
        contours = [cls(new_primitives)]
        return contours

    def cut_by_linesegments(self, lines: List[volmdlr.edges.LineSegment2D]):
        for c in lines:
            if not isinstance(c, volmdlr.edges.LineSegment2D):
                raise KeyError(
                    'contour must be a list of LineSegment2D object')

        cut_lines = []
        for p in lines:
            cut_lines.append(p.to_line())

        contour_to_cut = [self]
        for l in cut_lines:
            new_contour_to_cut = []
            for c in contour_to_cut:
                cs = c.cut_by_line(l)
                new_contour_to_cut.extend(cs)
            contour_to_cut.extend(new_contour_to_cut)

        p1 = Contour2D(lines).center_of_mass()
        dist_min = math.inf
        for c in contour_to_cut:
            if c.area() > 1e-10:
                p0 = c.center_of_mass()
                if p0.point_distance(p1) < dist_min:
                    c_opti = c
                    dist_min = p0.point_distance(p1)
        return c_opti

    def cut_by_line(self,
                    line: volmdlr.edges.Line2D) -> List['Contour2D']:
        """
        :param line: The line used to cut the contour
        :return: A list of resulting contours
        """
        intersections = self.line_crossings(line)
        if not intersections or len(intersections) < 2:
            return [self]
        if len(intersections) % 2 != 0:
            raise NotImplementedError(
                '{} intersections not supported yet'.format(
                    len(intersections)))

        points_intersections = [point for point, prim in intersections]
        sorted_points = line.sort_points_along_line(points_intersections)
        list_contours = []
        contour_to_cut = self
        cutting_points_counter = 0
        while cutting_points_counter != len(sorted_points):

            point1 = sorted_points[cutting_points_counter]
            point2 = sorted_points[cutting_points_counter + 1]
            closing_line = volmdlr.edges.LineSegment2D(point1, point2)
            closing_contour = Contour2D([closing_line])
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

    def get_pattern(self):
        """ A pattern is portion of the contour from which the contour can be
        reconstructed by rotations of this portion"""
        xmin, xmax, ymin, ymax = self.bounding_rectangle()

        # ax=plt.subplot()
        # line = Line2D(Point2D([xi, 0]),Point2D([xi,1]))
        line = volmdlr.edges.Line2D(volmdlr.Point2D([0, -0.17]),
                                    volmdlr.Point2D([0, 0.17]))
        line_2 = line.Rotation(self.center_of_mass(), 0.26)
        line_3 = line.Rotation(self.center_of_mass(), -0.26)

        intersections = []

        intersections += self.line_intersections(line_2)
        intersections += self.line_intersections(line_3)
        if isinstance(intersections[0][0], volmdlr.Point2D) and \
                isinstance(intersections[1][0], volmdlr.Point2D):
            ip1, ip2 = sorted([self.primitives.index(intersections[0][1]),
                               self.primitives.index(intersections[1][1])])

            ip3, ip4 = sorted([self.primitives.index(intersections[2][1]),
                               self.primitives.index(intersections[3][1])])

            sp11, sp12 = intersections[1][1].split(intersections[1][0])
            sp22, sp21 = intersections[2][1].split(intersections[2][0])

            primitives = []

            a = volmdlr.edges.Arc2D(sp12.end, sp12.interior, sp12.start)
            primitives.append(a)
            primitives.extend(self.primitives[:ip3])
            primitives.append(sp22)
            l = volmdlr.edges.LineSegment2D(sp22.start, sp12.end)
            interior = l.point_at_abscissa(l.Length() / 2)
            primitives.append(
                volmdlr.edges.Arc2D(sp22.start, interior, sp12.end))

        return Contour2D(primitives)

    def contour_from_pattern(self):
        # TODO: Delete this function?
        pattern = self.get_pattern()
        pattern_rotations = []
        # pattern_rotations.append(self)
        for k in range(1, 13):
            new_pattern = pattern.rotation(self.center_of_mass(),
                                           k * math.pi / 6)
            pattern_rotations.append(new_pattern)

        return pattern_rotations

    def simple_triangulation(self):
        lpp = len(self.polygon.points)
        if lpp == 3:
            return self.polygon.points, [(0, 1, 2)]
        elif lpp == 4:
            return self.polygon.points, [(0, 1, 2), (0, 2, 3)]

        # Use delaunay triangulation
        tri = Delaunay([p.vector for p in self.polygon.points])
        # indices = tri.simplices
        return self.polygon.points, tri.simplices

    def split_regularly(self, n):
        """
        Split in n slices
        """
        xmin, xmax, ymin, ymax = self.bounding_rectangle()
        cutted_contours = []
        iteration_contours = [self]
        for i in range(n - 1):
            xi = xmin + (i + 1) * (xmax - xmin) / n
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

    def to_polygon(self, angle_resolution):

        polygon_points = []
        # print([(line.start, line.end) for line in self.primitives])

        for primitive in self.primitives:
            polygon_points.extend(primitive.polygon_points()[:-1])
        #     print('1: ', primitive.polygon_points())
        #     print('2 :', primitive.polygon_points()[:-1])
        # print(polygon_points)
        return ClosedPolygon2D(polygon_points)

    def grid_triangulation(self, x_density: float = None,
                           y_density: float = None,
                           min_points_x: int = 20,
                           min_points_y: int = 20,
                           number_points_x: int = None,
                           number_points_y: int = None):
        """
        Use a n by m grid to triangulize the contour
        """
        bounding_rectangle = self.bounding_rectangle()
        # xmin, xmax, ymin, ymax = self.bounding_rectangle()
        dx = bounding_rectangle[1] - bounding_rectangle[0]  # xmax - xmin
        dy = bounding_rectangle[3] - bounding_rectangle[2]  # ymax - ymin
        if number_points_x is None:
            n = max(math.ceil(x_density * dx), min_points_x)
        else:
            n = number_points_x
        if number_points_y is None:
            m = max(math.ceil(y_density * dy), min_points_y)
        else:
            m = number_points_y

        x = [bounding_rectangle[0] + i * dx / n for i in range(n + 1)]
        y = [bounding_rectangle[2] + i * dy / m for i in range(m + 1)]

        point_index = {}
        ip = 0
        points = []
        triangles = []
        for xi in x:
            for yi in y:
                p = volmdlr.Point2D(xi, yi)
                if self.point_belongs(p):
                    point_index[p] = ip
                    points.append(p)
                    ip += 1

        for i in range(n):
            for j in range(m):
                p1 = volmdlr.Point2D(x[i], y[j])
                p2 = volmdlr.Point2D(x[i + 1], y[j])
                p3 = volmdlr.Point2D(x[i + 1], y[j + 1])
                p4 = volmdlr.Point2D(x[i], y[j + 1])
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

        return vmd.DisplayMesh2D(points, triangles)

    # def extract_contours(self, point1: volmdlr.Point2D, point2: volmdlr.Point2D):
    #     split_primitives  = []
    #     # primitives = [p for p in contour.primitives]
    #     primitives = self.primitives
    #     for point in [point1, point2]:
    #         dist_min = math.inf
    #         for primitive in primitives:
    #             # print(point)
    #             dist = primitive.point_distance(point)
    #             if dist < dist_min:
    #                 dist_min = dist
    #                 prim_opt = primitive
    #         split_primitives.append(prim_opt)
    #     print(len(split_primitives))
    #     return self.extract_primitives(point1, split_primitives[0], point2, split_primitives[1])

    def contour_intersections(self, contour2d):
        intersecting_points = []
        for primitive1 in self.primitives:
            for primitive2 in contour2d.primitives:
                line_intersection = primitive1.linesegment_intersections(
                    primitive2)
                if line_intersection:
                    if line_intersection[0] not in intersecting_points:
                        intersecting_points.extend(line_intersection)
                else:
                    point1, point2 = contour2d.primitives[0].start, \
                                     contour2d.primitives[-1].end
                    if primitive1.point_belongs(point1):
                        intersecting_points.append(point1)
                    if primitive1.point_belongs(point2):
                        intersecting_points.append(point2)
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
        if extracted_outerpoints_contour1.primitives[0].start == \
                closing_contour.primitives[0].start:
            cutting_contour_new = closing_contour.invert()
            primitives1 = cutting_contour_new.primitives + \
                extracted_outerpoints_contour1.primitives
        elif extracted_outerpoints_contour1.primitives[0].start == \
                closing_contour.primitives[-1].end:
            primitives1 = closing_contour.primitives + \
                          extracted_outerpoints_contour1.primitives

        if extracted_innerpoints_contour1.primitives[0].start == \
                closing_contour.primitives[0].start:
            cutting_contour_new = \
                closing_contour.invert()
            primitives2 = cutting_contour_new.primitives + \
                extracted_innerpoints_contour1.primitives
        elif extracted_innerpoints_contour1.primitives[
                0].start == closing_contour.primitives[-1].end:
            primitives2 = closing_contour.primitives + \
                          extracted_innerpoints_contour1.primitives
        contour1 = volmdlr.wires.Contour2D(primitives1)
        contour1.order_contour()
        contour2 = volmdlr.wires.Contour2D(primitives2)
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
            cutting_contour = contours[0]
            for base_contour in new_base_contours:
                cutting_points = []
                point1, point2 = [cutting_contour.primitives[0].start,
                                  cutting_contour.primitives[-1].end]
                if base_contour.point_over_contour(
                        point1) and base_contour.point_over_contour(point2):
                    cutting_points = [point1, point2]
                if cutting_points:
                    contour1, contour2 = base_contour.get_divided_contours(
                        cutting_points[0], cutting_points[1], cutting_contour,
                        inside)

                    new_base_contours.remove(base_contour)
                    for cntr in [contour1, contour2]:
                        valid_contour = True
                        for cut_contour in list_cutting_contours:
                            point_at_abs = cut_contour.point_at_abscissa(
                                cut_contour.length() / 2)
                            if cntr.point_belongs(
                                    point_at_abs) and (not cntr.point_over_contour(
                                    point_at_abs) or True not in [
                                    cntr.primitive_over_contour(
                                        prim) for prim in cut_contour.primitives]):
                                valid_contour = False

                        if valid_contour and cntr.area() != 0.0:
                            list_valid_contours.append(cntr)
                        else:

                            new_base_contours.append(cntr)
                    contours.remove(cutting_contour)
                    break
            if len(contours) == 0:
                finished = True
            counter += 1
            if counter >= 100 * len(list_contour):
                warnings.warn('There probably exists an open contour (two wires that could not be connected)')
                finished = True

        return list_valid_contours

    def discretized_contour(self, n: float):
        """
        discretize each contour's primitive and return a new contour with teses discretized primitives
        """

        contour = volmdlr.wires.Contour2D((self.discretized_primitives(n)))

        return contour.order_contour()

    @classmethod
    def from_bounding_rectangle(cls, xmin, xmax, ymin, ymax):
        """
        create a contour2d with bounding_box parameters, using linesegments2d
        """

        edge0 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(xmin, ymin), volmdlr.Point2D(xmax, ymin))
        edge1 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(xmax, ymin), volmdlr.Point2D(xmax, ymax))
        edge2 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(xmax, ymax), volmdlr.Point2D(xmin, ymax))
        edge3 = volmdlr.edges.LineSegment2D(volmdlr.Point2D(xmin, ymax), volmdlr.Point2D(xmin, ymin))

        edges = [edge0, edge1, edge2, edge3]

        return volmdlr.wires.Contour2D(edges)

    @classmethod
    def from_points(cls, points: List[volmdlr.Point2D]):
        '''
        create a contour2d from points with line_segments2D
        '''

        if len(points) < 3:
            raise ValueError('contour is defined at least with three points')
        else:
            edges = []
            for i in range(0, len(points) - 1):
                edges.append(volmdlr.edges.LineSegment2D(points[i], points[i + 1]))

            edges.append(volmdlr.edges.LineSegment2D(points[-1], points[0]))

            contour = cls(edges)

            return contour

    def cut_by_bspline_curve(self, bspline_curve2d: volmdlr.edges.BSplineCurve2D):
        '''
        cut a contou2d with bspline_curve2d to define two different contours
        '''
        # TODO: BsplineCurve is descretized and defined with a wire. To be improved!

        contours = self.cut_by_wire(bspline_curve2d.to_wire(20))

        return contours

    def clean_primitives(self):
        '''
        delete primitives with start=end, and return a new contour
        '''

        new_primitives = []
        for p in self.primitives:
            if p.start != p.end:
                new_primitives.append(p)

        return Contour2D(new_primitives)

    def merge_with(self, contour2d):
        '''
        merge two adjacent contours, sharing primitives, and returns one outer contour and inner contours (if there are any)
        '''

        if self.is_inside(contour2d):
            return [self]
        elif contour2d.is_inside(self):
            return [contour2d]
        merged_primitives = self.merge_primitives_with(contour2d)
        contours = volmdlr.wires.Contour2D.contours_from_edges(merged_primitives)
        contours = sorted(contours, key=lambda contour: contour.area(), reverse=True)

        return contours


class ClosedPolygon:

    def length(self):
        list_ = []
        for k in range(len(self.line_segments)):
            list_.append(self.line_segments[k].length())
        return sum(list_)

    def min_length(self):
        list_ = []
        for k in range(len(self.line_segments)):
            list_.append(self.line_segments[k].length())
        return min(list_)

    def max_length(self):
        list_ = []
        for k in range(len(self.line_segments)):
            list_.append(self.line_segments[k].length())
        return max(list_)

    def edge_statistics(self):
        distances = []
        for i, point in enumerate(self.points):
            if i != 0:
                distances.append(point.point_distance(self.points[i - 1]))
        mean_distance = mean(distances)
        std = npy.std(distances)
        return mean_distance, std

    def simplify_polygon(self, min_distance: float = 0.01,
                         max_distance: float = 0.05, angle: float = 20):
        points = [self.points[0]]
        previous_point = None
        for i, point in enumerate(self.points[1:]):
            distance = point.point_distance(points[-1])
            if distance > min_distance:
                if distance > max_distance:
                    number_segmnts = round(distance / max_distance) + 2
                    for n in range(number_segmnts):
                        new_point = points[-1] + (point - points[-1]) * (
                                n + 1) / number_segmnts
                        distance1 = new_point.point_distance(points[-1])
                        if distance1 > max_distance:
                            points.append(new_point)
                else:

                    if point not in points:
                        points.append(point)
            elif len(points) > 1:
                vector1 = points[-1] - points[-2]
                vector2 = point - points[-2]
                cos = vector1.dot(vector2) / (vector1.norm() * vector2.norm())
                cos = math.degrees(math.acos(round(cos, 6)))
                if abs(cos) > angle:
                    if previous_point not in points:
                        points.append(previous_point)
                    if point not in points:
                        points.append(point)
            if len(points) > 2:
                distance2 = points[-3].point_distance(points[-2])
                vector1 = points[-2] - points[-3]
                vector2 = points[-1] - points[-3]
                cos = vector1.dot(vector2) / (vector1.norm() * vector2.norm())
                cos = math.degrees(math.acos(round(cos, 6)))
                if distance2 < min_distance and cos < angle:
                    points = points[:-2] + [points[-1]]
            previous_point = point
        distance = points[0].point_distance(points[-1])
        if distance < min_distance:
            points.remove(points[-1])

        if volmdlr.wires.ClosedPolygon2D(points).area() == 0.0:
            return self

        return self.__class__(points)


class ClosedPolygon2D(Contour2D, ClosedPolygon):
    _non_serializable_attributes = ['line_segments', 'primitives',
                                    'basis_primitives']

    def __init__(self, points: List[volmdlr.Point2D], name: str = ''):
        self.points = points
        self.line_segments = self._line_segments()

        Contour2D.__init__(self, self.line_segments, name)

    def copy(self, *args, **kwargs):
        points = [p.copy() for p in self.points]
        return ClosedPolygon2D(points, self.name)

    def __hash__(self):
        return sum([hash(p) for p in self.points])

    def __eq__(self, other_):
        if not isinstance(other_, self.__class__):
            return False
        equal = True
        for point, other_point in zip(self.points, other_.points):
            equal = (equal and point == other_point)
        return equal

    def area(self):
        # TODO: perf: cache number of points
        if len(self.points) < 3:
            return 0.

        x = [point.x for point in self.points]
        y = [point.y for point in self.points]

        x1 = [x[-1]] + x[0:-1]
        y1 = [y[-1]] + y[0:-1]
        return 0.5 * abs(sum([i * j for i, j in zip(x, y1)])
                         - sum([i * j for i, j in zip(y, x1)]))
        # return 0.5 * npy.abs(
        #     npy.dot(x, npy.roll(y, 1)) - npy.dot(y, npy.roll(x, 1)))

    def center_of_mass(self):
        lp = len(self.points)
        if lp == 0:
            return volmdlr.O2D
        elif lp == 1:
            return self.points[0]
        elif lp == 2:
            return 0.5 * (self.points[0] + self.points[1])

        x = [point.x for point in self.points]
        y = [point.y for point in self.points]

        xi_xi1 = x + npy.roll(x, -1)
        yi_yi1 = y + npy.roll(y, -1)
        xi_yi1 = npy.multiply(x, npy.roll(y, -1))
        xi1_yi = npy.multiply(npy.roll(x, -1), y)

        a = 0.5 * npy.sum(xi_yi1 - xi1_yi)  # signed area!
        # print('a :', a)
        #        a=self.area()
        if not math.isclose(a, 0, abs_tol=1e-08):
            cx = npy.sum(npy.multiply(xi_xi1, (xi_yi1 - xi1_yi))) / 6. / a
            cy = npy.sum(npy.multiply(yi_yi1, (xi_yi1 - xi1_yi))) / 6. / a
            return volmdlr.Point2D(cx, cy)

        else:
            self.plot()
            raise NotImplementedError

    def barycenter(self):
        """
        calculates the geometric center of the polygon, which is the
        average position of all the points in it

        returns a Volmdlr.Point2D point
        """
        barycenter1_2d = self.points[0]
        for point in self.points[1:]:
            barycenter1_2d += point
        return barycenter1_2d / len(self.points)

    def point_belongs(self, point):
        """
        Ray casting algorithm copied from internet...
        """
        return polygon_point_belongs((point.x, point.y),
                                     [(p.x, p.y) for p in self.points])

    def second_moment_area(self, point):
        Ix, Iy, Ixy = 0., 0., 0.
        for pi, pj in zip(self.points, self.points[1:] + [self.points[0]]):
            xi, yi = (pi - point)
            xj, yj = (pj - point)
            Ix += (yi ** 2 + yi * yj + yj ** 2) * (xi * yj - xj * yi)
            Iy += (xi ** 2 + xi * xj + xj ** 2) * (xi * yj - xj * yi)
            Ixy += (xi * yj + 2 * xi * yi + 2 * xj * yj + xj * yi) * (
                    xi * yj - xj * yi)
        if Ix < 0:
            Ix = - Ix
            Iy = - Iy
            Ixy = - Ixy
        return Ix / 12., Iy / 12., Ixy / 24.

    def _line_segments(self):
        lines = []
        if len(self.points) > 1:
            for p1, p2 in zip(self.points,
                              list(self.points[1:]) + [self.points[0]]):
                if p1 != p2:
                    lines.append(volmdlr.edges.LineSegment2D(p1, p2))
        return lines

    def rotation(self, center, angle, copy=True):
        if copy:
            return ClosedPolygon2D(
                [p.rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.rotation(center, angle, copy=False)

    # @classmethod
    # def polygon_from_segments(cls, list_point_pairs):
    #     points = [list_point_pairs[0][0], list_point_pairs[0][1]]
    #     list_point_pairs.remove((list_point_pairs[0][0], list_point_pairs[0][1]))
    #     finished =  False

    #     while not finished:
    #         for p1, p2 in list_point_pairs:
    #             if p1 == points[-1]:
    #                 points.append(p2)
    #                 break
    #             elif p2 == points[-1]:
    #                 points.append(p1)
    #                 break
    #         list_point_pairs.remove((p1, p2))
    #         if len(list_point_pairs)==0:
    #             finished = True

    #     # for i, i_p1, i_p2 in enumerate(list_point_pairs):
    #     #     for j, j_p1, j_p2 in enumerate(list_point_pairs):
    #     #         if i != j:

    #     #             if p1 == points[-1]:
    #     #                 points.append(p2)
    #     #             elif p2 == points[-1]:
    #     #                 points.append(p1)
    #     # print('points : ', points)
    #     return cls(points)

    def translation(self, offset, copy=True):
        if copy:
            return ClosedPolygon2D(
                [p.translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.translation(offset, copy=False)

    def polygon_distance(self,
                         polygon: 'volmdlr.wires.ClosedPolygon2D'):
        p = self.points[0]
        d = []
        for point in polygon.points:
            d.append(p.point_distance(point))
        index = d.index(min(d))
        return d[index]

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
        xmin, xmax, ymin, ymax = self.bounding_rectangle()

        max_offset_len = min(xmax - xmin, ymax - ymin) / 2
        if offset <= -max_offset_len:
            print('Inadapted offset, '
                  'polygon might turn over. Offset must be greater than',
                  -max_offset_len)
            raise ValueError('inadapted offset')
        else:
            nb = len(self.points)
            vectors = []
            for i in range(nb - 1):
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

        for i in range(nb):

            # check = False
            ni = vectors[2 * i - 1] + vectors[2 * i]
            if ni == volmdlr.Vector2D(0, 0):
                ni = vectors[2 * i]
                ni = ni.normal_vector()
                offset_vectors.append(ni)
            else:
                ni.normalize()
                if ni.dot(vectors[2 * i - 1].normal_vector()) > 0:
                    ni = - ni
                    # check = True
                offset_vectors.append(ni)

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
        Compute the distance to the border distance of polygon
        Output is always positive, even if the point belongs to the polygon
        """
        d_min, other_point_min = self.line_segments[0].point_distance(
            point, return_other_point=True)
        for line in self.line_segments[1:]:
            d, other_point = line.point_distance(
                point, return_other_point=True)
            if d < d_min:
                d_min = d
                other_point_min = other_point
        if return_other_point:
            return d_min, other_point_min
        return d_min

    def to_polygon(self, angle_resolution=None):
        return self

    def self_intersects(self):
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

                        line1 = volmdlr.edges.LineSegment2D(
                            self.points[segment1[0]],
                            self.points[segment1[1]])
                        line2 = volmdlr.edges.LineSegment2D(
                            self.points[segment2[0]],
                            self.points[segment2[1]])

                        p, a, b = volmdlr.Point2D.line_intersection(line1,
                                                                    line2,
                                                                    True)
                        if p is not None:
                            if 0 + epsilon <= a <= 1 - epsilon \
                                    and 0 + epsilon <= b <= 1 - epsilon:
                                return True, line1, line2

        return False, None, None

    @classmethod
    def points_convex_hull(cls, points):
        if len(points) < 3:
            return
        ymax, pos_ymax = volmdlr.core.max_pos([pt.y for pt in points])
        point_start = points[pos_ymax]
        hull = [point_start]

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
            theta_i = -volmdlr.core.clockwise_angle(vec1, vec2)
            theta.append(theta_i)

        min_theta, posmin_theta = volmdlr.core.min_pos(theta)
        next_point = remaining_points[posmin_theta]
        hull.append(next_point)
        del remaining_points[posmin_theta]
        # Adding first point to close the loop at the end
        remaining_points.append(hull[0])

        initial_vector = vec1.copy()
        total_angle = 0
        while next_point != point_start:
            vec1 = next_point - hull[-2]
            theta = []
            for pt in remaining_points:
                vec2 = pt - next_point
                theta_i = -volmdlr.core.clockwise_angle(vec1, vec2)
                theta.append(theta_i)

            min_theta, posmin_theta = volmdlr.core.min_pos(theta)
            if math.isclose(min_theta, -2 * math.pi, abs_tol=1e-6) \
                    or math.isclose(min_theta, 0, abs_tol=1e-6):
                if remaining_points[posmin_theta] == point_start:
                    break

            else:
                next_point = remaining_points[posmin_theta]

                vec_next_point = next_point - barycenter
                total_angle += (2 * math.pi - volmdlr.core.clockwise_angle(initial_vector, vec_next_point))

                if total_angle > 2 * math.pi:
                    break
                else:
                    initial_vector = vec_next_point

                hull.append(next_point)

            del remaining_points[posmin_theta]

        hull.pop()

        return cls(hull)

    @classmethod
    def concave_hull(cls, points, concavity, scale_factor):
        """
        Calculates the concave hull from a cloud of points, i.e., it Unites all points under the smallest possible area.

        :param points: list of points corresponding to the cloud of points
        :type points: class: 'volmdlr.Point2D'
        :param concavity: Sets how sharp the concave angles can be. It goes from -1 (not concave at all. in fact,
                          the hull will be left convex) up to +1 (very sharp angles can occur. Setting concavity to +1 might
                          result in 0º angles!) concavity is defined as the cosine of the concave angles.
        :type concavity: float
        :param scale_factor: Sets how big is the area where concavities are going to be searched.
                             The bigger, the more sharp the angles can be. Setting it to a very high value might affect the performance of the program.
                             This value should be relative to how close to each other the points to be connected are.
        :type scale_factor: float

        """

        def get_nearby_points(line, points, scale_factor):
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
                for point in points:
                    if not ((
                                    point.x == line.start.x and point.y == line.start.y) or (
                                    point.x == line.end.x and point.y == line.end.y)):
                        point_x_rel_pos = int(point.x / scale_factor)
                        point_y_rel_pos = int(point.y / scale_factor)
                        if point_x_rel_pos >= boundary[
                                0] and point_x_rel_pos <= boundary[
                                1] and point_y_rel_pos >= boundary[
                                2] and point_y_rel_pos <= boundary[3]:
                            nearby_points.append(point)

                scale_factor *= 4 / 3
                tries += 1

            return nearby_points

        def line_colides_with_hull(line, concave_hull):
            for hull_line in concave_hull:
                if line.start != hull_line.start and line.start != hull_line.end and line.end != hull_line.start and line.end != hull_line.end:
                    if line.line_intersections(hull_line):
                        return True
            return False

        def get_divided_line(line, nearby_points, hull_concave_edges,
                             concavity):
            divided_line = []
            ok_middle_points = []
            list_cossines = []
            for middle_point in nearby_points:
                vect1 = line.start - middle_point
                vect2 = line.end - middle_point
                if middle_point in (line.start, line.end):
                    continue
                cos = round(vect1.dot(vect2) / (vect1.norm() * vect2.norm()),
                            4)
                if cos < concavity:
                    new_lineA = volmdlr.edges.LineSegment2D(start=line.start,
                                                            end=middle_point)
                    new_lineB = volmdlr.edges.LineSegment2D(start=middle_point,
                                                            end=line.end)
                    if not (line_colides_with_hull(line=new_lineA,
                                                   concave_hull=hull_concave_edges) and line_colides_with_hull(
                            line=new_lineB, concave_hull=hull_concave_edges)):
                        ok_middle_points.append(middle_point)
                        list_cossines.append(cos)
            if len(ok_middle_points) > 0:
                #  We want the middlepoint to be the one with widest angle (smallest cossine)
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
        hull_points = list(set(
            [pt for line in hull_concave_edges for pt in [line[0], line[1]]]))
        unused_points = []
        for point in points:
            if point not in hull_points:
                unused_points.append(point)

        aLineWasDividedInTheIteration = True
        while aLineWasDividedInTheIteration:
            aLineWasDividedInTheIteration = False
            for line_position_hull in range(len(hull_concave_edges)):

                line = hull_concave_edges[line_position_hull]
                nearby_points = get_nearby_points(line, unused_points,
                                                  scale_factor)
                divided_line = get_divided_line(line, nearby_points,
                                                hull_concave_edges, concavity)
                if len(divided_line) > 0:
                    aLineWasDividedInTheIteration = True
                    unused_points.remove(divided_line[0].end)
                    hull_concave_edges.remove(line)
                    hull_concave_edges.extend(divided_line)
                    break

            hull_concave_edges.sort(key=lambda x: x.length(), reverse=True)

        # line  = hull_concave_edges[0]
        # print('first line legth :', line.length())
        # nearby_points = get_nearby_points(line, unused_points, scale_factor)
        # print('points next the first line in the end: ', nearby_points)
        # divided_line = get_divided_line(line, nearby_points, hull_concave_edges, concavity)
        # print('len divided line :', len(divided_line))
        polygon_points = [(line.start, line.end) for line in hull_concave_edges]
        # polygon_points = [(line.start, line.end) for line in hull_concave_edges
        #                   if line.length() != 0]

        points = [polygon_points[0][0], polygon_points[0][1]]
        polygon_points.remove((polygon_points[0][0], polygon_points[0][1]))
        finished = False

        while not finished:
            for p1, p2 in polygon_points:
                if p1 == points[-1] and p2 not in points:
                    points.append(p2)
                    break
                elif p2 == points[-1] and p1 not in points:
                    points.append(p1)
                    break
            polygon_points.remove((p1, p2))
            if len(polygon_points) == 0:
                finished = True

        return cls(points)  # , nearby_points

    @classmethod
    def convex_hull_points(cls, points):
        """
        Uses the scipy method ConvexHull to calculate the convex hull from
        a cloud of points
        """
        numpy_points = npy.array([(p.x, p.y) for p in points])
        hull = ConvexHull(numpy_points)
        polygon_points = []
        for simplex in hull.simplices:
            polygon_points.append((points[simplex[0]], points[simplex[1]]))

        points = [polygon_points[0][0], polygon_points[0][1]]
        polygon_points.remove((polygon_points[0][0], polygon_points[0][1]))
        finished = False

        while not finished:
            for p1, p2 in polygon_points:
                if p1 == points[-1]:
                    points.append(p2)
                    break
                elif p2 == points[-1]:
                    points.append(p1)
                    break
            polygon_points.remove((p1, p2))
            if len(polygon_points) == 0:
                finished = True

        if points[0] == points[-1]:
            return cls(points[:-1])

        return cls(points)

    def to_3d(self, plane_origin, x, y):
        points3d = [point.to_3d(plane_origin, x, y) for point in self.points]
        return ClosedPolygon3D(points3d)

    def plot(self, ax=None, color='k', alpha=1,
             plot_points=False, point_numbering=False,
             fill=False, fill_color='w', equal_aspect=True):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')

        if fill:
            ax.fill([p[0] for p in self.points], [p[1] for p in self.points],
                    facecolor=fill_color)
        for ls in self.line_segments:
            ls.plot(ax=ax, color=color, alpha=alpha)

        if plot_points or point_numbering:
            for point in self.points:
                point.plot(ax=ax, color=color, alpha=alpha)

        if point_numbering:
            for ip, point in enumerate(self.points):
                ax.text(*point, 'point {}'.format(ip + 1),
                        ha='center', va='top')

        if equal_aspect:
            ax.set_aspect('equal')
        else:
            ax.set_aspect('auto')

        ax.margins(0.1)
        plt.show()

        return ax

    def triangulation(self):
        """
        Note: triangles have been inverted for a better rendering in babylonjs
        """
        # ear clipping
        points = self.points[:]
        initial_point_to_index = {p: i for i, p in enumerate(self.points)}
        triangles = []

        remaining_points = self.points[:]
        # ax = ClosedPolygon2D(remaining_points).plot()

        # inital_number_points = len(remaining_points)
        number_remaining_points = len(remaining_points)
        while number_remaining_points > 3:
            current_polygon = ClosedPolygon2D(remaining_points)
            # print('remaining_points')
            # print(len(remaining_points))
            # pl2 = ClosedPolygon2D(remaining_points[1:]+remaining_points[0:1])
            # pl3 = ClosedPolygon2D(remaining_points[2:]+remaining_points[0:2])
            # current_polygon.plot(ax = ax)
            # pl2.plot(point_numbering=True)
            # pl3.plot(point_numbering=True)

            found_ear = False
            for p1, p2, p3 in zip(remaining_points,
                                  remaining_points[1:] + remaining_points[0:1],
                                  remaining_points[2:] + remaining_points[
                                                         0:2]):
                # ax.text(*p2, '{}')
                # ax = current_polygon.plot(point_numbering=True)
                if p1 != p3:
                    line_segment = volmdlr.edges.LineSegment2D(p1, p3)
                # line_segment.plot(color='grey', ax=ax)

                # ax2 = p1.plot(color='r')
                # p2.plot(color='g', ax=ax2)
                # p3.plot(color='b', ax=ax2)

                # Checking if intersections does not contrain the verticies
                # of line_segment
                intersect = False
                intersections = current_polygon.linesegment_intersections(
                    line_segment)
                if intersections:
                    for inter in intersections:
                        if inter[0] not in [line_segment.start,
                                            line_segment.end]:
                            intersect = True
                            break

                if not intersect:
                    # if not current_polygon.linesegment_intersections(line_segment):
                    # May be an ear
                    # print('ear?')
                    # if current_polygon.point_belongs(line_segment.middle_point()):
                    #     line_segment.middle_point().plot(color='g', ax=ax)
                    # else:
                    #     line_segment.middle_point().plot(color='r', ax=ax)
                    # print(current_polygon.point_belongs(
                    #         line_segment.middle_point()))

                    if current_polygon.point_belongs(
                            line_segment.middle_point()):
                        # Confirmed as an ear
                        # print('ear!')

                        triangles.append((initial_point_to_index[p1],
                                          initial_point_to_index[p3],
                                          initial_point_to_index[p2]))
                        remaining_points.remove(p2)
                        # ax.text(*points[initial_point_to_index[p2]], str(number_remaining_points))
                        number_remaining_points -= 1
                        found_ear = True
                        break

            if not found_ear:
                remaining_polygon = ClosedPolygon2D(remaining_points)
                if remaining_polygon.area() > 0.:
                    # Searching for a flat ear
                    found_flat_ear = False
                    for p1, p2, p3 in zip(remaining_points,
                                          remaining_points[
                                              1:] + remaining_points[0:1],
                                          remaining_points[
                                              2:] + remaining_points[0:2]):
                        triangle = Triangle2D(p1, p2, p3)
                        if triangle.area() == 0:
                            remaining_points.remove(p2)
                            found_flat_ear = True
                            break

                    if not found_flat_ear:
                        # remaining_polygon.plot(point_numbering=True, plot_points=True)
                        # vmd.DisplayMesh2D(points, triangles).plot()
                        # print(remaining_points)
                        # raise ValueError('There are no ear in the polygon, it seems malformed')
                        print(
                            'Warning : There are no ear in the polygon, it seems malformed: skipping triangulation')
                        return vmd.DisplayMesh2D(points, triangles)
                else:
                    return vmd.DisplayMesh2D(points, triangles)

        if len(remaining_points) == 3:
            p1, p2, p3 = remaining_points
            triangles.append((initial_point_to_index[p1],
                              initial_point_to_index[p3],
                              initial_point_to_index[p2]))

        return vmd.DisplayMesh2D(points, triangles)

    def simplify(self, min_distance: float = 0.01, max_distance: float = 0.05):
        return ClosedPolygon2D(self.simplify_polygon(min_distance=min_distance,
                                                     max_distance=max_distance).points)

    def line_intersecting_closing_point(self, crossing_point):
        """
        finds closing point for the sewing method using intersection of lines
        drawn from the barycenter
        returns the closing point
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
        In case the barycenter of the polygon is outside, this method
        finds another point inside the polygon
        """
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
                        translation1 = (list(intersetions1.values())[
                                            i + 1][0] +
                                        list(intersetions1.values())[
                                            i + 1][1]) * 0.5
                        break
                if i % 2 != 0:
                    if len(list(intersetions1.values())[i - 1]) == 2:
                        translation1 = (list(intersetions1.values())[
                                            i - 1][0] +
                                        list(intersetions1.values())[
                                            i - 1][1]) * 0.5
                        break

        return translation1

    def repositioned_polygon(self, x, y):
        linex = volmdlr.edges.LineSegment2D(-x.to_2d(volmdlr.O2D, x, y),
                                            x.to_2d(volmdlr.O2D, x, y))
        way_back = volmdlr.O3D
        barycenter = self.barycenter()
        if not self.point_belongs(barycenter):
            barycenter1_2d = self.point_in_polygon()
            self.translation(-barycenter1_2d, False)
            way_back = barycenter1_2d.to_3d(volmdlr.O3D, x, y)
        else:
            inters = self.linesegment_intersections(linex)
            distance = inters[0][0].point_distance(inters[-1][0])
            if distance / 2 > 3 * min(
                    self.point_distance(inters[0][0]),
                    self.point_distance(inters[-1][0])):
                mid_point = (inters[0][0] + inters[-1][0]) * 0.5
                self.translation(-mid_point, False)
                way_back = mid_point.to_3d(volmdlr.O3D, x, y)

        return self, way_back

    def get_closing_point(self, polygon2_2d, primitive, ax=None):
        closing_point = volmdlr.O2D
        middle_point = primitive.middle_point()

        normal_vector = primitive.unit_normal_vector()
        line_segment1 = volmdlr.edges.LineSegment2D(
            middle_point - normal_vector, middle_point)
        line_segment2 = volmdlr.edges.LineSegment2D(middle_point,
                                                    middle_point + normal_vector)

        # print('angle in radians :', primitive.unit_direction_vector().dot(normal_vector))
        # print('vector1 :', primitive.unit_direction_vector())
        # print('normal: ', normal_vector)
        line_intersections = {line_segment1: [], line_segment2: []}
        for ls in [line_segment1, line_segment2]:
            inter_points = []
            for prim in polygon2_2d.line_segments + self.line_segments[
                                                    :self.line_segments.index(
                                                            primitive)] + self.line_segments[
                                                                           self.line_segments.index(
                                                                                   primitive) + 1:]:
                inters = prim.linesegment_intersections(ls)
                if inters:
                    line_intersections[ls].append((inters[0], prim))
                    inter_points.append(inters[0])
                elif ls.point_belongs(prim.start):
                    if prim.start not in inter_points:
                        line_intersections[ls].append((prim.start, prim))
                        inter_points.append(prim.start)
                elif ls.point_belongs(prim.end):
                    if prim.end not in inter_points:
                        line_intersections[ls].append((prim.end, prim))
                        inter_points.append(prim.end)
                elif prim.point_belongs(middle_point):
                    line_intersections[ls].append((prim.middle_point(), prim))
                    inter_points.append(prim.middle_point())
        if line_intersections[line_segment1]:
            # print('passing here, debug me')
            # inters_points = [intr_list[0] for intr_list in line_intersections[line_segment1]]
            distance = math.inf
            for intr_list in line_intersections[line_segment1]:
                # if intr_list[0] not in self.points:
                if intr_list[1] not in self.line_segments:
                    dist = intr_list[0].point_distance(line_segment1.start)
                    if dist < distance:
                        distance = dist
                        if intr_list[0].point_distance(intr_list[1].start) < \
                                intr_list[0].point_distance(intr_list[1].end):
                            closing_point = intr_list[1].start
                        else:
                            closing_point = intr_list[1].end
                elif intr_list[0] == middle_point and \
                        primitive.length() == intr_list[1].length():
                    closing_point = intr_list[1].start
        else:
            # print('last passed here, sure of it ')
            distance = math.inf
            for intr_list in line_intersections[line_segment2]:
                if intr_list[1] not in self.line_segments:
                    dist = intr_list[0].point_distance(line_segment1.start)
                    if dist < distance:
                        distance = dist
                        if intr_list[0].point_distance(intr_list[1].start) < \
                                intr_list[0].point_distance(intr_list[1].end):
                            closing_point = intr_list[1].start
                        else:
                            closing_point = intr_list[1].end
                elif intr_list[0] == middle_point and \
                        primitive.length() == intr_list[1].length():
                    closing_point = intr_list[1].start
        if ax is not None:
            middle_point.plot(ax=ax, color='r')
            line_segment1.plot(ax=ax, color='y')
            line_segment2.plot(ax=ax, color='b')
            closing_point.plot(ax=ax)

        return closing_point

    def validate_sewing_polygons(self, polygon2_2d):
        primitive1 = self.line_segments[0]
        closing_point = self.get_closing_point(polygon2_2d,
                                               primitive1)
        if closing_point != volmdlr.O2D:
            return True
        return False

    def get_valid_sewing_polygons(self, polygon2_2d):
        for i, primitive1 in enumerate(self.line_segments):
            closing_point = self.get_closing_point(polygon2_2d,
                                                   primitive1)
            if closing_point != volmdlr.O2D:
                if i == 0:
                    return self
                new_polygon_primitives = \
                    self.line_segments[self.line_segments.index(primitive1):] + \
                    self.line_segments[:self.line_segments.index(primitive1)]
                new_points = []
                for prim in new_polygon_primitives:
                    if prim.start not in new_points:
                        new_points.append(prim.start)
                    if prim.end not in new_points:
                        new_points.append(prim.end)
                return ClosedPolygon2D(new_points)
        raise NotImplementedError('make sure the two polygons '
                                  'you are trying to sew are valid ones')


class Triangle2D(ClosedPolygon2D):
    def __init__(self, point1: volmdlr.Point2D, point2: volmdlr.Point2D,
                 point3: volmdlr.Point2D, name: str = ''):
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
        self.name = name

        # ClosedPolygon2D.__init__(self, points=[point1, point2, point3],
        # name=name)

    def area(self):
        u = self.point2 - self.point1
        v = self.point3 - self.point1
        return abs(u.cross(v)) / 2

    def incircle_radius(self):
        a = self.point1.point_distance(self.point2)
        b = self.point1.point_distance(self.point3)
        c = self.point2.point_distance(self.point3)
        return 2 * self.area() / (a + b + c)

    def circumcircle_radius(self):
        a = self.point1.point_distance(self.point2)
        b = self.point1.point_distance(self.point3)
        c = self.point2.point_distance(self.point3)
        return a * b * c / (self.area() * 4.0)

    def ratio_circumr_length(self):
        return self.circumcircle_radius() / self.length()

    def ratio_incircler_length(self):
        return self.incircle_radius() / self.length()

    def aspect_ratio(self):
        a = self.point1.point_distance(self.point2)
        b = self.point1.point_distance(self.point3)
        c = self.point2.point_distance(self.point3)
        s = 0.5 * (a + b + c)
        try:
            return 0.125 * a * b * c / (s - a) / (s - b) / (s - c)
        except ZeroDivisionError:
            return 1000000.


class Circle2D(Contour2D):
    _non_serializable_attributes = ['internal_arcs', 'external_arcs',
                                    'polygon', 'straight_line_contour_polygon',
                                    'primitives', 'basis_primitives']

    def __init__(self, center: volmdlr.Point2D, radius: float, name: str = ''):
        self.center = center
        self.radius = radius
        self.angle = volmdlr.TWO_PI

        # self.points = self.tessellation_points()

        Contour2D.__init__(self, [self], name=name)  # !!! this is dangerous

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

    def to_polygon(self, angle_resolution: float):
        return ClosedPolygon2D(
            self.polygon_points(angle_resolution=angle_resolution))

    def tessellation_points(self, resolution=40):
        return [(self.center
                 + self.radius * math.cos(teta) * volmdlr.X2D
                 + self.radius * math.sin(teta) * volmdlr.Y2D)
                for teta in npy.linspace(0, volmdlr.TWO_PI, resolution + 1)][
               :-1]

    def point_belongs(self, point, tolerance=1e-9):
        return point.point_distance(self.center) <= self.radius + tolerance

    # def border_points(self):
    #     start = self.center - self.radius * volmdlr.Point2D(1, 0)
    #     end = self.center + self.radius * volmdlr.Point2D(1, 0)
    #     return [start, end]

    def bounding_rectangle(self):

        xmin = self.center.x - self.radius
        xmax = self.center.x + self.radius
        ymin = self.center.y - self.radius
        ymax = self.center.y + self.radius
        return xmin, xmax, ymin, ymax

    def line_intersections(self, line2d: volmdlr.edges.Line2D, tol=1e-9):
        # Duplicate from ffull arc
        Q = self.center
        if line2d.points[0] == self.center:
            P1 = line2d.points[1]
            V = line2d.points[0] - line2d.points[1]
        else:
            P1 = line2d.points[0]
            V = line2d.points[1] - line2d.points[0]
        a = V.dot(V)
        b = 2 * V.dot(P1 - Q)
        c = P1.dot(P1) + Q.dot(Q) - 2 * P1.dot(Q) - self.radius ** 2

        disc = b ** 2 - 4 * a * c
        if math.isclose(disc, 0., abs_tol=tol):
            t1 = -b / (2 * a)
            return [P1 + t1 * V]

        elif disc > 0:
            sqrt_disc = math.sqrt(disc)
            t1 = (-b + sqrt_disc) / (2 * a)
            t2 = (-b - sqrt_disc) / (2 * a)
            return [P1 + t1 * V,
                    P1 + t2 * V]
        else:
            return []

    def circle_intersections(self, circle: 'volmdlr.wires.Circle2D'):
        x0, y0 = self.center
        x1, y1 = circle.center
        r0 = self.radius
        r1 = circle.radius

        d = math.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)

        # non intersecting
        if d > r0 + r1:
            return []
        # One circle within other
        if d < abs(r0 - r1):
            return []
        # coincident circles
        if d == 0 and r0 == r1:
            return []
        else:
            a = (r0 ** 2 - r1 ** 2 + d ** 2) / (2 * d)
            h = math.sqrt(r0 ** 2 - a ** 2)
            x2 = x0 + a * (x1 - x0) / d
            y2 = y0 + a * (y1 - y0) / d
            x3 = x2 + h * (y1 - y0) / d
            y3 = y2 - h * (x1 - x0) / d

            x4 = x2 - h * (y1 - y0) / d
            y4 = y2 + h * (x1 - x0) / d

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
        return volmdlr.TWO_PI * self.radius

    def plot(self, ax=None, linestyle='-', color='k', linewidth=1, alpha=1.,
             equal_aspect=True):
        if ax is None:
            fig, ax = plt.subplots()
        # else:
        #     fig = ax.figure
        if self.radius > 0:
            ax.add_patch(matplotlib.patches.Arc((self.center.x, self.center.y),
                                                2 * self.radius,
                                                2 * self.radius,
                                                angle=0,
                                                theta1=0,
                                                theta2=360,
                                                color=color,
                                                alpha=alpha,
                                                linestyle=linestyle,
                                                linewidth=linewidth))
        if equal_aspect:
            ax.set_aspect('equal')
        return ax

    def to_3d(self, plane_origin, x, y):
        normal = x.cross(y)
        center3d = self.center.to_3d(plane_origin, x, y)
        return Circle3D(volmdlr.Frame3D(center3d, x, y, normal),
                        self.radius, self.name)

    def rotation(self, center, angle, copy=True):
        if copy:
            return Circle2D(self.center.rotation(center, angle, copy=True),
                            self.radius)
        else:
            self.center.rotation(center, angle, copy=False)

    def translation(self, offset, copy=True):
        if copy:
            return Circle2D(self.center.translation(offset, copy=True),
                            self.radius)
        else:
            self.center.translation(offset, copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return Circle2D(frame.old_coordinates(self.center),
                                self.radius)
            else:
                self.center = frame.old_coordinates(self.center)
        if side == 'new':
            if copy:
                return Circle2D(frame.new_coordinates(self.center),
                                self.radius)
            else:
                self.points = frame.new_coordinates(self.center)

    def area(self):
        return math.pi * self.radius ** 2

    def second_moment_area(self, point):
        """
        Second moment area of part of disk
        """
        I = math.pi * self.radius ** 4 / 4
        return volmdlr.geometry.huygens2d(I, I, 0, self.area(), self.center,
                                          point)

    def center_of_mass(self):
        return self.center

    def point_symmetric(self, point):
        center = 2 * point - self.center
        return Circle2D(center, self.radius)

    def plot_data(self, edge_style: plot_data.EdgeStyle = None,
                  surface_style: plot_data.SurfaceStyle = None):
        return plot_data.Circle2D(cx=self.center.x,
                                  cy=self.center.y,
                                  r=self.radius,
                                  edge_style=edge_style,
                                  surface_style=surface_style)

    def copy(self, *args, **kwargs):
        return Circle2D(self.center.copy(), self.radius)

    def point_at_abscissa(self, curvilinear_abscissa):
        start = self.center + self.radius * volmdlr.X3D
        return start.rotation(self.center,
                              curvilinear_abscissa / self.radius)

    # def triangulation(self, n=35):
    #     l = self.length()
    #     points = [self.point_at_abscissa(l * i / n) for i in range(n)]
    #     points.append(self.center)
    #     triangles = [(i, i + 1, n) for i in range(n - 1)] + [(n - 1, 0, n)]

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

    def discretise(self, n: float):
        # BUGGED: returns method
        circle_to_nodes = {}
        nodes = []
        if n * self.length() < 1:
            circle_to_nodes[self] = self.border_points
        else:
            n0 = int(math.ceil(n * self.length()))
            l0 = self.length() / n0

            for k in range(n0):
                node = self.point_at_abscissa(k * l0)

                nodes.append(node)

            circle_to_nodes[self] = nodes

        return circle_to_nodes[self]

    def polygon_points(self, angle_resolution=10):
        return volmdlr.edges.Arc2D.polygon_points(
            self, angle_resolution=angle_resolution)


class Contour3D(Contour, Wire3D):
    _non_serializable_attributes = ['points']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['points', 'name']
    _generic_eq = True
    """
    A collection of 3D primitives forming a closed wire3D
    """

    def __init__(self, primitives: List[volmdlr.core.Primitive3D],
                 name: str = ''):
        """

        """

        Wire3D.__init__(self, primitives=primitives, name=name)
        self._utd_edge_polygon = False

    def __hash__(self):
        return sum([hash(e) for e in self.primitives])

    def __eq__(self, other_):
        if self.__class__.__name__ != other_.__class__.__name__:
            return False
        equal = True
        for edge, other_edge in zip(self.primitives, other_.edges):
            equal = (equal and edge == other_edge)
        return equal

    @property
    def edge_polygon(self):
        if not self._utd_edge_polygon:
            self._edge_polygon = self._get_edge_polygon()
            self._utd_edge_polygon = True
        return self._edge_polygon

    def _get_edge_polygon(self):
        points = []
        for edge in self.primitives:
            if points:
                if edge.start != points[-1]:
                    points.append(edge.start)
            else:
                points.append(edge.start)
        return ClosedPolygon3D(points)

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]
        raw_edges = []
        # edge_ends = {}
        for ie, edge_id in enumerate(arguments[1]):
            edge = object_dict[int(edge_id[1:])]
            raw_edges.append(edge)

        if (len(raw_edges)) == 1:
            if isinstance(raw_edges[0], cls):
                # Case of a circle, ellipse...
                return raw_edges[0]
            else:
                return cls(raw_edges, name=name)

        distances = [raw_edges[0].end.point_distance(raw_edges[1].start),
                     raw_edges[0].start.point_distance(raw_edges[1].start),
                     raw_edges[0].end.point_distance(raw_edges[1].end),
                     raw_edges[0].start.point_distance(raw_edges[1].end)]
        index = distances.index(min(distances))
        if min(distances) > 5e-4:
            ax = raw_edges[0].plot()
            raw_edges[1].plot(ax=ax)
            deltax1 = abs(raw_edges[0].start.x - raw_edges[1].end.x)
            deltax2 = abs(raw_edges[0].end.x - raw_edges[1].end.x)
            deltay1 = abs(raw_edges[0].start.y - raw_edges[1].end.y)
            deltay2 = abs(raw_edges[0].end.y - raw_edges[1].end.y)
            deltaz1 = abs(raw_edges[0].start.z - raw_edges[1].end.z)
            deltaz2 = abs(raw_edges[0].end.z - raw_edges[1].end.z)
            raise NotImplementedError(
                'First 2 edges of contour not follwing each other',
                'delta = {}, {}, {}, {}, {}, {}'.format(deltax1, deltax2,
                                                        deltay1, deltay2,
                                                        deltaz1, deltaz2))
        # Making things right for first 2 primitives

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

        last_edge = edges[-1]
        for raw_edge in raw_edges[2:]:
            distances = [raw_edge.start.point_distance(last_edge.end),
                         raw_edge.end.point_distance(last_edge.end)]
            index = distances.index(min(distances))
            if min(distances) > 5e-4:
                ax = last_edge.plot(color='b')
                raw_edge.plot(ax=ax, color='r')
                last_edge.end.plot(ax=ax, color='b')
                raw_edges[0].plot(ax=ax, color='g')
                deltax1 = abs(raw_edge.start.x - last_edge.end.x)
                deltax2 = abs(raw_edge.end.x - last_edge.end.x)
                deltay1 = abs(raw_edge.start.y - last_edge.end.y)
                deltay2 = abs(raw_edge.end.y - last_edge.end.y)
                deltaz1 = abs(raw_edge.start.z - last_edge.end.z)
                deltaz2 = abs(raw_edge.end.z - last_edge.end.z)
                raise NotImplementedError(
                    'Edges of contour not follwing each other',
                    'delta = {}, {}, {}, {}, {}, {}'.format(deltax1, deltax2,
                                                            deltay1, deltay2,
                                                            deltaz1, deltaz2))
            if index == 0:
                last_edge = raw_edge
            elif index == 1:
                last_edge = raw_edge.reverse()

            edges.append(last_edge)
        return cls(edges, name=name)

    def to_step(self, current_id, surface_id=None):

        content = ''
        edge_ids = []
        for primitive in self.primitives:
            if isinstance(primitive, volmdlr.edges.BSplineCurve3D):
                continue
            primitive_content, primitive_ids = primitive.to_step(current_id)
            content += primitive_content
            current_id = primitive_ids[-1] + 1
            for primitive_id in primitive_ids:
                content += "#{} = ORIENTED_EDGE('{}',*,*,#{},.T.);\n".format(
                    current_id,
                    primitive.name,
                    primitive_id)
                edge_ids.append(current_id)

                current_id += 1

        content += "#{} = EDGE_LOOP('{}',({}));\n".format(
            current_id, self.name, volmdlr.core.step_ids_to_str(edge_ids))
        return content, current_id

    def average_center_point(self):
        nb = len(self.points)
        x = npy.sum([p[0] for p in self.points]) / nb
        y = npy.sum([p[1] for p in self.points]) / nb
        z = npy.sum([p[2] for p in self.points]) / nb

        return volmdlr.Point3D(x, y, z)

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            new_edges = [edge.rotation(center, axis, angle, copy=True) for edge
                         in self.primitives]
            # new_points = [p.rotation(center, axis, copy=True) for p in self.points]
            return Contour3D(new_edges, None, self.name)
        else:
            for edge in self.primitives:
                edge.rotation(center, axis, angle, copy=False)
            for point in self.tessel_points:
                point.rotation(center, axis, angle, copy=False)

    def translation(self, offset, copy=True):
        if copy:
            new_edges = [edge.translation(offset, copy=True) for edge in
                         self.primitives]
            # new_points = [p.translation(offset, copy=True) for p in self.points]
            return Contour3D(new_edges, self.name)
        else:
            for edge in self.primitives:
                edge.translation(offset, copy=False)
            for point in self.tessel_points:
                point.translation(offset, copy=False)

    def order_contour(self):
        # new_primitives = []
        # points = self.ordering_contour()
        # for p1, p2 in points:
        #     new_primitives.append(volmdlr.edges.LineSegment3D(p1, p2))
        # self.primitives = new_primitives

        initial_points = []
        for p in self.primitives:
            initial_points.append((p.start, p.end))

        new_primitives = []
        if self.is_ordered():
            return self
        points = self.ordering_contour()
        for p1, p2 in points:
            try:
                index = initial_points.index((p1, p2))
            except ValueError:
                index = initial_points.index((p2, p1))

            if isinstance(self.primitives[index], volmdlr.edges.LineSegment3D):
                new_primitives.append(volmdlr.edges.LineSegment3D(p1, p2))
            elif isinstance(self.primitives[index], volmdlr.edges.Arc3D):
                new_primitives.append(volmdlr.edges.Arc3D(p1, self.primitives[index].interior, p2))
            elif isinstance(self.primitives[index], volmdlr.edges.BSplineCurve3D):
                if (self.primitives[index].start == p1 and self.primitives[index].end == p2):
                    new_primitives.append(self.primitives[index])
                else:
                    new_primitives.append(self.primitives[index].reverse())

        self.primitives = new_primitives

        return self

    def point_over_contour(self, point, abs_tol=1e-7):
        belongs = False
        for primitive in self.primitives:
            if primitive.point_belongs(point, abs_tol):
                belongs = True
        return belongs

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

    def copy(self, deep=True, memo=None):
        new_edges = [edge.copy(deep=deep, memo=memo) for edge in self.primitives]
        if self.point_inside_contour is not None:
            new_point_inside_contour = self.point_inside_contour.copy()
        else:
            new_point_inside_contour = None
        return Contour3D(new_edges, new_point_inside_contour, self.name)

    def length(self):
        # TODO: this is duplicated code from Wire3D!
        length = 0.
        for edge in self.primitives:
            length += edge.length()
        return length

    def point_at_abscissa(self, curvilinear_abscissa):
        # TODO: this is duplicated code from Wire3D!
        length = 0.
        for primitive in self.primitives:
            primitive_length = primitive.length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.point_at_abscissa(
                    curvilinear_abscissa - length)
            length += primitive_length
        if math.isclose(curvilinear_abscissa, length, abs_tol=1e-6):
            return primitive.point_at_abscissa(primitive_length)
        raise ValueError('abscissa out of contour length')

    def plot(self, ax=None, color='k', alpha=1, edge_details=False):
        if ax is None:
            ax = Axes3D(plt.figure())

        for edge in self.primitives:
            edge.plot(ax=ax, color=color, alpha=alpha,
                      edge_ends=edge_details, edge_direction=edge_details)

        return ax

    def to_2d(self, plane_origin, x, y):
        z = x.cross(y)
        plane3d = volmdlr.faces.Plane3D(volmdlr.Frame3D(plane_origin, x, y, z))
        primitives2d = [plane3d.point3d_to_2d(p) for p in self.primitives]
        return Contour2D(primitives=primitives2d)

    def _bounding_box(self):
        """
        Flawed method, to be enforced by overloading
        """
        n = 50
        l = self.length()
        points = [self.point_at_abscissa(i / n * l)
                  for i in range(n)]
        return volmdlr.core.BoundingBox.from_points(points)

    @classmethod
    def extract_contours(cls, contour, point1: volmdlr.Point3D,
                         point2: volmdlr.Point3D, inside=False):

        new_primitives = contour.extract_with_points(point1, point2, inside)
        contours = [cls(new_primitives)]
        return contours

    def contour_intersection(self, contour3d):
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
        '''
        create a contour3d from points with line_segments3D
        '''

        if len(points) < 3:
            raise ValueError('contour is defined at least with three points')
        else:
            edges = []
            for i in range(0, len(points) - 1):
                edges.append(volmdlr.edges.LineSegment3D(points[i], points[i + 1]))

            edges.append(volmdlr.edges.LineSegment3D(points[len(points)], points[0]))

            contour = cls(edges)

            return contour

    def clean_primitives(self):
        '''
        delete primitives with start=end, and return a new contour
        '''

        new_primitives = []
        for p in self.primitives:
            if p.start != p.end:
                new_primitives.append(p)

        return Contour3D(new_primitives)

    def merge_with(self, contour3d):
        '''
        merge two adjacent contours, sharing primitives, and returns one outer contour and inner contours (if there are any)
        '''

        merged_primitives = self.merge_primitives_with(contour3d)
        contours = volmdlr.wires.Contour3D.contours_from_edges(merged_primitives, tol=3e-4)
        # contours = sorted(contours, key=lambda contour: contour.area(), reverse=True)

        return contours


class Circle3D(Contour3D):
    _non_serializable_attributes = ['point', 'edges', 'point_inside_contour']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, frame: volmdlr.Frame3D, radius: float,
                 name: str = ''):
        """
        frame.u, frame.v define the plane, frame.w the normal
        """
        self.radius = radius
        self.frame = frame
        self.angle = volmdlr.TWO_PI
        Contour3D.__init__(self, [self], name=name)

    @property
    def center(self):
        return self.frame.origin

    @property
    def normal(self):
        return self.frame.w

    def __hash__(self):
        return hash(self.frame.origin)

    def __eq__(self, other_circle):
        return self.frame.origin == other_circle.frame.origin \
               and self.frame.w.is_colinear(other_circle.frame.w) \
               and math.isclose(self.radius,
                                other_circle.radius, abs_tol=1e-06)

    def tessellation_points(self, resolution=20):

        tessellation_points_3d = [
                                     self.center + self.radius * math.cos(
                                         teta) * self.frame.u
                                     + self.radius * math.sin(
                                         teta) * self.frame.v
                                     for teta in
                                     npy.linspace(0, volmdlr.TWO_PI,
                                                  resolution + 1)][:-1]
        return tessellation_points_3d

    def length(self):
        return volmdlr.TWO_PI * self.radius

    def FreeCADExport(self, name, ndigits=3):
        xc, yc, zc = round(1000 * self.center, ndigits)
        xn, yn, zn = round(self.normal, ndigits)
        return '{} = Part.Circle(fc.Vector({},{},{}),fc.Vector({},{},{}),{})\n'.format(
            name, xc, yc, zc, xn, yn, zn, 1000 * self.radius)

    def rotation(self, rot_center, axis, angle, copy=True):
        new_center = self.center.rotation(rot_center, axis, angle, True)
        new_normal = self.normal.rotation(rot_center, axis, angle, True)
        if copy:
            return Circle3D(new_center, self.radius, new_normal, self.name)
        else:
            self.center = new_center
            self.normal = new_normal

    def translation(self, offset, copy=True):
        # new_frame = self.center.translation(offset, True)
        new_frame = self.frame.translation(offset, True)

        if copy:

            # return Circle3D(new_frame, self.radius, self.frame,
            #                 self.name)
            return Circle3D(new_frame, self.radius, self.name)

        else:
            self.frame = new_frame

    def plot(self, ax=None, color='k', alpha=1., edge_details=False):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        x = []
        y = []
        z = []
        for px, py, pz in self.tessellation_points():
            x.append(px)
            y.append(py)
            z.append(pz)
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])
        ax.plot(x, y, z, color=color, alpha=alpha)
        return ax

    def point_at_abscissa(self, curvilinear_abscissa):
        """
        start point is at intersection of frame.u axis
        """
        start = self.frame.origin + self.radius * self.frame.u
        return start.rotation(self.frame.origin, self.frame.w,
                              curvilinear_abscissa / self.radius,
                              copy=True)

    @classmethod
    def from_step(cls, arguments, object_dict):
        center = object_dict[arguments[1]].origin
        radius = float(arguments[2]) / 1000
        if object_dict[arguments[1]].u is not None:
            normal = object_dict[arguments[1]].u
            other_vec = object_dict[arguments[1]].v
            if other_vec is not None:
                other_vec.normalize()
        else:
            normal = object_dict[arguments[1]].v  # ou w
            other_vec = None
        normal.normalize()
        return cls.from_center_normal(center, normal, radius,
                                      arguments[0][1:-1])

    def to_step(self, current_id, surface_id=None):
        circle_frame = volmdlr.Frame3D(self.center, self.frame.w, self.frame.u,
                                       self.frame.v)
        content, frame_id = circle_frame.to_step(current_id)
        curve_id = frame_id + 1
        content += "#{} = CIRCLE('{}',#{},{});\n".format(
            curve_id, self.name, frame_id, round(self.radius * 1000, 3))

        if surface_id:
            content += "#{} = SURFACE_CURVE('',#{},(#{}),.PCURVE_S1.);\n".format(
                curve_id + 1, curve_id, surface_id)
            curve_id += 1

        p1 = self.frame.origin + self.frame.u * self.radius
        # p2 = self.frame.origin + self.frame.v*self.radius
        p3 = self.frame.origin - self.frame.u * self.radius
        # p4 = self.frame.origin - self.frame.v*self.radius

        p1_content, p1_id = p1.to_step(curve_id + 1, vertex=True)
        # p2_content, p2_id = p2.to_step(p1_id+1, vertex=True)
        p3_content, p3_id = p3.to_step(p1_id + 1, vertex=True)
        # p4_content, p4_id = p4.to_step(p3_id+1, vertex=True)
        content += p1_content + p3_content

        arc1_id = p3_id + 1
        content += "#{} = EDGE_CURVE('{}',#{},#{},#{},.T.);\n".format(
            arc1_id, self.name, p1_id, p3_id, curve_id)
        oriented_edge1_id = arc1_id + 1
        content += "#{} = ORIENTED_EDGE('',*,*,#{},.T.);\n".format(
            oriented_edge1_id, arc1_id)

        arc2_id = oriented_edge1_id + 1
        content += "#{} = EDGE_CURVE('{}',#{},#{},#{},.T.);\n".format(
            arc2_id, self.name, p3_id, p1_id, curve_id)
        oriented_edge2_id = arc2_id + 1
        content += "#{} = ORIENTED_EDGE('',*,*,#{},.T.);\n".format(
            oriented_edge2_id, arc2_id)

        current_id = oriented_edge2_id + 1
        content += "#{} = EDGE_LOOP('{}',(#{},#{}));\n".format(
            current_id, self.name, oriented_edge1_id, oriented_edge2_id)

        return content, current_id

    def _bounding_box(self):
        """
        """
        # u = self.normal.deterministic_unit_normal_vector()
        # v = self.normal.cross(u)
        points = [self.frame.origin + self.radius * v
                  for v in [self.frame.u, -self.frame.u,
                            self.frame.v, -self.frame.v]]
        return volmdlr.core.BoundingBox.from_points(points)

    def to_2d(self, plane_origin, x, y):
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
        u1 = (point2 - point1)
        u2 = (point2 - point3)
        try:
            u1.normalize()
            u2.normalize()
        except ZeroDivisionError:
            raise ValueError(
                'the 3 points must be distincts')

        normal = u2.cross(u1)
        normal.normalize()

        if u1 == u2:
            u2 = normal.cross(u1)
            u2.normalize()

        v1 = normal.cross(u1)  # v1 is normal, equal u2
        v2 = normal.cross(u2)  # equal -u1

        p11 = 0.5 * (point1 + point2)  # Mid point of segment s,m
        p21 = 0.5 * (point2 + point3)  # Mid point of segment s,m

        l1 = volmdlr.edges.Line3D(p11, p11 + v1)
        l2 = volmdlr.edges.Line3D(p21, p21 + v2)

        try:
            center, _ = l1.minimum_distance_points(l2)
        except ZeroDivisionError:
            raise ValueError(
                'Start, end and interior points  of an arc must be distincts')

        radius = (center - point1).norm()
        return cls(frame=volmdlr.Frame3D(center, u1, normal.cross(u1), normal),
                   radius=radius)

    def extrusion(self, extrusion_vector):

        if self.normal.is_colinear_to(extrusion_vector):
            u = self.normal.deterministic_unit_normal_vector()
            v = self.normal.cross(u)
            w = extrusion_vector.copy()
            w.normalize()
            cylinder = volmdlr.faces.CylindricalSurface3D(
                volmdlr.Frame3D(self.center, u, v, w), self.radius)
            return [cylinder.rectangular_cut(0, volmdlr.TWO_PI,
                                             0, extrusion_vector.norm())]
        else:
            raise NotImplementedError(
                'Extrusion along vector not colinar to normal for circle not handled yet: dot={}'.format(
                    self.normal.dot(extrusion_vector)))

    def revolution(self, axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                   angle: float):
        line3d = volmdlr.edges.Line3D(axis_point, axis_point + axis)
        tore_center, _ = line3d.point_projection(self.center)
        u = self.center - tore_center
        u.normalize()
        v = axis.cross(u)
        if not math.isclose(self.normal.dot(u), 0., abs_tol=1e-9):
            raise NotImplementedError(
                'Outside of plane revolution not supported')

        R = tore_center.point_distance(self.center)
        surface = volmdlr.faces.ToroidalSurface3D(
            volmdlr.Frame3D(tore_center, u, v, axis),
            R, self.radius)
        return [surface.rectangular_cut(0, angle, 0, volmdlr.TWO_PI)]

    def point_on_circle(self, point: volmdlr.Point3D):
        distance = point.point_distance(self.center)
        vec = volmdlr.Vector3D(*point - self.center)
        dot = self.normal.dot(vec)
        if math.isclose(distance, self.radius, abs_tol=1e-6)\
                and math.isclose(dot, 0, abs_tol=5e-6):
            return True
        return False

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D):
        if not self.point_on_circle(point1)\
                or not self.point_on_circle(point2):
            ax = self.plot()
            point1.plot(ax=ax, color='r')
            point2.plot(ax=ax, color='b')
            raise ValueError('Point not on circle for trim method')
        if point1 == point2:
            return volmdlr.edges.FullArc3D(self.frame.origin, point1,
                                           self.frame.w)
        interior = volmdlr.core.clockwise_interior_from_circle3d(
            point1, point2, self)
        return volmdlr.edges.Arc3D(point1, interior, point2)


class Ellipse3D(Contour3D):
    """
    :param major_axis: Largest radius of the ellipse
    :type major_axis: float
    :param minor_axis: Smallest radius of the ellipse
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
        Contour3D.__init__(self, [self], name=name)

    def tessellation_points(self, resolution=20):
        # plane = Plane3D.from_normal(self.center, self.normal)
        tessellation_points_3d = [
                                     self.center + self.major_axis * math.cos(
                                         teta) * self.major_dir
                                     + self.minor_axis * math.sin(
                                         teta) * self.major_dir.cross(
                                         self.normal) for teta in
                                     npy.linspace(0, volmdlr.TWO_PI,
                                                  resolution + 1)][:-1]
        return tessellation_points_3d

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D):
        # minor_dir = self.normal.cross(self.major_dir)
        # frame = volmdlr.Frame3D(self.center, self.major_dir,
        #                         minor_dir, self.normal)
        frame = volmdlr.Frame3D(self.center, self.major_dir,
                                self.normal.cross(self.major_dir), self.normal)

        # Positionnement des points dans leur frame
        p1_new, p2_new = frame.new_coordinates(
            point1), frame.new_coordinates(point2)

        # Angle pour le p1
        # u1, u2 = p1_new.x / self.major_axis, p1_new.y / self.minor_axis
        # theta1 = volmdlr.core.sin_cos_angle(u1, u2)
        theta1 = volmdlr.core.sin_cos_angle(p1_new.x / self.major_axis, p1_new.y / self.minor_axis)

        # Angle pour le p2
        # u3, u4 = p2_new.x / self.major_axis, p2_new.y / self.minor_axis
        # theta2 = volmdlr.core.sin_cos_angle(u3, u4)
        theta2 = volmdlr.core.sin_cos_angle(p2_new.x / self.major_axis, p2_new.y / self.minor_axis)

        if theta1 > theta2:  # sens trigo
            angle = math.pi + (theta1 + theta2) / 2
        else:
            angle = (theta1 + theta2) / 2

        # p_3 = volmdlr.Point3D(self.major_axis * math.cos(angle),
        #                       self.minor_axis * math.sin(angle), 0)
        # p3 = frame.old_coordinates(p_3)
        p3 = frame.old_coordinates(volmdlr.Point3D(self.major_axis * math.cos(angle),
                                                   self.minor_axis * math.sin(angle), 0))

        return volmdlr.edges.ArcEllipse3D(point1, p3, point2, self.center,
                                          self.major_dir)

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        xc, yc, zc = npy.round(1000 * self.center.vector, ndigits)
        major_vector = self.center + self.major_axis / 2 * self.major_dir
        xmaj, ymaj, zmaj = npy.round(1000 * major_vector.vector, ndigits)
        minor_vector = self.center + self.minor_axis / 2 * self.normal.cross(
            self.major_dir)
        xmin, ymin, zmin = npy.round(1000 * minor_vector.vector, ndigits)
        return '{} = Part.Ellipse(fc.Vector({},{},{}), fc.Vector({},{},{}), fc.Vector({},{},{}))\n'.format(
            name, xmaj, ymaj, zmaj, xmin, ymin, zmin, xc, yc, zc)

    def rotation(self, rot_center, axis, angle, copy=True):
        new_center = self.center.rotation(rot_center, axis, angle, True)
        new_normal = self.normal.rotation(rot_center, axis, angle, True)
        new_major_dir = self.major_dir.rotation(rot_center, axis, angle, True)
        if copy:
            return Ellipse3D(self.major_axis, self.minor_axis, new_center,
                             new_normal, new_major_dir, self.name)
        else:
            self.center = new_center
            self.normal = new_normal
            self.major_dir = new_major_dir

    def translation(self, offset, copy=True):
        new_center = self.center.translation(offset, True)
        new_normal = self.normal.translation(offset, True)
        new_major_dir = self.major_dir.translation(offset, True)
        if copy:
            return Ellipse3D(self.major_axis, self.minor_axis, new_center,
                             new_normal, new_major_dir, self.name)
        else:
            self.center = new_center
            self.normal = new_normal
            self.major_dir = new_major_dir

    def plot(self, ax=None, color='k', alpha=1, edge_details=False):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        x = []
        y = []
        z = []
        for px, py, pz in self.tessellation_points():
            x.append(px)
            y.append(py)
            z.append(pz)
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])
        ax.plot(x, y, z, color)
        return ax

    @classmethod
    def from_step(cls, arguments, object_dict):
        center = object_dict[arguments[1]].origin
        normal = object_dict[arguments[1]].u  # ancien w
        major_dir = object_dict[arguments[1]].v  # ancien u
        major_axis = float(arguments[2]) / 1000
        minor_axis = float(arguments[3]) / 1000
        return cls(major_axis, minor_axis, center, normal, major_dir,
                   arguments[0][1:-1])


class ClosedPolygon3D(Contour3D, ClosedPolygon):

    def __init__(self, points: List[volmdlr.Point3D], name: str = ''):
        self.points = points
        self.line_segments = self._line_segments()

        Contour3D.__init__(self, self.line_segments, name)

    def _line_segments(self):
        lines = []
        if len(self.points) > 1:
            for p1, p2 in zip(self.points,
                              list(self.points[1:]) + [self.points[0]]):
                lines.append(volmdlr.edges.LineSegment3D(p1, p2))
        return lines

    def copy(self, *args, **kwargs):
        points = [p.copy() for p in self.points]
        return ClosedPolygon2D(points, self.name)

    def __hash__(self):
        return sum([hash(p) for p in self.points])

    def __eq__(self, other_):
        if not isinstance(other_, self.__class__):
            return False
        equal = True
        for point, other_point in zip(self.points, other_.points):
            equal = (equal and point == other_point)
        return equal

    def plot(self, ax=None, color='k', alpha=1, edge_details=False):
        for line_segment in self.line_segments:
            ax = line_segment.plot(ax=ax, color=color, alpha=alpha,
                                   edge_ends=True, edge_direction=True)
        return ax

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            return ClosedPolygon3D(
                [p.rotation(center, axis, angle, copy=True) for p in
                 self.points])
        else:
            for p in self.points:
                p.rotation(center, axis, angle, copy=False)

    def translation(self, offset, copy=True):
        if copy:
            new_points = [point.translation(offset, copy=True) for point in
                          self.points]
            return ClosedPolygon3D(new_points, self.name)
        else:
            for point in self.points:
                point.translation(offset, copy=False)

    def to_2d(self, plane_origin, x, y):
        points2d = [point.to_2d(plane_origin, x, y) for point in self.points]
        return ClosedPolygon2D(points2d)

    def sewing_with(self, other_poly3d, x, y, normal, resolution=20):
        self_center, other_center = self.average_center_point(), \
                                    other_poly3d.average_center_point()

        self_poly2d, other_poly2d = self.to_2d(self_center, x, y), \
            other_poly3d.to_2d(other_center, x, y)
        self_center2d, other_center2d = self_poly2d.center_of_mass(), \
            other_poly2d.center_of_mass()
        self_poly2d.translation(-self_center2d, copy=False)
        other_poly2d.translation(-other_center2d, copy=False)

        bbox_self2d, bbox_other2d = self_poly2d.bounding_rectangle(), \
            other_poly2d.bounding_rectangle()
        position = [abs(value) for value in bbox_self2d] \
            + [abs(value) for value in bbox_other2d]
        max_scale = 2 * max(position)

        lines = [volmdlr.edges.LineSegment2D(volmdlr.O2D, max_scale * (
                volmdlr.X2D * math.sin(n * 2 * math.pi / resolution) +
                volmdlr.Y2D * math.cos(n * 2 * math.pi / resolution))
                                             ) for n in range(resolution)]

        self_new_points, other_new_points = [], []
        for l in lines:
            for self_line in self_poly2d.line_segments:
                intersect = l.linesegment_intersections(self_line)
                if intersect:
                    self_new_points.extend(intersect)
                    break

            for other_line in other_poly2d.line_segments:
                intersect = l.linesegment_intersections(other_line)
                if intersect:
                    other_new_points.extend(intersect)
                    break

        new_self_poly2d, new_other_poly2d = ClosedPolygon2D(
            self_new_points), ClosedPolygon2D(other_new_points)
        new_self_poly2d.translation(self_center2d, copy=False)
        new_other_poly2d.translation(other_center2d, copy=False)

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
        return ClosedPolygon3D(self.simplify_polygon(
            min_distance=min_distance, max_distance=max_distance).points)

    def sewing(self, polygon2, x, y):
        """
        x and y are used for plane projection to make sure it is being projected in the right plane
        """
        center1, center2 = self.average_center_point(), polygon2.average_center_point()
        center1_, center2_ = volmdlr.Point3D(center1.x, center1.y, 0), volmdlr.Point3D(center2.x, center2.y, 0)
        new_polygon1, new_polygon2 = self.translation(-center1_), polygon2.translation(-center2_)
        new_center1, new_center2 = new_polygon1.average_center_point(), new_polygon2.average_center_point()

        new_polygon1_2d, new_polygon2_2d = new_polygon1.to_2d(new_center1, x, y), new_polygon2.to_2d(new_center2, x, y)
        # barycenter1_2d = new_polygon1_2d.barycenter()
        # barycenter2_2d = new_polygon2_2d.barycenter()

        # # # ax2d = new_polygon1_2d.plot(color='r')
        # # # new_polygon2_2d.plot(ax=ax2d, color='g')
        # # barycenter1_2d = new_polygon1_2d.points[0]
        # # for point in new_polygon1_2d.points[1:]:
        # #     barycenter1_2d += point
        # # barycenter1_2d = barycenter1_2d / len(new_polygon1_2d.points)
        # # barycenter1_2d.plot(ax=ax2d, color='y')
        #
        # barycenter2_2d = new_polygon2_2d.points[0]
        # for point in new_polygon2_2d.points[1:]:
        #     barycenter2_2d += point
        # barycenter2_2d = barycenter2_2d / len(new_polygon2_2d.points)
        #
        # # barycenter2_2d.plot(ax=ax2d, color='r')
        #
        # # ax3d= new_polygon1.plot(color= 'r')
        # # new_polygon2.plot(ax=ax3d, color= 'g')
        # # volmdlr.Point3D(0,0, center1.z).plot(ax=ax3d)
        # # volmdlr.Point3D(0,0, center2.z).plot(ax=ax3d, color = 'r')

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
                        dict_closing_pairs[self.points[i - 1]] = (
                            previous_closing_point_index,
                            closing_point_index)

                if point_polygon1 == new_polygon1.points[0]:
                    if list(dict_closing_pairs.values())[-1][-1] != \
                            list(dict_closing_pairs.values())[0][0]:
                        dict_closing_pairs[self.points[0]] = (
                            list(dict_closing_pairs.values())[-1][-1],
                            list(dict_closing_pairs.values())[0][0])

                real_closing_point = polygon2.points[closing_point_index]

                if not volmdlr.edges.LineSegment3D(
                        self.points[new_polygon1.points.index(point_polygon1)],
                        real_closing_point).point_belongs(self.points[i - 1]):
                    face = volmdlr.faces.Triangle3D(self.points[new_polygon1.points.index(
                        point_polygon1)], self.points[i - 1],
                                      real_closing_point,
                                                    alpha=0.9,
                                                    color=(1, 0.1, 0.1))
                    triangles.append(face)

                list_closing_point_indexes.append(closing_point_index)
                previous_closing_point_index = closing_point_index
        triangles += polygon2.close_sewing(dict_closing_pairs)
        # for i, point_polygon2 in enumerate(
        #         new_polygon2.points + [new_polygon2.points[0]]):
        #     for j, index in enumerate(list(dict_closing_pairs.values())):
        #         if i != 0:
        #             if i - 1 >= index[0] and i <= index[1]:
        #                 # a = polygon2.points[i - 1]
        #                 # b = polygon2.points[
        #                 #     new_polygon2.points.index(point_polygon2)]
        #                 # c = list(dict_closing_pairs.keys())[j]
        #                 # if a != b and a != c and b != c and not volmdlr.edges.LineSegment3D(
        #                 #         a, c).point_belongs(b):
        #                 triangles.append([polygon2.points[i - 1],
        #                                   polygon2.points[
        #                                       new_polygon2.points.index(
        #                                           point_polygon2)],
        #                                   list(dict_closing_pairs.keys())[j]])
        #             elif index[0] > index[1]:
        #                 if (i - 1 <= index[0] and i <= index[1]) or (
        #                         (i - 1 >= index[0]) and i >= index[1]):
        #                     # a = polygon2.points[i - 1]
        #                     # b = polygon2.points[
        #                     #     new_polygon2.points.index(point_polygon2)]
        #                     # c = list(dict_closing_pairs.keys())[j]
        #                     # if a != b and a != c and b != c and not volmdlr.edges.LineSegment3D(
        #                     #         a, c).point_belongs(b):
        #                     triangles.append([polygon2.points[i - 1],
        #                                       polygon2.points[
        #                                           new_polygon2.points.index(
        #                                               point_polygon2)],
        #                                       list(dict_closing_pairs.keys())[j]])

        # print("list closing points indexes :", list_closing_point_indexes)
        # faces = []
        # coords = triangles
        # for trio in coords:
        #     faces.append(volmdlr.faces.Triangle3D(trio[0], trio[1], trio[2]))
        # volum = volmdlr.core.VolumeModel(faces)
        # volum.babylonjs()
        # ax = self.plot()
        # polygon2.plot(ax=ax, color='y')
        # for face in faces:
        #     volmdlr.edges.LineSegment3D(face.point1, face.point2).plot(ax=ax,
        #                                                                color='b')
        #     volmdlr.edges.LineSegment3D(face.point1, face.point3).plot(ax=ax,
        #                                                                color='b')
        #     volmdlr.edges.LineSegment3D(face.point2, face.point3).plot(ax=ax,
        #                                                                color='b')

        return triangles

    def close_sewing(self, dict_closing_pairs):
        triangles = []
        for i, point_polygon2 in enumerate(
                self.points + [self.points[0]]):
            for j, index in enumerate(list(dict_closing_pairs.values())):
                if i != 0:
                    if i - 1 >= index[0] and i <= index[1]:
                        face = volmdlr.faces.Triangle3D(self.points[i - 1],
                                                        point_polygon2,
                                                        list(
                                                            dict_closing_pairs.keys())[
                                                            j],
                                                        alpha=0.9,
                                                        color=(1, 0.1, 0.1))
                        triangles.append(face)
                    elif index[0] > index[1]:
                        if (i - 1 <= index[0] and i <= index[1]) or (
                                (i - 1 >= index[0]) and i >= index[1]):
                            face = volmdlr.faces.Triangle3D(
                                self.points[i - 1],
                                point_polygon2,
                                list(dict_closing_pairs.keys())[j],
                                alpha=0.9,
                                color=(1, 0.1, 0.1))
                            triangles.append(face)
        return triangles

    def check_sewing(self, polygon2, sewing_faces):
        if not len(self.line_segments) + len(polygon2.line_segments) == len(sewing_faces):
            return False
        return True

    @staticmethod
    def validate_closing_point(closing_point_index, list_closing_point_indexes,
                               passed_by_zero_index, ratio_denom):
        if closing_point_index == list_closing_point_indexes[-1]:
            return closing_point_index, [], passed_by_zero_index
        list_remove_closing_points = []
        ratio = (list_closing_point_indexes[-1] - closing_point_index) / ratio_denom
        new_list_closing_point = list(dict.fromkeys(list_closing_point_indexes))
        if closing_point_index < list_closing_point_indexes[-1]:
            if len(list_closing_point_indexes) > 2 and \
                    list_closing_point_indexes[
                        -2] < closing_point_index < list_closing_point_indexes[-1] - 1:
                list_remove_closing_points.append(list_closing_point_indexes[-1])
            elif len(new_list_closing_point) > 3 and \
                    new_list_closing_point[
                        -3] < closing_point_index < new_list_closing_point[-1] - 3:
                for idx in list_closing_point_indexes:
                    if idx in (new_list_closing_point[-1], new_list_closing_point[-2]):
                        list_remove_closing_points.append(idx)
            elif len(list_closing_point_indexes) > 5 and \
                    list_closing_point_indexes[
                        -5] < closing_point_index < \
                    list_closing_point_indexes[-1] - 5 and \
                    list_closing_point_indexes[-5] != \
                    list_closing_point_indexes[-6]:
                list_remove_closing_points.extend(
                    [list_closing_point_indexes[-1],
                     list_closing_point_indexes[-2],
                     list_closing_point_indexes[-3],
                     list_closing_point_indexes[-4]])
            elif closing_point_index in list_closing_point_indexes:
                closing_point_index = list_closing_point_indexes[-1]
            elif math.isclose(ratio, 0, abs_tol=0.3):
                closing_point_index = list_closing_point_indexes[-1]
            else:
                if passed_by_zero_index:
                    ratio = (list_closing_point_indexes[
                                 0] - closing_point_index) / ratio_denom
                    if math.isclose(ratio, 1, abs_tol=0.3):
                        closing_point_index = list_closing_point_indexes[0]
                else:
                    if closing_point_index > list_closing_point_indexes[0]:
                        ratio1 = (closing_point_index -
                                  list_closing_point_indexes[
                                      0]) / ratio_denom
                        if math.isclose(ratio1, 0, abs_tol=0.3) and\
                                math.isclose(ratio, 1, abs_tol=0.3):
                            passed_by_zero_index = True
                            closing_point_index = list_closing_point_indexes[0]
                        else:
                            closing_point_index = list_closing_point_indexes[-1]
                    else:
                        ratio1 = (closing_point_index -
                                  list_closing_point_indexes[
                                      -1]) / ratio_denom
                        if math.isclose(ratio1, -1, abs_tol=0.45):
                            if ratio_denom - list_closing_point_indexes[-1] < 10:
                                passed_by_zero_index = True
                            elif ratio_denom - list_closing_point_indexes[
                                    -1] >= 5:
                                closing_point_index = \
                                    list_closing_point_indexes[-1] + 5
                            else:
                                closing_point_index = ratio_denom - 2
                        else:
                            if list_closing_point_indexes.count(
                                    list_closing_point_indexes[-1]) < 5:
                                closing_point_index = \
                                    list_closing_point_indexes[-1]
                            elif ratio_denom - list_closing_point_indexes[-1] >= 5:
                                closing_point_index =\
                                    list_closing_point_indexes[-1] + 5
                            elif list_closing_point_indexes[-1] <= ratio_denom - 2:
                                closing_point_index = ratio_denom - 2

        elif closing_point_index in list_closing_point_indexes:
            closing_point_index = list_closing_point_indexes[-1]
        elif len(list_closing_point_indexes) > 2 and \
                list_closing_point_indexes[0] < closing_point_index < \
                list_closing_point_indexes[-1]:
            closing_point_index = list_closing_point_indexes[-1]
        elif passed_by_zero_index and closing_point_index > \
                list_closing_point_indexes[0]:
            closing_point_index = list_closing_point_indexes[-1]
        elif list_closing_point_indexes[0] == 0 and math.isclose(ratio, -1,
                                                                 abs_tol=0.3):
            closing_point_index = list_closing_point_indexes[-1]
        elif math.isclose(ratio, -1, abs_tol=0.3):
            closing_point_index = list_closing_point_indexes[-1]
        elif closing_point_index - list_closing_point_indexes[-1] > 5 and \
                list_closing_point_indexes[-1] + 4 <= ratio_denom - 1:
            closing_point_index = list_closing_point_indexes[-1] + 4

        return closing_point_index, list_remove_closing_points, passed_by_zero_index

    def sewing3(self, polygon2, x, y):
        polygon1_2d = self.to_2d(volmdlr.O2D, x, y)
        polygon2_2d = polygon2.to_2d(volmdlr.O2D, x, y)
        polygon1_3d = self
        polygon2_3d = polygon2
        if polygon2_2d.area() < polygon1_2d.area():
            polygon1_2d, polygon2_2d = polygon2_2d, polygon1_2d
            polygon1_3d = polygon2
            polygon2_3d = self

        valid_polygon = polygon1_2d.validate_sewing_polygons(polygon2_2d)
        if not valid_polygon:
            polygon1_2d = polygon1_2d.get_valid_sewing_polygons(polygon2_2d)
            polygon1_3d = polygon1_2d.to_3d(volmdlr.Point3D(0, 0,
                                                            polygon1_3d.points[0].z), x, y)
        # ax=polygon1_2d.plot()
        # polygon2_2d.plot(ax=ax, color='r')
        dict_closing_pairs = {}
        triangles = []
        list_closing_point_indexes = []
        passed_by_zero_index = False
        ratio_denom = len(polygon2_2d.points)
        previous_closing_point_index = None
        for i, primitive1 in enumerate(polygon1_2d.line_segments):
            list_remove_closing_points = []

            closing_point = polygon1_2d.get_closing_point(polygon2_2d,
                                                          primitive1)
            if closing_point == volmdlr.O2D:
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
                closing_point_index, list_remove_closing_points,\
                    passed_by_zero_index = self.validate_closing_point(
                        closing_point_index, list_closing_point_indexes,
                        passed_by_zero_index, ratio_denom)

            if list_remove_closing_points:
                new_list_closing_point_indexes = list(dict.fromkeys(list_closing_point_indexes))
                new_list_remove_closing_indexes = list(dict.fromkeys(list_remove_closing_points))
                for n in range(len(list_remove_closing_points)):
                    new_face = volmdlr.faces.Triangle3D(
                        triangles[-(n + 1)].point1,
                        triangles[-(n + 1)].point2,
                        polygon2_3d.points[closing_point_index],
                        alpha=0.9,
                        color=(1, 0.1, 0.1))
                    triangles[-(n + 1)] = new_face
                number_common_closing_point = 0
                for fc in triangles:
                    if fc.point3 == polygon2_3d.points[closing_point_index]:
                        number_common_closing_point += 1
                for n in range(0,
                               len(new_list_remove_closing_indexes) // 2 + 1):
                    del dict_closing_pairs[
                        list(dict_closing_pairs.keys())[-1]]

                dict_closing_pairs[triangles[-number_common_closing_point].point1] = (
                    new_list_closing_point_indexes[
                        -(len(new_list_remove_closing_indexes) + 1)],
                    closing_point_index)
                for pt_index in list_remove_closing_points:
                    list_closing_point_indexes.remove(pt_index)
                list_closing_point_indexes.append(closing_point_index)
            elif closing_point_index != previous_closing_point_index:
                dict_closing_pairs[polygon1_3d.line_segments[i].start] = (
                    previous_closing_point_index, closing_point_index)

            face = volmdlr.faces.Triangle3D(polygon1_3d.line_segments[i].start,
                                            polygon1_3d.line_segments[i].end,
                                            polygon2_3d.points[closing_point_index], alpha=0.9, color=(1, 0.1, 0.1))
            triangles.append(face)
            list_closing_point_indexes.append(closing_point_index)
            previous_closing_point_index = closing_point_index
            if primitive1 == polygon1_2d.line_segments[-1]:
                if list_closing_point_indexes[-1] != \
                        list_closing_point_indexes[0]:
                    ratio = (list_closing_point_indexes[-1] -
                             list_closing_point_indexes[0]) / len(
                        polygon2_2d.points)
                    if math.isclose(ratio, -1,
                                    abs_tol=0.2) and passed_by_zero_index:
                        dict_closing_pairs[
                            polygon1_3d.points[0]] = (
                            list_closing_point_indexes[-2],
                            list_closing_point_indexes[0])
                        new_face = volmdlr.faces.Triangle3D(
                            triangles[-1].point1,
                            triangles[-1].point2,
                            polygon2_3d.points[list_closing_point_indexes[-2]],
                            alpha=0.9,
                            color=(1, 0.1, 0.1))
                        triangles.remove(triangles[-1])
                        triangles.append(new_face)
                    else:
                        dict_closing_pairs[polygon1_3d.points[0]] = (
                            list(dict_closing_pairs.values())[-1][-1],
                            list(dict_closing_pairs.values())[0][0])
        triangles += polygon2_3d.close_sewing(dict_closing_pairs)

        # print('list closing indexes :', list_closing_point_indexes)
        # print('length polygon2 points: ', len(polygon2_3d.points))
        # print('dict_closing_pairs :', dict_closing_pairs)

        # volum = volmdlr.core.VolumeModel(triangles)
        # volum.babylonjs()
        # print('p1 3d points :', self.points)
        # print('p2 3d points :', polygon2.points)
        return triangles
