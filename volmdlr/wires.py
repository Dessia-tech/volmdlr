#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script checking offset and Curvilinear absissa of roundedline2D
"""

import math
import numpy as npy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import volmdlr
import volmdlr.core
from volmdlr.core_compiled import polygon_point_belongs


class Wire2D(volmdlr.core.CompositePrimitive2D):
    """
    A collection of simple primitives, following each other making a wire
    """

    def __init__(self, primitives, name=''):
        volmdlr.core.CompositePrimitive2D.__init__(self, primitives, name)

    # TODO: method to check if it is a wire

    def length(self):
        length = 0.
        for primitive in self.primitives:
            length += primitive.length()
        return length

    def point_at_abscissa(self, curvilinear_abscissa: float):
        length = 0.
        for primitive in self.primitives:
            primitive_length = primitive.length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.point_at_abscissa(
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

    def line_intersections(self, line: 'Line2D'):
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


class Wire3D(volmdlr.core.CompositePrimitive3D):
    """
    A collection of simple primitives, following each other making a wire
    """

    def __init__(self, primitives, name=''):
        volmdlr.core.CompositePrimitive3D.__init__(self, primitives, name)

    def length(self):
        length = 0.
        for primitive in self.primitives:
            length += primitive.length()
        return length

    def point_at_abscissa(self, curvilinear_abscissa):
        length = 0.
        for primitive in self.primitives:
            primitive_length = primitive.length()
            if length + primitive_length >= curvilinear_abscissa:
                return primitive.point_at_abscissa(
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

class Contour2D(Wire2D):
    """
    A collection of 2D primitives forming a closed wire2D
    TODO : center_of_mass and second_moment_area should be changed accordingly to
    area considering the triangle drawn by the arcs
    """
    _non_serializable_attributes = ['internal_arcs', 'external_arcs',
                                    'polygon', 'straight_line_contour_polygon']

    def __init__(self, primitives, name=''):
        Wire2D.__init__(self, primitives, name)
        self._utd_analysis = False
        # self.tessel_points = self.clean_points()

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

    def to_3d(self, plane_origin, x, y):
        p3d = []
        for edge in self.primitives:
            p3d.append(edge.to_3d(plane_origin, x, y))

        return Contour3D(p3d)

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
        return (volmdlr.Point2D((xmin, ymin)), volmdlr.Point2D((xmax, ymax)))


    # def To3D(self, plane_origin, x, y, name=None):
    #     if name is None:
    #         name = '3D of {}'.format(self.name)
    #     primitives3D = [p.To3D(plane_origin, x, y) for p in self.primitives]
    #     return Contour3D(primitives=primitives3D, name=name)

    def area(self):
        if len(self.primitives) == 1:
            return self.primitives[0].area()

        A = self.polygon.area()

        for arc in self.internal_arcs:
            triangle = Polygon2D([arc.start, arc.center, arc.end])
            A = A - arc.area() + triangle.area()
        for arc in self.external_arcs:
            triangle = Polygon2D([arc.start, arc.center, arc.end])
            A = A + arc.area() - triangle.area()

        return A

    def center_of_mass(self):
        if len(self.primitives) == 1:
            return self.primitives[0].center_of_mass()

        area = self.polygon.area()
        if area > 0.:
            c = area * self.polygon.center_of_mass()
        else:
            c = volmdlr.O2D

        for arc in self.internal_arcs:
            arc_area = arc.area()
            c -= arc_area * arc.center_of_mass()
            area -= arc_area
        for arc in self.external_arcs:
            arc_area = arc.area()
            c += arc_area * arc.center_of_mass()
            area += arc_area
        if area != 0:
            return c / area
        else:
            return False

    def second_moment_area(self, point):
        if len(self.primitives) == 1:
            return self.primitives[0].second_moment_area(point)

        A = self.polygon.second_moment_area(point)
        for arc in self.internal_arcs:
            A -= arc.second_moment_area(point)
        for arc in self.external_arcs:
            A += arc.second_moment_area(point)

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
        return volmdlr.Point2D((x, y))

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
                d1, d2 = (points_to_add[0] - points[0]).norm(), (
                            points_to_add[0] - points[-1]).norm()
                d3, d4 = (points_to_add[-1] - points[0]).norm(), (
                            points_to_add[-1] - points[-1]).norm()
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
            p = volmdlr.Point2D.random(xmin, xmax, ymin, ymax)
            if self.point_belongs(p):
                return p
    # def line_intersections(self, line:Line2D) -> List[Tuple[volmdlr.Point2D, Primitive2D]]:
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
            if isinstance(intersections[0][0], volmdlr.Point2D) and \
                    isinstance(intersections[1][0], volmdlr.Point2D):
                ip1, ip2 = sorted([self.primitives.index(intersections[0][1]),
                                   self.primitives.index(intersections[1][1])])

                sp11, sp12 = intersections[0][1].split(intersections[0][0])
                sp21, sp22 = intersections[1][1].split(intersections[1][0])

                primitives1 = self.primitives[:ip1]
                primitives1.append(sp11)
                primitives1.append(volmdlr.edges.LineSegment2D(intersections[0][0],
                                                 intersections[1][0]))
                primitives1.append(sp22)
                primitives1.extend(self.primitives[ip2 + 1:])

                primitives2 = self.primitives[ip1 + 1:ip2]
                primitives2.append(sp21)
                primitives2.append(volmdlr.edges.LineSegment2D(intersections[1][0],
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
            cut_line = volmdlr.edges.Line2D(volmdlr.Point2D(xi, 0),
                                            volmdlr.Point2D(xi, 1))

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
                p = volmdlr.Point2D((xi, yi))
                if self.point_belongs(p):
                    point_index[p] = ip
                    points.append(p)
                    ip += 1

        for i in range(n):
            for j in range(m):
                p1 = volmdlr.Point2D((x[i], y[j]))
                p2 = volmdlr.Point2D((x[i + 1], y[j]))
                p3 = volmdlr.Point2D((x[i + 1], y[j + 1]))
                p4 = volmdlr.Point2D((x[i], y[j + 1]))
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

        return volmdlr.display_mesh.DisplayMesh2D(points, triangles)

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

    def area(self):

        x = [point.x for point in self.points]
        y = [point.y for point in self.points]

        return 0.5 * npy.abs(
            npy.dot(x, npy.roll(y, 1)) - npy.dot(y, npy.roll(x, 1)))

    def center_of_mass(self):

        x = [point.vector[0] for point in self.points]
        y = [point.vector[1] for point in self.points]

        xi_xi1 = x + npy.roll(x, -1)
        yi_yi1 = y + npy.roll(y, -1)
        xi_yi1 = npy.multiply(x, npy.roll(y, -1))
        xi1_yi = npy.multiply(npy.roll(x, -1), y)

        a = 0.5 * npy.sum(xi_yi1 - xi1_yi)  # signed area!
        #        a=self.area()
        if not math.isclose(a, 0, abs_tol=1e-08):
            cx = npy.sum(npy.multiply(xi_xi1, (xi_yi1 - xi1_yi))) / 6. / a
            cy = npy.sum(npy.multiply(yi_yi1, (xi_yi1 - xi1_yi))) / 6. / a
            return volmdlr.Point2D((cx, cy))

        else:
            raise NotImplementedError

    def PointBelongs(self, point):
        """
        Ray casting algorithm copied from internet...
        """
        return polygon_point_belongs((point.x, point.y),
                                     [(p.x, p.y) for p in self.points])

    def second_moment_area(self, point):
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
            lines.append(volmdlr.edges.LineSegment2D(p1, p2))
        return lines

    def rotation(self, center, angle, copy=True):
        if copy:
            return Polygon2D(
                [p.rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.rotation(center, angle, copy=False)

    def translation(self, offset, copy=True):
        if copy:
            return Polygon2D(
                [p.translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.translation(offset, copy=False)

    def offset(self, offset):
        bound = self.bounding_rectangle2()
        max_offset_len = bound.min_length() / 2
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

            normal_vector1 = - vectors[2 * i - 1].normal_vector()
            normal_vector2 = vectors[2 * i].normal_vector()
            normal_vector1.normalize()
            normal_vector2.normalize()
            alpha = math.acos(normal_vector1.dot(normal_vector2))

            offset_point = self.points[i] + offset / math.cos(alpha / 2) * \
                offset_vectors[i]
            offset_points.append(offset_point)

        return self.__class__(offset_points)

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

                        line1 = volmdlr.edges.LineSegment2D(
                            volmdlr.Point2D(self.points[segment1[0]]),
                            volmdlr.Point2D(self.points[segment1[1]]))
                        line2 = volmdlr.edges.LineSegment2D(
                            volmdlr.Point2D(self.points[segment2[0]]),
                            volmdlr.Point2D(self.points[segment2[1]]))

                        p, a, b = volmdlr.Point2D.LinesIntersection(line1, line2, True)

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
        ymax, pos_ymax = volmdlr.max_pos([pt.vector[1] for pt in points])
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

    def plot(self, ax=None, color='k',
                plot_points=False, point_numbering=False,
                fill=False, fill_color='w'):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')

        if fill:
            ax.fill([p[0] for p in self.points], [p[1] for p in self.points],
                    facecolor=fill_color)
        for ls in self.line_segments:
            ls.plot(ax=ax ,color=color)

        if plot_points or point_numbering:
            for point in self.points:
                point.plot(ax=ax, color=color)

        if point_numbering:
            for ip, point in enumerate(self.points):
                ax.text(*point, 'point {}'.format(ip+1),
                        ha='center', va='top')

        ax.margins(0.1)
        plt.show()

        return ax


class Circle2D(Contour2D):
    _non_serializable_attributes = ['internal_arcs', 'external_arcs',
                                    'polygon', 'straight_line_contour_polygon',
                                    'primitives', 'basis_primitives']

    def __init__(self, center: volmdlr.Point2D, radius: float, name: str = ''):
        self.center = center
        self.radius = radius
        self.angle = volmdlr.TWO_PI
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
        return [(self.center
                 + self.radius * math.cos(teta) * volmdlr.X2D
                 + self.radius * math.sin(teta) * volmdlr.Y2D)\
                for teta in npy.linspace(0, volmdlr.TWO_PI, resolution + 1)][:-1]

    def point_belongs(self, point, tolerance=1e-9):
        return point.point_distance(self.center) <= self.radius + tolerance

    def line_intersections(self, line):
        Q = volmdlr.Vector2D(self.center.vector)
        if line.points[0].vector == self.center.vector:
            P1 = volmdlr.Vector2D(line.points[1].vector)
            V = volmdlr.Vector2D((line.points[0] - line.points[1]).vector)
        else:
            P1 = volmdlr.Vector2D(line.points[0].vector)
            V = volmdlr.Vector2D((line.points[1] - line.points[0]).vector)
        a = V.Dot(V)
        b = 2 * V.Dot(P1 - Q)
        c = P1.Dot(P1) + Q.Dot(Q) - 2 * P1.Dot(Q) - self.radius ** 2
        disc = b ** 2 - 4 * a * c
        if disc < 0:
            return []

        sqrt_disc = math.sqrt(disc)
        t1 = (-b + sqrt_disc) / (2 * a)
        t2 = (-b - sqrt_disc) / (2 * a)
        if line.__class__ is volmdlr.edges.Line2D:
            if t1 == t2:
                return [volmdlr.Point2D((P1 + t1 * V).vector)]
            else:
                return [volmdlr.Point2D((P1 + t1 * V).vector),
                        volmdlr.Point2D((P1 + t2 * V).vector)]
        else:
            if not 0 <= t1 <= 1 and not 0 <= t2 <= 1:
                return None
            elif 0 <= t1 <= 1 and not 0 <= t2 <= 1:
                return [volmdlr.Point2D((P1 + t1 * V).vector)]
            elif not 0 <= t1 <= 1 and 0 <= t2 <= 1:
                return [volmdlr.Point2D((P1 + t2 * V).vector)]
            else:
                return [volmdlr.Point2D((P1 + t1 * V).vector),
                        volmdlr.Point2D((P1 + t2 * V).vector)]

    def length(self):
        return volmdlr.TWO_PI * self.radius

    def plot(self, ax=None, linestyle='-', color='k', linewidth=1):
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
            self.utd_geo_points = False

    def translation(self, offset, copy=True):
        if copy:
            return Circle2D(self.center.translation(offset, copy=True),
                            self.radius)
        else:
            self.center.translation(offset, copy=False)
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

    def area(self):
        return math.pi * self.radius ** 2

    def second_moment_area(self, point):
        """
        Second moment area of part of disk
        """
        I = math.pi * self.radius ** 4 / 4
        Ic = npy.array([[I, 0], [0, I]])
        return geometry.Huygens2D(Ic, self.area(), self.center, point)

    def center_of_mass(self):
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

    def point_at_abscissa(self, curvilinear_abscissa):
        start = self.center + self.radius * volmdlr.X3D
        return start.rotation(self.center,
                              curvilinear_abscissa / self.radius)

    def triangulation(self, n=35):
        l = self.length()
        points = [self.point_at_abscissa(l * i / n) for i in range(n)]
        points.append(self.center)
        triangles = [(i, i + 1, n) for i in range(n - 1)] + [(n - 1, 0, n)]
        return volmdlr.display_mesh.DisplayMesh(points, triangles)

    def polygon_points(self, points_per_radian=10, min_x_density=None,
                       min_y_density=None):
        return volmdlr.edges.Arc2D.polygon_points(self, points_per_radian=points_per_radian,
                                    min_x_density=min_x_density,
                                    min_y_density=min_y_density)



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
                d1, d2 = (points_to_add[0] - points[0]).norm(), (
                            points_to_add[0] - points[-1]).norm()
                d3, d4 = (points_to_add[-1] - points[0]).norm(), (
                            points_to_add[-1] - points[-1]).norm()
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
            return Contour3D(new_edges, None, self.name)
        else:
            for edge in self.primitives:
                edge.translation(offset, copy=False)
            for point in self.tessel_points:
                point.translation(offset, copy=False)

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

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        for edge in self.primitives:
            edge.plot(ax=ax)

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
        points = [self.point_at_abscissa(i/n*l)\
                  for i in range(n)]
        return volmdlr.core.BoundingBox.from_points(points)

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

        tessellation_points_3D = [self.center
                                  + self.radius * math.cos(
            teta) * self.surface.frame.u
                                  + self.radius * math.sin(
            teta) * self.surface.frame.v \
                                  for teta in npy.linspace(0, TWO_PI,
                                                           resolution + 1)][
                                 :-1]
        return tessellation_points_3D

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
        new_frame = self.center.translation(offset, True)
        if copy:
            return Circle3D(new_frame, self.radius, self.frame,
                            self.name)
        else:
            self.frame = new_frame

    def plot(self, ax=None, color='k'):
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
            normal = object_dict[arguments[1]].v  ### ou w
            other_vec = None
        normal.normalize()

        return cls(center, radius, normal, arguments[0][1:-1], other_vec)

    def _bounding_box(self):
        """
        """
        u = self.normal.deterministic_unit_normal_vector()
        v = self.normal.cross(u)
        points = [self.frame.origin + self.radius * v \
                  for v in [self.frame.u,
                            -self.frame.u,
                            self.frame.v,
                            -self.frame.v]]
        return volmdlr.core.BoundingBox.from_points(points)

    def to_2d(self, plane_origin, x, y):
        z = x.cross(y)
        plane3d = volmdlr.faces.Plane3D(volmdlr.Frame3D(plane_origin, x, y, z))
        return Circle2D(plane3d.point3d_to_2d(self.center), self.radius)

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
            center, _ = l1.MinimumDistancePoints(l2)
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
            cylinder = volmdlr.faces.CylindricalSurface3D(volmdlr.Frame3D(self.center,
                                                    u,
                                                    v,
                                                    self.normal),
                                            self.radius
                                            )
            return cylinder.rectangular_cut(0, volmdlr.TWO_PI,
                                            0, extrusion_vector.norm())
        else:
            raise NotImplementedError('Elliptic faces not handled: dot={}'.format(
                self.normal.dot(extrusion_vector)
            ))

    def revolution(self, axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                   angle: float):
        line3d = volmdlr.edges.Line3D(axis_point, axis_point + axis)
        tore_center, _ = line3d.point_projection(self.center)
        u =  self.center - tore_center
        u.normalize()
        v = axis.cross(u)
        if not math.isclose(self.normal.dot(u), 0., abs_tol=1e-9):
            raise NotImplementedError(
                'Outside of plane revolution not supported')

        R = tore_center.point_distance(self.center)
        surface = volmdlr.faces.ToroidalSurface3D(volmdlr.Frame3D(tore_center, u, v, axis),
                                    R, self.radius)
        return surface.rectangular_cut(0, angle, 0, volmdlr.TWO_PI)


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

    def __init__(self, major_axis, minor_axis, center, normal, major_dir,
                 name=''):

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
        tessellation_points_3D = [self.center + self.major_axis * math.cos(
            teta) * self.major_dir + self.minor_axis * math.sin(
            teta) * self.major_dir.cross(self.normal) \
                                  for teta in npy.linspace(0, volmdlr.TWO_PI,
                                                           resolution + 1)][
                                 :-1]
        return tessellation_points_3D

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

    def plot(self, ax=None, color='k'):
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
