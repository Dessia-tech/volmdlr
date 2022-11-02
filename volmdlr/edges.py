#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from typing import List, Dict, Any
import math

import warnings
from packaging import version
import numpy as npy
import scipy as scp
import scipy.optimize

from geomdl import utilities, BSpline, fitting, operations
from geomdl.operations import length_curve, split_curve

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import __version__ as _mpl_version
import matplotlib.pyplot as plt
import matplotlib.patches

import plot_data.core as plot_data
import dessia_common as dc
import volmdlr.core_compiled
import volmdlr.core
import volmdlr.geometry


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


def insert_knots_and_mutiplicity(knots, knot_mutiplicities, knot_to_add, num):
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
    def __init__(self, start, end, name=''):
        self.start = start
        self.end = end
        self._length = None
        dc.DessiaObject.__init__(self, name=name)

    def __getitem__(self, key):
        if key == 0:
            return self.start
        if key == 1:
            return self.end
        raise IndexError

    def length(self):
        """
        Calculates the edge length
        :return: edges\' length
        """
        raise NotImplementedError(f'length method not implememented by {self.__class__.__name__}')

    def point_at_abscissa(self, abscissa):
        """
        Calcultes the point at given abscissa
        """
        raise NotImplementedError(f'point_at_absciss method not implememented by {self.__class__.__name__}')

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = None):
        """
        discretize a Edge to have "n" points
        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Usefull
            to mesh an arc
        :return: a list of sampled points
        """
        if number_points is None:
            number_points = 2
        step = self.length() / (number_points - 1)
        return [self.point_at_abscissa(i * step) for i in range(number_points)]

    def polygon_points(self, discretization_resolution: int):
        warnings.warn('polygon_points is deprecated,\
        please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(discretization_resolution)

    @classmethod
    def from_step(cls, arguments, object_dict):
        obj = object_dict[arguments[3]]
        p1 = object_dict[arguments[1]]
        p2 = object_dict[arguments[2]]
        if obj.__class__.__name__ == 'Line3D':
            return LineSegment3D(p1, p2, arguments[0][1:-1])

        else:
            if hasattr(obj, 'trim'):
                if obj.__class__.__name__ == 'Circle3D':
                    p1, p2 = p2, p1
                return obj.trim(p1, p2)

            else:
                raise NotImplementedError(f'Unsupported: {object_dict[arguments[3]]}')

    def normal_vector(self, abscissa):
        """
        Calculates the normal vector the edge at given abscissa
        :return: the normal vector
        """
        raise NotImplementedError('the normal_vector method must be'
                                  'overloaded by subclassing class')

    def unit_normal_vector(self, abscissa):
        """
        Calculates the unit normal vector the edge at given abscissa
        :param abscissa: edge abscissa
        :return: unit normal vector
        """
        raise NotImplementedError('the unit_normal_vector method must be'
                                  'overloaded by subclassing class')

    def direction_vector(self, abscissa):
        """
        Calculates the direction vector the edge at given abscissa
        :param abscissa: edge abscissa
        :return: direction vector
        """
        raise NotImplementedError('the direction_vector method must be'
                                  'overloaded by subclassing class')

    def unit_direction_vector(self, abscissa):
        """
        Calculates the unit direction vector the edge at given abscissa
        :param abscissa: edge abscissa
        :return: unit direction vector
        """
        raise NotImplementedError('the unit_direction_vector method must be'
                                  'overloaded by subclassing class')


class Line(dc.DessiaObject):
    """
    Abstract class
    """

    def __init__(self, point1, point2, name=''):
        self.point1 = point1
        self.point2 = point2
        dc.DessiaObject.__init__(self, name=name)

    def __getitem__(self, key):
        if key == 0:
            return self.point1
        elif key == 1:
            return self.point2
        else:
            raise IndexError

    def unit_direction_vector(self, *args, **kwargs):
        u = self.direction_vector()
        u.normalize()
        return u

    def direction_vector(self, *args, **kwargs):
        return self.point2 - self.point1

    def normal_vector(self, *args, **kwargs):
        return self.direction_vector().normal_vector()

    def unit_normal_vector(self, abscissa=0.):
        return self.unit_direction_vector().normal_vector()

    def point_projection(self, point):

        u = self.point2 - self.point1
        norm_u = u.norm()
        t = (point - self.point1).dot(u) / norm_u ** 2
        projection = self.point1 + t * u
        projection = projection.to_point()
        return projection, t * norm_u

    def abscissa(self, point):
        u = self.point2 - self.point1
        norm_u = u.norm()
        t = (point - self.point1).dot(u) / norm_u
        return t

    def split(self, split_point):
        return [self.__class__(self.point1, split_point),
                self.__class__(split_point, self.point2)]

    def is_between_points(self, point1: volmdlr.Point2D,
                          point2: volmdlr.Point2D):
        """
        Verifies if a line is between two points
        :param point1: first point
        :type point1: volmdlr.Point2D
        :param point2: second point
        :type point2: volmdlr.Point2D
        returns True is line is between the two given points or False if not
        """

        if point1 == point2:
            return False

        line_segment = LineSegment2D(point1, point2)
        if line_segment.line_intersections(self):
            return True
        return False


class LineSegment(Edge):
    """
    Abstract class
    """

    def length(self):
        if not self._length:
            self._length = self.end.point_distance(self.start)
        return self._length

    def abscissa(self, point, tol=1e-6):
        if point.point_distance(self.start) < tol:
            return 0
        if point.point_distance(self.end) < tol:
            return self.length()

        u = self.end - self.start
        length = u.norm()
        t = (point - self.start).dot(u) / length
        if t < -1e-9 or t > length + 1e-9:
            raise ValueError(f'Point is not on linesegment: abscissa={t}')
        return t

    def unit_direction_vector(self, abscissa=0.):
        """

        :param abscissa: defines where in the line_segement the unit
         direction vector is to be calculated
        :return: The unit direction vector of the LineSegement
        """
        direction_vector = self.direction_vector()
        direction_vector.normalize()
        return direction_vector

    def direction_vector(self, abscissa=0.):
        """
        :param abscissa: defines where in the line_segement
        direction vector is to be calculated
        :return: The direction vector of the LineSegement
        """
        return self.end - self.start

    def normal_vector(self, abscissa=0.):
        """
        :param abscissa: defines where in the line_segement
        normal vector is to be calculated
        :return: The normal vector of the LineSegement
        """
        return self.direction_vector(abscissa).normal_vector()

    def unit_normal_vector(self, abscissa=0.):
        """
        :param abscissa: defines where in the line_segement
        unit normal vector is to be calculated
        :return: The unit normal vector of the LineSegement
        """
        return self.unit_direction_vector(abscissa).normal_vector()

    def point_projection(self, point):
        p1, p2 = self.start, self.end
        u = p2 - p1
        norm_u = u.norm()
        t = (point - p1).dot(u) / norm_u ** 2
        projection = p1 + t * u

        return projection, t * norm_u

    def split(self, split_point):
        if split_point == self.start:
            return [None, self.copy()]
        elif split_point == self.end:
            return [self.copy(), None]
        else:
            return [self.__class__(self.start, split_point),
                    self.__class__(split_point, self.end)]

    def middle_point(self):
        return 0.5 * (self.start + self.end)

    def point_at_abscissa(self, abscissa):
        return self.start + self.unit_direction_vector() * abscissa


class BSplineCurve(Edge):
    _non_serializable_attributes = ['curve']

    def __init__(self,
                 degree: int,
                 control_points,
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
        self.name = name

        curve = BSpline.Curve()
        curve.degree = degree
        if weights is None:
            P = [[*point] for point in control_points]
            curve.ctrlpts = P
        else:
            Pw = [[*point * weights[i], weights[i]] for i, point in enumerate(control_points)]
            curve.ctrlptsw = Pw

        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot] * knot_multiplicities[i])
        curve.knotvector = knot_vector
        curve.delta = 0.01
        curve_points = curve.evalpts
        self.curve = curve

        self._length = None
        self.points = [getattr(volmdlr,
                               f'Point{self.__class__.__name__[-2::]}')(*p)
                       for p in curve_points]

        start = self.points[0]  # self.point_at_abscissa(0.)
        end = self.points[-1]  # self.point_at_abscissa(self.length())

        Edge.__init__(self, start, end, name=name)

    def reverse(self):
        '''
        reverse the bspline's direction by reversing its start and end points
        '''

        return self.__class__(degree=self.degree,
                              control_points=self.control_points[::-1],
                              knot_multiplicities=self.knot_multiplicities[::-1],
                              knots=self.knots[::-1],
                              weights=self.weights,
                              periodic=self.periodic)

    @classmethod
    def from_geomdl_curve(cls, curve):

        point_dimension = f'Point{cls.__name__[-2::]}'

        knots = list(sorted(set(curve.knotvector)))
        knot_multiplicities = [curve.knotvector.count(k) for k in knots]

        return cls(degree=curve.degree,
                   control_points=[getattr(volmdlr, point_dimension)(*p) for p in curve.ctrlpts],
                   knots=knots,
                   knot_multiplicities=knot_multiplicities)

    def length(self):
        if not self._length:
            self._length = length_curve(self.curve)
        return self._length

    def unit_direction_vector(self, abscissa: float):
        """
        :param abscissa: defines where in the BSplineCurve the
        unit direction vector is to be calculated
        :return: The unit direction vector of the BSplineCurve
        """

        direction_vector = self.direction_vector(abscissa)
        direction_vector.normalize()
        return direction_vector

    def middle_point(self):
        return self.point_at_abscissa(self.length() * 0.5)

    def abscissa(self, point, tol=1e-4):
        length = self.length()
        res = scp.optimize.least_squares(
            lambda u: (point - self.point_at_abscissa(u)).norm(),
            x0=npy.array(length / 2),
            bounds=([0], [length]),
            # ftol=tol / 10,
            # xtol=tol / 10,
            # loss='soft_l1'
        )

        if res.fun > tol:
            print('distance =', res.cost)
            print('res.fun:', res.fun)
            # ax = self.plot()
            # point.plot(ax=ax)
            # best_point = self.point_at_abscissa(res.x)
            # best_point.plot(ax=ax, color='r')
            raise ValueError('abscissa not found')
        return res.x[0]

    def split(self, point, tol=1e-5):
        if point.point_distance(self.start) < tol:
            return [None, self.copy()]
        elif point.point_distance(self.end) < tol:
            return [self.copy(), None]
        else:
            adim_abscissa = self.abscissa(point) / self.length()
            curve1, curve2 = split_curve(self.curve, adim_abscissa)

            return [self.__class__.from_geomdl_curve(curve1),
                    self.__class__.from_geomdl_curve(curve2)]

    def translation(self, offset):
        """
        BSplineCurve translation
        :param offset: translation vector (volmdlr.Vector2D/ volmdlr.Vector3D)
        :return: A new translated BSplineCurve
        """
        control_points = [point.translation(offset)
                          for point in self.control_points]
        return self.__class__(self.degree, control_points,
                              self.knot_multiplicities, self.knots,
                              self.weights, self.periodic)

    def translation_inplace(self, offset):
        """
        BSplineCurve translation. Object is updated inplace
        :param offset: translation vector (volmdlr.Vector2D/ volmdlr.Vector3D)
        """
        for point in self.control_points:
            point.translation_inplace(offset)

    def point_belongs(self, point, abs_tol=1e-10):
        '''
        check if a point belongs to the bspline_curve or not
        '''

        point_dimension = f'Point{self.__class__.__name__[-2::]}'

        def f(x):
            return (point - getattr(volmdlr, point_dimension)(*self.curve.evaluate_single(x))).norm()

        x = npy.linspace(0, 1, 5)
        x_init = []
        for xi in x:
            x_init.append(xi)

        for x0 in x_init:
            z = scp.optimize.least_squares(f, x0=x0, bounds=([0, 1]))
            if z.fun < abs_tol:
                return True
        return False

    def merge_with(self, bspline_curve):
        '''
        merge successives bspline_curves to define a new one
        '''

        point_dimension = f'Wire{self.__class__.__name__[-2::]}'
        wire = getattr(volmdlr.wires, point_dimension)(bspline_curve)
        ordered_wire = wire.order_wire()

        points, n = [], 10
        for primitive in ordered_wire.primitives:
            points.extend(primitive.polygon_points(n))
        points.pop(n + 1)

        return self.__class__.from_points_interpolation(points, min(self.degree, bspline_curve.degree))

    @classmethod
    def from_bsplines(cls, bsplines, discretization_points=10):
        '''
        define a bspline_curve using a list of bsplines
        '''
        point_dimension = f'Wire{cls.__name__[-2::]}'
        wire = getattr(volmdlr.wires, point_dimension)(bsplines)
        ordered_wire = wire.order_wire()

        points, degree = [], []
        for i, primitive in enumerate(ordered_wire.primitives):
            degree.append(primitive.degree)
            if i == 0:
                points.extend(primitive.polygon_points(discretization_points))
            else:
                points.extend(primitive.polygon_points(discretization_points)[1::])

        return cls.from_points_interpolation(points, min(degree))

    @classmethod
    def from_points_approximation(cls, points, degree, **kwargs):
        '''
        Bspline Curve approximation through points using least squares method
        It is better to specify the number of control points

        Parameters
        ----------
        points : volmdlr.Point
            data points
        degree: int
            degree of the output parametric curve

        Keyword Arguments:
            * ``centripetal``: activates centripetal parametrization method. *Default: False*
            * ``ctrlpts_size``: number of control points. *Default: len(points) - 1*

        Returns
        -------
        BSplineCurve

        '''

        curve = fitting.approximate_curve([[*point] for point in points], degree, **kwargs)
        return cls.from_geomdl_curve(curve)

    def tangent(self, position: float = 0.0):
        point, tangent = operations.tangent(self.curve, position,
                                            normalize=True)

        dimension = f'Vector{self.__class__.__name__[-2::]}'
        tangent = getattr(volmdlr, dimension)(*tangent)

        return tangent

    @classmethod
    def from_points_interpolation(cls, points, degree, periodic=False):

        curve = fitting.interpolate_curve([[*point] for point in points], degree)

        bsplinecurve = cls.from_geomdl_curve(curve)
        if not periodic:
            return bsplinecurve
        else:
            bsplinecurve.periodic = True
            return bsplinecurve


class Line2D(Line):
    """
    Define an infinite line given by two points.
    """

    def __init__(self, point1: volmdlr.Point2D,
                 point2: volmdlr.Point2D, *, name=''):
        self.points = [point1, point2]
        Line.__init__(self, point1, point2, name=name)

    def to_3d(self, plane_origin, x1, x2):
        p3D = [p.to_3d(plane_origin, x1, x2) for p in self.points]
        return Line2D(*p3D, self.name)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Line2D rotation
        :param center: rotation center
        :param angle: angle rotation
        :return: a new rotated Line2D
        """
        return Line2D(*[point.rotation(center, angle)
                        for point in self.points])

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        Line2D rotation. Object is updated inplace
        :param center: rotation center
        :param angle: rotation angle
        """
        for point in self.points:
            point.rotation_inplace(center, angle)

    def translation(self, offset: volmdlr.Vector2D):
        """
        Line2D translation
        :param offset: translation vector
        :return: A new translated Line2D
        """
        return Line2D(*[point.translation(offset) for point in self.points])

    def translation_inplace(self, offset: volmdlr.Vector2D):
        """
        Line2D translation. Object is updated inplace
        :param offset: translation vector
        """
        for point in self.points:
            point.translation_inplace(offset)

    def plot(self, ax=None, color='k', dashed=True):
        if ax is None:
            fig, ax = plt.subplots()

        if version.parse(_mpl_version) >= version.parse('3.3.2'):
            if dashed:
                ax.axline((self.point1.x, self.point1.y),
                          (self.point2.x, self.point2.y),
                          dashes=[30, 5, 10, 5],
                          color=color)
            else:
                ax.axline((self.point1.x, self.point1.y),
                          (self.point2.x, self.point2.y),
                          color=color)
        else:
            u = self.direction_vector()
            p3 = self.point1 - 3 * u
            p4 = self.point2 + 4 * u
            if dashed:
                ax.plot([p3[0], p4[0]], [p3[1], p4[1]], color=color,
                        dashes=[30, 5, 10, 5])
            else:
                ax.plot([p3[0], p4[0]], [p3[1], p4[1]], color=color)

        return ax

    def plot_data(self, edge_style=None):
        return plot_data.Line2D([self.point1.x, self.point1.y],
                                [self.point2.x, self.point2.y],
                                edge_style=edge_style)

    def line_intersections(self, line):

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
        else:
            return []

    def create_tangent_circle(self, point, other_line):
        """
        Computes the two circles that are tangent to 2 lines and intersect
        a point located on one of the two lines.
        """

        # point will be called I(x_I, y_I)
        # self will be (AB)
        # line will be (CD)

        if math.isclose(self.point_distance(point), 0, abs_tol=1e-10):
            I = volmdlr.Vector2D(point[0], point[1])
            A = volmdlr.Vector2D(self.points[0][0], self.points[0][1])
            B = volmdlr.Vector2D(self.points[1][0], self.points[1][1])
            C = volmdlr.Vector2D(other_line.points[0][0],
                                 other_line.points[0][1])
            D = volmdlr.Vector2D(other_line.points[1][0],
                                 other_line.points[1][1])

        elif math.isclose(other_line.point_distance(point), 0, abs_tol=1e-10):
            I = volmdlr.Vector2D(point[0], point[1])
            C = volmdlr.Vector2D(self.points[0][0], self.points[0][1])
            D = volmdlr.Vector2D(self.points[1][0], self.points[1][1])
            A = volmdlr.Vector2D(other_line.points[0][0],
                                 other_line.points[0][1])
            B = volmdlr.Vector2D(other_line.points[1][0],
                                 other_line.points[1][1])
        else:
            raise AttributeError("The point isn't on any of the two lines")

        # CHANGEMENT DE REPAIRE
        new_u = volmdlr.Vector2D((B - A))
        new_u.normalize()
        new_v = new_u.unit_normal_vector()
        new_basis = volmdlr.Frame2D(I, new_u, new_v)

        new_A = new_basis.new_coordinates(A)
        new_B = new_basis.new_coordinates(B)
        new_C = new_basis.new_coordinates(C)
        new_D = new_basis.new_coordinates(D)

        if new_C[1] == 0 and new_D[1] == 0:
            # Segments are on the same line: no solution
            return None, None

        elif math.isclose(self.unit_direction_vector().dot(
                other_line.unit_normal_vector()), 0, abs_tol=1e-06):
            # Parallel segments: one solution

            segments_distance = abs(new_C[1] - new_A[1])
            r = segments_distance / 2
            new_circle_center = volmdlr.Point2D(
                (0, npy.sign(new_C[1] - new_A[1]) * r))
            circle_center = new_basis.old_coordinates(new_circle_center)
            circle = volmdlr.wires.Circle2D(circle_center, r)

            return circle, None

        elif math.isclose(self.unit_direction_vector().dot(
                other_line.unit_direction_vector()), 0, abs_tol=1e-06):
            # Perpendicular segments: 2 solution
            line_AB = Line2D(volmdlr.Point2D(new_A), volmdlr.Point2D(new_B))
            line_CD = Line2D(volmdlr.Point2D(new_C), volmdlr.Point2D(new_D))
            new_pt_K = volmdlr.Point2D.line_intersection(line_AB, line_CD)

            r = abs(new_pt_K[0])
            new_circle_center1 = volmdlr.Point2D((0, r))
            new_circle_center2 = volmdlr.Point2D((0, -r))
            circle_center1 = new_basis.old_coordinates(new_circle_center1)
            circle_center2 = new_basis.old_coordinates(new_circle_center2)
            circle1 = volmdlr.wires.Circle2D(circle_center1, r)
            circle2 = volmdlr.wires.Circle2D(circle_center2, r)

            return circle1, circle2

        # =============================================================================
        # LES SEGMENTS SONT QUELCONQUES
        #   => 2 SOLUTIONS
        # =============================================================================
        else:

            line_AB = Line2D(volmdlr.Point2D(new_A), volmdlr.Point2D(new_B))
            line_CD = Line2D(volmdlr.Point2D(new_C), volmdlr.Point2D(new_D))
            new_pt_K = volmdlr.Point2D.line_intersection(line_AB, line_CD)
            pt_K = volmdlr.Point2D(new_basis.old_coordinates(new_pt_K))

            if pt_K == I:
                return None, None

            # CHANGEMENT DE REPERE:
            new_u2 = volmdlr.Vector2D(pt_K - I)
            new_u2.normalize()
            new_v2 = new_u2.normalVector(unit=True)
            new_basis2 = volmdlr.Frame2D(I, new_u2, new_v2)

            new_A = new_basis2.new_coordinates(A)
            new_B = new_basis2.new_coordinates(B)
            new_C = new_basis2.new_coordinates(C)
            new_D = new_basis2.new_coordinates(D)
            new_pt_K = new_basis2.new_coordinates(pt_K)

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

            new_circle_center1 = volmdlr.Point2D(0, -r1)
            new_circle_center2 = volmdlr.Point2D(0, r2)

            circle_center1 = new_basis2.old_coordinates(new_circle_center1)
            circle_center2 = new_basis2.old_coordinates(new_circle_center2)

            if new_basis.new_coordinates(circle_center1)[1] > 0:
                circle1 = volmdlr.wires.Circle2D(circle_center1, r1)
                circle2 = volmdlr.wires.Circle2D(circle_center2, r2)
            else:
                circle1 = volmdlr.wires.Circle2D(circle_center2, r2)
                circle2 = volmdlr.wires.Circle2D(circle_center1, r1)

            return circle1, circle2

    def cut_between_two_points(self, point1, point2):
        return LineSegment2D(point1, point2)

    def sort_points_along_line(self, points: List[volmdlr.Point2D]) -> List[
            volmdlr.Point2D]:
        most_distant_point = None
        farthest_distance = 0
        for i, point1 in enumerate(points):
            distances = []
            points_to_search = points[:i - 1] + points[i:]
            for point2 in points_to_search:
                distances.append(point1.point_distance(point2))
            max_point_distance = max(distances)
            farthest_point = points_to_search[
                distances.index(max_point_distance)]
            if max_point_distance > farthest_distance:
                most_distant_point = farthest_point
        list_points = [most_distant_point]
        points.remove(most_distant_point)
        distances_to_reference_point = {}
        for point in points:
            distances_to_reference_point[point] = \
                most_distant_point.point_distance(point)
        distances_to_reference_point = dict(
            sorted(distances_to_reference_point.items(),
                   key=lambda item: item[1]))
        list_points.extend(list(distances_to_reference_point.keys()))
        return list_points


class BSplineCurve2D(BSplineCurve):
    _non_serializable_attributes = ['curve']

    def __init__(self,
                 degree: int,
                 control_points: List[volmdlr.Point2D],
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

    def bounding_rectangle(self):
        points = self.discretization_points()
        points_x = [p.x for p in points]
        points_y = [p.y for p in points]

        return (min(points_x), max(points_x),
                min(points_y), max(points_y))

    def length(self):
        return length_curve(self.curve)

    def tangent(self, position: float = 0.0):
        _, tangent = operations.tangent(self.curve, position,
                                        normalize=True)
        tangent = volmdlr.Point2D(tangent[0], tangent[1])
        return tangent

    def point_at_abscissa(self, abscissa):
        l = self.length()
        adim_abs = max(min(abscissa / l, 1.), 0.)
        return volmdlr.Point2D(*self.curve.evaluate_single(adim_abs))

    def direction_vector(self, abscissa: float):
        """
        :param abscissa: defines where in the BSplineCurve2D the
        direction vector is to be calculated
        :return: The direection vector vector of the BSplineCurve2D
        """
        return self.tangent(abscissa)

    def normal_vector(self, abscissa: float):
        """
        :param abscissa: defines where in the BSplineCurve2D the
        normal vector is to be calculated
        :return: The normal vector of the BSplineCurve2D
        """
        tangent_vector = self.tangent(abscissa)
        normal_vector = tangent_vector.normal_vector()
        return normal_vector

    def unit_normal_vector(self, abscissa: float):
        """
        :param abscissa: defines where in the BSplineCurve2D the
        unit normal vector is to be calculated
        :return: The unit normal vector of the BSplineCurve2D
        """
        normal_vector = self.normal_vector(abscissa)
        normal_vector.normalize()
        return normal_vector

    def straight_line_area(self):
        points = self.discretization_points(number_points=100)
        x = [point.x for point in points]
        y = [point.y for point in points]
        x1 = [x[-1]] + x[0:-1]
        y1 = [y[-1]] + y[0:-1]
        return 0.5 * abs(sum(i * j for i, j in zip(x, y1))
                         - sum(i * j for i, j in zip(y, x1)))

    def straight_line_center_of_mass(self):
        polygon_points = self.discretization_points(number_points=100)
        cog = volmdlr.O2D
        for point in polygon_points:
            cog += point
        cog = cog / len(polygon_points)
        return cog

    def plot(self, ax=None, color='k', alpha=1, plot_points=False):
        if ax is None:
            _, ax = plt.subplots()

        # self.curve.delta = 0.01
        # points = [volmdlr.Point2D(px, py) for (px, py) in self.curve.evalpts]
        l = self.length()
        points = [self.point_at_abscissa(l * i / 50) for i in range(51)]

        xp = [p.x for p in points]
        yp = [p.y for p in points]
        ax.plot(xp, yp, color=color, alpha=alpha)

        return ax

    def to_3d(self, plane_origin, x1, x2):
        control_points3D = [p.to_3d(plane_origin, x1, x2) for p in
                            self.control_points]
        return BSplineCurve3D(self.degree, control_points3D,
                              self.knot_multiplicities, self.knots,
                              self.weights, self.periodic)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        BSplineCurve2D rotation
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
        BSplineCurve2D rotation. Object is updated inplace
        :param center: rotation center
        :param angle: rotation angle
        """
        for point in self.control_points:
            point.rotation_inplace(center, angle)

    def line_intersections(self, line2d: Line2D):
        polygon_points = self.discretization_points(number_points=200)
        list_intersections = []
        length = self.length()
        initial_abscissa = 0
        for point1, point2 in zip(polygon_points[:-1], polygon_points[1:]):
            linesegment = LineSegment2D(point1, point2)
            intersections = linesegment.line_intersections(line2d)
            initial_abscissa += linesegment.length()
            if intersections:
                if initial_abscissa < length * 0.1:
                    list_abcissas = [initial_abscissa * n for n in
                                     npy.linspace(0, 1, 100)]
                else:
                    list_abcissas = [initial_abscissa * n for n in
                                     npy.linspace(0.9, 1, 100)]
                distance = npy.inf
                for abscissa in list_abcissas:
                    point_in_curve = self.point_at_abscissa(abscissa)
                    dist = point_in_curve.point_distance(intersections[0])
                    if dist < distance:
                        distance = dist
                        intersection = point_in_curve
                list_intersections.append(intersection)
        return list_intersections

    def line_crossings(self, line2d: Line2D):
        polygon_points = self.discretization_points(number_points=50)
        crossings = []
        for p1, p2 in zip(polygon_points[:-1], polygon_points[1:]):
            l = LineSegment2D(p1, p2)
            crossings.extend(l.line_crossings(line2d))
        return crossings

    def to_wire(self, n: int):
        '''
        convert a bspline curve to a wire2d defined with 'n' line_segments
        '''

        u = npy.linspace(0, 1, num=n + 1).tolist()
        points = []
        for u0 in u:
            p = self.curve.evaluate_single(u0)
            points.append(volmdlr.Point2D(p[0], p[1]))

        return volmdlr.wires.Wire2D.from_points(points)

    def reverse(self):
        '''
        reverse the bspline's direction by reversing its start and end points
        '''

        return self.__class__(degree=self.degree,
                              control_points=self.control_points[::-1],
                              knot_multiplicities=self.knot_multiplicities[::-1],
                              knots=self.knots[::-1],
                              weights=self.weights,
                              periodic=self.periodic)

    def point_distance(self, point):
        distance = math.inf
        polygon_points = self.discretization_points(number_points=20)
        for p1, p2 in zip(polygon_points[:-1], polygon_points[1:]):
            line = LineSegment2D(p1, p2)
            dist = line.point_distance(point)
            if dist < distance:
                distance = dist
        return distance

    def nearest_point_to(self, point):
        '''
        find out the nearest point on the linesegment to point
        '''

        points = self.polygon_points(500)
        return point.nearest_point(points)

    def linesegment_intersections(self, linesegment):
        results = self.line_intersections(linesegment.to_line())
        intersections_points = []
        for result in results:
            if linesegment.point_belongs(result, 1e-6):
                intersections_points.append(result)
        return intersections_points

    def axial_symmetry(self, line):
        '''
        finds out the symmetric bsplinecurve2d according to a line
        '''

        points_symmetry = [point.axial_symmetry(line) for point in self.control_points]

        return self.__class__(degree=self.degree,
                              control_points=points_symmetry,
                              knot_multiplicities=self.knot_multiplicities[::-1],
                              knots=self.knots[::-1],
                              weights=self.weights,
                              periodic=self.periodic)


class BezierCurve2D(BSplineCurve2D):

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
    Define a line segment limited by two points
    """

    def __init__(self, start: volmdlr.Point2D, end: volmdlr.Point2D, *, name: str = ''):
        if start == end:
            raise NotImplementedError
        LineSegment.__init__(self, start, end, name=name)

    def __hash__(self):
        return self._data_hash()

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

    def point_belongs(self, point, abs_tol=1e-6):
        distance = self.start.point_distance(point) + self.end.point_distance(
            point)
        if math.isclose(distance, self.length(), abs_tol=abs_tol) and\
                math.isclose(self.point_distance(point), 0.0, abs_tol=abs_tol):
            return True
        return False

    def bounding_rectangle(self):
        return (min(self.start.x, self.end.x), max(self.start.x, self.end.x),
                min(self.start.y, self.end.y), max(self.start.y, self.end.y))

    def straight_line_area(self):
        return 0.

    def straight_line_second_moment_area(self, point: volmdlr.Point2D):
        return 0, 0, 0

    def straight_line_center_of_mass(self):
        return 0.5 * (self.start + self.end)

    def point_distance(self, point, return_other_point=False):
        """
        Computes the distance of a point to segment of line
        """
        distance, point = volmdlr.core_compiled.LineSegment2DPointDistance(
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
        point = volmdlr.Point2D.line_intersection(self, line)
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
            vector1 = self.start - line.point1
            vector2 = self.start - line.point2
            vector3 = self.end - line.point1
            vector4 = self.end - line.point2
            if math.isclose(vector1.cross(vector2), 0, abs_tol=1e-6):
                return [self.start]
            if math.isclose(vector3.cross(vector4), 0, abs_tol=1e-6):
                return [self.end]
        return []

    def linesegment_intersections(self, linesegment: 'LineSegment2D'):
        """
        touching linesegments does not intersect
        """
        point = volmdlr.Point2D.line_intersection(self, linesegment)
        # TODO: May be these commented conditions should be used for linesegment_crossings
        if point:  # and (point != self.start) and (point != self.end):
            point_projection1, _ = self.point_projection(point)
            if point_projection1 is None:
                return []

            point_projection2, _ = linesegment.point_projection(point)
            if point_projection2 is None:
                return []

            return [point_projection1]
        else:
            return []

    def line_crossings(self, line: 'Line2D'):
        if self.direction_vector().is_colinear_to(line.direction_vector()):
            return []
        else:
            line_intersection = self.line_intersections(line)
            if line_intersection and (line_intersection[0] == self.end or line_intersection[0] == self.start):
                return []
            return line_intersection

    def linesegment_crossings(self, linesegment: 'LineSegment2D'):
        if self.direction_vector().is_colinear_to(
                linesegment.direction_vector()):
            return []
        else:
            return self.linesegment_intersections(linesegment)

    def plot(self, ax=None, color='k', alpha=1, arrow=False, width=None,
             plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()

        p1, p2 = self.start, self.end
        if arrow:
            if plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color,
                        alpha=alpha, style='o-')
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color,
                        alpha=alpha)

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
                        marker='o', linewidth=width, alpha=alpha)
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color,
                        linewidth=width, alpha=alpha)
        return ax

    def to_3d(self, plane_origin, x1, x2):
        start = self.start.to_3d(plane_origin, x1, x2)
        end = self.end.to_3d(plane_origin, x1, x2)
        return LineSegment3D(start, end, name=self.name)

    def reverse(self):
        return LineSegment2D(self.end.copy(), self.start.copy())

    def to_line(self):
        return Line2D(self.start, self.end)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        LineSegment2D rotation
        :param center: rotation center
        :param angle: angle rotation
        :return: a new rotated LineSegment2D
        """
        return LineSegment2D(self.start.rotation(center, angle),
                             self.end.rotation(center, angle))

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        LineSegment2D rotation. Object is updated inplace
        :param center: rotation center
        :param angle: rotation angle
        """
        for point in [self.start, self.end]:
            point.rotation_inplace(center, angle)

    def translation(self, offset: volmdlr.Vector2D):
        """
        LineSegment2D translation
        :param offset: translation vector
        :return: A new translated LineSegment2D
        """
        return LineSegment2D(self.start.translation(offset),
                             self.end.translation(offset))

    def translation_inplace(self, offset: volmdlr.Vector2D):
        """
        LineSegment2D translation. Object is updated inplace
        :param offset: translation vector
        """
        for point in [self.start, self.end]:
            point.translation_inplace(offset)

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes vector frame_mapping and return a new LineSegment2D
        side = 'old' or 'new'
        """
        if side == 'old':
            new_start = frame.old_coordinates(self.start)
            new_end = frame.old_coordinates(self.end)
        elif side == 'new':
            new_start = frame.new_coordinates(self.start)
            new_end = frame.new_coordinates(self.end)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        return LineSegment2D(new_start, new_end)

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes vector frame_mapping and the object is updated inplace
        side = 'old' or 'new'
        """
        if side == 'old':
            new_start = frame.old_coordinates(self.start)
            new_end = frame.old_coordinates(self.end)
        elif side == 'new':
            new_start = frame.new_coordinates(self.start)
            new_end = frame.new_coordinates(self.end)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        self.start = new_start
        self.end = new_end

    def plot_data(self, edge_style: plot_data.EdgeStyle = None):
        return plot_data.LineSegment2D([self.start.x, self.start.y],
                                       [self.end.x, self.end.y],
                                       edge_style=edge_style)

    def CreateTangentCircle(self, point, other_line):
        circle1, circle2 = Line2D.CreateTangentCircle(other_line, point, self)
        if circle1 is not None:
            point_J1, curv_abs1 = Line2D.point_projection(self, circle1.center)
            if curv_abs1 < 0. or curv_abs1 > self.length():
                circle1 = None
        if circle2 is not None:
            point_J2, curv_abs2 = Line2D.point_projection(self, circle2.center)
            if curv_abs2 < 0. or curv_abs2 > self.length():
                circle2 = None
        return circle1, circle2

    def infinite_primitive(self, offset):
        n = self.normal_vector()
        offset_point_1 = self.start + offset * \
            n

        offset_point_2 = self.end + offset * \
            n

        return Line2D(offset_point_1, offset_point_2)

    def polygon_points(self, discretization_resolution: int):
        warnings.warn('polygon_points is deprecated,\
        please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(number_points=discretization_resolution)

    def to_wire(self, n: int):
        '''
        convert a linesegment2d to a wire2d defined with 'n' line_segments
        '''

        points = self.discretise(n)
        return volmdlr.wires.Wire2D.from_points(points)

    def nearest_point_to(self, point):
        '''
        find out the nearest point on the linesegment to point
        '''

        points = self.discretization_points(number_points=500)
        return point.nearest_point(points)

    def axial_symmetry(self, line):
        '''
        finds out the symmetric linesegment2d according to a line
        '''

        points_symmetry = [point.axial_symmetry(line) for point in [self.start, self.end]]

        return self.__class__(points_symmetry[0], points_symmetry[1])


class Arc(Edge):
    def __init__(self, start,
                 end,
                 interior,
                 name: str = ''):
        Edge.__init__(self, start=start, end=end, name=name)
        self.interior = interior
        self._utd_clockwise_and_trigowise_paths = False
        self._clockwise_and_trigowise_paths = None
        self._radius = None

    @property
    def center(self):
        """
        Gets the arc's center
        :return: The center of the arc
        """
        raise NotImplementedError(
            'the property method center must be overloaded by subclassing'
            'class if not a given parameter')

    @property
    def angle(self):
        """
        Gets the angle of the arc
        :return: The angle of the arc
        """
        return NotImplementedError(
            'the property method angle must be overloaded by subclassing'
            'class if not a given parameter')

    @property
    def is_trigo(self):
        """
        Verifies if arc is trigowise or clockwise
        :return: True if trigowise or False otherwise
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
        Calculates the length of the Arc, with its radius and it arc angle
        :return: the length fo the Arc
        """
        return self.radius * abs(self.angle)

    def point_at_abscissa(self, abscissa):
        if self.is_trigo:
            return self.start.rotation(self.center,
                                       abscissa / self.radius)
        else:
            return self.start.rotation(self.center,
                                       -abscissa / self.radius)

    @staticmethod
    def get_clockwise_and_trigowise_paths(radius_1, radius_2, radius_i):
        """
        :param radius_1: radius from center to start point
        :param radius_2: radius form center ro end point
        :param radius_i: radius from center to interior point
        :return: the clockwise and trigowise paths
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
        return self.point_at_abscissa(0.5 * self.length())

    def point_distance(self, point):
        points = self.discretization_points(angle_resolution=100)
        return point.point_distance(point.nearest_point(points))

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = None):
        """
        discretize a Edge to have "n" points
        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Usefull
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


class Arc2D(Arc):
    """
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

    @property
    def center(self):
        if not self._center:
            self._center = self.get_center()
        return self._center

    def get_center(self):
        xi, yi = self.interior.x, self.interior.y
        xe, ye = self.end.x, self.end.y
        xs, ys = self.start.x, self.start.y
        try:
            A = volmdlr.Matrix22(2 * (xs - xi), 2 * (ys - yi),
                                 2 * (xs - xe), 2 * (ys - ye))
            b = - volmdlr.Vector2D(xi ** 2 + yi ** 2 - xs ** 2 - ys ** 2,
                                   xe ** 2 + ye ** 2 - xs ** 2 - ys ** 2)
            inv_A = A.inverse()
            x = inv_A.vector_multiplication(b)
            center = volmdlr.Point2D(x.x, x.y)
        except ValueError:
            A = npy.array([[2 * (xs - xi), 2 * (ys - yi)],
                           [2 * (xs - xe), 2 * (ys - ye)]])
            b = - npy.array([xi ** 2 + yi ** 2 - xs ** 2 - ys ** 2,
                             xe ** 2 + ye ** 2 - xs ** 2 - ys ** 2])
            center = volmdlr.Point2D(*npy.linalg.solve(A, b))
        return center

    @property
    def is_trigo(self):
        if not self._is_trigo:
            self._is_trigo = self.get_arc_direction()
        return self._is_trigo

    @property
    def clockwise_and_trigowise_paths(self):
        if not self._clockwise_and_trigowise_paths:
            radius_1 = self.start - self.center
            radius_2 = self.end - self.center
            radius_i = self.interior - self.center
            self._clockwise_and_trigowise_paths =\
                self.get_clockwise_and_trigowise_paths(radius_1,
                                                       radius_2,
                                                       radius_i)
            self._utd_clockwise_and_trigowise_paths = True
        return self._clockwise_and_trigowise_paths

    def get_arc_direction(self):
        clockwise_path, trigowise_path =\
            self.clockwise_and_trigowise_paths
        if clockwise_path > trigowise_path:
            return True
        return False

    @property
    def angle(self):
        if not self._angle:
            self._angle = self.get_angle()
        return self._angle

    def get_angle(self):
        clockwise_path, trigowise_path = \
            self.clockwise_and_trigowise_paths
        if self.is_trigo:
            return trigowise_path
        return clockwise_path

    def _get_points(self):
        return [self.start, self.interior, self.end]

    points = property(_get_points)

    def point_distance(self, point):
        vector_start = self.start - self.center
        vector_point = point - self.center
        vector_end = self.end - self.center
        if self.is_trigo:
            vector_start, vector_end = vector_end, vector_start
        arc_angle = volmdlr.core.clockwise_angle(vector_start, vector_end)
        point_angle = volmdlr.core.clockwise_angle(vector_start, vector_point)
        if point_angle <= arc_angle:
            return abs(
                LineSegment2D(point, self.center).length() - self.radius)
        else:
            return min(LineSegment2D(point, self.start).length(),
                       LineSegment2D(point, self.end).length())

    def point_belongs(self, point2d, abs_tol=1e-10):
        """
        check if a Point2D belongs to the Arc2D
        """
        vector_start = self.start - self.center
        vector_end = self.end - self.center
        vector_point = point2d - self.center
        r1 = vector_start.norm()
        cp = vector_point.norm()
        if math.isclose(cp, r1, abs_tol=abs_tol):
            if self.get_arc_direction():
                arc_angle = - volmdlr.core.clockwise_angle(vector_start,
                                                           vector_end)
                point_angle = - volmdlr.core.clockwise_angle(vector_start,
                                                             vector_point)

            else:
                arc_angle = volmdlr.core.clockwise_angle(vector_start,
                                                         vector_end)
                point_angle = volmdlr.core.clockwise_angle(vector_start,
                                                           vector_point)
            if point_angle <= arc_angle:
                return True
        return False

    # def to_circle(self):
    #     return volmdlr.wires.Circle2D(self.center, self.radius)

    def to_full_arc_2d(self):
        return FullArc2D(center=self.center,
                         start_end=self.point_at_abscissa(0),
                         name=self.name)

    def line_intersections(self, line2d: Line2D):
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
        full_arc_2d = self.to_full_arc_2d()
        fa2d_intersection_points = full_arc_2d.linesegment_intersections(
            linesegment2d)
        intersection_points = []
        for pt in fa2d_intersection_points:
            if self.point_belongs(pt):
                intersection_points.append(pt)
        return intersection_points

    def abscissa(self, point2d: volmdlr.Point2D, tol=1e-9):
        if point2d.point_distance(self.start) < tol:
            return 0
        if point2d.point_distance(self.end) < tol:
            return self.length()

        p = point2d - self.center
        u = self.start - self.center
        u.normalize()
        if self.is_trigo:
            v = u.normal_vector()
        else:
            v = -u.normal_vector()

        x, y = p.dot(u), p.dot(v)
        theta = math.atan2(y, x)
        if theta < -tol or theta > self.angle + tol:
            raise ValueError('Point not in arc')

        if theta < 0:
            return 0.
        if theta > self.angle:
            return self.angle * self.radius

        return self.radius * theta

    def direction_vector(self, abscissa: float):
        """
        :param abscissa: defines where in the Arc2D the
        direction vector is to be calculated
        :return: The direction vector of the Arc2D
        """
        return -self.normal_vector(abscissa=abscissa).normal_vector()

    def unit_direction_vector(self, abscissa: float):
        """
        :param abscissa: defines where in the Arc2D the
        unit direction vector is to be calculated
        :return: The unit direction vector of the Arc2D
        """
        direction_vector = self.direction_vector(abscissa)
        direction_vector.normalize()
        return direction_vector

    def normal_vector(self, abscissa: float):
        """
        :param abscissa: defines where in the Arc2D the
        normal vector is to be calculated
        :return: The normal vector of the Arc2D
        """
        point = self.point_at_abscissa(abscissa)
        # if self.is_trigo:
        normal_vector = self.center - point
        # else:
        #     normal_vector = point - self.center
        return normal_vector

    def unit_normal_vector(self, abscissa: float):
        """
        :param abscissa: defines where in the Arc2D the
        unit normal vector is to be calculated
        :return: The unit normal vector of the Arc2D
        """
        normal_vector = self.normal_vector(abscissa)
        normal_vector.normalize()
        return normal_vector

    def area(self):
        return self.radius ** 2 * self.angle / 2

    def center_of_mass(self):
        #        u=self.middle.vector-self.center.vector
        u = self.middle_point() - self.center
        u.normalize()
        # alpha = abs(self.angle)
        return self.center + 4 / (3 * self.angle) * self.radius * math.sin(
            self.angle * 0.5) * u

    def bounding_rectangle(self):
        # TODO: Enhance this!!!
        return (self.center.x - self.radius, self.center.x + self.radius,
                self.center.y - self.radius, self.center.y + self.radius)

    def straight_line_area(self):
        if self.angle >= math.pi:
            angle = volmdlr.TWO_PI - self.angle
            area = math.pi * self.radius ** 2 - 0.5 * self.radius ** 2 * (
                angle - math.sin(angle))
        else:
            angle = self.angle
            area = 0.5 * self.radius ** 2 * (angle - math.sin(angle))

        if self.is_trigo:
            return area
        else:
            return -area

    def straight_line_second_moment_area(self, point: volmdlr.Point2D):

        if self.angle2 < self.angle1:
            angle2 = self.angle2 + volmdlr.TWO_PI

        else:
            angle2 = self.angle2
        angle1 = self.angle1

        # Full arc section
        Ix1 = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
            math.sin(2 * angle1) - math.sin(2 * angle2)))
        Iy1 = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
            math.sin(2 * angle2) - math.sin(2 * angle1)))
        Ixy1 = self.radius ** 4 / 8 * (
            math.cos(angle1) ** 2 - math.cos(angle2) ** 2)

        # Triangle
        xi, yi = (self.start - self.center)
        xj, yj = (self.end - self.center)
        Ix2 = (yi ** 2 + yi * yj + yj ** 2) * (xi * yj - xj * yi) / 12.
        Iy2 = (xi ** 2 + xi * xj + xj ** 2) * (xi * yj - xj * yi) / 12.
        Ixy2 = (xi * yj + 2 * xi * yi + 2 * xj * yj + xj * yi) * (
            xi * yj - xj * yi) / 24.
        if Ix2 < 0.:
            Ix2, Iy2, Ixy2 = -Ix2, -Iy2, -Ixy2
        if self.angle < math.pi:
            if self.is_trigo:
                Ix = Ix1 - Ix2
                Iy = Iy1 - Iy2
                Ixy = Ixy1 - Ixy2
            else:
                Ix = Ix2 - Ix1
                Iy = Iy2 - Iy1
                Ixy = Ixy2 - Ixy1
        else:
            # print('Ixy12', Ixy1, Ixy2)
            if self.is_trigo:
                Ix = Ix1 + Ix2
                Iy = Iy1 + Iy2
                Ixy = Ixy1 + Ixy2
            else:
                Ix = -Ix2 - Ix1
                Iy = -Iy2 - Iy1
                Ixy = -Ixy2 - Ixy1

        return volmdlr.geometry.huygens2d(Ix, Iy, Ixy,
                                          self.straight_line_area(),
                                          self.center,
                                          point)

    def straight_line_center_of_mass(self):
        if self.angle == math.pi:
            return self.center_of_mass()

        u = self.middle_point() - self.center
        u.normalize()
        if self.angle >= math.pi:
            u = -u
        bissec = Line2D(self.center, self.center + u)
        string = Line2D(self.start, self.end)
        p = volmdlr.Point2D.line_intersection(bissec, string)
        a = p.point_distance(self.start)
        h = p.point_distance(self.center)
        triangle_area = h * a
        # alpha = abs(self.angle)
        triangle_cog = self.center + 2 / 3. * h * u
        if self.angle < math.pi:
            cog = (
                self.center_of_mass() * self.area() - triangle_area * triangle_cog) / abs(
                self.straight_line_area())
        else:
            cog = (
                self.center_of_mass() * self.area() + triangle_area * triangle_cog) / abs(
                self.straight_line_area())

        # ax = self.plot()
        # bissec.plot(ax=ax, color='grey')
        # self.center.plot(ax=ax)
        # string.plot(ax=ax, color='grey')
        # triangle_cog.plot(ax=ax, color='green')
        # self.center_of_mass().plot(ax=ax, color='red')
        #
        # cog_line = Line2D(volmdlr.O2D, self.center_of_mass()*self.area()-triangle_area*triangle_cog)
        # cog_line.plot(ax=ax)
        #
        # cog.plot(ax=ax, color='b')
        # ax.set_aspect('equal')
        return cog

    def plot(self, ax=None, color='k', alpha=1, plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()

        if plot_points:
            for p in [self.center, self.start, self.interior, self.end]:
                p.plot(ax=ax, color=color, alpha=alpha)

        ax.add_patch(matplotlib.patches.Arc((self.center.x, self.center.y), 2 * self.radius,
                                            2 * self.radius, angle=0,
                                            theta1=self.angle1 * 0.5 / math.pi * 360,
                                            theta2=self.angle2 * 0.5 / math.pi * 360,
                                            color=color,
                                            alpha=alpha))
        return ax

    def to_3d(self, plane_origin, x, y):
        ps = self.start.to_3d(plane_origin, x, y)
        pi = self.interior.to_3d(plane_origin, x, y)
        pe = self.end.to_3d(plane_origin, x, y)

        return volmdlr.edges.Arc3D(ps, pi, pe, name=self.name)

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Arc2D rotation
        :param center: rotation center
        :param angle: angle rotation
        :return: a new rotated Arc2D
        """
        return Arc2D(*[point.rotation(center, angle,) for point in
                       [self.start, self.interior, self.end]])

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        Arc2D rotation. Object is updated inplace
        :param center: rotation center
        :param angle: rotation angle
        """
        self.start.rotation_inplace(center, angle)
        self.interior.rotation_inplace(center, angle)
        self.end.rotation_inplace(center, angle)
        self._angle = None
        self._is_trigo = None
        self._center = None
        self._clockwise_and_trigowise_paths = None

    def translation(self, offset: volmdlr.Vector2D):
        """
        Arc2D translation
        :param offset: translation vector
        :return: A new translated Arc2D
        """
        return Arc2D(*[point.translation(offset) for point in
                       [self.start, self.interior, self.end]])

    def translation_inplace(self, offset: volmdlr.Vector2D):
        """
        Arc2D translation. Object is updated inplace
        :param offset: translation vector
        """
        self.start.translation_inplace(offset)
        self.interior.translation_inplace(offset)
        self.end.translation_inplace(offset)
        self._angle = None
        self._is_trigo = None
        self._center = None
        self._clockwise_and_trigowise_paths = None

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes vector frame_mapping and return a new Arc2D
        side = 'old' or 'new'
        """
        return Arc2D(*[point.frame_mapping(frame, side) for point in
                       [self.start, self.interior, self.end]])

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes vector frame_mapping and the object is updated inplace
        side = 'old' or 'new'
        """
        self.__init__(*[point.frame_mapping(frame, side) for point in
                        [self.start, self.interior, self.end]])

    def second_moment_area(self, point):
        """
        Second moment area of part of disk
        """
        if self.angle2 < self.angle1:
            angle2 = self.angle2 + volmdlr.TWO_PI

        else:
            angle2 = self.angle2
        angle1 = self.angle1

        Ix = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
            math.sin(2 * angle1) - math.sin(2 * angle2)))
        Iy = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
            math.sin(2 * angle2) - math.sin(2 * angle1)))
        Ixy = self.radius ** 4 / 8 * (
            math.cos(angle1) ** 2 - math.cos(angle2) ** 2)
        # Ic = npy.array([[Ix, Ixy], [Ixy, Iy]])

        # Must be computed at center, so huygens related to center
        return volmdlr.geometry.huygens2d(Ix, Iy, Ixy, self.area(),
                                          self.center, point)

    def plot_data(self, edge_style: plot_data.EdgeStyle = None,
                  anticlockwise: bool = None):

        list_node = self.discretization_points()
        data = []
        for nd in list_node:
            data.append({'x': nd.x, 'y': nd.y})
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

    def split(self, split_point: volmdlr.Point2D):
        abscissa = self.abscissa(split_point)

        return [Arc2D(self.start,
                      self.point_at_abscissa(0.5 * abscissa),
                      split_point),
                Arc2D(split_point,
                      self.point_at_abscissa((self.abscissa(self.end)
                                              - abscissa) * 0.5 + abscissa),
                      self.end)
                ]

    def infinite_primitive(self, offset):

        if not self.is_trigo:
            radius = self.radius + offset
        else:
            radius = self.radius - offset

        return FullArc2D(self.center,
                         self.center + radius * volmdlr.Point2D(1, 0.))

    def complementary(self):

        interior = self.middle_point().rotation(self.center, math.pi)
        return Arc2D(self.start, interior, self.end)

    def to_wire(self, angle_resolution: float = 10.):
        '''
        convert an arc to a wire2d defined with line_segments
        '''

        return volmdlr.wires.Wire2D.from_points(self.polygon_points(angle_resolution))

    def axial_symmetry(self, line):
        '''
        finds out the symmetric arc2d according to a line
        '''

        points_symmetry = [point.axial_symmetry(line) for point in [self.start, self.interior, self.end]]

        return self.__class__(start=points_symmetry[0],
                              interior=points_symmetry[1],
                              end=points_symmetry[2])


class FullArc2D(Arc2D):
    """
    An edge that starts at start_end, ends at the same point after having described
    a circle
    """

    def __init__(self, center: volmdlr.Point2D, start_end: volmdlr.Point2D,
                 name: str = ''):
        self.__center = center
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
        dict_['name'] = self.name
        return dict_

    def copy(self, *args, **kwargs):
        return FullArc2D(self.center.copy(), self.start.copy())

    @classmethod
    def dict_to_object(cls, dict_, global_dict=None, pointers_memo: Dict[str, Any] = None, path: str = '#'):
        center = volmdlr.Point2D.dict_to_object(dict_['center'])
        start_end = volmdlr.Point2D.dict_to_object(dict_['start_end'])

        return cls(center, start_end, name=dict_['name'])

    def __hash__(self):
        return hash(self.radius)
        # return hash(self.center) + 5*hash(self.start)

    def __eq__(self, other_arc):
        if self.__class__.__name__ != other_arc.__class__.__name__:
            return False
        return (self.center == other_arc.center) \
            and (self.start_end == other_arc.start_end)

    def straight_line_area(self):
        area = self.area()
        return area

    def center_of_mass(self):
        return self.center

    def straight_line_center_of_mass(self):
        return self.center_of_mass()

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
        self._center.rotation(center, angle, False)
        self.start.rotation(center, angle, False)
        self.interior.rotation(center, angle, False)
        self.end.rotation(center, angle, False)

    def translation(self, offset: volmdlr.Vector2D):
        new_center = self._center.translation(offset)
        new_start_end = self.start.translation(offset)
        return FullArc2D(new_center, new_start_end)

    def translation_inplace(self, offset: volmdlr.Vector2D):
        self._center.translation_inplace(offset)
        self.start.translation_inplace(offset)
        self.end.translation_inplace(offset)
        self.interior.translation_inplace(offset)

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        side = 'old' or 'new'
        """
        return FullArc2D(*[point.frame_mapping(frame, side) for point in
                           [self._center, self.start]])

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        [p.frame_mapping_inplace(frame, side) for p in
         [self._center, self.start, self.end, self.interior]]

    def polygonization(self):
        return volmdlr.wires.ClosedPolygon2D(self.discretization_points(angle_resolution=15))

    def plot(self, ax=None, color='k', alpha=1, plot_points=False,
             linestyle='-', linewidth=1):
        if ax is None:
            fig, ax = plt.subplots()

        if self.radius > 0:
            ax.add_patch(matplotlib.patches.Arc((self.center.x, self.center.y),
                                                2 * self.radius,
                                                2 * self.radius,
                                                angle=0,
                                                theta1=0,
                                                theta2=360,
                                                color=color,
                                                linestyle=linestyle,
                                                linewidth=linewidth))
        if plot_points:
            ax.plot([self.start.x], [self.start.y], 'o',
                    color=color, alpha=alpha)
        return ax

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
            if line2d.start == self.center:
                pt1 = line2d.end
                vec = line2d.start - line2d.end
            else:
                pt1 = line2d.start
                vec = line2d.end - line2d.start
        except AttributeError:
            if line2d.point1 == self.center:
                pt1 = line2d.point2
                vec = line2d.point1 - line2d.point2
            else:
                pt1 = line2d.point1
                vec = line2d.point2 - line2d.point1
        a = vec.dot(vec)
        b = 2 * vec.dot(pt1 - self.center)
        c = pt1.dot(pt1) + self.center.dot(self.center) \
            - 2 * pt1.dot(self.center) - self.radius ** 2

        disc = b ** 2 - 4 * a * c
        if math.isclose(disc, 0., abs_tol=tol):
            t1 = -b / (2 * a)
            return [pt1 + t1 * vec]

        elif disc > 0:
            sqrt_disc = math.sqrt(disc)
            t1 = (-b + sqrt_disc) / (2 * a)
            t2 = (-b - sqrt_disc) / (2 * a)
            return [pt1 + t1 * vec,
                    pt1 + t2 * vec]

        return []

    def linesegment_intersections(self, linesegment2d: LineSegment2D,
                                  tol=1e-9):
        try:
            if linesegment2d.start == self.center:
                pt1 = linesegment2d.end
                vec = linesegment2d.start - linesegment2d.end
            else:
                pt1 = linesegment2d.start
                vec = linesegment2d.end - linesegment2d.start
        except AttributeError:
            if linesegment2d.point1 == self.center:
                pt1 = linesegment2d.point2
                vec = linesegment2d.point1 - linesegment2d.point2
            else:
                pt1 = linesegment2d.point1
                vec = linesegment2d.point2 - linesegment2d.point1
        a = vec.dot(vec)
        b = 2 * vec.dot(pt1 - self.center)
        c = pt1.dot(pt1) + self.center.dot(self.center) \
            - 2 * pt1.dot(self.center) - self.radius ** 2

        disc = b ** 2 - 4 * a * c
        if math.isclose(disc, 0., abs_tol=tol):
            t1 = -b / (2 * a)
            points = [pt1 + t1 * vec]
            if linesegment2d.point_belongs(points[0]):
                return points
            return []

        elif disc > 0:
            sqrt_disc = math.sqrt(disc)
            t1 = (-b + sqrt_disc) / (2 * a)
            t2 = (-b - sqrt_disc) / (2 * a)
            points = [pt1 + t1 * vec, pt1 + t2 * vec]
            valid_points = [pt for pt in points if
                            linesegment2d.point_belongs(pt)]
            return valid_points

        return []


class ArcEllipse2D(Edge):
    """

    """

    def __init__(self, start, interior, end, center, major_dir, name='',
                 extra=None):
        Edge.__init__(self, start, end, name)
        self.interior = interior
        self.center = center
        self.extra = extra
        self.major_dir = major_dir
        self.minor_dir = self.major_dir.deterministic_unit_normal_vector()

        frame = volmdlr.Frame2D(self.center, self.major_dir, self.minor_dir)
        start_new, end_new = frame.new_coordinates(
            self.start), frame.new_coordinates(self.end)
        interior_new, center_new = frame.new_coordinates(
            self.interior), frame.new_coordinates(self.center)

        def theta_A_B(s, i, e, c):
            """
            from : https://math.stackexchange.com/questions/339126/how-to-draw-an-ellipse-if-a-center-and-3-arbitrary-points-on-it-are-given
            theta=angle d'inclinaison ellipse par rapport à horizontal(sens horaire),A=demi grd axe, B=demi petit axe
            """
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
            extra_new = frame.new_coordinates(self.extra)
            theta, A, B = theta_A_B(start_new, extra_new, interior_new,
                                    center_new)
        else:
            theta, A, B = theta_A_B(start_new, interior_new, end_new,
                                    center_new)

        self.Gradius = A
        self.Sradius = B
        self.theta = theta

        # Angle pour start
        u1, u2 = start_new.x / self.Gradius, start_new.y / self.Sradius
        angle1 = volmdlr.core.sin_cos_angle(u1, u2)
        # Angle pour end
        u3, u4 = end_new.x / self.Gradius, end_new.y / self.Sradius
        angle2 = volmdlr.core.sin_cos_angle(u3, u4)
        # Angle pour interior
        u5, u6 = interior_new.x / self.Gradius, interior_new.y / self.Sradius
        anglei = volmdlr.core.sin_cos_angle(u5, u6)

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

        if self.start == self.end or self.angle == 0:
            self.angle = volmdlr.TWO_PI

        if self.is_trigo:  # sens trigo
            self.offset_angle = angle1
        else:
            self.offset_angle = angle2

    def _get_points(self):
        return self.discretization_points()

    points = property(_get_points)

    def discretization_points(self, *, number_points: int = None, angle_resolution: int = None):
        """
        discretize a Edge to have "n" points
        :param number_points: the number of points (including start and end points)
             if unset, only start and end will be returned
        :param angle_resolution: if set, the sampling will be adapted to have a controlled angular distance. Usefull
            to mesh an arc
        :return: a list of sampled points
        """
        if not number_points:
            if not angle_resolution:
                number_points = 2
            else:
                number_points = math.ceil(angle_resolution * abs(0.5 * self.angle / math.pi))

        step = self.length() / (number_points - 1)
        return [self.point_at_abscissa(i * step)
                for i in range(number_points)]

    def polygon_points(self, discretization_resolution: int):
        warnings.warn('polygon_points is deprecated,\
                please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(discretization_resolution)

    def to_3d(self, plane_origin, x, y):
        ps = self.start.to_3d(plane_origin, x, y)
        pi = self.interior.to_3d(plane_origin, x, y)
        pe = self.end.to_3d(plane_origin, x, y)
        pc = self.center.to_3d(plane_origin, x, y)

        major_dir = self.major_dir.to_3d(plane_origin, x, y)
        major_dir.normalize()

        return ArcEllipse3D(ps, pi, pe, pc, major_dir, name=self.name)

    def plot(self, ax=None, color='k', alpha=1):
        if ax is None:
            _, ax = plt.subplots()

        self.interior.plot(ax=ax, color='m')
        self.start.plot(ax=ax, color='r')
        self.end.plot(ax=ax, color='b')
        self.center.plot(ax=ax, color='y')

        x = []
        y = []
        for px, py in self.discretization_points():
            x.append(px)
            y.append(py)

        plt.plot(x, y, color=color, alpha=alpha)
        return ax

    def normal_vector(self, abscissa):
        raise NotImplementedError

    def unit_normal_vector(self, abscissa):
        raise NotImplementedError

    def direction_vector(self, abscissa):
        raise NotImplementedError

    def unit_direction_vector(self, abscissa):
        raise NotImplementedError


class Line3D(Line):
    _non_eq_attributes = ['name', 'basis_primitives', 'bounding_box']

    """
    Define an infinite line passing through the 2 points
    """

    def __init__(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D,
                 name: str = ''):
        Line.__init__(self, point1, point2, name=name)
        self.bounding_box = self._bounding_box()

    def _bounding_box(self):
        # points = [self.point1, self.point2]
        # xmin = min([pt[0] for pt in points])
        # xmax = max([pt[0] for pt in points])
        # ymin = min([pt[1] for pt in points])
        # ymax = max([pt[1] for pt in points])
        # zmin = min([pt[2] for pt in points])
        # zmax = max([pt[2] for pt in points])

        xmin = min([self.point1[0], self.point2[0]])
        xmax = max([self.point1[0], self.point2[0]])
        ymin = min([self.point1[1], self.point2[1]])
        ymax = max([self.point1[1], self.point2[1]])
        zmin = min([self.point1[2], self.point2[2]])
        zmax = max([self.point1[2], self.point2[2]])

        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def point_at_abscissa(self, abscissa):
        return self.point1 + (
            self.point2 - self.point1) * abscissa

    def point_belongs(self, point3d):
        if point3d == self.point1:
            return True
        return self.direction_vector().is_colinear_to(point3d - self.point1)

    def point_distance(self, point):
        vector1 = point - self.start
        vector1.to_vector()
        vector2 = self.end - self.start
        vector2.to_vector()
        return vector1.cross(vector2).norm() / vector2.norm()

    def plot(self, ax=None, color='k', alpha=1, dashed=True):
        if ax is None:
            ax = Axes3D(plt.figure())

        # Line segment
        x = [self.point1.x, self.point2.x]
        y = [self.point1.y, self.point2.y]
        z = [self.point1.z, self.point2.z]
        ax.plot(x, y, z, color=color, alpha=alpha)

        # Drawing 3 times length of segment on each side
        u = self.point2 - self.point1
        v1 = (self.point1 - 3 * u)
        x1, y1, z1 = v1.x, v1.y, v1.z
        v2 = (self.point2 - 3 * u)
        x2, y2, z2 = v2.x, v2.y, v2.z
        if dashed:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=color,
                    dashes=[30, 5, 10, 5])
        else:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=color)
        return ax

    def plane_projection2d(self, center, x, y):
        return Line2D(self.points[0].plane_projection2d(center, x, y),
                      self.point2.plane_projection2d(center, x, y))

    def minimum_distance_points(self, other_line):
        """
        Returns the points on this line and the other line that are the closest
        of lines
        """
        u = self.point2 - self.point1
        v = other_line.point2 - other_line.point1
        w = self.point1 - other_line.point1
        a = u.dot(u)
        b = u.dot(v)
        c = v.dot(v)
        d = u.dot(w)
        e = v.dot(w)

        s = (b * e - c * d) / (a * c - b ** 2)
        t = (a * e - b * d) / (a * c - b ** 2)
        p1 = self.point1 + s * u
        p2 = other_line.point1 + t * v
        return p1, p2

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Line3D rotation
        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Line3D
        """

        return Line3D(*[p.rotation(center, axis, angle) for p in
                        [self.point1, self.point2]])

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Line3D rotation. Object is updated inplace
        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        for p in [self.point1, self.point2]:
            p.rotation_inplace(center, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Line3D translation
        :param offset: translation vector
        :return: A new translated Line3D
        """
        return Line3D(*[point.translation(offset) for point in
                        [self.point1, self.point2]])

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Line3D translation. Object is updated inplace
        :param offset: translation vector
        """
        for point in [self.point1, self.point2]:
            point.translation_inplace(offset)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and return a new Line3D
        side = 'old' or 'new'
        """
        if side == 'old':
            new_start = frame.old_coordinates(self.point1)
            new_end = frame.old_coordinates(self.point2)
        elif side == 'new':
            new_start = frame.new_coordinates(self.point1)
            new_end = frame.new_coordinates(self.point2)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        return Line3D(new_start, new_end)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes Line3D frame_mapping and the object is updated inplace
        side = 'old' or 'new'
        """
        if side == 'old':
            new_start = frame.old_coordinates(self.point1)
            new_end = frame.old_coordinates(self.point2)
        elif side == 'new':
            new_start = frame.new_coordinates(self.point1)
            new_end = frame.new_coordinates(self.point2)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        self.point1 = new_start
        self.point2 = new_end
        self.bounding_box = self._bounding_box()

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D):
        if not self.point_belongs(point1) or not self.point_belongs(point2):
            raise ValueError('Point not on curve')
        return Line3D(point1, point2)

    def copy(self, *args, **kwargs):
        return Line3D(*[p.copy() for p in self.points])

    @classmethod
    def from_step(cls, arguments, object_dict):
        point1 = object_dict[arguments[1]]
        direction = object_dict[arguments[2]]
        point2 = point1 + direction
        return cls(point1, point2, arguments[0][1:-1])

    def intersection(self, line2):

        x1 = self.point1.x
        y1 = self.point1.y
        z1 = self.point1.z
        x2 = self.point2.x
        y2 = self.point2.y
        z2 = self.point2.z
        x3 = line2.point1.x
        y3 = line2.point1.y
        z3 = line2.point1.z
        x4 = line2.point2.x
        y4 = line2.point2.y
        z4 = line2.point2.z

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
                return volmdlr.Point3D(x1 + (x2 - x1) * t1,
                                       y1 + (y2 - y1) * t1,
                                       z1 + (z2 - z1) * t1)

        return None

    def to_step(self, current_id):
        p1_content, p1_id = self.point1.to_step(current_id)
        # p2_content, p2_id = self.point2.to_step(current_id+1)
        current_id = p1_id + 1
        u_content, u_id = volmdlr.Vector3D.to_step(
            self.unit_direction_vector(),
            current_id,
            vector=True)
        current_id = u_id + 1
        content = p1_content + u_content
        content += f"#{current_id} = LINE('{self.name}',#{p1_id},#{u_id});\n"
        return content, current_id


class LineSegment3D(LineSegment):
    """
    Define a line segment limited by two points
    """

    def __init__(self, start: volmdlr.Point3D, end: volmdlr.Point3D,
                 name: str = ''):
        if start == end:
            raise NotImplementedError
        self.points = [start, end]
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
        return 2 + hash(self.start) + hash(self.end)

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

    # def point_at_abscissa(self, abscissa):
    #     return self.start + abscissa * (
    #         self.end - self.start) / self.length()

    def point_belongs(self, point, abs_tol=1e-7):
        distance = self.start.point_distance(point) + self.end.point_distance(
            point)
        if math.isclose(distance, self.length(), abs_tol=abs_tol):
            return True
        return False

    def normal_vector(self, abscissa=0.):
        return None

    def unit_normal_vector(self, abscissa=0.):
        return None

    # def middle_point(self):
    #     return self.point_at_abscissa(0.5 * self.length())

    def point_distance(self, point):
        distance, point = volmdlr.core_compiled.LineSegment3DPointDistance(
            [(self.start.x, self.start.y, self.start.z),
             (self.end.x, self.end.y, self.end.z)],
            (point.x, point.y, point.z))
        return distance, volmdlr.Point3D(*point)

    def plane_projection2d(self, center, x, y):
        return LineSegment2D(self.start.plane_projection2d(center, x, y),
                             self.end.plane_projection2d(center, x, y))

    def intersection(self, segment2):
        x1 = self.start.x
        y1 = self.start.y
        z1 = self.start.z
        x2 = self.end.x
        y2 = self.end.y
        z2 = self.end.z
        x3 = segment2.start.x
        y3 = segment2.start.y
        z3 = segment2.start.z
        x4 = segment2.end.x
        y4 = segment2.end.y
        z4 = segment2.end.z

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
                    return volmdlr.Point3D(x1 + (x2 - x1) * t1,
                                           y1 + (y2 - y1) * t1,
                                           z1 + (z2 - z1) * t1)

        return None

    def linesegment_intersection(self, linesegment):
        intersection = self.intersection(linesegment)
        if intersection in [self.start, self.end] and intersection in [linesegment.start, linesegment.end]:
            return intersection
        if intersection is not None:
            if self.point_belongs(intersection) and linesegment.point_belongs(intersection):
                return intersection
            return None
        return None

    def rotation(self, center: volmdlr.Point3D,
                 axis: volmdlr.Vector3D, angle: float):
        """
        LineSegment3D rotation
        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated LineSegment3D
        """
        return LineSegment3D(
            *[point.rotation(center, axis, angle) for point in self.points])

    def rotation_inplace(self, center: volmdlr.Point3D,
                         axis: volmdlr.Vector3D, angle: float):
        """
        Line2D rotation. Object is updated inplace
        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        for point in self.points:
            point.rotation_inplace(center, axis, angle)
        self.bounding_box = self._bounding_box()

    def __contains__(self, point):

        point1, point2 = self.start, self.end
        axis = point2 - point1
        test = point.rotation(point1, axis, math.pi)
        if test == point:
            return True

        return False

    def translation(self, offset: volmdlr.Vector3D):
        """
        LineSegment3D translation
        :param offset: translation vector
        :return: A new translated LineSegment3D
        """
        return LineSegment3D(
            *[point.translation(offset) for point in self.points])

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        LineSegment3D translation. Object is updated inplace
        :param offset: translation vector
        """
        for point in self.points:
            point.translation_inplace(offset)
        self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes LineSegment3D frame_mapping and return a new LineSegment3D
        side = 'old' or 'new'
        """
        if side == 'old':
            return LineSegment3D(
                *[frame.old_coordinates(point) for point in self.points])
        elif side == 'new':
            return LineSegment3D(
                *[frame.new_coordinates(point) for point in self.points])
        raise ValueError('Please Enter a valid side: old or new')

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and the object is updated inplace
        side = 'old' or 'new'
        """
        if side == 'old':
            new_start = frame.old_coordinates(self.start)
            new_end = frame.old_coordinates(self.end)
        elif side == 'new':
            new_start = frame.new_coordinates(self.start)
            new_end = frame.new_coordinates(self.end)
        else:
            raise ValueError('Please Enter a valid side: old or new')
        self.start = new_start
        self.end = new_end
        self.bounding_box = self._bounding_box()

    def copy(self, *args, **kwargs):
        return LineSegment3D(self.start.copy(), self.end.copy())

    def plot(self, ax=None, color='k', alpha=1,
             edge_ends=False, edge_direction=False):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        points = [self.start, self.end]
        x = [p.x for p in points]
        y = [p.y for p in points]
        z = [p.z for p in points]
        if edge_ends:
            ax.plot(x, y, z, color=color, alpha=alpha, marker='o')
        else:
            ax.plot(x, y, z, color=color, alpha=alpha)
        if edge_direction:
            x, y, z = self.point_at_abscissa(0.5 * self.length())
            u, v, w = 0.05 * self.direction_vector()
            ax.quiver(x, y, z, u, v, w, length=self.length() / 100,
                      arrow_length_ratio=5, normalize=True,
                      pivot='tip', color=color)
        return ax

    def plot2d(self, x_3D, y_3D, ax=None, color='k', width=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        edge2D = self.plane_projection2d(volmdlr.O3D, x_3D, y_3D)
        edge2D.plot(ax=ax, color=color, width=width)
        return ax

    def plot_data(self, x_3D, y_3D, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        edge2D = self.plane_projection2d(volmdlr.O3D, x_3D, y_3D)
        return edge2D.plot_data(marker, color, stroke_width,
                                dash, opacity, arrow)

    def FreeCADExport(self, name, ndigits=6):
        name = 'primitive' + str(name)
        x1, y1, z1 = round(1000 * self.start, ndigits)
        x2, y2, z2 = round(1000 * self.end, ndigits)
        return '{} = Part.LineSegment(fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(
            name, x1, y1, z1, x2, y2, z2)

    def to_line(self):
        return Line3D(self.start, self.end)

    def babylon_script(self, color=(1, 1, 1), name='line', type_='line',
                       parent=None):
        if type_ in ['line', 'dashed']:
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
                s += f'var {name} = BABYLON.MeshBuilder.CreateDashedLines("lines", {{points: myPoints, dashNb:20}}, scene);'
            s += '{}.color = new BABYLON.Color3{};\n'.format(name, tuple(color))
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

    def to_2d(self, plane_origin, x1, x2):
        p2d = [p.to_2d(plane_origin, x1, x2) for p in (self.start, self.end)]
        if p2d[0] == p2d[1]:
            return None
        return LineSegment2D(*p2d, name=self.name)

    def reverse(self):
        return LineSegment3D(self.end.copy(), self.start.copy())

    def minimum_distance_points(self, other_line):
        """
        Returns the points on this line and the other line that are the closest
        of lines
        """
        u = self.end - self.start
        v = other_line.end - other_line.start
        w = self.start - other_line.start
        a = u.dot(u)
        b = u.dot(v)
        c = v.dot(v)
        d = u.dot(w)
        e = v.dot(w)
        if (a * c - b ** 2) != 0:
            s = (b * e - c * d) / (a * c - b ** 2)
            t = (a * e - b * d) / (a * c - b ** 2)
            p1 = self.start + s * u
            p2 = other_line.start + t * v
            return p1, p2
        else:
            return self.start, other_line.start

    def Matrix_distance(self, other_line):
        u = self.direction_vector()
        v = other_line.direction_vector()
        w = other_line.start - self.start

        a = u.dot(u)
        b = -u.dot(v)
        d = v.dot(v)

        e = w.dot(u)
        f = -w.dot(v)

        A = npy.array([[a, b],
                       [b, d]])
        B = npy.array([e, f])

        res = scp.optimize.lsq_linear(A, B, bounds=(0, 1))
        p1 = self.point_at_abscissa(res.x[0] * self.length())
        p2 = other_line.point_at_abscissa(
            res.x[1] * other_line.length())
        return p1, p2

    def parallele_distance(self, other_linesegment):
        ptA, ptB, ptC = self.start, self.end, other_linesegment.points[0]
        u = volmdlr.Vector3D((ptA - ptB).vector)
        u.normalize()
        plane1 = volmdlr.faces.Plane3D.from_3_points(ptA, ptB, ptC)
        v = u.cross(plane1.normal)  # distance vector
        # ptA = k*u + c*v + ptC
        res = (ptA - ptC).vector
        x, y, z = res[0], res[1], res[2]
        u1, u2, u3 = u.x, u.y, u.z
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
        else:
            return NotImplementedError

    def minimum_distance(self, element, return_points=False):
        if element.__class__ is Arc3D or element.__class__ is volmdlr.wires.Circle3D:
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

        elif element.__class__ is BSplineCurve3D:
            points = element.points
            lines = []
            dist_min = math.inf
            for p1, p2 in zip(points[0:-1], points[1:]):
                lines.append(LineSegment3D(p1, p2))
            for line in lines:
                p1, p2 = self.Matrix_distance(line)
                dist = p1.point_distance(p2)
                if dist < dist_min:
                    dist_min = dist
                    min_points = (p1, p2)
            if return_points:
                p1, p2 = min_points
                return dist_min, p1, p2
            else:
                return dist_min

        else:
            return NotImplementedError

    def extrusion(self, extrusion_vector):
        u = self.unit_direction_vector()
        v = extrusion_vector.copy()
        v.normalize()
        w = u.cross(v)
        l1 = self.length()
        l2 = extrusion_vector.norm()
        # outer_contour = Polygon2D([O2D, Point2D((l1, 0.)),
        #                            Point2D((l1, l2)), Point2D((0., l2))])
        plane = volmdlr.faces.Plane3D(volmdlr.Frame3D(self.start, u, v, w))
        return [plane.rectangular_cut(0, l1, 0, l2)]

    def revolution(self, axis_point, axis, angle):
        axis_line3d = Line3D(axis_point, axis_point + axis)
        if axis_line3d.point_belongs(self.start) and axis_line3d.point_belongs(
                self.end):
            return []

        p1_proj, _ = axis_line3d.point_projection(self.start)
        p2_proj, _ = axis_line3d.point_projection(self.end)
        d1 = self.start.point_distance(p1_proj)
        d2 = self.end.point_distance(p2_proj)
        if not math.isclose(d1, 0., abs_tol=1e-9):
            u = (self.start - p1_proj)  # Unit vector from p1_proj to p1
            u.normalize()
        elif not math.isclose(d2, 0., abs_tol=1e-9):
            u = (self.end - p2_proj)  # Unit vector from p1_proj to p1
            u.normalize()
        else:
            return []
        if u.is_colinear_to(self.direction_vector()):
            # Planar face
            v = axis.cross(u)
            surface = volmdlr.faces.Plane3D(
                volmdlr.Frame3D(p1_proj, u, v, axis))
            r, R = sorted([d1, d2])
            if angle == volmdlr.TWO_PI:
                # Only 2 circles as countours
                outer_contour2d = volmdlr.wires.Circle2D(volmdlr.O2D, R)
                if not math.isclose(r, 0, abs_tol=1e-9):
                    inner_contours2d = [volmdlr.wires.Circle2D(volmdlr.O2D, r)]
                else:
                    inner_contours2d = []
            else:
                inner_contours2d = []
                if math.isclose(r, 0, abs_tol=1e-9):
                    # One arc and 2 lines (pizza slice)
                    arc2_e = volmdlr.Point2D(R, 0)
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
                    arc1_s = volmdlr.Point2D(R, 0)
                    arc1_i = arc1_s.rotation(center=volmdlr.O2D,
                                             angle=0.5 * angle)
                    arc1_e = arc1_s.rotation(center=volmdlr.O2D, angle=angle)
                    arc1 = Arc2D(arc1_s, arc1_i, arc1_e)

                    arc2_e = volmdlr.Point2D(r, 0)
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

        elif not math.isclose(d1, d2, abs_tol=1e-9):
            # Conical
            v = axis.cross(u)
            dv = self.direction_vector()
            dv.normalize()

            semi_angle = math.atan2(dv.dot(u), dv.dot(axis))
            cone_origin = p1_proj - d1 / math.tan(semi_angle) * axis
            if semi_angle > 0.5 * math.pi:
                semi_angle = math.pi - semi_angle

                cone_frame = volmdlr.Frame3D(cone_origin, u, -v, -axis)
                angle2 = -angle
            else:
                angle2 = angle
                cone_frame = volmdlr.Frame3D(cone_origin, u, v, axis)

            surface = volmdlr.faces.ConicalSurface3D(cone_frame,
                                                     semi_angle)
            z1 = d1 / math.tan(semi_angle)
            z2 = d2 / math.tan(semi_angle)
            return [surface.rectangular_cut(0, angle2, z1, z2)]
        else:
            # Cylindrical face
            v = axis.cross(u)
            surface = volmdlr.faces.CylindricalSurface3D(
                volmdlr.Frame3D(p1_proj, u, v, axis), d1)
            return [surface.rectangular_cut(0, angle,
                                            0,
                                            (self.end - self.start).dot(axis))]

    def to_step(self, current_id, surface_id=None):
        line = self.to_line()
        content, line_id = line.to_step(current_id)

        if surface_id:
            content += "#{} = SURFACE_CURVE('',#{},(#{}),.PCURVE_S1.);\n".format(
                line_id + 1, line_id, surface_id)
            line_id += 1

        current_id = line_id + 1
        start_content, start_id = self.start.to_step(current_id, vertex=True)
        current_id = start_id + 1
        end_content, end_id = self.end.to_step(current_id + 1, vertex=True)
        content += start_content + end_content
        current_id = end_id + 1
        content += "#{} = EDGE_CURVE('{}',#{},#{},#{},.T.);\n".format(
            current_id, self.name,
            start_id, end_id, line_id)
        return content, [current_id]


class BSplineCurve3D(BSplineCurve, volmdlr.core.Primitive3D):
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
        volmdlr.core.Primitive3D.__init__(self, name=name)

        self.bounding_box = self._bounding_box()

    def _bounding_box(self):
        bbox = self.curve.bbox
        return volmdlr.core.BoundingBox(bbox[0][0], bbox[1][0],
                                        bbox[0][1], bbox[1][1],
                                        bbox[0][2], bbox[1][2])

    def look_up_table(self, resolution: int = 20, start_parameter: float = 0,
                      end_parameter: float = 1):
        """
        Creates a table of equivalence between the parameter t (evaluation
        of the BSplineCurve) and the cumulative distance.

        :param resolution: The precision of the table. Autoadjusted by the
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

    def point_at_abscissa(self, abscissa: float, resolution: int = 1000):
        """
        Returns the 3 dimensional point at a given curvilinear abscissa.
        This is an approximation. Resolution parameter can be increased
        for more accurate result.

        :param abscissa: The distance on the BSplineCurve3D from its start
        :type abscissa: float
        :param resolution: The precision of the approximation. Default value
            set to 1000
        :type resolution: int, optional
        :return: The Point3D at the given curvilinear abscissa.
        :rtype: :class:`volmdlr.Point3D`
        """
        if math.isclose(abscissa, 0, abs_tol=1e-10):
            return self.start
        elif math.isclose(abscissa, self.length(), abs_tol=1e-10):
            return self.end
        lut = self.look_up_table(resolution=resolution)
        if 0 < abscissa < self.length():
            last_param = 0
            for i, (t, dist) in enumerate(lut):
                if abscissa < dist:
                    t1 = last_param
                    t2 = t
                    return volmdlr.Point3D(
                        *self.curve.evaluate_single((t1 + t2) / 2))
                last_param = t
        else:
            raise ValueError('Curvilinear abscissa is bigger than length,'
                             ' or negative')

    def normal(self, position: float = 0.0):
        point, normal = operations.normal(self.curve, position, normalize=True)
        normal = volmdlr.Point3D(normal[0], normal[1], normal[2])
        return normal

    def direction_vector(self, abscissa=0.):
        l = self.length()
        if abscissa >= l:
            abscissa2 = l
            abscissa = abscissa2 - 0.001 * l

        else:
            abscissa2 = min(abscissa + 0.001 * l, l)

        tangent = self.point_at_abscissa(abscissa2) - self.point_at_abscissa(
            abscissa)
        return tangent

    def normal_vector(self, abscissa):
        return None

    def unit_normal_vector(self, abscissa):
        return None

    def point3d_to_parameter(self, point: volmdlr.Point3D):
        """
        Search for the value of the normalized evaluation parameter t
        (between 0 and 1) that would return the given point when the
        BSplineCurve3D is evaluated at the t value.
        """
        def f(param):
            p3d = volmdlr.Point3D(*self.curve.evaluate_single(param))
            return point.point_distance(p3d)
        res = scipy.optimize.minimize(fun=f, x0=(0.5), bounds=[(0, 1)],
                                      tol=1e-9)
        return res.x[0]

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        points = '['
        for i in range(len(self.control_points)):
            point = 'fc.Vector({},{},{}),'.format(self.control_points[i][0],
                                                  self.control_points[i][1],
                                                  self.control_points[i][2])
            points += point
        points = points[:-1]
        points += ']'
        # !!! : A QUOI SERT LE DERNIER ARG DE BSplineCurve (False)?
        # LA MULTIPLICITE EN 3e ARG ET LES KNOTS EN 2e ARG ?
        return '{} = Part.BSplineCurve({},{},{},{},{},{},{})\n'.format(name,
                                                                       points,
                                                                       self.knot_multiplicities,
                                                                       self.knots,
                                                                       self.periodic,
                                                                       self.degree,
                                                                       self.weights,
                                                                       False)

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]
        degree = int(arguments[1])
        points = [object_dict[int(i[1:])] for i in arguments[2]]
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

    def to_step(self, current_id, surface_id=None):

        points_ids = []
        content = ''
        for point in self.points:
            point_content, point_id = point.to_step(current_id,
                                                    vertex=True)
            content += point_content
            points_ids.append(point_id)

        curve_id = point_id + 1
        content += "#{} = B_SPLINE_CURVE_WITH_KNOTS('{}',{},({})," \
                   ".UNSPECIFIED.,.F.,.F.,({}),{}," \
                   ".PIECEWISE_BEZIER_KNOTS.);\n".format(curve_id,
                                                         self.name,
                                                         self.degree,
                                                         volmdlr.core.step_ids_to_str(
                                                             points_ids),
                                                         volmdlr.core.step_ids_to_str(
                                                             self.knot_multiplicities),
                                                         tuple(self.knots)
                                                         )

        if surface_id:
            content += "#{} = SURFACE_CURVE('',#{},(#{}),.PCURVE_S1.);\n".format(
                curve_id + 1, curve_id, surface_id)
            curve_id += 1

        current_id = curve_id + 1
        start_content, start_id = self.start.to_step(current_id, vertex=True)
        current_id = start_id + 1
        end_content, end_id = self.end.to_step(current_id + 1, vertex=True)
        content += start_content + end_content
        current_id = end_id + 1
        content += "#{} = EDGE_CURVE('{}',#{},#{},#{},.T.);\n".format(
            current_id, self.name,
            start_id, end_id, curve_id)
        return content, [current_id]

    def point_distance(self, pt1):
        distances = []
        for point in self.points:
            #            vmpt = Point3D((point[1], point[2], point[3]))
            distances.append(pt1.point_distance(point))
        return min(distances)

    # def point_belongs(self, point):
    #     polygon_points = self.polygon_points()
    #     for p1, p2 in zip(polygon_points[:-1], polygon_points[1:]):
    #         line = LineSegment3D(p1, p2)
    #         if line.point_belongs(point):
    #             return True
    #     return False

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        BSplineCurve3D rotation
        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated BSplineCurve3D
        """
        new_control_points = [p.rotation(center, axis, angle) for p in
                              self.control_points]
        new_bsplinecurve3d = BSplineCurve3D(self.degree, new_control_points,
                                            self.knot_multiplicities,
                                            self.knots, self.weights,
                                            self.periodic, self.name)
        return new_bsplinecurve3d

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        BSplineCurve3D rotation. Object is updated inplace
        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        new_control_points = [p.rotation(center, axis, angle) for p in
                              self.control_points]
        new_bsplinecurve3d = BSplineCurve3D(self.degree, new_control_points,
                                            self.knot_multiplicities,
                                            self.knots, self.weights,
                                            self.periodic, self.name)

        self.control_points = new_control_points
        self.curve = new_bsplinecurve3d.curve
        self.points = new_bsplinecurve3d.points

    def trim(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D):
        if (point1 == self.start and point2 == self.end) \
                or (point1 == self.end and point2 == self.start):
            return self

        elif point1 == self.start and point2 != self.end:
            parameter2 = self.point3d_to_parameter(point2)
            return self.cut_after(parameter2)

        elif point2 == self.start and point1 != self.end:
            parameter1 = self.point3d_to_parameter(point1)
            return self.cut_after(parameter1)

        elif point1 != self.start and point2 == self.end:
            parameter1 = self.point3d_to_parameter(point1)
            return self.cut_before(parameter1)

        elif point2 != self.start and point1 == self.end:
            parameter2 = self.point3d_to_parameter(point2)
            return self.cut_before(parameter2)

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
        elif math.isclose(parameter1, 0, abs_tol=1e-7):
            return self.cut_after(parameter2)
        elif math.isclose(parameter2, 1, abs_tol=1e-7):
            return self.cut_before(parameter1)

        # Cut before
        bspline_curve = self.insert_knot(parameter1, num=self.degree)
        if bspline_curve.weights is not None:
            raise NotImplementedError

        # Cut after
        bspline_curve = bspline_curve.insert_knot(parameter2, num=self.degree)
        if bspline_curve.weights is not None:
            raise NotImplementedError

        # Que faire quand on rajoute un noeud au milieu ?
        # plus simple de passer par cut_after cut_before
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
        # Is a value of parameter below 4e-3 a real need for precision ?
        if math.isclose(parameter, 0, abs_tol=4e-3):
            return self
        elif math.isclose(parameter, 1, abs_tol=4e-3):
            raise ValueError('Nothing will be left from the BSplineCurve3D')
        curves = operations.split_curve(self.curve, parameter)
        return self.from_geomdl_curve(curves[1])

    def cut_after(self, parameter: float):
        # Is a value of parameter below 4e-3 a real need for precision ?
        if math.isclose(parameter, 0, abs_tol=4e-3):
            raise ValueError('Nothing will be left from the BSplineCurve3D')
        if math.isclose(parameter, 1, abs_tol=4e-3):
            return self
        curves = operations.split_curve(self.curve, parameter)
        return self.from_geomdl_curve(curves[0])

    def insert_knot(self, knot: float, num: int = 1):
        """
        Returns a new BSplineCurve3D
        """
        curve_copy = self.curve.__deepcopy__({})
        modified_curve = operations.insert_knot(curve_copy, [knot], num=[num])
        return self.from_geomdl_curve(modified_curve)

    # Copy paste du LineSegment3D
    def plot(self, ax=None, edge_ends=False, color='k', alpha=1,
             edge_direction=False):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        x = [p.x for p in self.points]
        y = [p.y for p in self.points]
        z = [p.z for p in self.points]
        ax.plot(x, y, z, color=color, alpha=alpha)
        if edge_ends:
            ax.plot(x, y, z, 'o', color=color, alpha=alpha)
        return ax

    def to_2d(self, plane_origin, x1, x2):
        control_points2d = [p.to_2d(plane_origin, x1, x2) for p in
                            self.control_points]
        return BSplineCurve2D(self.degree, control_points2d,
                              self.knot_multiplicities, self.knots,
                              self.weights, self.periodic, self.name)

    def polygon_points(self):
        warnings.warn('polygon_points is deprecated,\
        please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points()

    def curvature(self, u: float, point_in_curve: bool = False):
        # u should be in the interval [0,1]
        curve = self.curve
        ders = curve.derivatives(u, 3)  # 3 first derivative
        c1, c2 = volmdlr.Point3D(*ders[1]), volmdlr.Point3D(*ders[2])
        denom = c1.cross(c2)
        if c1 == volmdlr.O3D or c2 == volmdlr.O3D or denom.norm() == 0.0:
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
        Returns the maximum curvature of a curve and the point where it is located
        """
        if point_in_curve:
            maximum_curvarture, point = max(self.global_maximum_curvature(nb_eval=21, point_in_curve=point_in_curve))
            return maximum_curvarture, point
        # print(self.global_maximum_curvature(point_in_curve))
        maximum_curvarture = max(self.global_maximum_curvature(nb_eval=21, point_in_curve=point_in_curve))
        return maximum_curvarture

    def minimum_radius(self, point_in_curve=False):
        """
        Returns the minimum curvature radius of a curve and the point where it is located
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

    def linesegment_intersection(self, linesegment: LineSegment3D):
        points = self.discretization_points()
        linesegments = [LineSegment3D(start, end) for start, end in zip(points[:-1], points[1:])]
        for line_segment in linesegments:
            intersection = line_segment.linesegment_intersection(linesegment)
            if intersection:
                return intersection
        return None


class BezierCurve3D(BSplineCurve3D):

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
    An arc is defined by a starting point, an end point and an interior point

    """

    def __init__(self, start, interior, end, name=''):
        """

        """
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
        # self.bounding_box = self._bounding_box()

    @property
    def bounding_box(self):
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def get_bounding_box(self):
        # TODO: implement exact calculation

        points = self.polygon_points()
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
                'Start, end and interior points of an arc must be distincts')

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
        u1 = self.interior - self.start
        u2 = self.interior - self.end
        if u1 == u2:
            u2 = self.normal.cross(u1)
            u2.normalize()

        v1 = self.normal.cross(u1)  # v1 is normal, equal u2
        v2 = self.normal.cross(u2)  # equal -u1

        p11 = 0.5 * (self.start + self.interior)  # Mid point of segment s,m
        p12 = p11 + v1
        p21 = 0.5 * (self.end + self.interior)  # Mid point of segment s,m
        p22 = p21 + v2

        l1 = Line3D(p11, p12)
        l2 = Line3D(p21, p22)

        try:
            center, _ = l1.minimum_distance_points(l2)
        except ZeroDivisionError:
            raise ValueError(
                'Start, end and interior points  of an arc must be distincts')

        return center

    @property
    def frame(self):
        if not self._utd_frame:
            self._frame = self.get_frame()
            self._utd_frame = True
        return self._frame

    def get_frame(self):
        vec1 = (self.start - self.center)
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
        Verifies if arc is clockwise of trigowise
        :return:
        """
        clockwise_path, trigowise_path = self.clockwise_and_trigowise_paths
        if clockwise_path > trigowise_path:
            return True
        return False

    @property
    def clockwise_and_trigowise_paths(self):
        """
        :return: clockwise path and trigonomectric path property
        """
        if not self._utd_clockwise_and_trigowise_paths:
            vec1 = (self.start - self.center)
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
        Arc angle property
        :return: arc angle
        """
        if not self._utd_angle:
            self._angle = self.get_angle()
            self._utd_angle = True
        return self._angle

    def get_angle(self):
        """
        Gets the arc angle
        :return: arc angle
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
        return self.__class__(self.end.copy(),
                              self.interior.copy(),
                              self.start.copy())

    def point_at_abscissa(self, abscissa):
        return self.start.rotation(self.center, self.normal,
                                   abscissa / self.radius)

    def normal_vector(self, abscissa):
        theta = abscissa / self.radius
        n_0 = self.center - self.start
        normal = n_0.rotation(self.center, self.normal, theta)
        return normal

    def unit_normal_vector(self, abscissa):
        normal_vector = self.normal_vector(abscissa)
        normal_vector.normalize()
        return normal_vector

    def direction_vector(self, abscissa):
        normal_vector = self.normal_vector(abscissa)
        tangent = normal_vector.cross(self.normal)
        return tangent

    def unit_direction_vector(self, abscissa):
        direction_vector = self.direction_vector(abscissa)
        direction_vector.normalize()
        return direction_vector

    def rotation(self, center: volmdlr.Point3D,
                 axis: volmdlr.Vector3D, angle: float):
        """
        Arc3D rotation
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
        Arc3D rotation. Object is updated inplace
        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        self.center.rotation_inplace(center, axis, angle)
        self.start.rotation_inplace(center, axis, angle)
        self.interior.rotation_inplace(center, axis, angle)
        self.end.rotation_inplace(center, axis, angle)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box
        [p.rotation_inplace(center, axis, angle) for p in self.primitives]

    def translation(self, offset: volmdlr.Vector3D):
        """
        Arc3D translation
        :param offset: translation vector
        :return: A new translated Arc3D
        """
        new_start = self.start.translation(offset)
        new_interior = self.interior.translation(offset)
        new_end = self.end.translation(offset)
        return Arc3D(new_start, new_interior, new_end, name=self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Arc3D translation. Object is updated inplace
        :param offset: translation vector
        """
        self.center.translation_inplace(offset)
        self.start.translation_inplace(offset)
        self.interior.translation_inplace(offset)
        self.end.translation_inplace(offset)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box
        [p.translation_inplace(offset) for p in self.primitives]

    def plot(self, ax=None, color='k', alpha=1,
             edge_ends=False, edge_direction=False):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None
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
        for px, py, pz in self.discretization_points():
            x.append(px)
            y.append(py)
            z.append(pz)

        ax.plot(x, y, z, color=color, alpha=alpha)
        if edge_ends:
            self.start.plot(ax=ax)
            self.end.plot(ax=ax)

        if edge_direction:
            x, y, z = self.point_at_abscissa(0.5 * self.length())
            u, v, w = 0.05 * self.unit_direction_vector(0.5 * self.length())
            ax.quiver(x, y, z, u, v, w, length=self.length() / 100,
                      arrow_length_ratio=5, normalize=True,
                      pivot='tip', color=color)
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
        l = self.length()
        x = []
        y = []
        for i in range(30):
            p = self.point_at_abscissa(i / (29.) * l)
            xi, yi = p.plane_projection2d(center, x3d, y3d)
            x.append(xi)
            y.append(yi)
        ax.plot(x, y, color=color)

        return ax

    def FreeCADExport(self, name, ndigits=6):
        xs, ys, zs = round(1000 * self.start, ndigits)
        xi, yi, zi = round(1000 * self.interior, ndigits)
        xe, ye, ze = round(1000 * self.end, ndigits)
        return '{} = Part.Arc(fc.Vector({},{},{}),fc.Vector({},{},{}),fc.Vector({},{},{}))\n' \
            .format(name, xs, ys, zs, xi, yi, zi, xe, ye, ze)

    def copy(self, *args, **kwargs):
        return Arc3D(self.start.copy(), self.interior.copy(), self.end.copy())

    def frame_mapping_parameters(self, frame: volmdlr.Frame3D, side: str):
        if side == 'old':
            new_start = frame.old_coordinates(self.start.copy())
            new_interior = frame.old_coordinates(self.interior.copy())
            new_end = frame.old_coordinates(self.end.copy())
        elif side == 'new':
            new_start = frame.new_coordinates(self.start.copy())
            new_interior = frame.new_coordinates(self.interior.copy())
            new_end = frame.new_coordinates(self.end.copy())
        else:
            raise ValueError('side value not valid, please specify'
                             'a correct value: \'old\' or \'new\'')
        return new_start, new_interior, new_end

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and return a new Arc3D
        side = 'old' or 'new'
        """
        new_start, new_interior, new_end =\
            self.frame_mapping_parameters(frame, side)

        return Arc3D(new_start, new_interior, new_end, name=self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes vector frame_mapping and the object is updated inplace
        side = 'old' or 'new'
        """
        new_start, new_interior, new_end = \
            self.frame_mapping_parameters(frame, side)
        self.start, self.interior, self.end = new_start, new_interior, new_end

    def abscissa(self, point3d: volmdlr.Point3D):
        x, y, z = self.frame.new_coordinates(point3d)
        u1 = x / self.radius
        u2 = y / self.radius
        theta = volmdlr.core.sin_cos_angle(u1, u2)

        return self.radius * abs(theta)

    def split(self, split_point: volmdlr.Point3D):
        abscissa = self.abscissa(split_point)

        return [Arc3D(self.start,
                      self.point_at_abscissa(0.5 * abscissa),
                      split_point),
                Arc3D(split_point,
                      self.point_at_abscissa(1.5 * abscissa),
                      self.end)
                ]

    def to_2d(self, plane_origin, x, y):
        ps = self.start.to_2d(plane_origin, x, y)
        pi = self.interior.to_2d(plane_origin, x, y)
        pe = self.end.to_2d(plane_origin, x, y)
        return Arc2D(ps, pi, pe, name=self.name)

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

        res1 = scp.optimize.least_squares(distance_squared, x01,
                                          bounds=[(0, 0), (
                                              self.angle, other_arc.angle)])

        p1 = self.point_at_abscissa(res1.x[0] * r1)
        p2 = other_arc.point_at_abscissa(res1.x[1] * r2)

        return p1, p2

    def minimum_distance_points_line(self, other_line):

        u = other_line.direction_vector()
        k = self.start - self.center
        k.normalize()
        w = self.center - other_line.start
        v = self.normal.cross(k)

        r = self.radius

        a = u.dot(u)
        b = u.dot(v)
        c = u.dot(k)
        d = v.dot(v)
        e = v.dot(k)
        f = k.dot(k)
        g = w.dot(u)
        h = w.dot(v)
        i = w.dot(k)
        j = w.dot(w)

        # x = (s, theta)
        def distance_squared(x):
            return (a * x[0] ** 2 + j + d * (
                    (math.sin(x[1])) ** 2) * r ** 2 + f * (
                (math.cos(x[1])) ** 2) * r ** 2
                - 2 * x[0] * g - 2 * x[0] * r * math.sin(x[1]) * b - 2 * x[
                        0] * r * math.cos(x[1]) * c
                + 2 * r * math.sin(x[1]) * h + 2 * r * math.cos(x[1]) * i
                + math.sin(2 * x[1]) * e * r ** 2)

        x01 = npy.array([0.5, self.angle / 2])
        x02 = npy.array([0.5, 0])
        x03 = npy.array([0.5, self.angle])

        res1 = scp.optimize.least_squares(distance_squared, x01,
                                          bounds=[(0, 0), (1, self.angle)])
        res2 = scp.optimize.least_squares(distance_squared, x02,
                                          bounds=[(0, 0), (1, self.angle)])
        res3 = scp.optimize.least_squares(distance_squared, x03,
                                          bounds=[(0, 0), (1, self.angle)])

        p1 = other_line.point_at_abscissa(
            res1.x[0] * other_line.length())
        p2 = self.point_at_abscissa(res1.x[1] * r)

        res = [res2, res3]
        for couple in res:
            ptest1 = other_line.point_at_abscissa(
                couple.x[0] * other_line.length())
            ptest2 = self.point_at_abscissa(couple.x[1] * r)
            dtest = ptest1.point_distance(ptest2)
            if dtest < d:
                p1, p2 = ptest1, ptest2

        return p1, p2

    def minimum_distance(self, element, return_points=False):
        if element.__class__ is Arc3D or element.__class__.__name__ == 'Circle3D':
            p1, p2 = self.minimum_distance_points_arc(element)
            if return_points:
                return p1.point_distance(p2), p1, p2
            else:
                return p1.point_distance(p2)

        elif element.__class__ is LineSegment3D:
            pt1, pt2 = self.minimum_distance_points_line(element)
            if return_points:
                return pt1.point_distance(pt2), pt1, pt2
            else:
                return pt1.point_distance(pt2)
        else:
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
            return [cylinder.rectangular_cut(angle1,
                                             angle2,
                                             0, extrusion_vector.norm())]
        else:
            raise NotImplementedError(
                'Elliptic faces not handled: dot={}'.format(
                    self.normal.dot(extrusion_vector)
                ))

    def revolution(self, axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                   angle: float):
        line3d = Line3D(axis_point, axis_point + axis)
        tore_center, _ = line3d.point_projection(self.center)
        if math.isclose(tore_center.point_distance(self.center), 0.,
                        abs_tol=1e-9):
            # Sphere
            start_p, _ = line3d.point_projection(self.start)
            u = self.start - start_p

            if math.isclose(u.norm(), 0, abs_tol=1e-9):
                end_p, _ = line3d.point_projection(self.end)
                u = self.end - end_p
                if math.isclose(u.norm(), 0, abs_tol=1e-9):
                    interior_p, _ = line3d.point_projection(self.interior)
                    u = self.interior - interior_p

            u.normalize()
            v = axis.cross(u)
            arc2d = self.to_2d(self.center, u, axis)

            surface = volmdlr.faces.SphericalSurface3D(
                volmdlr.Frame3D(self.center, u, v, axis), self.radius)

            return [surface.rectangular_cut(0, angle,
                                            arc2d.angle1, arc2d.angle2)]

        else:
            # Toroidal
            u = self.center - tore_center
            u.normalize()
            v = axis.cross(u)
            if not math.isclose(self.normal.dot(u), 0., abs_tol=1e-9):
                raise NotImplementedError(
                    'Outside of plane revolution not supported')

            R = tore_center.point_distance(self.center)
            surface = volmdlr.faces.ToroidalSurface3D(
                volmdlr.Frame3D(tore_center, u, v, axis), R,
                self.radius)
            arc2d = self.to_2d(tore_center, u, axis)
            return [surface.rectangular_cut(0, angle,
                                            arc2d.angle1, arc2d.angle2)]

    def to_step(self, current_id):
        if self.angle >= math.pi:
            l = self.length()
            arc1, arc2 = self.split(self.point_at_abscissa(0.33 * l))
            arc2, arc3 = arc2.split(self.point_at_abscissa(0.66 * l))
            content, arcs1_id = arc1.to_step_without_splitting(current_id)
            arc2_content, arcs2_id = arc2.to_step_without_splitting(
                arcs1_id[0] + 1)
            arc3_content, arcs3_id = arc3.to_step_without_splitting(
                arcs2_id[0] + 1)
            content += arc2_content + arc3_content
            return content, [arcs1_id[0], arcs2_id[0], arcs3_id[0]]
        else:
            return self.to_step_without_splitting(current_id)

    def to_step_without_splitting(self, current_id, surface_id=None):
        u = self.start - self.center
        u.normalize()
        v = self.normal.cross(u)
        frame = volmdlr.Frame3D(self.center, self.normal, u, v)

        content, frame_id = frame.to_step(current_id)
        curve_id = frame_id + 1
        content += "#{} = CIRCLE('{}', #{}, {:.6f});\n".format(curve_id,
                                                               self.name,
                                                               frame_id,
                                                               self.radius * 1000,
                                                               )

        if surface_id:
            content += "#{} = SURFACE_CURVE('',#{},(#{}),.PCURVE_S1.);\n".format(
                curve_id + 1, curve_id, surface_id)
            curve_id += 1

        current_id = curve_id + 1
        start_content, start_id = self.start.to_step(current_id, vertex=True)
        end_content, end_id = self.end.to_step(start_id + 1, vertex=True)
        content += start_content + end_content
        current_id = end_id + 1
        content += "#{} = EDGE_CURVE('{}',#{},#{},#{},.T.);\n".format(
            current_id, self.name,
            start_id, end_id, curve_id)
        return content, [current_id]

    def point_belongs(self, point3d, abs_tol=1e-10):
        '''
        check if a point3d belongs to the arc_3d or not
        '''
        def f(x):
            return (point3d - self.point_at_abscissa(x)).norm()
        length_ = self.length()
        x = npy.linspace(0, length_, 5)
        x_init = []
        for xi in x:
            x_init.append(xi)

        for x0 in x_init:
            z = scp.optimize.least_squares(f, x0=x0, bounds=([0, length_]))
            if z.fun < abs_tol:
                return True
        return False

    def triangulation(self):
        return None

    def middle_point(self):
        return self.point_at_abscissa(self.length() / 2)


class FullArc3D(Arc3D):
    """
    An edge that starts at start_end, ends at the same point after having described
    a circle
    """

    def __init__(self, center: volmdlr.Point3D, start_end: volmdlr.Point3D,
                 normal: volmdlr.Vector3D,
                 name: str = ''):
        self.__center = center
        self.__normal = normal
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

    def to_2d(self, plane_origin, x1, x2):
        center = self.center.to_2d(plane_origin, x1, x2)
        start_end = self.start.to_2d(plane_origin, x1, x2)
        return FullArc2D(center, start_end)

    def to_step(self, current_id, surface_id=None):
        # Not calling Circle3D.to_step because of circular imports
        u = self.start - self.center
        u.normalize()
        v = self.normal.cross(u)
        frame = volmdlr.Frame3D(self.center, self.normal, u, v)
        content, frame_id = frame.to_step(current_id)
        curve_id = frame_id + 1
        # Not calling Circle3D.to_step because of circular imports
        content += "#{} = CIRCLE('{}',#{},{:.6f});\n".format(curve_id,
                                                             self.name,
                                                             frame_id,
                                                             self.radius * 1000,
                                                             )

        if surface_id:
            content += "#{} = SURFACE_CURVE('',#{},(#{}),.PCURVE_S1.);\n".format(
                curve_id + 1, curve_id, surface_id)
            curve_id += 1

        p1 = (self.center + u * self.radius).to_point()
        # p2 = self.center + v*self.radius
        # p3 = self.center - u*self.radius
        # p4 = self.center - v*self.radius

        p1_content, p1_id = p1.to_step(curve_id + 1, vertex=True)
        content += p1_content
        # p2_content, p2_id = p2.to_step(p1_id+1, vertex=True)
        # p3_content, p3_id = p3.to_step(p2_id+1, vertex=True)
        # p4_content, p4_id = p4.to_step(p3_id+1, vertex=True)
        # content += p1_content + p2_content + p3_content + p4_content

        # arc1_id = p4_id + 1
        # content += "#{} = EDGE_CURVE('{}',#{},#{},#{},.T.);\n".format(arc1_id, self.name,
        #                                                             p1_id, p2_id,
        #                                                             circle_id)

        # arc2_id = arc1_id + 1
        # content += "#{} = EDGE_CURVE('{}',#{},#{},#{},.T.);\n".format(arc2_id, self.name,
        #                                                             p2_id, p3_id,
        #                                                             circle_id)

        # arc3_id = arc2_id + 1
        # content += "#{} = EDGE_CURVE('{}',#{},#{},#{},.T.);\n".format(arc3_id, self.name,
        #                                                             p3_id, p4_id,
        #                                                             circle_id)

        # arc4_id = arc3_id + 1
        # content += "#{} = EDGE_CURVE('{}',#{},#{},#{},.T.);\n".format(arc4_id, self.name,
        #                                                             p4_id, p1_id,
        #                                                             circle_id)

        edge_curve = p1_id + 1
        content += f"#{edge_curve} = EDGE_CURVE('{self.name}',#{p1_id},#{p1_id},#{curve_id},.T.);\n"
        curve_id += 1

        # return content, [arc1_id, arc2_id, arc3_id, arc4_id]
        return content, [edge_curve]

    def plot(self, ax=None, color='k', alpha=1., edge_ends=False,
             edge_direction=False):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)

        x = []
        y = []
        z = []
        for px, py, pz in self.discretization_points():
            x.append(px)
            y.append(py)
            z.append(pz)
        x.append(x[0])
        y.append(y[0])
        z.append(z[0])
        ax.plot(x, y, z, color=color, alpha=alpha)

        if edge_ends:
            self.start.plot(ax=ax)
            self.end.plot(ax=ax)

        if edge_direction:
            s = 0.5 * self.length()
            x, y, z = self.point_at_abscissa(s)
            tangent = self.unit_direction_vector(s)
            arrow_length = 0.15 * s
            ax.quiver(x, y, z, *arrow_length * tangent,
                      pivot='tip')

        return ax

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        new_start_end = self.start.rotation(center, axis, angle, True)
        new_center = self._center.rotation(center, axis, angle, True)
        new_normal = self._normal.rotation(center, axis, angle, True)
        return FullArc3D(new_center, new_start_end,
                         new_normal, name=self.name)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        self.start.rotation(center, axis, angle, False)
        self.end.rotation(center, axis, angle, False)
        self._center.rotation(center, axis, angle, False)
        self.interior.rotation(center, axis, angle, False)

    def translation(self, offset: volmdlr.Vector3D):
        new_start_end = self.start.translation(offset, True)
        new_center = self._center.translation(offset, True)
        new_normal = self._normal.translation(offset, True)
        return FullArc3D(new_center, new_start_end,
                         new_normal, name=self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        self.start.translation(offset, False)
        self.end.translation(offset, False)
        self._center.translation(offset, False)
        self.interior.translation(offset, False)


class ArcEllipse3D(Edge):
    """
    An arc is defined by a starting point, an end point and an interior point
    """

    def __init__(self, start, interior, end, center, major_dir,
                 name=''):  # , extra=None):
        # Extra is an additionnal point if start=end because you need 3 points on the arcellipse to define it
        Edge.__init__(self, start=start, end=end, name=name)
        self.interior = interior
        self.center = center
        major_dir.normalize()
        self.major_dir = major_dir  # Vector for Gradius
        # self.extra = extra

        u1 = (self.interior - self.start)
        u2 = (self.interior - self.end)
        u1.normalize()
        u2.normalize()

        if u1 == u2:
            u2 = (self.interior - self.extra)
            u2.normalize()

        # if normal is None:
        n = u2.cross(u1)
        n.normalize()
        self.normal = n
        # else:
        #     n = normal
        #     n.normalize()
        #     self.normal = normal

        self.minor_dir = self.normal.cross(self.major_dir)

        frame = volmdlr.Frame3D(self.center, self.major_dir, self.minor_dir,
                                self.normal)
        start_new, end_new = frame.new_coordinates(
            self.start), frame.new_coordinates(self.end)
        interior_new, center_new = frame.new_coordinates(
            self.interior), frame.new_coordinates(self.center)

        # from :
        # https://math.stackexchange.com/questions/339126/how-to-draw-an-ellipse-if-a-center-and-3-arbitrary-points-on-it-are-given
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
            theta = 0.5 * math.atan(2 * C[2] / (C[1] - C[0]))
            c1 = C[0] + C[1]
            c2 = (C[1] - C[0]) / math.cos(2 * theta)
            gdaxe = math.sqrt((2 / (c1 - c2)))
            ptax = math.sqrt((2 / (c1 + c2)))
            return theta, gdaxe, ptax

        if start == end:
            extra_new = frame.new_coordinates(self.extra)
            theta, A, B = theta_A_B(start_new, extra_new, interior_new,
                                    center_new)
        else:
            theta, A, B = theta_A_B(start_new, interior_new, end_new,
                                    center_new)

        self.Gradius = A
        self.Sradius = B
        self.theta = theta

        # Angle pour start
        u1, u2 = start_new.x / self.Gradius, start_new.y / self.Sradius
        angle1 = volmdlr.core.sin_cos_angle(u1, u2)
        # Angle pour end
        u3, u4 = end_new.x / self.Gradius, end_new.y / self.Sradius
        angle2 = volmdlr.core.sin_cos_angle(u3, u4)
        # Angle pour interior
        u5, u6 = interior_new.x / self.Gradius, interior_new.y / self.Sradius
        anglei = volmdlr.core.sin_cos_angle(u5, u6)

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

        if self.start == self.end:
            self.angle = volmdlr.TWO_PI

        if self.is_trigo:
            self.offset_angle = angle1
        else:
            self.offset_angle = angle2

        volmdlr.core.CompositePrimitive3D.__init__(self,
                                                   primitives=self.discretization_points(),
                                                   name=name)

    def polygon_points(self, discretization_resolution: int):
        warnings.warn('polygon_points is deprecated,\
        please use discretization_points instead',
                      DeprecationWarning)
        return self.discretization_points(discretization_resolution)

    def _get_points(self):
        return self.discretization_points()

    points = property(_get_points)

    def to_2d(self, plane_origin, x, y):
        ps = self.start.to_2d(plane_origin, x, y)
        pi = self.interior.to_2d(plane_origin, x, y)
        pe = self.end.to_2d(plane_origin, x, y)
        center = self.center.to_2d(plane_origin, x, y)

        maj_dir2d = self.major_dir.to_2d(plane_origin, x, y)
        maj_dir2d.normalize()
        return ArcEllipse2D(ps, pi, pe, center, maj_dir2d, name=self.name)

    def length(self):
        return self.angle * math.sqrt(
            (self.Gradius ** 2 + self.Sradius ** 2) / 2)

    def normal_vector(self, abscissa):
        raise NotImplementedError

    def unit_normal_vector(self, abscissa):
        raise NotImplementedError

    def direction_vector(self, abscissa):
        raise NotImplementedError

    def unit_direction_vector(self, abscissa):
        raise NotImplementedError

    def reverse(self):
        return self.__class__(self.end.copy(),
                              self.interior.copy(),
                              self.start.copy(),
                              self.center.copy(),
                              self.major_dir.copy(),
                              self.name)

    def plot(self, ax=None):
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
        for px, py, pz in self.discretization_points():
            x.append(px)
            y.append(py)
            z.append(pz)

        ax.plot(x, y, z, 'k')
        return ax

    def plot2d(self, x3d: volmdlr.Vector3D = volmdlr.X3D, y3d: volmdlr.Vector3D = volmdlr.Y3D,
               ax=None, color='k'):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        # TODO: Enhance this plot
        l = self.length()
        x = []
        y = []
        for i in range(30):
            p = self.point_at_abscissa(i / (29.) * l)
            xi, yi = p.plane_projection2d(x3d, y3d)
            x.append(xi)
            y.append(yi)
        ax.plot(x, y, color=color)
        return ax

    def FreeCADExport(self, name, ndigits=6):
        xs, ys, zs = round(1000 * self.start, ndigits).vector
        xi, yi, zi = round(1000 * self.interior, ndigits).vector
        xe, ye, ze = round(1000 * self.end, ndigits).vector
        return '{} = Part.Arc(fc.Vector({},{},{}),fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(
            name, xs, ys, zs, xi, yi, zi, xe, ye, ze)

    def triangulation(self):
        return None
