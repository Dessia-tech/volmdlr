#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

from typing import List
import math
import volmdlr
from volmdlr.primitives import RoundedLineSegments
import matplotlib.pyplot as plt


class Line2D(volmdlr.primitives.Line):
    """
    Define an infinite line given by two points.
    """

    def __init__(self, point1, point2, *, name=''):
        volmdlr.primitives.Line.__init__(self, point1, point2, name=name)


    def To3D(self, plane_origin, x1, x2):
        p3D = [p.To3D(plane_origin, x1, x2) for p in self.points]
        return Line2D(*p3D, self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Line2D(
                *[p.Rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Line2D(
                *[p.Translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset, copy=False)

    def MPLPlot(self, ax=None, color='k', dashed=True):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')

        p1, p2 = self.points

        if version.parse(_mpl_version) >= version.parse('3.3.2'):
            if dashed:
                ax.axline(p1.vector, p2.vector, dashes=[30, 5, 10, 5])
            else:
                ax.axline(p1.vector, p2.vector)
        else:
            u = p2 - p1
            p3 = p1 - 3 * u
            p4 = p2 + 4 * u
            if dashed:
                ax.plot([p3[0], p4[0]], [p3[1], p4[1]], color=color,
                        dashes=[30, 5, 10, 5])
            else:
                ax.plot([p3[0], p4[0]], [p3[1], p4[1]], color=color)

        return ax

    def plot_data(self, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        p1, p2 = self.points
        u = p2 - p1
        p3 = p1 - 3 * u
        p4 = p2 + 4 * u
        return {'type': 'line',
                'data': [p3[0], p3[1],
                         p4[0], p4[1]],
                'color': color,
                'marker': marker,
                'size': stroke_width,
                'dash': dash,
                'opacity': opacity,
                'arrow': arrow
                }

    def CreateTangentCircle(self, point, other_line):
        """
        Computes the two circles that are tangent to 2 lines and intersect
        a point located on one of the two lines.
        """

        # point will be called I(x_I, y_I)
        # self will be (AB)
        # line will be (CD)

        if math.isclose(self.point_distance(point), 0, abs_tol=1e-10):
            I = Vector2D((point[0], point[1]))
            A = Vector2D((self.points[0][0], self.points[0][1]))
            B = Vector2D((self.points[1][0], self.points[1][1]))
            C = Vector2D((other_line.points[0][0], other_line.points[0][1]))
            D = Vector2D((other_line.points[1][0], other_line.points[1][1]))

        elif math.isclose(other_line.point_distance(point), 0, abs_tol=1e-10):
            I = Vector2D((point[0], point[1]))
            C = Vector2D((self.points[0][0], self.points[0][1]))
            D = Vector2D((self.points[1][0], self.points[1][1]))
            A = Vector2D((other_line.points[0][0], other_line.points[0][1]))
            B = Vector2D((other_line.points[1][0], other_line.points[1][1]))
        else:
            raise AttributeError("The point isn't on any of the two lines")

        # CHANGEMENT DE REPAIRE
        new_u = Vector2D((B - A))
        new_u.Normalize()
        new_v = new_u.NormalVector(unit=True)
        new_basis = Frame2D(I, new_u, new_v)

        new_A = new_basis.NewCoordinates(A)
        new_B = new_basis.NewCoordinates(B)
        new_C = new_basis.NewCoordinates(C)
        new_D = new_basis.NewCoordinates(D)

        if new_C[1] == 0 and new_D[1] == 0:
            # Segments are on the same line: no solution
            return None, None

        elif math.isclose(self.DirectionVector(unit=True).Dot(
                other_line.NormalVector(unit=True)), 0, abs_tol=1e-06):
            # Parallel segments: one solution

            segments_distance = abs(new_C[1] - new_A[1])
            r = segments_distance / 2
            new_circle_center = Point2D((0, npy.sign(new_C[1] - new_A[1]) * r))
            circle_center = new_basis.OldCoordinates(new_circle_center)
            circle = Circle2D(circle_center, r)

            return circle, None

        elif math.isclose(self.DirectionVector(unit=True).Dot(
                other_line.DirectionVector(unit=True)), 0, abs_tol=1e-06):
            # Perpendicular segments: 2 solution
            line_AB = Line2D(Point2D(new_A), Point2D(new_B))
            line_CD = Line2D(Point2D(new_C), Point2D(new_D))
            new_pt_K = Point2D.LinesIntersection(line_AB, line_CD)

            r = abs(new_pt_K[0])
            new_circle_center1 = Point2D((0, r))
            new_circle_center2 = Point2D((0, -r))
            circle_center1 = new_basis.OldCoordinates(new_circle_center1)
            circle_center2 = new_basis.OldCoordinates(new_circle_center2)
            circle1 = Circle2D(circle_center1, r)
            circle2 = Circle2D(circle_center2, r)

            return circle1, circle2

        # =============================================================================
        # LES SEGMENTS SONT QUELCONQUES
        #   => 2 SOLUTIONS
        # =============================================================================
        else:

            line_AB = Line2D(Point2D(new_A), Point2D(new_B))
            line_CD = Line2D(Point2D(new_C), Point2D(new_D))
            new_pt_K = Point2D.LinesIntersection(line_AB, line_CD)
            pt_K = Point2D(new_basis.OldCoordinates(new_pt_K))

            if pt_K == I:
                return None, None

            # CHANGEMENT DE REPERE:
            new_u2 = Vector2D(pt_K - I)
            new_u2.Normalize()
            new_v2 = new_u2.NormalVector(unit=True)
            new_basis2 = Frame2D(I, new_u2, new_v2)

            new_A = new_basis2.NewCoordinates(A)
            new_B = new_basis2.NewCoordinates(B)
            new_C = new_basis2.NewCoordinates(C)
            new_D = new_basis2.NewCoordinates(D)
            new_pt_K = new_basis2.NewCoordinates(pt_K)

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

            new_circle_center1 = Point2D((0, -r1))
            new_circle_center2 = Point2D((0, r2))

            circle_center1 = new_basis2.OldCoordinates(new_circle_center1)
            circle_center2 = new_basis2.OldCoordinates(new_circle_center2)

            if new_basis.NewCoordinates(circle_center1)[1] > 0:
                circle1 = Circle2D(circle_center1, r1)
                circle2 = Circle2D(circle_center2, r2)
            else:
                circle1 = Circle2D(circle_center2, r2)
                circle2 = Circle2D(circle_center1, r1)

            return circle1, circle2


class BSplineCurve2D(volmdlr.Primitive2D):
    def __init__(self, degree, control_points, knot_multiplicities, knots,
                 weights=None, periodic=False, name=''):
        Primitive2D.__init__(self, name=name)
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
            P = [(control_points[i][0], control_points[i][1]) for i in
                 range(len(control_points))]
            curve.ctrlpts = P
        else:
            Pw = [(control_points[i][0] * weights[i],
                   control_points[i][1] * weights[i], weights[i]) for i in
                  range(len(control_points))]
            curve.ctrlptsw = Pw
        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot] * knot_multiplicities[i])
        curve.knotvector = knot_vector
        curve.delta = 0.1
        curve_points = curve.evalpts

        self.curve = curve
        self.points = [Point2D((p[0], p[1])) for p in curve_points]

    def Length(self):
        # Approximately
        length = 0
        for k in range(0, len(self.points) - 1):
            length += (self.points[k] - self.points[k + 1]).Norm()
        return length

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        # copy paste from wire3D
        length = 0.
        primitives = []
        for k in range(0, len(self.points) - 1):
            primitives.append(
                LineSegment2D(self.points[k], self.points[k + 1]))
        for primitive in primitives:
            primitive_length = primitive.Length()
            if length + primitive_length >= curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(
                    curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError

    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        x = [p.vector[0] for p in self.points]
        y = [p.vector[1] for p in self.points]
        ax.plot(x, y, 'o-k')
        return fig, ax

    def To3D(self, plane_origin, x1, x2):
        control_points3D = [p.To3D(plane_origin, x1, x2) for p in
                            self.control_points]
        return BSplineCurve3D(self.degree, control_points3D,
                              self.knot_multiplicities, self.knots,
                              self.weights, self.periodic, self.name)

    def tessellation_points(self):
        return self.points


class LineSegment2D(volmdlr.primitives.LineSegment):
    """
    Define a line segment limited by two points
    """

    def __init__(self, start, end, *, name=''):
        volmdlr.primitives.Edge.__init__(self, start, end, name=name)

    def Length(self):
        return self.end.point_distance(self.start)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.start + self.unit_direction_vector() * curvilinear_abscissa

    def point_distance(self, point, return_other_point=False):
        """
        Computes the distance of a point to segment of line
        """
        if self.point1 == self.point2:
            if return_other_point:
                return 0, Point2D(point)
            return 0
        distance, point = LineSegment2DPointDistance(
            [p.vector for p in self.points], point.vector)
        if return_other_point:
            return distance, Point2D(point)
        return distance

    def point_projection(self, point):
        """
        If the projection falls outside the LineSegment2D, returns None.
        """
        point, curv_abs = Line2D.point_projection(self, point)
        if curv_abs < 0 or curv_abs > self.Length():
            return None, curv_abs
        return point, curv_abs

    def line_intersections(self, line):
        point = Point2D.LinesIntersection(self, line)
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
            return []

    def MPLPlot(self, ax=None, color='k', arrow=False, width=None,
                plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            # ax.set_aspect('equal')
        # else:
        #     fig = ax.figure

        p1, p2 = self.start, self.end
        if arrow:
            if plot_points:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color,
                        style='o-')
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color)

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
                        marker='o', linewidth=width)
            else:
                ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color,
                        linewidth=width)
        return ax

    def To3D(self, plane_origin, x1, x2):
        start = self.start.To3D(plane_origin, x1, x2)
        end = self.end.To3D(plane_origin, x1, x2)
        return volmdlr.primitives3D.LineSegment3D(start, end, self.name)

    def reverse(self):
        return LineSegment2D(self.end.copy(), self.points[0].copy())

    def to_line(self):
        return Line2D(*self.points)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return LineSegment2D(
                *[p.Rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return LineSegment2D(
                *[p.Translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset, copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return LineSegment2D(
                    *[frame.OldCoordinates(p) for p in self.points])
            else:
                self.points = [frame.OldCoordinates(p) for p in self.points]
        if side == 'new':
            if copy:
                return LineSegment2D(
                    *[frame.NewCoordinates(p) for p in self.points])
            else:
                self.points = [frame.NewCoordinates(p) for p in self.points]

    def plot_data(self, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        return {'type': 'line',
                'data': [self.points[0].vector[0], self.points[0].vector[1],
                         self.end.vector[0], self.end.vector[1]],
                'color': color,
                'marker': marker,
                'size': stroke_width,
                'dash': dash,
                'opacity': opacity,
                'arrow': arrow
                }

    def CreateTangentCircle(self, point, other_line):
        circle1, circle2 = Line2D.CreateTangentCircle(other_line, point, self)
        if circle1 is not None:
            point_J1, curv_abs1 = Line2D.point_projection(self, circle1.center)
            if curv_abs1 < 0. or curv_abs1 > self.Length():
                circle1 = None
        if circle2 is not None:
            point_J2, curv_abs2 = Line2D.point_projection(self, circle2.center)
            if curv_abs2 < 0. or curv_abs2 > self.Length():
                circle2 = None
        return circle1, circle2

    def tessellation_points(self):
        return [self.start, self.end]

    def polygon_points(self, min_x_density=None, min_y_density=None):
        n = 0# Number of points to insert between start and end
        if min_x_density:
            dx = abs(self.start[0]-self.end[0])
            n = max(n, math.floor(dx*min_x_density))
        if min_y_density:
            dy = abs(self.start[1]-self.end[1])
            n = max(n, math.floor(dy*min_y_density))

        if n:
            l = self.Length()
            return [self.PointAtCurvilinearAbscissa(i*l/(n+1)) for i in range(n+2)]
        else:
            return [self.start, self.end]



class Arc2D(volmdlr.Primitive2D):
    """
    angle: the angle measure always >= 0
    """

    def __init__(self, start, interior, end, name=''):
        volmdlr.Primitive2D.__init__(self, name)
        self.interior = interior
        self.start = start
        self.end = end
        xi, yi = interior.vector
        xe, ye = end.vector
        xs, ys = start.vector
        try:
            A = volmdlr.Matrix22(2 * (xs - xi), 2 * (ys - yi),
                         2 * (xs - xe), 2 * (ys - ye))
            b = - volmdlr.Vector2D((xi ** 2 + yi ** 2 - xs ** 2 - ys ** 2,
                            xe ** 2 + ye ** 2 - xs ** 2 - ys ** 2))
            inv_A = A.inverse()
            x = inv_A.vector_multiplication(b)
            self.center = volmdlr.Point2D(x.vector)
        except ValueError:
            A = npy.array([[2 * (xs - xi), 2 * (ys - yi)],
                           [2 * (xs - xe), 2 * (ys - ye)]])
            b = - npy.array([xi ** 2 + yi ** 2 - xs ** 2 - ys ** 2,
                             xe ** 2 + ye ** 2 - xs ** 2 - ys ** 2])
            self.center = Point2D(solve(A, b))

        r1 = self.start - self.center
        r2 = self.end - self.center
        ri = self.interior - self.center

        self.radius = r1.Norm()
        angle1 = math.atan2(r1.vector[1], r1.vector[0])
        anglei = math.atan2(ri.vector[1], ri.vector[0])
        angle2 = math.atan2(r2.vector[1], r2.vector[0])

        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei + volmdlr.two_pi) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + volmdlr.two_pi

        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2 + volmdlr.two_pi) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + volmdlr.two_pi

        if clockwise_path > trigowise_path:
            self.is_trigo = True
            self.angle1 = angle1
            self.angle2 = angle2
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle1 = angle2
            self.angle2 = angle1
            self.angle = clockwise_path

    def _get_points(self):
        return [self.start, self.interior, self.end]

    points = property(_get_points)


    def tessellation_points(self, resolution_for_circle=40):
        number_points_tesselation = math.ceil(
            resolution_for_circle * abs(self.angle) / 2 / math.pi)
        number_points_tesselation = max(number_points_tesselation, 5)
        l = self.Length()
        return [self.PointAtCurvilinearAbscissa(
            i / (number_points_tesselation - 1) * l) for i in
                range(number_points_tesselation)]

    def point_belongs(self, point):
        """
        Computes if the point belongs to the pizza slice drawn by the arc and its center
        """
        circle = Circle2D(self.center, self.radius)
        if not circle.point_belongs(point):
            return False
        vector_start = self.start - self.center
        vector_point = point - self.center
        vector_end = self.end - self.center
        if self.is_trigo:
            vector_start, vector_end = vector_end, vector_start
        arc_angle = clockwise_angle(vector_start, vector_end)
        point_angle = clockwise_angle(vector_start, vector_point)
        if point_angle <= arc_angle:
            return True

    def point_distance(self, point):
        vector_start = self.start - self.center
        vector_point = point - self.center
        vector_end = self.end - self.center
        if self.is_trigo:
            vector_start, vector_end = vector_end, vector_start
        arc_angle = clockwise_angle(vector_start, vector_end)
        point_angle = clockwise_angle(vector_start, vector_point)
        if point_angle <= arc_angle:
            return abs(
                LineSegment2D(point, self.center).Length() - self.radius)
        else:
            return min(LineSegment2D(point, self.start).Length(),
                       LineSegment2D(point, self.end).Length())

    def line_intersections(self, line):
        circle = Circle2D(self.center, self.radius)
        circle_intersection_points = circle.line_intersections(line)

        if circle_intersection_points is None:
            return None

        intersection_points = []
        for pt in circle_intersection_points:
            if self.point_belongs(pt):
                intersection_points.append(pt)
        return intersection_points

    def Length(self):
        return self.radius * abs(self.angle)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        if self.is_trigo:
            return self.start.Rotation(self.center,
                                       curvilinear_abscissa / self.radius)
            # return self.start.Rotation(self.center, curvilinear_abscissa*self.angle)
        else:
            return self.start.Rotation(self.center,
                                       -curvilinear_abscissa / self.radius)
            # return self.start.Rotation(self.center, -curvilinear_abscissa*self.angle)

    def MiddlePoint(self):
        l = self.Length()
        return self.PointAtCurvilinearAbscissa(0.5 * l)

    def Area(self):
        if self.angle2 < self.angle1:
            angle = self.angle2 + volmdlr.two_pi - self.angle1
        else:
            angle = self.angle2 - self.angle1
        return self.radius ** 2 * angle / 2

    def CenterOfMass(self):
        #        u=self.middle.vector-self.center.vector
        u = self.MiddlePoint() - self.center
        u.Normalize()
        alpha = abs(self.angle)
        return self.center + 4 / (3 * alpha) * self.radius * math.sin(
            alpha * 0.5) * u

    def MPLPlot(self, ax=None, color='k', plot_points=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')

        if plot_points:
            for p in [self.center, self.start, self.interior, self.end]:
                p.MPLPlot(ax=ax)

        pc = self.center.vector
        ax.add_patch(volmdlr.Arc(pc, 2 * self.radius, 2 * self.radius, angle=0,
                             theta1=self.angle1 * 0.5 / math.pi * 360,
                           theta2=self.angle2 * 0.5 / math.pi * 360,
                            color=color))

        return ax

    def To3D(self, plane_origin, x, y):
        ps = self.start.To3D(plane_origin, x, y)
        pi = self.interior.To3D(plane_origin, x, y)
        pe = self.end.To3D(plane_origin, x, y)

        return volmdlr.primitives3D.Arc3D(ps, pi, pe, name=self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Arc2D(*[p.Rotation(center, angle, copy=True) for p in
                           [self.start, self.interior, self.end]])
        else:
            self.__init__(*[p.Rotation(center, angle, copy=True) for p in
                            [self.start, self.interior, self.end]])


    def Translation(self, offset, copy=True):
        if copy:
            return Arc2D(*[p.Translation(offset, copy=True) for p in
                           [self.start, self.interior, self.end]])
        else:
            self.__init__(*[p.Translation(offset, copy=True) for p in
                            [self.start, self.interior, self.end]])


    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            return Arc2D(*[p.frame_mapping(frame, side, copy=True) for p in
                           [self.start, self.interior, self.end]])
        else:
            self.__init__(*[p.frame_mapping(frame, side, copy=True) for p in
                            [self.start, self.interior, self.end]])


    def SecondMomentArea(self, point):
        """
        Second moment area of part of disk
        """
        if self.angle2 < self.angle1:
            angle2 = self.angle2 + two_pi

        else:
            angle2 = self.angle2
        angle1 = self.angle1

        Ix = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                    math.sin(2 * angle1) - math.sin(2 * angle2)))
        Iy = self.radius ** 4 / 8 * (angle2 - angle1 + 0.5 * (
                    math.sin(2 * angle2) - math.sin(2 * angle1)))
        Ixy = self.radius ** 4 / 8 * (
                    math.cos(angle1) ** 2 - math.cos(angle2) ** 2)
        Ic = npy.array([[Ix, Ixy], [Ixy, Iy]])
        return geometry.Huygens2D(Ic, self.Area(), self.center, point)

    def Discretise(self, num=10):
        list_node = []
        if (self.angle1 < 0) and (self.angle2 > 0):
            delta_angle = -self.angle1 + self.angle2
        elif (self.angle1 > 0) and (self.angle2 < 0):
            delta_angle = (2 * npy.pi + self.angle2) - self.angle1
        else:
            delta_angle = self.angle2 - self.angle1
        for angle in npy.arange(self.angle1, self.angle1 + delta_angle,
                                delta_angle / (num * 1.)):
            list_node.append(Point2D(self.center + self.radius * Vector2D(
                (npy.cos(angle), npy.sin(angle)))))
        list_node.append(Point2D(self.center + self.radius * Vector2D((npy.cos(
            self.angle1 + delta_angle), npy.sin(self.angle1 + delta_angle)))))
        if list_node[0] == self.start:
            return list_node
        else:
            return list_node[::-1]

    def plot_data(self, marker=None, color='black', stroke_width=1, dash=False,
                  opacity=1):
        list_node = self.Discretise()
        data = []
        for nd in list_node:
            data.append({'x': nd.vector[0], 'y': nd.vector[1]})
        return {'type': 'arc',
                'cx': self.center.vector[0],
                'cy': self.center.vector[1],
                'data': data,
                'r': self.radius,
                'color': color,
                'opacity': opacity,
                'size': stroke_width,
                'dash': None,
                'marker': marker,
                'angle1': self.angle1,
                'angle2': self.angle2, }

    def copy(self):
        return Arc2D(self.start.copy(),
                     self.interior.copy(),
                     self.end.copy())

    def split(self, split_point: volmdlr.Point2D):
        raise NotImplementedError
        return [Arc2D(self.start, self.split_point)]

    def polygon_points(self, points_per_radian=10, min_x_density=None,
                       min_y_density=None):

        number_points = math.ceil(self.angle*points_per_radian)
        densities = []
        for d in [min_x_density, min_y_density]:
            if d:
                densities.append(d)
        if densities:
            number_points = max(number_points,
                                min(densities)*self.angle*self.radius)
        l = self.Length()
        return [self.PointAtCurvilinearAbscissa(i*l/number_points)\
                for i in range(number_points+1)]


class ArcEllipse2D(volmdlr.Primitive2D):
    """

    """

    def __init__(self, start, interior, end, center, major_dir, name='',
                 extra=None):
        self.start = start
        self.interior = interior
        self.end = end
        self.center = center
        self.extra = extra
        self.major_dir = major_dir
        self.minor_dir = self.major_dir.deterministic_unit_normal_vector()

        frame = Frame2D(self.center, self.major_dir, self.minor_dir)
        start_new, end_new = frame.NewCoordinates(
            self.start), frame.NewCoordinates(self.end)
        interior_new, center_new = frame.NewCoordinates(
            self.interior), frame.NewCoordinates(self.center)

        #### from : https://math.stackexchange.com/questions/339126/how-to-draw-an-ellipse-if-a-center-and-3-arbitrary-points-on-it-are-given
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
            theta = 0
            c1 = C[0] + C[1]
            c2 = (C[1] - C[0]) / math.cos(2 * theta)
            gdaxe = math.sqrt((2 / (c1 - c2)))
            ptax = math.sqrt((2 / (c1 + c2)))
            return theta, gdaxe, ptax

        if start == end:
            extra_new = frame.NewCoordinates(self.extra)
            theta, A, B = theta_A_B(start_new, extra_new, interior_new,
                                    center_new)
        else:
            theta, A, B = theta_A_B(start_new, interior_new, end_new,
                                    center_new)

        self.Gradius = A
        self.Sradius = B
        self.theta = theta

        # Angle pour start
        u1, u2 = start_new.vector[0] / self.Gradius, start_new.vector[
            1] / self.Sradius
        angle1 = sin_cos_angle(u1, u2)
        # Angle pour end
        u3, u4 = end_new.vector[0] / self.Gradius, end_new.vector[
            1] / self.Sradius
        angle2 = sin_cos_angle(u3, u4)
        # Angle pour interior
        u5, u6 = interior_new.vector[0] / self.Gradius, interior_new.vector[
            1] / self.Sradius
        anglei = sin_cos_angle(u5, u6)

        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei + two_pi) - angle1
            clockwise_path = angle1 - anglei
        else:
            trigowise_path = anglei - angle1
            clockwise_path = angle1 - anglei + two_pi

        # Going trigo wise from interior to interior
        if angle2 < anglei:
            trigowise_path += (angle2 + two_pi) - anglei
            clockwise_path += anglei - angle2
        else:
            trigowise_path += angle2 - anglei
            clockwise_path += anglei - angle2 + two_pi

        if clockwise_path > trigowise_path:
            self.is_trigo = True
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle = clockwise_path

        if self.start == self.end or self.angle == 0:
            self.angle = two_pi

        if self.is_trigo:  # sens trigo
            self.offset_angle = angle1
        else:
            self.offset_angle = angle2

        Primitive2D.__init__(self, name=name)

    def _get_points(self):
        return self.tessellation_points()

    points = property(_get_points)

    def tessellation_points(self, resolution_for_ellipse=40):
        number_points_tesselation = math.ceil(
            resolution_for_ellipse * abs(0.5 * self.angle / math.pi))

        frame2d = Frame2D(self.center, self.major_dir, self.minor_dir)

        tessellation_points_2D = [(Point2D((self.Gradius * math.cos(
            self.offset_angle + self.angle * i / (number_points_tesselation)),
                                            self.Sradius * math.sin(
                                                self.offset_angle + self.angle * i / (
                                                    number_points_tesselation)))))
                                  for i in
                                  range(number_points_tesselation + 1)]

        global_points = []
        for pt in tessellation_points_2D:
            global_points.append(frame2d.OldCoordinates(pt))

        return global_points

    def To3D(self, plane_origin, x, y):
        ps = self.start.To3D(plane_origin, x, y)
        pi = self.interior.To3D(plane_origin, x, y)
        pe = self.end.To3D(plane_origin, x, y)
        pc = self.center.To3D(plane_origin, x, y)
        if self.extra is None:
            pextra = None
        else:
            pextra = self.extra.To3D(plane_origin, x, y)
        if ps == pe:
            p3 = pextra
        else:
            p3 = pe
        plane = volmdlr.surfaces3d.Plane3D.from_3_points(ps, pi, p3)
        n = plane.normal
        major_dir = self.major_dir.To3D(plane_origin, x, y)
        major_dir.Normalize()

        return ArcEllipse3D(ps, pi, pe, pc, major_dir, normal=n,
                            name=self.name, extra=pextra)

    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = ax.figure

        self.interior.MPLPlot(ax=ax, color='m')
        self.start.MPLPlot(ax=ax, color='r')
        self.end.MPLPlot(ax=ax, color='b')
        self.center.MPLPlot(ax=ax, color='y')

        x = []
        y = []
        for px, py in self.tessellation_points():
            x.append(px)
            y.append(py)

        plt.plot(x, y, 'k')
        return fig, ax
        


class Wire2D(volmdlr.CompositePrimitive2D):
    """
    A collection of simple primitives, following each other making a wire
    """

    def __init__(self, primitives, name=''):
        volmdlr.CompositePrimitive2D.__init__(self, primitives, name)

    # TODO: method to check if it is a wire

    def Length(self):
        length = 0.
        for primitive in self.primitives:
            length += primitive.Length()
        return length

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa: float):
        length = 0.
        for primitive in self.primitives:
            primitive_length = primitive.Length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(
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


class OpenedRoundedLineSegments2D(RoundedLineSegments, Wire2D):
    closed = False

    def __init__(self, points, radius, adapt_radius=False, name=''):
        primitives = RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.LineSegment2D,
                                                  volmdlr.Arc2D,
                                                  closed=False,
                                                  adapt_radius=adapt_radius,
                                                  name='')

        volmdlr.Wire2D.__init__(self, primitives, name)

    def ArcFeatures(self, ipoint):
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
        dist1 = (pt1 - pti).Norm()
        dist2 = (pt2 - pti).Norm()
        dist3 = (pt1 - pt2).Norm()
        alpha = math.acos(
            -(dist3 ** 2 - dist1 ** 2 - dist2 ** 2) / (2 * dist1 * dist2)) / 2.
        dist = radius / math.tan(alpha)

        u1 = (pt1 - pti) / dist1
        u2 = (pt2 - pti) / dist2

        p3 = pti + u1 * dist
        p4 = pti + u2 * dist

        w = (u1 + u2)
        if w != volmdlr.Vector2D((0, 0)):
            w.Normalize()

        v1 = u1.deterministic_unit_normal_vector()
        if v1.Dot(w) < 0:
            v1 = -v1

        pc = p3 + v1 * radius
        pm = pc - radius * w

        return p3, pm, p4, dist, alpha

    def Rotation(self, center, angle, copy=True):
        if copy:
            return self.__class__([p.Rotation(center, angle, copy=True) \
                                   for p in self.points],
                                  self.radius,
                                  adapt_radius=self.adapt_radius,
                                  name=self.name)
        else:
            self.__init__(
                [p.Rotation(center, angle, copy=True) for p in self.points],
                self.radius,
                adapt_radius=self.adapt_radius, name=self.name)

    def Translation(self, offset, copy=True):
        if copy:
            return self.__class__(
                [p.Translation(offset, copy=True) for p in self.points],
                self.radius, adapt_radius=self.adapt_radius, name=self.name)
        else:
            self.__init__(
                [p.Translation(offset, copy=True) for p in self.points],
                self.radius, adapt_radius=self.adapt_radius, name=self.name)

    def Offset(self, offset):
        nb = len(self.points)
        vectors = []
        for i in range(nb - 1):
            v1 = volmdlr.Vector2D((self.points[i + 1] - self.points[i]))
            v2 = volmdlr.Vector2D((self.points[i] - self.points[i + 1]))
            v1.Normalize()
            v2.Normalize()
            vectors.append(v1)
            vectors.append(v2)

        if self.closed:
            v1 = volmdlr.Vector2D((self.points[0] - self.points[-1]))
            v2 = volmdlr.Vector2D((self.points[-1] - self.points[0]))
            v1.Normalize()
            v2.Normalize()
            vectors.append(v1)
            vectors.append(v2)

        offset_vectors = []
        new_radii = {}
        offset_points = []

        for i in range((not self.closed), nb - (not self.closed)):

            check = False
            ni = vectors[2 * i - 1] + vectors[2 * i]
            if ni == volmdlr.Vector2D((0, 0)):
                ni = vectors[2 * i]
                ni = ni.NormalVector()
                offset_vectors.append(ni)
            else:
                ni.Normalize()
                if ni.Dot(vectors[2 * i - 1].NormalVector()) > 0:
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

            normal_vector1 = - vectors[2 * i - 1].NormalVector()
            normal_vector2 = vectors[2 * i].NormalVector()
            normal_vector1.Normalize()
            normal_vector2.Normalize()
            alpha = math.acos(normal_vector1.Dot(normal_vector2))

            offset_point = self.points[i] + offset / math.cos(alpha / 2) * \
                           offset_vectors[i - (not self.closed)]
            offset_points.append(offset_point)

        if not self.closed:
            n1 = vectors[0].NormalVector(unit=True)
            offset_vectors.insert(0, n1)
            offset_points.insert(0,
                                 self.points[0] + offset * offset_vectors[0])

            n_last = vectors[-1].NormalVector(unit=True)
            n_last = - n_last
            offset_vectors.append(n_last)
            offset_points.append(self.points[-1] + offset * offset_vectors[-1])

        return self.__class__(offset_points, new_radii,
                              adapt_radius=self.adapt_radius)

    def OffsetSingleLine(self, line_index, offset):
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
                    dir_vec_1.Normalize()
                    dir_vec_2 = dir_vec_1
                    dont_add_last_point = True
                # The offset line is the first one
                elif i == 0:
                    dir_vec_2 = volmdlr.Vector2D(
                        self.points[i + 1] - self.points[i + 2])
                    dir_vec_2.Normalize()
                    if not self.closed:
                        dir_vec_1 = dir_vec_2
                    else:
                        dir_vec_1 = volmdlr.Vector2D(
                            point - self.points[i - 1])
                        dir_vec_1.Normalize()
                # Closed RLS2D and the offset line is the last one
                elif i == len(self.points) - 1:
                    dir_vec_1 = volmdlr.Vector2D(point - self.points[i - 1])
                    dir_vec_1.Normalize()
                    dir_vec_2 = volmdlr.Vector2D(
                        self.points[0] - self.points[1])
                    dir_vec_2.Normalize()
                    dont_add_last_point = True
                else:
                    dir_vec_1 = volmdlr.Vector2D(point - self.points[i - 1])
                    dir_vec_1.Normalize()
                    dir_vec_2 = volmdlr.Vector2D(
                        self.points[i + 1] - self.points[i + 2])
                    dir_vec_2.Normalize()

                if self.closed and line_index == len(self.points) - 1:
                    normal_vector = volmdlr.Vector2D(
                        self.points[0] - point).NormalVector(unit=True)
                else:
                    normal_vector = volmdlr.Vector2D(
                        self.points[i + 1] - point).NormalVector(unit=True)

                alpha1 = math.acos(dir_vec_1.Dot(normal_vector))
                alpha2 = math.acos(dir_vec_2.Dot(normal_vector))

                # If 3 segments are aligned and the middle one have to be offset
                if math.isclose(math.cos(alpha1), 0,
                                abs_tol=1e-9) or math.isclose(math.cos(alpha2),
                                                              0, abs_tol=1e-9):
                    print('ca sort direct')
                    print('direction vector', dir_vec_1, dir_vec_2)
                    print('normal vector', normal_vector)
                    print('alpha', alpha1 * 180 / math.pi,
                          alpha2 * 180 / math.pi)
                    print(point, self.points[i + 1])
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

    def OffsetLines(self, line_indexes, offset):
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
            dir_vec_1 = volmdlr.Vector2D((self.points[line_indexes[0]] -
                                          self.points[line_indexes[0] - 1]))

        if line_indexes[-1] == len(self.points) - (2 - self.closed):
            if not self.closed:
                pass
            else:
                dir_vec_2 = volmdlr.Vector2D((self.points[0] - self.points[1]))
        elif self.closed and line_indexes[-1] == len(self.points) - 2:
            dir_vec_2 = volmdlr.Vector2D(
                (self.points[line_indexes[-1] + 1] - self.points[0]))
        else:
            dir_vec_2 = volmdlr.Vector2D((self.points[line_indexes[-1] + 1] -
                                          self.points[line_indexes[-1] + 2]))

        if dir_vec_1 is None:
            dir_vec_1 = dir_vec_2
        if dir_vec_2 is None:
            dir_vec_2 = dir_vec_1

        dir_vec_1.Normalize()
        dir_vec_2.Normalize()

        # =============================================================================
        # COMPUTES THE ANGLE BETWEEN THE NORMAL VECTOR OF THE SURFACE TO OFFSET AND
        # THE DIRECTIVE VECTOR IN ORDER TO SET THE NEW POINT AT THE RIGHT DISTANCE
        # =============================================================================
        normal_vectors = []
        for index in line_indexes:
            if index == len(self.points) - 1:
                normal_vectors.append(volmdlr.Vector2D(
                    self.points[0] - self.points[index]).NormalVector(
                    unit=True))
            else:
                normal_vectors.append(volmdlr.Vector2D(
                    self.points[index + 1] - self.points[index]).NormalVector(
                    unit=True))

        dot1 = dir_vec_1.Dot(normal_vectors[0])
        dot2 = dir_vec_2.Dot(normal_vectors[-1])

        if math.isclose(dot1, 0, abs_tol=1e-9):
            # call function considering the line before, because the latter and
            # the first offset segment are parallel
            return self.OffsetLines([line_indexes[0] - 1] + line_indexes,
                                    offset)
        if math.isclose(dot2, 0, abs_tol=1e-9):
            # call function considering the line after, because the latter and
            # the last offset segment are parallel
            return self.OffsetLines(line_indexes + [line_indexes[-1] + 1],
                                    offset)

        distance_dir1 = offset / dot1
        distance_dir2 = offset / dot2

        if len(line_indexes) > 1:
            intersection = volmdlr.Point2D.LinesIntersection(
                volmdlr.Line2D(self.points[line_indexes[0]],
                               self.points[line_indexes[0]] + dir_vec_1),
                volmdlr.Line2D(self.points[line_indexes[-1] + 1],
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
            if dir_vec_1.Dot(normal_vectors[nb + 1]) < 0:
                coeff_vec_1 = - coeff_vec_1
            if dir_vec_2.Dot(normal_vectors[nb + 1]) < 0:
                coeff_vec_2 = - coeff_vec_2
            index_dir_vector = coeff_vec_1 * vec1 + coeff_vec_2 * vec2
            index_dot = index_dir_vector.Dot(normal_vectors[nb + 1])
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

class Contour2D(Wire2D):
    """
    A collection of 2D primitives forming a closed wire2D
    TODO : CenterOfMass and SecondMomentArea should be changed accordingly to
    Area considering the triangle drawn by the arcs
    """
    _non_serializable_attributes = ['internal_arcs', 'external_arcs',
                                    'polygon', 'straight_line_contour_polygon']

    def __init__(self, primitives, name=''):
        Wire2D.__init__(self, primitives, name)
        self._utd_analysis = False
        self.tessel_points = self.clean_points()

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
        return (Point2D((xmin, ymin)), Point2D((xmax, ymax)))

    def to_2d(self, plane3d, name=None):
        primitives3D = [p.to_2d(plane3d) for p in self.primitives]
        return Contour3D(edges=primitives3D, name=name)

    def To3D(self, plane_origin, x, y, name=None):
        if name is None:
            name = '3D of {}'.format(self.name)
        primitives3D = [p.To3D(plane_origin, x, y) for p in self.primitives]
        return Contour3D(primitives=primitives3D, name=name)

    def Area(self):
        if len(self.primitives) == 1:
            return self.primitives[0].Area()

        A = self.polygon.Area()

        for arc in self.internal_arcs:
            triangle = Polygon2D([arc.start, arc.center, arc.end])
            A = A - arc.Area() + triangle.Area()
        for arc in self.external_arcs:
            triangle = Polygon2D([arc.start, arc.center, arc.end])
            A = A + arc.Area() - triangle.Area()

        return A

    def CenterOfMass(self):
        if len(self.primitives) == 1:
            return self.primitives[0].CenterOfMass()

        area = self.polygon.Area()
        if area > 0.:
            c = area * self.polygon.CenterOfMass()
        else:
            c = O2D

        for arc in self.internal_arcs:
            arc_area = arc.Area()
            c -= arc_area * arc.CenterOfMass()
            area -= arc_area
        for arc in self.external_arcs:
            arc_area = arc.Area()
            c += arc_area * arc.CenterOfMass()
            area += arc_area
        if area != 0:
            return c / area
        else:
            return False

    def SecondMomentArea(self, point):
        if len(self.primitives) == 1:
            return self.primitives[0].SecondMomentArea(point)

        A = self.polygon.SecondMomentArea(point)
        for arc in self.internal_arcs:
            A -= arc.SecondMomentArea(point)
        for arc in self.external_arcs:
            A += arc.SecondMomentArea(point)

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
        return Point2D((x, y))

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
                d1, d2 = (points_to_add[0] - points[0]).Norm(), (
                            points_to_add[0] - points[-1]).Norm()
                d3, d4 = (points_to_add[-1] - points[0]).Norm(), (
                            points_to_add[-1] - points[-1]).Norm()
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
            p = Point2D.random(xmin, xmax, ymin, ymax)
            if self.point_belongs(p):
                return p
    # def line_intersections(self, line:Line2D) -> List[Tuple[Point2D, Primitive2D]]:
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
            if intersections[0][0].__class__.__name__ == 'Point2D' and \
                    intersections[1][0].__class__.__name__ == 'Point2D':
                ip1, ip2 = sorted([self.primitives.index(intersections[0][1]),
                                   self.primitives.index(intersections[1][1])])

                sp11, sp12 = intersections[0][1].split(intersections[0][0])
                sp21, sp22 = intersections[1][1].split(intersections[1][0])

                primitives1 = self.primitives[:ip1]
                primitives1.append(sp11)
                primitives1.append(LineSegment2D(intersections[0][0],
                                                 intersections[1][0]))
                primitives1.append(sp22)
                primitives1.extend(self.primitives[ip2 + 1:])

                primitives2 = self.primitives[ip1 + 1:ip2]
                primitives2.append(sp21)
                primitives2.append(LineSegment2D(intersections[1][0],
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
            cut_line = Line2D(Point2D((xi, 0)), Point2D((xi, 1)))

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
                p = Point2D((xi, yi))
                if self.point_belongs(p):
                    point_index[p] = ip
                    points.append(p)
                    ip += 1

        for i in range(n):
            for j in range(m):
                p1 = Point2D((x[i], y[j]))
                p2 = Point2D((x[i + 1], y[j]))
                p3 = Point2D((x[i + 1], y[j + 1]))
                p4 = Point2D((x[i], y[j + 1]))
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

        return DisplayMesh2D(points, triangles)

    
class ClosedRoundedLineSegments2D(OpenedRoundedLineSegments2D, Contour2D):
    """
    :param points: Points used to draw the wire 
    :type points: List of Point2D
    :param radius: Radius used to connect different parts of the wire
    :type radius: {position1(n): float which is the radius linked the n-1 and n+1 points, position2(n+1):...}
    """
    closed = True
    def __init__(self, points, radius, adapt_radius=False, name=''):
        primitives = RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.primitives2D.LineSegment2D,
                                                  volmdlr.primitives2D.Arc2D,
                                                  closed=True,
                                                  adapt_radius=adapt_radius, name='')

        volmdlr.primitives2D.Contour2D.__init__(self, primitives, name)

class Measure2D(LineSegment2D):
    def __init__(self, point1, point2, label='', unit='mm', type_='distance'):
        """
        :param unit: 'mm', 'm' or None. If None, the distance won't be in the label

        """
        # TODO: offset parameter
        LineSegment2D.__init__(self, point1, point2)
        self.label = label
        self.unit = unit
        self.type_ = type_

    def MPLPlot(self, ax, ndigits=6):
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
            arrow = FancyArrowPatch((x1, y1), (x2, y2),
                                    arrowstyle='<|-|>,head_length=10,head_width=5',
                                    shrinkA=0, shrinkB=0,
                                    color='k')
        elif self.type_ == 'radius':
            arrow = FancyArrowPatch((x1, y1), (x2, y2),
                                    arrowstyle='-|>,head_length=10,head_width=5',
                                    shrinkA=0, shrinkB=0,
                                    color='k')

        ax.add_patch(arrow)
        if x2 - x1 == 0.:
            theta = 90.
        else:
            theta = math.degrees(math.atan((y2 - y1) / (x2 - x1)))
        ax.text(xm, ym, label, va='bottom', ha='center', rotation=theta)


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

    def Area(self):

        x = [point.vector[0] for point in self.points]
        y = [point.vector[1] for point in self.points]

        return 0.5 * npy.abs(
            npy.dot(x, npy.roll(y, 1)) - npy.dot(y, npy.roll(x, 1)))

    def CenterOfMass(self):

        x = [point.vector[0] for point in self.points]
        y = [point.vector[1] for point in self.points]

        xi_xi1 = x + npy.roll(x, -1)
        yi_yi1 = y + npy.roll(y, -1)
        xi_yi1 = npy.multiply(x, npy.roll(y, -1))
        xi1_yi = npy.multiply(npy.roll(x, -1), y)

        a = 0.5 * npy.sum(xi_yi1 - xi1_yi)  # signed area!
        #        a=self.Area()
        if not math.isclose(a, 0, abs_tol=1e-08):
            cx = npy.sum(npy.multiply(xi_xi1, (xi_yi1 - xi1_yi))) / 6. / a
            cy = npy.sum(npy.multiply(yi_yi1, (xi_yi1 - xi1_yi))) / 6. / a
            return Point2D((cx, cy))

        else:
            raise NotImplementedError

    def PointBelongs(self, point):
        """
        Ray casting algorithm copied from internet...
        """
        return PolygonPointBelongs(point.vector,
                                   [p.vector for p in self.points])

    def SecondMomentArea(self, point):
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
            lines.append(LineSegment2D(p1, p2))
        return lines

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Polygon2D(
                [p.Rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Polygon2D(
                [p.Translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset, copy=False)

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

                        line1 = LineSegment2D(
                            Point2D(self.points[segment1[0]]),
                            Point2D(self.points[segment1[1]]))
                        line2 = LineSegment2D(
                            Point2D(self.points[segment2[0]]),
                            Point2D(self.points[segment2[1]]))

                        p, a, b = Point2D.LinesIntersection(line1, line2, True)

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
        ymax, pos_ymax = max_pos([pt.vector[1] for pt in points])
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

    def MPLPlot(self, ax=None, color='k',
                plot_points=False, point_numbering=False):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')

        for ls in self.line_segments:
            ls.MPLPlot(ax=ax ,color=color)

        if plot_points or point_numbering:
            for point in self.points:
                point.MPLPlot(ax=ax, color=color)

        if point_numbering:
            for ip, point in enumerate(self.points):
                ax.text(*point, 'point {}'.format(ip+1),
                        ha='center', va='top')

        ax.margins(0.1)
        plt.show()

        return ax

class Surface2D(volmdlr.Primitive2D):
    """
    A surface bounded by an outer contour
    """
    def __init__(self, outer_contour: Contour2D,
                 inner_contours: List[Contour2D],
                 name:str='name'):
        self.outer_contour = outer_contour
        self.inner_contours = inner_contours

        Primitive2D.__init__(self, name=name)

    def triangulation(self, min_x_density=None, min_y_density=None):
        outer_polygon = self.outer_contour.polygonization(min_x_density=15, min_y_density=12)
        # ax2 = outer_polygon.MPLPlot(color='r', point_numbering=True)
        points = outer_polygon.points
        vertices = [p.vector for p in points]
        n = len(outer_polygon.points)
        segments = [(i, i+1) for i in range(n-1)]
        segments.append((n-1, 0))
        point_index = {p:i for i,p in enumerate(points)}
        holes = []

        for inner_contour in self.inner_contours:
            inner_polygon = inner_contour.polygonization()
            # inner_polygon.MPLPlot(ax=ax2)
            for point in inner_polygon.points:
                if not point in point_index:
                    points.append(point)
                    vertices.append(point.vector)
                    point_index[point] = n
                    n += 1
            for point1, point2 in zip(inner_polygon.points[:-1],
                                      inner_polygon.points[1:]):
                segments.append((point_index[point1],
                                 point_index[point2]))
            segments.append((point_index[inner_polygon.points[-1]],
                             point_index[inner_polygon.points[0]]))
            holes.append(inner_contour.random_point_inside().vector)


        tri = {'vertices': npy.array(vertices).reshape((-1, 2)),
               'segments': npy.array(segments).reshape((-1, 2)),
               }
        if holes:
            tri['holes'] = npy.array(holes).reshape((-1, 2))

        t = triangle.triangulate(tri, 'p')
        triangles = t['triangles'].tolist()

        return DisplayMesh2D(points, triangles=triangles, edges=None)

    def plot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        self.outer_contour.MPLPlot(ax=ax)
        for inner_contour in self.inner_contours:
            inner_contour.MPLPlot(ax=ax)

        ax.set_aspect('equal')
        ax.margins(0.1)
        return ax




class Circle2D(Contour2D):
    _non_serializable_attributes = ['internal_arcs', 'external_arcs',
                                    'polygon', 'straight_line_contour_polygon',
                                    'primitives', 'basis_primitives']

    def __init__(self, center: volmdlr.Point2D, radius: float, name: str = ''):
        self.center = center
        self.radius = radius
        self.angle = two_pi
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
        return [self.center + self.radius * math.cos(teta) * Vector2D(
            (1, 0)) + self.radius * math.sin(teta) * Vector2D((0, 1)) \
                for teta in npy.linspace(0, two_pi, resolution + 1)][:-1]

    def point_belongs(self, point, tolerance=1e-9):
        return point.point_distance(self.center) <= self.radius + tolerance

    def line_intersections(self, line):
        V = Vector2D((line.points[1] - line.points[0]).vector)
        Q = Vector2D(self.center.vector)
        P1 = Vector2D(line.points[0].vector)

        a = V.Dot(V)
        b = 2 * V.Dot(P1 - Q)
        c = P1.Dot(P1) + Q.Dot(Q) - 2 * P1.Dot(Q) - self.radius ** 2

        disc = b ** 2 - 4 * a * c
        if disc < 0:
            return []

        if math.isclose(disc, 0, abs_tol=1e-8):
            t = -b / (2 * a)
            if line.__class__ is Line2D:
                return [Point2D((P1 + t * V).vector)]
            else:
                if 0 <= t <= 1:
                    return [Point2D((P1 + t * V).vector)]
                else:
                    return []

        sqrt_disc = math.sqrt(disc)
        t1 = (-b + sqrt_disc) / (2 * a)
        t2 = (-b - sqrt_disc) / (2 * a)
        if line.__class__ is Line2D:
            return [Point2D((P1 + t1 * V).vector),
                    Point2D((P1 + t2 * V).vector)]
        else:
            if not (0 <= t1 <= 1 or 0 <= t2 <= 1):
                return []
            elif 0 <= t1 <= 1 and not 0 <= t2 <= 1:
                return [Point2D((P1 + t1 * V).vector)]
            elif not 0 <= t1 <= 1 and 0 <= t2 <= 1:
                return [Point2D((P1 + t2 * V).vector)]
            else:
                [Point2D((P1 + t1 * V).vector), Point2D((P1 + t2 * V).vector)]

    def Length(self):
        return two_pi * self.radius

    def MPLPlot(self, ax=None, linestyle='-', color='k', linewidth=1):
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

    def To3D(self, plane_origin, x, y):
        normal = x.Cross(y)
        center3d = self.center.To3D(plane_origin, x, y)
        return Circle3D(Frame3D(center3d, x, y, normal),
                        self.radius, self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Circle2D(self.center.Rotation(center, angle, copy=True),
                            self.radius)
        else:
            self.center.Rotation(center, angle, copy=False)
            self.utd_geo_points = False

    def Translation(self, offset, copy=True):
        if copy:
            return Circle2D(self.center.Translation(offset, copy=True),
                            self.radius)
        else:
            self.center.Translation(offset, copy=False)
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

    def Area(self):
        return math.pi * self.radius ** 2

    def SecondMomentArea(self, point):
        """
        Second moment area of part of disk
        """
        I = math.pi * self.radius ** 4 / 4
        Ic = npy.array([[I, 0], [0, I]])
        return geometry.Huygens2D(Ic, self.Area(), self.center, point)

    def CenterOfMass(self):
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

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        start = self.center + self.radius * X3D
        return start.Rotation(self.center,
                              curvilinear_abscissa / self.radius)

    def triangulation(self, n=35):
        l = self.Length()
        points = [self.PointAtCurvilinearAbscissa(l * i / n) for i in range(n)]
        points.append(self.center)
        triangles = [(i, i + 1, n) for i in range(n - 1)] + [(n - 1, 0, n)]
        return DisplayMesh(points, triangles)

    def polygon_points(self, points_per_radian=10, min_x_density=None,
                       min_y_density=None):
        return Arc2D.polygon_points(self, points_per_radian=points_per_radian,
                                    min_x_density=min_x_density,
                                    min_y_density=min_y_density)
