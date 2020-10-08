#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common primitives 3D
"""

import math

import numpy as npy
npy.seterr(divide='raise')

import volmdlr
import volmdlr.primitives
import volmdlr.surfaces3d
from typing import Tuple



import matplotlib.pyplot as plt


class Line3D(volmdlr.primitives.Line):
    _non_eq_attributes = ['name', 'basis_primitives', 'bounding_box']

    """
    Define an infinite line passing through the 2 points
    """

    def __init__(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D,
                 name: str=''):
        volmdlr.primitives.Line.__init__(self, point1, point2, name=name)
        self.bounding_box = self._bounding_box()


    def _bounding_box(self):
        points = [self.point1, self.point2]

        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])

        return volmdlr.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.points[0] + (
                self.point2 - self.points[0]) * curvilinear_abscissa

    def MPLPlot(self, ax=None, color='k', dashed=True):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = ax.figure

        # Line segment
        x = [p.vector[0] for p in self.points]
        y = [p.vector[1] for p in self.points]
        z = [p.vector[2] for p in self.points]
        ax.plot(x, y, z, 'ok')

        # Drawing 3 times length of segment on each side
        u = self.point2 - self.points[0]
        x1, y1, z1 = (self.points[0] - 3 * u).vector
        x2, y2, z2 = (self.point2 + 3 * u).vector
        if dashed:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=color,
                    dashes=[30, 5, 10, 5])
        else:
            ax.plot([x1, x2], [y1, y2], [z1, z2], color=color)
        return ax

    def PlaneProjection2D(self, center, x, y):
        return Line2D(self.points[0].PlaneProjection2D(center, x, y),
                      self.point2.PlaneProjection2D(center, x, y))

    def MinimumDistancePoints(self, other_line):
        """
        Returns the points on this line and the other line that are the closest
        of lines
        """
        u = self.point2 - self.point1
        v = other_line.point2 - other_line.point1
        w = self.point1 - other_line.point1
        a = u.Dot(u)
        b = u.Dot(v)
        c = v.Dot(v)
        d = u.Dot(w)
        e = v.Dot(w)

        s = (b * e - c * d) / (a * c - b ** 2)
        t = (a * e - b * d) / (a * c - b ** 2)
        p1 = self.point1 + s * u
        p2 = other_line.point1 + t * v
        return p1, p2

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            return Line3D(*[p.Rotation(center, axis, angle, copy=True) for p in
                            self.points])
        else:
            for p in self.points:
                p.Rotation(center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Line3D(
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
                return Line3D(*[frame.OldCoordinates(p) for p in self.points])
            else:
                for p in self.points:
                    self.points = [frame.OldCoordinates(p) for p in
                                   self.points]
        if side == 'new':
            if copy:
                return Line3D(*[frame.NewCoordinates(p) for p in self.points])
            else:
                for p in self.points:
                    self.points = [frame.NewCoordinates(p) for p in
                                   self.points]

    def copy(self):
        return Line3D(*[p.copy() for p in self.points])

    @classmethod
    def from_step(cls, arguments, object_dict):
        point1 = object_dict[arguments[1]]
        direction = object_dict[arguments[2]]
        point2 = point1 + direction
        return cls(point1, point2, arguments[0][1:-1])

    def Intersection(self, line2):

        x1 = self.point1.vector[0]
        y1 = self.point1.vector[1]
        z1 = self.point1.vector[2]
        x2 = self.point2.vector[0]
        y2 = self.point2.vector[1]
        z2 = self.point2.vector[2]
        x3 = line2.points[0].vector[0]
        y3 = line2.points[0].vector[1]
        z3 = line2.points[0].vector[2]
        x4 = line2.point2.vector[0]
        y4 = line2.point2.vector[1]
        z4 = line2.point2.vector[2]

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
                return Point3D([x1 + (x2 - x1) * t1, y1 + (y2 - y1) * t1,
                                z1 + (z2 - z1) * t1])

        return None

class LineSegment3D(volmdlr.primitives.LineSegment):
    """
    Define a line segment limited by two points
    """

    def __init__(self, start, end, name=''):
        volmdlr.primitives.LineSegment.__init__(self, start, end, name)
        self.bounding_box = self._bounding_box()

    def __hash__(self):
        return 2 + hash(self.start) + hash(self.end)

    def __eq__(self, other_linesegment3d):
        if other_linesegment3d.__class__ != self.__class__:
            return False
        return

    def _bounding_box(self):
        points = [self.start, self.end]

        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])

        return volmdlr.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)


    def Length(self):
        return self.end.point_distance(self.start)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.start + curvilinear_abscissa * (
                    self.end - self.start) / self.Length()

    def frenet(self, curvilinear_abscissa):
        return self.unit_direction_vector(), None

    def middle_point(self):
        l = self.Length()
        return self.PointAtCurvilinearAbscissa(0.5 * l)

    def PlaneProjection2D(self, center, x, y):
        return LineSegment2D(self.start.PlaneProjection2D(center, x, y),
                             self.point2.PlaneProjection2D(center, x, y))

    def Intersection(self, segment2):
        x1 = self.start.vector[0]
        y1 = self.start.vector[1]
        z1 = self.start.vector[2]
        x2 = self.point2.vector[0]
        y2 = self.point2.vector[1]
        z2 = self.point2.vector[2]
        x3 = segment2.start.vector[0]
        y3 = segment2.start.vector[1]
        z3 = segment2.start.vector[2]
        x4 = segment2.end_point.vector[0]
        y4 = segment2.end_point.vector[1]
        z4 = segment2.end_point.vector[2]

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
                    return Point3D([x1 + (x2 - x1) * t1, y1 + (y2 - y1) * t1,
                                    z1 + (z2 - z1) * t1])

        return None

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            return LineSegment3D(
                *[p.Rotation(center, axis, angle, copy=True) for p in
                  self.points])
        else:
            Edge.Rotation(self, center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def __contains__(self, point):
        point1, point2 = self.start, self.end
        axis = Vector3D(point2 - point1)
        test = point.Rotation(point1, axis, math.pi)
        if test == point:
            return True
        else:
            return False

    def Translation(self, offset, copy=True):
        if copy:
            return LineSegment3D(
                *[p.Translation(offset, copy=True) for p in self.points])
        else:
            Edge.Translation(self, offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return LineSegment3D(
                    *[frame.OldCoordinates(p) for p in self.points])
            else:
                Edge.frame_mapping(self, frame, side, copy=False)
                self.bounding_box = self._bounding_box()
        if side == 'new':
            if copy:
                return LineSegment3D(
                    *[frame.NewCoordinates(p) for p in self.points])
            else:
                Edge.frame_mapping(self, frame, side, copy=False)
                self.bounding_box = self._bounding_box()

    def copy(self):
        return LineSegment3D(*[p.copy() for p in self.points])

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        points = [self.start, self.end]
        x = [p.vector[0] for p in points]
        y = [p.vector[1] for p in points]
        z = [p.vector[2] for p in points]
        ax.plot(x, y, z, 'o-k')
        return ax

    def MPLPlot2D(self, x_3D, y_3D, ax=None, color='k', width=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        edge2D = self.PlaneProjection2D(O3D, x_3D, y_3D)
        edge2D.MPLPlot(ax=ax, color=color, width=width)
        return ax

    def plot_data(self, x_3D, y_3D, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        edge2D = self.PlaneProjection2D(O3D, x_3D, y_3D)
        return edge2D.plot_data(marker, color, stroke_width,
                                dash, opacity, arrow)

    def FreeCADExport(self, name, ndigits=6):
        name = 'primitive' + str(name)
        x1, y1, z1 = round(1000 * self.start, ndigits).vector
        x2, y2, z2 = round(1000 * self.end, ndigits).vector
        return '{} = Part.LineSegment(fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(
            name, x1, y1, z1, x2, y2, z2)

    def to_line(self):
        return Line3D(*self.points)

    def babylon_script(self, color=(1, 1, 1), name='line', type_='line',
                       parent=None):
        if type_ == 'line' or type_ == 'dashed':
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
                s += 'var {} = BABYLON.MeshBuilder.CreateDashedLines("lines", {{points: myPoints, dashNb:20}}, scene);'.format(
                    name)
            s += '{}.color = new BABYLON.Color3{};\n'.format(name,
                                                             tuple(color))
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

    def To2D(self, plane_origin, x1, x2):
        p2D = [p.To2D(plane_origin, x1, x2) for p in (self.start, self.end)]
        return LineSegment2D(*p2D, name=self.name)

    def reverse(self):
        return LineSegment3D(self.end.copy(), self.start.copy())

    def MinimumDistancePoints(self, other_line):
        """
        Returns the points on this line and the other line that are the closest
        of lines
        """
        u = self.end - self.start
        v = other_line.end_point - other_line.start
        w = self.start - other_line.start
        a = u.Dot(u)
        b = u.Dot(v)
        c = v.Dot(v)
        d = u.Dot(w)
        e = v.Dot(w)
        if (a * c - b ** 2) != 0:
            s = (b * e - c * d) / (a * c - b ** 2)
            t = (a * e - b * d) / (a * c - b ** 2)
            p1 = self.start + s * u
            p2 = other_line.start + t * v
            return p1, p2, s, t
        else:
            return None, None, -1, -1

    def Matrix_distance(self, other_line):
        u = self.direction_vector()
        v = other_line.direction_vector()
        w = other_line.start - self.start

        a = u.Dot(u)
        b = -u.Dot(v)
        d = v.Dot(v)

        e = w.Dot(u)
        f = -w.Dot(v)

        A = npy.array([[a, b],
                       [b, d]])
        B = npy.array([e, f])

        res = scp.optimize.lsq_linear(A, B, bounds=(0, 1))
        p1 = self.PointAtCurvilinearAbscissa(res.x[0] * self.Length())
        p2 = other_line.PointAtCurvilinearAbscissa(
            res.x[1] * other_line.Length())
        return p1, p2

    def parallele_distance(self, other_linesegment):
        ptA, ptB, ptC = self.start, self.end, other_linesegment.points[0]
        u = Vector3D((ptA - ptB).vector)
        u.Normalize()
        plane1 = volmdlr.faces.Plane3D.from_3_points(ptA, ptB, ptC)
        v = u.Cross(plane1.normal)  # distance vector
        # ptA = k*u + c*v + ptC
        res = (ptA - ptC).vector
        x, y, z = res[0], res[1], res[2]
        u1, u2, u3 = u.vector[0], u.vector[1], u.vector[2]
        v1, v2, v3 = v.vector[0], v.vector[1], v.vector[2]

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
        if element.__class__ is Arc3D or element.__class__ is Circle3D:
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

        else:
            return NotImplementedError

    def extrusion(self, extrusion_vector):
        u = self.unit_direction_vector()
        v = extrusion_vector.copy()
        v.Normalize()
        w = u.Cross(v)
        l1 = self.Length()
        l2 = extrusion_vector.Norm()
        # outer_contour = Polygon2D([O2D, Point2D((l1, 0.)),
        #                            Point2D((l1, l2)), Point2D((0., l2))])
        plane = volmdlr.surfaces3d.Plane3D(Frame3D(self.start, u, v, w))
        return plane.rectangular_cut(0, l1, 0, l2)

    def revolution(self, axis_point, axis, angle):
        axis_line3d = Line3D(axis_point, axis_point + axis)

        p1_proj, _ = axis_line3d.point_projection(self.start)
        p2_proj, _ = axis_line3d.point_projection(self.end)
        d1 = self.start.point_distance(p1_proj)
        d2 = self.end.point_distance(p2_proj)
        u1 = (self.start - p1_proj)  # Unit vector from p1_proj to p1
        u1.Normalize()
        if u1.is_colinear_to(self.direction_vector()):
            # Planar face
            surface = volmdlr.faces.Plane3D(Frame3D(p1_proj, u1, axis.Cross(u1), axis))
            if angle == two_pi:
                # Only 2 circles as countours
                r, R = sorted([d1, d2])
                v1 = axis.Cross(u1)
                outer_contour3d = Circle3D(Frame3D(p1_proj, u1, v1, axis),
                                           R)
                if not math.isclose(r, 0, abs_tol=1e-9):
                    inner_contours3d = [Circle3D(Frame3D(p1_proj, u1, v1, axis)
                                                 , r)]
                else:
                    inner_contours3d = []
            else:
                # Two arcs and lines
                arc1_i = self.end.Rotation(axis=axis,
                                                 axis_point=axis_point,
                                                 angle=0.5 * angle)
                arc1_e = self.end.Rotation(axis=axis,
                                                 axis_point=axis_point,
                                                 angle=angle)
                arc2_i = self.start.Rotation(axis=axis,
                                                 axis_point=axis_point,
                                                 angle=0.5 * angle)
                arc2_s = self.start.Rotation(axis=axis,
                                                 axis_point=axis_point,
                                                 angle=angle)

                arc1 = Arc3D(self.end, arc1_i, arc1_e)
                arc2 = Arc3D(arc2_s, arc1_i, self.start)
                line2 = LineSegment3D(arc1_e, arc2_s)
                outer_contour3d = Contour3D([self, arc1, line2, arc2])
                inner_contours3d = []
            return volmdlr.faces.PlaneFace3D.from_contours3d(
                outer_contour3d=outer_contour3d,
                inner_contours3d=inner_contours3d)

        elif d1 != d2:
            # Conical
            v = axis.Cross(u1)
            w = axis.Cross(v)
            u1 = self.direction_vector()
            semi_angle = math.asin(u1.Cross(axis).Norm())
            surface = volmdlr.surfaces.ConicalSurface3D(Frame3D(p1_proj, axis, v, w),
                                       semi_angle)
            return surface.rectangular_cut(0, self.Length(), 0, angle)
        else:
            v = axis.Cross(u1)
            surface = volmdlr.surfaces.CylindricalSurface3D(Frame3D(p1_proj, u1, v, axis), d1)
            return surface.rectangular_cut(0, angle,
                                           0, (self.end-self.start).Dot(axis))


class BSplineCurve3D(volmdlr.Primitive3D):
    def __init__(self, degree, control_points, knot_multiplicities, knots,
                 weights=None, periodic=False, name=''):
        Primitive3D.__init__(self, basis_primitives=control_points, name=name)
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
            P = [(control_points[i][0], control_points[i][1],
                  control_points[i][2]) for i in range(len(control_points))]
            curve.ctrlpts = P
        else:
            Pw = [(control_points[i][0] * weights[i],
                   control_points[i][1] * weights[i],
                   control_points[i][2] * weights[i], weights[i]) for i in
                  range(len(control_points))]
            curve.ctrlptsw = Pw
        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot] * knot_multiplicities[i])
        curve.knotvector = knot_vector
        curve.delta = 0.1
        curve_points = curve.evalpts

        self.curve = curve
        self.points = [Point3D((p[0], p[1], p[2])) for p in curve_points]

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
                LineSegment3D(self.points[k], self.points[k + 1]))
        for primitive in primitives:
            primitive_length = primitive.Length()
            if length + primitive_length >= curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(
                    curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError

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
        closed_curve = False
        return cls(degree, points, knot_multiplicities, knots, weight_data,
                   closed_curve, name)

    def point_distance(self, pt1):
        distances = []
        for point in self.points:
            #            vmpt = Point3D((point[1], point[2], point[3]))
            distances.append(pt1.point_distance(point))
        return min(distances)

    def Rotation(self, center, axis, angle, copy=True):
        new_control_points = [p.Rotation(center, axis, angle, True) for p in
                              self.control_points]
        new_BSplineCurve3D = BSplineCurve3D(self.degree, new_control_points,
                                            self.knot_multiplicities,
                                            self.knots, self.weights,
                                            self.periodic, self.name)
        if copy:
            return new_BSplineCurve3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineCurve3D.curve
            self.points = new_BSplineCurve3D.points

    def Translation(self, offset, copy=True):
        new_control_points = [p.Translation(offset, True) for p in
                              self.control_points]
        new_BSplineCurve3D = BSplineCurve3D(self.degree, new_control_points,
                                            self.knot_multiplicities,
                                            self.knots, self.weights,
                                            self.periodic, self.name)
        if copy:
            return new_BSplineCurve3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineCurve3D.curve
            self.points = new_BSplineCurve3D.points

    # Copy paste du LineSegment3D
    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        x = [p.vector[0] for p in self.points]
        y = [p.vector[1] for p in self.points]
        z = [p.vector[2] for p in self.points]
        ax.plot(x, y, z, 'o-k')
        return ax

    def To2D(self, plane_origin, x1, x2):
        control_points2D = [p.To2D(plane_origin, x1, x2) for p in
                            self.control_points]
        return BSplineCurve2D(self.degree, control_points2D,
                              self.knot_multiplicities, self.knots,
                              self.weights, self.periodic, self.name)

    def tessellation_points(self):
        return self.points


class Arc3D(volmdlr.Primitive3D):
    """
    An arc is defined by a starting point, an end point and an interior point

    """

    def __init__(self, start, interior, end, name='',
                 other_vec=None):
        """

        """
        self.start = start
        self.interior = interior
        self.end = end
        self.other_vec = other_vec
        self.setup_arc(start, interior, end, name=name)

    @classmethod
    def from_angle(cls, start: volmdlr.Point3D, angle: float,
                   axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D):
        start_gen = start
        int_gen = start_gen.Rotation(axis_point, axis, angle / 2, copy=True)
        end_gen = start_gen.Rotation(axis_point, axis, angle, copy=True)
        if angle == two_pi:
            line = Line3D(axis_point, axis_point + axis)
            center, _ = line.point_projection(start)
            radius = center.point_distance(start)
            u = start - center
            v = axis.Cross(u)
            return Circle3D(Frame3D(center, u, v, axis), radius)
        return cls(start_gen, int_gen, end_gen, axis)

    def setup_arc(self, start, interior, end, name=''):
        u1 = (self.interior - self.start)
        u2 = (self.interior - self.end)
        try:
            u1.Normalize()
            u2.Normalize()
        except ZeroDivisionError:
            raise ValueError(
                'Start, end and interior points of an arc must be distincts')

        self.normal = u2.Cross(u1)
        self.normal.Normalize()

        if u1 == u2:
            u2 = self.normal.Cross(u1)
            u2.Normalize()

        v1 = self.normal.Cross(u1)  # v1 is normal, equal u2
        v2 = self.normal.Cross(u2)  # equal -u1

        p11 = 0.5 * (start + interior)  # Mid point of segment s,m
        p12 = p11 + v1
        p21 = 0.5 * (end + interior)  # Mid point of segment s,m
        p22 = p21 + v2

        l1 = Line3D(p11, p12)
        l2 = Line3D(p21, p22)

        try:
            c1, _ = l1.MinimumDistancePoints(l2)
        except ZeroDivisionError:
            raise ValueError(
                'Start, end and interior points  of an arc must be distincts')

        self.center = c1
        self.radius = (self.center - self.start).Norm()

        # Determining angle

        if self.other_vec is None:
            vec1 = (self.start - self.center)
        else:
            vec1 = self.other_vec
        vec1.Normalize()
        vec2 = self.normal.Cross(vec1)

        r1 = self.start.To2D(self.center, vec1, vec2)
        r2 = self.end.To2D(self.center, vec1, vec2)
        ri = self.interior.To2D(self.center, vec1, vec2)

        angle1 = math.atan2(r1.vector[1], r1.vector[0])
        anglei = math.atan2(ri.vector[1], ri.vector[0])
        angle2 = math.atan2(r2.vector[1], r2.vector[0])

        # Going trigo/clock wise from start to interior
        if anglei < angle1:
            trigowise_path = (anglei + two_pi) - angle1
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
            self.angle = trigowise_path
        else:
            # Clock wise
            self.is_trigo = False
            self.angle = clockwise_path

        if self.angle > math.pi:
            # Inverting normal to be sure to have a right defined normal for rotation
            self.normal = -self.normal
        volmdlr.Primitive3D.__init__(self, name=name)

    @property
    def points(self):
        return [self.start, self.interior, self.end]

    def tessellation_points(self, resolution_for_circle=40):
        number_points_tesselation = resolution_for_circle
        l = self.Length()
        tessellation_points_3D = [self.PointAtCurvilinearAbscissa(
            l * i / (number_points_tesselation)) for i in
                                  range(number_points_tesselation + 1)]
        return tessellation_points_3D

    def Length(self):
        return self.radius * abs(self.angle)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.start.Rotation(self.center, self.normal,
                                   curvilinear_abscissa / self.radius,
                                   copy=True)

    def frenet(self, curvilinear_abscissa):
        theta = curvilinear_abscissa/self.radius
        t0 = self.normal.Cross(self.start-self.center)
        t0.Normalize()
        tangent = t0.Rotation(self.center, self.normal, theta, copy=True)
        normal = -self.normal.Cross(tangent)
        return tangent, normal

    def Rotation(self, rot_center, axis, angle, copy=True):
        if copy:
            new_start = self.start.Rotation(rot_center, axis, angle, True)
            new_interior = self.interior.Rotation(rot_center, axis, angle,
                                                  True)
            new_end = self.end.Rotation(rot_center, axis, angle, True)
            return Arc3D(new_start, new_interior, new_end, name=self.name)
        else:
            self.center.Rotation(rot_center, axis, angle, False)
            self.start.Rotation(rot_center, axis, angle, False)
            self.interior.Rotation(rot_center, axis, angle, False)
            self.end.Rotation(rot_center, axis, angle, False)
            [p.Rotation(rot_center, axis, angle, False) for p in
             self.primitives]

    def Translation(self, offset, copy=True):
        if copy:
            new_start = self.start.Translation(offset, True)
            new_interior = self.interior.Translation(offset, True)
            new_end = self.end.Translation(offset, True)
            return Arc3D(new_start, new_interior, new_end, name=self.name)
        else:
            self.center.Translation(offset, False)
            self.start.Translation(offset, False)
            self.interior.Translation(offset, False)
            self.end.Translation(offset, False)
            [p.Translation(offset, False) for p in self.primitives]

    def MPLPlot(self, ax=None):
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
        for px, py, pz in self.tessellation_points():
            x.append(px)
            y.append(py)
            z.append(pz)

        ax.plot(x, y, z, 'k')
        return ax

    def MPLPlot2D(self, center=volmdlr.O3D,
                  x3d=volmdlr.X3D, y3D=volmdlr.Y3D,
                  ax=None, color='k'):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        # TODO: Enhance this plot
        l = self.Length()
        x = []
        y = []
        for i in range(30):
            p = self.PointAtCurvilinearAbscissa(i / (29.) * l)
            xi, yi = p.PlaneProjection2D(center, X3D, Y3D)
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

    def copy(self):
        return Arc3D(self.start.copy(), self.interior.copy(), self.end.copy())

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            new_start = frame.OldCoordinates(self.start.copy())
            new_interior = frame.OldCoordinates(self.interior.copy())
            new_end = frame.OldCoordinates(self.end.copy())
            if copy:
                return Arc3D(new_start, new_interior, new_end, normal=None,
                             name=self.name)
            else:
                self.start, self.interior, self.end = new_start, new_interior, new_end
                self.setup_arc(self.start, self.interior, self.end)

        if side == 'new':
            new_start = frame.NewCoordinates(self.start.copy())
            new_interior = frame.NewCoordinates(self.interior.copy())
            new_end = frame.NewCoordinates(self.end.copy())
            if copy:
                return Arc3D(new_start, new_interior, new_end, normal=None,
                             name=self.name)
            else:
                self.start, self.interior, self.end = new_start, new_interior, new_end
                self.setup_arc(self.start, self.interior, self.end)

    def To2D(self, plane_origin, x, y):
        ps = self.start.To2D(plane_origin, x, y)
        pi = self.interior.To2D(plane_origin, x, y)
        pe = self.end.To2D(plane_origin, x, y)
        return volmdlr.primitives2D.Arc2D(ps, pi, pe, name=self.name)

    def minimum_distance_points_arc(self, other_arc):

        u1 = self.start - self.center
        u1.Normalize()
        u2 = self.normal.Cross(u1)

        w = other_arc.center - self.center

        u3 = other_arc.start - other_arc.center
        u3.Normalize()
        u4 = other_arc.normal.Cross(u3)

        r1, r2 = self.radius, other_arc.radius

        a, b, c, d = u1.Dot(u1), u1.Dot(u2), u1.Dot(u3), u1.Dot(u4)
        e, f, g = u2.Dot(u2), u2.Dot(u3), u2.Dot(u4)
        h, i = u3.Dot(u3), u3.Dot(u4)
        j = u4.Dot(u4)
        k, l, m, n, o = w.Dot(u1), w.Dot(u2), w.Dot(u3), w.Dot(u4), w.Dot(w)

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

        p1 = self.PointAtCurvilinearAbscissa(res1.x[0] * r1)
        p2 = other_arc.PointAtCurvilinearAbscissa(res1.x[1] * r2)

        return p1, p2

    def minimum_distance_points_line(self, other_line):

        u = other_line.DirectionVector()
        k = self.start - self.center
        k.Normalize()
        w = self.center - other_line.points[0]
        v = self.normal.Cross(k)

        r = self.radius

        a = u.Dot(u)
        b = u.Dot(v)
        c = u.Dot(k)
        d = v.Dot(v)
        e = v.Dot(k)
        f = k.Dot(k)
        g = w.Dot(u)
        h = w.Dot(v)
        i = w.Dot(k)
        j = w.Dot(w)

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

        p1 = other_line.PointAtCurvilinearAbscissa(
            res1.x[0] * other_line.Length())
        p2 = self.PointAtCurvilinearAbscissa(res1.x[1] * r)

        res = [res2, res3]
        for couple in res:
            ptest1 = other_line.PointAtCurvilinearAbscissa(
                couple.x[0] * other_line.Length())
            ptest2 = self.PointAtCurvilinearAbscissa(couple.x[1] * r)
            dtest = ptest1.point_distance(ptest2)
            if dtest < d:
                p1, p2 = ptest1, ptest2

        return p1, p2

    def minimum_distance(self, element, return_points=False):
        if element.__class__ is Arc3D or element.__class__ is Circle3D:
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
            u.Normalize()
            w = extrusion_vector.copy()
            w.Normalize()
            v = w.Cross(u)
            arc2d = self.To2D(self.center, u, v)
            cylinder = volmdlr.surfaces3d.CylindricalSurface3D(volmdlr.Frame3D(self.center,
                                                    u,
                                                    v,
                                                    w),
                                            self.radius
                                            )
            return cylinder.rectangular_cut(arc2d.angle1,
                                            arc2d.angle2,
                                            0, extrusion_vector.Norm())
        else:
            raise NotImplementedError('Elliptic faces not handled: dot={}'.format(
                self.normal.Dot(extrusion_vector)
            ))


    def revolution(self, axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                   angle: float):
        line3d = Line3D(axis_point, axis_point + axis)
        tore_center, _ = line3d.point_projection(self.center)
        u =  self.center - tore_center
        u.Normalize()
        v = axis.Cross(u)
        if not math.isclose(self.normal.Dot(u), 0., abs_tol=1e-9):
            raise NotImplementedError(
                'Outside of plane revolution not supported')

        R = tore_center.point_distance(self.center)
        surface = surfaces3d.ToroidalSurface3D(Frame3D(tore_center, u, v, axis), R,
                                    self.radius)
        arc2d = self.To2D(tore_center, u, axis)
        return surface.rectangular_cut(0, angle,
                                       arc2d.angle1, arc2d.angle2)



class ArcEllipse3D(volmdlr.Primitive3D):
    """
    An arc is defined by a starting point, an end point and an interior point

    """

    def __init__(self, start, interior, end, center, major_dir, normal=None,
                 name='', extra=None):
        # Extra is an additionnal point if start=end because you need 3 points on the arcellipse to define it
        self.start = start
        self.interior = interior
        self.end = end
        self.center = center
        major_dir.Normalize()
        self.major_dir = major_dir  # Vector for Gradius
        self.extra = extra

        u1 = (self.interior - self.start)
        u2 = (self.interior - self.end)
        u1.Normalize()
        u2.Normalize()

        if u1 == u2:
            u2 = (self.interior - self.extra)
            u2.Normalize()

        if normal is None:
            n = u2.Cross(u1)
            n.Normalize()
            self.normal = n
        else:
            n = normal
            n.Normalize()
            self.normal = normal

        self.minor_dir = self.normal.Cross(self.major_dir)

        frame = Frame3D(self.center, self.major_dir, self.minor_dir,
                        self.normal)
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
            theta = 0.5 * math.atan(2 * C[2] / (C[1] - C[0]))
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

        if self.start == self.end:
            self.angle = two_pi

        if self.is_trigo:
            self.offset_angle = angle1
        else:
            self.offset_angle = angle2

        Primitive3D.__init__(self, basis_primitives=self.tessellation_points(),
                             name=name)

    def _get_points(self):
        return self.tessellation_points()

    points = property(_get_points)

    def tessellation_points(self, resolution_for_ellipse=40):
        number_points_tesselation = math.ceil(
            resolution_for_ellipse * abs(0.5 * self.angle / math.pi))

        plane3d = volmdlr.faces.Plane3D(self.center, self.major_dir, self.minor_dir,
                          self.normal)
        frame3d = Frame3D(self.center, plane3d.vectors[0], plane3d.vectors[1],
                          plane3d.normal)

        tessellation_points_3D = [Point3D((self.Gradius * math.cos(
            self.offset_angle + self.angle * i / (number_points_tesselation)),
                                           self.Sradius * math.sin(
                                               self.offset_angle + self.angle * i / (
                                                   number_points_tesselation)),
                                           0)) for i in
                                  range(number_points_tesselation + 1)]

        global_points = []
        for pt in tessellation_points_3D:
            global_points.append(frame3d.OldCoordinates(pt))

        return global_points

    def To2D(self, plane_origin, x, y):
        ps = self.start.To2D(plane_origin, x, y)
        pi = self.interior.To2D(plane_origin, x, y)
        pe = self.end.To2D(plane_origin, x, y)
        center = self.center.To2D(plane_origin, x, y)

        if self.extra is None:
            pextra = None
        else:
            pextra = self.extra.To2D(plane_origin, x, y)

        maj_dir2d = self.major_dir.To2D(plane_origin, x, y)
        maj_dir2d.Normalize()
        return ArcEllipse2D(ps, pi, pe, center, maj_dir2d, name=self.name,
                            extra=pextra)

    def Length(self):
        return self.angle * math.sqrt(
            (self.Gradius ** 2 + self.Sradius ** 2) / 2)

    def MPLPlot(self, ax=None):
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
        for px, py, pz in self.tessellation_points():
            x.append(px)
            y.append(py)
            z.append(pz)

        ax.plot(x, y, z, 'k')
        return ax

    def MPLPlot2D(self, x3d, y3D, ax, color='k'):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        # TODO: Enhance this plot
        l = self.Length()
        x = []
        y = []
        for i in range(30):
            p = self.PointAtCurvilinearAbscissa(i / (29.) * l)
            xi, yi = p.PlaneProjection2D(X3D, Y3D)
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


class OpenedRoundedLineSegments3D(volmdlr.Wire3D, volmdlr.primitives.RoundedLineSegments):
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    
    def __init__(self, points, radius, adapt_radius=False, name=''):
        primitives = RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.LineSegment3D,
                                                  volmdlr.Arc3D,
                                                  closed=False,
                                                  adapt_radius=adapt_radius,
                                                  name='')

        volmdlr.Wire3D.__init__(self, primitives, name)

    def ArcFeatures(self, ipoint):
        radius = self.radius[ipoint]
        if self.closed:
            if ipoint == 0:
                pt1 = self.points[-1]
            else:
                pt1 = self.points[ipoint -1]
            pti = self.points[ipoint]
            if ipoint < self.npoints-1:
                pt2 = self.points[ipoint+1]
            else:
                pt2 = self.points[0]
        else:
            pt1 = self.points[ipoint - 1]
            pti = self.points[ipoint]
            pt2 = self.points[ipoint + 1]

        dist1 = (pt1 - pti).Norm()
        dist2 = (pt2 - pti).Norm()
        dist3 = (pt1 - pt2).Norm()
        alpha = math.acos(-(dist3**2-dist1**2-dist2**2)/(2*dist1*dist2))/2.
        dist = radius/math.tan(alpha)

        u1 = (pt1 - pti) / dist1
        u2 = (pt2 - pti) / dist2

        p3 = pti + u1*dist
        p4 = pti + u2*dist

        n = u1.Cross(u2)
        n /= n.Norm()
        v1 = u1.Cross(n)
        v2 = u2.Cross(n)

        l1 = volmdlr.Line3D(p3, p3+v1)
        l2 = volmdlr.Line3D(p4, p4+v2)
        c, _ = l1.MinimumDistancePoints(l2)

        u3 = u1 + u2# mean of v1 and v2
        u3 /= u3.Norm()

        interior = c - u3 * radius
        return p3, interior, p4, dist, alpha


    def Rotation(self, center, angle, copy=True):
        if copy:
            return self.__class__([p.Rotation(center, angle, copy=True)\
                                          for p in self.points],
                                         self.radius, self.closed, self.name)
        else:
            self.__init__([p.Rotation(center, angle, copy=True)\
                           for p in self.points],
                          self.radius, self.closed, self.name)

    def Translation(self, offset, copy=True):
        if copy:
            return self.__class__([p.Translation(offset, copy=True)\
                                          for p in self.points],
                                         self.radius, self.closed, self.name)
        else:
            self.__init__([p.Translation(offset, copy=True)\
                           for p in self.points],
                          self.radius, self.closed, self.name)


class ClosedRoundedLineSegments3D(volmdlr.Contour3D, OpenedRoundedLineSegments3D):
    """
    :param points: Points used to draw the wire 
    :type points: List of Point3D
    :param radius: Radius used to connect different parts of the wire
    :type radius: {position1(n): float which is the radius linked the n-1 and n+1 points, position2(n+1):...}
    """
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    
    def __init__(self, points, radius, adapt_radius=False, name=''):
        primitives = RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.LineSegment3D,
                                                  volmdlr.Arc3D,
                                                  closed=True,
                                                  adapt_radius=adapt_radius,
                                                  name='')

        volmdlr.Contour3D.__init__(self, primitives, name=name)



class Block(volmdlr.Shell3D):
    _standalone_in_db = True
    _generic_eq = True
    _non_serializable_attributes  = ['size']
    _non_eq_attributes = ['name', 'color', 'alpha', 'size', 'bounding_box', 'faces', 'contours',
                          'plane', 'points', 'polygon2D']
    _non_hash_attributes = []

    """
    Creates a block
    :param frame: a frame 3D. The origin of the frame is the center of the block,
     the 3 vectors are defining the edges. The frame has not to be orthogonal
    """
    def __init__(self, frame:volmdlr.Frame3D, *,
                 color:Tuple[float, float, float]=None, alpha:float=1.,
                 name:str=''):
        self.frame = frame
        self.size = (self.frame.u.Norm(), self.frame.v.Norm(), self.frame.w.Norm())

        faces = self.shell_faces()
        volmdlr.Shell3D.__init__(self, faces,  color=color, alpha=alpha, name=name)

    # def __hash__(self):
    #     return hash(self.frame)

    def Vertices(self):
        return [self.frame.origin - 0.5*self.frame.u - 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin - 0.5*self.frame.u + 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u + 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u - 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin - 0.5*self.frame.u - 0.5*self.frame.v + 0.5*self.frame.w,
                self.frame.origin - 0.5*self.frame.u + 0.5*self.frame.v + 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u + 0.5*self.frame.v + 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u - 0.5*self.frame.v + 0.5*self.frame.w]

    def Edges(self):
        p1, p2, p3, p4, p5, p6, p7, p8 = self.Vertices()
        return [volmdlr.LineSegment3D(p1.copy(), p2.copy()),
                volmdlr.LineSegment3D(p2.copy(), p3.copy()),
                volmdlr.LineSegment3D(p3.copy(), p4.copy()),
                volmdlr.LineSegment3D(p4.copy(), p1.copy()),
                volmdlr.LineSegment3D(p5.copy(), p6.copy()),
                volmdlr.LineSegment3D(p6.copy(), p7.copy()),
                volmdlr.LineSegment3D(p7.copy(), p8.copy()),
                volmdlr.LineSegment3D(p8.copy(), p5.copy()),
                volmdlr.LineSegment3D(p1.copy(), p5.copy()),
                volmdlr.LineSegment3D(p2.copy(), p6.copy()),
                volmdlr.LineSegment3D(p3.copy(), p7.copy()),
                volmdlr.LineSegment3D(p4.copy(), p8.copy())]

    def face_contours(self):
        e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12 = self.Edges()
        e5_switch = e5.reverse()
        e6_switch = e6.reverse()
        e7_switch = e7.reverse()
        e8_switch = e8.reverse()
        e9_switch = e9.reverse()
        e10_switch = e10.reverse()
        e11_switch = e11.reverse()
        e12_switch = e12.reverse()
        return [volmdlr.Contour3D([e1.copy(), e2.copy(), e3.copy(), e4.copy()]),
                volmdlr.Contour3D([e5.copy(), e6.copy(), e7.copy(), e8.copy()]),
                # volmdlr.Contour3D([e1.copy(), e9.copy(), e5.copy(), e10.copy()]),
                volmdlr.Contour3D([e1.copy(), e10.copy(), e5_switch.copy(), e9_switch.copy()]),
                # volmdlr.Contour3D([e2.copy(), e10.copy(), e6.copy(), e11.copy()]),
                volmdlr.Contour3D([e2.copy(), e11.copy(), e6_switch.copy(), e10_switch.copy()]),
                # volmdlr.Contour3D([e3.copy(), e11.copy(), e7.copy(), e12.copy()]),
                volmdlr.Contour3D([e3.copy(), e12.copy(), e7_switch.copy(), e11_switch.copy()]),
                # volmdlr.Contour3D([e4.copy(), e12.copy(), e8.copy(), e9.copy()])]
                volmdlr.Contour3D([e4.copy(), e9.copy(), e8_switch.copy(), e12_switch.copy()])]

    def shell_faces(self):
        c1, c2, c3, c4, c5, c6 = self.face_contours()
        return [volmdlr.PlaneFace3D.from_contours3d(c1),
                volmdlr.PlaneFace3D.from_contours3d(c2),
                volmdlr.PlaneFace3D.from_contours3d(c3),
                volmdlr.PlaneFace3D.from_contours3d(c4),
                volmdlr.PlaneFace3D.from_contours3d(c5),
                volmdlr.PlaneFace3D.from_contours3d(c6)]

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_frame = self.frame.Rotation(center, axis, angle, copy=True)
            return Block(new_frame, color=self.color, alpha=self.alpha, name=self.name)
        else:
            self.frame.Rotation(center, axis, angle, copy=False)
            volmdlr.Shell3D.Rotation(self, center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            new_frame = self.frame.Translation(offset, copy=True)
            return Block(new_frame, color=self.color, alpha=self.alpha, name=self.name)
        else:
            self.frame.Translation(offset, copy=False)
            volmdlr.Shell3D.Translation(self, offset, copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.frame.origin)
            new_u = basis.NewCoordinates(self.frame.u)
            new_v = basis.NewCoordinates(self.frame.v)
            new_w = basis.NewCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return Block(new_frame, color=self.color, alpha=self.alpha, name=self.name)
            else:
                self.frame = new_frame
                volmdlr.Shell3D.frame_mapping(self, frame, side, copy=False)

        if side == 'old':
            new_origin = frame.OldCoordinates(self.frame.origin)
            new_u = basis.OldCoordinates(self.frame.u)
            new_v = basis.OldCoordinates(self.frame.v)
            new_w = basis.OldCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return Block(new_frame, color=self.color, alpha=self.alpha, name=self.name)
            else:
                self.frame = new_frame
                volmdlr.Shell3D.frame_mapping(self, frame, side, copy=False)

    def copy(self):
        new_origin = self.frame.origin.Copy()
        new_u = self.frame.u.Copy()
        new_v = self.frame.v.Copy()
        new_w = self.frame.w.Copy()
        new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
        return Block(new_frame, color=self.color, alpha=self.alpha, name=self.name)

    def plot_data(self, x3D, y3D, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        lines = []
        for edge3D in self.Edges():
            lines.append(edge3D.plot_data(x3D, y3D, marker, color, stroke_width,
                         dash, opacity, arrow))

        return lines

    def MPLPlot2D(self, x3D, y3D, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = None

        for edge3D in self.Edges():
#            edge2D = edge3D.PlaneProjection2D()
            edge3D.MPLPlot2D(x3D, y3D, ax)

        return fig, ax



class Cone(volmdlr.Primitive3D):
    def __init__(self, position, axis, radius, length, name=''):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.position = position
        axis.Normalize()
        self.axis = axis
        self.radius = radius
        self.length = length
        self.bounding_box = self._bounding_box()

    def _bounding_box(self):
        """
        A is the point at the basis
        B is the top
        """
        pointA = self.position - self.length/2 * self.axis
        pointB = self.position + self.length/2 * self.axis

        dx2 = (pointA[0]-pointB[0])**2
        dy2 = (pointA[1]-pointB[1])**2
        dz2 = (pointA[2]-pointB[2])**2

        kx = ((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5
        ky = ((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5
        kz = ((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5

        x_bound = (pointA[0] - kx * self.radius, pointA[0] + kx * self.radius, pointB[0])
        xmin = min(x_bound)
        xmax = max(x_bound)

        y_bound = (pointA[1] - ky * self.radius, pointA[1] + ky * self.radius, pointB[1])
        ymin = min(y_bound)
        ymax = max(y_bound)

        z_bound = (pointA[2] - kz * self.radius, pointA[2] + kz * self.radius, pointB[2])
        zmin = min(z_bound)
        zmax = max(z_bound)

        return volmdlr.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def Volume(self):
        return self.length * math.pi * self.radius**2 / 3

    def babylon_script(self):
        new_axis = volmdlr.Vector3D((self.axis[0], self.axis[1], self.axis[2]))
        normal_vector1 = new_axis.RandomUnitNormalVector()
        normal_vector2 = new_axis.Cross(normal_vector1)
        x, y, z = self.position
        s = 'var cone = BABYLON.MeshBuilder.CreateCylinder("cone", {{diameterTop:0, diameterBottom:{}, height: {}, tessellation: 100}}, scene);\n'.format(2*self.radius, self.length)
        s += 'cone.position = new BABYLON.Vector3({},{},{});\n;'.format(x,y,z)
        s += 'var axis1 = new BABYLON.Vector3({},{},{});\n'.format(new_axis[0], new_axis[1], new_axis[2])
        s += 'var axis2 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector1[0], normal_vector1[1], normal_vector1[2])
        s += 'var axis3 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector2[0], normal_vector2[1], normal_vector2[2])
        s += 'cone.rotation = BABYLON.Vector3.RotationFromAxis(axis3, axis1, axis2);\n'
        return s


class ExtrudedProfile(volmdlr.Shell3D):
    """
    
    """
    _non_serializable_attributes  = ['faces', 'inner_contours3d', 'outer_contour3d']
    def __init__(self, plane_origin, x, y, outer_contour2d, inner_contours2d,
                 extrusion_vector, color=None, alpha=1., name=''):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.plane_origin = plane_origin
        
        self.outer_contour2d = outer_contour2d
        self.outer_contour3d = outer_contour2d.To3D(plane_origin, x, y)
        
        self.inner_contours2d = inner_contours2d
        self.extrusion_vector = extrusion_vector
        self.inner_contours3d = []
        self.x = x
        self.y = y
        self.color = color

        bool_areas = []
        for contour in inner_contours2d:
            self.inner_contours3d.append(contour.To3D(plane_origin, x, y))
            if contour.Area() > outer_contour2d.Area():
                bool_areas.append(True)
            else:
                bool_areas.append(False)
        if any(bool_areas):
            raise ValueError('At least one inner contour is not contained in outer_contour.')

        faces = self.shell_faces()
        volmdlr.Shell3D.__init__(self, faces, color=color, alpha=alpha, name=name)

    def shell_faces(self):
        lower_face = volmdlr.surfaces3d.PlaneFace3D.from_contours3d(self.outer_contour3d,
                                                         self.inner_contours3d)

        upper_face = lower_face.Translation(self.extrusion_vector)
        lateral_faces = [p.extrusion(self.extrusion_vector)
                         for p in self.outer_contour3d.primitives]

        for inner_contour in self.inner_contours3d:
            lateral_faces.extend([p.extrusion(self.extrusion_vector)
                                  for p in inner_contour.primitives])

        return [lower_face]+[upper_face]+lateral_faces

    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        for contour in [self.outer_contour2d]+self.inner_contours2d:
            for primitive in contour.primitives:
                primitive.MPLPlot(ax)
        ax.margins(0.1)
        return ax

    def FreeCADExport(self, ip):
        name='primitive'+str(ip)
        s = 'Wo = []\n'
        s += 'Eo = []\n'
        for ip, primitive in enumerate(self.outer_contour3d.edges):
            s += primitive.FreeCADExport('L{}'.format(ip))
            s += 'Eo.append(Part.Edge(L{}))\n'.format(ip)
        s += 'Wo.append(Part.Wire(Eo[:]))\n'
        s += 'Fo = Part.Face(Wo)\n'

        s += 'Fi = []\n'
        s += 'W = []\n'
        for ic,contour in enumerate(self.inner_contours3d):
            s+='E = []\n'
            for ip, primitive in enumerate(contour.edges):
                s += primitive.FreeCADExport('L{}_{}'.format(ic, ip))
                s += 'E.append(Part.Edge(L{}_{}))\n'.format(ic, ip)
            s += 'Wi = Part.Wire(E[:])\n'
            s += 'Fi.append(Part.Face(Wi))\n'

        if len(self.inner_contours3d) != 0:
            s += 'Fo = Fo.cut(Fi)\n'
        e1, e2, e3 = round(1000*self.extrusion_vector, 6)

        s+='{} = Fo.extrude(fc.Vector({}, {}, {}))\n'.format(name,e1,e2,e3)
        return s

    def Area(self):
        areas = self.outer_contour2d.Area()
        areas -= sum([c.Area() for c in self.inner_contours2d])
        # sic=list(npy.argsort(areas))[::-1]# sorted indices of contours
        # area=areas[sic[0]]

        # for i in sic[1:]:
        #     area-=self.contours2D[i].Area()
        return areas

    def Volume(self):
        z = self.x.Cross(self.y)
        z.Normalize()
        coeff = npy.dot(self.extrusion_vector, z)
        return self.Area()*coeff


    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'old':
            extrusion_vector = basis.OldCoordinates(self.extrusion_vector)
            x = basis.OldCoordinates(self.x)
            y = basis.OldCoordinates(self.y)
        elif side == 'new':
            extrusion_vector = basis.NewCoordinates(self.extrusion_vector)
            x = basis.NewCoordinates(self.x)
            y = basis.NewCoordinates(self.y)
        else:
            raise ValueError('side must be either old or new')

        if copy:
            return ExtrudedProfile(self.plane_origin.frame_mapping(frame, side, copy),
                                   x,
                                   y,
                                   self.outer_contour2d,
                                   self.inner_contours2d,
                                   extrusion_vector)
        else:
            self.__init__(self.plane_origin.frame_mapping(frame, side, copy),
                          x,
                          y,
                          self.outer_contour2d,
                          self.inner_contours2d,
                          extrusion_vector)


class RevolvedProfile(volmdlr.Shell3D):
    """

    """
    _non_serializable_attributes  = ['faces', 'contour3D']

    def __init__(self, plane_origin, x, y, contour2d, axis_point,
                 axis, angle=2*math.pi, *, color=None, alpha=1, name=''):

        self.contour2d = contour2d
        self.axis_point = axis_point
        self.axis = axis
        self.angle = angle
        self.plane_origin = plane_origin
        self.x = x
        self.y = y
        self.contour3d = self.contour2d.To3D(plane_origin, x, y)

        faces = self.shell_faces()
        volmdlr.Shell3D.__init__(self, faces, color=color,
                                 alpha=alpha, name=name)

    def shell_faces(self):
        faces = []
                        
        for edge in self.contour3d.primitives:
            faces.append(edge.revolution(self.axis_point,
                                         self.axis, self.angle))
                    
        return faces

    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        for contour in self.contours3D:
            contour.MPLPlot(ax)

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive'+str(ip)
        s = 'W=[]\n'
#        for ic, contour in enumerate(self.contours3D):
        s += 'L=[]\n'
        for ibp, basis_primitive in enumerate(self.contour3D.edges):
            s += basis_primitive.FreeCADExport('L{}_{}'.format(1, ibp), 8)
            s += 'L.append(L{}_{})\n'.format(1,ibp)
        s += 'S = Part.Shape(L)\n'
        s += 'W.append(Part.Wire(S.Edges))\n'
        s += 'F=Part.Face(W)\n'
        a1, a2, a3 = self.axis.vector
        ap1, ap2, ap3 = self.axis_point.vector
        ap1 = round(ap1*1000, ndigits)
        ap2 = round(ap2*1000, ndigits)
        ap3 = round(ap3*1000, ndigits)
        angle = self.angle/math.pi*180
        s += '{} = F.revolve(fc.Vector({},{},{}), fc.Vector({},{},{}),{})\n'.format(name, ap1,ap2,ap3,a1,a2,a3,angle)
        return s

    def Volume(self):
        p1 = self.axis_point.PlaneProjection3D(self.plane_origin,self.x,self.y)
        p1_2D = p1.To2D(self.axis_point,self.x,self.y)
        p2_3D = self.axis_point+volmdlr.Point3D(self.axis.vector)
        p2_2D = p2_3D.To2D(self.plane_origin,self.x,self.y)
        axis_2D = volmdlr.Line2D(p1_2D,p2_2D)
        com = self.contour2D.CenterOfMass()
        if com is not False:
            rg = axis_2D.point_distance(com)
            return self.angle*rg*self.contour2D.Area()
        else:
            return 0

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'old':
            axis = basis.OldCoordinates(self.axis)
            x = basis.OldCoordinates(self.x)
            y = basis.OldCoordinates(self.y)
        elif side == 'new':
            axis = basis.NewCoordinates(self.axis)
            x = basis.NewCoordinates(self.x)
            y = basis.NewCoordinates(self.y)
        else:
            raise ValueError('side must be either old or new')

        if copy:

            return RevolvedProfile(self.plane_origin.frame_mapping(frame, side, copy),
                                   x,
                                   y,
                                   self.contour2D,
                                   self.axis_point.frame_mapping(frame, side, copy),
                                   axis,
                                   self.angle)
        else:
            self.__init__(self.plane_origin.frame_mapping(frame, side, copy),
                          x,
                          y,
                          self.contour2D,
                          self.axis_point.frame_mapping(frame, side, copy),
                          axis,
                          self.angle)

class Cylinder(RevolvedProfile):
    """
    Creates a full cylinder with the position, the axis of revolution, the radius and the length.
    """
    def __init__(self, position, axis, radius, length, color=None, alpha=1., name=''):
        self.position = position
        axis.Normalize()
        self.axis = axis
        self.radius = radius
        self.length = length
        self.bounding_box = self._bounding_box()

        # Revolved Profile
        p1 = volmdlr.Point2D((-0.5*self.length, 0))
        p2 = volmdlr.Point2D((0.5*self.length, 0))
        p3 = volmdlr.Point2D((0.5*self.length, self.radius))
        p4 = volmdlr.Point2D((-0.5*self.length, self.radius))
        l1 = volmdlr.LineSegment2D(p1, p2)
        l2 = volmdlr.LineSegment2D(p2, p3)
        l3 = volmdlr.LineSegment2D(p3, p4)
        l4 = volmdlr.LineSegment2D(p4, p1)
        contour = volmdlr.Contour2D([l1, l2, l3, l4])
        y = axis.RandomUnitNormalVector()
        RevolvedProfile.__init__(self, position, axis, y, contour, position, axis,
                                 color=color, alpha=alpha, name=name)


    def _bounding_box(self):
        
        if hasattr(self, 'radius'):
            radius = self.radius
        elif hasattr(self, 'outer_radius'):
            radius = self.outer_radius

            
        pointA = self.position - self.length/2 * self.axis
        pointB = self.position + self.length/2 * self.axis

        dx2 = (pointA[0]-pointB[0])**2
        dy2 = (pointA[1]-pointB[1])**2
        dz2 = (pointA[2]-pointB[2])**2

        kx = ((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5
        ky = ((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5
        kz = ((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5

        if pointA[0] > pointB[0]:
            pointA, pointB = pointB, pointA
        xmin = pointA[0] - kx * radius
        xmax = pointB[0] + kx * radius

        if pointA[1] > pointB[1]:
            pointA, pointB = pointB, pointA
        ymin = pointA[1] - ky * radius
        ymax = pointB[1] + ky * radius

        if pointA[2] > pointB[2]:
            pointA, pointB = pointB, pointA
        zmin = pointA[2] - kz * radius
        zmax = pointB[2] + kz * radius

        return volmdlr.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def Volume(self):
        return self.length * math.pi * self.radius**2

    def FreeCADExport(self, ip):
        if self.radius > 0:
            name = 'primitive'+str(ip)
            e = str(1000*self.length)
            r = str(1000*self.radius)
            position = 1000*(self.position - self.axis*self.length/2.)
            x, y, z = position
            x = str(x)
            y = str(y)
            z = str(z)

            ax, ay, az = self.axis
            ax = str(ax)
            ay = str(ay)
            az = str(az)
            return name+'=Part.makeCylinder('+r+','+e+',fc.Vector('+x+','+y+','+z+'),fc.Vector('+ax+','+ay+','+az+'),360)\n'
        else:
            return ''

    def babylon_script(self, name='primitive_mesh'):
        normal_vector1 = self.axis.RandomUnitNormalVector()
#        normal_vector2 = new_axis.Cross(normal_vector1)
#        x, y, z = self.position
#        s='var {} = BABYLON.Mesh.CreateCylinder("{}", {}, {}, {}, 30, 1, scene,false, BABYLON.Mesh.DEFAULTSIDE);'.format(name, self.name,self.length,2*self.radius,2*self.radius)
#        s += '{}.position = new BABYLON.Vector3({},{},{});\n;'.format(name, x,y,z)
#        s += 'var axis1 = new BABYLON.Vector3({},{},{});\n'.format(new_axis[0], new_axis[1], new_axis[2])
#        s += 'var axis2 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector1[0], normal_vector1[1], normal_vector1[2])
#        s += 'var axis3 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector2[0], normal_vector2[1], normal_vector2[2])
#        s += '{}.rotation = BABYLON.Vector3.RotationFromAxis(axis3, axis1, axis2);\n'.format(name)
        p1 = volmdlr.Point2D((-0.5*self.length, self.radius))
        p2 = volmdlr.Point2D((0.5*self.length, self.radius))
        p3 = volmdlr.Point2D((0.5*self.length, 0.))
        p4 = volmdlr.Point2D((-0.5*self.length, 0.))
        l1 = volmdlr.LineSegment2D(p1, p2)
        l2 = volmdlr.LineSegment2D(p2, p3)
        l3 = volmdlr.LineSegment2D(p3, p4)
        l4 = volmdlr.LineSegment2D(p4, p1)
        extruded_profile = RevolvedProfile(self.position, self.axis, normal_vector1,
                                           volmdlr.Contour2D([l1, l2, l3, l4]), self.position, self.axis, name=self.name)
        return extruded_profile.babylon_script(name=name)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'old':
            axis = basis.OldCoordinates(self.axis)
        elif side == 'new':
            axis = basis.NewCoordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')

        if copy:
            return Cylinder(self.position.frame_mapping(frame, side, copy),
                            axis,
                            self.radius, self.length, color=self.color,
                            alpha=self.alpha)
        else:
            self.position.frame_mapping(frame, side, copy)
            self.axis = axis
            Cylinder.__init__(self, self.position, self.axis, self.radius, 
                              self.length, color=self.color, alpha=self.alpha)
            
    def copy(self):
        new_position = self.position.copy()
        new_axis = self.axis.copy()
        return Cylinder(new_position, new_axis, self.radius, self.length,
                        color=self.color, alpha=self.alpha, name=self.name)
        

class HollowCylinder(Cylinder):
    def __init__(self, position, axis, inner_radius, outer_radius, length,
                 color=None, alpha=1, name=''):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.position = position
        axis.Normalize()
        self.axis = axis
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.length = length
        
        # Revolved Profile
        p1 = volmdlr.Point2D((-0.5*self.length, self.inner_radius))
        p2 = volmdlr.Point2D((0.5*self.length, self.inner_radius))
        p3 = volmdlr.Point2D((0.5*self.length, self.outer_radius))
        p4 = volmdlr.Point2D((-0.5*self.length, self.outer_radius))
        l1 = volmdlr.LineSegment2D(p1, p2)
        l2 = volmdlr.LineSegment2D(p2, p3)
        l3 = volmdlr.LineSegment2D(p3, p4)
        l4 = volmdlr.LineSegment2D(p4, p1)
        contour = volmdlr.Contour2D([l1, l2, l3, l4])
        y = axis.RandomUnitNormalVector()
        # contour.MPLPlot()
        RevolvedProfile.__init__(self, position, axis, y, contour, position, axis,
                                 color=color, alpha=alpha, name=name)

        

    def Volume(self):
        return self.length * math.pi* (self.outer_radius**2 - self.inner_radius**2)


    def FreeCADExport(self, ip):
        if self.outer_radius > 0.:
            name = 'primitive'+str(ip)
            re = round(1000*self.outer_radius, 6)
            ri = round(1000*self.inner_radius, 6)
            x, y, z = round((1000*(self.position - self.axis*self.length/2)), 6)
            ax, ay, az = npy.round(self.axis.vector, 6)

            s='C2 = Part.makeCircle({}, fc.Vector({}, {}, {}),fc.Vector({}, {}, {}))\n'.format(re, x, y, z, ax, ay, az)
            s+='W2 = Part.Wire(C2.Edges)\n'
            s+='F2 = Part.Face(W2)\n'

            if self.inner_radius!=0.:
                s+='C1 = Part.makeCircle({}, fc.Vector({}, {}, {}),fc.Vector({}, {}, {}))\n'.format(ri, x, y, z, ax, ay, az)
                s+='W1 = Part.Wire(C1.Edges)\n'
                s+='F1 = Part.Face(W1)\n'
                s+='F2 = F2.cut(F1)\n'

            vx, vy, vz = round(self.axis*self.length*1000, 6)

            s += '{} = F2.extrude(fc.Vector({}, {}, {}))\n'.format(name, vx, vy, vz)
            return s

        else:
            return ''


    def babylon_script(self, name='primitive_mesh'):
        normal_vector1 = self.axis.RandomUnitNormalVector()
        p1 = volmdlr.Point2D((-0.5*self.length, self.outer_radius))
        p2 = volmdlr.Point2D((0.5*self.length, self.outer_radius))
        p3 = volmdlr.Point2D((0.5*self.length, self.inner_radius))
        p4 = volmdlr.Point2D((-0.5*self.length, self.inner_radius))
        l1 = volmdlr.LineSegment2D(p1, p2)
        l2 = volmdlr.LineSegment2D(p2, p3)
        l3 = volmdlr.LineSegment2D(p3, p4)
        l4 = volmdlr.LineSegment2D(p4, p1)
        extruded_profile = RevolvedProfile(self.position, self.axis, normal_vector1,
                                           volmdlr.Contour2D([l1, l2, l3, l4]),
                                           self.position, self.axis, name=self.name)
        return extruded_profile.babylon_script(name=name)


    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'old':
            axis = basis.OldCoordinates(self.axis)
        elif side == 'new':
            axis = basis.NewCoordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')

        if copy:
            return HollowCylinder(position=self.position.frame_mapping(frame, side, copy),
                                  axis=axis,
                                  inner_radius=self.inner_radius, 
                                  outer_radius=self.outer_radius,
                                  length=self.length)
        else:
            self.position.frame_mapping(frame, side, copy)
            self.axis = axis


class HelicalExtrudedProfile(volmdlr.Primitive3D):
    """

    """
    def __init__(self, plane_origin, x, y, axis_point,axis, pitch,
                 outer_contour2d, inner_contours2d=None, name=''):
        volmdlr.Primitive3D.__init__(self, name)
        if inner_contours2d is not None:
            self.inner_contours2d = inner_contours2d
        else:
            self.inner_contours2d = []
            
        self.outer_contour2d = outer_contour2d
        self.axis_point=axis_point
        self.axis=axis
        self.pitch=pitch

        self.inner_contours3d=[c.To3D(plane_origin,x,y) for c in self.inner_contours2d]
        self.outer_contour3d=outer_contour2d.To3D(plane_origin,x,y)


    def FreeCADExport(self,ip,ndigits=3):
        name='primitive{}'.format(ip)
        s="E = []\n"
        for icontour, contour in enumerate(self.outer_contour3d.edges):
            s += contour.FreeCADExport('L_{}'.format(icontour))
            s += 'E.append(Part.Edge(L_{}))\n'.format(icontour)
        s += 'W = Part.Wire(E[:])\n'

#        a1,a2,a3=self.axis
        ap1, ap2, ap3 = self.axis_point
        ap1 = round(ap1*1000, ndigits)
        ap2 = round(ap2*1000, ndigits)
        ap3 = round(ap3*1000, ndigits)

        width = self.axis.Norm()*1000
        direction = bool(self.pitch < 0)
        pitch = round(abs(self.pitch)*1000, ndigits)
        s += "helix = Part.makeHelix({}, {}, 50., 0, {})\n".format(pitch, width, direction)
        s += "helix.translate(fc.Vector({},{},{}))\n".format(ap1, ap2, ap3)

        s += '{} = helix.makePipeShell([W],True,True)\n'.format(name)
        for ic, contour in enumerate(self.inner_contours3d):
            s += "Ei=[]\n"
            s += "helix2 = Part.makeHelix({}, {}, 50., 0, {})\n".format(pitch,1.01*width+pitch,direction)
            s += "helix2.translate(fc.Vector({},{},{}))\n".format(ap1,ap2,ap3-pitch)
            for ip, primitive in enumerate(contour.basis_primitives):
                s += primitive.FreeCADExport('L_{}_{}'.format(ic,ip))
                s += 'Ei.append(Part.Edge(L_{}_{}))\n'.format(ic,ip)
            s+= 'Wi = Part.Wire(Ei[:])\n'
            s+= "{} = {}.cut(helix2.makePipeShell([Wi],True,True))\n".format(name,name)

        return s


class Sweep(volmdlr.Shell3D):
    """
    Sweep a 2D contour along a Wire3D
    """

    def __init__(self, contour2d, wire3d, *, color=None, alpha=1, name=''):
        self.contour2d = contour2d
        self.wire3d = wire3d
        self.frames = []
        
        faces = self.shell_faces()
        volmdlr.Shell3D.__init__(self, faces, color=color, alpha=alpha, name=name)


    def shell_faces(self):
        """
        For now it does not take into account rotation of sections
        """
        faces = []
        last_n = None
        for wire_primitive in self.wire3d.primitives :
            tangent, normal = wire_primitive.frenet(0.)
            if normal is None:
                normal = tangent.deterministic_unit_normal_vector()
            n2 = tangent.Cross(normal)
            contour3d = self.contour2d.To3D(wire_primitive.start,
                                            normal,
                                            n2)

            if wire_primitive.__class__ is volmdlr.LineSegment3D:
                for contour_primitive in contour3d.primitives:
                    faces.append(contour_primitive.extrusion(
                        wire_primitive.direction_vector()))
            elif wire_primitive.__class__ is volmdlr.Arc3D:
                for contour_primitive in contour3d.primitives:
                    faces.append(contour_primitive.revolution(
                        wire_primitive.center,
                        wire_primitive.normal,
                        wire_primitive.angle))
            elif wire_primitive.__class__ is volmdlr.Arc3D:
                for contour_primitive in contour3d.primitives:
                    faces.append(contour_primitive.revolution(
                        wire_primitive.center,
                        wire_primitive.normal,
                        volmdlr.two_pi))
            else:
                raise NotImplementedError(
                    'Unimplemented primitive for sweep: {}'\
                        .format(wire_primitive.__class__.__name__)
                )

        return faces

    # def FreeCADExport(self, ip, ndigits=3):
    #     name = 'primitive{}'.format(ip)
    #     s = "E = []\n"
    #     # for icontour, contour in enumerate(self.contour3d.edges):
    #     #     s += contour.FreeCADExport('L_{}'.format(icontour))
    #     #     s += 'E.append(Part.Edge(L_{}))\n'.format(icontour)
    #     # s += 'contour = Part.Wire(E[:])\n'

    #     # s += "E=[]\n"
    #     # for iwire, wire in enumerate(self.wire3d.edges):
    #     #     s += wire.FreeCADExport('L_{}'.format(iwire))
    #     #     s += 'E.append(Part.Edge(L_{}))\n'.format(iwire)
    #     # s += 'wire = Part.Wire(E[:])\n'

    #     # s += '{} = wire.makePipeShell([contour],True, True)\n'.format(name)
    #     return s

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_wire = self.wire3d.frame_mapping(frame, side, copy)
            return Sweep(self.contour2d, new_wire, color=self.color, alpha=self.alpha, name=self.name)
        else:
            self.wire3d.frame_mapping(frame, side, copy=False)
            for face in self.faces :
                face.frame_mapping(frame, side, copy=False)
            

    def copy(self):
        new_contour2d = self.contour2d.copy()
        new_wire3d = self.wire3d.copy()
        return Sweep(new_contour2d, new_wire3d, color=self.color, alpha=self.alpha, name=self.name)

class Cut(volmdlr.Primitive3D):
    """
    Cut primitive 1 by primitive 2
    """
    def __init__(self,primitive,cut_primitives,name=''):
        volmdlr.Primitive3D.__init__(self,name)
        self.primitive=primitive
        self.cut_primitives = cut_primitives


    def FreeCADExport(self,ip):
        name = 'primitive{}'.format(ip)

        s = self.primitive.FreeCADExport('{}'.format(ip))
        for icp, cut_primitive in enumerate(self.cut_primitives):
            s += cut_primitive.FreeCADExport('{}_{}'.format(ip, icp))
            s += "{} = {}.cut({}_{})\n".format(name, name, name, icp)

        return s

class Fuse(volmdlr.Primitive3D):
    """
    Fuse primitives
    """
    def __init__(self, primitives, name=''):
        volmdlr.Primitive3D.__init__(self, name)
        self.primitives = primitives


    def FreeCADExport(self,ip):
        name = 'primitive{}'.format(ip)


        s = self.primitives[0].FreeCADExport(ip)
        for primitive in self.primitives[1:]:
            s += primitive.FreeCADExport('{}_0'.format(ip))
            s += "{} = {}.fuse({}_0)\n".format(name,name,name)

        return s
    
    
# class Sphere(volmdlr.Primitive3D):
class Sphere(RevolvedProfile):
    def __init__(self, center, radius, color=None, alpha=1., name=''):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.center = center
        self.radius = radius
        self.position = center
        
        # Revolved Profile for complete sphere
        s = volmdlr.Point2D((-self.radius, 0.01*self.radius))
        i = volmdlr.Point2D((0, 1.01*self.radius))
        e = volmdlr.Point2D((self.radius, 0.01*self.radius)) #Not coherent but it works at first, to change !!
        
        # s = volmdlr.Point2D((-self.radius, 0))
        # i = volmdlr.Point2D(((math.sqrt(2)/2)*self.radius,(math.sqrt(2)/2)*self.radius))
        # e = volmdlr.Point2D(((-math.sqrt(2)/2)*self.radius,(-math.sqrt(2)/2)*self.radius)) 
        
        contour = volmdlr.Contour2D([volmdlr.Arc2D(s, i , e), volmdlr.LineSegment2D(s, e)])
        # fig, ax = plt.subplots()
        # c.MPLPlot(ax=ax)
        
        # contour = volmdlr.Contour2D([c])
        axis = volmdlr.X3D
        y = axis.RandomUnitNormalVector()
        RevolvedProfile.__init__(self, center, axis, y, contour, center, axis,
                                 color=color, alpha=alpha, name=name)

    def Volume(self):
        return 4/3*math.pi*self.radius**3

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive'+str(ip)
        r = 1000*self.radius
        x, y, z = round(1000*self.center, ndigits)
        return '{} = Part.makeSphere({}, fc.Vector({}, {}, {}))\n'.format(name,r,x,y,z)

    def babylon_script(self, name='primitive_mesh'):
        p1 = volmdlr.Point2D((-self.radius, 0))
        p2 = volmdlr.Point2D((0, self.radius))
        p3 = volmdlr.Point2D((self.radius, 0))
        line = volmdlr.LineSegment2D(p1, p3)
        arc = volmdlr.Arc2D(p1, p2, p3)
        extruded_profile = RevolvedProfile(self.position, volmdlr.X3D, volmdlr.Y3D,
                                           volmdlr.Contour2D([line, arc]), self.position, volmdlr.X3D, name=self.name)
        return extruded_profile.babylon_script(name=name)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            return Sphere(self.center.frame_mapping(frame, side, copy),
                          self.radius)
        else:
            self.center.frame_mapping(frame, side, copy)


class Circle3D(volmdlr.Contour3D):
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
        self.angle = two_pi
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
                                  for teta in npy.linspace(0, two_pi,
                                                           resolution + 1)][
                                 :-1]
        return tessellation_points_3D

    def Length(self):
        return two_pi * self.radius

    def FreeCADExport(self, name, ndigits=3):
        xc, yc, zc = round(1000 * self.center, ndigits)
        xn, yn, zn = round(self.normal, ndigits)
        return '{} = Part.Circle(fc.Vector({},{},{}),fc.Vector({},{},{}),{})\n'.format(
            name, xc, yc, zc, xn, yn, zn, 1000 * self.radius)

    def Rotation(self, rot_center, axis, angle, copy=True):
        new_center = self.center.Rotation(rot_center, axis, angle, True)
        new_normal = self.normal.Rotation(rot_center, axis, angle, True)
        if copy:
            return Circle3D(new_center, self.radius, new_normal, self.name)
        else:
            self.center = new_center
            self.normal = new_normal

    def Translation(self, offset, copy=True):
        new_frame = self.center.Translation(offset, True)
        if copy:
            return Circle3D(new_frame, self.radius, self.frame,
                            self.name)
        else:
            self.frame = new_frame

    def MPLPlot(self, ax=None, color='k'):
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

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        """
        start point is at intersection of frame.u axis
        """
        start = self.frame.origin + self.radius * self.frame.u
        return start.Rotation(self.frame.origin, self.frame.w,
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
                other_vec.Normalize()
        else:
            normal = object_dict[arguments[1]].v  ### ou w
            other_vec = None
        normal.Normalize()

        return cls(center, radius, normal, arguments[0][1:-1], other_vec)

    def _bounding_box(self):
        """
        """
        u = self.normal.deterministic_unit_normal_vector()
        v = self.normal.Cross(u)
        points = [self.frame.origin + self.radius * v \
                  for v in [self.frame.u,
                            -self.frame.u,
                            self.frame.v,
                            -self.frame.v]]
        return BoundingBox.from_points(points)

    def to_2d(self, plane3d, name=None):
        return Circle2D(plane3d.point3d_to_2d(self.center), self.radius)

    @classmethod
    def from_3_points(cls, point1, point2, point3):
        u1 = (point2 - point1)
        u2 = (point2 - point3)
        try:
            u1.Normalize()
            u2.Normalize()
        except ZeroDivisionError:
            raise ValueError(
                'the 3 points must be distincts')

        normal = u2.Cross(u1)
        normal.Normalize()

        if u1 == u2:
            u2 = normal.Cross(u1)
            u2.Normalize()

        v1 = normal.Cross(u1)  # v1 is normal, equal u2
        v2 = normal.Cross(u2)  # equal -u1

        p11 = 0.5 * (point1 + point2)  # Mid point of segment s,m
        p21 = 0.5 * (point2 + point3)  # Mid point of segment s,m

        l1 = Line3D(p11, p11 + v1)
        l2 = Line3D(p21, p21 + v2)

        try:
            center, _ = l1.MinimumDistancePoints(l2)
        except ZeroDivisionError:
            raise ValueError(
                'Start, end and interior points  of an arc must be distincts')

        radius = (center - point1).Norm()
        return cls(frame=Frame3D(center, u1, normal.Cross(u1), normal),
                   radius=radius)

    def extrusion(self, extrusion_vector):
        if self.normal.is_colinear_to(extrusion_vector):
            u = self.normal.deterministic_unit_normal_vector()
            v = self.normal.Cross(u)
            cylinder = volmdlr.surfaces.CylindricalSurface3D(Frame3D(self.center,
                                                    u,
                                                    v,
                                                    self.normal),
                                            self.radius
                                            )
            return cylinder.rectangular_cut(0, two_pi,
                                            0, extrusion_vector.Norm())
        else:
            raise NotImplementedError('Elliptic faces not handled: dot={}'.format(
                self.normal.Dot(extrusion_vector)
            ))

    def revolution(self, axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                   angle: float):
        line3d = Line3D(axis_point, axis_point + axis)
        tore_center, _ = line3d.point_projection(self.center)
        u =  self.center - tore_center
        u.Normalize()
        v = axis.Cross(u)
        if not math.isclose(self.normal.Dot(u), 0., abs_tol=1e-9):
            raise NotImplementedError(
                'Outside of plane revolution not supported')

        R = tore_center.point_distance(self.center)
        surface = volmdlr.surfaces.ToroidalSurface3D(volmdlr.Frame3D(tore_center, u, v, axis),
                                    R, self.radius)
        return surface.rectangular_cut(0, angle, 0, volmdlr.two_pi)


class Ellipse3D(volmdlr.Contour3D):
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
        normal.Normalize()
        self.normal = normal
        major_dir.Normalize()
        self.major_dir = major_dir
        Contour3D.__init__(self, [self], name=name)

    def tessellation_points(self, resolution=20):
        # plane = Plane3D.from_normal(self.center, self.normal)
        tessellation_points_3D = [self.center + self.major_axis * math.cos(
            teta) * self.major_dir + self.minor_axis * math.sin(
            teta) * self.major_dir.Cross(self.normal) \
                                  for teta in npy.linspace(0, two_pi,
                                                           resolution + 1)][
                                 :-1]
        return tessellation_points_3D

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        xc, yc, zc = npy.round(1000 * self.center.vector, ndigits)
        major_vector = self.center + self.major_axis / 2 * self.major_dir
        xmaj, ymaj, zmaj = npy.round(1000 * major_vector.vector, ndigits)
        minor_vector = self.center + self.minor_axis / 2 * self.normal.Cross(
            self.major_dir)
        xmin, ymin, zmin = npy.round(1000 * minor_vector.vector, ndigits)
        return '{} = Part.Ellipse(fc.Vector({},{},{}), fc.Vector({},{},{}), fc.Vector({},{},{}))\n'.format(
            name, xmaj, ymaj, zmaj, xmin, ymin, zmin, xc, yc, zc)

    def Rotation(self, rot_center, axis, angle, copy=True):
        new_center = self.center.Rotation(rot_center, axis, angle, True)
        new_normal = self.normal.Rotation(rot_center, axis, angle, True)
        new_major_dir = self.major_dir.Rotation(rot_center, axis, angle, True)
        if copy:
            return Ellipse3D(self.major_axis, self.minor_axis, new_center,
                             new_normal, new_major_dir, self.name)
        else:
            self.center = new_center
            self.normal = new_normal
            self.major_dir = new_major_dir

    def Translation(self, offset, copy=True):
        new_center = self.center.Translation(offset, True)
        new_normal = self.normal.Translation(offset, True)
        new_major_dir = self.major_dir.Translation(offset, True)
        if copy:
            return Ellipse3D(self.major_axis, self.minor_axis, new_center,
                             new_normal, new_major_dir, self.name)
        else:
            self.center = new_center
            self.normal = new_normal
            self.major_dir = new_major_dir

    def MPLPlot(self, ax=None, color='k'):
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



class Measure3D(Line3D):
    def __init__(self, point1, point2, color=(1, 0, 0)):
        self.points = [point1, point2]
        self.color = color
        self.distance = Vector3D(self.points[0] - self.points[1]).Norm()
        self.bounding_box = self._bounding_box()

    # !!! no eq defined!
    def __hash__(self):
        return sum([hash(p) for p in self.points])

    def babylon_script(self):
        s = 'var myPoints = [];\n'
        s += 'var point1 = new BABYLON.Vector3({},{},{});\n'.format(
            self.points[0][0], self.points[0][1], self.points[0][2])
        s += 'myPoints.push(point1);\n'
        s += 'var point2 = new BABYLON.Vector3({},{},{});\n'.format(
            self.points[1][0], self.points[1][1], self.points[1][2])
        s += 'myPoints.push(point2);\n'
        s += 'var line = BABYLON.MeshBuilder.CreateLines("lines", {points: myPoints}, scene);\n'
        s += 'line.color = new BABYLON.Color3({}, {}, {});\n'.format(
            self.color[0], self.color[1], self.color[2])
        return s


class BSplineExtrusion(volmdlr.Primitive3D):

    def __init__(self, obj, vectorextru, name=''):
        self.obj = obj
        vectorextru.Normalize()
        self.vectorextru = vectorextru
        if obj.__class__ is Ellipse3D:
            self.points = obj.tessel_points
        else:
            self.points = obj.points

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]
        if object_dict[arguments[1]].__class__ is Ellipse3D:
            ell = object_dict[arguments[1]]
            vectextru = -object_dict[arguments[2]]
            return cls(ell, vectextru, name)

        elif object_dict[arguments[1]].__class__ is BSplineCurve3D:
            bsplinecurve = object_dict[arguments[1]]
            vectextru = object_dict[arguments[2]]
            return cls(bsplinecurve, vectextru, name)
        else:
            raise NotImplementedError  ## a adapter pour les bpsline
