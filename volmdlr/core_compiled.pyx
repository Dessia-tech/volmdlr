#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#cython: language_level=3
"""

Cython functions

"""
# from __future__ import annotations
from typing import TypeVar, List, Tuple, Text, Any, Dict
import math
import warnings
import random

import numpy as npy
from mpl_toolkits.mplot3d import proj3d
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, FancyArrowPatch

from dessia_common import DessiaObject
import plot_data

# =============================================================================

cdef(double, double) Csub2D(double u1, double u2,
                            double v1, double v2):
    return (u1 - v1, u2 - v2)


# =============================================================================

cdef(double, double) Cadd2D(double u1, double u2,
                            double v1, double v2,):
    return (u1 + v1, u2 + v2)


# =============================================================================

cdef(double, double) Cmul2D(double u1, double u2, double value):
    return (u1 * value, u2 * value)

# def mul2D(vector, value):
#    return Cmul2D(vector.x, vector.y, value)

# =============================================================================

cdef double CVector2DDot(double u1, double u2,
                         double v1, double v2):
    return u1 * v1 + u2 * v2

# def Vector2DDot(vector1, vector2):
#    return CVector2DDot(vector1.x, vector1.y,
#                        vector2.x, vector2.y)

# =============================================================================

cdef double CVector2Dnorm(double u1, double u2):
    return (u1 * u1 + u2 * u2)**0.5

# def Vector2Dnorm(vector):
#    return CVector2Dnorm(vector.x, vector.y)

# =============================================================================

cdef(double, double, double) Csub3D(double u1, double u2, double u3,
                                    double v1, double v2, double v3):
    return (u1 - v1, u2 - v2, u3 - v3)

# def sub3D(vector1, vector2):
#    return Csub3D(vector1.x, vector1.y, vector1.z,
#                  vector2.x, vector2.y, vector2.z)

# =============================================================================

cdef(double, double, double) Cadd3D(double u1, double u2, double u3,
                                    double v1, double v2, double v3):
    return (u1 + v1, u2 + v2, u3 + v3)

# def add3D(vector1, vector2):
#    return Cadd3D(vector1.x, vector1.y, vector1.z,
#                  vector2.x, vector2.y, vector2.z)

# =============================================================================

cdef(double, double, double) Cmul3D(double u1, double u2, double u3,
                                    double value):
    return (u1 * value, u2 * value, u3 * value)

# def mul3D(vector, value):
#    return Cmul3D(vector.x, vector.y, vector.z, value)

# =============================================================================

cdef double CVector3DDot(double u1, double u2, double u3,
                         double v1, double v2, double v3):
    return u1 * v1 + u2 * v2 + u3 * v3

# def Vector3DDot(vector1, vector2):
#    return CVector3DDot(vector1.x, vector1.y, vector1.z,
#                        vector2.x, vector2.y, vector2.z)

# =============================================================================

cdef double CVector3Dnorm(double u1, double u2, double u3):
    return (u1 * u1 + u2 * u2 + u3 * u3)**0.5

# def Vector3Dnorm(vector):
#    return CVector3Dnorm(vector.x, vector.y, vector.z)


# =============================================================================

cdef(double, double, double) CVector3D_cross(double u1, double u2, double u3,
                                             double v1, double v2, double v3):
    return (u2 * v3 - u3 * v2, u3 * v1 - u1 * v3, u1 * v2 - u2 * v1)

# def vector3D_cross(vector1, vector2):
#    return C_vector3D_cross(vector1.x, vector1.y, vector1.z,
#                            vector2.x, vector2.y, vector2.z)


# =============================================================================

cdef(double, double, double) C_vector3D_rotation(double vx, double vy, double vz,
                                                 double center_x, double center_y, double center_z,
                                                 double axis_x, double axis_y, double axis_z,
                                                 double angle):

    cdef double ux = vx - center_x
    cdef double uy = vy - center_y
    cdef double uz = vz - center_z

    cdef double cos_angle = math.cos(angle)
    cdef double sin_angle = math.sin(angle)

    cdef double rv1_x = cos_angle * ux
    cdef double rv1_y = cos_angle * uy
    cdef double rv1_z = cos_angle * uz

    rv2_x, rv2_y, rv2_z = Cmul3D(axis_x, axis_y, axis_z,
                                 (1 - cos_angle) * CVector3DDot(
                                         ux, uy, uz,
                                         axis_x, axis_y, axis_z)
                                 )

    rv3_x, rv3_y, rv3_z = CVector3D_cross(axis_x, axis_y, axis_z,
                                          ux, uy, uz)

    return (rv1_x + rv2_x + rv3_x * sin_angle + center_x,
            rv1_y + rv2_y + rv3_y * sin_angle + center_y,
            rv1_z + rv2_z + rv3_z * sin_angle + center_z)


def vector3D_rotation(vector, center, axis, angle):
    return C_vector3D_rotation(vector.x, vector.y, vector.z,
                               center.x, center.y, center.z,
                               axis.x, axis.y, axis.z,
                               angle)


cdef(double, double, double) C_matrix_vector_multiplication3(double M11, double M12, double M13,
                                                             double M21, double M22, double M23,
                                                             double M31, double M32, double M33,
                                                             double v1, double v2, double v3):

    return (M11 * v1 + M12 * v2 + M13 * v3,
            M21 * v1 + M22 * v2 + M23 * v3,
            M31 * v1 + M32 * v2 + M33 * v3)

cdef(double, double) C_matrix_vector_multiplication2(double M11, double M12,
                                                     double M21, double M22,
                                                     double v1, double v2):

    return (M11 * v1 + M12 * v2,
            M21 * v1 + M22 * v2)


cdef(double, double, double,
     double, double, double,
     double, double, double) Cmatrix_multiplication3(double A11, double A12, double A13,
                                                     double A21, double A22, double A23,
                                                     double A31, double A32, double A33,
                                                     double B11, double B12, double B13,
                                                     double B21, double B22, double B23,
                                                     double B31, double B32, double B33):

    return (A11 * B11 + A12 * B21 + A13 * B31,
            A11 * B12 + A12 * B22 + A13 * B32,
            A11 * B13 + A12 * B23 + A13 * B33,
            A21 * B11 + A22 * B21 + A23 * B31,
            A21 * B12 + A22 * B22 + A23 * B32,
            A21 * B13 + A22 * B23 + A23 * B33,
            A31 * B11 + A32 * B21 + A33 * B31,
            A31 * B12 + A32 * B22 + A33 * B32,
            A31 * B13 + A32 * B23 + A33 * B33)


# =============================================================================

def polygon_point_belongs(point, points):

    cdef int i
    cdef int n = len(points)
    cdef bint inside = False
    cdef float x, y, p1x, p1y, p2x, p2y, xints
    x, y = point
    p1x, p1y = points[0]

    for i in range(n + 1):
        p2x, p2y = points[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside

# =============================================================================


cdef(double, (double, double)) CLineSegment2DPointDistance((double, double) p1, (double, double) p2, (double, double) point):
    cdef double t

    ux, uy = Csub2D(p2[0], p2[1], p1[0], p1[1])
    ppx, ppy = Csub2D(point[0], point[1], p1[0], p1[1])

    t = max(0, min(1, CVector2DDot(ppx, ppy, ux, uy) / CVector2Dnorm(ux, uy)**2))
    vx, vy = Cmul2D(ux, uy, t)

    projection = Cadd2D(p1[0], p1[1], vx, vy)
    ppx, ppy = projection[0] - point[0], projection[1] - point[1]
    return CVector2Dnorm(ppx, ppy), projection


def LineSegment2DPointDistance(points, point):
    return CLineSegment2DPointDistance(tuple(points[0]), tuple(points[1]), tuple(point))

# =============================================================================

cdef (double, (double, double, double)) CLineSegment3DPointDistance((double, double, double) p1, (double, double, double) p2, (double, double, double) point):
    cdef double t

    ux, uy, uz = Csub3D(p2[0], p2[1], p2[2], p1[0], p1[1], p1[2])
    ppx, ppy, ppz = Csub3D(point[0], point[1], point[2], p1[0], p1[1], p1[2])
    t = max(0, min(1, CVector3DDot(ppx, ppy, ppz, ux, uy, uz) / CVector3Dnorm(ux, uy, uz)**2))
    vx, vy, vz = Cmul3D(ux, uy, uz, t)
    projection = Cadd3D(p1[0], p1[1], p1[2], vx, vy, vz)
    ppx, ppy, ppz = projection[0]-point[0], projection[1]-point[1], projection[2]-point[2]
    return CVector3Dnorm(ppx, ppy, ppz), projection

def LineSegment3DPointDistance(points, point):
    return CLineSegment3DPointDistance(tuple(points[0]), tuple(points[1]), tuple(point))

# =============================================================================
#  Points, Vectors
# =============================================================================

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def plot2d(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
        FancyArrowPatch.draw(self, renderer)

    def plot(self, ax=None, color='b'):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        points = [self.start, self.end]
        x = [p.x for p in points]
        y = [p.y for p in points]
        z = [p.z for p in points]
        ax.plot(x, y, z, 'o-k')
        return ax


class Vector(DessiaObject):
    """
    Abstract class of vector
    """

    def __radd__(self, other_vector):
        return self + other_vector

    def __rsub__(self, other_vector):
        return self - other_vector

    def __rmul__(self, value):
        return self * value

    def __rtruediv__(self, value):
        return self / value

    def __lt__(self, other_vector):
        return self.norm() < other_vector.norm()

    def __le__(self, other_vector):
        return self.norm() <= other_vector.norm()

    def is_colinear_to(self, other_vector):
        try:
            return math.isclose(abs(self.dot(other_vector)) / self.norm() / other_vector.norm(),
                                1,
                                abs_tol=1e-5)

        except ZeroDivisionError:
            return False

    @classmethod
    def mean_point(cls, points):
        n = 1
        point = points[0].copy()
        for point2 in points[1:]:
            point += point2
            n += 1
        point /= n
        return point

    @classmethod
    def remove_duplicate(cls, points):
        dict_ = {p.approx_hash(): p for p in points}
        return list(dict_.values())


class Vector2D(Vector):
    def __init__(self, x: float, y: float, name=''):
        self.x = x
        self.y = y
        self.name = name

    def __repr__(self):
        return '{}: [{}, {}]'.format(self.__class__.__name__, self.x, self.y)

    def __setitem__(self, key, item):
        if key == 0:
            self.x = item
        elif key == 1:
            self.y = item
        else:
            raise IndexError

    def __getitem__(self, key):
        if key == 0:
            return self.x
        elif key == 1:
            return self.y
        else:
            raise IndexError

    def __add__(self, other_vector):
        return Vector2D(*Cadd2D(self.x, self.y,
                                other_vector.x, other_vector.y))

    def __neg__(self):
        return Vector2D(-self.x, -self.y)

    def __sub__(self, other_vector):
        return Vector2D(*Csub2D(self.x, self.y,
                                other_vector.x, other_vector.y))

    def __mul__(self, value: float):
        return Vector2D(*Cmul2D(self.x, self.y, value))

    def __truediv__(self, value: float):
        if value == 0:
            raise ZeroDivisionError
        return Vector2D(self.x / value,
                        self.y / value)

    def __round__(self, ndigits: int = 6):
        return self.__class__(round(self.x, ndigits),
                              round(self.y, ndigits))

    def __hash__(self):
        """
        hash returns 0 because points are difficult to hash if they are meant
        to be equalized at a given tolerance
        """
        return 0

    def __eq__(self, other_vector):
        return self.is_close(other_vector)

    def is_close(self, other_vector, tol=1e-6):
        if other_vector.__class__.__name__ not in ['Vector2D', 'Point2D']:
            return False
        # return math.isclose(self.x, other_vector.x, abs_tol=tol) \
        # and math.isclose(self.y, other_vector.y, abs_tol=tol)
        return math.isclose(self.point_distance(other_vector), 0, abs_tol=tol)

    def approx_hash(self):
        return round(1e6 * (self.x + self.y))

    def to_dict(self, *args, **kwargs):
        return {'object_class': 'volmdlr.Vector2D',
                'x': self.x, 'y': self.y,
                'name': self.name}

    def copy(self, deep=True, memo=None):
        return self.__class__(self.x, self.y)

    def norm(self):
        """
        :returns: norm of vector
        """
        return CVector2Dnorm(self.x, self.y)

    def normalize(self):
        """
        normalize the vector modifying its coordinates in place
        """
        n = self.norm()
        if math.isclose(n, 0, abs_tol=1e-9):
            raise ZeroDivisionError

        self.x /= n
        self.y /= n

    def dot(self, other_vector):
        return CVector2DDot(self.x,
                            self.y,
                            other_vector.x,
                            other_vector.y)

    def cross(self, other_vector):
        return self.x * other_vector.y - self.y * other_vector.x

    def point_distance(self, other_vector):
        return (self - other_vector).norm()

    def rotation_parameters(self, center: 'Point2D', angle: float):
        """Calculates the parameters to be used in rotation methods"""
        u = self - center
        v2x = math.cos(angle) * u[0] - math.sin(angle) * u[1] + center[0]
        v2y = math.sin(angle) * u[0] + math.cos(angle) * u[1] + center[1]
        return v2x, v2y

    def rotation(self, center: 'Point2D', angle: float):
        """Rotates the vector and returns a new rotated vector"""
        v2x, v2y = self.rotation_parameters(center, angle)
        return self.__class__(v2x, v2y)

    def rotation_inplace(self, center: 'Point2D', angle: float):
        """Rotates the vector and changes its values inplace"""
        v2x, v2y = self.rotation_parameters(center, angle)
        self.x = v2x
        self.y = v2y

    def translation(self, offset: 'Vector2D'):
        """
        Translates the vector and returns a new translated vector
        :param offset: an other Vector2D
        """
        v2x = self.x + offset[0]
        v2y = self.y + offset[1]
        return self.__class__(v2x, v2y)

    def translation_inplace(self, offset: 'Vector2D'):
        """Translates the vector and changes its values inplace"""
        v2x = self.x + offset[0]
        v2y = self.y + offset[1]
        self.x = v2x
        self.y = v2y

    def frame_mapping(self, frame: 'Frame2D', side: str):
        """
        Changes vector frame_mapping and return a new vector
        side = 'old' or 'new'
        """
        if side == 'old':
            new_vector = frame.old_coordinates(self)
        if side == 'new':
            new_vector = frame.new_coordinates(self)
        return new_vector

    def frame_mapping_inplace(self, frame: 'Frame2D', side: str):
        """
        Changes vector frame_mapping in place
        side = 'old' or 'new'
        """
        if side == 'old':
            new_vector = frame.old_coordinates(self)
        if side == 'new':
            new_vector = frame.new_coordinates(self)
        self.x = new_vector.x
        self.y = new_vector.y

    def to_3d(self, plane_origin, vx, vy):
        return Vector3D(plane_origin.x + vx.x * self.x + vy.x * self.y,
                        plane_origin.y + vx.y * self.x + vy.y * self.y,
                        plane_origin.z + vx.z * self.x + vy.z * self.y,
                        )

    def to_point(self):
        return Point2D(self.x, self.y)

    def normal_vector(self):
        n = Vector2D(-self.y, self.x)
        return n

    def unit_normal_vector(self):
        n = self.normal_vector()
        n.normalize()
        return n

    def deterministic_unit_normal_vector(self):
        return self.unit_normal_vector()

    @classmethod
    def random(cls, xmin, xmax, ymin, ymax):
        return cls(random.uniform(xmin, xmax),
                   random.uniform(ymin, ymax))

    def plot(self, amplitude=0.5, width=None, head_width=None, origin=None,
             ax=None, color='k', line=False, label=None, normalize=False):
        if origin is None:
            origin = Vector2D(0., 0.)

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        if math.isclose(self.norm(), 0, abs_tol=1e-9):
            point = origin.copy()
            point.plot(ax=ax, color=color)
            return ax

        if width is None:
            width = 0.001 * 5 * amplitude
        if head_width is None:
            head_width = 0.3 * amplitude

        if not normalize:
            ax.add_patch(FancyArrow(origin[0], origin[1],
                                    self.x * amplitude, self.y * amplitude,
                                    width=width,
                                    head_width=head_width,
                                    length_includes_head=True,
                                    color=color))
        else:
            normalized_vector = self.copy()
            normalized_vector.normalize()
            ax.add_patch(FancyArrow(origin[0], origin[1],
                                    normalized_vector.x * amplitude,
                                    normalized_vector.y * amplitude,
                                    width=width,
                                    head_width=head_width,
                                    length_includes_head=True,
                                    color=color))

        if line:
            style = '-' + color
            linestyle = '-.'
            origin = Point2D(origin)
            p1, p2 = origin, origin + self
            u = p2 - p1
            p3 = p1 - 3 * u
            p4 = p2 + 4 * u
            ax.plot([p3[0], p4[0]], [p3[1], p4[1]], style, linestyle=linestyle)

        if label is not None:
            ax.text(*(origin + self * amplitude), label)

        return ax


X2D = Vector2D(1, 0)
Y2D = Vector2D(0, 1)


class Point2D(Vector2D):
    def __init__(self, x: float, y: float, name: Text = ''):
        Vector2D.__init__(self, x=x, y=y, name=name)

    def __add__(self, other_vector):
        return Point2D(*Cadd2D(self.x, self.y, other_vector.x, other_vector.y))

    def __neg__(self):
        return Point2D(-self.x, -self.y)

    def __sub__(self, other_vector):
        return Point2D(*Csub2D(self.x, self.y,
                               other_vector.x, other_vector.y))

    def __mul__(self, value: float):
        return Point2D(*Cmul2D(self.x, self.y, value))

    def __truediv__(self, value: float):
        if value == 0:
            raise ZeroDivisionError
        return Point2D(self.x / value,
                       self.y / value)

    def to_dict(self, *args, **kwargs):
        return {'object_class': 'volmdlr.Point2D',
                'x': self.x, 'y': self.y,
                'name': self.name}

    def to_3d(self, plane_origin, vx, vy):
        return Point3D(plane_origin.x + vx.x * self.x + vy.x * self.y,
                       plane_origin.y + vx.y * self.x + vy.y * self.y,
                       plane_origin.z + vx.z * self.x + vy.z * self.y)

    def to_vector(self):
        return Vector2D(self.x, self.y)

    def plot(self, ax=None, color='k', alpha=1, plot_points=True):
        if ax is None:
            fig, ax = plt.subplots()

        ax.plot([self.x], [self.y], color=color, alpha=alpha, marker='o')
        return ax

    def point_distance(self, other_point: 'Point2D'):
        return (self - other_point).norm()

    @classmethod
    def line_intersection(cls, line1, line2, curvilinear_abscissa=False):
        #        point11, point12 = line1
        #        point21, point22 = line2
        (x1, y1), (x2, y2) = line1
        (x3, y3), (x4, y4) = line2

#        x1 = line1.points[0][0]
#        y1 = line1.points[0][1]
#        x2 = line1.points[1][0]
#        y2 = line1.points[1][1]
#        x3 = line2.points[0][0]
#        y3 = line2.points[0][1]
#        x4 = line2.points[1][0]
#        y4 = line2.points[1][1]

        denominateur = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if math.isclose(denominateur, 0, abs_tol=1e-6):
            if not curvilinear_abscissa:
                return None
            else:
                return None, None, None
        else:
            x = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)
            x = x / denominateur
            y = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)
            y = y / denominateur
            if not curvilinear_abscissa:
                return cls(x, y)
            else:
                t = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)
                t = t / denominateur
                u = (x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)
                u = -u / denominateur
                return (cls(x, y), t, u)

    @classmethod
    def segment_intersection(cls, segment1, segment2,
                             curvilinear_abscissa=False):
        x1 = segment1.points[0].x
        y1 = segment1.points[0].y
        x2 = segment1.points[1].x
        y2 = segment1.points[1].y
        x3 = segment2.points[0].x
        y3 = segment2.points[0].y
        x4 = segment2.points[1].x
        y4 = segment2.points[1].y

        denominateur = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if math.isclose(denominateur, 0, abs_tol=1e-6):
            if not curvilinear_abscissa:
                return None
            else:
                return None, None, None

        t = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)
        t = t / denominateur
        u = (x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)
        u = -u / denominateur
        if (0 <= t and t <= 1) or (0 <= u and u <= 1):
            x = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)
            x = x / denominateur
            y = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)
            y = y / denominateur
            if not curvilinear_abscissa:
                return cls((x, y))
            else:
                return (cls((x, y)), t, u)
        else:
            if not curvilinear_abscissa:
                return None
            else:
                return None, None, None

    def plot_data(self, marker=None, color='black', size=1,
                  opacity=1, arrow=False, stroke_width=None):
        return plot_data.Point2D(self.x, self.y)

    @classmethod
    def middle_point(cls, point1, point2):
        return (point1 + point2) * 0.5

    @classmethod
    def line_projection(cls, point, line):
        p1, p2 = line[0], line[1]
        n = line.unit_normal_vector()
        pp1 = point - p1
        return pp1 - pp1.dot(n) * n + p1

    def nearest_point(self, points):
        distances = []
        for p in points:
            distances.append(self.point_distance(p))
        return points[distances.index(min(distances))]


    def axial_symmetry(self, line):
        '''
        finds out the symmetric point according to a line
        '''

        point_projection = line.point_projection(self)[0]
        point_symmetry = point_projection + (point_projection - self)

        return point_symmetry

    def coordinates(self):
        '''
        gets x,y coordinates of a point2d
        '''

        return (self.x, self.y)


O2D = Point2D(0, 0)


class Vector3D(Vector):

    def __init__(self, x: float, y: float, z: float, name: Text = ''):
        self.x = x
        self.y = y
        self.z = z
        self.name = name

    def __repr__(self):
        return '{}: [{}, {}, {}]'.format(self.__class__.__name__, self.x, self.y, self.z)

    def __setitem__(self, key, item):
        if key == 0:
            self.x = item
        elif key == 1:
            self.y = item
        elif key == 2:
            self.z = item
        else:
            raise IndexError

    def __getitem__(self, key):
        if key == 0:
            return self.x
        elif key == 1:
            return self.y
        elif key == 2:
            return self.z
        else:
            raise IndexError

    def __add__(self, other_vector):
        return Vector3D(*Cadd3D(self.x, self.y, self.z,
                                other_vector.x,
                                other_vector.y,
                                other_vector.z))

    def __neg__(self):
        return Vector3D(-self.x, -self.y, -self.z)

    def __sub__(self, other_vector):
        return Vector3D(*Csub3D(self.x, self.y, self.z,
                                other_vector.x,
                                other_vector.y,
                                other_vector.z))

    def __mul__(self, value):
        return Vector3D(*Cmul3D(self.x, self.y, self.z, value))

    def __truediv__(self, value):
        if value == 0:
            raise ZeroDivisionError
        return Vector3D(self.x / value,
                        self.y / value,
                        self.z / value)

    def __round__(self, ndigits: int = 6):
        return self.__class__(round(self.x, ndigits),
                              round(self.y, ndigits),
                              round(self.z, ndigits))

    def __hash__(self):
        """
        hash returns 0 because points are difficult to hash if they are meant
        to be equalized at a given tolerance
        """

        return 0

    def __eq__(self, other_vector: 'Vector3D'):
        return self.is_close(other_vector)

    def is_close(self, other_vector, tol=1e-6):
        if other_vector.__class__.__name__ not in ['Vector3D', 'Point3D']:
            return False
        # return math.isclose(self.x, other_vector.x, abs_tol=tol) \
        # and math.isclose(self.y, other_vector.y, abs_tol=tol) \
        # and math.isclose(self.z, other_vector.z, abs_tol=tol)
        return math.isclose(self.point_distance(other_vector), 0, abs_tol=tol)

    def approx_hash(self):
        return round(1e6 * (self.x + self.y + self.z))

    def to_dict(self, *args, **kwargs):
        return {'object_class': 'volmdlr.Vector3D',
                'x': self.x, 'y': self.y, 'z': self.z,
                'name': self.name}

    @classmethod
    def dict_to_object(cls, dict_, global_dict=None, pointers_memo: Dict[str, Any] = None, path: str = '#'):
        return Vector3D(dict_['x'], dict_['y'], dict_['z'], dict_.get('name', ''))

    def dot(self, other_vector):
        return CVector3DDot(self.x, self.y, self.z,
                            other_vector.x, other_vector.y, other_vector.z)

    def cross(self, other_vector: 'Vector3D') -> 'Vector3D':
        return self.__class__(*CVector3D_cross(self.x, self.y, self.z,
                                               other_vector.x,
                                               other_vector.y,
                                               other_vector.z))

    def norm(self) -> float:
        return CVector3Dnorm(self.x, self.y, self.z)

    def normalize(self) -> None:
        """
        normalize the vector modifying its coordinates
        """
        n = self.norm()
        if n == 0:
            raise ZeroDivisionError

        self.x /= n
        self.y /= n
        self.z /= n

    def point_distance(self, point2: 'Vector3D') -> float:
        return (self - point2).norm()

    def rotation(self, center: 'Point3D', axis: 'Vector3D', angle: float):

        """
        Rotation of angle around axis. Returns a new rotated vector
        Used Rodrigues Formula:
            https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        """
        vector2 = vector3D_rotation(self, center, axis, angle)
        return self.__class__(*vector2)

    def rotation_inplace(self, center: 'Point3D', axis: 'Vector3D', angle: float):
        """Rotation of angle around axis.
        The vector parameters are changed inplace
        Used Rodrigues Formula:
            https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula"""
        vector2 = vector3D_rotation(self, center, axis, angle)
        self.x = vector2[0]
        self.y = vector2[1]
        self.z = vector2[2]

    @staticmethod
    def axis_rotation_parameters(axis1_value, axis2_value, angle):

        """
        Calcules new axis1 and axis2 new values
        after vector rotation
        """
        cos_angle = math.cos(angle)
        sin_angle = math.sin(angle)

        axis1 = cos_angle * axis1_value + sin_angle * axis2_value
        axis2 = -sin_angle * axis1_value + cos_angle * axis2_value
        return axis1, axis2

    def x_rotation(self, angle:float):
        """
        rotation of angle around X axis and returns a new vector as a result
        """
        y1, z1 = self.axis_rotation_parameters(self.y, self.z, angle)

        return Point3D(self.x, y1, z1)

    def x_rotation_inplace(self, angle: float):
        """
        Rotation of angle around X axis and
        changes the vector parameters inplace
        """
        y1, z1 = self.axis_rotation_parameters(self.y, self.z, angle)
        self.y = y1
        self.z = z1

    def y_rotation(self, angle:float):
        """
        rotation of angle around Y axis and
        returns a new vector as result.
        """
        z1, x1 = self.axis_rotation_parameters(self.z, self.x, angle)
        return Point3D(x1, self.y, z1)

    def y_rotation_inplace(self, angle):
        """
        Rotation of vector around the Y axis and
        changes its parameters inplace
        """
        z1, x1 = self.axis_rotation_parameters(self.z, self.x, angle)
        self.x = x1
        self.z = z1

    def z_rotation(self, angle:float):
        """
        rrotation of angle around Z axis and
        returns a new vector as result.
        """
        x1, y1 = self.axis_rotation_parameters(self.x, self.y, angle)
        return Point3D(x1, y1, self.z)

    def z_rotation_inplace(self, angle: float):
        """Rotation of vector around the Z axis and
        changes its parameters inplace"""
        x1, y1 = self.axis_rotation_parameters(self.x, self.y, angle)
        self.x = x1
        self.y = y1


    def translation(self, offset):
        """
        Translates the vector and returns a new translated vector
        :param offset: an other Vector3D
        """
        return self + offset

    def translation_inplace(self, offset):
        """
        Translates the vector and changes its values inplace
        :param offset: an other Vector3D
        """
        self.x += offset[0]
        self.y += offset[1]
        self.z += offset[2]

    def frame_mapping(self, frame: 'Frame3D', side: str):
        """
        Changes vector frame_mapping and return a new vector
        side = 'old' or 'new'
        """
        if side == 'old':
            new_vector = frame.old_coordinates(self)

        if side == 'new':
            new_vector = frame.new_coordinates(self)
        return new_vector

    def frame_mapping_inplace(self, frame: 'Frame3D', side: str):
        """
        Changes vector frame_mapping in place
        side = 'old' or 'new'
        """
        if side == 'old':
            new_vector = frame.old_coordinates(self)

        if side == 'new':
            new_vector = frame.new_coordinates(self)
        self.x = new_vector.x
        self.y = new_vector.y
        self.z = new_vector.z

    def plane_projection3d(self, plane_origin, x, y):
        z = x.cross(y)
        z.normalize()
        return self - z.dot(self - plane_origin) * z

    def plane_projection2d(self, plane_origin, x, y):
        z = x.cross(y)
        z.normalize()
        p3d = self - (self - plane_origin).dot(z) * z
        u1 = p3d.dot(x)
        u2 = p3d.dot(y)
        return Point2D(u1, u2)

    def to_2d(self, plane_origin, x, y):
        x2d = self.dot(x) - plane_origin.dot(x)
        y2d = self.dot(y) - plane_origin.dot(y)
        return Point2D(x2d, y2d)

    def random_unit_normal_vector(self):
        """
        Returns a random normal vector
        """
        v = Vector3D.random(0, 1, 0, 1, 0, 1)

        v = v - v.dot(self) * self / (self.norm()**2)
        v.normalize()
        return v

    def deterministic_unit_normal_vector(self):
        """
        Returns a deterministic normal vector
        """
        v = X3D
        if not math.isclose(self.y, 0, abs_tol=1e-7) \
                or not math.isclose(self.z, 0, abs_tol=1e-7):
            v = X3D
        else:
            v = Y3D
        v = v - v.dot(self) * self / (self.norm()**2)
        v.normalize()
        return v

    def copy(self, deep=True, memo=None):
        return Vector3D(self.x, self.y, self.z)

    @classmethod
    def random(cls, xmin, xmax, ymin, ymax, zmin, zmax):
        return cls(random.uniform(xmin, xmax),
                   random.uniform(ymin, ymax),
                   random.uniform(zmin, zmax))

    def to_point(self):
        return Point3D(self.x, self.y, self.z)

    @classmethod
    def from_step(cls, arguments, object_dict):
        if type(arguments[1]) is int:
            # VECTOR
            return cls(*object_dict[arguments[1]], arguments[0][1:-1])
        else:
            # DIRECTION
            # return cls(*[float(i)/1000 for i in arguments[1][1:-1].split(",")],
            #             arguments[0][1:-1])
            return cls(*[float(i) for i in arguments[1][1:-1].split(",")],
                       arguments[0][1:-1])

    def to_step(self, current_id, vector=False, vertex=False):
        if vertex:
            return self.to_point().to_step(current_id=current_id, vertex=True)
        content = "#{} = DIRECTION('{}',({:.6f},{:.6f},{:.6f}));\n"\
            .format(current_id, self.name,
                    self.x, self.y, self.z)
        if vector:
            content += "#{} = VECTOR('{}',#{},1.);\n".format(current_id + 1,
                                                             self.name,
                                                             current_id)
            current_id += 1
        return content, current_id

    def plot(self, ax=None, starting_point=None, color=''):
        if starting_point is None:
            starting_point = Point3D(0, 0, 0)
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        xs = [starting_point[0], self.x + starting_point[0]]
        ys = [starting_point[1], self.y + starting_point[1]]
        zs = [starting_point[2], self.z + starting_point[2]]
        if color:
            a = Arrow3D(xs, ys, zs, mutation_scale=10, lw=3, arrowstyle="-|>", color=color)
        else:
            a = Arrow3D(xs, ys, zs, mutation_scale=10, lw=3, arrowstyle="-|>")
        ax.add_artist(a)
        return ax


X3D = Vector3D(1, 0, 0)
Y3D = Vector3D(0, 1, 0)
Z3D = Vector3D(0, 0, 1)


class Point3D(Vector3D):
    _standalone_in_db = False

    def __init__(self, x: float, y: float, z: float, name: Text = ''):
        Vector3D.__init__(self, x, y, z, name)

    def __add__(self, other_vector):
        return Point3D(*Cadd3D(self.x, self.y, self.z,
                               other_vector.x,
                               other_vector.y,
                               other_vector.z))

    def __neg__(self):
        return Point3D(-self.x, -self.y, -self.z)

    def __sub__(self, other_vector):
        return Point3D(*Csub3D(self.x, self.y, self.z,
                               other_vector.x, other_vector.y, other_vector.z))

    def __mul__(self, value):
        return Point3D(*Cmul3D(self.x, self.y, self.z, value))

    def __truediv__(self, value):
        if value == 0:
            raise ZeroDivisionError
        return Point3D(self.x / value,
                       self.y / value,
                       self.z / value)

    def copy(self, deep=True, memo=None):
        return Point3D(self.x, self.y, self.z)

    def to_dict(self, *args, **kwargs):
        return {'object_class': 'volmdlr.Point3D',
                'x': self.x, 'y': self.y, 'z': self.z,
                'name': self.name}

    @classmethod
    def dict_to_object(cls, dict_, global_dict=None, pointers_memo: Dict[str, Any] = None, path: str = '#'):
        return Point3D(dict_['x'], dict_['y'], dict_['z'], dict_.get('name', ''))

    def plot(self, ax=None, color='k', alpha=1, marker='o'):

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        ax.plot([self.x], [self.y], [self.z], color=color, alpha=alpha,
                marker=marker)
        return ax

    def to_2d(self, plane_origin, x, y):
        x2d = self.dot(x) - plane_origin.dot(x)
        y2d = self.dot(y) - plane_origin.dot(y)
        return Point2D(x2d, y2d)

    @classmethod
    def from_step(cls, arguments, object_dict):
        return cls(*[float(i) / 1000 for i in arguments[1][1:-1].split(",")],
                   arguments[0][1:-1])

    def to_vector(self):
        return Vector3D(self.x, self.y, self.z)

    def point_distance(self, point2: 'Point3D') -> float:
        return (self - point2).norm()

    @classmethod
    def middle_point(cls, point1, point2):
        return (point1 + point2) * 0.5

    def to_step(self, current_id, vertex=False):
        content = "#{} = CARTESIAN_POINT('{}',({:.6f},{:.6f},{:.6f}));\n"\
                        .format(current_id, self.name,
                                1000. * self.x,
                                1000. * self.y,
                                1000. * self.z)
        if vertex:
            content += "#{} = VERTEX_POINT('{}',#{});\n".format(current_id + 1,
                                                                self.name,
                                                                current_id)
            current_id += 1

        return content, current_id

    def babylon_script(self):
        s = 'var sphere = BABYLON.MeshBuilder.CreateSphere("point", {diameter: 0.05}, scene);\n'
        s += "sphere.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n".format(self.x, self.y, self.z)
        s += 'var mat = new BABYLON.StandardMaterial("mat", scene);\n'
        s += 'mat.diffuseColor = new BABYLON.Color3(1, 0, 0);\n'
        s += 'sphere.material = mat;\n'
        return s

    def nearest_point(self, points):
        distances = []
        for p in points:
            distances.append(self.point_distance(p))
        return points[distances.index(min(distances))]

    def coordinates(self):
        '''
        gets x,y,z coordinates of a point3d
        '''

        return (self.x, self.y, self.z)


O3D = Point3D(0, 0, 0)


# =============================================================================
#  Basis, Frames
# =============================================================================

class Matrix22:
    def __init__(self, M11: float, M12: float, M21: float, M22: float):
        self.M11 = M11
        self.M12 = M12
        self.M21 = M21
        self.M22 = M22

    def __add__(self, other_matrix):
        return Matrix22(self.M11 + other_matrix.M11,
                        self.M12 + other_matrix.M12,
                        self.M21 + other_matrix.M21,
                        self.M22 + other_matrix.M22,
                        )

    def __mul__(self, other_matrix):
        return Matrix22(self.M11 * other_matrix.M11 + self.M12 * other_matrix.M21,
                        self.M11 * other_matrix.M12 + self.M12 * other_matrix.M22,
                        self.M21 * other_matrix.M11 + self.M22 * other_matrix.M21,
                        self.M21 * other_matrix.M12 + self.M22 * other_matrix.M22)

    def vector_multiplication(self, vector):
        u1, u2 = C_matrix_vector_multiplication2(self.M11, self.M12,
                                                 self.M21, self.M22,
                                                 vector.x, vector.y)

        return vector.__class__(u1, u2)

    def determinent(self):
        return self.M11 * self.M22 - self.M12 * self.M21

    def inverse(self):
        det = self.determinent()
        if not math.isclose(det, 0, abs_tol=1e-10):
            det_inv = 1 / self.determinent()
            return Matrix22(det_inv * self.M22, -det_inv * self.M12,
                            -det_inv * self.M21, det_inv * self.M11)
        else:
            raise ValueError('The matrix is singular')

    def vector_multiplication(self, vector):
        return vector.__class__(self.M11 * vector.x + self.M12 * vector.y,
                                self.M21 * vector.x + self.M22 * vector.y)


class Matrix33:
    def __init__(self, M11: float, M12: float, M13: float,
                 M21: float, M22: float, M23: float,
                 M31: float, M32: float, M33: float):
        self.M11 = M11
        self.M12 = M12
        self.M13 = M13
        self.M21 = M21
        self.M22 = M22
        self.M23 = M23
        self.M31 = M31
        self.M32 = M32
        self.M33 = M33

    def __add__(self, other_matrix):
        return Matrix33(self.M11 + other_matrix.M11,
                        self.M12 + other_matrix.M12,
                        self.M13 + other_matrix.M13,
                        self.M21 + other_matrix.M21,
                        self.M22 + other_matrix.M22,
                        self.M23 + other_matrix.M23,
                        self.M31 + other_matrix.M31,
                        self.M32 + other_matrix.M32,
                        self.M33 + other_matrix.M33)

    def __mul__(self, other_matrix):
        M11, M12, M13, M21, M22, M23, M31, M32, M33 = Cmatrix_multiplication3(self.M11, self.M12, self.M13,
                                                                              self.M21, self.M22, self.M23,
                                                                              self.M31, self.M32, self.M33,
                                                                              other_matrix.M11, other_matrix.M12, other_matrix.M13,
                                                                              other_matrix.M21, other_matrix.M22, other_matrix.M23,
                                                                              other_matrix.M31, other_matrix.M32, other_matrix.M33)

        return Matrix33(M11, M12, M13, M21, M22, M23, M31, M32, M33)

    def __repr__(self):
        s = '[{} {} {}]\n[{} {} {}]\n[{} {} {}]\n'.format(self.M11, self.M12, self.M13,
                                                          self.M21, self.M22, self.M23,
                                                          self.M31, self.M32, self.M33)
        return s

    def float_multiplication(self, float_value):
        return Matrix33(self.M11 * float_value, self.M12 * float_value, self.M13 * float_value,
                        self.M21 * float_value, self.M22 * float_value, self.M23 * float_value,
                        self.M31 * float_value, self.M32 * float_value, self.M33 * float_value)

    def vector_multiplication(self, vector):
        u1, u2, u3 = C_matrix_vector_multiplication3(self.M11, self.M12, self.M13,
                                                     self.M21, self.M22, self.M23,
                                                     self.M31, self.M32, self.M33,
                                                     vector.x, vector.y, vector.z)

        return vector.__class__(u1, u2, u3)

    def determinent(self):
        det = self.M11 * self.M22 * self.M33 + self.M12 * self.M23 * self.M31 \
            + self.M13 * self.M21 * self.M32 - self.M13 * self.M22 * self.M31 \
            - self.M23 * self.M32 * self.M11 - self.M33 * self.M12 * self.M21
        return det

    def inverse(self):
        det = self.determinent()

        if not math.isclose(det, 0, abs_tol=1e-10):
            det_inv = 1 / det
            return Matrix33(det_inv * (self.M22 * self.M33 - self.M23 * self.M32),  # a22a33−a23a32
                            det_inv * (self.M13 * self.M32 - self.M12 * self.M33),  # a13a32−a12a33
                            det_inv * (self.M12 * self.M23 - self.M13 * self.M22),  # a12a23−a13a22
                            det_inv * (self.M23 * self.M31 - self.M21 * self.M33),  # a23a31−a21a33
                            det_inv * (self.M11 * self.M33 - self.M13 * self.M31),  # a11a33−a31a13
                            det_inv * (self.M21 * self.M13 - self.M23 * self.M11),  # a13a21−a23a11
                            det_inv * (self.M21 * self.M32 - self.M31 * self.M22),  # a21a32−a31a22
                            det_inv * (self.M12 * self.M31 - self.M32 * self.M11),  # a12a31−a32a11
                            det_inv * (self.M11 * self.M22 - self.M21 * self.M12)  # a11a22−a21a12
                            )
        else:
            raise ValueError('The matrix is singular')

    @classmethod
    def random_matrix(cls, minimum=0, maximum=1):
        range_ = maximum - minimum
        return cls(*[minimum + range_ * random.random() for _ in range(9)])

    def to_numpy(self):
        return npy.array([[self.M11, self.M12, self.M13],
                          [self.M21, self.M22, self.M23],
                          [self.M31, self.M32, self.M33]])


class Basis(DessiaObject):
    """
    Abstract class of a basis
    """

    def __contains__(self, vector):
        return vector in self.vectors

    def __hash__(self):
        """
        hash returns 0 because points are difficult to hash if they are meant
        to be equalized at a given tolerance
        """
        0

    def copy(self, deep=True, memo=None):
        return self.__class__(*self.vectors)


class Basis2D(Basis):
    """
    Defines a 2D basis
    :param u:Vector2D: first vector of the basis
    :param v:Vector2D: second vector of the basis
    """

    def __init__(self, u: Vector2D, v: Vector2D, name: Text = ''):
        self.u = u
        self.v = v
        self.name = name

    def __eq__(self, other_basis):
        if other_basis.__class__.__name__ != self.__class__.__name__:
            return False
        all_equal = all([other_vector == vector
                         for other_vector, vector
                         in zip([other_basis.u, other_basis.v], [self.u, self.v])])
        return all_equal

    def __neg__(self):
        Pinv = self.inverse_transfer_matrix()
        return Basis2D(Vector3D(Pinv[:, 0]),
                       Vector3D(Pinv[:, 1]))

    def __repr__(self):
        return '{}: U={}, V={}'.format(self.__class__.__name__, *self.vectors)

    def _get_vectors(self):
        return (self.u, self.v)

    vectors = property(_get_vectors)

    def to_dict(self, *args, **kwargs):
        return {'object_class': 'volmdlr.Basis2D',
                'name': self.name,
                'u': self.u.to_dict(),
                'v': self.v.to_dict()
                }

    def to_frame(self, origin: Point3D) -> 'Frame3D':
        return Frame2D(origin, self.u, self.v)

    def transfer_matrix(self):
        return Matrix22(self.u.x, self.v.x,
                        self.u.y, self.v.y)
        # return npy.array([[self.u[0], self.v[0]],
        #                   [self.u[1], self.v[1]]])

    # def inverse_transfer_matrix(self):
    #     det = self.u[0]*self.v[1] - self.v[0]*self.u[1]
    #     if not math.isclose(det, 0, abs_tol=1e-10):
    #         return 1/det * npy.array([[self.v[1], -self.v[0]],
    #                                  [-self.u[1], self.u[0]]])
    #     else:
    #         raise ZeroDivisionError
    def inverse_transfer_matrix(self):
        return self.transfer_matrix().inverse()

    def new_coordinates(self, vector):
        matrix = self.inverse_transfer_matrix()
        return matrix.vector_multiplication(vector)

    def old_coordinates(self, point):
        matrix = self.transfer_matrix()
        return matrix.vector_multiplication(point)

    # def new_coordinates(self, vector):
    #     matrix = self.inverse_transfer_matrix()
    #     return Point2D((matrix[0][0]*vector.x + matrix[0][1]*vector.y,
    #                      matrix[1][0]*vector.x + matrix[1][1]*vector.y))
    #
    #
    # def old_coordinates(self, vector):
    #     matrix = self.transfer_matrix()
    #     return Point2D(matrix[0][0]*vector.x + matrix[0][1]*vector.y,
    #                    matrix[1][0]*vector.x + matrix[1][1]*vector.y)

    def rotation(self, angle: float):
        """Rotates the base and returns a new rotated one"""

        center = O2D
        new_u = self.u.rotation(center, angle)
        new_v = self.v.rotation(center, angle)
        return Basis2D(new_u, new_v)

    def roation_inplace(self, angle: float):
        """Rotates the base and changes its parameters inplace"""
        center = O2D
        new_u = self.u.rotation(center, angle)
        new_v = self.v.rotation(center, angle)
        self.u = new_u
        self.v = new_v

    def copy(self, deep=True, memo=None):
        return Basis2D(self.u, self.v)

    def normalize(self):
        """
        normalize the basis modifying its coordinates in place
        """
        self.u.normalize()
        self.v.normalize()


XY = Basis2D(X2D, Y2D)


class Basis3D(Basis):
    """
    Defines a 3D basis

    :param u:Vector3D: first vector of the basis
    :param v:Vector3D: second vector of the basis
    :param w:Vector3D: third vector of the basis
    """
    _standalone_in_db = False

    # TODO: create a Basis and Frame class to mutualize between 2D and 2D
    def __init__(self, u: Vector3D, v: Vector3D, w: Vector3D, name: Text = ''):
        self.u = u
        self.v = v
        self.w = w
        self.name = name

    def __eq__(self, other_basis):
        if other_basis.__class__.__name__ != self.__class__.__name__:
            return False

        for other_vector, vector in zip([other_basis.u,
                                         other_basis.v, other_basis.w],
                                        [self.u, self.v, self.w]):
            if other_vector != vector:
                return False
        return True

    def __hash__(self):
        """
        hash returns 0 because points are difficult to hash if they are meant
        to be equalized at a given tolerance
        """
        return 0

    def __add__(self, other_basis):
        M = self.transfer_matrix() * other_basis.transfer_matrix()
        return Basis3D(Vector3D(M.M11, M.M21, M.M31),
                       Vector3D(M.M12, M.M22, M.M32),
                       Vector3D(M.M13, M.M23, M.M33))

    def __neg__(self):
        M = self.inverse_transfer_matrix()
        return Basis3D(Vector3D(M.M11, M.M21, M.M31),
                       Vector3D(M.M12, M.M22, M.M32),
                       Vector3D(M.M13, M.M23, M.M33))

    def __sub__(self, other_frame):
        P1inv = other_frame.inverse_transfer_matrix()
        P2 = self.transfer_matrix()
        M = P1inv * P2
        return Basis3D(Vector3D(M.M11, M.M21, M.M31),
                       Vector3D(M.M12, M.M22, M.M32),
                       Vector3D(M.M13, M.M23, M.M33))

    def __round__(self, ndigits: int = 6):
        return self.__class__((round(self.u, ndigits),
                               round(self.v, ndigits),
                               round(self.w, ndigits)))

    def __repr__(self):
        return '{}: U={}, V={}, W={}'.format(self.__class__.__name__, *self.vectors)

    def _get_vectors(self):
        return (self.u, self.v, self.w)

    def to_dict(self, *args, **kwargs):
        return {'object_class': 'volmdlr.Basis3D',
                'name': self.name,
                'u': self.u.to_dict(),
                'v': self.v.to_dict(),
                'w': self.w.to_dict()
                }

    vectors = property(_get_vectors)

    # TODO: transform to annotation when available
    @classmethod
    def from_two_vectors(cls, vector1: Vector3D, vector2: Vector3D) -> 'Basis3D':
        """
        Create a basis with first vector1 adimensionned, as u, v is the vector2 substracted of u component,
        w is the cross product of u and v
        """
        u = vector1.copy()
        u.normalize()
        v = vector2 - vector2.dot(vector1) * vector1
        v.normalize()
        w = u.cross(v)

        return Basis3D(u, v, w)

    def to_frame(self, origin):
        return Frame3D(origin, self.u, self.v, self.w)

    def rotation(self, axis:Vector3D, angle:float):
        """Rotates the Basis and return new Basis as a result"""

        center = O3D
        new_u = self.u.rotation(center, axis, angle)
        new_v = self.v.rotation(center, axis, angle)
        new_w = self.w.rotation(center, axis, angle)

        return Basis3D(new_u, new_v, new_w, self.name)

    def rotation_inplace(self, axis: Vector3D, angle: float):
        """Rotates the basis and changes is parameters inplace"""
        center = O3D
        new_u = self.u.rotation(center, axis, angle)
        new_v = self.v.rotation(center, axis, angle)
        new_w = self.w.rotation(center, axis, angle)


        self.u = new_u
        self.v = new_v
        self.w = new_w

    def x_rotation(self, angle:float):
        """
        Basis Rotation around the X axis and a new
        Basis is returned as results
        :param angle: rotation angle
        """
        new_u = self.u.x_rotation(angle)
        new_v = self.v.x_rotation(angle)
        new_w = self.w.x_rotation(angle)
        return Basis3D(new_u, new_v, new_w, self.name)


    def x_rotation_inplace(self, angle: float):
        """
        Basis Rotation around the X axis and its
        parameters are changed inplace
        :param angle: rotation angle
        """
        self.u = self.u.x_rotation(angle)
        self.v = self.v.x_rotation(angle)
        self.w = self.w.x_rotation(angle)

    def y_rotation(self, angle:float):
        """
        Basis Rotation around the Y axis and a new
        Basis is returned as results
        :param angle: rotation angle
        """
        new_u = self.u.y_rotation(angle)
        new_v = self.v.y_rotation(angle)
        new_w = self.w.y_rotation(angle)
        return Basis3D(new_u, new_v, new_w, self.name)

    def y_rotation_inplace(self, angle):
        """
        Basis Rotation around the Y axis and its
        parameters are changed inplace
        :param angle: rotation angle
        """
        self.u = self.u.y_rotation(angle)
        self.v = self.v.y_rotation(angle)
        self.w = self.w.y_rotation(angle)

    def z_rotation(self, angle:float):
        """
        Basis Rotation around the Z axis and a new
        Basis is returned as results
        :param angle: rotation angle
        """
        new_u = self.u.z_rotation(angle)
        new_v = self.v.z_rotation(angle)
        new_w = self.w.z_rotation(angle)
        return Basis3D(new_u, new_v, new_w, self.name)

    def z_rotation_inplace(self, angle: float):
        """
        Basis Rotation around the Z axis and a its
        parameters are changed inplace
        :param angle: rotation angle
        """
        self.u = self.u.z_rotation(angle)
        self.v = self.v.z_rotation(angle)
        self.w = self.w.z_rotation(angle)

    # def Eulerrotation(self, angles:Tuple[float, float, float], copy:bool=True):
    #     psi, theta, phi = angles
    #     center = O3D
    #
    #     vect_u = self.u.copy()
    #     vect_v = self.v.copy()
    #     vect_w = self.w.copy()
    #
    #     # rotation around w
    #     vect_u.rotation(center, vect_w, psi, False)
    #     vect_v.rotation(center, vect_w, psi, False)
    #
    #     # rotation around v
    #     vect_v.rotation(center, vect_u, theta, False)
    #     vect_w.rotation(center, vect_u, theta, False)
    #
    #     # rotation around w
    #     vect_u.rotation(center, vect_w, phi, False)
    #     vect_v.rotation(center, vect_w, phi, False)
    #
    #     if copy:
    #         return Basis3D(vect_u, vect_v, vect_w)
    #     self.u = vect_u
    #     self.v = vect_v
    #     self.w = vect_w

    def euler_rotation_parameters(self, angles: Tuple[float, float, float]):

        psi, theta, phi = angles
        center = O3D

        vect_u = self.u.copy()
        vect_v = self.v.copy()
        vect_w = self.w.copy()

        # rotation around w
        vect_u.rotation_inplace(center, vect_w, psi)
        vect_v.rotation_inplace(center, vect_w, psi)

        # rotation around v
        vect_v.rotation_inplace(center, vect_u, theta)
        vect_w.rotation_inplace(center, vect_u, theta)

        # rotation around w
        vect_u.rotation_inplace(center, vect_w, phi)
        vect_v.rotation_inplace(center, vect_w, phi)
        return vect_u, vect_v, vect_w

    def euler_rotation(self, angles:Tuple[float, float, float]):
        """
        Basis rotation using euler rotation and
        returns a new basis as a result
        :param angles: psi, theta, phi
        """
        vect_u, vect_v, vect_w = self.euler_rotation_parameters(angles)
        return Basis3D(vect_u, vect_v, vect_w)

    def euler_rotation_inplace(self, angles: Tuple[float, float, float]):
        """
        Basis rotation using euler rotation and
        its parameters are changed inplace
        :param angles: psi, theta, phi
        """
        vect_u, vect_v, vect_w = self.euler_rotation_parameters(angles)
        self.u = vect_u
        self.v = vect_v
        self.w = vect_w

    def transfer_matrix(self):
        return Matrix33(self.u.x, self.v.x, self.w.x,
                        self.u.y, self.v.y, self.w.y,
                        self.u.z, self.v.z, self.w.z)

    def inverse_transfer_matrix(self):
        return self.transfer_matrix().inverse()

    def new_coordinates(self, vector):
        matrix = self.inverse_transfer_matrix()
        return matrix.vector_multiplication(vector)

    def old_coordinates(self, point):
        matrix = self.transfer_matrix()
        return matrix.vector_multiplication(point)

    def copy(self, deep=True, memo=None):
        return Basis3D(self.u, self.v, self.w)

    def normalize(self):
        """
        normalize the basis modifying its coordinates in place
        """
        self.u.normalize()
        self.v.normalize()
        self.w.normalize()


class Frame2D(Basis2D):
    """
    Defines a 2D basis
    :param origin:Point2D: origin of the basis
    :param u:Vector2D: first vector of the basis
    :param v:Vector2D: second vector of the basis
    """

    def __init__(self, origin: Point2D, u: Vector2D, v: Vector2D, name: Text = ''):
        self.origin = origin
        Basis2D.__init__(self, u, v, name=name)

    def __repr__(self):
        return '{}: O={} U={}, V={}'.format(self.__class__.__name__, self.origin, self.u, self.v)

    def __neg__(self):
        Pinv = self.inverse_transfer_matrix()
        new_origin = Point2D(npy.dot(Pinv, self.origin))
        return Frame2D(new_origin,
                       Vector2D(Pinv[:, 0]),
                       Vector2D(Pinv[:, 1]))

    def __add__(self, other_frame):
        P1 = self.transfer_matrix()
        new_origin = P1.vector_multiplication(other_frame.origin) + self.origin
        M = P1 * other_frame.transfer_matrix()
        return Frame2D(new_origin,
                       Vector2D(M.M11, M.M21),
                       Vector2D(M.M12, M.M22))
        # new_origin = Point2D(npy.dot(P1, other_frame.origin) + self.origin)
        # M = npy.dot(P1, other_frame.transfer_matrix())
        # return Frame2D(new_origin,
        #                Vector2D(M[:, 0]),
        #                Vector2D(M[:, 1]))

    def __sub__(self, other_frame):
        P1inv = other_frame.inverse_transfer_matrix()
        P2 = self.transfer_matrix()
        new_origin = Point2D(npy.dot(P1inv, (self.origin - other_frame.origin)))
        M = npy.dot(P1inv, P2)
        return Frame2D(new_origin,
                       Vector2D(M[:, 0]),
                       Vector2D(M[:, 1]))

    def to_dict(self, *args, **kwargs):
        return {'object_class': 'volmdlr.Frame2D',
                'name': self.name,
                'origin': self.origin.to_dict(),
                'u': self.u.to_dict(),
                'v': self.v.to_dict()
                }

    def basis(self):
        return Basis2D(self.u, self.v)

    def new_coordinates(self, vector):
        return Basis2D.new_coordinates(self, vector - self.origin)

    def old_coordinates(self, vector):
        return Basis2D.old_coordinates(self, vector) + self.origin

    def translation(self, vector):
        """
        Frame Translation
        :param vector: translation vector
        :return: a new translated Frame
        """
        new_origin = self.origin.translation(vector)
        return Frame2D(new_origin, self.u, self.v)

    def translation_inplace(self, vector):
        """
        Frame Translation. Object updated inplace
        :param vector: translation vector
        """
        self.origin = self.origin.translation(vector)

    def rotation(self, angle):
        """
        Frame rotation
        :param angle: rotation angle
        :return: New rotated Frame
        """
        new_base = Basis2D.rotation(self, angle)
        new_frame = Frame2D(self.origin, new_base.u, new_base.v)
        return new_frame

    def rotation_inplace(self, angle: float):
        """
        Frame rotation. Object updated inplace
        :param angle: rotation angle
        """
        new_base = Basis2D.rotation(self, angle)
        self.u = new_base.u
        self.v = new_base.v

    def Draw(self, ax=None, style='ok'):
        if ax is None:
            fig, ax = plt.subplots()

        ax.plot(*self.origin, style)
        self.u.plot(origin=self.origin, ax=ax, color='r')
        self.v.plot(origin=self.origin, ax=ax, color='g')
        ax.axis('equal')

    def copy(self, deep=True, memo=None):
        return Frame2D(self.origin, self.u, self.v)


OXY = Frame2D(O2D, X2D, Y2D)


class Frame3D(Basis3D):
    """
    Defines a 3D frame
    :param origin:Point3D: origin of the basis
    :param u:Vector3D: first vector of the basis
    :param v:Vector3D: second vector of the basis
    :param w:Vector3D: third vector of the basis
    """

    def __init__(self, origin: Point3D, u: Vector3D, v: Vector3D, w: Vector3D, name: Text = ''):
        self.origin = origin
        Basis3D.__init__(self, u, v, w)
        self.name = name

    def __repr__(self):
        return '{}: O={} U={}, V={}, W={}'.format(self.__class__.__name__,
                                                  self.origin,
                                                  self.u, self.v, self.w)

    def __hash__(self):
        """
        hash returns 0 because points are difficult to hash if they are meant
        to be equalized at a given tolerance
        """
        return 0

    def __eq__(self, other_frame):
        if other_frame.__class__.__name__ != self.__class__.__name__:
            return False

        for other_vector, vector in zip([other_frame.origin, other_frame.u,
                                         other_frame.v, other_frame.w],
                                        [self.origin, self.u, self.v, self.w]):
            if other_vector != vector:
                return False
        return True

    def __neg__(self):
        M = self.inverse_transfer_matrix()
        new_origin = M.vector_multiplication(self.origin)
        return Frame3D(new_origin,
                       Vector3D(M.M11, M.M21, M.M31),
                       Vector3D(M.M12, M.M22, M.M32),
                       Vector3D(M.M13, M.M23, M.M33))

    def __add__(self, other_frame):
        P1 = self.transfer_matrix()
        new_origin = P1.vector_multiplication(other_frame.origin) + self.origin

        M = P1 * other_frame.transfer_matrix()
        return Frame3D(new_origin,
                       Vector3D(M.M11, M.M21, M.M31),
                       Vector3D(M.M12, M.M22, M.M32),
                       Vector3D(M.M13, M.M23, M.M33))

    def __sub__(self, other_frame):
        P1inv = other_frame.inverse_transfer_matrix()
        P2 = self.transfer_matrix()
        new_origin = P1inv.vector_multiplication(self.origin - other_frame.origin)
        M = P1inv * P2
        return Frame3D(new_origin,
                       Vector3D(M.M11, M.M21, M.M31),
                       Vector3D(M.M12, M.M22, M.M32),
                       Vector3D(M.M13, M.M23, M.M33))

    def __round__(self, ndigits=6):
        return self.__class__(round(self.origin, ndigits),
                              round(self.u, ndigits),
                              round(self.v, ndigits),
                              round(self.w, ndigits))

    def __hash__(self):
        """
        hash returns 0 because points are difficult to hash if they are meant
        to be equalized at a given tolerance
        """
        return 0

    def to_dict(self, *args, **kwargs):
        return {'object_class': 'volmdlr.Frame3D',
                'name': self.name,
                'origin': self.origin.to_dict(),
                'u': self.u.to_dict(),
                'v': self.v.to_dict(),
                'w': self.w.to_dict()
                }

    # @classmethod
    # def dict_to_object(cls, dict_, global_dict=None,
    #                    pointers_memo: Dict[str, Any] = None, path: str = '#'):
    #     return Frame3D(Point3D.dict_to_object(dict_['origin']),
    #                    Vector3D.dict_to_object(dict_['u']),
    #                    Vector3D.dict_to_object(dict_['v']),
    #                    Vector3D.dict_to_object(dict_['w']),
    #                    dict_.get('name', ''))

    def basis(self):
        return Basis3D(self.u, self.v, self.w)

    def new_coordinates(self, vector):
        """ You have to give coordinates in the global landmark """
        return Basis3D.new_coordinates(self, vector - self.origin)

    def old_coordinates(self, vector):
        """ You have to give coordinates in the local landmark """
        return Basis3D.old_coordinates(self, vector) + self.origin

    def rotation(self, center: Point3D, axis: Vector3D, angle: float):
        """
        Rotate the center as a point and vectors as directions
        (calling Basis), and returns a new Frame
        """
        new_base = Basis3D.rotation(self, axis, angle)
        new_origin = self.origin.rotation(center, axis, angle)
        return Frame3D(new_origin,
                       new_base.u, new_base.v, new_base.w,
                       self.name)

    def rotation_inplace(self, center: Point3D,  axis: Vector3D, angle: float):
        """
        Rotate the center as a point and vectors as directions (calling Basis).
        Object is updated inplace
        """
        new_base = Basis3D.rotation(self, axis, angle)
        new_origin = self.origin.rotation(center, axis, angle)
        self.origin = new_origin
        self.u = new_base.u
        self.v = new_base.v
        self.w = new_base.w

    def translation(self, offset):
        """
        Frame translation
        :param offset: translation vector
        :return: new translated frame
        """
        return Frame3D(self.origin.translation(offset),
                       self.u, self.v, self.w, self.name)

    def translation_inplace(self, offset: Vector3D):
        """
        Frame translation. Object is updated is inplace
        :param offset: translation vector
        """
        self.origin.translation_inplace(offset)

    def copy(self, deep=True, memo=None):
        return Frame3D(self.origin.copy(),
                       self.u.copy(), self.v.copy(), self.w.copy())

    def to_step(self, current_id):

        content, origin_id = self.origin.to_point().to_step(current_id)
        current_id = origin_id + 1
        u_content, u_id = Vector3D.to_step(self.u, current_id)
        current_id = u_id + 1
        v_content, v_id = Vector3D.to_step(self.v, current_id)
        current_id = v_id + 1
        content += u_content + v_content
        content += "#{} = AXIS2_PLACEMENT_3D('{}',#{},#{},#{});\n"\
            .format(current_id, self.name, origin_id, u_id, v_id)
        return content, current_id

    def plot2d(self, x=X3D, y=Y3D, ax=None, color='k'):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        origin2d = self.origin.to_2d(O3D, x, y)

        for iv, vector in enumerate(self.vectors):
            vector2D = vector.to_2d(O3D, x, y)
            if vector2D.norm() > 1e-8:
                vector2D.plot(origin=origin2d, ax=ax, color=color, label=str(iv + 1))

        return fig, ax

    def plot(self, ax=None, color='b', alpha=1., plot_points=True,
             ratio=1):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        x1 = [p.x for p in (self.origin, self.origin + self.u * ratio)]
        y1 = [p.y for p in (self.origin, self.origin + self.u * ratio)]
        z1 = [p.z for p in (self.origin, self.origin + self.u * ratio)]
        ax.plot(x1, y1, z1, 'r')

        x2 = [p.x for p in (self.origin, self.origin + self.v * ratio)]
        y2 = [p.y for p in (self.origin, self.origin + self.v * ratio)]
        z2 = [p.z for p in (self.origin, self.origin + self.v * ratio)]
        ax.plot(x2, y2, z2, 'g')

        x3 = [p.x for p in (self.origin, self.origin + self.w * ratio)]
        y3 = [p.y for p in (self.origin, self.origin + self.w * ratio)]
        z3 = [p.z for p in (self.origin, self.origin + self.w * ratio)]
        ax.plot(x3, y3, z3, 'b')
        return ax

    @classmethod
    def from_step(cls, arguments, object_dict):
        origin = object_dict[arguments[1]]
        if arguments[2] == '$':
            u = None
        else:
            u = object_dict[arguments[2]]
        if arguments[3] == '$':
            v = None
        else:
            v = object_dict[arguments[3]]
        if u is None or v is None:
            w = None
        else:
            w = u.cross(v)

        return cls(origin, u, v, w, arguments[0][1:-1])

    @classmethod
    def from_point_and_vector(cls, point: Point3D, vector: Vector3D, main_axis: Vector3D = X3D):
        """
        Create a new frame from a point and vector by rotating the global frame.
        Global frame rotate in order to have 'vector' and 'main_axis' collinear.
        This method is very useful to compute a local frame of an object.

        :param point: origin of the new frame
        :param vector: vector used to define one of the main axis (by default X-axis) of the local frame
        :param main_axis: the axis of global frame you want to match 'vector' (can be X3D, Y3D or Z3D).
        :return: created local frame
        """
        if main_axis not in [X3D, Y3D, Z3D]:
            raise ValueError("main_axis must be X, Y or Z of the global frame")

        vector.normalize()

        if vector == main_axis:
            # The local frame is oriented like the global frame
            return cls(O3D, X3D, Y3D, Z3D)

        if vector == -main_axis:
            return cls(O3D, -X3D, -Y3D, -Z3D)

        # The local frame is oriented differently from the global frame
        # Rotation angle
        dot = main_axis.dot(vector)
        rot_angle = math.acos(dot / (vector.norm() * main_axis.norm()))

        # Rotation axis
        vector2 = vector - main_axis
        rot_axis = main_axis.cross(vector2)
        rot_axis.normalize()

        u = X3D.rotation(O3D, rot_axis, rot_angle)
        v = Y3D.rotation(O3D, rot_axis, rot_angle)
        w = Z3D.rotation(O3D, rot_axis, rot_angle)

        return cls(point, u, v, w)

    def babylonjs(self, size=0.1, parent=None):
        s = 'var origin = new BABYLON.Vector3({},{},{});\n'.format(*self.origin)
        s += 'var o_u = new BABYLON.Vector3({}, {}, {});\n'.format(*(size * self.u + self.origin))
        s += 'var o_v = new BABYLON.Vector3({}, {}, {});\n'.format(*(size * self.v + self.origin))
        s += 'var o_w = new BABYLON.Vector3({}, {}, {});\n'.format(*(size * self.w + self.origin))
        s += 'var line1 = BABYLON.MeshBuilder.CreateTube("frame_U", {{path: [origin, o_u], radius: {}}}, scene);'.format(
            0.03 * size)
        s += 'line1.material = red_material;\n'
        s += 'var line2 = BABYLON.MeshBuilder.CreateTube("frame_V", {{path: [origin, o_v], radius: {}}}, scene);'.format(
            0.03 * size)
        s += 'line2.material = green_material;\n'
        s += 'var line3 = BABYLON.MeshBuilder.CreateTube("frame_W", {{path: [origin, o_w], radius: {}}}, scene);'.format(
            0.03 * size)
        s += 'line3.material = blue_material;\n'
        if parent is not None:
            s += 'line1.parent = {};\n'.format(parent)
            s += 'line2.parent = {};\n'.format(parent)
            s += 'line3.parent = {};\n'.format(parent)

        return s
