#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# cython: language_level=3
# cython: c_string_type=str, c_string_encoding=ascii
"""

Cython functions

"""
import cython
import cython.cimports.libc.math as math_c
from cython.parallel import prange
import math
import random
import sys
import warnings
from typing import List, Text, Tuple

import matplotlib.pyplot as plt
import numpy as npy
cimport numpy as np
import plot_data
import volmdlr
from dessia_common.core import DessiaObject
from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d


# C functions
# =============================================================================
cdef(double, double) c_sub_2d(double u1, double u2, double v1, double v2):
    return (u1 - v1, u2 - v2)


cdef(double, double) c_add_2d(double u1, double u2, double v1, double v2,):
    return (u1 + v1, u2 + v2)


cdef(double, double) c_mul_2d(double u1, double u2, double value):
    return (u1 * value, u2 * value)


cdef double c_vector2d_dot(double u1, double u2, double v1, double v2):
    return u1 * v1 + u2 * v2


cdef double c_vector2d_norm(double u1, double u2):
    return (u1 * u1 + u2 * u2)**0.5


cdef(double, double, double) c_sub_3d(double u1, double u2, double u3,
                                      double v1, double v2, double v3):
    return (u1 - v1, u2 - v2, u3 - v3)


cdef(double, double, double) c_add_3d(double u1, double u2, double u3,
                                      double v1, double v2, double v3):
    return (u1 + v1, u2 + v2, u3 + v3)


cdef(double, double, double) c_mul_3d(double u1, double u2, double u3, double value):
    return (u1 * value, u2 * value, u3 * value)


cdef double c_vector3d_dot(double u1, double u2, double u3,
                           double v1, double v2, double v3):
    return u1 * v1 + u2 * v2 + u3 * v3


cdef double c_vector3d_norm(double u1, double u2, double u3):
    return (u1 * u1 + u2 * u2 + u3 * u3)**0.5


cdef(double, double, double) c_vector3d_cross(double u1, double u2, double u3,
                                              double v1, double v2, double v3):
    return (u2 * v3 - u3 * v2, u3 * v1 - u1 * v3, u1 * v2 - u2 * v1)


cdef(double, double, double) c_vector3d_rotation(double vx, double vy, double vz,
                                                 double center_x, double center_y, double center_z,
                                                 double axis_x, double axis_y, double axis_z,
                                                 double angle):

    cdef double ux = vx - center_x
    cdef double uy = vy - center_y
    cdef double uz = vz - center_z

    cdef double cos_angle = math_c.cos(angle)
    cdef double sin_angle = math_c.sin(angle)

    cdef double rv1_x = cos_angle * ux
    cdef double rv1_y = cos_angle * uy
    cdef double rv1_z = cos_angle * uz
    cdef double rv2_x, rv2_y, rv2_z, rv3_x, rv3_y, rv3_z
    rv2_x, rv2_y, rv2_z = c_mul_3d(axis_x, axis_y, axis_z,
                                   (1 - cos_angle) * c_vector3d_dot(ux, uy, uz, axis_x, axis_y, axis_z))

    rv3_x, rv3_y, rv3_z = c_vector3d_cross(axis_x, axis_y, axis_z, ux, uy, uz)

    return (rv1_x + rv2_x + rv3_x * sin_angle + center_x,
            rv1_y + rv2_y + rv3_y * sin_angle + center_y,
            rv1_z + rv2_z + rv3_z * sin_angle + center_z)


cdef(double, double, double) c_matrix_vector_multiplication3(double M11, double M12, double M13,
                                                             double M21, double M22, double M23,
                                                             double M31, double M32, double M33,
                                                             double v1, double v2, double v3):

    return (M11 * v1 + M12 * v2 + M13 * v3,
            M21 * v1 + M22 * v2 + M23 * v3,
            M31 * v1 + M32 * v2 + M33 * v3)

cdef(double, double) c_matrix_vector_multiplication2(double M11, double M12,
                                                     double M21, double M22,
                                                     double v1, double v2):

    return (M11 * v1 + M12 * v2,
            M21 * v1 + M22 * v2)


cdef(double, double, double,
     double, double, double,
     double, double, double) c_matrix_multiplication3(double A11, double A12, double A13,
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


cdef(double, (double, double)) c_linesegment2d_point_distance((double, double) p1,
                                                              (double, double) p2, (double, double) point):
    cdef double t, ux, uy, ppx, ppy, vx, vy
    cdef (double, double) projection

    ux, uy = c_sub_2d(p2[0], p2[1], p1[0], p1[1])
    ppx, ppy = c_sub_2d(point[0], point[1], p1[0], p1[1])

    t = max(0, min(1, c_vector2d_dot(ppx, ppy, ux, uy) / c_vector2d_norm(ux, uy)**2))
    vx, vy = c_mul_2d(ux, uy, t)

    projection = c_add_2d(p1[0], p1[1], vx, vy)
    ppx, ppy = projection[0] - point[0], projection[1] - point[1]
    return c_vector2d_norm(ppx, ppy), projection


cdef (double, (double, double, double)) c_linesegment3d_point_distance((double, double, double) p1,
                                                                       (double, double, double) p2,
                                                                       (double, double, double) point):
    cdef double t, ux, uy, uz, ppx, ppy, ppz, vx, vy, vz
    cdef (double, double, double) projection
    ux, uy, uz = c_sub_3d(p2[0], p2[1], p2[2], p1[0], p1[1], p1[2])
    ppx, ppy, ppz = c_sub_3d(point[0], point[1], point[2], p1[0], p1[1], p1[2])
    t = max(0, min(1, c_vector3d_dot(ppx, ppy, ppz, ux, uy, uz) / c_vector3d_norm(ux, uy, uz)**2))
    vx, vy, vz = c_mul_3d(ux, uy, uz, t)
    projection = c_add_3d(p1[0], p1[1], p1[2], vx, vy, vz)
    ppx, ppy, ppz = projection[0]-point[0], projection[1]-point[1], projection[2]-point[2]
    return c_vector3d_norm(ppx, ppy, ppz), projection


# =============================================================================

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef bint polygon_point_belongs(double[:, ::1] polygon, double[:] point,
                                 bint include_edge_points=False, double tol= 1e-6):
    cdef int i
    cdef int n = polygon.shape[0]
    cdef bint inside = False
    cdef double x, y, p1x, p1y, p2x, p2y, xints, dot_product, length_squared, t, distance_projection_to_point
    cdef double[2] u, v, projection_vector, projection_point
    x, y = point
    for i in range(n):
        p1x = polygon[i][0]
        p1y = polygon[i][1]
        p2x = polygon[(i + 1) % n][0]
        p2y = polygon[(i + 1) % n][1]
        v = [p2x - p1x, p2y - p1y]
        u = [x - p1x, y - p1y]
        dot_product = u[0] * v[0] + u[1] * v[1]
        length_squared = v[0] ** 2 + v[1] ** 2
        t = (dot_product / length_squared)
        if 0.0 <= t <= 1.0:
            projection_vector = [v[0] * t, v[1] * t]
            projection_point = [p1x + projection_vector[0], p1y + projection_vector[1]]
            distance_projection_to_point = math_c.sqrt((projection_point[0] - x) ** 2 + (projection_point[1] - y) ** 2)
            if distance_projection_to_point <= tol:
                if include_edge_points:
                    return True
                return False
        xints = math_c.HUGE_VAL
        if min(p1y, p2y) <= y <= max(p1y, p2y) and min(p1x, p2x) <= x <= max(p1x, p2x):
            if p1y != p2y:
                xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
            if p1y == p2y or x == xints:
                if include_edge_points:
                    return True
                return False
        if min(p1y, p2y) < y <= max(p1y, p2y) and x <= max(p1x, p2x):
            if p1y != p2y:
                xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
            if p1x == p2x or x < xints:
                inside = not inside
    return inside


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef np.ndarray[np.uint8_t, ndim = 1] points_in_polygon(double[:, ::1] polygon, double[:, ::1] points,
                                                         bint include_edge_points = False, double tol = 1e-6):
    cdef int n = polygon.shape[0]
    cdef int m = points.shape[0]
    cdef int i, j
    cdef double x, y, p1x, p1y, p2x, p2y, xints, dot_product, length_squared, t, distance_projection_to_point
    cdef double[2] u, v, projection_vector, projection_point
    cdef np.ndarray[np.uint8_t, ndim = 1] results = npy.zeros(m, dtype=npy.uint8)
    cdef bint inside

    for i in prange(m, nogil=True):
        x = points[i][0]
        y = points[i][1]
        inside = False
        for j in range(n):
            p1x = polygon[j][0]
            p1y = polygon[j][1]
            p2x = polygon[(j + 1) % n][0]
            p2y = polygon[(j + 1) % n][1]
            v[0] = p2x - p1x
            v[1] = p2y - p1y
            u[0] = x - p1x
            u[1] = y - p1y
            dot_product = u[0] * v[0] + u[1] * v[1]
            length_squared = v[0] * v[0] + v[1] * v[1]
            t = dot_product / length_squared
            if 0.0 <= t <= 1.0:
                projection_vector[0] = v[0] * t
                projection_vector[1] = v[1] * t
                projection_point[0] = p1x + projection_vector[0]
                projection_point[1] = p1y + projection_vector[1]
                distance_projection_to_point = math_c.sqrt((projection_point[0] - x) ** 2 + (projection_point[1]
                                                                                             - y) ** 2)
                if distance_projection_to_point <= tol:
                    if include_edge_points:
                        results[i] = True
                        break
                    else:
                        results[i] = False
                        break
            xints = math_c.HUGE_VAL
            if min(p1y, p2y) <= y <= max(p1y, p2y) and min(p1x, p2x) <= x <= max(p1x, p2x):
                if p1y != p2y:
                    xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                if p1y == p2y or x == xints:
                    if include_edge_points:
                        results[i] = True
                        break
                    else:
                        results[i] = False
                        break
            if min(p1y, p2y) < y <= max(p1y, p2y) and x <= max(p1x, p2x):
                if p1y != p2y:
                    xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                if p1x == p2x or x < xints:
                    inside = not inside
        results[i] = inside

    return results


# =============================================================================
def bbox_is_intersecting(bbox1, bbox2, tol):
    """Verifies if the two bounding boxes are intersecting, or touching."""
    cdef double x1_min, x1_max, y1_min, y1_max, z1_min, z1_max, x2_min, x2_max, y2_min, y2_max, z2_min, z2_max
    x1_min = bbox1.xmin - tol
    x1_max = bbox1.xmax + tol
    y1_min = bbox1.ymin - tol
    y1_max = bbox1.ymax + tol
    z1_min = bbox1.zmin - tol
    z1_max = bbox1.zmax + tol

    x2_min = bbox2.xmin - tol
    x2_max = bbox2.xmax + tol
    y2_min = bbox2.ymin - tol
    y2_max = bbox2.ymax + tol
    z2_min = bbox2.zmin - tol
    z2_max = bbox2.zmax + tol

    # Check for non-intersection cases
    if (x1_max < x2_min or x1_min > x2_max or
            y1_max < y2_min or y1_min > y2_max or
            z1_max < z2_min or z1_min > z2_max):
        return False

    return True


def linesegment2d_point_distance(tuple start, tuple end, tuple point):
    return c_linesegment2d_point_distance(start, end, point)


def linesegment3d_point_distance(tuple start, tuple end, tuple point):
    return c_linesegment3d_point_distance(start, end, point)


cdef (double, double) CLineSegment3DDistance((double, double, double) u,
                                             (double, double, double) v,
                                             (double, double, double) w):
    cdef double a = c_vector3d_dot(u[0], u[1], u[2], u[0], u[1], u[2])
    cdef double b = c_vector3d_dot(u[0], u[1], u[2], v[0], v[1], v[2])
    cdef double c = c_vector3d_dot(v[0], v[1], v[2], v[0], v[1], v[2])
    cdef double d = c_vector3d_dot(u[0], u[1], u[2], w[0], w[1], w[2])
    cdef double e = c_vector3d_dot(v[0], v[1], v[2], w[0], w[1], w[2])
    cdef double determinant = a * c - b * c
    cdef double s_parameter
    cdef double t_parameter
    if determinant > - 1e-6:
        b_times_e = b * e
        c_times_d = c * d
        if b_times_e <= c_times_d:
            s_parameter = 0.0
            if e <= 0.0:
                t_parameter = 0.0
                negative_d = -d
                if negative_d >= a:
                    s_parameter = 1.0
                elif negative_d > 0.0:
                    s_parameter = negative_d / a
            elif e < c:
                t_parameter = e / c
            else:
                t_parameter = 1.0
                b_minus_d = b - d
                if b_minus_d >= a:
                    s_parameter = 1.0
                elif b_minus_d > 0.0:
                    s_parameter = b_minus_d / a
        else:
            s_parameter = b_times_e - c_times_d
            if s_parameter >= determinant:
                s_parameter = 1.0
                b_plus_e = b + e
                if b_plus_e <= 0.0:
                    t_parameter = 0.0
                    negative_d = -d
                    if negative_d <= 0.0:
                        s_parameter = 0.0
                    elif negative_d < a:
                        s_parameter = negative_d / a
                elif b_plus_e < c:
                    t_parameter = b_plus_e / c
                else:
                    t_parameter = 1.0
                    b_minus_d = b - d
                    if b_minus_d <= 0.0:
                        s_parameter = 0.0
                    elif b_minus_d < a:
                        s_parameter = b_minus_d / a
            else:
                a_times_e = a * e
                b_times_d = a * d
                if a_times_e <= b_times_d:
                    t_parameter = 0.0
                    negative_d = -d
                    if negative_d <= 0.0:
                        s_parameter = 0.0
                    elif negative_d >= a:
                        s_parameter = 1.0
                    else:
                        s_parameter = negative_d / a
                else:
                    t_parameter = a_times_e - b_times_d
                    if t_parameter >= determinant:
                        t_parameter = 1.0
                        b_minus_d = b - d
                        if b_minus_d <= 0.0:
                            s_parameter = 0.0
                        elif b_minus_d >= a:
                            s_parameter = 1.0
                        else:
                            s_parameter = b_minus_d / a
                    else:
                        s_parameter /= determinant
                        t_parameter /= determinant
    else:
        if e <= 0.0:
            t_parameter = 0.0
            negative_d = -d
            if negative_d <= 0.0:
                s_parameter = 0.0
            elif negative_d >= a:
                s_parameter = 1.0
            else:
                s_parameter = negative_d / a
        elif e >= c:
            t_parameter = 1.0
            b_minus_d = b - d
            if b_minus_d <= 0.0:
                s_parameter = 0.0
            elif b_minus_d >= a:
                s_parameter = 1.0
            else:
                s_parameter = b_minus_d / a
        else:
            s_parameter = 0.0
            t_parameter = e / c
    return s_parameter, t_parameter


def LineSegment3DDistance(points_linesegment1, points_linesegment2):
    u = c_sub_3d(points_linesegment1[1].x, points_linesegment1[1].y, points_linesegment1[1].z,
                 points_linesegment1[0].x, points_linesegment1[0].y, points_linesegment1[0].z)
    v = c_sub_3d(points_linesegment2[1].x, points_linesegment2[1].y, points_linesegment2[1].z,
                 points_linesegment2[0].x, points_linesegment2[0].y, points_linesegment2[0].z)
    w = c_sub_3d(points_linesegment1[0].x, points_linesegment1[0].y, points_linesegment1[0].z,
                 points_linesegment2[0].x, points_linesegment2[0].y, points_linesegment2[0].z)
    s_parameter, t_parameter = CLineSegment3DDistance(u, v, w)
    point1 = points_linesegment1[0] + Vector3D(*u) * s_parameter
    point2 = points_linesegment2[0] + Vector3D(*v) * t_parameter
    return point1, point2
# =============================================================================
#  Points, Vectors
# =============================================================================


class Arrow3D(FancyArrowPatch):
    def __init__(self, x, y, z, starting_point=None, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
        if starting_point is None:
            starting_point = (0., 0., 0.)

        self.starting_point = starting_point
        self._xyz = x, y, z

    def draw(self, renderer):
        """
        Draw the arrow by overloading a Matplotlib method. Do not rename this method!

        """
        x1, y1, z1 = self.starting_point
        dx, dy, dz = self._xyz
        x2, y2, z2 = (dx + x1, dy + y1, dz + z1)
        xn, yn, zn = proj3d.proj_transform((x1, x2), (y1, y2), (z1, z2), self.axes.M)
        self.set_positions((xn[0], yn[0]), (xn[1], yn[1]))
        FancyArrowPatch.draw(self, renderer)

    def plot(self, ax=None, color="b"):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")

        points = [self.start, self.end]
        x = [p.x for p in points]
        y = [p.y for p in points]
        z = [p.z for p in points]
        ax.plot(x, y, z, "o-k")
        return ax


cdef class Vector:
    """
    Abstract class of vector
    """
    def __init__(self, name = ""):
        self.name = name

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

    def is_colinear_to(self, other_vector: Vector, abs_tol: float = 1e-6):
        """
        Checks if two vectors are colinear.
        The two vectors should be of same dimension.

        :param other_vector: A vector-like object
        :type other_vector: :class:`volmdlr.Vector`
        :param abs_tol: Absolute tolerance to consider colinear
        :type abs_tol: float
        :return: `True` if the two vectors are colinear, `False` otherwise
        :rtype: bool
        """
        try:
            return math.isclose(abs(self.dot(other_vector)) / self.norm() / other_vector.norm(),
                                1,
                                abs_tol=abs_tol)

        except ZeroDivisionError:
            return False

    def is_perpendicular_to(self, other_vector: Vector, abs_tol: float = 1e-5):
        """
        Checks if two vectors are perpendicular.
        The two vectors should be of same dimension.

        :param other_vector: A vector-like object
        :type other_vector: :class:`volmdlr.Vector`
        :param abs_tol: Absolute tolerance to consider perpendicular
        :type abs_tol: float
        :return: `True` if the two vectors are perpendicular, `False` otherwise
        :rtype: bool
        """
        return math.isclose(abs(self.dot(other_vector)), 0, abs_tol=abs_tol)

    @classmethod
    def mean_point(cls, points: List["Vector"], name = ""):
        """
        Find the mean point from a list of points. All the objects of this list
        should be of same dimension.

        :param points: A list of vector-like objects
        :type points: List[:class:`volmdlr.Vector`]
        :param name: object's name.
        :return: The mean point or vector
        :rtype: :class:`volmdlr.Vector`
        """
        n = 1
        point = points[0].copy()
        for point2 in points[1:]:
            point += point2
            n += 1
        point /= n
        point.name = name
        return point

    def vector_projection(self, other_vector):
        """
        Projects the vector onto other_vector.

        :param other_vector: Vector to project self.
        """
        return (self.dot(other_vector) / other_vector.dot(other_vector)) * other_vector

    @classmethod
    def remove_duplicate(cls, points: List["Vector"]):
        """
        An approximative method to remove duplicated points from a list.
        All the objects of this list should be of same dimension.

        :param points: A list of vector-like objects with potential duplicates
        :type points: List[:class:`volmdlr.Vector`]
        :return: The new list of vector-like objects without duplicates&
        :rtype: List[:class:`volmdlr.Vector`]
        """
        dict_ = {p.approx_hash(): p for p in points}
        return list(dict_.values())

    def to_vector(self):
        return self


cdef class Vector2D(Vector):
    """
    Class representing a 2-dimensional vector.

    :param x: The vector's abscissa
    :type x: float
    :param y: The vector's ordinate
    :type y: float

    """
    def __init__(self, x: float, y: float, name: str = ""):
        self.x = x
        self.y = y
        Vector.__init__(self, name=name)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.x}, {self.y})"

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
        return Vector2D(*c_add_2d(self.x, self.y, other_vector.x, other_vector.y))

    def __neg__(self):
        return Vector2D(-self.x, -self.y)

    def __sub__(self, other_vector):
        return Vector2D(*c_sub_2d(self.x, self.y, other_vector.x, other_vector.y))

    def __mul__(self, value: float):
        return Vector2D(*c_mul_2d(self.x, self.y, value))

    def __truediv__(self, value: float):
        if value == 0:
            raise ZeroDivisionError
        return Vector2D(self.x / value,
                        self.y / value)

    def __round__(self, ndigits: int = 6):
        return self.__class__(round(self.x, ndigits), round(self.y, ndigits))

    def __hash__(self):
        """Return a hash value for the vector."""
        return hash(("vector", self.x, self.y))

    def __eq__(self, other):
        """Return True if the other point has the same x and y coordinates, False otherwise."""
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        return False

    def _data_eq(self, other):
        return self == other

    def is_close(self, other_vector: Vector2D, tol: float = 1e-6):
        """
        Checks if two vectors are close to each other considering the
        Euclidean distance. The tolerance can be modified. The two vectors
        should be of same dimension.

        :param other_vector: A Vector2D-like object
        :type other_vector: :class:`volmdlr.Vector2D`
        :param tol: The tolerance under which the euclidean distance is
            considered equal to 0
        :type tol: float
        :return: `True` if the two Vector2D-like objects are close enough
            to each other, `False` otherwise
        :rtype: bool
        """
        if other_vector.__class__.__name__ not in ["Vector2D", "Point2D", "Node2D"]:
            return False
        return math.isclose(self.point_distance(other_vector), 0, abs_tol=tol)

    def approx_hash(self):
        """
        Computes an approximative hash value based on the coordinates.

        :return: An approximative hash value
        :rtype: int
        """
        return round(1e6 * (self.x + self.y))

    def to_dict(self, *args, **kwargs):
        """
        Serializes a 2-dimensional vector into a dictionary.

        :return: A serialized version of the Vector2D
        :rtype: dict

        .. seealso::
            How `serialization and deserialization`_ works in dessia_common

        .. _serialization and deserialization:
            https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method
        """
        return {"object_class": "volmdlr.Vector2D",
                "x": self.x, "y": self.y,
                "name": self.name}

    @classmethod
    def dict_to_object(cls, dict_, *args, **kwargs):
        """
        Deserializes a dictionary to a 3-dimensional point.

        :param dict_: The dictionary of a serialized Point3D
        :type dict_: dict
        :return:
        :rtype: :class:`volmdlr.Point3D`

        .. seealso::
            How `serialization and deserialization`_ works in dessia_common

        .. _serialization and deserialization:
            https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method
        """
        return cls(dict_["x"], dict_["y"], dict_.get("name", ""))

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of a 2-dimensional vector.

        :param deep: *not used*
        :param memo: *not used*
        :return: A copy of the Vector2D-like object
        :rtype: :class:`volmdlr.Vector2D`
        """
        return self.__class__(self.x, self.y)

    def norm(self):
        """
        Computes the euclidiean norm of a 2-dimensional vector.

        :return: Norm of the Vector2D-like object
        :rtype: float
        """
        return c_vector2d_norm(self.x, self.y)

    def unit_vector(self):
        """Calculates the unit vector."""
        n = self.norm()
        if n == 0:
            raise ZeroDivisionError
        return Vector2D(self.x / n, self.y / n)

    def dot(self, other_vector: "Vector2D"):
        """
        Computes the dot product (scalar product) of two 2-dimensional vectors.

        :param other_vector: A Vector2D-like object
        :type other_vector: :class:`volmdlr.Vector2D`
        :return: A scalar, result of the dot product
        :rtype: float
        """
        return c_vector2d_dot(self.x, self.y, other_vector.x, other_vector.y)

    def cross(self, other_vector: "Vector2D"):
        """
        Computes the cross product of two 2-dimensional vectors.

        :param other_vector: A Vector2D-like object
        :type other_vector: :class:`volmdlr.Vector2D`
        :return: A scalar, result of the cross product
        :rtype: float
        """
        return self.x * other_vector.y - self.y * other_vector.x

    def point_distance(self, other_vector: "Vector2D"):
        """
        Computes the euclidiean distance between two Vector2D objects.

        :param other_vector: A Vector2D object
        :type other_vector: :class:`volmdlr.Vector2D`
        :return: The euclidiean distance
        :rtype: float
        """
        return (self - other_vector).norm()

    def rotation_parameters(self, center: "Point2D", angle: float):
        """
        Calculates the parameters to be used in rotation methods

        :param center: The center of rotation
        :type center: :class:`volmdlr.Point2D`
        :param angle: The angle of the rotation in radian
        :type angle: float
        :return: The abscissa and ordinate of the rotated vector
        :rtype: tuple
        """
        u = self - center
        v2x = math.cos(angle) * u.x - math.sin(angle) * u.y + center.x
        v2y = math.sin(angle) * u.x + math.cos(angle) * u.y + center.y
        return v2x, v2y

    def rotation(self, center: "Point2D", angle: float):
        """
        Rotates the 2-dimensional vector and returns a new rotated vector

        :param center: The center of rotation
        :type center: :class:`volmdlr.Point2D`
        :param angle: The angle of the rotation in radian
        :type angle: float
        :return: A rotated Vector2D-like object
        :rtype: :class:`volmdlr.Vector2D`
        """
        v2x, v2y = self.rotation_parameters(center, angle)
        return self.__class__(v2x, v2y)

    def translation(self, offset: "Vector2D"):
        """
        Translates the 2-dimensional vector and returns a new translated vector

        :param offset: The offset vector of the translation
        :type offset: :class:`volmdlr.Vector2D`
        :return: A translated Vector2D-like object
        :rtype: :class:`volmdlr.Vector2D`
        """
        v2x = self.x + offset.x
        v2y = self.y + offset.y
        return self.__class__(v2x, v2y)

    def frame_mapping(self, frame: "Frame2D", side: str):
        """
        # TODO: Needs correction. Add an example ?
        Transforms a 2-dimensional vector from the current reference frame to a
        new one. Choose side equals to 'old' if the current reference frame is
        the old one ; choose side equals to 'new' if the input reference frame
        is the new one. This way, choosing 'old' will return the frame mapped
        vector of the input reference frame.

        :param frame: The input reference frame
        :type frame: :class:`volmdlr.Frame2D`
        :param side: Choose between 'old' and 'new'
        :type side: str
        :return: A frame mapped Vector2D-like object
        :rtype: :class:`volmdlr.Vector2D`
        """
        if side == "old":
            new_vector = frame.local_to_global_coordinates(self)
        if side == "new":
            new_vector = frame.global_to_local_coordinates(self)
        return new_vector

    def to_3d(self, plane_origin: Point3D, vx: Vector3D, vy: Vector3D):
        """
        Returns the 3-dimensional vector corresponding to the 2-dimensional
        vector placed on the 3-dimensional plane (XY) of the 3-dimensional
        frame (centered on `plane_origin`, having for basis (`vx`, `vy`, vz),
        vz being the cross product of `vx` and `vy`).

        :param plane_origin: The origin of the plane, on which lies the
            Vector2D
        :type plane_origin: :class:`volmdlr.Vector3D`
        :param vx: The first direction of the plane
        :type vx: :class:`volmdlr.Vector3D`
        :param vy: The second direction of the plane
        :type vy: :class:`volmdlr.Vector3D`
        :return: The Vector3D from the Vector2D set in the 3-dimensional space
        :rtype: :class:`volmdlr.Vector3D`
        """
        return Vector3D(plane_origin.x + vx.x * self.x + vy.x * self.y,
                        plane_origin.y + vx.y * self.x + vy.y * self.y,
                        plane_origin.z + vx.z * self.x + vy.z * self.y)

    def to_point(self):
        """
        Transforms a Vector2D into a Point2D and returns it.

        :return: A Point2D
        :rtype: :class:`volmdlr.Point2D`
        """
        return Point2D(self.x, self.y)

    def normal_vector(self):
        """
        Returns the normal vector located pi/2 (counterclockwise) to the
        2-dimensional vector.

        :return: A normal Vector2D
        :rtype: :class:`volmdlr.Vector2D`
        """
        return Vector2D(-self.y, self.x)

    def unit_normal_vector(self):
        """
        Returns the unit normal vector located pi/2 (counterclockwise) to the
        2-dimensional vector.

        :return: A unit normal Vector2D
        :rtype: :class:`volmdlr.Vector2D`
        """
        n = self.normal_vector()
        return n.unit_vector()

    def deterministic_unit_normal_vector(self):
        """
        # TODO: to be deleted ?
        # TODO: Or unit_normal_vector should be renamed deterministic_unit_normal_vector ?
        """
        return self.unit_normal_vector()

    @classmethod
    def random(cls, xmin: float, xmax: float, ymin: float, ymax: float, name: str = ""):
        """
        Returns a random 2-dimensional point.

        :param xmin: The minimal abscissa
        :type xmin: float
        :param xmax: The maximal abscissa
        :type xmax: float
        :param ymin: The minimal ordinate
        :type ymin: float
        :param ymax: The maximal ordinate
        :type ymax: float
        :param name: object's name.
        :return: A random Vector2D
        :rtype: :class:`volmdlr.Vector2D`
        """
        return cls(random.uniform(xmin, xmax),
                   random.uniform(ymin, ymax), name=name)

    def plot(self,
             head_width: float = 3, origin: "Vector2D" = None,
             ax: "matplotlib.axes.Axes" = None,
             color: str = "k", label: str = None):
        """
        Plots the 2-dimensional vector. If the vector has a norm greater than
        1e-9, it will be plotted with an arrow, else it will be plotted with
        a point.

        :param head_width: The width of the head of the arrow
        :type head_width: float, optional
        :param origin: The starting point of the tail of the arrow
        :type origin: :class:`volmdlr.Vector2D`, optional
        :param ax: The Axes on which the Vector2D will be drawn
        :type ax: :class:`matplotlib.axes.Axes`, optional
        :param color: The color of the arrow
        :type color: str, optional
        :param label: The text you want to display
        :type label: str, optional
        :return: A matplotlib Axes object on which the Vector2D have been
            plotted
        :rtype: :class:`matplotlib.axes.Axes`
        """
        if origin is None:
            origin = Vector2D(0., 0.)

        if ax is None:
            fig, ax = plt.subplots()

        if math.isclose(self.norm(), 0, abs_tol=1e-9):
            point = origin.copy()
            point.plot(ax=ax, color=color)
            return ax

        ax.quiver(origin[0], origin[1], self[0], self[1], angles="xy", scale_units="xy",
                  scale=1, color=color)

        ax.set(xlim=sorted([self.x, origin.x]), ylim=sorted([self.y, origin.y]))

        if label is not None:
            ax.text(*(origin + self * 0.5), label)

        return ax

    def to_step(self, current_id, vector=False, vertex=False):
        """
        Write a step primitive from a 2-dimensional vector.

        :param current_id: The id of the last written primitive
        :type current_id: int
        :param vector: If 'True' creates a step VECTOR primitive. Otherwise,
            only a DIRECTION primitive will be created. Default value set to
            'False'
        :type vector: bool, optional
        :param vertex: If 'True' calls the to_step method of Point3D. Default
            value set to 'False'
        :type vertex: bool, optional
        :return: A tuple containing the string representing the step primitive
            and the new current id
        :rtype: tuple
        """
        if vertex:
            return self.to_point().to_step(current_id=current_id, vertex=True)
        content = f"#{current_id} = DIRECTION('{self.name}',({self.x:.6f},{self.y:.6f}));\n"
        if vector:
            content += f"#{current_id + 1} = VECTOR('{self.name}',#{current_id},1.);\n"
            current_id += 1
        return content, current_id


X2D = Vector2D(1, 0)
Y2D = Vector2D(0, 1)


cdef class Point2D(Vector2D):
    """
    Class representing a 2-dimensional point.

    :param x: The vector's abscissa
    :type x: float
    :param y: The vector's ordinate
    :type y: float
    :param name: The vector's name
    :type name: str
    """

    def __init__(self, x, y, name = ""):
        self.x = x
        self.y = y
        self.name = name

    def __add__(self, other_vector):
        return Point2D(*c_add_2d(self.x, self.y, other_vector.x, other_vector.y))

    def __neg__(self):
        return Point2D(-self.x, -self.y)

    def __sub__(self, other_vector):
        return Point2D(*c_sub_2d(self.x, self.y, other_vector.x, other_vector.y))

    def __mul__(self, value: float):
        return Point2D(*c_mul_2d(self.x, self.y, value))

    def __truediv__(self, value: float):
        if value == 0:
            raise ZeroDivisionError
        return Point2D(self.x / value,
                       self.y / value)

    def __hash__(self):
        """Return a hash value for the point 2d."""
        return hash(("point", self.x, self.y))

    def __eq__(self, other):
        """Return True if the other point has the same x and y coordinates, False otherwise."""
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y
        return False

    def _data_eq(self, other):
        return self == other

    def to_dict(self, *args, **kwargs):
        """
        Serializes a 2-dimensional point into a dictionary.

        :return: A serialized version of the Point2D
        :rtype: dict

        .. seealso::
            How `serialization and deserialization`_ works in dessia_common

        .. _serialization and deserialization:
            https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method
        """
        dict_ = {"object_class": "volmdlr.Point2D",
                 "x": self.x, "y": self.y}
        if self.name:
            dict_["name"] = self.name
        return dict_

    def to_3d(self, plane_origin: "Vector3D", vx: "Vector3D", vy: "Vector3D"):
        """
        Returns the 3-dimensional point corresponding to the 2-dimensional
        point placed on the 3-dimensional plane (XY) of the 3-dimensional
        frame (centered on `plane_origin`, having for basis (`vx`, `vy`, vz),
        vz being the cross product of `vx` and `vy`).

        :param plane_origin: The origin of the plane, on which lies the
            Vector2D
        :type plane_origin: :class:`volmdlr.Vector3D`
        :param vx: The first direction of the plane
        :type vx: :class:`volmdlr.Vector3D`
        :param vy: The second direction of the plane
        :type vy: :class:`volmdlr.Vector3D`
        :return: The Point3D from the Point2D set in the 3-dimensional space
        :rtype: :class:`volmdlr.Point3D`
        """
        return Point3D(round(plane_origin.x + vx.x * self.x + vy.x * self.y, 12),
                       round(plane_origin.y + vx.y * self.x + vy.y * self.y, 12),
                       round(plane_origin.z + vx.z * self.x + vy.z * self.y, 12))

    def to_vector(self):
        """
        Transforms a Point2D into a Vector2D and returns it.

        :return: A Vector2D
        :rtype: :class:`volmdlr.Vector2D`
        """
        return Vector2D(self.x, self.y)

    def to_step(self, current_id, vertex=False):
        content = "#{} = CARTESIAN_POINT('{}',({:.6f},{:.6f}));\n"\
                        .format(current_id, self.name,
                                1000.*self.x,
                                1000.*self.y)
        if vertex:
            content += "#{} = VERTEX_POINT('{}',#{});\n".format(current_id+1,
                                                                self.name,
                                                                current_id)
            current_id += 1

        return content, current_id

    def plot(self, ax=None, color="k", alpha=1, plot_points=True):
        """
        Plots the 2-dimensional point as a dot.

        :param ax: The Axes on which the Vector2D will be drawn
        :type ax: :class:`matplotlib.axes.Axes`, optional
        :param color: The color of the arrow
        :type color: str, optional
        :param alpha: The transparency of the point from 0 to 1. 0 being
            fully transparent
        :type alpha: float, optional
        :param plot_points: # TODO: delete this attribute
        :type plot_points: bool, optional
        :return: A matplotlib Axes object on which the Point2D have been plotted
        :rtype: :class:`matplotlib.axes.Axes`
        """
        if ax is None:
            fig, ax = plt.subplots()

        ax.plot([self.x], [self.y], color=color, alpha=alpha, marker="o")
        return ax

    def point_distance(self, other_point: "Point2D"):
        """
        Computes the euclidiean distance between two Point2D objects.

        :param other_point: A Point2D object
        :type other_point: :class:`volmdlr.Point2D`
        :return: The euclidiean distance
        :rtype: float
        """
        return (self - other_point).norm()

    @classmethod
    def line_intersection(cls, line1: "volmdlr.edges.Line2D",
                          line2: "volmdlr.edges.Line2D",
                          curvilinear_abscissa: bool = False, name: str = ""):
        """
        Returns a Point2D based on the intersection between two infinite lines.

        :param line1: The first line
        :type line1: :class:`volmdlr.edges.Line2D`
        :param line2: The second line
        :type line2: :class:`volmdlr.edges.Line2D`
        :param curvilinear_abscissa: `True` will return, in addition to the
            intersection point, the curvilinear abscissa of the point on the
            first line and on the second line. Otherwise, only the point will
            be returned
        :type curvilinear_abscissa: bool, optional
        :param name: Object's name.
        :return: The two-dimensional point at the intersection of the two lines
        :rtype: :class:`volmdlr.Point2D`
        """
        (x1, y1), (x2, y2) = line1
        (x3, y3), (x4, y4) = line2

        denominateur = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if math.isclose(denominateur, 0, abs_tol=1e-15):
            if not curvilinear_abscissa:
                return None
            else:
                return None, None, None
        x = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)
        x = x / denominateur
        y = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)
        y = y / denominateur
        if not curvilinear_abscissa:
            return cls(x, y)
        t = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)
        t = t / denominateur
        u = (x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)
        u = -u / denominateur
        return cls(x, y), t, u

    @classmethod
    def segment_intersection(cls, segment1: "volmdlr.edges.LineSegment2D",
                             segment2: "volmdlr.edges.LineSegment2D",
                             curvilinear_abscissa: bool = False, name: str = ""):
        """
        Returns a Point2D based on the intersection between two finite lines.

        :param segment1: The first line segment
        :type segment1: :class:`volmdlr.edges.LineSegment2D`
        :param segment2: The second line segment
        :type segment2: :class:`volmdlr.edges.LineSegment2D`
        :param curvilinear_abscissa: `True` will return, in addition to the
            intersection point, the curvilinear abscissa of the point on the
            first line segment and on the second line segment. Otherwise, only
            the point will be returned
        :type curvilinear_abscissa: bool, optional
        :param name: object's name.
        :return: The two-dimensional point at the intersection of the two lines
            segments
        :rtype: :class:`volmdlr.Point2D`
        """
        (x1, y1), (x2, y2) = segment1
        (x3, y3), (x4, y4) = segment2

        denominateur = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
        if math.isclose(denominateur, 0, abs_tol=1e-15):
            if not curvilinear_abscissa:
                return None
            else:
                return None, None, None

        t = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)
        t = t / denominateur
        u = (x1 - x2) * (y1 - y3) - (y1 - y2) * (x1 - x3)
        u = -u / denominateur
        if (0 <= t <= 1) or (0 <= u <= 1):
            x = (x1 * y2 - y1 * x2) * (x3 - x4) - (x1 - x2) * (x3 * y4 - y3 * x4)
            x = x / denominateur
            y = (x1 * y2 - y1 * x2) * (y3 - y4) - (y1 - y2) * (x3 * y4 - y3 * x4)
            y = y / denominateur
            if not curvilinear_abscissa:
                return cls(x, y, name=name)
            else:
                return cls(x, y, name=name), t, u
        else:
            if not curvilinear_abscissa:
                return None
            else:
                return None, None, None

    def plot_data(self, marker=None, color="black", size=1,
                  opacity=1, arrow=False, stroke_width=None):
        """
        Transforms the two-dimensional point into a plot_data twe-dimensional point.

        :param marker: # TODO: unused parameter
        :type marker: str, optional
        :param color: # TODO: unused parameter
        :type color: str, optional
        :param size: # TODO: unused parameter
        :type size: float, optional
        :param opacity: # TODO: unused parameter
        :type opacity: float, optional
        :param arrow: # TODO: unused parameter
        :type arrow: bool, optional
        :param stroke_width: # TODO: unused parameter
        :type stroke_width: float, optional
        :return: a plot_data two-dimensional point
        :rtype: :class:`plot_data.Point2D`
        """
        return plot_data.Point2D(self.x, self.y)

    @classmethod
    def middle_point(cls, point1: Vector2D, point2: Vector2D, name = ""):
        """
        Computes the middle point between two two-dimensional vector-like objects.

        :param point1: the first point
        :type point1: :class:`volmdlr.Vector2D`
        :param point2: the second point
        :type point2: :class:`volmdlr.Vector2D`.
        :param name: object's name.
        :return: the middle point
        :rtype: :class:`volmdlr.Point2D`
        """
        middle_point = (point1 + point2) * 0.5
        middle_point.name = name
        return middle_point

    @classmethod
    def line_projection(cls, point: Vector2D,
                        line: "volmdlr.edges.Line2D", name = ""):
        """
        Computes the projection of a two-dimensional vector-like object on an
        infinite two-dimensional line

        :param point: the point to be projected
        :type point: :class:`volmdlr.Vector2D`
        :param line: the infinite line
        :type line: :class:`volmdlr.edges.Line2D`.
        :param name: object's name.
        :return: the projected point
        :rtype: :class:`volmdlr.Point2D`
        """
        p1, p2 = line[0], line[1]
        n = line.unit_normal_vector()
        pp1 = point - p1
        point = pp1 - pp1.dot(n) * n + p1
        point.name = name
        return point

    def nearest_point(self, points: List[Vector2D]):
        """
        Finds the nearest point out of a list of two-dimensional vector-like
        objects.

        :param points: a list of points
        :type points: List[:class:`volmdlr.Vector2D`]
        :return: the nearest point out of the list
        :rtype: :class:`volmdlr.Vector2D`
        """
        min_distance = self.point_distance(points[0])
        min_point = points[0]
        for point in points:
            pd = self.point_distance(point)
            if pd < min_distance:
                min_distance, min_point = pd, point
        return min_point

    def axial_symmetry(self, line: "volmdlr.curves.Line2D"):
        """
        Returns the symmetric two-dimensional point according to a line.

        :param line: the line used for axial symmetry
        :type line: :class:`volmdlr.edges.Line2D`
        :return: the symmetrical point
        :rtype: :class:`volmdlr.Point2D`
        """

        point_projection = line.point_projection(self)[0]
        point_symmetry = point_projection + (point_projection - self)

        return point_symmetry

    def coordinates(self):
        """
        Gets x,y coordinates of a point2d.
        """

        return (self.x, self.y)

    def get_geo_lines(self, tag: int, point_mesh_size: float = None):
        """
        Gets the lines that define a Point2D in a .geo file.

        :param tag: The point index
        :type tag: int
        :param mesh_size: The target mesh size close to the point, defaults to None
        :type mesh_size: float, optional

        :return: A line
        :rtype: str
        """

        if point_mesh_size:
            return "Point("+str(tag)+") = {"+str([*self, 0])[1:-1]+", "+str(point_mesh_size)+"};"
        else:
            return "Point("+str(tag)+") = {"+str([*self, 0])[1:-1]+"};"


O2D = Point2D(0, 0)


cdef class Vector3D(Vector):
    """
    Class representing a 3-dimensional vector.

    :param x: The vector's abscissa
    :type x: float
    :param y: The vector's ordinate
    :type y: float
    :param Z: The vector's applicate
    :type z: float
    :param name: The vector's name
    :type name: str
    """

    def __init__(self, double x, double y, double z, name: str = ""):
        self.x = x
        self.y = y
        self.z = z
        Vector.__init__(self, name=name)

    def __repr__(self):
        return f"{self.__class__.__name__}({self.x}, {self.y}, {self.z})"

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
        return Vector3D(*c_add_3d(self.x, self.y, self.z, other_vector.x, other_vector.y, other_vector.z))

    def __neg__(self):
        return Vector3D(-self.x, -self.y, -self.z)

    def __sub__(self, other_vector):
        return Vector3D(*c_sub_3d(self.x, self.y, self.z, other_vector.x, other_vector.y, other_vector.z))

    def __mul__(self, value):
        return Vector3D(*c_mul_3d(self.x, self.y, self.z, value))

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
        """Return a hash value for the vector 3d."""
        return hash(("vector", self.x, self.y, self.z))

    def __eq__(self, other):
        """Return True if the other point has the same x, y and z coordinates, False otherwise."""
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y and self.z == other.z
        return False

    def _data_eq(self, other):
        return self == other

    def is_close(self, other_vector, tol=1e-6):
        """
        Checks if two vectors are close to each other considering the
        Euclidean distance. The tolerance can be modified. The two vectors
        should be of same dimension.

        :param other_vector: A Vector3D-like object
        :type other_vector: :class:`volmdlr.Vector3D`
        :param tol: The tolerance under which the euclidean distance is
            considered equal to 0
        :type tol: float
        :return: `True` if the two Vector3D-like objects are close enough
            to each other, `False` otherwise
        :rtype: bool
        """
        if other_vector.__class__.__name__ not in ["Vector3D", "Point3D", "Node3D"]:
            return False
        return math.isclose(self.point_distance(other_vector), 0, abs_tol=tol)

    def approx_hash(self):
        """
        Computes an approximative hash value based on the coordinates.

        :return: An approximative hash value
        :rtype: int
        """
        return round(1e6 * (self.x + self.y + self.z))

    def to_dict(self, *args, **kwargs):
        """
        Serializes a 3-dimensional vector into a dictionary.

        :return: A serialized version of the Vector3D
        :rtype: dict

        .. seealso::
            How `serialization and deserialization`_ works in dessia_common

        .. _serialization and deserialization:
            https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method
        """
        dict_ = {"object_class": "volmdlr.Vector3D",
                 "x": self.x, "y": self.y, "z": self.z}

        if self.name:
            dict_["name"] = self.name
        return dict_

    @classmethod
    def dict_to_object(cls, dict_, *args, **kwargs):
        """
        Deserializes a dictionary to a 3-dimensional vector.

        :param dict_: The dictionary of a serialized Vector3D
        :type dict_: dict
        :param global_dict: The global dictionary. Default value is None
        :type global_dict: dict, optional
        :param pointers_memo: A dictionary from path to python object of
            already serialized values. Default value is None
        :type pointers_memo: dict, optional
        :param path: The path in the global object. In most cases, append
            ‘/attribute_name’ to given path for your attributes.
            Default value is '#'
        :type path: str, optional
        :return:
        :rtype: :class:`volmdlr.Vector3D`

        .. seealso::
            How `serialization and deserialization`_ works in dessia_common

        .. _serialization and deserialization:
            https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method
        """

        return cls(dict_["x"], dict_["y"], dict_["z"], dict_.get("name", ""))

    def dot(self, other_vector):
        """
        Computes the dot product between two 3-dimensional vectors.

        :param other_vector: A Vector3D-like object
        :type other_vector: :class:`volmdlr.Vector3D`
        :return: Value of the dot product
        :rtype: float
        """
        return c_vector3d_dot(self.x, self.y, self.z,
                              other_vector.x, other_vector.y, other_vector.z)

    def cross(self, other_vector: "Vector3D") -> "Vector3D":
        """
        Computes the cross product between two 3-dimensional vectors.

        :param other_vector: A Vector3D-like object
        :type other_vector: :class:`volmdlr.Vector3D`
        :return: Value of the cross product
        :rtype: float
        """
        return self.__class__(*c_vector3d_cross(self.x, self.y, self.z,
                                                other_vector.x, other_vector.y, other_vector.z))

    def norm(self) -> float:
        """
        Computes the euclidiean norm of a 3-dimensional vector.

        :return: Norm of the Vector3D-like object
        :rtype: float
        """
        return c_vector3d_norm(self.x, self.y, self.z)

    def unit_vector(self):
        """Calculates the unit vector."""
        n = self.norm()
        if n == 0:
            raise ZeroDivisionError
        return Vector3D(self.x / n, self.y / n, self.z / n)

    def point_distance(self, point2: "Vector3D") -> float:
        """
        Computes the euclidiean distance between two Vector3D objects.

        :param other_vector: A Vector3D object
        :type other_vector: :class:`volmdlr.Vector3D`
        :return: The euclidiean distance
        :rtype: float
        """
        return (self - point2).norm()

    def rotation(self, center: "Point3D", axis: "Vector3D", angle: float):
        """
        Rotates of angle around axis the 2-dimensional vector and returns
        a new rotated vector.
        Using Rodrigues Formula:
            https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula.

        :param center: The center of rotation
        :type center: :class:`volmdlr.Point3D`
        :param axis: The axis of rotation
        :type axis: :class:`volmdlr.Vector3D`
        :param angle: The angle of the rotation in radian
        :type angle: float
        :return: A rotated Vector3D-like object
        :rtype: :class:`volmdlr.Vector3D`
        """
        vector2 = c_vector3d_rotation(self.x, self.y, self.z,
                                      center.x, center.y, center.z,
                                      axis.x, axis.y, axis.z, angle)
        return self.__class__(*vector2)

    @staticmethod
    def axis_rotation_parameters(axis1_value, axis2_value, angle):
        """
        # TODO: to be completed
        Calcules new axis1 and axis2 new values after vector rotation.

        :param axis1_value:
        :type axis1_value:
        :param axis2_value:
        :type axis2_value:
        :param angle:
        :type angle:
        :return:
        :rtype: tuple
        """
        cos_angle = math.cos(angle)
        sin_angle = math.sin(angle)

        axis1 = cos_angle * axis1_value + sin_angle * axis2_value
        axis2 = -sin_angle * axis1_value + cos_angle * axis2_value
        return axis1, axis2

    def x_rotation(self, angle: float):
        """
        Rotation of angle around X axis and returns a new vector as a result.

        :param angle: Value of the angle
        :type angle: float
        :return: A 3-dimensional point
        :rtype: :class:`volmdlr.Point3D`
        """
        y1, z1 = self.axis_rotation_parameters(self.y, self.z, angle)

        return Point3D(self.x, y1, z1)

    def y_rotation(self, angle: float):
        """
        Rotation of angle around Y axis and returns a new vector as result.

        :param angle: Value of the angle
        :type angle: float
        :return: A 3-dimensional point
        :rtype: :class:`volmdlr.Point3D`
        """
        z1, x1 = self.axis_rotation_parameters(self.z, self.x, angle)
        return Point3D(x1, self.y, z1)

    def z_rotation(self, angle: float):
        """
        rrotation of angle around Z axis and returns a new vector as result.

        :param angle: Value of the angle
        :type angle: float
        :return: A 3-dimensional point
        :rtype: :class:`volmdlr.Point3D`
        """
        x1, y1 = self.axis_rotation_parameters(self.x, self.y, angle)
        return Point3D(x1, y1, self.z)

    def translation(self, offset: "Vector3D"):
        """
        Translates the vector and returns a new translated vector

        :param offset: A Vector3D-like object used for offsetting
        :type offset: :class:`volmdlr.Vector3D`
        :return: A translated Vector3D-like object
        :rtype: :class:`volmdlr.Vector3D`
        """
        return self + offset

    def frame_mapping(self, frame: "Frame3D", side: str):
        """
        # TODO: Needs correction. Add an example ?
        Transforms a 3-dimensional vector from the current reference frame to a
        new one. Choose side equals to 'old' if the current reference frame is
        the old one ; choose side equals to 'new' if the input reference frame
        is the new one. This way, choosing 'old' will return the frame mapped
        vector of the input reference frame.

        :param frame: The input reference frame
        :type frame: :class:`volmdlr.Frame3D`
        :param side: Choose between 'old' and 'new'
        :type side: str
        :return: A frame mapped Vector3D-like object
        :rtype: :class:`volmdlr.Vector3D`
        """
        if side == "old":
            new_vector = frame.local_to_global_coordinates(self)

        if side == "new":
            new_vector = frame.global_to_local_coordinates(self)
        return new_vector

    def plane_projection3d(self, plane_origin: Vector3D, x: Vector3D, y: Vector3D):
        """
        Projects a Vector3D-like object on a 3D plane.

        :param plane_origin: The origin of the 3D projection plane
        :type plane_origin: :class:`volmdlr.Vector3D`
        :param x: The X axis of the 3D plane
        :type x: :class:`volmdlr.Vector3D`
        :param y: The Y axis of the 3D plane
        :type y: :class:`volmdlr.Vector3D`
        :return: The projection on the 3D plane
        :rtype: :class:`volmdlr.Vector3D`
        """
        z = x.cross(y)
        z = z.unit_vector()
        return self - z.dot(self - plane_origin) * z

    def plane_projection2d(self, plane_origin: Vector3D, x: Vector3D, y: Vector3D):
        """
        Projects a Vector3D-like object on a 2D plane.

        :param plane_origin: The 3D origin of the 2D projection plane
        :type plane_origin: :class:`volmdlr.Vector3D`
        :param x: The 3D X axis of the 2D plane
        :type x: :class:`volmdlr.Vector3D`
        :param y: The 3D Y axis of the 2D plane
        :type y: :class:`volmdlr.Vector3D`
        :return: The projection on the 2D plane
        :rtype: :class:`volmdlr.Point2D`
        """
        p3d = self.plane_projection3d(plane_origin, x, y)
        u1 = p3d.dot(x)
        u2 = p3d.dot(y)
        return Point2D(u1, u2)

    def to_2d(self, plane_origin: Point3D, x: Vector3D, y: Vector3D):
        """
        # TODO: difference with plane_projection2d needs details
        Transforms a Vector3D-like object to a Point2D.

        :param plane_origin: The 3D origin of the 2D projection plane
        :type plane_origin: :class:`volmdlr.Vector3D`
        :param x: The 3D X axis of the 2D plane
        :type x: :class:`volmdlr.Vector3D`
        :param y: The 3D Y axis of the 2D plane
        :type y: :class:`volmdlr.Vector3D`
        :return: The transformed Point2D
        :rtype: :class:`volmdlr.Point2D`
        """
        x2d = self.dot(x) - plane_origin.dot(x)
        y2d = self.dot(y) - plane_origin.dot(y)
        class_name = self.__class__.__name__[:-2] + "2D"
        if class_name in ("Vector2D", "Point2D"):
            return getattr(sys.modules["volmdlr.core_compiled"], class_name)(x2d, y2d)
        return getattr(sys.modules["volmdlr.display"], class_name)(x2d, y2d)

    def random_unit_normal_vector(self):
        """
        Returns a random normal 3-dimensional vector.

        :return: A normal Vector3D
        :rtype: :class:`volmdlr.Vector3D`
        """
        v = Vector3D.random(0, 1, 0, 1, 0, 1)

        v = v - v.dot(self) * self / (self.norm()**2)
        return v.unit_vector()

    def deterministic_normal_vector(self):
        """
        Returns a deterministic normal 3-dimensional vector.

        :return: A normal Vector3D
        :rtype: :class:`volmdlr.Vector3D`
        """
        if not math.isclose(self.y, 0, abs_tol=1e-7) \
                or not math.isclose(self.z, 0, abs_tol=1e-7):
            v = X3D
        else:
            v = Y3D
        v = v - v.dot(self) * self / (self.norm()**2)
        return v

    def deterministic_unit_normal_vector(self):
        """
        Returns a deterministic unit normal 3-dimensional vector.

        :return: A normal Vector3D
        :rtype: :class:`volmdlr.Vector3D`
        """
        normal_vector = self.deterministic_normal_vector()
        return normal_vector.unit_vector()

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of a 2-dimensional vector.

        :param deep: *not used*
        :param memo: *not used*
        :return: A copy of the Vector2D-like object
        :rtype: :class:`volmdlr.Vector2D`
        """
        return self.__class__(self.x, self.y, self.z)

    @classmethod
    def random(cls, xmin: float, xmax: float, ymin: float, ymax: float, zmin: float, zmax: float, name: str = ""):
        """
        Returns a random 2-dimensional point.

        :param xmin: The minimal abscissa
        :type xmin: float
        :param xmax: The maximal abscissa
        :type xmax: float
        :param ymin: The minimal ordinate
        :type ymin: float
        :param ymax: The maximal ordinate
        :type ymax: float
        :param zmin: The minimal applicate
        :type zmin: float
        :param zmax: The maximal applicate
        :type zmax: float.
        :param name: object's name.
        :return: A random Vector3D
        :rtype: :class:`volmdlr.Vector3D`
        """
        return cls(random.uniform(xmin, xmax),
                   random.uniform(ymin, ymax),
                   random.uniform(zmin, zmax))

    def to_point(self):
        """
        Converts a Vector3D object to a Point3D object.

        :return: A Point3D
        :rtype: :class:`volmdlr.Point3D`
        """
        return Point3D(self.x, self.y, self.z)

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive from a 3-dimensional vector to a Vector3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated
        :type object_dict: dict
        :return: The corresponding Vector3D object
        :rtype: :class:`volmdlr.Vector3D`
        """
        length_conversion_factor = kwargs.get("length_conversion_factor", 1)

        if type(arguments[1]) is int:
            # VECTOR
            new_vector = length_conversion_factor * float(arguments[2])*object_dict[arguments[1]]
            new_vector.name = arguments[0][1:-1]
            return new_vector
        else:
            # DIRECTION
            # return cls(*[float(i)/1000 for i in arguments[1][1:-1].split(",")],
            #             arguments[0][1:-1])
            return cls(*[float(i) for i in arguments[1][1:-1].split(",")],
                       arguments[0][1:-1])

    def to_step(self, current_id, vector=False, vertex=False):
        """
        Write a step primitive from a 3-dimensional vector.

        :param current_id: The id of the last written primitive
        :type current_id: int
        :param vector: If 'True' creates a step VECTOR primitive. Otherwise,
            only a DIRECTION primitive will be created. Default value set to
            'False'
        :type vector: bool, optional
        :param vertex: If 'True' calls the to_step method of Point3D. Default
            value set to 'False'
        :type vertex: bool, optional
        :return: A tuple containing the string representing the step primitive
            and the new current id
        :rtype: tuple
        """
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

    def plot(self, ax=None, starting_point=None, color="k"):
        """
        Plots the 3-dimensional vector.

        :param ax: The Axes on which the Vector2D will be drawn
        :type ax: :class:`matplotlib.axes.Axes`, optional
        :param starting_point: The location of the origin of the vector.
            Default value is None, corresponding to (0, 0, 0)
        :type starting_point: :class:`volmdlr.Vector3D`, optional
        :param color: The color of the drawn vector. Default value is empty
            string for black
        :type color: str, optional
        :return: A matplotlib Axes object on which the Vector3D have been
            plotted
        :rtype: :class:`matplotlib.axes.Axes`
        """
        if starting_point is None:
            starting_point = Point3D(0, 0, 0)
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")
        # Change for head length
        arrow = Arrow3D(self.x, self.y, self.z, starting_point=starting_point, mutation_scale=20,
                        arrowstyle="-|>", color=color)

        ax.add_artist(arrow)
        return ax


X3D = Vector3D(1, 0, 0)
Y3D = Vector3D(0, 1, 0)
Z3D = Vector3D(0, 0, 1)


cdef class Point3D(Vector3D):
    """
    Class representing a 3-dimensional point.

    :param x: The vector's abscissa
    :type x: float
    :param y: The vector's ordinate
    :type y: float
    :param z: The vector's applicate
    :type y: float
    :param name: The vector's name
    :type name: str
    """

    def __init__(self, x, y, z, name = ""):
        self.x = x
        self.y = y
        self.z = z
        self.name = name

    def __add__(self, other_vector):
        return Point3D(*c_add_3d(self.x, self.y, self.z, other_vector.x, other_vector.y, other_vector.z))

    def __neg__(self):
        return Point3D(-self.x, -self.y, -self.z)

    def __sub__(self, other_vector):
        return Point3D(*c_sub_3d(self.x, self.y, self.z, other_vector.x, other_vector.y, other_vector.z))

    def __mul__(self, value):
        return Point3D(*c_mul_3d(self.x, self.y, self.z, value))

    def __truediv__(self, value):
        if value == 0:
            raise ZeroDivisionError
        return Point3D(self.x / value,
                       self.y / value,
                       self.z / value)

    def __hash__(self):
        """Return a hash value for the point 3d."""
        return hash(("point", self.x, self.y, self.z))

    def __eq__(self, other):
        """Return True if the other point has the same x, y and z coordinates, False otherwise."""
        if isinstance(other, self.__class__):
            return self.x == other.x and self.y == other.y and self.z == other.z
        return False

    def _data_eq(self, other):
        return self == other

    def to_dict(self, *args, **kwargs):
        """
        Serializes a 3-dimensional point into a dictionary.

        :return: A serialized version of the Point3D
        :rtype: dict

        .. seealso::
            How `serialization and deserialization`_ works in dessia_common

        .. _serialization and deserialization:
            https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method
        """
        dict_ = {"object_class": "volmdlr.Point3D",
                 "x": self.x, "y": self.y, "z": self.z}
        if self.name:
            dict_["name"] = self.name
        return dict_

    def plot(self, ax=None, color="k", alpha=1, marker="o"):
        """
        Plots the 3-dimensional point.

        :param ax: The Axes on which the Point3D will be drawn. Default value
            is None, creating a new drawing figure
        :type ax: :class:`matplotlib.axes.Axes`, optional
        :param color: The color of the point. Default value is 'k', for black
        :type color: str, optional
        :param alpha: The transparency of the Point3D. Default value is 1, for
            full opacity
        :type alpha: float, optional
        :param marker: The shape of the Point3D. Default value is 'o', for a
            round marker
        :type marker: str, optional
        :return: A matplotlib Axes object on which the Point3D have been
            plotted
        :rtype: :class:`matplotlib.axes.Axes`
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")

        ax.plot([self.x], [self.y], [self.z], color=color, alpha=alpha,
                marker=marker)
        return ax

    # def to_2d(self, plane_origin, x, y):
    #     """
    #     Using Vector3D.to_2d
    #     """
    #     x2d = self.dot(x) - plane_origin.dot(x)
    #     y2d = self.dot(y) - plane_origin.dot(y)
    #     return Point2D(x2d, y2d)

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive from a 3-dimensional point to a Point3D.

        :param arguments: The arguments of the step primitive
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated
        :type object_dict: dict
        :return: The corresponding Point3D object
        :rtype: :class:`volmdlr.Point3D`
        """
        length_conversion_factor = kwargs.get("length_conversion_factor", 1)

        return cls(*[float(i) * length_conversion_factor for i in arguments[1][1:-1].split(",")],
                   arguments[0][1:-1])

    def to_vector(self):
        """
        Converts a Point3D object to a Vector3D object.

        :return: A Vector3D
        :rtype: :class:`volmdlr.Vector3D`
        """
        return Vector3D(self.x, self.y, self.z)

    def point_distance(self, point2: "Point3D") -> float:
        """
        Computes the euclidean distance between two 3-dimensional points.

        :param point2: The other 3-dimensional point
        :type point2: :class:`volmdlr.Point3D`
        :return: The euclidean distance
        :rtype: float
        """
        return (self - point2).norm()

    @classmethod
    def middle_point(cls, point1: "Point3D", point2: "Point3D", name = ""):
        """
        Computes the middle point between two 3-dimensional points.

        :param point1: The first 3-dimensional point
        :type point1: :class:`volmdlr.Point3D`
        :param point2: The second 3-dimensional point
        :type point2: :class:`volmdlr.Point3D`.
        :param name: object's name.
        :return: The middle point
        :rtype: :class:`volmdlr.Point3D`
        """
        middle_point = (point1 + point2) * 0.5
        middle_point.name = name
        return middle_point

    def to_step(self, current_id, vertex=False):
        """
        Writes a step primitive from a 3-dimensional point.

        :param current_id: The id of the last written primitive
        :type current_id: int
        :param vertex: If 'True', adds a VERTEX_POINT step primitive on top of
            the CARTESIAN_POINT step primitive. Default value set to 'False'
        :type vertex: bool, optional
        :return: A tuple containing the string representing the step primitive
            and the new current id
        :rtype: tuple
        """
        current_id += 1
        content = "#{} = CARTESIAN_POINT('{}',({:.6f},{:.6f},{:.6f}));\n".format(current_id, self.name,
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
        """
        # TODO: to be deleted ?
        Returns the babylonjs script for 3D display in browser.

        :return: A babylonjs script
        :rtype: str
        """
        s = 'var sphere = BABYLON.MeshBuilder.CreateSphere("point", {diameter: 0.05}, scene);\n'
        s += "sphere.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n".format(self.x, self.y, self.z)
        s += 'var mat = new BABYLON.StandardMaterial("mat", scene);\n'
        s += "mat.diffuseColor = new BABYLON.Color3(1, 0, 0);\n"
        s += "sphere.material = mat;\n"
        return s

    def nearest_point(self, points: List["Point3D"]):
        """
        Returns the nearest 3-dimensional point out of the list.

        :param points: A list of 3-dimensional points
        :type points: List[:class:`volmdlr.Point3D`]
        :return: The closest point
        :rtype: :class:`volmdlr.Point3D`
        """
        min_distance, closest_point = math.inf, None
        for point in points:
            distance = self.point_distance(point)
            if distance < min_distance:
                closest_point = point
                min_distance = distance
        return closest_point

    def coordinates(self):
        """
        Returns the coordinates of a Point3D as a tuple of values.

        :return: A tuple containing the abscissan, the ordiante and the
            applicate of the Point3D
        :rtype: tuple
        """
        return self.x, self.y, self.z

    def get_geo_lines(self, tag: int, point_mesh_size: float = None):
        """
        Gets the lines that define a Point3D in a .geo file.

        :param tag: The point index
        :type tag: int
        :param mesh_size: The target mesh size close to the point, defaults to None
        :type mesh_size: float, optional

        :return: A line
        :rtype: str
        """

        if point_mesh_size:
            return "Point("+str(tag)+") = {"+str([*self, 0])[1:-1]+", "+str(point_mesh_size)+"};"
        else:
            return "Point("+str(tag)+") = {"+str([*self, 0])[1:-1]+"};"


O3D = Point3D(0, 0, 0)


# =============================================================================
#  Basis, Frames
# =============================================================================

class Matrix22:
    """
    Class representing a 2x2 matrix.

    :param M11: The first line, first column value
    :type M11: float
    :param M12: The first line, second column value
    :type M12: float
    :param M21: The second line, first column value
    :type M21: float
    :param M22: The second line, second column value
    :type M22: float
    """

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
        """
        Multiplies the matrix by a 2-dimensional vector.

        :param vector: A Vector2D-like object
        :type vector: :class:`volmdlr.Vector2D`
        :return: A Vector2D-like object
        :rtype: :class:`volmdlr.Vector2D`
        """
        u1, u2 = c_matrix_vector_multiplication2(self.M11, self.M12,
                                                 self.M21, self.M22,
                                                 vector.x, vector.y)

        return vector.__class__(u1, u2)

    def determinent(self):
        """
        Computes the determinent of the matrix.

        :return: The determinent of the matrix
        :rtype: float
        """
        return self.M11 * self.M22 - self.M12 * self.M21

    def inverse(self):
        """
        Computes the invert matrix.

        :return: The inverse of the matrix
        :rtype: :class:`volmdlr.Matrix22`
        """
        det = self.determinent()
        if not math.isclose(det, 0, abs_tol=1e-10):
            det_inv = 1 / self.determinent()
            return Matrix22(det_inv * self.M22, -det_inv * self.M12,
                            -det_inv * self.M21, det_inv * self.M11)
        else:
            raise ValueError("The matrix is singular")

    # def vector_multiplication(self, vector):
    #     return vector.__class__(self.M11 * vector.x + self.M12 * vector.y,
    #                             self.M21 * vector.x + self.M22 * vector.y)


class Matrix33:
    """
        Class representing a 3x3 matrix.

        :param M11: The first line, first column value
        :type M11: float
        :param M12: The first line, second column value
        :type M12: float
        :param M13: The first line, third column value
        :type M13: float
        :param M21: The second line, first column value
        :type M21: float
        :param M22: The second line, second column value
        :type M22: float
        :param M23: The second line, third column value
        :type M23: float
        :param M31: The third line, first column value
        :type M31: float
        :param M32: The third line, second column value
        :type M32: float
        :param M33: The third line, third column value
        :type M33: float
        """

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
        (M11, M12, M13,
         M21, M22, M23,
         M31, M32, M33) = c_matrix_multiplication3(self.M11, self.M12, self.M13,
                                                   self.M21, self.M22, self.M23,
                                                   self.M31, self.M32, self.M33,
                                                   other_matrix.M11, other_matrix.M12, other_matrix.M13,
                                                   other_matrix.M21, other_matrix.M22, other_matrix.M23,
                                                   other_matrix.M31, other_matrix.M32, other_matrix.M33)

        return Matrix33(M11, M12, M13, M21, M22, M23, M31, M32, M33)

    def __repr__(self):
        s = (f"[{self.M11} {self.M12} {self.M13}]\n"
             f"[{self.M21} {self.M22} {self.M23}]\n"
             f"[{self.M31} {self.M32} {self.M33}]\n")
        return s

    def float_multiplication(self, float_value: float):
        """
        Multiplies the whole matrix by a scalar value.

        :param float_value: The value of the scalar
        :type float_value: float
        :return: The new matrix after multiplication
        :rtype: :class:`volmdlr.Matrix33`
        """
        return Matrix33(self.M11 * float_value, self.M12 * float_value, self.M13 * float_value,
                        self.M21 * float_value, self.M22 * float_value, self.M23 * float_value,
                        self.M31 * float_value, self.M32 * float_value, self.M33 * float_value)

    def vector_multiplication(self, vector):
        """
       Multiplies the matrix by a 3-dimensional vector.

       :param vector: A Vector3D-like object
       :type vector: :class:`volmdlr.Vector3D`
       :return: A Vector3D-like object
       :rtype: :class:`volmdlr.Vector3D`
       """
        u1, u2, u3 = c_matrix_vector_multiplication3(self.M11, self.M12, self.M13,
                                                     self.M21, self.M22, self.M23,
                                                     self.M31, self.M32, self.M33,
                                                     vector.x, vector.y, vector.z)
        if abs(u1) < 1e-12:
            u1 = 0.
        if abs(u2) < 1e-12:
            u2 = 0.
        if abs(u3) < 1e-12:
            u3 = 0.
        return vector.__class__(u1, u2, u3)

    def determinent(self):
        """
        Computes the determinent of the matrix.

        :return: The determinent of the matrix
        :rtype: float
        """
        det = self.M11 * self.M22 * self.M33 + self.M12 * self.M23 * self.M31 \
            + self.M13 * self.M21 * self.M32 - self.M13 * self.M22 * self.M31 \
            - self.M23 * self.M32 * self.M11 - self.M33 * self.M12 * self.M21
        return det

    def inverse(self):
        """
        Computes the invert matrix.

        :return: The inverse of the matrix
        :rtype: :class:`volmdlr.Matrix33`
        """
        det = self.determinent()

        if not abs(det) <= 1e-12:
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
            raise ValueError("The matrix is singular")

    @classmethod
    def random_matrix(cls, minimum: float = 0., maximum: float = 1., name: str = ""):
        """
        Creates a random matrix with values between bounds.

        :param minimum: Minimum possible value of matrix coefficients. Default
            value is 0
        :type minimum: float, optional
        :param maximum: Maximum possible value of matrix coefficients. Default
            value is 1
        :type maximum: float, optional.
        :param name: object's name.
        :return: A random matrix
        :rtype: :class:`volmdlr.Matrix33`
        """
        range_ = maximum - minimum
        return cls(*[minimum + range_ * random.random() for _ in range(9)])

    def to_numpy(self):
        """
        Returns the numpy array corresponding to the matrix.

        :return: A numpy array of the matrix
        :rtype: :class:`numpy.array`
        """
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
        return 0

    def copy(self, deep=True, memo=None):
        return self.__class__(*self.vectors)


class Basis2D(Basis):
    """
    Defines a 2D basis.

    :param u: First vector of the basis
    :type u: :class:`volmdlr.Vector2D`
    :param v: Second vector of the basis
    :type v: :class:`volmdlr.Vector2D`
    """

    def __init__(self, u: Vector2D, v: Vector2D, name: Text = ""):
        self.u = u
        self.v = v
        self.name = name

    def __eq__(self, other_basis):
        if other_basis.__class__.__name__ != self.__class__.__name__:
            return False
        return all([other_vector == vector for other_vector, vector in zip([other_basis.u, other_basis.v],
                                                                           [self.u, self.v])])

    def __neg__(self):
        p_inv = self.inverse_transfer_matrix()
        return Basis2D(Vector3D(p_inv[:, 0]),
                       Vector3D(p_inv[:, 1]))

    def __repr__(self):
        return "{}: U={}, V={}".format(self.__class__.__name__, *self.vectors)

    def _get_vectors(self):
        return (self.u, self.v)

    vectors = property(_get_vectors)

    def to_dict(self, *args, **kwargs):
        """
        Serializes a 2-dimensional basis into a dictionary.

        :return: A serialized version of the Basis2D
        :rtype: dict

        .. seealso::
            How `serialization and deserialization`_ works in dessia_common

        .. _serialization and deserialization:
            https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method
        """
        return {"object_class": "volmdlr.Basis2D",
                "name": self.name,
                "u": self.u.to_dict(),
                "v": self.v.to_dict()
                }

    def to_frame(self, origin: Point2D) -> "Frame2D":
        """
        Returns the 2-dimensional frame oriented the same way as the Basis2D
        and having for origin the given 2-dimensional point.

        :param origin: The origin of the 2-dimensional frame
        :type origin: :class:`volmdlr.Point2D`
        :return: A 2-dimensional frame
        :rtype: :class:`volmdlr.Frame2D`
        """
        return Frame2D(origin, self.u, self.v)

    def transfer_matrix(self):
        """
        Computes the transfer matrix of the 2-dimensional basis.

        :return: The 2x2 transfer matrix
        :rtype: :class:`volmdlr.Matrix22`
        """
        return Matrix22(self.u.x, self.v.x,
                        self.u.y, self.v.y)

    def inverse_transfer_matrix(self):
        """
        Computes the inverse transfer matrix of the 2-dimensional basis.

        :return: The 2x2 inverse transfer matrix
        :rtype: :class:`volmdlr.Matrix22`
        """
        return self.transfer_matrix().inverse()

    def new_coordinates(self, vector: Vector2D) -> Vector2D:
        """
        Convert the given vector's coordinates from the global landmark to the local landmark of this Basis2D.

        :param vector: The vector to convert, given in global coordinates.
        :type vector: :class:`volmdlr.Vector2D`
        :return: The converted vector, in local coordinates.
        :rtype: :class:`volmdlr.Vector2D`

        .. deprecated:: Use global_to_local_coordinates instead.
        """
        warnings.warn(
            "new_coordinates is deprecated. Use global_to_local_coordinates instead.",
            DeprecationWarning,
        )
        return self.global_to_local_coordinates(vector)

    def old_coordinates(self, vector: Vector2D) -> Vector2D:
        """
        Convert the given vector's coordinates from the local landmark of this Basis2D to the global landmark.

        :param vector: The vector to convert, given in local coordinates.
        :type vector: :class:`volmdlr.Vector2D`
        :return: The converted vector, in global coordinates.
        :rtype: :class:`volmdlr.Vector2D`

        .. deprecated:: Use local_to_global_coordinates instead.
        """
        warnings.warn(
            "old_coordinates is deprecated. Use local_to_global_coordinates instead.",
            DeprecationWarning,
        )
        return self.local_to_global_coordinates(vector)

    def global_to_local_coordinates(self, vector: Vector2D) -> Vector2D:
        """
        Convert the given vector's coordinates from the global landmark to the local landmark of this Basis2D.

        :param vector: The vector to convert, given in global coordinates.
        :type vector: :class:`volmdlr.Vector2D`
        :return: The converted vector, in local coordinates.
        :rtype: :class:`volmdlr.Vector2D`
        """
        matrix = self.inverse_transfer_matrix()
        return matrix.vector_multiplication(vector)

    def local_to_global_coordinates(self, vector: Vector2D) -> Vector2D:
        """
        Convert the given vector's coordinates from the local landmark of this Basis2D to the global landmark.

        :param vector: The vector to convert, given in local coordinates.
        :type vector: :class:`volmdlr.Vector2D`
        :return: The converted vector, in global coordinates.
        :rtype: :class:`volmdlr.Vector2D`
        """
        matrix = self.transfer_matrix()
        return matrix.vector_multiplication(vector)

    def rotation(self, angle: float):
        """
        Rotates the 2-dimensional basis and returns a new rotated one.

        :param angle: The angle of rotation in rad
        :type angle: float
        :return: The rotated Basis2D
        :rtype: :class:`volmdlr.Basis2D`
        """
        center = O2D
        new_u = self.u.rotation(center, angle)
        new_v = self.v.rotation(center, angle)
        return Basis2D(new_u, new_v)

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of a 2-dimensional basis.

        :param deep: *not used*
        :param memo: *not used*
        :return: A copy of the Basis2D
        :rtype: :class:`volmdlr.Basis2D`
        """
        return Basis2D(self.u, self.v)

    def normalize(self):
        """
        Normalizes the basis, and return a new object.

        :return: A new normalized basis
        :rtype: Basis2D
        """
        return Basis2D(self.u.unit_vector(), self.v.unit_vector())


XY = Basis2D(X2D, Y2D)


class Basis3D(Basis):
    """
    Defines a 3D basis.

    :param u: First vector of the basis
    :type u: :class:`volmdlr.Vector3D`
    :param v: Second vector of the basis
    :type v: :class:`volmdlr.Vector3D`
    :param w: Third vector of the basis
    :type w: :class:`volmdlr.Vector3D`
    """
    _standalone_in_db = False

    # TODO: create a Basis and Frame class to mutualize between 2D and 2D
    def __init__(self, u: Vector3D, v: Vector3D, w: Vector3D, name: Text = ""):
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
        return "{}: U={}, V={}, W={}".format(self.__class__.__name__, *self.vectors)

    def _get_vectors(self):
        return (self.u, self.v, self.w)

    def to_dict(self, *args, **kwargs):
        """
        Serializes a 3-dimensional basis into a dictionary.

        :return: A serialized version of the Basis3D
        :rtype: dict

        .. seealso::
            How `serialization and deserialization`_ works in dessia_common

        .. _serialization and deserialization:
            https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method
        """
        return {"object_class": "volmdlr.Basis3D",
                "name": self.name,
                "u": self.u.to_dict(),
                "v": self.v.to_dict(),
                "w": self.w.to_dict()
                }

    vectors = property(_get_vectors)

    # TODO: transform to annotation when available
    @classmethod
    def from_two_vectors(cls, vector1: Vector3D, vector2: Vector3D, name: str = "") -> "Basis3D":
        """
        Creates a basis with first vector1 adimensionned, as u, v is the
        vector2 subtracted of u component, w is the cross product of u and v.

        :param vector1: The first vector of the Basis3D
        :type vector1: :class:`volmdlr.Vector3D`.
        :param name: object's name.
        :param vector2: The second vector of the Basis3D
        :type vector2: :class:`volmdlr.Vector3D`
        """
        u = vector1.copy()
        u = u.unit_vector()
        v = vector2 - vector2.dot(vector1) * vector1
        v = v.unit_vector()
        w = u.cross(v)
        return Basis3D(u, v, w, name=name)

    def to_frame(self, origin):
        """
        Returns the 3-dimensional frame oriented the same way as the Basis3D
        and having for origin the given 3-dimensional point.

        :param origin: The origin of the 3-dimensional frame
        :type origin: :class:`volmdlr.Point3D`
        :return: A 3-dimensional frame
        :rtype: :class:`volmdlr.Frame3D`
        """
        return Frame3D(origin, self.u, self.v, self.w)

    def rotation(self, axis: Vector3D, angle: float):
        """
        Rotates the 3-dimensional basis and returns a new rotated one.

        :param axis: The axis around which the rotation is made
        :type axis: :class:`volmdlr.Vector3D`
        :param angle: The angle of rotation in rad
        :type angle: float
        :return: The rotated Basis3D
        :rtype: :class:`volmdlr.Basis3D`
        """
        center = O3D
        new_u = self.u.rotation(center, axis, angle)
        new_v = self.v.rotation(center, axis, angle)
        new_w = self.w.rotation(center, axis, angle)
        return Basis3D(new_u, new_v, new_w, self.name)

    def x_rotation(self, angle: float):
        """
        Rotates the basis around the X axis and a new basis is returned
        as a result.

        :param angle: The rotation angle
        :type angle: float
        :return: The rotated Basis3D
        :rtype: :class:`volmdlr.Basis3D`
        """
        new_u = self.u.x_rotation(angle)
        new_v = self.v.x_rotation(angle)
        new_w = self.w.x_rotation(angle)
        return Basis3D(new_u, new_v, new_w, self.name)

    def y_rotation(self, angle: float):
        """
        Rotates the basis around the Y axis and a new basis is returned
        as a result.

        :param angle: The rotation angle
        :type angle: float
        :return: The rotated Basis3D
        :rtype: :class:`volmdlr.Basis3D`
        """
        new_u = self.u.y_rotation(angle)
        new_v = self.v.y_rotation(angle)
        new_w = self.w.y_rotation(angle)
        return Basis3D(new_u, new_v, new_w, self.name)

    def z_rotation(self, angle: float):
        """
        Rotates the basis around the Z axis and a new basis is returned
        as a result.

        :param angle: The rotation angle
        :type angle: float
        :return: The rotated Basis3D
        :rtype: :class:`volmdlr.Basis3D`
        """
        new_u = self.u.z_rotation(angle)
        new_v = self.v.z_rotation(angle)
        new_w = self.w.z_rotation(angle)
        return Basis3D(new_u, new_v, new_w, self.name)

    def euler_rotation_parameters(self, angles: Tuple[float, float, float]):
        """
        Computes the new basis' parameter after rotation of the basis using
        euler rotation.

        :param angles: Three angles corresponding to psi, theta, phi in rad
        :type angles: tuple
        :return: The new vectors of rotated basis
        :rtype: tuple
        """

        psi, theta, phi = angles
        center = O3D

        # rotation around w
        vect_u = self.u.rotation(center, self.w, psi)
        vect_v = self.v.rotation(center, self.w, psi)

        # rotation around v
        vect_v = vect_v.rotation(center, vect_u, theta)
        vect_w = self.w.rotation(center, vect_u, theta)

        # rotation around w
        vect_u = vect_u.rotation(center, vect_w, phi)
        vect_v = vect_v.rotation(center, vect_w, phi)
        return vect_u, vect_v, vect_w

    def euler_rotation(self, angles: Tuple[float, float, float]):
        """
        Rotates the 3-dimensional basis using euler rotation and
        returns a new basis as a result.

        :param angles: Three angles corresponding to psi, theta, phi in rad
        :type angles: tuple
        :return: The rotated basis
        :rtype: :class:`volmdlr.Basis3D`
        """
        vect_u, vect_v, vect_w = self.euler_rotation_parameters(angles)
        return Basis3D(vect_u, vect_v, vect_w)

    def transfer_matrix(self):
        """
        Computes the transfer matrix of the 3-dimensional basis.

        :return: The 3x3 transfer matrix
        :rtype: :class:`volmdlr.Matrix33`
        """
        return Matrix33(self.u.x, self.v.x, self.w.x,
                        self.u.y, self.v.y, self.w.y,
                        self.u.z, self.v.z, self.w.z)

    def inverse_transfer_matrix(self):
        """
        Computes the inverse transfer matrix of the 3-dimensional basis.

        :return: The 3x3 inverse transfer matrix
        :rtype: :class:`volmdlr.Matrix33`
        """
        return self.transfer_matrix().inverse()

    def new_coordinates(self, vector: Vector3D) -> Vector3D:
        """
        Convert the given vector's coordinates from the global landmark to the local landmark of this Basis3D.

        :param vector: The vector to convert, given in global coordinates.
        :type vector: :class:`volmdlr.Vector3D`
        :return: The converted vector, in local coordinates.
        :rtype: :class:`volmdlr.Vector3D`

        .. deprecated:: Use global_to_local_coordinates instead.
        """
        warnings.warn(
            "new_coordinates is deprecated. Use global_to_local_coordinates instead.",
            DeprecationWarning,
        )
        return self.global_to_local_coordinates(vector)

    def old_coordinates(self, vector: Vector3D) -> Vector3D:
        """
        Convert the given vector's coordinates from the local landmark of this Basis3D to the global landmark.

        :param vector: The vector to convert, given in local coordinates.
        :type vector: :class:`volmdlr.Vector3D`
        :return: The converted vector, in global coordinates.
        :rtype: :class:`volmdlr.Vector3D`

        .. deprecated:: Use local_to_global_coordinates instead.
        """
        warnings.warn(
            "old_coordinates is deprecated. Use local_to_global_coordinates instead.",
            DeprecationWarning,
        )
        return self.local_to_global_coordinates(vector)

    def global_to_local_coordinates(self, vector: Vector3D) -> Vector3D:
        """
        Convert the given vector's coordinates from the global landmark to the local landmark of this Basis3D.

        :param vector: The vector to convert, given in global coordinates.
        :type vector: :class:`volmdlr.Vector3D`
        :return: The converted vector, in local coordinates.
        :rtype: :class:`volmdlr.Vector3D`
        """
        matrix = self.inverse_transfer_matrix()
        return matrix.vector_multiplication(vector)

    def local_to_global_coordinates(self, vector: Vector3D) -> Vector3D:
        """
        Convert the given vector's coordinates from the local landmark of this Basis3D to the global landmark.

        :param vector: The vector to convert, given in local coordinates.
        :type vector: :class:`volmdlr.Vector3D`
        :return: The converted vector, in global coordinates.
        :rtype: :class:`volmdlr.Vector3D`
        """
        matrix = self.transfer_matrix()
        return matrix.vector_multiplication(vector)

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of a 3-dimensional basis.

        :param deep: *not used*
        :param memo: *not used*
        :return: A copy of the Basis3D
        :rtype: :class:`volmdlr.Basis3D`
        """
        return Basis3D(self.u, self.v, self.w)

    def normalize(self):
        """
        Normalizes the basis, modifying its coordinates in place.

        :return: New normalized basis
        :rtype: Basis3D
        """
        u, v, w = self.u, self.v, self.w
        if not math.isclose(self.u.norm(), 0.0, abs_tol=1e-10):
            u = self.u.unit_vector()
        if not math.isclose(self.v.norm(), 0.0, abs_tol=1e-10):
            v = self.v.unit_vector()
        if not math.isclose(self.w.norm(), 0.0, abs_tol=1e-10):
            w = self.w.unit_vector()
        return Basis3D(u, v, w)

    def is_orthogonal(self, tol: float = 1e-6):
        """
        Check if the basis vectors are orthogonal to each other.

        :param tol: Tolerance for considering a dot product as zero.
        :type tol: float
        :return: True if the basis vectors are orthogonal, False otherwise.
        :rtype: bool
        """
        dot_uv = self.u.dot(self.v)
        dot_uw = self.u.dot(self.w)
        dot_vw = self.v.dot(self.w)

        return abs(dot_uv) < tol and abs(dot_uw) < tol and abs(dot_vw) < tol

    def is_normalized(self, tol: float = 1e-6):
        """
        Check if the basis vectors are normal to each other.

        :param tol: Tolerance for considering a unit vector.
        :type tol: float
        :return: True if the basis vectors are normalized, False otherwise.
        :rtype: bool
        """
        return (math.isclose(self.u.norm(), 1.0, abs_tol=tol) and math.isclose(self.v.norm(), 1.0, abs_tol=tol) and
                math.isclose(self.w.norm(), 1.0, abs_tol=tol))

    def is_orthonormal(self, tol: float = 1e-6):
        """
        Check if the basis vectors are orthonormal to each other.

        :param tol: Tolerance for considering an orthonormal basis.
        :type tol: float
        :return: True if the basis vectors are orthonormal, False otherwise.
        :rtype: bool
        """
        return self.is_orthogonal(tol) and self.is_normalized(tol)


class Frame2D(Basis2D):
    """
    Defines a 2D basis
    :param origin:Point2D: origin of the basis
    :param u:Vector2D: first vector of the basis
    :param v:Vector2D: second vector of the basis
    """

    def __init__(self, origin: Point2D, u: Vector2D, v: Vector2D, name: Text = ""):
        self.origin = origin
        Basis2D.__init__(self, u, v, name=name)

    def __repr__(self):
        return "{}: O={} U={}, V={}".format(self.__class__.__name__, self.origin, self.u, self.v)

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

    def __hash__(self):
        """
        Hash returns 0 because points are difficult to hash if they are meant
        to be equalized at a given tolerance
        """
        return 0

    def to_dict(self, *args, **kwargs):
        """
        Serializes a 2-dimensional frame into a dictionary.

        :return: A serialized version of the Frame2D
        :rtype: dict

        .. seealso::
            How `serialization and deserialization`_ works in dessia_common

        .. _serialization and deserialization:
            https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method
        """
        return {"object_class": "volmdlr.Frame2D",
                "name": self.name,
                "origin": self.origin.to_dict(),
                "u": self.u.to_dict(),
                "v": self.v.to_dict()
                }

    def normalize(self):
        """
        Normalizes the Frame, and return a new object.

        :return: A new normalized basis
        :rtype: Frame2D
        """
        return Frame2D(self.origin, self.u.unit_vector(), self.v.unit_vector())

    def basis(self):
        """
        Returns the 2-dimensional basis oriented the same way as the Frame2D.

        :return: A 2-dimensional basis
        :rtype: :class:`volmdlr.Basis2D`
        """
        return Basis2D(self.u, self.v)

    def new_coordinates(self, vector: Vector2D) -> Vector2D:
        """
        Convert the given vector's coordinates from the global landmark to the local landmark of this Frame2D.

        :param vector: The vector to convert, given in global coordinates.
        :type vector: :class:`volmdlr.Vector2D`
        :return: The converted vector, in local coordinates.
        :rtype: :class:`volmdlr.Vector2D`

        .. deprecated:: Use global_to_local_coordinates instead.
        """
        warnings.warn(
            "new_coordinates is deprecated. Use global_to_local_coordinates instead.",
            DeprecationWarning,
        )
        return self.global_to_local_coordinates(vector)

    def old_coordinates(self, vector: Vector2D) -> Vector2D:
        """
        Convert the given vector's coordinates from the local landmark of this Frame2D to the global landmark.

        :param vector: The vector to convert, given in local coordinates.
        :type vector: :class:`volmdlr.Vector2D`
        :return: The converted vector, in global coordinates.
        :rtype: :class:`volmdlr.Vector2D`

        .. deprecated:: Use local_to_global_coordinates instead.
        """
        warnings.warn(
            "old_coordinates is deprecated. Use local_to_global_coordinates instead.",
            DeprecationWarning,
        )
        return self.local_to_global_coordinates(vector)

    def global_to_local_coordinates(self, vector: Vector2D) -> Vector2D:
        """
        Convert the given vector's coordinates from the global landmark to the local landmark of this Frame2D.

        :param vector: The vector to convert, given in global coordinates.
        :type vector: :class:`volmdlr.Vector2D`
        :return: The converted vector, in local coordinates.
        :rtype: :class:`volmdlr.Vector2D`
        """
        return Basis2D.global_to_local_coordinates(self, vector - self.origin)

    def local_to_global_coordinates(self, vector: Vector2D) -> Vector2D:
        """
        Convert the given vector's coordinates from the local landmark of this Frame2D to the global landmark.

        :param vector: The vector to convert, given in local coordinates.
        :type vector: :class:`volmdlr.Vector2D`
        :return: The converted vector, in global coordinates.
        :rtype: :class:`volmdlr.Vector2D`
        """
        return Basis2D.local_to_global_coordinates(self, vector) + self.origin

    def frame_mapping(self, frame: "Frame2D", side: str):
        basis = frame.basis()
        if side == "new":
            new_origin = frame.global_to_local_coordinates(self.origin)
            new_u = basis.global_to_local_coordinates(self.u)
            new_v = basis.global_to_local_coordinates(self.v)
        elif side == "old":
            new_origin = frame.local_to_global_coordinates(self.origin)
            new_u = basis.local_to_global_coordinates(self.u)
            new_v = basis.local_to_global_coordinates(self.v)
        else:
            raise ValueError("side value not valid, please specify a correct value: \'old\' or \'new\'")
        return Frame2D(new_origin, new_u, new_v)

    def translation(self, vector):
        """
        Returns a translated 2-dimensional frame.

        :param vector: The translation vector
        :type vector: :class:`volmdlr.Vector2D`
        :return: A new translated 2-dimensional frame
        :rtype: :class:`volmdlr.Frame2D`
        """
        new_origin = self.origin.translation(vector)
        return Frame2D(new_origin, self.u, self.v)

    def rotation(self, angle):
        """
        Returns a rotated 2-dimensional frame.

        :param angle: The rotation angle
        :type angle: float
        :return: New rotated frame
        :rtype: :class:`volmdlr.Frame2D`
        """
        new_base = Basis2D.rotation(self, angle)
        return Frame2D(self.origin, new_base.u, new_base.v)

    def Draw(self, ax=None, style="ok"):
        """
        # TODO : unused ? to be deleted ?

        :param ax:
        :param style:
        :return:
        """
        if ax is None:
            fig, ax = plt.subplots()

        ax.plot(*self.origin, style)
        self.u.plot(origin=self.origin, ax=ax, color="r")
        self.v.plot(origin=self.origin, ax=ax, color="g")
        ax.axis("equal")

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of a 2-dimensional frame.

        :param deep: *not used*
        :param memo: *not used*
        :return: A copy of the Frame2D
        :rtype: :class:`volmdlr.Frame2D`
        """
        return Frame2D(self.origin, self.u, self.v)

    def plot(self, ax=None, color="b", alpha=1., plot_points=True,
             ratio=1.):
        """
        Plots the 2-dimensional frame.

        :param ax: The Axes on which the Point2D will be drawn. Default value
            is None, creating a new drawing figure
        :type ax: :class:`matplotlib.axes.Axes`, optional
        :param color: *not used*
        :type color: str, optional
        :param alpha: *not used*
        :type alpha: float, optional
        :param plot_points: *not used*
        :type plot_points: bool, optional
        :param ratio: A ratio for controlling the size of the 2-dimensional
            frame. Default value is 1
        :type ratio: float, optional
        :return: A matplotlib Axes object on which the Point2D have been
            plotted
        :rtype: :class:`matplotlib.axes.Axes`
        """
        if ax is None:
            fig, ax = plt.subplots()

        x1 = [p.x for p in (self.origin, self.origin + self.u * ratio)]
        y1 = [p.y for p in (self.origin, self.origin + self.u * ratio)]
        ax.plot(x1, y1, "r")

        x2 = [p.x for p in (self.origin, self.origin + self.v * ratio)]
        y2 = [p.y for p in (self.origin, self.origin + self.v * ratio)]
        ax.plot(x2, y2, "g")
        ax.axline((self.origin.x, self.origin.y), (self.u.x, self.u.y),
                  color="black", linewidth=0.8, linestyle="--")
        ax.axline((self.origin.x, self.origin.y), (self.v.x, self.v.y),
                  color="black", linewidth=0.8, linestyle="--")
        return ax


OXY = Frame2D(O2D, X2D, Y2D)


class Frame3D(Basis3D):
    """
    Defines a 3D frame
    :param origin:Point3D: origin of the basis
    :param u:Vector3D: first vector of the basis
    :param v:Vector3D: second vector of the basis
    :param w:Vector3D: third vector of the basis
    """

    def __init__(self, origin: Point3D, u: Vector3D, v: Vector3D, w: Vector3D, name: Text = ""):
        self.origin = origin
        Basis3D.__init__(self, u, v, w)
        self.name = name

    def __repr__(self):
        return f"{self.__class__.__name__}(origin={self.origin}, u={self.u}, v={self.v}, w={self.w})"

    def __hash__(self):
        """
        hash returns 0 because points are difficult to hash if they are meant
        to be equalized at a given tolerance
        """
        return hash((self.origin, self.u, self.v, self.w))

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

    def to_dict(self, *args, **kwargs):
        """
        Serializes a 3-dimensional frame into a dictionary.

        :return: A serialized version of the Frame3D
        :rtype: dict

        .. seealso::
            How `serialization and deserialization`_ works in dessia_common

        .. _serialization and deserialization:
            https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method
        """
        return {"object_class": "volmdlr.Frame3D",
                "name": self.name,
                "origin": self.origin.to_dict(),
                "u": self.u.to_dict(),
                "v": self.v.to_dict(),
                "w": self.w.to_dict()
                }

    def normalize(self):
        """
        Normalizes the Frame, and return a new object.

        :return: A new normalized basis
        :rtype: Frame2D
        """
        return Frame3D(self.origin, self.u.unit_vector(), self.v.unit_vector(), self.w.unit_vector())

    def basis(self):
        """
        Returns the 3-dimensional basis oriented the same way as the Frame3D.

        :return: A 3-dimensional basis
        :rtype: :class:`volmdlr.Basis3D`
        """
        return Basis3D(self.u, self.v, self.w)

    def new_coordinates(self, vector: Vector3D) -> Vector3D:
        """
        Convert the given vector's coordinates from the global landmark to the local landmark of this Frame3D.

        :param vector: The vector to convert, given in global coordinates.
        :type vector: :class:`volmdlr.Vector3D`
        :return: The converted vector, in local coordinates.
        :rtype: :class:`volmdlr.Vector3D`

        .. deprecated:: Use global_to_local_coordinates instead.
        """
        warnings.warn(
            "new_coordinates is deprecated. Use global_to_local_coordinates instead.",
            DeprecationWarning,
        )
        return self.global_to_local_coordinates(vector)

    def old_coordinates(self, vector: Vector3D) -> Vector3D:
        """
        Convert the given vector's coordinates from the local landmark of this Frame3D to the global landmark.

        :param vector: The vector to convert, given in local coordinates.
        :type vector: :class:`volmdlr.Vector3D`
        :return: The converted vector, in global coordinates.
        :rtype: :class:`volmdlr.Vector3D`

        .. deprecated:: Use local_to_global_coordinates instead.
        """
        warnings.warn(
            "old_coordinates is deprecated. Use local_to_global_coordinates instead.",
            DeprecationWarning,
        )
        return self.local_to_global_coordinates(vector)

    def global_to_local_coordinates(self, vector: Vector3D) -> Vector3D:
        """
        Convert the given vector's coordinates from the global landmark to the local landmark of this Frame3D.

        :param vector: The vector to convert, given in global coordinates.
        :type vector: :class:`volmdlr.Vector3D`
        :return: The converted vector, in local coordinates.
        :rtype: :class:`volmdlr.Vector3D`
        """
        return Basis3D.global_to_local_coordinates(self, vector - self.origin)

    def local_to_global_coordinates(self, vector: Vector3D) -> Vector3D:
        """
        Convert the given vector's coordinates from the local landmark of this Frame3D to the global landmark.

        :param vector: The vector to convert, given in local coordinates.
        :type vector: :class:`volmdlr.Vector3D`
        :return: The converted vector, in global coordinates.
        :rtype: :class:`volmdlr.Vector3D`
        """
        return Basis3D.local_to_global_coordinates(self, vector) + self.origin

    def frame_mapping(self, frame: "Frame3D", side: str):
        basis = frame.basis()
        if side == "new":
            new_origin = frame.global_to_local_coordinates(self.origin)
            new_u = basis.global_to_local_coordinates(self.u)
            new_v = basis.global_to_local_coordinates(self.v)
            new_w = basis.global_to_local_coordinates(self.w)

        elif side == "old":
            new_origin = frame.local_to_global_coordinates(self.origin)
            new_u = basis.local_to_global_coordinates(self.u)
            new_v = basis.local_to_global_coordinates(self.v)
            new_w = basis.local_to_global_coordinates(self.w)
        else:
            raise ValueError("side value not valid, please specify"
                             'a correct value: \'old\' or \'new\'')
        return Frame3D(new_origin, new_u, new_v, new_w)

    def rotation(self, center: Point3D, axis: Vector3D, angle: float):
        """
        Rotates the center as a point and vectors as directions
        (calling Basis), and returns a new 3-dimensional frame.

        :param center: The center of rotation
        :type center: :class:`volmdlr.Point3D`
        :param axis: The axis around which the rotation will be made
        :type axis: :class:`volmdlr.Vector3D`
        :param angle: The rotation angle
        :type angle: float
        :return: New rotated frame
        :rtype: :class:`volmdlr.Frame3D`
        """
        new_base = Basis3D.rotation(self, axis, angle)
        new_origin = self.origin.rotation(center, axis, angle)
        return Frame3D(new_origin,
                       new_base.u, new_base.v, new_base.w,
                       self.name)

    def translation(self, offset: Vector3D):
        """
        Translates a 3-dimensional frame.

        :param offset: The translation vector
        :type offset: :class:`volmdlr.Vector3D`
        :return: new translated frame
        :rtype: Frame3D
        """
        return Frame3D(self.origin.translation(offset),
                       self.u, self.v, self.w, self.name)

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of a 3-dimensional frame.

        :param deep: *not used*
        :param memo: *not used*
        :return: A copy of the Frame3D
        :rtype: :class:`volmdlr.Frame3D`
        """
        return Frame3D(self.origin.copy(),
                       self.u.copy(), self.v.copy(), self.w.copy())

    def to_step(self, current_id):
        """
        Writes a step primitive from a 3-dimensional frame.

        :param current_id: The id of the last written primitive
        :type current_id: int
        :return: A tuple containing the string representing the step primitive
            and the new current id
        :rtype: tuple
        """
        content, origin_id = self.origin.to_point().to_step(current_id)
        current_id = origin_id + 1
        w_content, w_id = Vector3D.to_step(self.w, current_id)
        current_id = w_id + 1
        u_content, u_id = Vector3D.to_step(self.u, current_id)
        current_id = u_id + 1
        content += w_content + u_content
        content += f"#{current_id} = AXIS2_PLACEMENT_3D('{self.name}',#{origin_id},#{w_id},#{u_id});\n"
        return content, current_id

    def plot2d(self, x=X3D, y=Y3D, ax=None, color="k"):
        """
        Plots the 3-dimensional frame on a 2-dimensional surface given
        by (x, y).

        :param x: The first 3-dimensional vector of the 2-dimensional surface.
            Default value is X3D, the vector (1, 0, 0)
        :type x: :class:`volmdlr.Vector3D`, optional
        :param y: The second 3-dimensional vector of the 2-dimensional surface.
            Default value is Y3D, the vector (0, 1, 0)
        :type y: :class:`volmdlr.Vector3D`, optional
        :param ax: The Axes on which the Frame3D will be drawn. Default value
            is None, creating a new drawing figure
        :type ax: :class:`matplotlib.axes.Axes`, optional
        :param color: The color of the frame. Default value is 'k', for black
        :type color: str, optional
        :return: A matplotlib Axes object on which the 2-dimensional
            representation of the Frame3D have been plotted
        :rtype: :class:`matplotlib.axes.Axes`
        """
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

    def plot(self, ax=None, color="b", alpha=1., plot_points=True,
             ratio=1.):
        """
        Plots the 3-dimensional frame.

        :param ax: The Axes on which the Point3D will be drawn. Default value
            is None, creating a new drawing figure
        :type ax: :class:`matplotlib.axes.Axes`, optional
        :param color: *not used*
        :type color: str, optional
        :param alpha: *not used*
        :type alpha: float, optional
        :param plot_points: *not used*
        :type plot_points: bool, optional
        :param ratio: A ratio for controlling the size of the 3-dimensional
            frame. Default value is 1
        :type ratio: float, optional
        :return: A matplotlib Axes object on which the Point3D have been
            plotted
        :rtype: :class:`matplotlib.axes.Axes`
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection="3d")

        x1 = [p.x for p in (self.origin, self.origin + self.u * ratio)]
        y1 = [p.y for p in (self.origin, self.origin + self.u * ratio)]
        z1 = [p.z for p in (self.origin, self.origin + self.u * ratio)]
        ax.plot(x1, y1, z1, "r")

        x2 = [p.x for p in (self.origin, self.origin + self.v * ratio)]
        y2 = [p.y for p in (self.origin, self.origin + self.v * ratio)]
        z2 = [p.z for p in (self.origin, self.origin + self.v * ratio)]
        ax.plot(x2, y2, z2, "g")

        x3 = [p.x for p in (self.origin, self.origin + self.w * ratio)]
        y3 = [p.y for p in (self.origin, self.origin + self.w * ratio)]
        z3 = [p.z for p in (self.origin, self.origin + self.w * ratio)]
        ax.plot(x3, y3, z3, "b")
        return ax

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive from a 3-dimensional point to a Frame3D.

        :param arguments: The arguments of the step primitive. The last element represents the unit_conversion_factor.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated
        :type object_dict: dict
        :return: The corresponding Frame3D object
        :rtype: :class:`volmdlr.Frame3D`
        """
        origin = object_dict[arguments[1]]
        if arguments[2] == "$" and arguments[3] == "$":
            return cls(origin, volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D, arguments[0][1:-1])
        if arguments[2] == "$":
            return cls.from_point_and_vector(origin, object_dict[arguments[3]], main_axis=X3D,
                                             name= arguments[0][1:-1])
        if arguments[3] == "$":
            return cls.from_point_and_vector(origin, object_dict[arguments[2]], main_axis=Z3D,
                                             name= arguments[0][1:-1])
        w = object_dict[arguments[2]]
        u = object_dict[arguments[3]]
        u = u - u.dot(w) * w
        u = u.unit_vector()
        v = w.cross(u)
        return cls(origin, u, v, w, arguments[0][1:-1])

    @classmethod
    def from_point_and_vector(cls, point: Point3D, vector: Vector3D,
                              main_axis: Vector3D = X3D, name: str = ""):
        """
        Creates a new frame from a point and vector by rotating the global
        frame. Global frame rotates in order to have 'vector' and 'main_axis'
        collinear. This method is very useful to compute a local frame of
        an object.

        :param point: The origin of the new frame
        :type point: :class:`volmdlr.Point3D`
        :param vector: The vector used to define one of the main axis
            (by default X-axis) of the local frame
        :type vector: :class:`volmdlr.Vector3D`
        :param main_axis: The axis of global frame you want to match 'vector'
            (can be X3D, Y3D or Z3D). Default value is X3D,
            the vector (1, 0, 0)
        :type main_axis: :class:`volmdlr.Vector3D`, optional
        :param name: Frame's name.
        :type name: str
        :return: The created local frame
        :rtype: :class:`volmdlr.Frame3D`
        """
        if main_axis not in [X3D, Y3D, Z3D]:
            raise ValueError("main_axis must be X, Y or Z of the global frame")

        vector = vector.unit_vector()

        if vector == main_axis:
            # The local frame is oriented like the global frame
            return cls(point, X3D, Y3D, Z3D)

        if vector == -main_axis:
            return cls(point, -X3D, -Y3D, -Z3D)

        # The local frame is oriented differently from the global frame
        # Rotation angle
        dot = main_axis.dot(vector)
        rot_angle = math.acos(dot / (vector.norm() * main_axis.norm()))

        # Rotation axis
        vector2 = vector - main_axis
        rot_axis = main_axis.cross(vector2)
        rot_axis = rot_axis.unit_vector()

        u = X3D.rotation(O3D, rot_axis, rot_angle)
        v = Y3D.rotation(O3D, rot_axis, rot_angle)
        w = Z3D.rotation(O3D, rot_axis, rot_angle)

        return cls(point, u, v, w, name=name)

    @classmethod
    def from_3_points(cls, point1, point2, point3, name: str = ""):
        """
        Creates a frame 3d from 3 points.

        :param point1: point 1.
        :param point2: point 2.
        :param point3: point 3.
        :param name: object's name.
        :return: a frame 3d.
        """
        vector1 = point2 - point1
        vector2 = point3 - point1
        vector1 = vector1.to_vector().unit_vector()
        vector2 = vector2.to_vector().unit_vector()
        normal = vector1.cross(vector2)
        normal = normal.unit_vector()
        return cls(point1, vector1, normal.cross(vector1), normal, name=name)

    @classmethod
    def from_point_and_normal(cls, origin, normal, name: str = ""):
        """Creates a frame 3D from a point and a normal vector."""
        u_vector = normal.deterministic_unit_normal_vector()
        v_vector = normal.cross(u_vector)
        return cls(origin, u_vector, v_vector, normal, name=name)

    # def babylonjs(self, size=0.1, parent=None):
    #     """
    #     # TODO: to be deleted ?
    #     Returns the babylonjs script for 3D display in browser.

    #     :param size: The adjustable size of the 3-dimensional frame. Default
    #         value is 0.1
    #     :type size: float, optional
    #     :param parent:
    #     :type parent:
    #     :return: A babylonjs script
    #     :rtype: str
    #     """
    #     s = "var origin = new BABYLON.Vector3({},{},{});\n".format(*self.origin)
    #     s += "var o_u = new BABYLON.Vector3({}, {}, {});\n".format(*(size * self.u + self.origin))
    #     s += "var o_v = new BABYLON.Vector3({}, {}, {});\n".format(*(size * self.v + self.origin))
    #     s += "var o_w = new BABYLON.Vector3({}, {}, {});\n".format(*(size * self.w + self.origin))
    #     s += 'var line1 = BABYLON.MeshBuilder.CreateTube("frame_U",{{path:[origin, o_u], radius:{}}},scene);'.format(
    #         0.03 * size)
    #     s += "line1.material = red_material;\n"
    #     s += 'var line2 = BABYLON.MeshBuilder.CreateTube("frame_V",{{path:[origin, o_v], radius:{}}},scene);'.format(
    #         0.03 * size)
    #     s += "line2.material = green_material;\n"
    #     s += 'var line3 = BABYLON.MeshBuilder.CreateTube("frame_W",{{path:[origin, o_w], radius:{}}},scene);'.format(
    #         0.03 * size)
    #     s += "line3.material = blue_material;\n"
    #     if parent is not None:
    #         s += "line1.parent = {};\n".format(parent)
    #         s += "line2.parent = {};\n".format(parent)
    #         s += "line3.parent = {};\n".format(parent)

    #     return s
