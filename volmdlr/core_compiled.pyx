#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#cython: language_level=3
"""

Cython functions

"""
# from __future__ import annotations
from typing import TypeVar, List, Tuple
import math
from dessia_common import DessiaObject
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow, FancyArrowPatch
import warnings
import random
import numpy as npy
from mpl_toolkits.mplot3d import proj3d

# =============================================================================

cdef (double, double) Csub2D(double u1, double u2,
                             double v1, double v2):
    return (u1-v1, u2-v2)

def sub2D(vector1, vector2):
    return Csub2D(vector1[0], vector1[1],
                  vector2[0], vector2[1])

# =============================================================================

cdef (double, double) Cadd2D(double u1, double u2,
                             double v1, double v2,):
    return (u1+v1, u2+v2)

def add2D(vector1, vector2):
    return Cadd2D(vector1[0], vector1[1],
                  vector2[0], vector2[1])

# =============================================================================

cdef (double, double) Cmul2D(double u1, double u2, double value):
    return (u1*value, u2*value)

def mul2D(vector, value):
    return Cmul2D(vector[0], vector[1], value)

# =============================================================================

cdef double CVector2DDot(double u1, double u2,
                         double v1, double v2):
    return u1*v1 + u2*v2

def Vector2DDot(vector1, vector2):
    return CVector2DDot(vector1[0], vector1[1],
                        vector2[0], vector2[1])

# =============================================================================

cdef double CVector2DNorm(double u1, double u2):
    return (u1*u1 + u2*u2)**0.5

def Vector2DNorm(vector):
    return CVector2DNorm(vector[0], vector[1])

# =============================================================================

cdef (double, double, double) Csub3D(double u1, double u2, double u3,
                                     double v1, double v2, double v3):
    return (u1-v1, u2-v2, u3-v3)

def sub3D(vector1, vector2):
    return Csub3D(vector1[0], vector1[1], vector1[2],
                  vector2[0], vector2[1], vector2[2])

# =============================================================================

cdef (double, double, double) Cadd3D(double u1, double u2, double u3,
                                     double v1, double v2, double v3):
    return (u1+v1, u2+v2, u3+v3)

def add3D(vector1, vector2):
    return Cadd3D(vector1[0], vector1[1], vector1[2],
                  vector2[0], vector2[1], vector2[2])

# =============================================================================

cdef (double, double, double) Cmul3D(double u1, double u2, double u3,
                                     double value):
    return (u1*value, u2*value, u3*value)

def mul3D(vector, value):
    return Cmul3D(vector[0], vector[1], vector[2], value)

# =============================================================================

cdef double CVector3DDot(double u1, double u2, double u3,
                         double v1, double v2, double v3):
    return u1*v1 + u2*v2 + u3*v3

def Vector3DDot(vector1, vector2):
    return CVector3DDot(vector1[0], vector1[1], vector1[2],
                        vector2[0], vector2[1], vector2[2])

# =============================================================================

cdef double CVector3DNorm(double u1, double u2, double u3):
    return (u1*u1 + u2*u2 + u3*u3)**0.5

def Vector3DNorm(vector):
    return CVector3DNorm(vector[0], vector[1], vector[2])


# =============================================================================

cdef (double, double, double) C_vector3D_cross(double u1, double u2, double u3,
                                               double v1, double v2, double v3):
    return (u2*v3 - u3*v2, u3*v1 - u1*v3, u1*v2 - u2*v1)

def vector3D_cross(vector1, vector2):
    return C_vector3D_cross(vector1[0], vector1[1], vector1[2],
                            vector2[0], vector2[1], vector2[2])


# =============================================================================

cdef (double, double, double) C_vector3D_rotation(double vx, double vy, double vz,
                                                  double center_x, double center_y, double center_z,
                                                  double axis_x, double axis_y, double axis_z,
                                                  double angle):
    
    cdef double ux = vx - center_x
    cdef double uy = vy - center_y
    cdef double uz = vz - center_z
    
    cdef double cos_angle = math.cos(angle)
    cdef double sin_angle = math.sin(angle)
    
    cdef double rv1_x = cos_angle*ux
    cdef double rv1_y = cos_angle*uy
    cdef double rv1_z = cos_angle*uz
    
    rv2_x, rv2_y, rv2_z = Cmul3D(axis_x, axis_y, axis_z,
                                 (1-cos_angle)*CVector3DDot(
                                         ux, uy, uz,
                                         axis_x, axis_y, axis_z)
                                 )
    
    rv3_x, rv3_y, rv3_z = C_vector3D_cross(axis_x, axis_y, axis_z,
                                           ux, uy, uz)
    
    return (rv1_x + rv2_x + rv3_x*sin_angle + center_x,
            rv1_y + rv2_y + rv3_y*sin_angle + center_y,
            rv1_z + rv2_z + rv3_z*sin_angle + center_z)

def vector3D_rotation(vector, center, axis, angle):
        return C_vector3D_rotation(vector[0], vector[1], vector[2],
                                   center[0], center[1], center[2],
                                   axis[0], axis[1], axis[2],
                                   angle)
        
    
    
cdef (double, double, double) C_matrix_vector_multiplication3(double M11, double M12, double M13,
                                                              double M21, double M22, double M23,
                                                              double M31, double M32, double M33,
                                                              double v1, double v2, double v3):

    return (M11*v1 + M12*v2 + M13*v3,
            M21*v1 + M22*v2 + M23*v3,
            M31*v1 + M32*v2 + M33*v3)


cdef (double, double, double,
      double, double, double,
      double, double, double) C_matrix_multiplication3(double A11, double A12, double A13,
                                                       double A21, double A22, double A23,
                                                       double A31, double A32, double A33,
                                                       double B11, double B12, double B13,
                                                       double B21, double B22, double B23,
                                                       double B31, double B32, double B33):

    return (A11*B11 + A12*B21 + A13*B31,
            A11*B12 + A12*B22 + A13*B32,
            A11*B13 + A12*B23 + A13*B33,
            A21*B11 + A22*B21 + A23*B31,
            A21*B12 + A22*B22 + A23*B32,
            A21*B13 + A22*B23 + A23*B33,
            A31*B11 + A32*B21 + A33*B31,
            A31*B12 + A32*B22 + A33*B32,
            A31*B13 + A32*B23 + A33*B33)
    

# =============================================================================

def PolygonPointBelongs(point, points):

    cdef int i
    cdef int n = len(points)
    cdef bint inside = False
    cdef float x, y, p1x, p1y, p2x, p2y, xints
    x,y=point
    p1x,p1y = points[0]

    for i in range(n+1):
        p2x, p2y = points[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside

# =============================================================================

cdef (double, (double, double)) CLineSegment2DPointDistance((double, double) p1, (double, double) p2, (double, double) point):
    cdef double t
    cdef (double, double) u, projection

    u = sub2D(p2, p1)
    t = max(0, min(1, Vector2DDot(sub2D(point, p1), u) / Vector2DNorm(u)**2))
    projection = add2D(p1, mul2D(u, t))
    return Vector2DNorm(sub2D(projection, point)), projection

def LineSegment2DPointDistance(points, point):
    return CLineSegment2DPointDistance(tuple(points[0]), tuple(points[1]), tuple(point))


# =============================================================================
#  Points, Vectors
# =============================================================================

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

class Vector(DessiaObject):
    """
    Abstract class of vector
    """
    def __setitem__(self, key, item):
        self.vector[key] = item

    def __getitem__(self, key):
        return self.vector[key]

    def __repr__(self):
        return '{}: {}'.format(self.__class__.__name__, self.vector)

    def __radd__(self, other_vector):
        return self + other_vector

    def __rsub__(self, other_vector):
        return self - other_vector

    def __rmul__(self, value):
        return self * value

    def __rtruediv__(self, value):
        return self / value

    def __lt__(self, other_vector):
        return self.Norm() < other_vector.Norm()

    def __le__(self, other_vector):
        return self.Norm() <= other_vector.Norm()

    # def __ne__(self, other_vector):
    #     return not npy.allclose(self.vector, other_vector.vector)

    def to_numpy(self):
        return npy.array(self.vector)

    def copy(self):
        return self.__class__(self.vector)

#    def Dict(self):
#        d = {'vector': [float(i) for i in self.vector],
#             'object_class' : self.full_classname}
#        return d

    @classmethod
    def mean_point(cls, points):
        n = 1
        point = points[0].copy()
        for point2 in points[1:]:
            point += point2
            n += 1
        point /= n
        return point

class Vector2D(Vector):
    def __init__(self, vector:Tuple[float, float], name=''):
        # TODO: change this list to 2 values vx and vy
        self.vector = [0, 0]
#        self.vector = npy.zeros(2)
        self.vector[0] = vector[0]
        self.vector[1] = vector[1]
        self.name = name

    def __add__(self, other_vector):
        return Vector2D(add2D(self.vector, other_vector.vector))

    def __neg__(self):
        return Vector2D((-self.vector[0], -self.vector[1]))

    def __sub__(self, other_vector):
        return Vector2D(sub2D(self.vector, other_vector.vector))

    def __mul__(self, value:float):
        return Vector2D(mul2D(self.vector, value))

    def __truediv__(self, value:float):
        if value == 0:
            raise ZeroDivisionError
        return Vector2D((self.vector[0] / value,
                         self.vector[1] / value))

    def __round__(self, ndigits:int=6):
        return self.__class__((round(self.vector[0], ndigits),
                               round(self.vector[1], ndigits)))

    def __hash__(self):
#        return int(1000*(self.vector[0]+self.vector[1]))
        return int(round(1e6*(self.vector[0]+self.vector[1])))

    def __eq__(self, other_vector):

        if other_vector.__class__.__name__ not in ['Vector2D', 'Point2D']:
            return False
        return math.isclose(self.vector[0], other_vector.vector[0], abs_tol=1e-06) \
        and math.isclose(self.vector[1], other_vector.vector[1], abs_tol=1e-06)

    def Norm(self):
        """
        :returns: norm of vector
        """
        return Vector2DNorm(self.vector)

    def Normalize(self):
        """
        Normalize the vector modifying its coordinates
        """
        n = self.Norm()
        if math.isclose(n, 0, abs_tol=1e-9):
            raise ZeroDivisionError

        self.vector[0] /= n
        self.vector[1] /= n

    def Dot(self, other_vector):
        v1, v2 = self.vector
        ov1, ov2 = other_vector.vector
        return CVector2DDot(v1, v2, ov1, ov2)

    def Cross(self, other_vector):
        u1, u2 = self.vector
        v1, v2 = other_vector.vector
        return u1*v2 - u2*v1

    def point_distance(self, point2):
        return (self-point2).Norm()


    def Rotation(self, center, angle, copy=True):
        u = self - center
        vector2 = [math.cos(angle)*u[0] - math.sin(angle)*u[1] + center[0],
                   math.sin(angle)*u[0] + math.cos(angle)*u[1] + center[1]]
#        vector2 = (npy.dot(npy.array([[math.cos(angle), -math.sin(angle)],
#                                      [math.sin(angle), math.cos(angle)]]),
#                           (self.vector-center.vector))
#                   + center.vector)
        if copy:
            return self.__class__(vector2)
        else:
            self.vector = vector2

    def Translation(self, offset, copy=True):
        """
        :param offset: an other Vector2D
        """
        vector2 = [self.vector[0] + offset[0],
                   self.vector[1] + offset[1]]
        if copy:
            return self.__class__(vector2)
        else:
            self.vector = vector2

    def frame_mapping(self, frame, side, copy=True):
        # """
        # side = 'old' or 'new'
        # """
        if side == 'old':
            new_vector = frame.OldCoordinates(self)
            if copy:
                return new_vector
            else:
                self.vector = new_vector.vector

        if side == 'new':
            new_vector = frame.NewCoordinates(self)
            if copy:
                return new_vector
            else:
                self.vector = new_vector.vector

    def To3D(self, plane_origin, vx, vy):
        return Vector3D([plane_origin.vector[0] + vx.vector[0]*self.vector[0] + vy.vector[0]*self.vector[1],
                         plane_origin.vector[1] + vx.vector[1]*self.vector[0] + vy.vector[1]*self.vector[1],
                         plane_origin.vector[2] + vx.vector[2]*self.vector[0] + vy.vector[2]*self.vector[1],
                         ])

    def NormalVector(self, unit=False):
        n = Vector2D((-self.vector[1], self.vector[0]))
        if unit:
            n.Normalize()
        return n
    
    def deterministic_unit_normal_vector(self):
        return self.NormalVector(unit=True)

    @classmethod
    def random(cls, xmin, xmax, ymin, ymax):
        return cls((random.uniform(xmin, xmax),
                    random.uniform(ymin, ymax)))

    def Draw(self):
        warnings.warn(
            "Draw is deprecated and will be removed in next versions, use plot() instead",
            DeprecationWarning
        )
        self.plot()

    def plot(self, amplitude=0.5, width=None, head_width=None, origin=None, ax=None, color='k', line=False, label=None, normalize=False):
        if origin is None:
            origin = Vector2D((0., 0.))

        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

#        if self.vector == [0., 0.]:
        if math.isclose(self.Norm(), 0, abs_tol=1e-9):
            point = Point2D(origin.vector)
            point.MPLPlot(ax=ax, color=color)
            return fig, ax
        
        if width is None:
            width = 0.001*5*amplitude
        if head_width is None:
            head_width = 0.3*amplitude
        
        if not normalize:
            ax.add_patch(FancyArrow(origin[0], origin[1],
                                    self.vector[0]*amplitude, self.vector[1]*amplitude,
                                    width=width,
                                    head_width=head_width,
                                    length_includes_head=True,
                                    color=color))
        else:
            normalized_vector = self.copy()
            normalized_vector.Normalize()
            ax.add_patch(FancyArrow(origin[0], origin[1],
                                    normalized_vector[0]*amplitude, normalized_vector[1]*amplitude,
                                    width=width,
                                    head_width=head_width,
                                    length_includes_head=True,
                                    color=color))
            
        if line:
            style='-'+color
            linestyle = '-.'
            origin = Point2D(origin)
            p1, p2 = origin, origin+self
            u = p2 - p1
#            plt.plot([p1[0], p2[0]], [p1[1], p2[1]], style)
            p3 = p1 - 3*u
            p4 = p2 + 4*u
            ax.plot([p3[0], p4[0]], [p3[1], p4[1]], style, linestyle=linestyle)

        if label is not None:
            ax.text(*(origin+self*amplitude).vector, label)

        return fig, ax
    


X2D = Vector2D((1, 0))
Y2D = Vector2D((0, 1))


class Point2D(Vector2D):
    def __init__(self, vector:Tuple[float, float], name=''):
        Vector2D.__init__(self, vector, name)

    def __add__(self, other_vector):
        return Point2D(add2D(self.vector, other_vector.vector))

    def __neg__(self):
        return Point2D((-self.vector[0], -self.vector[1]))

    def __sub__(self, other_vector):
        return Point2D(sub2D(self.vector, other_vector.vector))

    def __mul__(self, value:float):
        return Point2D(mul2D(self.vector, value))

    def __truediv__(self, value:float):
        if value == 0:
            raise ZeroDivisionError
        return Point2D((self.vector[0] / value,
                        self.vector[1] / value))

    def To3D(self, plane_origin, vx, vy):
        return Point3D([plane_origin.vector[0] + vx.vector[0]*self.vector[0] + vy.vector[0]*self.vector[1],
                        plane_origin.vector[1] + vx.vector[1]*self.vector[0] + vy.vector[1]*self.vector[1],
                        plane_origin.vector[2] + vx.vector[2]*self.vector[0] + vy.vector[2]*self.vector[1],
                        ])

    def MPLPlot(self, ax=None, color='k'):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        x1, y1 = self.vector
        ax.plot([x1], [y1], color=color, marker='o')
        return fig, ax

    def point_distance(self, point2):
        return (self-point2).Norm()


    @classmethod
    def LinesIntersection(cls, line1, line2, curvilinear_abscissa=False):
        x1 = line1.points[0].vector[0]
        y1 = line1.points[0].vector[1]
        x2 = line1.points[1].vector[0]
        y2 = line1.points[1].vector[1]
        x3 = line2.points[0].vector[0]
        y3 = line2.points[0].vector[1]
        x4 = line2.points[1].vector[0]
        y4 = line2.points[1].vector[1]

        denominateur = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
        if math.isclose(denominateur, 0, abs_tol=1e-6):
            if not curvilinear_abscissa:
                return None
            else:
                return None, None, None
        else:
            x = (x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4)
            x = x / denominateur
            y = (x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4)
            y = y / denominateur
            if not curvilinear_abscissa:
                return cls((x,y))
            else:
                t = (x1-x3)*(y3-y4)-(y1-y3)*(x3-x4)
                t = t / denominateur
                u = (x1-x2)*(y1-y3)-(y1-y2)*(x1-x3)
                u = -u / denominateur
                return (cls((x,y)), t, u)
            
    @classmethod
    def SegmentsIntersection(cls, segment1, segment2, curvilinear_abscissa=False):
        x1 = segment1.points[0].vector[0]
        y1 = segment1.points[0].vector[1]
        x2 = segment1.points[1].vector[0]
        y2 = segment1.points[1].vector[1]
        x3 = segment2.points[0].vector[0]
        y3 = segment2.points[0].vector[1]
        x4 = segment2.points[1].vector[0]
        y4 = segment2.points[1].vector[1]

        denominateur = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)
        if math.isclose(denominateur, 0, abs_tol=1e-6):
            if not curvilinear_abscissa:
                return None
            else:
                return None, None, None
            
        t = (x1-x3)*(y3-y4)-(y1-y3)*(x3-x4)
        t = t / denominateur
        u = (x1-x2)*(y1-y3)-(y1-y2)*(x1-x3)
        u = -u / denominateur
        if (0 <= t and t <= 1) or (0 <= u and u <= 1):
            x = (x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4)
            x = x / denominateur
            y = (x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4)
            y = y / denominateur
            if not curvilinear_abscissa:
                return cls((x,y))
            else:
                return (cls((x,y)), t, u)
        else:
            if not curvilinear_abscissa:
                return None
            else:
                return None, None, None

    def plot_data(self, marker=None, color='black', size=1,
                  opacity=1, arrow=False, stroke_width=None):
        return {'type' : 'point',
                'data' : [self.vector[0], self.vector[1]],
                'color' : color,
                'marker' : marker,
                'size' : size,
                'opacity' : opacity
                }

    @classmethod
    def MiddlePoint(cls, point1, point2):
        return (point1 + point2)*0.5

    @classmethod
    def LineProjection(cls, point, line):
        p1, p2 = line.points
        n = line.NormalVector(unit=True)
        pp1 = point - p1
        return  pp1 - pp1.Dot(n)*n + p1

O2D = Point2D((0, 0))

class Vector3D(Vector):

    # _jsonschema = {
    #     "definitions": {},
    #     "$schema": "http://json-schema.org/draft-07/schema#",
    #     "type": "object",
    #     "title": "powerpack.mechanical.Vector3D Base Schema",
    #     "required": ["vector"],
    #     "properties": {
    #         'vector' : {
    #             "type" : "array",
    #             "order" : 0,
    #             "items" : {
    #                 "type" : "number",
    #                 "step" : 1,
    #                 "minimum" : -1,
    #                 "maximum" : 1
    #                 },
    #             "minItems": 3,
    #             "maxItems": 3,
    #             "examples": [[1, 0, 0]],
    #             "editable" : True,
    #             "description" : "Vector array"
    #             }
    #         }
    #     }

    def __init__(self, vector:Tuple[float, float, float], name:str=''):
        self.vector = [0, 0, 0]
#        self.vector = npy.zeros(3)
        self.vector[0] = vector[0]
        self.vector[1] = vector[1]
        self.vector[2] = vector[2]
        self.name = name

    def __add__(self, other_vector):
        return Vector3D(add3D(self.vector, other_vector.vector))

    def __neg__(self):
        return Vector3D((-self.vector[0], -self.vector[1], -self.vector[2]))

    def __sub__(self, other_vector):
        return Vector3D(sub3D(self.vector, other_vector.vector))

    def __mul__(self, value):
        return Vector3D(mul3D(self.vector, value))

    def __truediv__(self, value):
        if value == 0:
            raise ZeroDivisionError
        return Vector3D((self.vector[0] / value,
                         self.vector[1] / value,
                         self.vector[2] / value))

    def __round__(self, ndigits:int=6):
        return self.__class__((round(self.vector[0], ndigits),
                               round(self.vector[1], ndigits),
                               round(self.vector[2], ndigits)))

    def __hash__(self):
#        return int(1000*(self.vector[0]+self.vector[1]+self.vector[2]))
        return int(round(1e6*(self.vector[0]+self.vector[1]+self.vector[2])))

    def __eq__(self, other_vector:'Vector3D'):
        if other_vector.__class__.__name__ not in ['Vector3D', 'Point3D']:
            return False
        return math.isclose(self.vector[0], other_vector.vector[0], abs_tol=1e-06) \
        and math.isclose(self.vector[1], other_vector.vector[1], abs_tol=1e-06) \
        and math.isclose(self.vector[2], other_vector.vector[2], abs_tol=1e-06)

    def Dot(self, other_vector):
        v1, v2, v3 = self.vector
        ov1, ov2, ov3 = other_vector.vector
        return CVector3DDot(v1, v2, v3, ov1, ov2, ov3)

    def Cross(self, other_vector:'Vector3D') -> 'Vector3D':
        return self.__class__(vector3D_cross(self.vector, other_vector.vector))

    def Norm(self) -> float:
        vx, vy, vz = self.vector
        return CVector3DNorm(vx, vy, vz)


    def Normalize(self) -> 'Vector3D':
        """
        Normalize the vector modifying its coordinates
        """
        n = self.Norm()
        if n == 0:
            raise ZeroDivisionError

        self.vector[0] /= n
        self.vector[1] /= n
        self.vector[2] /= n

    def point_distance(self, point2:'Vector3D') -> float:
        return (self-point2).Norm()

    def Rotation(self, center:'Point3D', axis:'Vector3D', angle:float, copy:bool=True):
        """
        Rotation of angle around axis.
        Used Rodrigues Formula:
            https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
        """
#        u = self - center
#        vector2 = (math.cos(angle)*u
#                   + (1-math.cos(angle))*(u.Dot(axis))*axis
#                   + math.sin(angle)*axis.Cross(u)
#                   + center)
        vector2 = vector3D_rotation(self.vector, center.vector, axis.vector, angle)

        if copy:
            return self.__class__(vector2)
            # return Point3D(vector2)
        else:
            self.vector = list(vector2)

    def x_rotation(self, angle:float, copy:bool=True):
        """
        Rotation of angle around X axis.
        """
        cos_angle = math.cos(angle)
        sin_angle = math.sin(angle)

        y1 = cos_angle * self.vector[1] + sin_angle * self.vector[2]
        z1 = -sin_angle * self.vector[1] + cos_angle * self.vector[2]


        if copy:
            return Point3D([self.vector[0], y1, z1])
        else:
            self.vector[1] = y1
            self.vector[2] = z1

    def y_rotation(self, angle:float, copy:bool=True):
        """
        Rotation of angle around Y axis.
        """
        cos_angle = math.cos(angle)
        sin_angle = math.sin(angle)

        z1 = cos_angle * self.vector[2] + sin_angle * self.vector[0]
        x1 = -sin_angle * self.vector[2] + cos_angle * self.vector[0]


        if copy:
            return Point3D([x1, self.vector[1], z1])
        else:
            self.vector[0] = x1
            self.vector[2] = z1

    def z_rotation(self, angle:float, copy:bool=True):
        """
        Rotation of angle around Z axis.
        """
        cos_angle = math.cos(angle)
        sin_angle = math.sin(angle)

        x1 = cos_angle * self.vector[0] + sin_angle * self.vector[1]
        y1 = -sin_angle * self.vector[0] + cos_angle * self.vector[1]


        if copy:
            return Point3D([x1, y1, self.vector[2]])
        else:
            self.vector[0] = x1
            self.vector[1] = y1

    def Translation(self, offset, copy=True):
        if copy:
            return self + offset
        else:
            self.vector = [self.vector[0] + offset[0],
                           self.vector[1] + offset[1],
                           self.vector[2] + offset[2]]

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            new_vector = frame.OldCoordinates(self)
            if copy:
                return new_vector
            else:
                self.vector = new_vector.vector#.copy() #copy was adding to avoid connection point in RoundedLS

        if side == 'new':
            new_vector = frame.NewCoordinates(self)
            if copy:
                return new_vector
            else:
                self.vector = new_vector.vector#.copy()

    def PlaneProjection3D(self, plane_origin, x, y):
        z = x.Cross(y)
        z.Normalize()
        return self - z.Dot(self-plane_origin)*z

    def PlaneProjection2D(self, plane_origin, x, y):
        z = x.Cross(y)
        z.Normalize()
        p3d = self - (self-plane_origin).Dot(z)*z
        u1 = p3d.Dot(x)
        u2 = p3d.Dot(y)
        return Point2D((u1, u2))


    def To2D(self, plane_origin, x, y):
        # print(self.Dot(x))
        # print(plane_origin.Dot(x))
        # print(self.Dot(y))
        # print(plane_origin.Dot(y))
        x2d = self.Dot(x) - plane_origin.Dot(x)
        y2d = self.Dot(y) - plane_origin.Dot(y)
        return Point2D((x2d,y2d))

    def RandomUnitNormalVector(self):
        """
        Returns a random normal vector
        """
        v = Vector3D(npy.random.random(3))

        v = v - v.Dot(self)*self/(self.Norm()**2)
        v.Normalize()
        return v
    
    def deterministic_unit_normal_vector(self):
        """
        Returns a deterministic normal vector
        """
        v = X3D
        if not math.isclose(self.vector[1], 0, abs_tol=1e-7) \
        or not math.isclose(self.vector[2], 0, abs_tol=1e-7):
            v = X3D
        else:
            v = Y3D
        v = v - v.Dot(self)*self/(self.Norm()**2)
        v.Normalize()
        return v

    def Copy(self):
        return Vector3D(self.vector)
    
    @classmethod
    def random(cls, xmin, xmax, ymin, ymax, zmin, zmax):
        return cls((random.uniform(xmin, xmax),
                   random.uniform(ymin, ymax),
                   random.uniform(zmin, zmax)))

#    @classmethod
#    def DictToObject(cls, dict_):
#        return cls(dict_['vector'])

    @classmethod
    def from_step(cls, arguments, object_dict):
        if type(arguments[1]) is int:
        # VECTOR
            return cls(object_dict[arguments[1]], arguments[0][1:-1])
        else:
        # DIRECTION
            return cls([float(i)/1000 for i in arguments[1][1:-1].split(",")],
                        arguments[0][1:-1])

    def MPLPlot(self, ax=None, starting_point=None, color=''):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = ax.figure

        if starting_point is None:
            starting_point = Point3D((0,0,0))
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        xs = [starting_point[0], self.vector[0]+starting_point[0]]
        ys = [starting_point[1], self.vector[1]+starting_point[1]]
        zs = [starting_point[2], self.vector[2]+starting_point[2]]
        if color:
            a = Arrow3D(xs, ys, zs, mutation_scale=10, lw=3, arrowstyle="-|>", color=color)
        else:
            a = Arrow3D(xs, ys, zs, mutation_scale=10, lw=3, arrowstyle="-|>")
        ax.add_artist(a)
        return fig, ax


X3D = Vector3D((1, 0, 0))
Y3D = Vector3D((0, 1, 0))
Z3D = Vector3D((0, 0, 1))


class Point3D(Vector3D):
    _standalone_in_db = False
    # _jsonschema = {
    #     "definitions": {},
    #     "$schema": "http://json-schema.org/draft-07/schema#",
    #     "type": "object",
    #     "title": "powerpack.mechanical.Point3D Base Schema",
    #     "required": ["vector"],
    #     "properties": {
    #         'vector' : {
    #             "type" : "object",
    #             "order" : 0,
    #             "classes" : ["volmdlr.core.Vector3D"],
    #             "editable" : True,
    #             "description" : "Vector array"
    #             }
    #         }
    #     }

    def __init__(self, vector:Tuple[float, float, float], name:str=''):
        Vector3D.__init__(self, vector, name)

    def __add__(self, other_vector):
        return Point3D(add3D(self.vector, other_vector.vector))

    def __neg__(self):
        return Point3D((-self.vector[0], -self.vector[1], -self.vector[2]))

    def __sub__(self, other_vector):
        return Point3D(sub3D(self.vector, other_vector.vector))

    def __mul__(self, value):
        return Point3D(mul3D(self.vector, value))

    def __truediv__(self, value):
        if value == 0:
            raise ZeroDivisionError
        return Point3D((self.vector[0] / value,
                        self.vector[1] / value,
                        self.vector[2] / value))

    def Copy(self):
        return Point3D(self.vector)

    def copy(self):
        return Point3D(self.vector)

    def MPLPlot(self, ax=None, color='k'):

        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        ax.scatter(*self.vector, color=color)
        return fig, ax


    def To2D(self, plane_origin, x, y):
        x2d = self.Dot(x) - plane_origin.Dot(x)
        y2d = self.Dot(y) - plane_origin.Dot(y)
        return Point2D((x2d,y2d))

    @classmethod
    def from_step(cls, arguments, object_dict):
        return cls([float(i)/1000 for i in arguments[1][1:-1].split(",")],
                    arguments[0][1:-1])

    def babylon_script(self):
        s = 'var sphere = BABYLON.MeshBuilder.CreateSphere("point", {diameter: 0.05}, scene);\n'
        s += "sphere.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n".format(self.vector[0],self.vector[1],self.vector[2])
        s += 'var mat = new BABYLON.StandardMaterial("mat", scene);\n'
        s += 'mat.diffuseColor = new BABYLON.Color3(1, 0, 0);\n'
        s += 'sphere.material = mat;\n'
        return s

O3D = Point3D((0, 0, 0))


# =============================================================================
#  Basis, Frames
# =============================================================================
   
    
class Matrix22:
    def __init__(self, M11:float, M12:float, M21:float, M22:float):
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
        return Matrix22(self.M11*other_matrix.M11 + self.M12*other_matrix.M21,
                        self.M11*other_matrix.M12 + self.M12*other_matrix.M22,
                        self.M21*other_matrix.M11 + self.M22*other_matrix.M21,
                        self.M21*other_matrix.M12 + self.M22*other_matrix.M22)
    
    def determinent(self):
        return self.M11 * self.M22 - self.M12 * self.M21
    
    def inverse(self):
        det = self.determinent()
        if not math.isclose(det, 0, abs_tol=1e-10):
            det_inv = 1/self.determinent()
            return Matrix22(det_inv*self.M22, -det_inv*self.M12,
                            -det_inv*self.M21, det_inv*self.M11)
        else:
            raise ValueError('The matrix is singular')
    
    def vector_multiplication(self, vector):
        return vector.__class__((self.M11*vector[0] + self.M12*vector[1],
                                 self.M21*vector[0] + self.M22*vector[1]))


class Matrix33:
    def __init__(self, M11:float, M12:float, M13:float,
                       M21:float, M22:float, M23:float,
                       M31:float, M32:float, M33:float):
        self.M11 = M11
        self.M12 = M12
        self.M13 = M13
        self.M21 = M21
        self.M22 = M22
        self.M23 = M23
        self.M31 = M31
        self.M32 = M32
        self.M33 = M33

#    def __getitem__(self, key):
#        return self.vector[key]

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
        M11, M12, M13, M21, M22, M23, M31, M32, M33 = C_matrix_multiplication3(self.M11, self.M12, self.M13,
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
        return Matrix33(self.M11*float_value, self.M12*float_value, self.M13*float_value,
                        self.M21*float_value, self.M22*float_value, self.M23*float_value,
                        self.M31*float_value, self.M32*float_value, self.M33*float_value)


    def vector_multiplication(self, vector):
        v1, v2, v3 = vector.vector
        u1, u2, u3 = C_matrix_vector_multiplication3(self.M11, self.M12, self.M13,
                                                     self.M21, self.M22, self.M23,
                                                     self.M31, self.M32, self.M33,
                                                     v1, v2, v3)
        
        return vector.__class__((u1, u2, u3))

    def determinent(self):
        det = self.M11*self.M22*self.M33 + self.M12*self.M23*self.M31 \
            + self.M13*self.M21*self.M32 - self.M13*self.M22*self.M31 \
            - self.M23*self.M32*self.M11 - self.M33*self.M12*self.M21
        return det

    def inverse(self):
        det = self.determinent()

        if not math.isclose(det, 0, abs_tol=1e-10):
            det_inv = 1/det
            return Matrix33(det_inv*(self.M22*self.M33 - self.M23*self.M32),# a22a33−a23a32
                            det_inv*(self.M13*self.M32 - self.M12*self.M33),# a13a32−a12a33
                            det_inv*(self.M12*self.M23 - self.M13*self.M22),# a12a23−a13a22
                            det_inv*(self.M23*self.M31 - self.M21*self.M33),# a23a31−a21a33
                            det_inv*(self.M11*self.M33 - self.M13*self.M31),# a11a33−a31a13
                            det_inv*(self.M21*self.M13 - self.M23*self.M11),# a13a21−a23a11
                            det_inv*(self.M21*self.M32 - self.M31*self.M22),# a21a32−a31a22
                            det_inv*(self.M12*self.M31 - self.M32*self.M11),# a12a31−a32a11
                            det_inv*(self.M11*self.M22 - self.M21*self.M12) # a11a22−a21a12
                            )
        else:
            # print(self.__dict__, det)
            raise ValueError('The matrix is singular')

    @classmethod
    def random_matrix(cls, minimum=0, maximum=1):
        range_ = maximum - minimum
        return cls(*[minimum + range_*random.random() for _ in range(9)])

    def to_numpy(self):
        return npy.array([[self.M11, self.M12, self.M13],
                          [self.M21, self.M22, self.M23],
                          [self.M31, self.M32, self.M33]])

class Basis(DessiaObject):
    """
    Abstract class of a basis
    """
    def __getitem__(self, key):
        return self.vectors[key]

    def __setitem__(self, key, item):
        self.vectors[key] = item

    def __contains__(self, vector):
        return vector in self.vectors

    def __eq__(self, other_basis):
        all_equal = all([other_vector == vector\
                         for other_vector, vector\
                         in zip(other_basis.vectors, self.vectors)])
        return all_equal

    def __hash__(self):
        return hash(self.vectors)

    def copy(self):
        return self.__class__(*self.vectors)

class Basis2D(Basis):
    """
    Defines a 2D basis
    :param u:Vector2D: first vector of the basis
    :param v:Vector2D: second vector of the basis
    """
    def __init__(self, u:Vector2D, v:Vector2D, name:str=''):
        self.u = u
        self.v = v
        self.name = name

    def __neg__(self):
        Pinv = self.InverseTransferMatrix()
        return Basis2D(Vector3D(Pinv[:, 0]),
                       Vector3D(Pinv[:, 1]))

    def __repr__(self):
        return '{}: U={}, V={}'.format(self.__class__.__name__, *self.vectors)

    def _get_vectors(self):
        return (self.u, self.v)

    vectors = property(_get_vectors)


    def to_frame(self, origin:Point3D) -> 'Frame3D':
        return Frame2D(origin, self.u, self.v)


    def TransferMatrix(self):
        return npy.array([[self.u[0], self.v[0]],
                          [self.u[1], self.v[1]]])

    def InverseTransferMatrix(self):
        det = self.u[0]*self.v[1] - self.v[0]*self.u[1]
        if not math.isclose(det, 0, abs_tol=1e-10):
            return 1/det * npy.array([[self.v[1], -self.v[0]],
                                     [-self.u[1], self.u[0]]])
        else:
            raise ZeroDivisionError

    def NewCoordinates(self, vector):
        matrix = self.InverseTransferMatrix()
        return Point2D((matrix[0][0]*vector[0] + matrix[0][1]*vector[1],
                         matrix[1][0]*vector[0] + matrix[1][1]*vector[1]))

    def OldCoordinates(self, vector):
        matrix = self.TransferMatrix()
        return Point2D((matrix[0][0]*vector[0] + matrix[0][1]*vector[1],
                         matrix[1][0]*vector[0] + matrix[1][1]*vector[1]))

    def Rotation(self, angle:float, copy=True):
        center = O2D
        new_u = self.u.Rotation(center, angle, True)
        new_v = self.v.Rotation(center, angle, True)

        if copy:
            return Basis2D(new_u, new_v)
        self.u = new_u
        self.v = new_v

    def Copy(self):
        return Basis2D(self.u, self.v)

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
    def __init__(self, u:Vector3D, v:Vector3D, w:Vector3D, name:str=''):
        self.u = u
        self.v = v
        self.w = w
        self.name = name

    def __hash__(self):
        return hash(self.u) + hash(self.v) + hash(self.w)

    def __add__(self, other_basis):
        M = self.TransferMatrix()*other_basis.TransferMatrix()
        return Basis3D(Vector3D((M.M11, M.M21, M.M31)),
                       Vector3D((M.M12, M.M22, M.M32)),
                       Vector3D((M.M13, M.M23, M.M33)))


    def __neg__(self):
        M = self.InverseTransferMatrix()
        return Basis3D(Vector3D((M.M11, M.M21, M.M31)),
                       Vector3D((M.M12, M.M22, M.M32)),
                       Vector3D((M.M13, M.M23, M.M33)))

    def __sub__(self, other_frame):
        P1inv = other_frame.InverseTransferMatrix()
        P2 = self.TransferMatrix()
        M = P1inv * P2
        return Basis3D(Vector3D((M.M11, M.M21, M.M31)),
                       Vector3D((M.M12, M.M22, M.M32)),
                       Vector3D((M.M13, M.M23, M.M33)))

    def __round__(self, ndigits:int=6):
        return self.__class__((round(self.u, ndigits),
                               round(self.v, ndigits),
                               round(self.w, ndigits)))

    def __repr__(self):
        return '{}: U={}, V={}, W={}'.format(self.__class__.__name__, *self.vectors)

    def _get_vectors(self):
        return (self.u, self.v, self.w)

    vectors = property(_get_vectors)

    # TODO: transform to annotation when available
    @classmethod
    def from_two_vectors(cls, vector1:Vector3D, vector2:Vector3D) -> 'Basis3D':
        """
        Create a basis with first vector1 adimensionned, as u, v is the vector2 substracted of u component,
        w is the cross product of u and v
        """
        u = vector1.copy()
        u.Normalize()
        v = vector2 - vector2.Dot(vector1)*vector1
        v.Normalize()
        w = u.Cross(v)

        return Basis3D(u, v, w)

    def to_frame(self, origin):
        return Frame3D(origin, self.u, self.v, self.w)

    def Rotation(self, axis:Vector3D, angle:float, copy:bool=True):
        center = O3D
        new_u = self.u.Rotation(center, axis, angle, True)
        new_v = self.v.Rotation(center, axis, angle, True)
        new_w = self.w.Rotation(center, axis, angle, True)

        if copy:
            return Basis3D(new_u, new_v, new_w, self.name)
        else:
            self.u = new_u
            self.v = new_v
            self.w = new_w

    def x_rotation(self, angle:float, copy:bool=True):
        new_u = self.u.x_rotation(angle, True)
        new_v = self.v.x_rotation(angle, True)
        new_w = self.w.x_rotation(angle, True)

        if copy:
            return Basis3D(new_u, new_v, new_w, self.name)
        else:
            self.u = new_u
            self.v = new_v
            self.w = new_w

    def y_rotation(self, angle:float, copy:bool=True):
        new_u = self.u.y_rotation(angle, True)
        new_v = self.v.y_rotation(angle, True)
        new_w = self.w.y_rotation(angle, True)

        if copy:
            return Basis3D(new_u, new_v, new_w, self.name)
        else:
            self.u = new_u
            self.v = new_v
            self.w = new_w

    def z_rotation(self, angle:float, copy:bool=True):
        new_u = self.u.z_rotation(angle, True)
        new_v = self.v.z_rotation(angle, True)
        new_w = self.w.z_rotation(angle, True)

        if copy:
            return Basis3D(new_u, new_v, new_w, self.name)
        else:
            self.u = new_u
            self.v = new_v
            self.w = new_w

    def EulerRotation(self, angles:Tuple[float, float, float], copy:bool=True):
        psi, theta, phi = angles
        center = O3D

        vect_u = self.u.Copy()
        vect_v = self.v.Copy()
        vect_w = self.w.Copy()

        # Rotation around w
        vect_u.Rotation(center, vect_w, psi, False)
        vect_v.Rotation(center, vect_w, psi, False)

        # Rotation around v
        vect_v.Rotation(center, vect_u, theta, False)
        vect_w.Rotation(center, vect_u, theta, False)

        # Rotation around w
        vect_u.Rotation(center, vect_w, phi, False)
        vect_v.Rotation(center, vect_w, phi, False)

        if copy:
            return Basis3D(vect_u, vect_v, vect_w)
        self.u = vect_u
        self.v = vect_v
        self.w = vect_w

    def TransferMatrix(self):
        return Matrix33(self.u.vector[0], self.v.vector[0], self.w.vector[0],
                        self.u.vector[1], self.v.vector[1], self.w.vector[1],
                        self.u.vector[2], self.v.vector[2], self.w.vector[2])

    def InverseTransferMatrix(self):
        return self.TransferMatrix().inverse()

    def NewCoordinates(self, vector):
        matrix = self.InverseTransferMatrix()
        return matrix.vector_multiplication(vector)

    def OldCoordinates(self, point):
        matrix = self.TransferMatrix()
        return matrix.vector_multiplication(point)

    def copy(self):
        return Basis3D(self.u, self.v, self.w)
    
    
class Frame2D(Basis2D):
    """
    Defines a 2D basis
    :param origin:Point2D: origin of the basis
    :param u:Vector2D: first vector of the basis
    :param v:Vector2D: second vector of the basis
    """
    def __init__(self, origin:Point2D, u:Vector2D, v:Vector2D, name:str=''):
        self.origin = origin
        Basis2D.__init__(self, u, v, name=name)

    def __repr__(self):
        return '{}: O={} U={}, V={}'.format(self.__class__.__name__, self.origin, self.u, self.v)

    def __neg__(self):
        Pinv = self.InverseTransferMatrix()
        new_origin = Point2D(npy.dot(Pinv, self.origin.vector))
        return Frame2D(new_origin,
                       Vector2D(Pinv[:, 0]),
                       Vector2D(Pinv[:, 1]))


    def __add__(self, other_frame):
        P1 = self.TransferMatrix()
        new_origin = Point2D(npy.dot(P1, other_frame.origin.vector) + self.origin.vector)
        M = npy.dot(P1, other_frame.TransferMatrix())
        return Frame2D(new_origin,
                       Vector2D(M[:, 0]),
                       Vector2D(M[:, 1]))


    def __sub__(self, other_frame):
        P1inv = other_frame.InverseTransferMatrix()
        P2 = self.TransferMatrix()
        new_origin = Point2D(npy.dot(P1inv, (self.origin - other_frame.origin).vector))
        M = npy.dot(P1inv, P2)
        return Frame2D(new_origin,
                       Vector2D(M[:, 0]),
                       Vector2D(M[:, 1]))

    def Basis(self):
        return Basis2D(self.u, self.v)

    def NewCoordinates(self, vector):
        return Basis2D.NewCoordinates(self, vector - self.origin)

    def OldCoordinates(self, vector):
        return Basis2D.OldCoordinates(self, vector) + self.origin

    def Translation(self, vector, copy=True):
        new_origin = self.origin.Translation(vector, True)
        if copy:
            return Frame2D(new_origin, self.u, self.v)
        self.origin = new_origin


    def Rotation(self, angle, copy=True):
        new_base = Basis2D.Rotation(self, angle, True)
        if copy:
            new_frame = Frame2D(self.origin, new_base.u, new_base.v)
            return new_frame
        self.u = new_base.u
        self.v = new_base.v

    def Draw(self, ax=None, style='ok'):
        if ax is None:
            fig, ax = plt.subplots()

        ax.plot(*self.origin.vector, style)
        self.u.plot(origin=self.origin, ax=ax, color='r')
        self.v.plot(origin=self.origin, ax=ax, color='g')
        ax.axis('equal')

    def Copy(self):
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
    def __init__(self, origin:Point3D, u:Vector3D, v:Vector3D, w:Vector3D, name:str=''):
        self.origin = origin
        Basis3D.__init__(self, u, v, w)
        self.name = name

    def __repr__(self):
        return '{}: O={} U={}, V={}, W={}'.format(self.__class__.__name__,
                                                  self.origin,
                                                  self.u, self.v, self.w)


    def __neg__(self):
        M = self.InverseTransferMatrix()
        new_origin = M.vector_multiplication(self.origin)
        return Frame3D(new_origin,
                       Vector3D((M.M11, M.M21, M.M31)),
                       Vector3D((M.M12, M.M22, M.M32)),
                       Vector3D((M.M13, M.M23, M.M33)))


    def __add__(self, other_frame):
        P1 = self.TransferMatrix()
        new_origin = P1.vector_multiplication(other_frame.origin) + self.origin


        M = P1 * other_frame.TransferMatrix()
        return Frame3D(new_origin,
                       Vector3D((M.M11, M.M21, M.M31)),
                       Vector3D((M.M12, M.M22, M.M32)),
                       Vector3D((M.M13, M.M23, M.M33)))


    def __sub__(self, other_frame):
        P1inv = other_frame.InverseTransferMatrix()
        P2 = self.TransferMatrix()
        new_origin = P1inv.vector_multiplication(self.origin - other_frame.origin)
        M = P1inv * P2
        return Frame3D(new_origin,
                       Vector3D((M.M11, M.M21, M.M31)),
                       Vector3D((M.M12, M.M22, M.M32)),
                       Vector3D((M.M13, M.M23, M.M33)))

    def __round__(self, ndigits=6):
        return self.__class__(round(self.origin, ndigits),
                              round(self.u, ndigits),
                              round(self.v, ndigits),
                              round(self.w, ndigits))

    def __hash__(self):
        return hash(self.u) + hash(self.v) + hash(self.w) + hash(self.origin)

    def Basis(self):
        return Basis3D(self.u, self.v, self.w)

    def NewCoordinates(self, vector):
        """ You have to give coordinates in the global landmark """
        return Basis3D.NewCoordinates(self, vector - self.origin)

    def OldCoordinates(self, vector): 
        """ You have to give coordinates in the local landmark """
        return Basis3D.OldCoordinates(self, vector) + self.origin

    def Rotation(self, axis, angle, copy=True):
        new_base = Basis3D.Rotation(self, axis, angle, copy=True)
        if copy:
            new_frame = Frame3D(self.origin.copy(), new_base.u, new_base.v, new_base.w, self.name)
            return new_frame
        self.u = new_base.u
        self.v = new_base.v
        self.w = new_base.w

    def Translation(self, offset, copy=True):
        if copy:
            return Frame3D(self.origin.Translation(offset, copy=True), self.u, self.v, self.w, self.name)
        self.origin.Translation(offset, copy=False)

    def copy(self):
        return Frame3D(self.origin.Copy(), self.u.Copy(), self.v.Copy(), self.w.Copy())

    def plot2d(self, x=X3D, y=Y3D, ax=None, color='k'):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure

        origin2d = self.origin.To2D(O3D, x, y)

        for iv, vector in enumerate(self.vectors):
            vector2D = vector.To2D(O3D, x, y)
            if vector2D.Norm() > 1e-8:
                vector2D.plot(origin=origin2d, ax=ax, color=color, label=str(iv+1))

        return fig, ax


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
            w = u.Cross(v)
        return cls(origin, u, v, w, arguments[0][1:-1])
    
    # @classmethod
    # def from_step(cls, arguments, object_dict):
    #     origin = object_dict[arguments[1]]
    #     if arguments[2] == '$':
    #         w = None
    #     else:
    #         w = object_dict[arguments[2]]
    #     if arguments[3] == '$':
    #         u = None
    #     else:
    #         u = object_dict[arguments[3]]
    #     if u is None or w is None:
    #         v = None
    #     else:
    #         v = w.Cross(u)
    #     return cls(origin, u, v, w, arguments[0][1:-1])

    def babylonjs(self, size=0.1, parent=None):
        s = 'var origin = new BABYLON.Vector3({},{},{});\n'.format(*self.origin)
        s += 'var o_u = new BABYLON.Vector3({}, {}, {});\n'.format(*(size*self.u+self.origin))
        s += 'var o_v = new BABYLON.Vector3({}, {}, {});\n'.format(*(size*self.v+self.origin))
        s += 'var o_w = new BABYLON.Vector3({}, {}, {});\n'.format(*(size*self.w+self.origin))
        s += 'var line1 = BABYLON.MeshBuilder.CreateTube("frame_U", {{path: [origin, o_u], radius: {}}}, scene);'.format(0.03*size)
        s += 'line1.material = red_material;\n'
        s += 'var line2 = BABYLON.MeshBuilder.CreateTube("frame_V", {{path: [origin, o_v], radius: {}}}, scene);'.format(0.03*size)
        s += 'line2.material = green_material;\n'
        s += 'var line3 = BABYLON.MeshBuilder.CreateTube("frame_W", {{path: [origin, o_w], radius: {}}}, scene);'.format(0.03*size)
        s += 'line3.material = blue_material;\n'
        if parent is not None:
            s += 'line1.parent = {};\n'.format(parent)
            s += 'line2.parent = {};\n'.format(parent)
            s += 'line3.parent = {};\n'.format(parent)

        return s
    








