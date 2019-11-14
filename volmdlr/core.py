#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:07:37 2017

@author: steven
"""


#import bezier
from geomdl import NURBS



import warnings
import math
import numpy as npy
npy.seterr(divide='raise')
#from itertools import permutations

import matplotlib.pyplot as plt
import  mpl_toolkits
from matplotlib.patches import Arc, FancyArrow
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import proj3d
from matplotlib.patches import FancyArrowPatch

import networkx as nx

from .vmcy import (sub2D, add2D, mul2D, Vector2DNorm, Vector2DDot, 
                   sub3D, add3D, mul3D, Vector3DNorm, Vector3DDot, 
                   vector3D_cross, vector3D_rotation,
                   LineSegment2DPointDistance, PolygonPointBelongs)

from scipy.linalg import solve, LinAlgError

import volmdlr.geometry as geometry
from volmdlr import plot_data
from volmdlr import triangulation as tri

from jinja2 import Environment, PackageLoader, select_autoescape

import webbrowser
import os

import tempfile
import subprocess

import random



def standardize_knot_vector(knot_vector):
    u0 = knot_vector[0]
    u1 = knot_vector[-1]
    standard_u_knots = []
    if u0 != 0 or u1 != 1:
        x = 1/(u1-u0)
        y = u0/(u0-u1)
        for u in knot_vector:
            standard_u_knots.append(u*x+y)
        return standard_u_knots
    else:
        return knot_vector


def find_and_replace(string, find, replace):
    """
    Finds a string in a string and replace it
    """
    index = string.find(find)
    if index != -1:
        try:
            # verifie si il reste pas des chiffre apres
            int(string[index+len(find)])
        except (ValueError, IndexError):
            # on remplace
            return string[:index]+replace+string[index+len(find):]
        else:
            return string[:index]+find_and_replace(string[index+len(find)], find, replace)
    return string

def step_split_arguments(function_arg):
    """
    Split the arguments of a function that doesn't start with '(' but end with ')'
    ex: IN: '#123,#124,#125)'
       OUT: ['#123', '#124', '#125']
    """
    if len(function_arg) > 0 and function_arg[-1] != ')':
            function_arg += ')'
    arguments = []
    argument = ""
    parenthesis = 1
    for char in function_arg:
        if char == "(":
            parenthesis += 1

        if char != "," or parenthesis > 1:
            argument += char
        else:
            arguments.append(argument)
            argument = ""

        if char == ")":
            parenthesis -= 1
            if parenthesis == 0:
                arguments.append(argument[:-1])
                argument = ""
                break
    return arguments

def set_to_list(step_set):
    char_list = step_set.split(',')
    char_list[0] = char_list[0][1:]
    char_list[-1] = char_list[-1][:-1]
    return [elem for elem in char_list]

def delete_node_and_predecessors(graph, node):
    predecessors = list(graph.predecessors(node))
    graph.remove_node(node)
    for predecessor in predecessors:
        delete_node_and_predecessors(graph, predecessor)
        
def delete_node_and_successors(graph, node):
    successors = list(graph.successors(node))
    graph.remove_node(node)
    for successor in successors:
        delete_node_and_successors(graph, successor)

class Matrix22:
    def __init__(self, M11, M12, M21, M22):
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


class Matrix33:
    def __init__(self, M11, M12, M13, M21, M22, M23, M31, M32, M33):
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
        return Matrix33(self.M11*other_matrix.M11 + self.M12*other_matrix.M21 + self.M13*other_matrix.M31,
                        self.M11*other_matrix.M12 + self.M12*other_matrix.M22 + self.M13*other_matrix.M32,
                        self.M11*other_matrix.M13 + self.M12*other_matrix.M23 + self.M13*other_matrix.M33,
                        self.M21*other_matrix.M11 + self.M22*other_matrix.M21 + self.M23*other_matrix.M31,
                        self.M21*other_matrix.M12 + self.M22*other_matrix.M22 + self.M23*other_matrix.M32,
                        self.M21*other_matrix.M13 + self.M22*other_matrix.M23 + self.M23*other_matrix.M33,
                        self.M31*other_matrix.M11 + self.M32*other_matrix.M21 + self.M33*other_matrix.M31,
                        self.M31*other_matrix.M12 + self.M32*other_matrix.M22 + self.M33*other_matrix.M32,
                        self.M31*other_matrix.M13 + self.M32*other_matrix.M23 + self.M33*other_matrix.M33)

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
        return vector.__class__((self.M11*vector[0] + self.M12*vector[1] + self.M13*vector[2],
                                 self.M21*vector[0] + self.M22*vector[1] + self.M23*vector[2],
                                 self.M31*vector[0] + self.M32*vector[1] + self.M33*vector[2]))

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
            print(self.__dict__, det)
            raise ValueError

    @classmethod
    def random_matrix(cls, minimum=0, maximum=1):
        range_ = maximum - minimum
        return cls(*[minimum + range_*random.random() for _ in range(9)])

    def to_numpy(self):
        return npy.array([[self.M11, self.M12, self.M13],
                          [self.M21, self.M22, self.M23],
                          [self.M31, self.M32, self.M33]])

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)


class Vector:
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

    def __ne__(self, other_vector):
        return not npy.allclose(self.vector, other_vector.vector)

    def to_numpy(self):
        return npy.array(self.vector)

    def copy(self):
        return self.__class__(self.vector)

    def Dict(self):
        d = {'vector': [float(i) for i in self.vector]}
        return d

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
    def __init__(self, vector, name=''):
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
        
    def __mul__(self, value):
        return Vector2D(mul2D(self.vector, value))

    def __truediv__(self, value):
        if value == 0:
            raise ZeroDivisionError
        return Vector2D((self.vector[0] / value,
                         self.vector[1] / value))

    def __round__(self, ndigits=6):
        return self.__class__((round(self.vector[0], ndigits),
                               round(self.vector[1], ndigits)))
    
    def __hash__(self):
        return int(1000*(self.vector[0]+self.vector[1]))
    
    def __eq__(self, other_vector):
        return math.isclose(self.vector[0], other_vector.vector[0], abs_tol=1e-08) \
        and math.isclose(self.vector[1], other_vector.vector[1], abs_tol=1e-08)

    def Norm(self):
        """
        :returns: norm of vector
        """
        return Vector2DNorm(self.vector)
    
    def Normalize(self):
        """
        Normalize the vector modifying it's coordinate
        """
        n = self.Norm()
        if math.isclose(n, 0, abs_tol=1e-9):
            raise ZeroDivisionError

        self.vector[0] /= n
        self.vector[1] /= n

    def Dot(self, other_vector):
        return Vector2DDot(self.vector, other_vector.vector)

    def Cross(self, other_vector):
        u1, u2 = self.vector
        v1, v2 = other_vector.vector
        return u1*v2 - u2*v1


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
    
    def Draw(self):
        warnings.warn(
            "Draw is deprecated and will be removed in next versions, use plot() instead",
            DeprecationWarning
        )
        self.plot()

    def plot(self, amplitude=0.5, origin=None, ax=None, color='k', line=False, label=None):
        if origin is None:
            origin = Vector2D((0., 0.))
        
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
        ax.add_patch(FancyArrow(origin[0], origin[1],
                                self.vector[0]*amplitude, self.vector[1]*amplitude,
                                width=0.001,
                                head_width=0.01,
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

#        print(origin, self, origin+self)
        if label is not None:
            ax.text(*(origin+self*amplitude).vector, label)

        return fig, ax

    @classmethod
    def DictToObject(cls, dict_):
        return cls(dict_['vector'])


x2D = Vector2D((1, 0))
y2D = Vector2D((0, 1))

X2D = Vector2D((1, 0))
Y2D = Vector2D((0, 1))

class Point2D(Vector2D):
    def __init__(self, vector, name=''):
        Vector2D.__init__(self, vector)
        self.name = name

    def __add__(self, other_vector):
        return Point2D(add2D(self.vector, other_vector.vector))

    def __neg__(self):
        return Point2D((-self.vector[0], -self.vector[1]))

    def __sub__(self, other_vector):
        return Point2D(sub2D(self.vector, other_vector.vector))

    def __mul__(self, value):
        return Point2D(mul2D(self.vector, value))

    def __truediv__(self, value):
        if value == 0:
            raise ZeroDivisionError
        return Point2D((self.vector[0] / value,
                        self.vector[1] / value))

    def To3D(self, plane_origin, vx, vy):
        return Point3D([plane_origin.vector[0] + vx.vector[0]*self.vector[0] + vy.vector[0]*self.vector[1],
                        plane_origin.vector[1] + vx.vector[1]*self.vector[0] + vy.vector[1]*self.vector[1],
                        plane_origin.vector[2] + vx.vector[2]*self.vector[0] + vy.vector[2]*self.vector[1],
                        ])

    def MPLPlot(self, ax=None, style='ob'):
        if ax is None:
            fig, ax = plt.subplots()
        x1 = self.vector
        ax.plot([x1[0]], [x1[1]], style)
        return ax

    def PointDistance(self, point2):
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
#        d = (p2-p1) / p2.PointDistance(p1)
#        n = d.Rotation(Point2D((0, 0)), math.pi/2).vector
        n = line.NormalVector(unit=True)
        pp1 = point - p1
        return  pp1 - pp1.Dot(n)*n + p1

o2D = Point2D((0, 0))
O2D = Point2D((0, 0))

class Basis:
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

    def Dict(self):
        d = {'vectors' : [vector.Dict() for vector in self.vectors]}
        return d

    def copy(self):
        return self.__class__(*self.vectors)

class Basis2D(Basis):
    """
    Defines a 2D basis
    :param u: first vector of the basis
    :param v: second vector of the basis
    """
    def __init__(self, u, v):
        self.u = u
        self.v = v
    
    def __neg__(self):
        Pinv = self.InverseTransferMatrix()
        return Basis2D(Vector3D(Pinv[:, 0]),
                       Vector3D(Pinv[:, 1]))

    def __repr__(self):
        return '{}: U={}, V={}'.format(self.__class__.__name__, *self.vectors)

    def _get_vectors(self):
        return (self.u, self.v)

    vectors = property(_get_vectors)
    

    def to_frame(self, origin):
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

    def Rotation(self, angle, copy=True):
        center = o2D
        new_u = self.u.Rotation(center, angle, True)
        new_v = self.v.Rotation(center, angle, True)

        if copy:
            return Basis2D(new_u, new_v)
        self.u = new_u
        self.v = new_v

    def Copy(self):
        return Basis2D(self.u, self.v)

    @classmethod
    def DictToObject(cls, dict_):
        vectors = [Vector2D.DictToObject(vector_dict) for vector_dict in dict_['vectors']]
        return cls(*vectors)

xy = Basis2D(x2D, y2D)

class Frame2D(Basis2D):
    """
    Defines a 2D basis
    :param origin: origin of the basis
    :param u: first vector of the basis
    :param v: second vector of the basis
    """
    def __init__(self, origin, u, v):
        self.origin = origin
        Basis2D.__init__(self, u, v)

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

oxy = Frame2D(o2D, x2D, y2D)

class Primitive2D:
    def __init__(self, name=''):
        self.name = name

class CompositePrimitive2D(Primitive2D):
    """
    A collection of simple primitives
    """
    def __init__(self, primitives, name=''):
        Primitive2D.__init__(self, name)
        self.primitives = primitives
        self.UpdateBasisPrimitives()

    def UpdateBasisPrimitives(self):
        basis_primitives = []
        for primitive in self.primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.basis_primitives)
            else:
                basis_primitives.append(primitive)

        self.basis_primitives = basis_primitives


    def Rotation(self, center, angle, copy=True):
        if copy:
            return self.__class__([p.Rotation(center, angle, copy=True)\
                                   for p in self.primitives])
        else:
            for p in self.basis_primitives:
                p.Rotation(center, angle, copy=False)
            self.UpdateBasisPrimitives()

    def Translation(self, offset, copy=True):
        if copy:
            return self.__class__([p.Translation(offset, copy=True)\
                                   for p in self.primitives])
        else:
            for p in self.basis_primitives:
                p.Translation(offset, copy=False)
            self.UpdateBasisPrimitives()

    def To3D(self, plane_origin, x, y, name=None):
        if name is None:
            name = '3D of {}'.format(self.name)
        primitives3D = [p.To3D(plane_origin, x, y) for p in self.primitives]
        return CompositePrimitive3D(primitives3D, name)

    # TODO: change style to color!
    def MPLPlot(self, ax=None, color='k', arrow=False, width=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = None

        for element in self.basis_primitives:
            if element.__class__.__name__ == 'LineSegment2D':
                element.MPLPlot(ax, color, arrow, width)
            else:
                element.MPLPlot(ax, color=color)

        ax.margins(0.1)
        plt.show()

        return fig, ax
    
    def plot_data(self, name, fill=None, color='black', stroke_width=0.2, opacity=1):
        plot_data = {}
        plot_data['fill'] = fill
        plot_data['name'] = name
        plot_data['type'] = 'contour'
        plot_data['plot_data'] = []
        for item in self.basis_primitives:
            plot_data['plot_data'].append(item.plot_data(color=color,
                                                        stroke_width=stroke_width,
                                                        opacity=opacity))
        return plot_data


class Wire2D(CompositePrimitive2D):
    """
    A collection of simple primitives, following each other making a wire
    """
    def __init__(self, primitives, name=''):
        CompositePrimitive2D.__init__(self, primitives, name)

    # TODO: method to check if it is a wire

    def Length(self):
        length = 0.
        for primitive in self.basis_primitives:
            length += primitive.Length()
        return length

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        length = 0.
        for primitive in self.basis_primitives:
            primitive_length = primitive.Length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(curvilinear_abscissa - length)
            length += primitive_length
        return ValueError
    
    def plot_data(self, name='', fill=None, color='black', stroke_width=1, opacity=1):
        plot_data = {}
        plot_data['name'] = name
        plot_data['type'] = 'wire'
        plot_data['plot_data'] = []
        for item in self.basis_primitives:
            plot_data['plot_data'].append(item.plot_data(color=color,
                                                        stroke_width=stroke_width,
                                                        opacity=opacity))
        return plot_data


class Contour2D(Wire2D):
    """
    A collection of 2D primitives forming a closed wire2D
    """
    def __init__(self, primitives, name=''):
        Wire2D.__init__(self, primitives, name)

    def To3D(self, plane_origin, x, y, name=None):
        if name is None:
            name = '3D of {}'.format(self.name)
        primitives3D = [p.To3D(plane_origin, x, y) for p in self.primitives]
        return Contour3D(primitives3D, name=name)

    def Area(self):
        if len(self.basis_primitives) == 1:
            return self.basis_primitives[0].Area()
        arcs = []
        points_polygon = []
        for primitive in self.basis_primitives:
            if primitive.__class__.__name__ == 'LineSegment2D':
                points_polygon.extend(primitive.points)
            elif primitive.__class__.__name__ == 'Arc2D':
                points_polygon.append(primitive.center)
                arcs.append(primitive)
        polygon = Polygon2D(points_polygon)
        A = polygon.Area()

        for arc in arcs:
            if polygon.PointBelongs(arc.interior):
                A -= arc.Area()
            else:
                A += arc.Area()

        return A

    def CenterOfMass(self):
        if len(self.basis_primitives) == 1:
            return self.basis_primitives[0].CenterOfMass()

        arcs = []
        points_polygon = []
        for primitive in self.basis_primitives:
            if primitive.__class__.__name__ == 'LineSegment2D':
                points_polygon.extend(primitive.points)
            elif primitive.__class__.__name__ == 'Arc2D':
                points_polygon.append(primitive.center)
                arcs.append(primitive)
        polygon = Polygon2D(points_polygon)

        area = polygon.Area()
        if area > 0.:
            c = area*polygon.CenterOfMass()
        else:
            c = o2D

        for arc in arcs:
            arc_area = arc.Area()
            if polygon.PointBelongs(arc.interior):
                c -= arc_area*arc.CenterOfMass()
                area -= arc_area
            else:
                c += arc_area*arc.CenterOfMass()
                area += arc_area
        return c/area

    def SecondMomentArea(self, point):
        if len(self.primitives) == 1:
            return self.primitives[0].SecondMomentArea(point)

        arcs = []
        points_polygon = []
        for primitive in self.primitives:
            if primitive.__class__.__name__ == 'Line2D':
                points_polygon.extend(primitive.points)
            elif primitive.__class__.__name__ == 'Arc2D':
                points_polygon.append(primitive.center)
                arcs.append(primitive)
        polygon = Polygon2D(points_polygon)
        A = polygon.SecondMomentArea(point)
        for arc in arcs:
            if polygon.PointBelongs(arc.middle):
                A -= arc.SecondMomentArea(point)
            else:
                A += arc.SecondMomentArea(point)
        return A

    def plot_data(self, name='', fill=None, color='black', stroke_width=1, opacity=1):
        plot_data = {}
        plot_data['fill'] = fill
        plot_data['name'] = name
        plot_data['type'] = 'contour'
        plot_data['plot_data'] = []
        for item in self.basis_primitives:
            print(item)
            plot_data['plot_data'].append(item.plot_data(color=color,
                                                        stroke_width=stroke_width,
                                                        opacity=opacity))
        return plot_data

class Mesh2D:
    def __init__(self, contours, points_densities, default_density):
        self.contours = contours
        self.points_densities = points_densities
        self.default_density = default_density

    def GeoScript(self, filepath=''):
        s = ''
        ipt = 1# point index
        ipr = 1# primitive index
        points_index = {}
        #assigning an index to point
        for contour in self.contours:
            for primitive in contour.primitives:
                for point in primitive.geo_points:
                    try:
                        points_index[point]
                    except KeyError:
                        points_index[point]=ipt
                        try:
                            d=self.points_densities[point]
                        except KeyError:
                            d=self.default_density
                        s += 'Point({})={{{},{},0.,{}}};\n'.format(ipt,*point.vector,d)
                        ipt += 1
        contours_indices = []
        for contour in self.contours:
            contour_iprs = []
            for primitive in contour.primitives:
                spr,ipr2=primitive.GeoScript(ipr,[points_index[p] for p in primitive.geo_points])
                s+=spr
                contour_iprs.extend(range(ipr,ipr2))
                ipr=ipr2
            s+='Line Loop({}) = {{{}}};\n'.format(ipr,str(contour_iprs)[1:-1])
            contours_indices.append(ipr)
            ipr+=1
        s+='Plane Surface({}) = {{{}}};\n'.format(ipr,str(contours_indices)[1:-1])
        # Saving to file if required
        if filepath!='':
            with open(filepath, 'w') as file:
                file.write(s)
        return s

class Line:
    def __neg__(self):
        return self.__class__(self.points[::-1])

    def DirectionVector(self, unit=False):
        u = self.points[1] - self.points[0]
        if unit:
            u.Normalize()
        return u

    def NormalVector(self, unit=False):
        return self.DirectionVector(unit).NormalVector()

class Line2D(Primitive2D, Line):
    """
    Define an infinte line given by two points.
    """
    def __init__(self, point1, point2, name=''):
        Primitive2D.__init__(self, name=name)
        self.points=[point1, point2]

    def To3D(self, plane_origin, x1, x2):
        p3D = [p.To3D(plane_origin, x1, x2) for p in self.points]
        return Line2D(*p3D, self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Line2D(*[p.Rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Line2D(*[p.Translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset, copy=False)

    def PointDistance(self, point):
        """
        Computes the distance of a point to line
        """
        p1, p2 = self.points
        u = p2 - p1
        t = (point-p1).Dot(u) / u.Norm()**2
        projection = p1 + t * u # Projection falls on the segment
        return (point-projection).Norm()

    def PointProjection(self, point, curvilinear_abscissa=False):
        p1, p2 = self.points
        t = (point - p1).Dot(p2 - p1) / (p2-p1).Norm()**2
        projection = p1 + t * (p2-p1)
        if curvilinear_abscissa:
            return projection,t
        return projection

    def MPLPlot(self, ax=None, color='k', linestyle = '-.'):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = ax.figure

        p1, p2 = self.points
        u = p2 - p1
        plt.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color)
        p3 = p1 - 3*u
        p4 = p2 + 4*u
        ax.plot([p3[0], p4[0]], [p3[1], p4[1]], color=color, linestyle = linestyle)
        return fig ,ax

    def CreateTangentCircle(self, point, other_line):
        """
        Computes the two circles that are tangent to 2 lines and intersect
        a point located on one of the two lines.
        """

        # point will be called I(x_I, y_I)
        # self will be (AB)
        # line will be (CD)

        if math.isclose(self.PointDistance(point), 0, abs_tol=1e-10):
            I = Vector2D((point[0], point[1]))
            A = Vector2D((self.points[0][0], self.points[0][1]))
            B = Vector2D((self.points[1][0], self.points[1][1]))
            C = Vector2D((other_line.points[0][0], other_line.points[0][1]))
            D = Vector2D((other_line.points[1][0], other_line.points[1][1]))

        elif math.isclose(other_line.PointDistance(point), 0, abs_tol=1e-10):
            I = Vector2D((point[0], point[1]))
            C = Vector2D((self.points[0][0], self.points[0][1]))
            D = Vector2D((self.points[1][0], self.points[1][1]))
            A = Vector2D((other_line.points[0][0], other_line.points[0][1]))
            B = Vector2D((other_line.points[1][0], other_line.points[1][1]))
        else:
            raise AttributeError("The point isn't on any of the two lines")

        # CHANGEMENT DE REPAIRE
        new_u = Vector2D((B-A))
        new_u.Normalize()
        new_v = new_u.NormalVector(unit=True)
        new_basis = Frame2D(I, new_u, new_v)

        new_A = new_basis.NewCoordinates(A)
        new_B = new_basis.NewCoordinates(B)
        new_C = new_basis.NewCoordinates(C)
        new_D = new_basis.NewCoordinates(D)

# =============================================================================
# LES SEGMENTS DECRIVENT UNE SEULE ET MEME DROITE
#   => AUCUNE SOLUTION
# =============================================================================
        if new_C[1] == 0 and new_D[1] == 0:

            return None, None

# =============================================================================
# LES SEGMENTS SONT PARALLELES
#   => 1 SOLUTION
# =============================================================================
        elif math.isclose(self.DirectionVector(unit=True).Dot(other_line.NormalVector(unit=True)), 0, abs_tol=1e-06):

            segments_distance = abs(new_C[1] - new_A[1])
            r = segments_distance / 2
            new_circle_center = Point2D((0, npy.sign(new_C[1] - new_A[1])*r))
            circle_center = new_basis.OldCoordinates(new_circle_center)
            circle = Circle2D(circle_center, r)

            return circle, None

# =============================================================================
# LES SEGMENTS SONT PERPENDICULAIRES
#   => 2 SOLUTIONS
# =============================================================================
        elif math.isclose(self.DirectionVector(unit=True).Dot(other_line.DirectionVector(unit=True)), 0, abs_tol=1e-06):

            line_AB = Line2D(Point2D(new_A), Point2D(new_B))
            line_CD = Line2D(Point2D(new_C), Point2D(new_D))
            new_pt_K = Point2D.LinesIntersection(line_AB ,line_CD)

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
            new_pt_K = Point2D.LinesIntersection(line_AB ,line_CD)
            pt_K = Point2D(new_basis.OldCoordinates(new_pt_K))

            if pt_K == I:
                return None, None

            # CHANGEMENT DE REPERE:
            new_u2 = Vector2D(pt_K-I)
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
                if math.isclose(teta1, math.pi, abs_tol=1e-08) or math.isclose(teta1, 0., abs_tol=1e-08):
                    teta = teta2
                elif math.isclose(teta2, math.pi, abs_tol=1e-08) or math.isclose(teta2, 0., abs_tol=1e-08):
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

class LineSegment2D(Line2D):
    """
    Define a line segment limited by two points
    """
    def __init__(self,point1, point2, name=''):
        Line2D.__init__(self, point1, point2, name = name)

    def _get_geo_points(self):
        return self.points

    geo_points=property(_get_geo_points)

    def Length(self):
        return self.points[1].PointDistance(self.points[0])

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.points[0] + self.DirectionVector(unit=True) * curvilinear_abscissa

    def PointDistance(self, point, return_other_point=False):
        """
        Computes the distance of a point to segment of line
        """
        if self.points[0] == self.points[1]:
            return 0, Point2D(point)
        distance, point = LineSegment2DPointDistance([p.vector for p in self.points], point.vector)
        if return_other_point:
            return distance, Point2D(point)
        return distance

    def PointProjection(self, point, curvilinear_abscissa=False):
        point, curv_abs = Line2D.PointProjection(self, point, True)
        if curv_abs <= 0.:
            point = self.points[0]
            curv_abs = 0.
        elif curv_abs >= 1.:
            point = self.points[1]
            curv_abs = 1.

        if curvilinear_abscissa:
            return point, curv_abs
        else:
            return point

    def MPLPlot(self, ax=None, color='k', arrow=False, width=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = None

        p1, p2 = self.points
        if arrow:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color)
            length = ((p1[0] - p2[0])**2 + (p1[1] - p2[1])**2)**0.5
            if width is None:
                width = length / 1000.
                head_length = length/20.
                head_width = head_length/2.
            else:
                head_width = 2*width
                head_length = head_width
            ax.arrow(p1[0], p1[1], (p2[0] - p1[0])/length*(length - head_length),
                     (p2[1] - p1[1])/length*(length - head_length),
                     head_width = head_width, fc = 'b', linewidth = 0,
                     head_length = head_length, width = width, alpha = 0.3)
        else:
            ax.plot([p1[0], p2[0]], [p1[1], p2[1]], color=color)
        return fig, ax

    def To3D(self, plane_origin, x1, x2):
        p3D=[p.To3D(plane_origin,x1,x2) for p in self.points]
        return LineSegment3D(*p3D,self.name)

    def to_line(self):
        return Line2D(*self.points)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return LineSegment2D(*[p.Rotation(center,angle,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center,angle,copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return LineSegment2D(*[p.Translation(offset,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset,copy=False)

    def GeoScript(self, primitive_index, points_indices):
        s='Line({}) = {{{}, {}}};\n'.format(primitive_index,*points_indices)
        return s,primitive_index+1

    def plot_data(self, marker=None, color='black', stroke_width=1,
                 dash=False, opacity=1, arrow=False):
        return {'type' : 'line',
                'data' : [self.points[0].vector[0], self.points[0].vector[1],
                          self.points[1].vector[0], self.points[1].vector[1]],
                'color' : color,
                'marker' : marker,
                'size' : stroke_width,
                'dash' : dash,
                'opacity' : opacity,
                'arrow': arrow
                }

    def CreateTangentCircle(self, point, other_line):
        circle1, circle2 = Line2D.CreateTangentCircle(other_line, point, self)
        if circle1 is not None:
            point_J1, curv_abs1 = Line2D.PointProjection(self, circle1.center, True)
            if curv_abs1 < 0. or curv_abs1 > 1.:
                circle1 = None
        if circle2 is not None:
            point_J2, curv_abs2 = Line2D.PointProjection(self, circle2.center, True)
            if curv_abs2 < 0. or curv_abs2 > 1.:
                circle2 = None
        return circle1, circle2
    

class Arc2D(Primitive2D):
    def __init__(self, start, interior, end, name=''):
        Primitive2D.__init__(self, name)
        self.interior = interior
        self.start = start
        self.end = end
        xi, yi = interior.vector
        xe, ye = end.vector
        xs, ys = start.vector
        A = npy.array([[2*(xs-xi), 2*(ys-yi)],
                       [2*(xs-xe), 2*(ys-ye)]])
        b = - npy.array([xi**2 + yi**2 - xs**2 - ys**2,
                         xe**2 + ye**2 - xs**2 - ys**2])
        self.center = Point2D(solve(A,b))
        r1 = self.start - self.center
        r2 = self.end - self.center
        ri = self.interior - self.center

        self.radius = r1.Norm()
        angle1 = npy.arctan2(r1.vector[1], r1.vector[0])
        angle2 = npy.arctan2(r2.vector[1], r2.vector[0])

        anglei = npy.arctan2(ri.vector[1], ri.vector[0])
        order = [y for x, y in sorted(zip([angle1, anglei, angle2], [0, 1, 2]))]
        order = order*2
        i = order.index(0)
        if order[i+1] == 1:
            # Trigo wise angle should be plus
            self.is_trigo = True
            self.angle1 = angle1
            self.angle2 = angle2
            if angle1 > angle2:
                self.angle = angle2 - angle1 + 2 * math.pi
            else:
                self.angle = angle2 - angle1
#            self.angle = abs(self.angle2 - self.angle1)
        else:
            # Clock wise
            self.is_trigo = False
            self.angle1 = angle2
            self.angle2 = angle1
#            self.angle = -abs(self.angle2 - self.angle1)
            if angle1 < angle2:
                self.angle = angle2 - angle1 + 2 * math.pi
            else:
                self.angle = angle2 - angle1
                
    def _get_points(self):
        return [self.start,self.interior,self.end]

    points=property(_get_points)


    def _get_geo_points(self):
        return [self.start,self.interior,self.end]

    geo_points=property(_get_geo_points)
    
    def tessellation_points(self, resolution=20):
#        return [self.center + self.radius*math.cos(teta)*Vector2D((1,0)) + self.radius*math.sin(teta)*Vector2D((0,1)) \
#                for teta in npy.linspace(0, self.angle, resolution)]
        points = []
        if not self.is_trigo:
            delta_angle = -abs(self.angle1-self.angle2)/(resolution-1)
        else:
            delta_angle = abs(self.angle2-self.angle1)/(resolution-1)
        points.append(self.start)
        for i in range(resolution-2):
            point_to_add = points[-1].Rotation(self.center, delta_angle)
            points.append(point_to_add)
        points.append(self.end)
        return points
        

    def Length(self):
        return self.radius * abs(self.angle)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        if self.angle>0:
            return self.start.Rotation(self.center, curvilinear_abscissa/self.radius)
        else:
            return self.start.Rotation(self.center, -curvilinear_abscissa/self.radius)

    def MiddlePoint(self):
        l = self.Length()
        return self.PointAtCurvilinearAbscissa(0.5*l)

    def GeoScript(self, primitive_index, points_indices):
        s='Circle({}) = {{{}, {}, {}}};\n'.format(primitive_index,*points_indices)
        return s,primitive_index+1

    def Area(self):
        if self.angle2<self.angle1:
            angle=self.angle2+2*math.pi-self.angle1
        else:
            angle=self.angle2-self.angle1
        return self.radius**2*angle/2

    def CenterOfMass(self):
#        u=self.middle.vector-self.center.vector
        u = self.MiddlePoint() - self.center
        u.Normalize()
        alpha = abs(self.angle)
        return self.center + 4/(3*alpha)*self.radius*math.sin(alpha*0.5)*u

    def MPLPlot(self, ax, color='k'):
        pc = self.center.vector
#        ax.plot([pc[0]], [pc[1]], 'or')
#        ax.plot([self.interior[0]], [self.interior[1]], 'ob')
        ax.add_patch(Arc(pc, 2*self.radius, 2*self.radius, angle=0,
                    theta1=self.angle1*0.5/math.pi*360,
                    theta2=self.angle2*0.5/math.pi*360,
                    color=color))

    def To3D(self,plane_origin, x, y):
        ps = self.start.To3D(plane_origin, x, y)
        pi = self.interior.To3D(plane_origin, x, y)
        pe = self.end.To3D(plane_origin, x, y)

        return Arc3D(ps, pi, pe, self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Arc2D(*[p.Rotation(center,angle,copy=True) for p in [self.start,self.interior,self.end]])
        else:
            self.__init__(*[p.Rotation(center,angle,copy=True) for p in [self.start,self.interior,self.end]])

    def Translation(self, offset, copy=True):
        if copy:
            return Arc2D(*[p.Translation(offset,copy=True) for p in [self.start,self.interior,self.end]])
        else:
            self.__init__(*[p.Translation(offset,copy=True) for p in [self.start,self.interior,self.end]])

    def SecondMomentArea(self, point):
        """
        Second moment area of part of disk
        """
        if self.angle2<self.angle1:
            angle2 = self.angle2+2*math.pi

        else:
            angle2 = self.angle2
        angle1 = self.angle1

        Ix = self.radius**4/8*(angle2-angle1+0.5*(math.sin(2*angle1)-math.sin(2*angle2)))
        Iy = self.radius**4/8*(angle2-angle1+0.5*(math.sin(2*angle2)-math.sin(2*angle1)))
        Ixy = self.radius**4/8*(math.cos(angle1)**2-math.cos(angle2)**2)
        Ic = npy.array([[Ix, Ixy], [Ixy, Iy]])
        return geometry.Huygens2D(Ic, self.Area(), self.center, point)

    def Discretise(self, num=10):
        list_node = []
        if (self.angle1 < 0) and (self.angle2 > 0):
            delta_angle = -self.angle1 + self.angle2
        elif (self.angle1 > 0) and (self.angle2 < 0):
            delta_angle =  (2*npy.pi + self.angle2) - self.angle1
        else:
            delta_angle = self.angle2 - self.angle1
        for angle in npy.arange(self.angle1, self.angle1 + delta_angle, delta_angle/(num*1.)):
            list_node.append(Point2D(self.center + self.radius*Vector2D((npy.cos(angle), npy.sin(angle)))))
        list_node.append(Point2D(self.center + self.radius*Vector2D((npy.cos(self.angle1 + delta_angle), npy.sin(self.angle1 + delta_angle)))))
        if list_node[0] == self.start:
            return list_node
        else:
            return list_node[::-1]

    def plot_data(self, marker=None, color='black', stroke_width=1, opacity=1):
        list_node = self.Discretise()
        data = []
        for nd in list_node:
            data.append({'x': nd.vector[0], 'y': nd.vector[1]})
        return {'type' : 'arc',
                    'cx' : self.center.vector[0],
                    'cy' : self.center.vector[1],
                    'data' : data,
                    'r' : self.radius,
                    'color' : color,
                    'opacity' : opacity,
                    'size' : stroke_width,
                    'dash' : None,
                    'marker' : marker,
                    'angle1' : self.angle1,
                    'angle2' : self.angle2, }

class Circle2D(Primitive2D):
    def __init__(self,center,radius,name=''):
        Primitive2D.__init__(self,name)
        self.center = center
        self.radius = radius
        self.utd_geo_points = False

    def _get_geo_points(self):
        if not self.utd_geo_points:
            self._geo_start = self.center+self.radius*Point2D((1,0))
            self.utd_geo_points = True
        return [self._geo_start, self.center, self._geo_start]

    geo_points = property(_get_geo_points)
    
    def tessellation_points(self, resolution=40):
        return [self.center + self.radius*math.cos(teta)*Vector2D((1,0)) + self.radius*math.sin(teta)*Vector2D((0,1)) \
                for teta in npy.linspace(0, 2*math.pi, resolution+1)]

    def Length(self):
        return 2* math.pi * self.radius

    def GeoScript(self, primitive_index, points_indices):
        s = 'Circle({}) = {{{}, {}, {}}};\n'.format(primitive_index,*points_indices)
        return s, primitive_index+1

    def MPLPlot(self, ax, linestyle='-', color='k', linewidth=1):
        pc = self.center.vector
        if self.radius > 0:
            ax.add_patch(Arc(pc,
                             2*self.radius,
                             2*self.radius,
                             angle=0,
                             theta1=0,
                             theta2=360,
                             color=color,
                             linestyle=linestyle,
                             linewidth=linewidth))

    def To3D(self, plane_origin, x, y):
        normal = Vector3D(npy.cross(x.vector, y.vector))
        pc=self.center.To3D(plane_origin, x, y)
        return Circle3D(pc,self.radius,normal, self.name)

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Circle2D(self.center.Rotation(center,angle,copy=True),self.radius)
        else:
            self.center.Rotation(center,angle,copy=False)
            self.utd_geo_points=False

    def Translation(self,offset,copy=True):
        if copy:
            return Circle2D(self.center.Translation(offset,copy=True),self.radius)
        else:
            self.center.Translation(offset,copy=False)
            self.utd_geo_points=False

    def Area(self):
        return math.pi*self.radius**2

    def SecondMomentArea(self, point):
        """
        Second moment area of part of disk
        """
        I = math.pi*self.radius**4/4
        Ic = npy.array([[I,0],[0,I]])
        return geometry.Huygens2D(Ic,self.Area(),self.center,point)

    def CenterOfMass(self):
        return self.center

    def plot_data(self, marker=None, color='black', stroke_width=1, opacity=1, fill=None):
        return {'type' : 'circle',
                'cx' : self.center.vector[0],
                'cy' : self.center.vector[1],
                'r' : self.radius,
                'color' : color,
                'opacity' : opacity,
                'size' : stroke_width,
                'dash' : None,
                'fill' : fill}

class Polygon2D(CompositePrimitive2D):
    # TODO: inherit from contour?
    def __init__(self,points, name=''):
        self.points = points
        primitives = []
        for p1,p2 in zip(points,points[1:]+[points[0]]):
            primitives.append(LineSegment2D(p1,p2))

        self.line_segments = self._LineSegments()

        CompositePrimitive2D.__init__(self, primitives, name)


    def Area(self):

        x=[point.vector[0]for point in self.points]
        y=[point.vector[1]for point in self.points]

        return 0.5*npy.abs(npy.dot(x,npy.roll(y,1))-npy.dot(y,npy.roll(x,1)))

    def CenterOfMass(self):

        x = [point.vector[0] for point in self.points]
        y = [point.vector[1] for point in self.points]


        xi_xi1 = x+npy.roll(x,-1)
        yi_yi1 = y+npy.roll(y,-1)
        xi_yi1 = npy.multiply(x,npy.roll(y,-1))
        xi1_yi = npy.multiply(npy.roll(x,-1),y)

        a = 0.5*npy.sum(xi_yi1-xi1_yi)# signed area!
#        a=self.Area()
        if not math.isclose(a, 0, abs_tol=1e-08):
            cx = npy.sum(npy.multiply(xi_xi1,(xi_yi1-xi1_yi)))/6./a
            cy = npy.sum(npy.multiply(yi_yi1,(xi_yi1-xi1_yi)))/6./a
            return Point2D((cx, cy))

        else:
            raise NotImplementedError

    def PointBelongs(self, point):
        """
        Ray casting algorithm copied from internet...
        """
        return PolygonPointBelongs(point.vector, [p.vector for p in self.points])

    def SecondMomentArea(self, point):
        Ix, Iy, Ixy = 0, 0, 0
        for pi, pj in zip(self.points,self.points[1:]+[self.points[0]]):
            xi, yi = pi.vector-point.vector
            xj, yj = pj.vector-point.vector
            Ix += (yi**2 + yi*yj + yj**2)*(xi*yj - xj*yi)
            Iy += (xi**2 + xi*xj + xj**2)*(xi*yj - xj*yi)
            Ixy += (xi*yj + 2*xi*yi + 2*xj*yj + xj*yi)*(xi*yj - xj*yi)
        if Ix < 0:
            Ix =- Ix
            Iy =- Iy
            Ixy =- Ixy
        return npy.array([[Ix/12., Ixy/24.], [Ixy/24., Iy/12.]])

    def _LineSegments(self):
        lines=[]
        for p1,p2 in zip(self.points,self.points[1:]+[self.points[0]]):
            lines.append(LineSegment2D(p1,p2))
        return lines

    def Rotation(self, center, angle, copy=True):
        if copy:
            return Polygon2D([p.Rotation(center,angle,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center,angle,copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Polygon2D([p.Translation(offset,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset,copy=False)

    def PointBorderDistance(self, point, return_other_point=False):
        """
        Compute the distance to the border distance of polygon
        Output is always positive, even if the point belongs to the polygon
        """
        d_min, other_point_min = self.line_segments[0].PointDistance(point, return_other_point=True)
        for line in self.line_segments[1:]:
            d, other_point = line.PointDistance(point, return_other_point=True)
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
        sorted_index = sorted(range(len(self.points)), key=lambda p: (self.points[p][0], self.points[p][1]))
        nb = len(sorted_index)
        segments = []
        deleted = []

        while len(sorted_index) != 0: # While all the points haven't been swept
            # Stock the segments between 2 consecutive edges
            # Ex: for the ABCDE polygon, if Sweep Line is on C, the segments
            #   will be (C,B) and (C,D)
            if sorted_index[0]-1 < 0:
                segments.append((sorted_index[0], nb-1))
            else:
                segments.append((sorted_index[0], sorted_index[0]-1))
            if sorted_index[0] >= len(self.points)-1:
                segments.append((sorted_index[0], 0))
            else:
                segments.append((sorted_index[0], sorted_index[0]+1))


            # Once two edges linked by a segment have been swept, delete the
            # segment from the list
            to_del = []
            for index in deleted:
                if abs(index-sorted_index[0]) == 1 or abs(index-sorted_index[0]) == nb-1:
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
                    if segment1[0] != segment2[0] and segment1[1] != segment2[1] and segment1[0] != segment2[1] and segment1[1] != segment2[0]:

                        line1 = LineSegment2D(Point2D(self.points[segment1[0]]), Point2D(self.points[segment1[1]]))
                        line2 = LineSegment2D(Point2D(self.points[segment2[0]]), Point2D(self.points[segment2[1]]))

                        p, a, b = Point2D.LinesIntersection(line1, line2, True)

                        if p is not None:
                            if a >= 0+epsilon and a <= 1-epsilon and b >= 0+epsilon and b <= 1-epsilon:
                                return True, line1, line2

        return False, None, None


    def Dict(self):
        d = {'points': [point.Dict() for point in self.points], 'name':self.name}
        return d

    @classmethod
    def DictToObject(cls, dict_):
        return cls([Point2D.DictToObject(p) for p in dict_['points']], name=dict_['name'])

    def plot_data(self, marker=None, color='black', stroke_width=1, opacity=1):
        data = []
        for nd in self.points:
            data.append({'x': nd.vector[0], 'y': nd.vector[1]})
        return {'type' : 'path',
                    'data' : data,
                    'color' : color,
                    'size' : stroke_width,
                    'dash' : None,
                    'marker' : marker,
                    'opacity' : opacity}
        
    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
            
        ax.plot([p[0] for p in self.points]+[self.points[0][0]], [p[1] for p in self.points]+[self.points[0][1]], '-')
        return ax

class Primitive3D:
    def __init__(self, basis_primitives=None, name=''):
        self.name = name
        self.basis_primitives = basis_primitives # une liste
        if basis_primitives is None:
            self.basis_primitives = []

class Vector3D(Vector):
    _standalone_in_db = False

    _jsonschema = {
        "definitions": {},
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "title": "powerpack.mechanical.Vector3D Base Schema",
        "required": ["vector"],
        "properties": {
            'vector' : {
                "type" : "array",
                "order" : 0,
                "items" : {
                    "type" : "number",
                    "step" : 1,
                    "minimum" : -1,
                    "maximum" : 1
                    },
                "minItems": 3,
                "maxItems": 3,
                "examples": [[1, 0, 0]],
                "editable" : True,
                "description" : "Vector array"
                }
            }
        }    

    def __init__(self, vector, name=''):
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

    def __round__(self, ndigits=6):
        return self.__class__((round(self.vector[0], ndigits),
                               round(self.vector[1], ndigits),
                               round(self.vector[2], ndigits)))
        
    def __hash__(self):
        return int(1000*(self.vector[0]+self.vector[1]+self.vector[2]))

    def __eq__(self, other_vector):
        return math.isclose(self.vector[0], other_vector.vector[0], abs_tol=1e-08) \
        and math.isclose(self.vector[1], other_vector.vector[1], abs_tol=1e-08) \
        and math.isclose(self.vector[2], other_vector.vector[2], abs_tol=1e-08) 
    
    def Dot(self, other_vector):
        return Vector3DDot(self.vector, other_vector.vector)

    def Cross(self, other_vector):
        return self.__class__(vector3D_cross(self.vector, other_vector.vector))

    def Norm(self):
        return Vector3DNorm(self.vector)


    def Normalize(self):
        """
        Normalize the vector modifying it's coordinate
        """
        n = self.Norm()
        if n == 0:
            raise ZeroDivisionError

        self.vector[0] /= n
        self.vector[1] /= n
        self.vector[2] /= n


    def Rotation(self, center, axis, angle, copy=True):
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
            return Point3D(vector2)
        else:
            self.vector = list(vector2)
            
    def x_rotation(self, angle, copy=True):
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
            
    def y_rotation(self, angle, copy=True):
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
            
    def z_rotation(self, angle, copy=True):
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
                self.vector = new_vector.vector
                
        if side == 'new':
            new_vector = frame.NewCoordinates(self)
            if copy:
                return new_vector
            else:
                self.vector = new_vector.vector

    def To2D(self, plane_origin, x, y):
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

    def Copy(self):
        return Vector3D(self.vector)

    @classmethod
    def DictToObject(cls, dict_):
        return cls(dict_['vector'])

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
        return ax


x3D = Vector3D((1, 0, 0))
y3D = Vector3D((0, 1, 0))
z3D = Vector3D((0, 0, 1))

X3D = Vector3D((1, 0, 0))
Y3D = Vector3D((0, 1, 0))
Z3D = Vector3D((0, 0, 1))


class Point3D(Vector3D):
    _standalone_in_db = False
    _jsonschema = {
        "definitions": {},
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "title": "powerpack.mechanical.Point3D Base Schema",
        "required": ["vector"],
        "properties": {
            'vector' : {
                "type" : "object",
                "order" : 0,
                "classes" : ["volmdlr.core.Vector3D"],
                "editable" : True,
                "description" : "Vector array"
                }
            }
        }

    def __init__(self, vector, name=''):
        Vector3D.__init__(self, vector)
        self.name=name

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

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        ax.scatter(*self.vector)
        return ax

    def PlaneProjection3D(self, plane_origin, x, y):
        z = x.Cross(y)
        z /= z.Norm()
        return self - z.Dot(self-plane_origin)*z

    def PlaneProjection2D(self, x, y):
        z = npy.cross(x.vector,y.vector)
        z = x.Cross(y)
        z.Normalize()
        p3d = self - self.Dot(z)*z
        u1 = p3d.Dot(x)
        u2 = p3d.Dot(y)
        return Point2D((u1, u2))


    def To2D(self, plane_origin, x, y):
        x2d = self.Dot(x) - plane_origin.Dot(x)
        y2d = self.Dot(y) - plane_origin.Dot(y)
        return Point2D((x2d,y2d))

    def PointDistance(self, point2):
        return (self-point2).Norm()

    @classmethod
    def from_step(cls, arguments, object_dict):
        return cls([float(i)/1000 for i in arguments[1][1:-1].split(",")],
                    arguments[0][1:-1])
        
    def Babylon(self):
        s = 'var sphere = BABYLON.MeshBuilder.CreateSphere("point", {diameter: 0.05}, scene);\n'
        s += "sphere.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n".format(self.vector[0],self.vector[1],self.vector[2])
        s += 'var mat = new BABYLON.StandardMaterial("mat", scene);\n'
        s += 'mat.diffuseColor = new BABYLON.Color3(1, 0, 0);\n'
        s += 'sphere.material = mat;\n'
        return s
        

o3D = Point3D((0, 0, 0))
O3D = Point3D((0, 0, 0))


class Plane3D:
    def __init__(self, origin, vector1, vector2, name=''):
        self.origin = Point3D(origin.vector)
        vector1 = Vector3D(vector1.vector)
        vector1.Normalize()
        vector2 = Vector3D(vector2.vector)
        vector2.Normalize()
        self.vectors = [vector1, vector2]
        self.name = name
        self.normal = self.vectors[0].Cross(self.vectors[1])
        self.normal.Normalize()
                
    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        origin = frame3d.origin
        vector1 = frame3d.u
        vector2 = frame3d.v
        # TRANSFORMER EN 3D TOUS LES OBJETS LIES AU PLAN
        return cls(origin, vector1, vector2, arguments[0][1:-1])
    
    @classmethod
    def from_3_points(cls, point1, point2, point3):
        vector1 = point2 - point1
        vector2 = point3 - point1
        vector1.Normalize()
        vector2.Normalize()
        normal = vector1.Cross(vector2)
        normal.Normalize()
        vector = normal.Cross(vector1)
        return cls(point1.copy(), vector1.copy(), vector.copy())
    
    @classmethod
    def from_normal(cls, point, normal):
        v1 = normal.RandomUnitNormalVector()
        v2 = v1.Cross(normal)
        return cls(point, v1+point, v2+point)
    
    @classmethod
    def from_points(cls, points):
        if len(points) < 3:
            raise ValueError
        elif len(points) == 3:
            return cls.from_3_points(Point3D(points[0].vector), Vector3D(points[1].vector), Vector3D(points[2].vector))
        else:
            points = [p.copy() for p in points]
            indexes_to_del = []
            for i, point in enumerate(points[1:]):
                if point == points[0]:
                    indexes_to_del.append(i)
            for index in indexes_to_del[::-1]:
                del points[index+1]

            origin = Point3D(points[0].vector)
            vector1 = Vector3D(points[1]-origin)
            vector1.Normalize()
            vector2_min = Vector3D(points[2]-origin)
            vector2_min.Normalize()
            dot_min = abs(vector1.Dot(vector2_min))
            for point in points[3:]:
                vector2 = Vector3D(point-origin)
                vector2.Normalize()
                dot = abs(vector1.Dot(vector2))
                if dot < dot_min:
                    vector2_min = vector2
                    dot_min = dot
            return cls.from_3_points(origin, vector1+origin, vector2_min+origin)
        
    def point_on_plane(self, point):
        if math.isclose(self.normal.Dot(point-self.origin), 0, abs_tol=1e-6):
            return True
        return False

    def line_intersection(self, line):
        u = line.points[1] - line.points[0]
        w = line.points[0] - self.origin
        if math.isclose(self.normal.Dot(u), 0, abs_tol=1e-08):
            return None
        intersection_abscissea = - self.normal.Dot(w) / self.normal.Dot(u)
        return line.points[0] + intersection_abscissea * u

    def linesegment_intersection(self, linesegment, abscissea=False):
        u = linesegment.points[1] - linesegment.points[0]
        w = linesegment.points[0] - self.origin
        normalDotu = self.normal.Dot(u)
        if math.isclose(normalDotu, 0, abs_tol=1e-08):
            if abscissea:
                return None, None
            return None
        intersection_abscissea = - self.normal.Dot(w) / normalDotu
        if intersection_abscissea < 0 or intersection_abscissea > 1:
            if abscissea:
                return None, None
            return None
        if abscissea:
            return linesegment.points[0] + intersection_abscissea * u, intersection_abscissea
        return linesegment.points[0] + intersection_abscissea * u

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_origin = self.origin.Rotation(center, axis, angle, True)
            new_vector1 = self.vectors[0].Rotation(center, axis, angle, True)
            new_vector2 = self.vectors[1].Rotation(center, axis, angle, True)
            return Plane3D(new_origin, new_vector1, new_vector2, self.name)
        else:
            self.origin.Rotation(center, axis, angle, True)
            self.vectors[0].Rotation(center, axis, angle, True)
            self.vectors[1].Rotation(center, axis, angle, True)

    def Translation(self, offset, copy=True):
        if copy:
            new_origin = self.origin.Translation(offset, True)
            return Plane3D(new_origin, self.vectors[0], self.vectors[1], self.name)
        else:
            self.origin.Translation(offset, False)
    
    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            new_origin = frame.OldCoordinates(self.origin)
            new_vector1 = frame.Basis().OldCoordinates(self.vectors[0])
            new_vector2 = frame.Basis().OldCoordinates(self.vectors[1])
            if copy:
                return Plane3D(new_origin, new_vector1, new_vector2, self.name)
            else:
                self.origin = new_origin
                self.vectors = [new_vector1, new_vector2]
                self.normal = frame.Basis().OldCoordinates(self.normal)
                self.normal.Normalize()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.origin)
            new_vector1 = frame.Basis().NewCoordinates(self.vectors[0])
            new_vector2 = frame.Basis().NewCoordinates(self.vectors[1])
            if copy:
                return Plane3D(new_origin, new_vector1, new_vector2, self.name)
            else:
                self.origin = new_origin
                self.vectors = [new_vector1, new_vector2]
                self.normal = frame.Basis().NewCoordinates(self.normal)
                self.normal.Normalize()
                
    def copy(self):
        new_origin = self.origin.Copy()
        new_vector1 = self.vectors[0].Copy()
        new_vector2 = self.vectors[1].Copy()
        return Plane3D(new_origin, new_vector1, new_vector2, self.name)
        
    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        self.origin.MPLPlot(ax)
        self.vectors[0].MPLPlot(ax, starting_point=self.origin, color='r')
        self.vectors[1].MPLPlot(ax, starting_point=self.origin, color='g')
        return ax
    
    def Babylon(self):
        s = 'var myPlane = BABYLON.MeshBuilder.CreatePlane("myPlane", {width: 0.5, height: 0.5, sideOrientation: BABYLON.Mesh.DOUBLESIDE}, scene);\n'
        s += 'myPlane.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n'.format(self.origin[0], self.origin[1], self.origin[2])

        s += 'var axis1 = new BABYLON.Vector3({}, {}, {});\n'.format(self.vectors[0][0], self.vectors[0][1], self.vectors[0][2])
        s += 'var axis2 = new BABYLON.Vector3({}, {}, {});\n'.format(self.vectors[1][0], self.vectors[1][1], self.vectors[1][2])
        s += 'var axis3 = new BABYLON.Vector3({}, {}, {});\n'.format(self.normal[0], self.normal[1], self.normal[2])
        s += 'var orientation = BABYLON.Vector3.RotationFromAxis(axis1, axis2, axis3);\n'
        s += 'myPlane.rotation = orientation;\n'

        s += 'var planemat = new BABYLON.StandardMaterial("planemat", scene);\n'
        s += 'planemat.alpha = 0.4;\n'
        s += 'myPlane.material = planemat;\n'
        
#        print(self.__dict__)
        
        return s


class Basis3D(Basis):
    """
    Defines a 3D basis
    :param u: first vector of the basis
    :param v: second vector of the basis
    :param w: third vector of the basis
    """
    _standalone_in_db = False
    _jsonschema = {
        "definitions": {},
        "$schema": "http://json-schema.org/draft-07/schema#",
        "type": "object",
        "title": "powerpack.mechanical.Basis3D Base Schema",
        "required": ['vectors'],
        "properties": {
            'vectors' : {
                'type' : 'array',
                'items' : {
                    'type' : 'object',
                    "editable" : True,
                    'classes' : ['volmdlr.core.Vector3D']
                    },
                'order' : 0,
                'editable' : True
                }
            }
        }

    # TODO: create a Basis and Frame class to mutualize between 2D and 2D
    def __init__(self, u, v, w, name=''):
        self.u = u
        self.v = v
        self.w = w
        self.name = name
        
        
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

    def __round__(self, ndigits=6):
        return self.__class__((round(self.u, ndigits),
                               round(self.v, ndigits),
                               round(self.w, ndigits)))

    def __repr__(self):
        return '{}: U={}, V={}, W={}'.format(self.__class__.__name__, *self.vectors)
    
    def _get_vectors(self):
        return (self.u, self.v, self.w)

    vectors = property(_get_vectors)

    @classmethod
    def from_two_vectors(cls, vector1, vector2):
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

    def Rotation(self, axis, angle, copy=True):
        center = o3D
        new_u = self.u.Rotation(center, axis, angle, True)
        new_v = self.v.Rotation(center, axis, angle, True)
        new_w = self.w.Rotation(center, axis, angle, True)

        if copy:
            return Basis3D(new_u, new_v, new_w, self.name)
        else:
            self.u = new_u
            self.v = new_v
            self.w = new_w

    def x_rotation(self, angle, copy=True):
        new_u = self.u.x_rotation(angle, True)
        new_v = self.v.x_rotation(angle, True)
        new_w = self.w.x_rotation(angle, True)

        if copy:
            return Basis3D(new_u, new_v, new_w, self.name)
        else:
            self.u = new_u
            self.v = new_v
            self.w = new_w
        
    def y_rotation(self, angle, copy=True):
        new_u = self.u.y_rotation(angle, True)
        new_v = self.v.y_rotation(angle, True)
        new_w = self.w.y_rotation(angle, True)

        if copy:
            return Basis3D(new_u, new_v, new_w, self.name)
        else:
            self.u = new_u
            self.v = new_v
            self.w = new_w
        
    def z_rotation(self, angle, copy=True):
        new_u = self.u.z_rotation(angle, True)
        new_v = self.v.z_rotation(angle, True)
        new_w = self.w.z_rotation(angle, True)

        if copy:
            return Basis3D(new_u, new_v, new_w, self.name)
        else:
            self.u = new_u
            self.v = new_v
            self.w = new_w
        
    def EulerRotation(self, angles, copy=True):
        psi, theta, phi = angles
        center = o3D

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
        return Matrix33(self.u[0], self.v[0], self.w[0],
                        self.u[1], self.v[1], self.w[1],
                        self.u[2], self.v[2], self.w[2])
#        return npy.array([[self.u[0], self.v[0], self.w[0]],
#                          [self.u[1], self.v[1], self.w[1]],
#                          [self.u[2], self.v[2], self.w[2]]])

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

    @classmethod
    def DictToObject(cls, dict_):
        vectors = [Vector3D.DictToObject(vector_dict) for vector_dict in dict_['vectors']]
        return cls(*vectors)


xyz = Basis3D(x3D, y3D, z3D)
XYZ = Basis3D(x3D, y3D, z3D)
YZX = Basis3D(y3D, z3D, x3D)
ZXY = Basis3D(z3D, x3D, y3D)

class Frame3D(Basis3D):
    """
    Defines a 3D frame
    :param origin: origin of the basis
    :param u: first vector of the basis
    :param v: second vector of the basis
    :param w: third vector of the basis
    """
    def __init__(self, origin, u, v, w, name=''):
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
        

    def Basis(self):
        return Basis3D(self.u, self.v, self.w)

    def NewCoordinates(self, vector):
        return Basis3D.NewCoordinates(self, vector - self.origin)

    def OldCoordinates(self, vector):
        return Basis3D.OldCoordinates(self, vector) + self.origin

    def Rotation(self, axis, angle, copy=True):
        new_base = Basis3D.Rotation(self, axis, angle, True)
        if copy:
            new_frame = Frame3D(self.origin, new_base.u, new_base.v, new_base.w, self.name)
            return new_frame
        self.u = new_base.u
        self.v = new_base.v
        self.w = new_base.w

    def Translation(self, offset, copy=True):
        if copy:
            return Frame3D(self.origin.Translation(offset), self.u, self.v, self.w, self.name)
        self.origin.Translation(offset, copy=False)

    def copy(self):
        return Frame3D(self.origin, self.u, self.v, self.w)

    def plot2d(self, x=x3D, y=y3D, ax=None, color='k'):
        if ax is None:
            fig, ax = plt.subplots()
        else:
            fig = ax.figure
            
        origin2d = self.origin.To2D(o3D, x, y)
        
        for iv, vector in enumerate(self.vectors):
            vector2D = vector.To2D(o3D, x, y)
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
    

oxyz = Frame3D(o3D, x3D, y3D, z3D)
OXYZ = Frame3D(O3D, x3D, y3D, z3D)


class Line3D(Primitive3D, Line):
    """
    Define an infinite line passing through the 2 points
    """
    def __init__(self, point1, point2, name=''):
        Primitive3D.__init__(self, basis_primitives=[point1, point2], name=name)
        self.points = [point1, point2]
        self.bounding_box = self._bounding_box()
        
    def _bounding_box(self):
        points = self.points
        
        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])
        
        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.points[0] + (self.points[1]-self.points[0]) * curvilinear_abscissa

    def MPLPlot(self, ax):
        # Line segment
        x = [p.vector[0] for p in self.points]
        y = [p.vector[1] for p in self.points]
        z = [p.vector[2] for p in self.points]
        ax.plot(x,y,z, 'ok')

        # Drawing 3 times length of segment on each side
        u = self.points[1] - self.points[0]
        x1, y1, z1 = (self.points[0] - 3*u).vector
        x2, y2, z2 = (self.points[1] + 3*u).vector
        ax.plot([x1, x2], [y1, y2], [z1, z2], color='k')

    def MinimumDistancePoints(self, other_line):
        """
        Returns the points on this line and the other line that are the closest
        of lines
        """
        u = self.points[1] - self.points[0]
        v = other_line.points[1] - other_line.points[0]
        w = self.points[0] - other_line.points[0]
        a = u.Dot(u)
        b = u.Dot(v)
        c = v.Dot(v)
        d = u.Dot(w)
        e = v.Dot(w)

        s = (b*e -c*d) / (a*c - b**2)
        t = (a*e -b*d) / (a*c - b**2)
        p1 = self.points[0] + s*u
        p2 = other_line.points[0] + t*v
        return p1, p2

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            return Line3D(*[p.Rotation(center, axis, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            return Line3D(*[p.Translation(offset, copy=True) for p in self.points])
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
                    self.points = [frame.OldCoordinates(p) for p in self.points]
        if side == 'new':
            if copy:
                return Line3D(*[frame.NewCoordinates(p) for p in self.points])
            else:
                for p in self.points:
                    self.points = [frame.NewCoordinates(p) for p in self.points]
                    
    def copy(self):
        return Line3D(*[p.Copy() for p in self.points])
        
    @classmethod
    def from_step(cls, arguments, object_dict):
        point1 = object_dict[arguments[1]]
        direction = object_dict[arguments[2]]
        point2 = point1 + direction
        return cls(point1, point2, arguments[0][1:-1])


class BSplineCurve3D(Primitive3D):
    def __init__(self, degree, control_points, knot_multiplicities, knots, weights=None, periodic=False, name=''):
        Primitive3D.__init__(self, basis_primitives=control_points, name=name)
        self.control_points = control_points
        self.degree = degree
        knots = standardize_knot_vector(knots)
        self.knots = knots
        self.knot_multiplicities = knot_multiplicities
        self.weights = weights
        self.periodic = periodic
        self.name = name

        curve = NURBS.Curve()
        curve.degree = degree
        if weights is None:
            P = [(control_points[i][0], control_points[i][1], control_points[i][2]) for i in range(len(control_points))]
            curve.ctrlpts = P
        else:
            Pw = [(control_points[i][0]*weights[i], control_points[i][1]*weights[i], control_points[i][2]*weights[i], weights[i]) for i in range(len(control_points))]
            curve.ctrlptsw = Pw
        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot]*knot_multiplicities[i])
        curve.knotvector = knot_vector
        curve.delta = 0.01
        curve_points = curve.evalpts

        self.curve = curve
        self.points = [Point3D((p[0], p[1], p[2])) for p in curve_points]

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        points = '['
        for i in range(len(self.control_points)):
            point = 'fc.Vector({},{},{}),'.format(self.control_points[i][0],self.control_points[i][1],self.control_points[i][2])
            points += point
        points = points[:-1]
        points += ']'
        # !!! : A QUOI SERT LE DERNIER ARG DE BSplineCurve (False)?
        # LA MULTIPLICITE EN 3e ARG ET LES KNOTS EN 2e ARG ?
        return '{} = Part.BSplineCurve({},{},{},{},{},{},{})\n'.format(name,points,self.knot_multiplicities,self.knots,self.periodic,self.degree,self.weights,False)

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]
        degree = int(arguments[1])
        points = [object_dict[int(i[1:])] for i in arguments[2]]
        curve_form = arguments[3]
        if arguments[4] == '.F.':
            closed_curve = False
        elif arguments[4] == '.T.':
            closed_curve = True
        else:
            raise ValueError
        self_intersect = arguments[5]
        knot_multiplicities = [int(i) for i in arguments[6][1:-1].split(",")]
        knots = [float(i) for i in arguments[7][1:-1].split(",")]
        knot_spec = arguments[8]
        knot_vector = []
        for i, knot in enumerate(knots):
            knot_vector.extend([knot]*knot_multiplicities[i])

        if 9 in range(len(arguments)):
            weight_data = [float(i) for i in arguments[9][1:-1].split(",")]
        else:
            weight_data = None

        # FORCING CLOSED_CURVE = FALSE:
        closed_curve = False
        return cls(degree, points, knot_multiplicities, knots, weight_data, closed_curve, name)

    def point_distance(self, pt1):
        distances = []
        for point in self.points:
#            vmpt = Point3D((point[1], point[2], point[3]))
            distances.append(pt1.PointDistance(point))
        return min(distances)

    def Rotation(self, center, axis, angle, copy=True):
        new_control_points = [p.Rotation(center, axis, angle, True) for p in self.control_points]
        new_BSplineCurve3D = BSplineCurve3D(self.degree, new_control_points, self.knot_multiplicities, self.knots, self.weights, self.periodic, self.name)
        if copy:
            return new_BSplineCurve3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineCurve3D.curve
            self.points = new_BSplineCurve3D.points

    def Translation(self, offset, copy=True):
        new_control_points = [p.Translation(offset, True) for p in self.control_points]
        new_BSplineCurve3D = BSplineCurve3D(self.degree, new_control_points, self.knot_multiplicities, self.knots, self.weights, self.periodic, self.name)
        if copy:
            return new_BSplineCurve3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineCurve3D.curve
            self.points = new_BSplineCurve3D.points

class Circle3D(Primitive3D):
    def __init__(self, center, radius, normal, name=''):
        self.center = center
        self.radius = radius
        self.normal = normal
        Primitive3D.__init__(self, basis_primitives=self.tessellation_points(), name=name)
        
    def tessellation_points(self, resolution=20):
        plane = Plane3D.from_normal(self.center, self.normal)
        center_2D = self.center.To2D(plane.origin, plane.vectors[0], plane.vectors[1])
        circle2D = Circle2D(center_2D, self.radius)
        tessellation_points_2D = circle2D.tessellation_points()
        tessellation_points_3D = [p.To3D(plane.origin, x3D, y3D) for p in tessellation_points_2D]
        return tessellation_points_3D
    
    def Length(self):
        return 2* math.pi * self.radius

    def FreeCADExport(self,ip,ndigits=3):
        name = 'primitive{}'.format(ip)
        xc,yc,zc = npy.round(1000*self.center.vector,ndigits)
        xn,yn,zn = npy.round(self.normal.vector,ndigits)
        return '{} = Part.Circle(fc.Vector({},{},{}),fc.Vector({},{},{}),{})\n'.format(name,xc,yc,zc,xn,yn,zn,1000*self.radius)

    def Rotation(self, rot_center, axis, angle, copy=True):
        new_center = self.center.Rotation(rot_center, axis, angle, True)
        new_normal = self.normal.Rotation(rot_center, axis, angle, True)
        if copy:
            return Circle3D(new_center, self.radius, new_normal, self.name)
        else:
            self.center = new_center
            self.normal = new_normal

    def Translation(self, offset, copy=True):
        new_center = self.center.Translation(offset, True)
        new_normal = self.normal.Translation(offset, True)
        if copy:
            return Circle3D(new_center, self.radius, new_normal, self.name)
        else:
            self.center = new_center
            self.normal = new_normal

    @classmethod
    def from_step(cls, arguments, object_dict):
        center = object_dict[arguments[1]].origin
        radius = float(arguments[2])/1000
        normal = object_dict[arguments[1]].w
        return cls(center, radius, normal, arguments[0][1:-1])


class Ellipse3D(Primitive3D):
    def __init__(self, major_axis, minor_axis, center, normal, major_dir, name=''):
        Primitive3D.__init__(self, basis_primitives=None, name=name)
        self.major_axis = major_axis
        self.minor_axis = minor_axis
        self.center = center
        self.normal = normal
        major_dir.Normalize()
        self.major_dir = major_dir

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        xc, yc, zc = npy.round(1000*self.center.vector, ndigits)
        major_vector = self.center + self.major_axis/2 * self.major_dir
        xmaj, ymaj, zmaj = npy.round(1000*major_vector.vector, ndigits)
        minor_vector = self.center + self.minor_axis/2 * self.normal.Cross(self.major_dir)
        xmin, ymin, zmin = npy.round(1000*minor_vector.vector, ndigits)
        return '{} = Part.Ellipse(fc.Vector({},{},{}), fc.Vector({},{},{}), fc.Vector({},{},{}))\n'.format(name,xmaj,ymaj,zmaj,xmin,ymin,zmin,xc,yc,zc)

    def Rotation(self, rot_center, axis, angle, copy=True):
        new_center = self.center.Rotation(rot_center, axis, angle, True)
        new_normal = self.normal.Rotation(rot_center, axis, angle, True)
        new_major_dir = self.major_dir.Rotation(rot_center, axis, angle, True)
        if copy:
            return Ellipse3D(self.major_axis, self.minor_axis, new_center, new_normal, new_major_dir, self.name)
        else:
            self.center = new_center
            self.normal = new_normal
            self.major_dir = new_major_dir

    def Translation(self, offset, copy=True):
        new_center = self.center.Translation(offset, True)
        new_normal = self.normal.Translation(offset, True)
        new_major_dir = self.major_dir.Translation(offset, True)
        if copy:
            return Ellipse3D(self.major_axis, self.minor_axis, new_center, new_normal, new_major_dir, self.name)
        else:
            self.center = new_center
            self.normal = new_normal
            self.major_dir = new_major_dir

    @classmethod
    def from_step(cls, arguments, object_dict):
        center = object_dict[arguments[1]].origin
        normal = object_dict[arguments[1]].w
        major_dir = object_dict[arguments[1]].u
        major_axis = float(arguments[2])/1000
        minor_axis = float(arguments[3])/1000
        return cls(major_axis, minor_axis, center, normal, major_dir, arguments[0][1:-1])


class Arc3D(Primitive3D):
    """
    An arc is defined by a starting point, an end point and an interior point
    """
    def __init__(self, start, interior, end, name=''):
        self.start = start
        self.interior = interior
        self.end = end
        Primitive3D.__init__(self, basis_primitives=self.tessellation_points(), name=name)
        
        u1 = (self.interior - self.start)
        u2 = (self.interior - self.end)
        u1.Normalize()
        u2.Normalize()

        n = u2.Cross(u1)
        n.Normalize()
        self.normal = n
        v1 = n.Cross(u1)# v1 is normal
        v2 = n.Cross(u2)

        p11 = 0.5 * (start + interior)# Mid point of segment s,m
        p12 = p11 + v1
        p21 = 0.5 * (end + interior)# Mid point of segment s,m
        p22 = p21 + v2

        l1 = Line3D(p11, p12)
        l2 = Line3D(p21, p22)

        c1, _ = l1.MinimumDistancePoints(l2)
        self.center = c1
        self.radius = (self.center - self.start).Norm()

        # Determining angle

        r1 = (self.start).To2D(self.center, u1, v1)
        r2 = (self.end).To2D(self.center, u1, v1)
        ri = (self.interior).To2D(self.center, u1, v1)

        angle1 = npy.arctan2(r1.vector[1], r1.vector[0])
        angle2 = npy.arctan2(r2.vector[1], r2.vector[0])

        anglei = npy.arctan2(ri.vector[1], ri.vector[0])
        order = [y for x, y in sorted(zip([angle1, anglei, angle2], [0, 1, 2]))]
        order = order*2
        i = order.index(0)
        if order[i+1] == 1:
            # Trigo wise angle should be plus
            self.angle1 = angle1
            self.angle2 = angle2
            if angle1 > angle2:
                self.angle = angle2 - angle1 + 2 * math.pi
            else:
                self.angle = angle2 - angle1
#            self.angle = abs(self.angle2 - self.angle1)
        else:
            # Clock wise
            self.angle1 = angle2
            self.angle2 = angle1
            if angle1 < angle2:
                self.angle = -(angle2 - angle1 + 2 * math.pi)
            else:
                self.angle = -(angle2 - angle1)

    def _get_points(self):
        return [self.start,self.interior,self.end]

    points=property(_get_points)
    
    def tessellation_points(self, resolution=10):
        plane = Plane3D.from_3_points(self.interior, self.start, self.end)
        interior_2D = self.interior.To2D(plane.origin, plane.vectors[0], plane.vectors[1])
        start_2D = self.start.To2D(plane.origin, plane.vectors[0], plane.vectors[1])
        end_2D = self.end.To2D(plane.origin, plane.vectors[0], plane.vectors[1])
        arc2D = Arc2D(start_2D, interior_2D, end_2D)
        tessellation_points_2D = arc2D.tessellation_points()
#        ax = interior_2D.MPLPlot()
#        start_2D.MPLPlot(ax=ax)
#        end_2D.MPLPlot(ax=ax)
#        for pt in tessellation_points_2D:
#            pt.MPLPlot(ax=ax)
#        ax.set_aspect('equal')
#        tessellation_points_3D = [p.To3D(self.interior, self.start, self.end) for p in tessellation_points_2D]
        tessellation_points_3D = [p.To3D(plane.origin, plane.vectors[0], plane.vectors[1]) for p in tessellation_points_2D]
        return tessellation_points_3D

    def Length(self):
        return self.radius * abs(self.angle)

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.start.Rotation(self.center, self.normal, curvilinear_abscissa/self.radius)

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = Axes3D(fig)
        else:
            fig = None

        ax.plot(*self.center.vector,color='b')
        ax.plot(*self.start.vector,c='r')
        ax.plot(*self.end.vector,c='r')
        ax.plot(*self.interior.vector,c='g')
        x = []
        y = []
        z = []
        l = self.Length()
        for i in range(31):
            p = self.PointAtCurvilinearAbscissa(l*(i)/30)
            x.append(p[0])
            y.append(p[1])
            z.append(p[2])

        ax.plot(x, y, z, 'k')

    def MPLPlot2D(self, x3d, y3D, ax, color='k'):
        # TODO: Enhance this plot
        l = self.Length()
        x = []
        y = []
        for i in range(30):
            p = self.PointAtCurvilinearAbscissa(i/(29.)*l)
            xi, yi = p.PlaneProjection2D(x3D, y3D)
            x.append(xi)
            y.append(yi)
        ax.plot(x, y, color=color)

    def FreeCADExport(self, name, ndigits=6):
        xs, ys, zs = npy.round(1000*self.start.vector, ndigits)
        xm, ym, zm = npy.round(1000*self.interior.vector, ndigits)
        xe, ye, ze = npy.round(1000*self.end.vector, ndigits)
        return '{} = Part.Arc(fc.Vector({},{},{}),fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,xs,ys,zs,xm,ym,zm,xe,ye,ze)


class BSplineSurface3D(Primitive3D):
    def __init__(self, degree_u, degree_v, control_points, nb_u, nb_v, u_multiplicities, v_multiplicities, u_knots, v_knots, weights=None, name=''):
        Primitive3D.__init__(self, basis_primitives=control_points, name=name)
        self.control_points = control_points
        self.degree_u = degree_u
        self.degree_v = degree_v
        self.nb_u = nb_u
        self.nb_v = nb_v
        u_knots = standardize_knot_vector(u_knots)
        v_knots = standardize_knot_vector(v_knots)
        self.u_knots = u_knots
        self.v_knots = v_knots
        self.u_multiplicities = u_multiplicities
        self.v_multiplicities = v_multiplicities
        self.weights = weights


        self.control_points_table = []
        points_row = []
        i = 1
        for pt in control_points:
            points_row.append(pt)
            if i == nb_v:
                self.control_points_table.append(points_row)
                points_row = []
                i = 1
            else:
                i += 1
        # TRANSPOSE THE LIST OF LISTS
#        self.control_points_table = list(map(list, zip(*self.control_points_table)))

        surface = NURBS.Surface()
        surface.degree_u = degree_u
        surface.degree_v = degree_v
        if weights is None:
            P = [(control_points[i][0], control_points[i][1], control_points[i][2]) for i in range(len(control_points))]
            surface.set_ctrlpts(P, nb_u, nb_v)
        else:
            Pw = [(control_points[i][0]*weights[i], control_points[i][1]*weights[i], control_points[i][2]*weights[i], weights[i]) for i in range(len(control_points))]
            surface.set_ctrlpts(Pw, nb_u, nb_v)
        knot_vector_u = []
        for i, u_knot in enumerate(u_knots):
            knot_vector_u.extend([u_knot]*u_multiplicities[i])
        knot_vector_v = []
        for i, v_knot in enumerate(v_knots):
            knot_vector_v.extend([v_knot]*v_multiplicities[i])
        surface.knotvector_u = knot_vector_u
        surface.knotvector_v = knot_vector_v
        surface.delta = 0.01
        surface_points = surface.evalpts

        self.surface = surface
        self.points = [Point3D((p[0], p[1], p[2])) for p in surface_points]


    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        script = ""
        points = '['
        for i, pts_row in enumerate(self.control_points_table):
            pts = '['
            for j, pt in enumerate(pts_row):
                point = 'fc.Vector({},{},{}),'.format(pt[0],pt[1],pt[2])
                pts += point
            pts = pts[:-1] + '],'
            points += pts
        points = points[:-1] + ']'

        script += '{} = Part.BSplineSurface()\n'.format(name)
        if self.weights is None:
            script += '{}.buildFromPolesMultsKnots({},{},{},udegree={},vdegree={},uknots={},vknots={})\n'.format(name,points,self.u_multiplicities,self.v_multiplicities,self.degree_u,self.degree_v,self.u_knots,self.v_knots)
        else:
            script += '{}.buildFromPolesMultsKnots({},{},{},udegree={},vdegree={},uknots={},vknots={},weights={})\n'.format(name,points,self.u_multiplicities,self.v_multiplicities,self.degree_u,self.degree_v,self.u_knots,self.v_knots,self.weights)

        return script

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]
        degree_u = int(arguments[1])
        degree_v = int(arguments[2])
        points_sets = arguments[3][1:-1].split("),")
        points_sets = [elem+")" for elem in points_sets[:-1]]+[points_sets[-1]]
        control_points = []
        for points_set in points_sets:
            points = [object_dict[int(i[1:])] for i in points_set[1:-1].split(",")]
            nb_v = len(points)
            control_points.extend(points)
        nb_u = int(len(control_points) / nb_v)
        surface_form = arguments[4]
        if arguments[5] == '.F.':
            u_closed = False
        elif arguments[5] == '.T.':
            u_closed = True
        else:
            raise ValueError
        if arguments[6] == '.F.':
            v_closed = False
        elif arguments[6] == '.T.':
            v_closed = True
        else:
            raise ValueError
        self_intersect = arguments[7]
        u_multiplicities = [int(i) for i in arguments[8][1:-1].split(",")]
        v_multiplicities = [int(i) for i in arguments[9][1:-1].split(",")]
        u_knots = [float(i) for i in arguments[10][1:-1].split(",")]
        v_knots = [float(i) for i in arguments[11][1:-1].split(",")]
        knot_spec = arguments[12]

        if 13 in range(len(arguments)):
            weight_data = [float(i) for i in arguments[13][1:-1].replace("(", "").replace(")", "").split(",")]
        else:
            weight_data = None

        return cls(degree_u, degree_v, control_points, nb_u, nb_v, u_multiplicities, v_multiplicities, u_knots, v_knots, weight_data, name)

    def Rotation(self, center, axis, angle, copy=True):
        new_control_points = [p.Rotation(center, axis, angle, True) for p in self.control_points]
        new_BSplineSurface3D = BSplineSurface3D(self.degree_u, self.degree_v, new_control_points, self.nb_u, self.nb_v, self.u_multiplicities, self.v_multiplicities, self.u_knots, self.v_knots, self.weights, self.name)
        if copy:
            return new_BSplineSurface3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineSurface3D.curve
            self.points = new_BSplineSurface3D.points

    def Translation(self, offset, copy=True):
        new_control_points = [p.Translation(offset, True) for p in self.control_points]
        new_BSplineSurface3D = BSplineSurface3D(self.degree_u, self.degree_v, new_control_points, self.nb_u, self.nb_v, self.u_multiplicities, self.v_multiplicities, self.u_knots, self.v_knots, self.weights, self.name)
        if copy:
            return new_BSplineSurface3D
        else:
            self.control_points = new_control_points
            self.curve = new_BSplineSurface3D.curve
            self.points = new_BSplineSurface3D.points

class CompositePrimitive3D(Primitive3D):
    """
    A collection of simple primitives3D
    """
    def __init__(self, primitives, name=''):
        self.primitives = primitives
        basis_primitives=[]
        for primitive in primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.basis_primitives)
            else:
                basis_primitives.append(primitive)
        
        Primitive3D.__init__(self, basis_primitives=basis_primitives, name=name)

    def UpdateBasisPrimitives(self):
        # TODO: This is a copy/paste from CompositePrimitive2D, in the future make a Common abstract class
        basis_primitives=[]
        for primitive in self.primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.basis_primitives)
            else:
                basis_primitives.append(primitive)

        self.basis_primitives = basis_primitives

    def MPLPlot(self, ax = None):
        if ax is None:
            fig = plt.figure()
#            ax = fig.add_subplot(111, projection='3d', adjustable='box')
            ax = Axes3D(fig)
        else:
            fig = None

        for primitive in self.edges:
            primitive.MPLPlot(ax)

        ax.set_aspect('equal')

        return fig, ax


class Wire3D(CompositePrimitive3D):
    """
    A collection of simple primitives, following each other making a wire
    """
    def __init__(self, primitives, name=''):
        CompositePrimitive3D.__init__(self, primitives, name)

    def Length(self):
        length = 0.
        for primitive in self.basis_primitives:
            length += primitive.Length()
        return length

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        length = 0.
        for primitive in self.basis_primitives:
            primitive_length = primitive.Length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError


    # TODO: method to check if it is a wire
    def FreeCADExport(self, ip):
        name='primitive'+str(ip)

        s = 'E = []\n'
        for ip, primitive in enumerate(self.basis_primitives):
            s += primitive.FreeCADExport('L{}'.format(ip))
            s += 'E.append(Part.Edge(L{}))\n'.format(ip)
        s += '{} = Part.Wire(E[:])\n'.format(name)

        return s
            

class Edge3D(Primitive3D):
    def __init__(self, edge_start, edge_end, name=''):
        Primitive3D.__init__(self, basis_primitives=[edge_start, edge_end], name=name)
        self.points = [edge_start, edge_end]

    @classmethod
    def from_step(cls, arguments, object_dict):
        return LineSegment3D(object_dict[arguments[1]], object_dict[arguments[2]], arguments[0][1:-1])

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_edge_start = self.points[0].Rotation(center, axis, angle, copy=True)
            new_edge_end = self.points[1].Rotation(center, axis, angle, copy=True)
            return Edge3D(new_edge_start, new_edge_end)
        else:
            self.points[0].Rotation(center, axis, angle, copy=False)
            self.points[1].Rotation(center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            new_edge_start = self.points[0].Translation(offset, copy=True)
            new_edge_end = self.points[1].Translation(offset, copy=True)
            return Edge3D(new_edge_start, new_edge_end)
        else:
            self.points[0].Translation(offset, copy=False)
            self.points[1].Translation(offset, copy=False)
            
    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_edge_start = self.points[0].frame_mapping(frame, side, copy=True)
            new_edge_end = self.points[1].frame_mapping(frame, side, copy=True)
            return Edge3D(new_edge_start, new_edge_end)
        else:
            self.points[0].frame_mapping(frame, side, copy=False)
            self.points[1].frame_mapping(frame, side, copy=False)
            
    def copy(self):
        new_edge_start = self.points[0].copy()
        new_edge_end = self.points[1].copy()
        return Edge3D(new_edge_start, new_edge_end)


class LineSegment3D(Edge3D):
    """
    Define a line segment limited by two points
    """
    def __init__(self, point1, point2, name=''):
        Edge3D.__init__(self, point1, point2, name='')
        self.bounding_box = self._bounding_box()
        
    def _bounding_box(self):
        points = self.points
        
        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])
        
        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def Length(self):
        return self.points[1].PointDistance(self.points[0])

    def PlaneProjection2D(self, x, y):
        return LineSegment2D(self.points[0].PlaneProjection2D(x, y),
                             self.points[1].PlaneProjection2D(x, y))

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            return LineSegment3D(*[p.Rotation(center, axis, angle, copy=True) for p in self.points])
        else:
            Edge3D.Rotation(self, center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def Translation(self, offset, copy=True):
        if copy:
            return LineSegment3D(*[p.Translation(offset, copy=True) for p in self.points])
        else:
            Edge3D.Translation(self, offset, copy=False)
            self.bounding_box = self._bounding_box()
        
    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            if copy:
                return LineSegment3D(*[frame.OldCoordinates(p) for p in self.points])
            else:
                Edge3D.frame_mapping(self, frame, side, copy=False)
                self.bounding_box = self._bounding_box()
        if side == 'new':
            if copy:
                return LineSegment3D(*[frame.NewCoordinates(p) for p in self.points])
            else:
                Edge3D.frame_mapping(self, frame, side, copy=False)
                self.bounding_box = self._bounding_box()
                
    def copy(self):
        return LineSegment3D(*[p.copy() for p in self.points])

    def MPLPlot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        x=[p.vector[0] for p in self.points]
        y=[p.vector[1] for p in self.points]
        z=[p.vector[2] for p in self.points]
        ax.plot(x,y,z, 'o-k')
        return ax

    def MPLPlot2D(self, x_3D, y_3D, ax, color='k'):
        edge2D =  self.PlaneProjection2D(x_3D, y_3D)
        edge2D.MPLPlot(ax=ax, color=color)

    def plot_data(self, x_3D, y_3D, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        edge2D =  self.PlaneProjection2D(x_3D, y_3D)
        return edge2D.plot_data(marker, color, stroke_width,
                         dash, opacity, arrow)

    def FreeCADExport(self, name, ndigits=6):
        x1, y1, z1 = npy.round(1000*self.points[0].vector, ndigits)
        x2, y2, z2 = npy.round(1000*self.points[1].vector, ndigits)
        return '{} = Part.LineSegment(fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,x1,y1,z1,x2,y2,z2)

    def to_line(self):
        return Line3D(*self.points)        
    
    def Babylon(self, color=(1, 1, 1), name='line',  type_='line', parent=None):
        if type_ == 'line':
            s = 'var myPoints = [];\n'
            s += 'var point1 = new BABYLON.Vector3({},{},{});\n'.format(*self.points[0])
            s += 'myPoints.push(point1);\n'
            s += 'var point2 = new BABYLON.Vector3({},{},{});\n'.format(*self.points[1])
            s += 'myPoints.push(point2);\n'
            s += 'var {} = BABYLON.MeshBuilder.CreateLines("lines", {{points: myPoints}}, scene);\n'.format(name)
            s += '{}.color = new BABYLON.Color3{};\n'.format(name, tuple(color))
        elif type_ == 'tube':
            radius = 0.03*self.points[0].PointDistance(self.points[1])
            s = 'var points = [new BABYLON.Vector3({},{},{}), new BABYLON.Vector3({},{},{})];\n'.format(*self.points[0], *self.points[1])
            s += 'var {} = BABYLON.MeshBuilder.CreateTube("frame_U", {{path: points, radius: {}}}, {});'.format(name, radius, parent)
#            s += 'line.material = red_material;\n'

        else:
            raise NotImplementedError

        if parent is not None:
            s += '{}.parent = {};\n'.format(name, parent)

        return s 


class Contour3D(Wire3D):
    """
    A collection of 3D primitives forming a closed wire3D
    """
    def __init__(self, edges, points=None, name=''):
        """
        Faire un choix : soit edges c'est un CompositePrimitives3D
        ou alors un ensemble de primitives
        ou alors un ensemble de basis_primtives (qui sont des points pour le moment)
        """
        self.name = name
        self.points = points
        
        edges_primitives = []
        for edge in edges:
            if edge.__class__ == CompositePrimitive3D:
                edges_primitives.extend(edge.primitives)
            else:
                edges_primitives.append(edge)
        self.edges = edges_primitives
                    
        self.points = self.clean_points()


    @classmethod
    def from_step(cls, arguments, object_dict):
        edges = []
        for edge in arguments[1]:
            edges.append(object_dict[int(edge[1:])])

#        points = edges[0].points[:]                
#        for i, edge in enumerate(edges[1:-1]):
#            if edge.points[0] in points[-2:]:
#                points.append(edge.points[1])
#            elif edge.points[1] in points[-2:]:
#                points.append(edge.points[0])
#            else:
#                raise NotImplementedError
#        contour_points = [p.copy() for p in points]

        return cls(edges, points=None, name=arguments[0][1:-1])
    
    def clean_points(self):
        points = self.edges[0].basis_primitives[::]
        last_points_added = points
        for edge in self.edges[1:]:
            if hasattr(edge, 'basis_primitives'):
                points_to_add = edge.basis_primitives
                if points_to_add[0] in [last_points_added[0], last_points_added[-1]]:
                    points.extend(points_to_add[1:])
                elif points_to_add[-1] in [last_points_added[0], last_points_added[-1]]:
                    points.extend(points_to_add[:-1])
                else:
                    self.MPLPlot()
                    raise NotImplementedError
                last_points_added = points_to_add
            else:
                raise NotImplementedError
        return points
        

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_edges = [edge.Rotation(center, axis, angle, copy=True) for edge in self.edges]
            new_points = [p.Rotation(center, axis, copy=True) for p in self.points]
            return Contour3D(new_edges, new_points, self.name)
        else:
            for edge in self.edges:
                edge.Rotation(center, axis, angle, copy=False)
            for point in self.points:
                point.Rotation(center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            new_edges = [edge.Translation(offset, copy=True) for edge in self.edges]
            new_points = [p.Translation(offset, copy=True) for p in self.points]
            return Contour3D(new_edges, new_points, self.name)
        else:
            for edge in self.edges:
                edge.Translation(offset, copy=False)
            for point in self.points:
                point.Translation(offset, copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_edges = [edge.frame_mapping(frame, side, copy=True) for edge in self.edges]
            new_points = [p.frame_mapping(frame, side, copy=True) for p in self.points]
            return Contour3D(new_edges, new_points, self.name)
        else:
            for edge in self.edges:
                edge.frame_mapping(frame, side, copy=False)
            for point in self.points:
                point.frame_mapping(frame, side, copy=False)
            
    def copy(self):
        new_edges = [edge.copy() for edge in self.edges]
        new_points = [p.copy() for p in self.points]
        return Contour3D(new_edges, new_points, self.name)


class Face3D(Primitive3D):
    def __init__(self, contours, plane=None, points=None, polygon2D=None, name=''):
#        Primitive3D.__init__(self, name=name)
        self.contours = contours
        self.plane = plane
        self.points = points
        self.polygon2D = polygon2D
        
        self.name = name
        
        contour_points = [p.copy() for p in self.contours[0].points[:]]
        if plane is None:
#            print('face3D recalcule self.plane')
            self.plane = Plane3D.from_points(contour_points)
        
        if points is None or polygon2D is None:
#            print('face3D recalcule self.points et self.polygon2D')
            self.points, self.polygon2D = self._repair_points_and_polygon2d(contour_points, self.plane)
            self.contours[0].points = [p.copy() for p in self.points]
                        
        self.bounding_box = self._bounding_box()
                        
        # CHECK #
        for pt in self.points:
            if not self.plane.point_on_plane(pt):
                print('WARNING', pt, 'not on', self.plane)
                print('dot =', self.plane.normal.Dot(pt-self.plane.origin))
                raise ValueError
        
    @classmethod
    def from_step(cls, arguments, object_dict):                    
        contours = []
        contours.append(object_dict[int(arguments[1][0][1:])])
        
        plane = Plane3D.from_points(contours[0].points)
        contours[0].points, polygon2D = cls._repair_points_and_polygon2d(contours[0].points, plane)
        points = [p.copy() for p in contours[0].points[:]]
        
        return cls(contours, plane=plane, points=points, polygon2D=polygon2D, name=arguments[0][1:-1])
    
    @classmethod
    def _repair_points_and_polygon2d(cls, points, plane):
        if points[0] == points[-1]:
            points = points[:-1]
        polygon_points = [p.To2D(plane.origin, plane.vectors[0], plane.vectors[1]) for p in points]
        repaired_points = [p.copy() for p in points]
        polygon2D = Polygon2D(polygon_points)
        if polygon2D.SelfIntersect()[0]:
            repaired_points = [repaired_points[1]]+[repaired_points[0]]+repaired_points[2:]
            polygon_points = [polygon_points[1]]+[polygon_points[0]]+polygon_points[2:]
            if polygon_points[0] == polygon_points[-1]:
                repaired_points = repaired_points[:-1]
                polygon_points = polygon_points[:-1]
            polygon2D = Polygon2D(polygon_points)
        return repaired_points, polygon2D

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_contour = [subcontour.Rotation(center, axis, angle, copy=True) for subcontour in self.contour]
            new_plane = self.plane.Rotation(center, axis, angle, copy=True)
            new_points = [p.Rotation(center, axis, angle, copy=True) for p in self.points]
            return Face3D(new_contour, new_plane, new_points, self.polygon2D, self.name)
        else:
            for contour in self.contours:
                contour.Rotation(center, axis, angle, copy=False)
            for point in self.points:
                point.Rotation(center, axis, angle, copy=False)
            self.plane.Rotation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def Translation(self, offset, copy=True):
        if copy:
            new_contour = [subcontour.Translation(offset, copy=True) for subcontour in self.contours]
            new_plane = self.plane.Translation(offset, copy=True)
            new_points = [p.Translation(offset, copy=True) for p in self.points]
            return Face3D(new_contour, new_plane, new_points, self.polygon2D, self.name)
        else:
            for contour in self.contours:
                contour.Translation(offset, copy=False)
            for point in self.points:
                point.Translation(offset, copy=False)
            self.plane.Translation(offset, copy=False)
            self.bounding_box = self._bounding_box()
            
    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_contour = [subcontour.frame_mapping(frame, side, copy=True) for subcontour in self.contours]
            new_plane = self.plane.frame_mapping(frame, side, copy=True)
            new_points = [p.frame_mapping(frame, side, copy=True) for p in self.points]
            return Face3D(new_contour, new_plane, new_points, None, self.name)
        else:
            for contour in self.contours:
                contour.frame_mapping(frame, side, copy=False)
            for point in self.points:
                point.frame_mapping(frame, side, copy=False)
            self.plane.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()
            
    def copy(self):
        new_contour = [contour.copy() for contour in self.contours]
        new_plane = self.plane.copy()
        new_points = [p.copy() for p in self.points]
        return Face3D(new_contour, new_plane, new_points, self.polygon2D, self.name)
        
    def average_center_point(self):
        """
        excluding holes
        """
        points = self.points
        nb = len(points)
        x = npy.sum([p[0] for p in points]) / nb
        y = npy.sum([p[1] for p in points]) / nb    
        z = npy.sum([p[2] for p in points]) / nb
        return Point3D((x,y,z))
    
    def triangulation(self):
        points = [tuple(p.vector) for p in self.polygon2D.points]
        points_dict = {p: i for i, p in enumerate(points)}
        points_3D = [Point2D(p).To3D(self.plane.origin, self.plane.vectors[0], self.plane.vectors[1]) for p in points]
        triangles = tri.earclip(points)
        triangles_indexes = []
        for triangle in triangles:
            triangle_indexes = []
            for point in triangle:
                triangle_indexes.append(points_dict[point])
            triangles_indexes.append(triangle_indexes)
        return points_3D, triangles_indexes
        
    def _bounding_box(self):
        points = self.points
        
        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])
                
        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def distance_to_point(self, point, return_other_point=False):
        """
        Only works if the surface is planar
        TODO : this function does not take into account if Face has holes
        """
        # On projette le point sur la surface plane
        # Si le point est à l'intérieur de la face, on retourne la distance de projection
        # Si le point est à l'extérieur, on projette le point sur le plan
        # On calcule en 2D la distance entre la projection et le polygone contour
        # On utilise le theroeme de Pytagore pour calculer la distance minimale entre le point et le contour

        projected_pt = point.PlaneProjection3D(self.plane.origin, self.plane.vectors[0], self.plane.vectors[1])
        projection_distance = point.PointDistance(projected_pt)

        if self.point_on_face(projected_pt):
            if return_other_point:
                return projection_distance, projected_pt
            return projection_distance
            
        point_2D = point.To2D(self.plane.origin, self.plane.vectors[0], self.plane.vectors[1])

        border_distance, other_point = self.polygon2D.PointBorderDistance(point_2D, return_other_point=True)

        other_point = other_point.To3D(self.plane.origin , self.plane.vectors[0], self.plane.vectors[1])        
        
        if return_other_point:
            return (projection_distance**2 + border_distance**2)**0.5, other_point
        return (projection_distance**2 + border_distance**2)**0.5

    def distance_to_face(self, face2, return_points=False):
        """
        Only works if the surface is planar
        TODO : this function does not take into account if Face has holes
        TODO : TRAITER LE CAS OU LA DISTANCE LA PLUS COURTE N'EST PAS D'UN SOMMET
        """
        # On calcule la distance entre la face 1 et chaque point de la face 2
        # On calcule la distance entre la face 2 et chaque point de la face 1
        
        if self.face_intersection(face2) is not None:
            return 0, None, None
            
        polygon1_points_3D = [Point3D(p.vector) for p in self.points]
        polygon2_points_3D = [Point3D(p.vector) for p in face2.points]

        distances = []
        if not return_points:
            d_min = face2.distance_to_point(polygon1_points_3D[0])
            for point1 in polygon1_points_3D[1:]:
                d = face2.distance_to_point(point1)
                if d < d_min:
                    d_min = d
            for point2 in polygon2_points_3D:
                d = self.distance_to_point(point2)
                if d < d_min:
                    d_min = d
            return d_min
        
        else:
            for point1 in polygon1_points_3D:
                d, other_point = face2.distance_to_point(point1, return_other_point=True)
                distances.append((d, point1, other_point))
            for point2 in polygon2_points_3D:
                d, other_point = self.distance_to_point(point2, return_other_point=True)
                distances.append((d, point2, other_point))
            
        d_min, point_min, other_point_min = distances[0]
        for distance in distances[1:]:
            if distance[0] < d_min:
                d_min = distance[0]
                point_min = distance[1]
                other_point_min = distance[2]
        
        return d_min, point_min, other_point_min

    def point_on_face(self, point):
        """
        Only works if the surface is planar
        TODO : this function does not take into account if Face has holes

        Tells you if a point is on the 3D face and inside its contour
        """
        
        point_on_plane = self.plane.point_on_plane(point)
        # The point is not in the same plane
        if not point_on_plane:
            print('point not on plane so not on face')
            return False
        
        point_2D = point.To2D(self.plane.origin, self.plane.vectors[0], self.plane.vectors[1])
        
        ### PLOT ###
#        ax = self.polygon2D.MPLPlot()
#        point_2D.MPLPlot(ax)
        ############

        if not self.polygon2D.PointBelongs(point_2D):
            ### PLOT ###
#            ax = self.polygon2D.MPLPlot()
#            point_2D.MPLPlot(ax)
#            ax.set_title('DEHORS')
#            ax.set_aspect('equal')
            ############
            return False
        ### PLOT ###
#        ax = self.polygon2D.MPLPlot()
#        point_2D.MPLPlot(ax)
#        ax.set_title('DEDANS')
#        ax.set_aspect('equal')
        ############
        return True

    def edge_intersection(self, edge):
        linesegment = LineSegment3D(*edge.points)
        intersection_point = self.plane.linesegment_intersection(linesegment)

        if intersection_point is None:
            return None

        point_on_face_boo = self.point_on_face(intersection_point)
        if not point_on_face_boo:
            return None

        return intersection_point

    def linesegment_intersection(self, linesegment, abscissea=False):
        if abscissea:
            intersection_point, intersection_abscissea = self.plane.linesegment_intersection(linesegment, True)
        else:
            intersection_point = self.plane.linesegment_intersection(linesegment)
            
        if intersection_point is None:
            if abscissea:
                return None, None
            return None
        point_on_face_boo = self.point_on_face(intersection_point)
        ### PLOT ###
#        if point_on_face_boo:
#            ax = self.plot()
#            linesegment.MPLPlot(ax)
#            Point3D(intersection_point.vector).MPLPlot(ax)
#            self.plane.MPLPlot(ax)
#            ax.set_aspect('equal')
        
#            print('=>', self.plane.normal.Dot(intersection_point-self.plane.origin))
#            print('point_on_face_boo', point_on_face_boo)
        ############
        if not point_on_face_boo:
            if abscissea:
                return None, None
            return None
        
        if abscissea:
            return intersection_point, intersection_abscissea
        return intersection_point

    def face_intersection(self, face2):
        """
        Only works if the surface is planar
        TODO : this function does not take into account if Face has holes
        """
        bbox1 = self.bounding_box
        bbox2 = face2.bounding_box
        if not bbox1.bbox_intersection(bbox2):
            return None
        
        intersection_points = []

        for edge2 in face2.contours[0].edges:
            intersection_point = self.edge_intersection(edge2)
            if intersection_point is not None:
                intersection_points.append(intersection_point)

        for edge1 in self.contours[0].edges:
            intersection_point = face2.edge_intersection(edge1)
            if intersection_point is not None:
                intersection_points.append(intersection_point)

        if not intersection_points:
            return None
        
        return intersection_points
    
    def plot(self, ax=None):
        fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')

#        x = [p[0] for edge in self.contours[0].edges for p in edge.points]
#        y = [p[1] for edge in self.contours[0].edges for p in edge.points]
#        z = [p[2] for edge in self.contours[0].edges for p in edge.points]
#        print(x,y,z)
        x = [p[0] for p in self.contours[0].points]
        y = [p[1] for p in self.contours[0].points]
        z = [p[2] for p in self.contours[0].points]
        
        ax.scatter(x, y, z)
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        for edge in self.contours[0].edges:
            for point1, point2 in (edge.points, edge.points[1:]+[edge.points[0]]):
                xs = [point1[0], point2[0]]
                ys = [point1[1], point2[1]]
                zs = [point1[2], point2[2]]
                line = mpl_toolkits.mplot3d.art3d.Line3D(xs, ys, zs)
                ax.add_line(line)

        plt.show()
        return ax
    

class Shell3D(CompositePrimitive3D):
    def __init__(self, faces, name='', color=None):
        self.faces = faces
        self.name = name
        self.color = color
        self.bounding_box = self._bounding_box()

    @classmethod
    def from_step(cls, arguments, object_dict):
        faces = []
        for face in arguments[1]:
            faces.append(object_dict[int(face[1:])])
        return cls(faces, arguments[0][1:-1])

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_faces = [face.Rotation(center, axis, angle, copy=True) for face in self.faces]
            return Shell3D(new_faces, self.name)
        else:
            for face in self.faces:
                face.Rotation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def Translation(self, offset, copy=True):
        if copy:
            new_faces = [face.Translation(offset, copy=True) for face in self.faces]
            return Shell3D(new_faces, self.name)
        else:
            for face in self.faces:
                face.Translation(offset, copy=False)
            self.bounding_box = self._bounding_box()
    
    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_faces = [face.frame_mapping(frame, side, copy=True) for face in self.faces]
            return Shell3D(new_faces, self.name)
        else:
            for face in self.faces:
                face.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()
            
    def copy(self):
        new_faces = [face.copy() for face in self.faces]
        return Shell3D(new_faces, self.name)
    
    def union(self, shell2):
        new_faces = [face for face in self.faces+shell2.faces] 
        new_name = self.name+' union '+shell2.name
        new_color = self.color
        return Shell3D(new_faces, new_name, new_color)

    def _bounding_box(self):
        """
        Returns the boundary box
        """
        points = []
        for face in self.faces:
            points.extend(face.bounding_box.points)
                
        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])
        
        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax, self.name)

    def point_belongs(self, point):
        """
        Ray Casting algorithm
        Returns True if the point is inside the Shell, False otherwise
        """
        epsilon = 1

        bbox = self.bounding_box
        if point[0] < bbox.xmin or point[0] > bbox.xmax:
            return False
        if point[1] < bbox.ymin or point[1] > bbox.ymax:
            return False
        if point[2] < bbox.zmin or point[2] > bbox.zmax:
            return False
        
        rays = []
        rays.append(LineSegment3D(point, Point3D((bbox.xmin-random.uniform(0, 1)*epsilon, bbox.ymin-random.uniform(0, 1)*epsilon, bbox.zmin-random.uniform(0, 1)*epsilon))))
        rays.append(LineSegment3D(point, Point3D((bbox.xmax+random.uniform(0, 1)*epsilon, bbox.ymin-random.uniform(0, 1)*epsilon, bbox.zmin-random.uniform(0, 1)*epsilon))))
        rays.append(LineSegment3D(point, Point3D((bbox.xmin-random.uniform(0, 1)*epsilon, bbox.ymax+random.uniform(0, 1)*epsilon, bbox.zmin-random.uniform(0, 1)*epsilon))))
        rays.append(LineSegment3D(point, Point3D((bbox.xmin-random.uniform(0, 1)*epsilon, bbox.ymin-random.uniform(0, 1)*epsilon, bbox.zmax+random.uniform(0, 1)*epsilon))))
        rays.append(LineSegment3D(point, Point3D((bbox.xmax+random.uniform(0, 1)*epsilon, bbox.ymax+random.uniform(0, 1)*epsilon, bbox.zmax+random.uniform(0, 1)*epsilon))))
        rays.append(LineSegment3D(point, Point3D((bbox.xmax+random.uniform(0, 1)*epsilon, bbox.ymax+random.uniform(0, 1)*epsilon, bbox.zmin-random.uniform(0, 1)*epsilon))))
        rays.append(LineSegment3D(point, Point3D((bbox.xmax+random.uniform(0, 1)*epsilon, bbox.ymin-random.uniform(0, 1)*epsilon, bbox.zmax+random.uniform(0, 1)*epsilon))))
        rays.append(LineSegment3D(point, Point3D((bbox.xmin-random.uniform(0, 1)*epsilon, bbox.ymax+random.uniform(0, 1)*epsilon, bbox.zmax+random.uniform(0, 1)*epsilon))))
        
        rays = sorted(rays, key=lambda ray: ray.Length())
        
        rays_intersections = []
        tests = []
        for ray in rays[:3]:
            count = 0
            ray_intersection = []
            is_inside = True
            for face in self.faces:
                intersection_point = face.linesegment_intersection(ray)
                if intersection_point is not None:
                    ray_intersection.append(intersection_point)
                    count += 1
            if count%2 == 0:
                is_inside = False
            tests.append(is_inside)
            rays_intersections.append(ray_intersection)
            
        if sum(tests) == 0 or sum(tests) == 3:
            return tests[0]
        else:
            print('PROBLEME')
#            print(point)
#            raise NotImplementedError
#            if tuple(point.vector) == (-0.3, 0.25, -0.25):
            ### BABYLON ###
#                print('------------------------')
#                point = Point3D(point.vector)
#                print(point.Babylon())
#                for ray in rays[:3]:
#                    print(ray.Babylon())
#                for ray in rays_intersections:
#                    for point in ray:
#                        point = Point3D(point.vector)
#                        print(point.Babylon())
#                print('------------------------')
            ###############
            
        return sum(tests) > 1
                
    def is_inside_shell(self, shell2):
        """
        Returns True if all the points of self are inside shell2 and no face \
        are intersecting
        """
        bbox1 = self.bounding_box
        bbox2 = shell2.bounding_box
        if not bbox1.is_inside_bbox(bbox2):
            return False

        points = []
        for face in self.faces:
            points.extend(face.contours[0].points)
                    
        for point in points:
            if not shell2.point_belongs(point):
                print('point belongs', point)
                return False
            
        # Check if any faces are intersecting
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
#                    print('Two faces are intersecting :', face1, face2)
                    return False
                
        return True
        
    def shell_intersection(self, shell2):
        """
        Return None if disjointed
        Return (1, 0) or (0, 1) if one is inside the other
        Return (n1, n2) if intersection
        
        4 cases :
            (n1, n2) with face intersection             => (n1, n2)
            (0, 0) with face intersection               => (0, 0)
            (0, 0) with no face intersection            => None
            (1, 0) or (0, 1) with no face intersection  => 1
        """
        # Check if boundary boxes don't intersect
        bbox1 = self.bounding_box
        bbox2 = shell2.bounding_box
        if not bbox1.bbox_intersection(bbox2):
#            print("No intersection of shells' BBox")
            return None
        
        # Check if any point of the first shell is in the second shell
        points1 = []
        for face in self.faces:
            points1.extend(face.contours[0].points)
        points2 = []
        for face in shell2.faces:
            points2.extend(face.contours[0].points)
        
        nb_pts1 = len(points1)
        nb_pts2 = len(points2)
        compteur1 = 0
        compteur2 = 0
        for point1 in points1:
            if shell2.point_belongs(point1):
                compteur1 += 1
        for point2 in points2:
            if self.point_belongs(point2):
                compteur2 += 1
                
        inter1 = compteur1/nb_pts1
        inter2 = compteur2/nb_pts2
#        print('shell intersection')
#        print('shell1 intersecte shell2 à', inter1*100, '%')
#        print('shell2 intersecte shell1 à', inter2*100, '%')
        
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
#                    print('Two faces are intersecting :', face1, face2)
#                    ax = face1.plot()
#                    face2.plot(ax)
                    return inter1, inter2
        if (inter1, inter2) == (0, 0):
            return None
        return 1
        
    def distance_to_shell(self, shell2, add_to_volumemodel=None):
        """
        Returns a Mesure object if the distance is not zero, otherwise returns None
        """
        
        if self.shell_intersection(shell2) is not None and self.shell_intersection(shell2) != 1:
            return None

        distance_min, point1_min, point2_min = self.faces[0].distance_to_face(shell2.faces[0], return_points=True)
        for face1 in self.faces:
            bbox1 = face1.bounding_box
            for face2 in shell2.faces:
                bbox2 = face2.bounding_box
                bbox_distance = bbox1.distance_to_bbox(bbox2)
                if bbox_distance < distance_min:
                    distance, point1, point2 = face1.distance_to_face(face2, return_points=True)
                    if distance == 0:
                        return None
                    elif distance < distance_min:
                        distance_min, point1_min, point2_min = distance, point1, point2
                        
        mesure = Mesure(point1_min, point2_min)
                        
        if add_to_volumemodel is not None:
            add_to_volumemodel.primitives.append(mesure)
                            
        return mesure

    def distance_to_point(self, point, add_to_volumemodel=None):
        """
        Computes the distance of a point to a Shell3D, whether it is inside or outside the Shell3D
        """
        distance_min, point1_min, point2_min = self.faces[0].distance_to_point(point, return_other_point=True)
        for face in self.faces[1:]:
            bbox_distance = self.bounding_box.distance_to_point(point)
            if bbox_distance < distance_min:
                distance, point1, point2 = face.distance_to_point(point, return_other_point=True)
                if distance < distance_min:
                    distance_min, point1_min, point2_min = distance, point1, point2
                
        mesure = Mesure(point1_min, point2_min)
                        
        if add_to_volumemodel is not None:
            add_to_volumemodel.primitives.append(mesure)
                            
        return mesure
    
    def intersection_internal_aabb_volume(self, shell2):
        """
        aabb made of the intersection points and the points of self internal to shell2
        """
        intersections_points = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    intersections_points.extend(intersection_points)
       
        shell1_points_inside_shell2 = []
        for face in self.faces:
            for point in face.points:
                if shell2.point_belongs(point):
                    shell1_points_inside_shell2.append(point)
                            
        if len(intersections_points+shell1_points_inside_shell2) == 0:
            return 0
        bbox = BoundingBox.from_points(intersections_points+shell1_points_inside_shell2)
        return bbox.volume()
    
    def intersection_external_aabb_volume(self, shell2):
        """
        aabb made of the intersection points and the points of self external to shell2
        """
        intersections_points = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    intersections_points.extend(intersection_points)
        
        shell1_points_outside_shell2 = []
        for face in self.faces:
            for point in face.points:
                if not shell2.point_belongs(point):
                    shell1_points_outside_shell2.append(point)
        
        if len(intersections_points+shell1_points_outside_shell2) == 0:
            return 0
        bbox = BoundingBox.from_points(intersections_points+shell1_points_outside_shell2)
        return bbox.volume()
    
    def Babylon(self):
        s = 'var customMesh = new BABYLON.Mesh("custom", scene);\n'
        
        positions = ''
        indices = ''
        
        nb_points = 0
        for i, face in enumerate(self.faces):
            points_3D, triangles_indexes = face.triangulation()
            
            for point in points_3D:
                positions += '{},{},{},'.format(round(point[0],6), round(point[1],6), round(point[2],6))
            
            for j, indexes in enumerate(triangles_indexes):                    
                indices += '{},{},{},'.format(indexes[0]+nb_points, indexes[1]+nb_points, indexes[2]+nb_points)
            nb_points += len(points_3D)
            
#            s += face.plane.Babylon()
                    
        positions = positions[:-1]
        indices = indices[:-1]
        
        s += 'var positions = [{}];\n'.format(positions)
        s += 'var indices = [{}];\n'.format(indices)
        s += 'var normals = [];\n'
        s += 'var vertexData = new BABYLON.VertexData();\n'
        s += 'BABYLON.VertexData.ComputeNormals(positions, indices, normals);\n'
        s += 'vertexData.positions = positions;\n'
        s += 'vertexData.indices = indices;\n'
        s += 'vertexData.normals = normals;\n'
        s += 'vertexData.applyToMesh(customMesh);\n'
        s += 'customMesh.enableEdgesRendering(0.95);\n'
        s += 'customMesh.edgesWidth = 0.1;\n'
        s += 'customMesh.edgesColor = new BABYLON.Color4(0, 0, 0, 0.6);\n'
        s += 'var mat = new BABYLON.StandardMaterial("mat", scene);\n'
#        s += 'mat.diffuseColor = BABYLON.Color3.Green();\n'
#        s += 'mat.specularColor = new BABYLON.Color3(0.5, 0.6, 0.87);\n'
#        s += 'mat.emissiveColor = new BABYLON.Color3(1, 1, 1);\n'
#        s += 'mat.ambientColor = new BABYLON.Color3(0.23, 0.98, 0.53);\n'
        s += 'mat.backFaceCulling = false;\n'
        s += 'customMesh.material = mat;\n'
        if self.color is not None:
            s += 'mat.diffuseColor = new BABYLON.Color3({}, {}, {});\n'.format(self.color[0], self.color[1], self.color[2])

        return s


class BoundingBox:
    """
    An axis aligned boundary box
    """
    def __init__(self, xmin, xmax, ymin, ymax, zmin, zmax, name=''):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax
        self.points = [Point3D((self.xmin, self.ymin, self.zmin)), \
                       Point3D((self.xmax, self.ymin, self.zmin)), \
                       Point3D((self.xmax, self.ymax, self.zmin)), \
                       Point3D((self.xmin, self.ymax, self.zmin)), \
                       Point3D((self.xmin, self.ymin, self.zmax)), \
                       Point3D((self.xmax, self.ymin, self.zmax)), \
                       Point3D((self.xmax, self.ymax, self.zmax)), \
                       Point3D((self.xmin, self.ymax, self.zmax))]
        self.center = (self.points[0]+self.points[-2])/2
        self.name = name
    
    def plot(self, ax=None, color=''):
        fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')
            
        bbox_edges =  [[self.points[0], self.points[1]], [self.points[0], self.points[3]], \
                       [self.points[0], self.points[4]], [self.points[1], self.points[2]], \
                       [self.points[1], self.points[5]], [self.points[2], self.points[3]], \
                       [self.points[2], self.points[6]], [self.points[3], self.points[7]], \
                       [self.points[4], self.points[5]], [self.points[5], self.points[6]], \
                       [self.points[6], self.points[7]], [self.points[7], self.points[4]]]

        x = [p[0] for p in self.points]
        y = [p[1] for p in self.points]
        z = [p[2] for p in self.points]
        ax.scatter(x, y, z)
        for edge in bbox_edges:
            ax.plot3D([edge[0][0], edge[1][0]],
                      [edge[0][1], edge[1][1]],
                      [edge[0][2], edge[1][2]],
                      'gray')
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        plt.show()
        return ax
    
    @classmethod
    def from_points(cls, points):
        xmin = min([pt[0] for pt in points])
        xmax = max([pt[0] for pt in points])
        ymin = min([pt[1] for pt in points])
        ymax = max([pt[1] for pt in points])
        zmin = min([pt[2] for pt in points])
        zmax = max([pt[2] for pt in points])
        return cls(xmin, xmax, ymin, ymax, zmin, zmax)
        
    def volume(self):
        return (self.xmax-self.xmin)*(self.ymax-self.ymin)*(self.zmax-self.zmin)
    
    def bbox_intersection(self, bbox2):
        return (self.xmin < bbox2.xmax and self.xmax > bbox2.xmin \
                and self.ymin < bbox2.ymax and self.ymax > bbox2.ymin \
                and self.zmin < bbox2.zmax and self.zmax > bbox2.zmin)
        
    def is_inside_bbox(self, bbox2):
        return (self.xmin > bbox2.xmin and self.xmax < bbox2.xmax \
                and self.ymin > bbox2.ymin and self.ymax < bbox2.ymax \
                and self.zmin > bbox2.zmin and self.zmax < bbox2.zmax)
        
    def intersection_volume(self, bbox2):
        if not self.bbox_intersection(bbox2):
            return 0
        
        permute_bbox1 = self
        permute_bbox2 = bbox2
        
        if permute_bbox2.xmin < permute_bbox1.xmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        lx = permute_bbox1.xmax - permute_bbox2.xmin
            
        if permute_bbox2.ymin < permute_bbox1.ymin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        ly = permute_bbox1.ymax - permute_bbox2.ymin
        
        if permute_bbox2.zmin < permute_bbox1.zmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        lz = permute_bbox1.zmax - permute_bbox2.zmin
        
        return lx*ly*lz
        
    def distance_to_bbox(self, bbox2):
        if self.bbox_intersection(bbox2):
            return 0
        
        permute_bbox1 = self
        permute_bbox2 = bbox2
        
        if permute_bbox2.xmin < permute_bbox1.xmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dx = permute_bbox2.xmin - permute_bbox1.xmax
        if dx < 0:
            dx = 0
        
        if permute_bbox2.ymin < permute_bbox1.ymin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dy = permute_bbox2.ymin - permute_bbox1.ymax
        if dy < 0:
            dy = 0
        
        if permute_bbox2.zmin < permute_bbox1.zmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dz = permute_bbox2.zmin - permute_bbox1.zmax
        if dz < 0:
            dz = 0
        
        return (dx**2 + dy**2 + dz**2)**0.5
    
    def point_belongs(self, point):
        return self.xmin < point[0] and point[0] < self.xmax \
        and self.ymin < point[1] and point[1] < self.ymax \
        and self.zmin < point[2] and point[2] < self.zmax 
    
    def distance_to_point(self, point):
        if self.point_belongs(point):
            return min([self.xmax-point[0], point[0]-self.xmin,
                        self.ymax-point[1], point[1]-self.ymin,
                        self.zmax-point[2], point[2]-self.zmin])
        else:
            if point[0] < self.xmin:
                dx = self.xmin - point[0]
            elif self.xmax < point[0]:
                dx = point[0] - self.xmax
            else:
                dx = 0
                
            if point[1] < self.ymin:
                dy = self.ymin - point[1]
            elif self.ymax < point[1]:
                dy = point[1] - self.ymax
            else:
                dy = 0
                
            if point[2] < self.zmin:
                dz = self.zmin - point[2]
            elif self.zmax < point[2]:
                dz = point[2] - self.zmax
            else:
                dz = 0
        return (dx**2 + dy**2 + dz**2)**0.5
    
    def distance_between_two_points_on_bbox(self, point1, point2):
        
        if   math.isclose(point1[0], self.xmin, abs_tol=1e-8):
            face_point1 = 5
        elif math.isclose(point1[0], self.xmax, abs_tol=1e-8):
            face_point1 = 3
        elif math.isclose(point1[1], self.ymin, abs_tol=1e-8):
            face_point1 = 4
        elif math.isclose(point1[1], self.ymax, abs_tol=1e-8):
            face_point1 = 2
        elif math.isclose(point1[2], self.zmin, abs_tol=1e-8):
            face_point1 = 6
        elif math.isclose(point1[2], self.zmax, abs_tol=1e-8):
            face_point1 = 1
        else:
            raise NotImplementedError
            
        if   math.isclose(point2[0], self.xmin, abs_tol=1e-8):
            face_point2 = 5
        elif math.isclose(point2[0], self.xmax, abs_tol=1e-8):
            face_point2 = 3
        elif math.isclose(point2[1], self.ymin, abs_tol=1e-8):
            face_point2 = 4
        elif math.isclose(point2[1], self.ymax, abs_tol=1e-8):
            face_point2 = 2
        elif math.isclose(point2[2], self.zmin, abs_tol=1e-8):
            face_point2 = 6
        elif math.isclose(point2[2], self.zmax, abs_tol=1e-8):
            face_point2 = 1
        else:
            raise NotImplementedError
        
        point1_copy = point1.copy()
        point2_copy = point2.copy()
        print(face_point1, face_point2)
        if face_point1 > face_point2:
            print('inversion')
            point1, point2 = point2, point1
            face_point1, face_point2 = face_point2, face_point1
            
        # The points are on the same face
        if face_point1 == face_point2:
            return point1.PointDistance(point2)

        deltax = self.xmax - self.xmin
        deltay = self.ymax - self.ymin
        deltaz = self.zmax - self.zmin
        
        point1_2d_coordinate_dict = {1: Point2D((point1[0]-self.xmin-deltax/2, point1[1]-self.ymin-deltay/2)),
                                     2: Point2D((point1[2]-self.zmin-deltaz/2, point1[0]-self.xmin-deltax/2)),
                                     3: Point2D((point1[1]-self.ymin-deltay/2, point1[2]-self.zmin-deltaz/2)),
                                     4: Point2D((point1[0]-self.xmin-deltax/2, point1[2]-self.zmin-deltaz/2)),
                                     5: Point2D((point1[2]-self.zmin-deltaz/2, point1[1]-self.ymin-deltay/2)),
                                     6: Point2D((point1[1]-self.ymin-deltay/2, point1[0]-self.xmin-deltax/2))}
        
        point2_2d_coordinate_dict = {1: Point2D((point2[0]-self.xmin-deltax/2, point2[1]-self.ymin-deltay/2)),
                                     2: Point2D((point2[2]-self.zmin-deltaz/2, point2[0]-self.xmin-deltax/2)),
                                     3: Point2D((point2[1]-self.ymin-deltay/2, point2[2]-self.zmin-deltaz/2)),
                                     4: Point2D((point2[0]-self.xmin-deltax/2, point2[2]-self.zmin-deltaz/2)),
                                     5: Point2D((point2[2]-self.zmin-deltaz/2, point2[1]-self.ymin-deltay/2)),
                                     6: Point2D((point2[1]-self.ymin-deltay/2, point2[0]-self.xmin-deltax/2))}
        
        vertex_2d_coordinate_dict = {1: [Point2D((self.xmin-self.xmin-deltax/2, self.ymin-self.ymin-deltay/2)), Point2D((self.xmin-self.xmin-deltax/2, self.ymax-self.ymin-deltay/2)), Point2D((self.xmax-self.xmin-deltax/2, self.ymax-self.ymin-deltay/2)), Point2D((self.xmax-self.xmin-deltax/2, self.ymin-self.ymin-deltay/2))],
                                     2: [Point2D((self.zmin-self.zmin-deltaz/2, self.xmin-self.xmin-deltax/2)), Point2D((self.zmin-self.zmin-deltaz/2, self.xmax-self.xmin-deltax/2)), Point2D((self.zmax-self.zmin-deltaz/2, self.xmax-self.xmin-deltax/2)), Point2D((self.zmax-self.zmin-deltaz/2, self.xmin-self.xmin-deltax/2))],
                                     3: [Point2D((self.ymin-self.ymin-deltay/2, self.zmin-self.zmin-deltaz/2)), Point2D((self.ymin-self.ymin-deltay/2, self.zmax-self.zmin-deltaz/2)), Point2D((self.ymax-self.ymin-deltay/2, self.zmax-self.zmin-deltaz/2)), Point2D((self.ymax-self.ymin-deltay/2, self.zmin-self.zmin-deltaz/2))],
                                     4: [Point2D((self.xmin-self.xmin-deltax/2, self.zmin-self.zmin-deltaz/2)), Point2D((self.xmin-self.xmin-deltax/2, self.zmax-self.zmin-deltaz/2)), Point2D((self.xmax-self.xmin-deltax/2, self.zmax-self.zmin-deltaz/2)), Point2D((self.xmax-self.xmin-deltax/2, self.zmin-self.zmin-deltaz/2))],
                                     5: [Point2D((self.zmin-self.zmin-deltaz/2, self.ymin-self.ymin-deltay/2)), Point2D((self.zmin-self.zmin-deltaz/2, self.ymax-self.ymin-deltay/2)), Point2D((self.zmax-self.zmin-deltaz/2, self.ymax-self.ymin-deltay/2)), Point2D((self.zmax-self.zmin-deltaz/2, self.ymin-self.ymin-deltay/2))],
                                     6: [Point2D((self.ymin-self.ymin-deltay/2, self.xmin-self.xmin-deltax/2)), Point2D((self.ymin-self.ymin-deltay/2, self.xmax-self.xmin-deltax/2)), Point2D((self.ymax-self.ymin-deltay/2, self.xmax-self.xmin-deltax/2)), Point2D((self.ymax-self.ymin-deltay/2, self.xmin-self.xmin-deltax/2))],}
        
        vertex_2d_coordinate_dict2 = {1: [Point2D((self.xmin, self.ymin)), Point2D((self.xmin, self.ymax)), Point2D((self.xmax, self.ymax)), Point2D((self.xmax, self.ymin))],
                                    2: [Point2D((self.zmin, self.xmin)), Point2D((self.zmin, self.xmax)), Point2D((self.zmax, self.xmax)), Point2D((self.zmax, self.xmin))],
                                     3: [Point2D((self.ymin, self.zmin)), Point2D((self.ymin, self.zmax)), Point2D((self.ymax, self.zmax)), Point2D((self.ymax, self.zmin))],
                                     4: [Point2D((self.xmin, self.zmin)), Point2D((self.xmin, self.zmax)), Point2D((self.xmax, self.zmax)), Point2D((self.xmax, self.zmin))],
                                     5: [Point2D((self.zmin, self.ymin)), Point2D((self.zmin, self.ymax)), Point2D((self.zmax, self.ymax)), Point2D((self.zmax, self.ymin))],
                                     6: [Point2D((self.ymin, self.xmin)), Point2D((self.ymin, self.xmax)), Point2D((self.ymax, self.xmax)), Point2D((self.ymax, self.xmin))],}
        
        
        vertex_to_3d_dict = {1: (2, self.zmax, 0, 1),
                             2: (1, self.ymax, 2, 0),
                             3: (0, self.xmax, 1, 2),
                             4: (1, self.ymin, 0, 2),
                             5: (0, self.xmin, 2, 1),
                             6: (2, self.zmin, 1, 0)}
        
        offset_dict = {0: self.xmin+deltax/2,
                       1: self.ymin+deltay/2,
                       2: self.zmin+deltaz/2}
        
        opposite_face_dict = {1: 6, 2: 4, 3: 5, 4: 2, 5: 3, 6: 1}
        
        combination_dict = {(1, 2): Frame2D(Point2D((0,  deltay/2+deltaz/2)), Vector2D((0,-1)), Vector2D(( 1, 0))),
                            (2, 1): Frame2D(Point2D(( deltay/2+deltaz/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (1, 3): Frame2D(Point2D(( deltax/2+deltaz/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (3, 1): Frame2D(Point2D((0,  deltax/2+deltaz/2)), Vector2D((0,-1)), Vector2D(( 1, 0))),
                            (1, 4): Frame2D(Point2D((0, -deltay/2-deltaz/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (4, 1): Frame2D(Point2D((-deltay/2-deltaz/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (1, 5): Frame2D(Point2D((-deltax/2-deltaz/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (5, 1): Frame2D(Point2D((0, -deltax/2-deltaz/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (2, 3): Frame2D(Point2D((0,  deltax/2+deltay/2)), Vector2D((0,-1)), Vector2D(( 1, 0))),
                            (3, 2): Frame2D(Point2D(( deltax/2+deltay/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (2, 5): Frame2D(Point2D((0, -deltax/2-deltay/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (5, 2): Frame2D(Point2D((-deltax/2-deltay/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (2, 6): Frame2D(Point2D((-deltaz/2-deltay/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (6, 2): Frame2D(Point2D((0, -deltaz/2-deltay/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (3, 4): Frame2D(Point2D((-deltay/2-deltax/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (4, 3): Frame2D(Point2D((0, -deltay/2-deltax/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (3, 6): Frame2D(Point2D((0, -deltaz/2-deltax/2)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (6, 3): Frame2D(Point2D((-deltaz/2-deltax/2, 0)), Vector2D((1, 0)), Vector2D(( 0, 1))),
                            (4, 5): Frame2D(Point2D((-deltax/2-deltay/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (5, 4): Frame2D(Point2D((0, -deltax/2-deltay/2)), Vector2D((0,-1)), Vector2D(( 1, 0))),
                            (4, 6): Frame2D(Point2D((0, -deltaz/2-deltay/2)), Vector2D((0,-1)), Vector2D(( 1, 0))),
                            (6, 4): Frame2D(Point2D((-deltaz/2-deltay/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (5, 6): Frame2D(Point2D((-deltaz/2-deltax/2, 0)), Vector2D((0, 1)), Vector2D((-1, 0))),
                            (6, 5): Frame2D(Point2D((0, -deltaz/2-deltax/2)), Vector2D((0,-1)), Vector2D(( 1, 0)))}
        
        point1_2d = point1_2d_coordinate_dict[face_point1]
        point2_2d = point2_2d_coordinate_dict[face_point2]
        
        # The points are on adjacent faces
        if opposite_face_dict[face_point1] != face_point2:
            frame =  combination_dict[(face_point1, face_point2)]
            net_point2 = frame.OldCoordinates(point2_2d)
            
            # Computes the 3D intersection between the net_line and the edges of the face_point1
            net_line = LineSegment2D(point1_2d, net_point2)
            vertex_points = vertex_2d_coordinate_dict[face_point1]
            edge_lines = [LineSegment2D(p1, p2) for p1, p2 in zip(vertex_points, vertex_points[1:]+[vertex_points[0]])]
            for line in edge_lines:
                edge_intersection_point, a, b = Point2D.LinesIntersection(net_line, line, curvilinear_abscissa=True)
                if edge_intersection_point is not None \
                and a > 0 and a < 1 and b > 0 and b < 1:
                    break
            offset_indice, offset, indice1, indice2 = vertex_to_3d_dict[face_point1]
            disordered_coordinate = [(indice1, edge_intersection_point[0]+offset_dict[indice1]),
                                     (indice2, edge_intersection_point[1]+offset_dict[indice2]),
                                     (offset_indice, offset)]
            disordered_coordinate = sorted(disordered_coordinate, key=lambda a: a[0])
            intersection_point_3d = Point3D(tuple([p[1] for p in disordered_coordinate]))
            
            mesures = [Mesure(point1_copy, intersection_point_3d), Mesure(intersection_point_3d, point2_copy)]
            
            return mesures
        
        # The points are on opposite faces
        else:
            net_points2_and_frame = []
            
            faces_number = [1, 2, 3, 4, 5, 6]
            faces_number.remove(face_point1)
            faces_number.remove(face_point2)
            pathes = []
            for face_nb in faces_number:
                path = [(face_point1, face_nb), (face_nb, face_point2)]
                pathes.append(path)

            for path in pathes:
                frame1 = combination_dict[(path[0][0], path[0][1])]
                frame2 = combination_dict[(path[1][0], path[1][1])]
                frame = frame1 + frame2
                net_points2_and_frame.append((Point2D(frame.OldCoordinates(point2_2d).vector), frame))
            net_point2, frame = min(net_points2_and_frame, key=lambda pt: pt[0].PointDistance(point1_2d))
            net_line = LineSegment2D(point1_2d, net_point2)

            # Computes the 3D intersection between the net_line and the edges of the face_point1
            vertex_points = vertex_2d_coordinate_dict[face_point1]
            edge_lines = [LineSegment2D(p1, p2) for p1, p2 in zip(vertex_points, vertex_points[1:]+[vertex_points[0]])]
            for line in edge_lines:
                edge_intersection_point1, a, b = Point2D.LinesIntersection(net_line, line, curvilinear_abscissa=True)                
                if edge_intersection_point1 is not None \
                and a > 0 and a < 1 and b > 0 and b < 1:
                    break
            offset_indice, offset, indice1, indice2 = vertex_to_3d_dict[face_point1]
            disordered_coordinate = [(indice1, edge_intersection_point1[0]+offset_dict[indice1]), 
                                     (indice2, edge_intersection_point1[1]+offset_dict[indice2]),
                                     (offset_indice, offset)]
            disordered_coordinate = sorted(disordered_coordinate, key=lambda a: a[0])
            intersection_point1_3d = Point3D(tuple([p[1] for p in disordered_coordinate]))

            # Computes the 3D intersection between the net_line and the edges of the face_point2            
            vertex_points = [frame.OldCoordinates(p) for p in vertex_2d_coordinate_dict[face_point2]]
            edge_lines = [LineSegment2D(p1, p2) for p1, p2 in zip(vertex_points, vertex_points[1:]+[vertex_points[0]])]
            for line in edge_lines:
                edge_intersection_point2, a, b = Point2D.LinesIntersection(net_line, line, curvilinear_abscissa=True)
                if edge_intersection_point2 is not None \
                and a > 0 and a < 1 and b > 0 and b < 1:
                    break
            edge_intersection_point2 = Point2D(frame.NewCoordinates(edge_intersection_point2))
            offset_indice, offset, indice1, indice2 = vertex_to_3d_dict[face_point2]
            disordered_coordinate = [(indice1, edge_intersection_point2[0]+offset_dict[indice1]),
                                     (indice2, edge_intersection_point2[1]+offset_dict[indice2]),
                                     (offset_indice, offset)]
            disordered_coordinate = sorted(disordered_coordinate, key=lambda a: a[0])
            intersection_point2_3d = Point3D(tuple([p[1] for p in disordered_coordinate]))
            
            if point1 == point1_copy:
                mesures = [Mesure(point1, intersection_point1_3d), Mesure(intersection_point1_3d, intersection_point2_3d), Mesure(intersection_point2_3d, point2)]
            else:
                mesures = [Mesure(point2, intersection_point2_3d), Mesure(intersection_point2_3d, intersection_point1_3d), Mesure(intersection_point1_3d, point1)]
            return mesures

    def Babylon(self):
        height = self.ymax-self.ymin
        width = self.xmax-self.xmin
        depth = self.zmax-self.zmin
        s = 'var box = BABYLON.MeshBuilder.CreateBox("box", {{height: {}, width: {}, depth: {}}}, scene);\n'.format(height, width, depth)
        s += 'box.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n'.format(self.center[0], self.center[1], self.center[2])
        s += 'var bboxmat = new BABYLON.StandardMaterial("bboxmat", scene);\n'
        s += 'bboxmat.alpha = 0.4;\n'
        s += 'var DTWidth = {};\n'.format(width*60)
        s += 'var DTHeight = {};\n'.format(height*60)
        s += 'var font_type = "Arial";\n'
        s += 'var text = "{}";\n'.format(self.name)
        s += 'var dynamicTexture = new BABYLON.DynamicTexture("DynamicTexture", {width:DTWidth, height:DTHeight}, scene);\n'
        s += 'var ctx = dynamicTexture.getContext();\n'
        s += 'var size = 0.8;\n'
        s += 'ctx.font = size + "px " + font_type;\n'
        s += 'var textWidth = ctx.measureText(text).width;\n'
        s += 'var ratio = textWidth/size;\n'
        s += 'var font_size = Math.floor(DTWidth / ratio);\n'
        s += 'var font = font_size + "px " + font_type;\n'
        s += 'dynamicTexture.drawText(text, null, null, font, "#000000", "#ffffff", false);\n'
        s += 'bboxmat.diffuseTexture = dynamicTexture;\n'
        s += 'box.material = bboxmat;\n'
        return s


class Mesure(Line3D):
    def __init__(self, point1, point2, color=(1,0,0)):
        self.points = [point1, point2]
        self.color = color
        self.distance = Vector3D(self.points[0]-self.points[1]).Norm()
        self.bounding_box = self._bounding_box()
        
    def Babylon(self):
        s = 'var myPoints = [];\n'
        s += 'var point1 = new BABYLON.Vector3({},{},{});\n'.format(self.points[0][0],self.points[0][1],self.points[0][2])
        s += 'myPoints.push(point1);\n'
        s += 'var point2 = new BABYLON.Vector3({},{},{});\n'.format(self.points[1][0],self.points[1][1],self.points[1][2])
        s += 'myPoints.push(point2);\n'
        s += 'var line = BABYLON.MeshBuilder.CreateLines("lines", {points: myPoints}, scene);\n'
        s += 'line.color = new BABYLON.Color3({}, {}, {});\n'.format(self.color[0], self.color[1], self.color[2])
        return s 
    

class Group:
    def __init__(self, primitives, name):
        self.primitives = primitives
        self.name = name


class StepFunction:
    def __init__(self, function_id, function_name, function_arg):
        self.id = function_id
        self.name = function_name
        self.arg = function_arg

        if self.name == "":
            if self.arg[1][0] == 'B_SPLINE_SURFACE':
                self.simplify('B_SPLINE_SURFACE')
            if self.arg[1][0] == 'B_SPLINE_CURVE':
                self.simplify('B_SPLINE_CURVE')

    def simplify(self, new_name):
        # ITERATE ON SUBFUNCTIONS
        args = [subfun[1] for (i, subfun) in enumerate(self.arg) if (len(subfun[1]) != 0 or i == 0)]
        arguments = []
        for arg in args:
            if arg == []:
                arguments.append('')
            else:
                arguments.extend(arg)
        arguments.pop() # DELETE REPRESENTATION_ITEM('')

        self.name = new_name
        self.arg = arguments
        

class Step:
    def __init__(self, stepfile):
        self.stepfile = stepfile

        self.functions, self.all_connections = self.read_functions()
        self.graph = self.create_graph()

    def read_functions(self):
        f = open(self.stepfile, "r", encoding = "ISO-8859-1")

        all_connections = []

        previous_line = ""
        functions = {}

        for line in f:

            line = line.replace(" ", "")
            line = line.replace("\n", "")

            # SKIP EMPTY LINE
            if not line:
                continue

            # ASSEMBLE LINES IF THEY ARE SEPARATED
            if line[-1] != ';':
                previous_line = previous_line + line
                continue

            line = previous_line + line

            # SKIP HEADER
            if line[0] != "#":
                previous_line = str()
                continue

            function = line.split("=")
            function_id = int(function[0][1:])
            function_name_arg = function[1].split("(", 1)
            function_name = function_name_arg[0]
            function_arg = function_name_arg[1].split("#")
            function_connections = []
            for connec in function_arg[1:]:
                connec = connec.split(",")
                connec = connec[0].split(")")
                function_connection = int(connec[0])
                function_connections.append((function_id, function_connection))
                
            all_connections.extend(function_connections)

            previous_line = str()

            # FUNCTION ARGUMENTS
            function_arg = function_name_arg[1]
            arguments = step_split_arguments(function_arg)
            if function_name == "":
                arguments = self.step_subfunctions(arguments)

            for i, argument in enumerate(arguments):
                if argument[:2] == '(#' and argument[-1] == ')':
                    arg_list = set_to_list(argument)
                    arguments[i] = arg_list

            function = StepFunction(function_id, function_name, arguments)
            functions[function_id] = function

        f.close()

        return functions, all_connections

    def create_graph(self, draw=False, html=False):

        G = nx.Graph()
        F = nx.DiGraph()
        labels = {}

        for function in self.functions.values():
            if function.name in step_to_volmdlr_primitive:
                G.add_node(function.id)
                F.add_node(function.id)
                labels[function.id] = str(function.id)+' '+function.name

        # Delete connection if node not found
        node_list = list(F.nodes())
        delete_connection = []
        for connection in self.all_connections:
            if connection[0] not in node_list or connection[1] not in node_list:
                delete_connection.append(connection)
        for delete in delete_connection:
            self.all_connections.remove(delete)

        # Create graph connections
        G.add_edges_from(self.all_connections)
        F.add_edges_from(self.all_connections)

        # Remove single nodes
        delete_nodes = []
        for node in F.nodes:
            if F.degree(node) == 0:
                delete_nodes.append(node)
        for node in delete_nodes:
            F.remove_node(node)
            G.remove_node(node)

        if draw:
            # ----------------PLOT----------------
            pos = nx.kamada_kawai_layout(G)
            plt.figure()
            nx.draw_networkx_nodes(F, pos)
            nx.draw_networkx_edges(F, pos)
            nx.draw_networkx_labels(F, pos, labels)
            # ------------------------------------

        if html:

            env = Environment(loader=PackageLoader('powertransmission', 'templates'),
                              autoescape=select_autoescape(['html', 'xml']))
            template = env.get_template('graph_visJS.html')

            nodes = []
            edges = []
            for label in list(labels.values()):
                nodes.append({'name': label, 'shape': 'circular'})

            for edge in G.edges:
                edge_dict = {}
                edge_dict['inode1'] = int(edge[0])-1
                edge_dict['inode2'] = int(edge[1])-1
                edges.append(edge_dict)

            options = {}
            s = template.render(
                name=self.stepfile,
                nodes=nodes,
                edges=edges,
                options=options)

            with open('graph_visJS.html', 'wb') as file:
                file.write(s.encode('utf-8'))

            webbrowser.open('file://' + os.path.realpath('graph_visJS.html'))

        return F

    def draw_graph(self):
        labels = {}
        for id_nb, function in self.functions.items():
            labels[id_nb] = str(id_nb)+' '+function.name
        print(labels)
        print()
        pos = nx.kamada_kawai_layout(self.graph)
        print(pos)
        plt.figure()
        nx.draw_networkx_nodes(self.graph, pos)
        nx.draw_networkx_edges(self.graph, pos)
        nx.draw_networkx_labels(self.graph, pos, labels)

    def step_subfunctions(self, subfunctions):
        subfunctions = subfunctions[0]
        parenthesis_count = 0
        subfunction_names = []
        subfunction_args = []
        subfunction_name = ""
        subfunction_arg = ""
        for char in subfunctions:

            if char == "(":
                parenthesis_count += 1
                if parenthesis_count == 1:
                    subfunction_names.append(subfunction_name)
                    subfunction_name = ""
                else:
                    subfunction_arg += char
                continue
            if char == ")":
                parenthesis_count -= 1
                if parenthesis_count == 0:
                    subfunction_args.append(subfunction_arg)
                    subfunction_arg = ""
                else:
                    subfunction_arg += char
                continue

                subfunction_name += char
            else:
                subfunction_arg += char
        return [(subfunction_names[i], step_split_arguments(subfunction_args[i])) for i in range(len(subfunction_names))]

    def instanciate(self, instanciate_id, object_dict, primitives):
        """
        Returns None if the object was instanciate
        """

        name = self.functions[instanciate_id].name
        arguments = self.functions[instanciate_id].arg[:]

        for i, arg in enumerate(arguments):
            if type(arg) == str and arg[0] == '#':
                arguments[i] = int(arg[1:])
            
        if name == 'VERTEX_POINT':
            object_dict[instanciate_id] = object_dict[arguments[1]]
            
        elif name == 'LINE':
            pass
            
        elif name == 'ORIENTED_EDGE':
            object_dict[instanciate_id] = object_dict[arguments[3]]
            
        elif name == 'FACE_OUTER_BOUND':
            object_dict[instanciate_id] = object_dict[arguments[1]]
            
        elif name == 'FACE_BOUND':
            object_dict[instanciate_id] = object_dict[arguments[1]]
            
        elif name == 'SURFACE_CURVE':
            pass
        
        
        elif name in step_to_volmdlr_primitive and hasattr(step_to_volmdlr_primitive[name], "from_step"):
            volmdlr_object = step_to_volmdlr_primitive[name].from_step(arguments, object_dict)

            object_dict[instanciate_id] = volmdlr_object
            if hasattr(volmdlr_object, "primitive"):
                primitives.append(volmdlr_object.primitive)
        
        return None, object_dict, primitives
    
    def to_shells3d(self, name):
        object_dict = {}
        primitives = []
                                    
        self.graph.add_node("#0")
        for node in self.graph.nodes:
            if node != '#0' and self.functions[node].name == "CLOSED_SHELL":
                self.graph.add_edge("#0", node)
        
        edges = list(nx.algorithms.traversal.breadth_first_search.bfs_edges(self.graph, "#0"))[::-1]        

        for edge_nb, edge in enumerate(edges):
            instanciate_id = edge[1]
            res, object_dict, primitives = self.instanciate(instanciate_id, object_dict, primitives)
            if res is not None:
                raise NotImplementedError

        shells = []
        for node in list(self.graph.nodes):
            if node != '#0' and self.functions[node].name == 'CLOSED_SHELL':
                shells.append(object_dict[node])
        print(shells)
        return shells
        

class VolumeModel:
    """
    :param groups: A list of two element tuple. The first element is a string naming the group and the second element is a list of primitives of the group
    """
    def __init__(self, primitives, name=''):
        self.primitives = primitives
        self.name = name
        if self.primitives:
            self.shells = self._extract_shells()
        if self.shells:
            self.bounding_box = self._bounding_box()
            
    def _extract_shells(self):
        shells = []
        for primitive in self.primitives:
            if isinstance(primitive, Shell3D):
                shells.append(primitive)
        return shells

    def Volume(self):
        volume=0
        for primitive in self.primitives:
            volume+=primitive.Volume()
        return volume

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_primitives = [primitive.Rotation(center, axis, angle, copy=True) for primitive in self.primitives]
            return VolumeModel(new_primitives, self.name)
        else:
            for primitives in self.primitives:
                primitives.Translation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def Translation(self, offset, copy=True):
        if copy:
            new_primitives = [primitive.Translation(offset, copy=True) for primitive in self.primitives]
            return VolumeModel(new_primitives, self.name)
        else:
            for primitives in self.primitives:
                primitives.Translation(offset, copy=False)
            self.bounding_box = self._bounding_box()
            
        
    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_primitives = [primitive.frame_mapping(frame, side, copy=True) for primitive in self.primitives]
            return VolumeModel(new_primitives, self.name)
        else:
            for primitives in self.primitives:
                primitives.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()
            
    def copy(self):
        new_primitives = [primitive.copy() for primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)
            
    def _bounding_box(self):
        bboxes = []
        for primitive in self.primitives:
            if hasattr(primitive, 'bounding_box'):
                bboxes.append(primitive.bounding_box)

        xmin = min([box.xmin for box in bboxes])
        xmax = max([box.xmax for box in bboxes])
        ymin = min([box.ymin for box in bboxes])
        ymax = max([box.ymax for box in bboxes])
        zmin = min([box.zmin for box in bboxes])
        zmax = max([box.zmax for box in bboxes])
        
        return BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)
            
    def plot(self, ax=None, color=None):
        fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')
        
        for i, shell in enumerate(self.shells):
            bbox = shell.bbox()
            bbox.plot(ax, color[i])
            
        return ax

    def MPLPlot(self):
        """
        Matplotlib plot of model.
        To use for debug.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', adjustable='box')
#        ax.set_aspect('equal')
        for primitive in self.primitives:
            primitive.MPLPlot(ax)
        ax.set_aspect('equal')
        return fig, ax

    def FreeCADScript(self, fcstd_filepath,
                      freecad_lib_path='/usr/lib/freecad/lib',
                      export_types=('fcstd',),
                      save_to = '',
                      tolerance=0.0001):
        """
        Generate python a FreeCAD definition of model
        :param fcstd_filename: a filename without extension to give the name at the fcstd part written in python code
        :type fcstd_filename:str
        """
        fcstd_filepath = os.path.abspath(fcstd_filepath)
        fcstd_filepath = fcstd_filepath.replace('\\','\\\\')
        freecad_lib_path = freecad_lib_path.replace('\\','\\\\')

        s='# -*- coding: utf-8 -*-\n'
        if freecad_lib_path != '':
            s+="import sys\nsys.path.append('"+freecad_lib_path+"')\n"

        s+="import math\nimport FreeCAD as fc\nimport Part\n\ndoc=fc.newDocument('doc')\n\n"
        for ip, primitive in enumerate(self.primitives):
            if primitive.name == '':
                primitive_name = 'Primitive_{}'.format(ip)
            else:
                primitive_name = 'Primitive_{}_{}'.format(ip, primitive_name)
            s += "part = doc.addObject('App::Part','{}')\n".format(primitive_name)
            if hasattr(primitive, 'FreeCADExport'):
                sp = primitive.FreeCADExport(ip)
                if sp != '':
#                        s += (sp+'\n')
                    s += (sp)
                    s += 'shapeobj = doc.addObject("Part::Feature","{}")\n'.format(primitive_name)
                    if isinstance(primitive, BSplineCurve3D) \
                    or isinstance(primitive, BSplineSurface3D) \
                    or isinstance(primitive, Circle3D) \
                    or isinstance(primitive, Ellipse3D):
#                            print(primitive)
#                            s += 'S = Part.Shape([primitive{}])\n'.format(ip)
#                            s += 'shapeobj.Shape = S\n'
                        s += 'shapeobj.Shape = primitive{}.toShape()\n'.format(ip)
                    else:
                        s += "shapeobj.Shape = primitive{}\n".format(ip)
                    s += 'part.addObject(shapeobj)\n\n'.format(ip, primitive.name)
            # --------------------DEBUG-------------------
#                else:
#                    raise NotImplementedError
            # ---------------------------------------------

        s+='doc.recompute()\n'
        if 'fcstd' in export_types:
            s+="doc.saveAs('"+fcstd_filepath+".fcstd')\n\n"
        if 'stl' in export_types:
            s+="import Mesh\nMesh.export(doc.Objects,'{}.stl', tolerance={})\n".format(fcstd_filepath, tolerance)
        if 'step' in export_types:
            s+="Part.export(doc.Objects,'{}.step')\n".format(fcstd_filepath)


        if save_to != '':
            with open(os.path.abspath(save_to),'w') as file:
                file.write(s)
        return s

    def FreeCADExport(self,fcstd_filepath,
                      python_path='python',
                      freecad_lib_path='/usr/lib/freecad/lib',
                      export_types=('fcstd',),
                      tolerance=0.0001):
        """
        Export model to .fcstd FreeCAD standard

        :param python_path: path of python binded to freecad

            * on windows: something like C:\\\\Program Files\\\\FreeCAD X.XX\\\\bin\\\\python
            * on linux: python if installed by a dstribution package
        :param filepath: path of fcstd file (without extension)
        :param freecad_lib_path: FreeCAD.so lib path (/usr/lib/freecad/lib in general)
        :param tolerance: the tolerance of tesselation for mesh exports

        """
        fcstd_filepath=os.path.abspath(fcstd_filepath)
        s=self.FreeCADScript(fcstd_filepath,
                             freecad_lib_path = freecad_lib_path,
                             export_types = export_types,
                             tolerance=tolerance)
        with tempfile.NamedTemporaryFile(suffix=".py",delete=False) as f:
            f.write(bytes(s,'utf8'))

        arg=f.name
        output=subprocess.call([python_path, arg])

        f.close()
        os.remove(f.name)
        return output

    def BabylonScript(self, use_cdn=True):

        env = Environment(loader=PackageLoader('volmdlr', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))

        template = env.get_template('babylon.html')
        
        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])

        primitives_strings=[]
        for primitive in self.primitives:
            if hasattr(primitive, 'Babylon'):
                primitives_strings.append(primitive.Babylon())
        return template.render(name=self.name,
                               center=tuple(center),
                               length=2*max_length,
                               primitives_strings=primitives_strings,
                               use_cdn=use_cdn)

    def BabylonShow(self,page='vm_babylonjs', use_cdn=True):
        page+='.html'
        with open(page,'w') as file:
            file.write(self.BabylonScript(use_cdn=use_cdn))

        webbrowser.open('file://' + os.path.realpath(page))
        
        

class Routing:
    def __init__(self, point1, point2, volumemodel):
        self.points = [point1, point2]
        self.volumemodel = volumemodel
        
    def straight_line(self):
        """
        Returns 2 distances : 
            - no collision distance
            - collision distance
        """
        line = LineSegment3D(self.points[0], self.points[1])
        
        intersection_points = []
        for shell in self.volumemodel.shells:
            for face in shell.faces:
                intersection_point, intersection_abscissea = face.linesegment_intersection(line, abscissea=True)
                if intersection_point is not None and intersection_abscissea != 0 and intersection_abscissea != 1:
                    intersection_points.append((intersection_point, intersection_abscissea))
        
        if len(intersection_points)%2 != 0:
            raise NotImplementedError

        intersection_points = sorted(intersection_points, key=lambda abscissea: abscissea[1])
        all_points_abscissea = [(self.points[0], 0)] + intersection_points[:] + [(self.points[1], 1)]
        all_points = [p[0] for p in all_points_abscissea]
        
        no_collision_mesures = []
        collision_mesures = []
        i = 0
        for pt1, pt2 in zip(all_points[:-1], all_points[1:]):
            if i%2 == 0:
                no_collision_mesures.append(Mesure(pt1, pt2, color=(0,0,1)))
            else:
                collision_mesures.append(Mesure(pt1, pt2, color=(1,0,0)))
            i += 1
            
        return no_collision_mesures, collision_mesures
    
    def straight_line2(self):
        """
        Returns the distance of the line going around each shell's bbox encountered along the path
        """
        line = LineSegment3D(self.points[0], self.points[1])
        
        all_mesures_abscissea = []
        intersection_points = []
        for shell in self.volumemodel.shells:
            shell_intersection_points = []
            bbox = shell.bounding_box
            for face in shell.faces:
                intersection_point, intersection_abscissea = face.linesegment_intersection(line, abscissea=True)
                if intersection_point is not None and intersection_abscissea != 0 and intersection_abscissea != 1:
                    intersection_points.append((intersection_point, intersection_abscissea))
                    shell_intersection_points.append((intersection_point, intersection_abscissea))
                    
            if len(shell_intersection_points) == 2:
                shell_intersection_points = sorted(shell_intersection_points, key=lambda abscissea: abscissea[1])
                abscissea1 = shell_intersection_points[0][1]
                abscissea2 = shell_intersection_points[1][1]
                shell_intersection_points = [p[0] for p in shell_intersection_points]
                around_bbox_mesures = bbox.distance_between_two_points_on_bbox(shell_intersection_points[0], shell_intersection_points[1])
                all_mesures_abscissea.append((around_bbox_mesures, abscissea1, abscissea2))
            elif len(shell_intersection_points) > 2 or len(shell_intersection_points) == 1:
                raise NotImplementedError
        
        intersection_points = sorted(intersection_points, key=lambda abscissea: abscissea[1])
        all_mesures_abscissea = sorted(all_mesures_abscissea, key=lambda abscissea: abscissea[1])
        all_points_abscissea = [(self.points[0], 0)] + intersection_points[:] + [(self.points[1], 1)]
        all_points = [p[0] for p in all_points_abscissea]
                    
        no_collision_mesures = []
        i = 0
        for pt1, pt2 in zip(all_points[:-1], all_points[1:]):
            if i%2 == 0:
                no_collision_mesures.append(Mesure(pt1, pt2, color=(0,0,1)))
            else:
                no_collision_mesures.extend(all_mesures_abscissea[i//2][0])
            i += 1
            
        return no_collision_mesures


class ViewIso:
    def __init__(self, component, frame, size):
        self.component = component
        self.frame = frame
        self.size = size
        self.plot_datas = self.plot_data()
        
    def plot_data(self, detail=True):
        wide = min(self.size)/2
        plot_datas = []
        plot_datas.extend(self.component.plot_data(self.frame, detail=detail))
        plot_datas.extend(self.component.plot_data(Frame3D(self.frame.origin + Point3D((0, self.size[1]/2 + self.size[2]/2 + wide, 0)), self.frame.u, self.frame.w, self.frame.v), detail = detail))   
        plot_datas.extend(self.component.plot_data(Frame3D(self.frame.origin + Point3D((self.size[0]/2 + self.size[2]/2 + wide, 0, 0)), self.frame.w, self.frame.v, self.frame.u), detail = detail))
        return plot_datas
        
    def plot(self):
        plot_data.plot(self.plot_datas)



step_to_volmdlr_primitive = {
        # GEOMETRICAL ENTITIES
        'CARTESIAN_POINT': Point3D,
        'DIRECTION': Vector3D,
        'VECTOR': Vector3D,

        'AXIS1_PLACEMENT': None,
        'AXIS2_PLACEMENT_2D': None,
        'AXIS2_PLACEMENT_3D': Frame3D,

        'LINE': None,
        'CIRCLE': Circle3D,
        'ELLIPSE': Ellipse3D,
        'PARABOLA': None,
        'HYPERBOLA': None,
        'PCURVE': None,
        'CURVE_REPLICA': None,
        'OFFSET_CURVE_3D': None,
        'TRIMMED_CURVE': None, # BSplineCurve3D cannot be trimmed on FreeCAD
        'B_SPLINE_CURVE': BSplineCurve3D,
        'B_SPLINE_CURVE_WITH_KNOTS': BSplineCurve3D,
        'BEZIER_CURVE': BSplineCurve3D,
        'RATIONAL_B_SPLINE_CURVE': BSplineCurve3D,
        'UNIFORM_CURVE': BSplineCurve3D,
        'QUASI_UNIFORM_CURVE': BSplineCurve3D,
        'SURFACE_CURVE': None, # TOPOLOGICAL EDGE
        'SEAM_CURVE': LineSegment3D, # TOPOLOGICAL EDGE
        'COMPOSITE_CURVE_SEGMENT': None, # TOPOLOGICAL EDGE
        'COMPOSITE_CURVE': Wire3D, # TOPOLOGICAL WIRE
        'COMPOSITE_CURVE_ON_SURFACE': Wire3D, # TOPOLOGICAL WIRE
        'BOUNDARY_CURVE': Wire3D, # TOPOLOGICAL WIRE

        'PLANE': Plane3D,
        'CYLINDRICAL_SURFACE': None,
        'CONICAL_SURFACE': None,
        'SPHERICAL_SURFACE': None,
        'TOROIDAL_SURFACE': None,
        'DEGENERATE_TOROIDAL_SURFACE': None,
        'B_SPLINE_SURFACE_WITH_KNOTS': BSplineSurface3D,
        'B_SPLINE_SURFACE': BSplineSurface3D,
        'BEZIER_SURFACE': BSplineSurface3D,
        'OFFSET_SURFACE': None,
        'SURFACE_REPLICA': None,
        'RATIONAL_B_SPLINE_SURFACE': BSplineSurface3D,
        'RECTANGULAR_TRIMMED_SURFACE': None,
        'SURFACE_OF_LINEAR_EXTRUSION': None, # CAN BE A BSplineSurface3D
        'SURFACE_OF_REVOLUTION': None,
        'UNIFORM_SURFACE': BSplineSurface3D,
        'QUASI_UNIFORM_SURFACE': BSplineSurface3D,
        'RECTANGULAR_COMPOSITE_SURFACE': Face3D, # TOPOLOGICAL FACES
        'CURVE_BOUNDED_SURFACE': Face3D, # TOPOLOGICAL FACE


        # TOPOLOGICAL ENTITIES
        'VERTEX_POINT': None,

        'EDGE_CURVE': LineSegment3D, # TOPOLOGICAL EDGE
        'ORIENTED_EDGE': None, # TOPOLOGICAL EDGE
        # The one above can influence the direction with their last argument
        # TODO : maybe take them into consideration

        'FACE_BOUND': None, # TOPOLOGICAL WIRE
        'FACE_OUTER_BOUND': None, # TOPOLOGICAL WIRE
        # Both above can influence the direction with their last argument
        # TODO : maybe take them into consideration
        'EDGE_LOOP': Contour3D, # TOPOLOGICAL WIRE
        'POLY_LOOP': Contour3D, # TOPOLOGICAL WIRE
        'VERTEX_LOOP': Contour3D, # TOPOLOGICAL WIRE

        'ADVANCED_FACE': Face3D,
        'FACE_SURFACE': Face3D,

        'CLOSED_SHELL': Shell3D,
        'OPEN_SHELL': Shell3D,
#        'ORIENTED_CLOSED_SHELL': None,
        'CONNECTED_FACE_SET': Shell3D,

        }
