#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module containing mesh and relative objects
"""

from typing import List  # TypeVar, Tuple, Dict
from itertools import combinations
import math
import matplotlib.pyplot as plt
import numpy as npy
from dessia_common.core import DessiaObject
# import volmdlr.core_compiled
import volmdlr as vm
import volmdlr.wires as vmw
import volmdlr.edges as vme
# from volmdlr.core_compiled import Matrix33

# from itertools import combinations
# import numpy as npy
# import volmdlr.wires
# import volmdlr.faces
# from volmdlr.core_compiled import Matrix33
# import matplotlib
# import random
# from itertools import product
# from matplotlib.colors import LinearSegmentedColormap

# cdict = {'red':  [(0.0, 0.0, 0.0),
#                    (1.0, 1.0, 1.0)],
#          'green': [(0.0, 0.0, 0.0),
#                    (1.0, 0.0, 0.0)],
#          'blue':  [(0.0, 1.0, 1.0),
#                    (1.0, 0.0, 0.0)]}
# blue_red = LinearSegmentedColormap('BLueRed', cdict)


class FlatElementError(Exception):
    pass

# def find_duplicate_linear_element(linear_elements1, linear_elements2):
#     duplicates = []
#     for linear_element in linear_elements1:
#         if linear_element in linear_elements2 and linear_element not in duplicates:
#             duplicates.append(linear_element)
#     return duplicates


class Node2D(vm.Point2D):
    def __hash__(self):
        return int(1e6 * (self.x + self.y))

    def __eq__(self, other_node: 'Node2D'):
        if other_node.__class__.__name__ not in ['Vector2D', 'Point2D',
                                                 'Node2D']:
            return False
        return math.isclose(self.x, other_node.x, abs_tol=1e-12) \
            and math.isclose(self.y, other_node.y, abs_tol=1e-12)

    @classmethod
    def from_point(cls, point2d):
        return cls(point2d.x, point2d.y)


class Node3D(vm.Point3D):
    def __hash__(self):
        return int(1e6 * (self.x + self.y + self.z))

    def __eq__(self, other_node: 'Node3D'):
        if other_node.__class__.__name__ not in ['Vector3D', 'Point3D',
                                                 'Node3D']:
            return False
        return math.isclose(self.x, other_node.x, abs_tol=1e-12) \
            and math.isclose(self.y, other_node.y, abs_tol=1e-12)

    @classmethod
    def from_point(cls, point3d):
        return cls(point3d.x, point3d.y, point3d.z)


class LinearElement(vme.LineSegment2D):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, start: vm.Point2D, end: vm.Point2D,
                 interior_normal: vm.Vector2D, name: str = ''):
        self.points = [start, end]
        self.interior_normal = interior_normal

        vme.LineSegment2D.__init__(self, start=start, end=end, name=name)

    # def __hash__(self):
    #     return self.start.__hash__() + self.end.__hash__()

    # def __eq__(self, other_linear_element):
    #     if self.__class__ != other_linear_element.__class__:
    #         return False
    #     return (self.start == other_linear_element.start and self.end == other_linear_element.end) \
    #         or (self.start== other_linear_element.end and self.end == other_linear_element.start)

    # def plot(self, ax=None, color='k', width=None, plot_points=False):
    #     if ax is None:
    #         fig, ax = plt.subplots()
    #         ax.set_aspect('equal')
    #     if width is None:
    #         width=1
    #     if plot_points:
    #         ax.plot([self.start.x, self.end.x], [self.start.y, self.end.y], color=color, marker='o', linewidth=width)
    #     else:
    #         ax.plot([self.start.x, self.end.x], [self.start.y, self.end.y], color=color, linewidth=width)
    #     return ax


class TriangularElement(vmw.Triangle):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points):
        self.points = points
        # self.linear_elements = self._to_linear_elements()
        # self.form_functions = self._form_functions()
        # self.line_segments = self.line_segments()
        self.center = (self.points[0] + self.points[1] + self.points[2]) / 3

        self._line_segments = None

        # self.area = self._area()

        # vmw.Triangle.__init__(self, points)

    def _to_linear_elements(self):
        vec1 = vm.Vector2D(self.points[1].x - self.points[0].x,
                           self.points[1].y - self.points[0].y)
        vec2 = vm.Vector2D(self.points[2].x - self.points[1].x,
                           self.points[2].y - self.points[1].y)
        vec3 = vm.Vector2D(self.points[0].x - self.points[2].x,
                           self.points[0].y - self.points[2].y)
        normal1 = vm.Vector2D(-vec1.y, vec1.x)
        normal2 = vm.Vector2D(-vec2.y, vec2.x)
        normal3 = vm.Vector2D(-vec3.y, vec3.x)
        normal1.normalize()
        normal2.normalize()
        normal3.normalize()
        if normal1.dot(vec2) < 0:
            normal1 = - normal1
        if normal2.dot(vec3) < 0:
            normal2 = - normal2
        if normal3.dot(vec1) < 0:
            normal3 = - normal3
        linear_element_1 = LinearElement(self.points[0], self.points[1],
                                         normal1)
        linear_element_2 = LinearElement(self.points[1], self.points[2],
                                         normal2)
        linear_element_3 = LinearElement(self.points[2], self.points[0],
                                         normal3)
        return [linear_element_1, linear_element_2, linear_element_3]

    def _form_functions(self):
        a_matrix = vm.Matrix33(1, self.points[0].x, self.points[0].y,
                               1, self.points[1].x, self.points[1].y,
                               1, self.points[2].x, self.points[2].y)
        try:
            inv_a = a_matrix.inverse()
        except ValueError as expt:
            # self.plot()
            print('buggy element area', self._area())
            raise FlatElementError('form function bug') from expt
        x_1 = inv_a.vector_multiplication(vm.X3D)
        x_2 = inv_a.vector_multiplication(vm.Y3D)
        x_3 = inv_a.vector_multiplication(vm.Z3D)
        return x_1, x_2, x_3

    # def _quadratic_form_functions(self):
    #     a = [[1, self.points[0][0], self.points[0][1],self.points[0][0]**2,
    #           self.points[0][0]*self.points[0][1],self.points[0][1]**2],
    #           [1, self.points[1][0], self.points[1][1],self.points[1][0]**2,
    #            self.points[1][0]*self.points[1][1],self.points[1][1]**2],
    #           [1, self.points[2][0], self.points[2][1],self.points[2][0]**2,
    #            self.points[2][0]*self.points[2][1],self.points[2][1]**2],
    #           [1, self.points[3][0], self.points[3][1],self.points[3][0]**2,
    #            self.points[3][0]*self.points[3][1],self.points[3][1]**2],
    #           [1, self.points[4][0], self.points[4][1],self.points[4][0]**2,
    #            self.points[4][0]*self.points[4][1],self.points[4][1]**2],
    #           [1, self.points[5][0], self.points[5][1],self.points[5][0]**2,
    #            self.points[5][0]*self.points[5][1],self.points[5][1]**2]]

    #     try :
    #         inv_a = a.inverse()
    #     except ValueError:
    #         self.plot()
    #         print(self._area())
    #         raise FlatElementError('form function bug')
    #     x1 = inv_a.dot([1,0,0,0,0,0])
    #     x2 = inv_a.dot([1,0,0,0,0,0])
    #     x3 = inv_a.dot([1,0,0,0,0,0])
    #     x4=inv_a.dot([1,0,0,0,0,0])

    #     return x1, x2, x3

    def _area(self):
        u_vect = self.points[1] - self.points[0]
        v_vect = self.points[2] - self.points[0]
        return abs(u_vect.cross(v_vect)) / 2

    # def point_belongs(self, point):
    #     polygon = volmdlr.wires.ClosedPolygon2D(self.points)
    #     point_belongs = polygon.point_belongs(point)
    #     return point_belongs

    # def rotation(self, center, angle, copy=True):
    #     if copy:
    #         return TriangularElement([pt.rotation(center, angle, copy=True)
    #                                   for pt in self.points])
    #     else:
    #         for pt in self.points:
    #             pt.Rotation(center, angle, copy=False)

    # def translation(self, offset, copy=True):
    #     if copy:
    #         return TriangularElement([pt.translation(offset, copy=True)
    #                                   for pt in self.points])
    #     else:
    #         for pt in self.points:
    #             pt.translation(offset, copy=False)

    def axial_symmetry(self, line):
        new_points = []
        for point in self.points:
            new_points.append(point.axial_symmetry(line))
        return self.__class__(new_points)

    # def axial_symmetry(self, line, copy=True):
    #     p1, p2 = line.points
    #     symmetric_points = []
    #     for point in self.points:
    #         u = p2 - p1
    #         t = (point-p1).dot(u) / u.norm()**2
    #         projection = p1 + t * u
    #         symmetric_point = volmdlr.Point2D(*(2 * projection - point))
    #         symmetric_points.append(symmetric_point)
    #     if copy:
    #         return TriangularElement(symmetric_points)
    #     else:
    #         for i, point in enumerate(self.points):
    #             point = symmetric_points[i]

    # def triangle_to_polygon(self):
    #     points = self.points
    #     return volmdlr.wires.ClosedPolygon2D(points)


class TriangularElement2D(TriangularElement, vmw.ClosedPolygon2D):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points, name: str = ''):
        self.points = points
        self.name = name
        self.linear_elements = self._to_linear_elements()
        self.form_functions = self._form_functions()
        # self._line_segments = None

        # self.line_segments = self.get_line_segments()
        self.center = (self.points[0] + self.points[1] + self.points[2]) / 3

        self.area = self._area()
        # vmw.Triangle.__init__(self, points)

    def _to_linear_elements(self):
        vec1 = vm.Vector2D(self.points[1].x - self.points[0].x,
                           self.points[1].y - self.points[0].y)
        vec2 = vm.Vector2D(self.points[2].x - self.points[1].x,
                           self.points[2].y - self.points[1].y)
        vec3 = vm.Vector2D(self.points[0].x - self.points[2].x,
                           self.points[0].y - self.points[2].y)
        normal1 = vm.Vector2D(-vec1.y, vec1.x)
        normal2 = vm.Vector2D(-vec2.y, vec2.x)
        normal3 = vm.Vector2D(-vec3.y, vec3.x)
        normal1.normalize()
        normal2.normalize()
        normal3.normalize()
        if normal1.dot(vec2) < 0:
            normal1 = - normal1
        if normal2.dot(vec3) < 0:
            normal2 = - normal2
        if normal3.dot(vec1) < 0:
            normal3 = - normal3
        linear_element_1 = LinearElement(self.points[0], self.points[1],
                                         normal1)
        linear_element_2 = LinearElement(self.points[1], self.points[2],
                                         normal2)
        linear_element_3 = LinearElement(self.points[2], self.points[0],
                                         normal3)
        return [linear_element_1, linear_element_2, linear_element_3]

    def _form_functions(self):
        a_matrix = vm.Matrix33(1, self.points[0].x, self.points[0].y,
                               1, self.points[1].x, self.points[1].y,
                               1, self.points[2].x, self.points[2].y)
        try:
            inv_a = a_matrix.inverse()
        except ValueError as expt:
            self.plot()
            print('buggy element area', self.area)
            raise FlatElementError('form function bug') from expt
        x_1 = inv_a.vector_multiplication(vm.X3D)
        x_2 = inv_a.vector_multiplication(vm.Y3D)
        x_3 = inv_a.vector_multiplication(vm.Z3D)
        return x_1, x_2, x_3

    # def _quadratic_form_functions(self):
    #     a = [[1, self.points[0][0], self.points[0][1],self.points[0][0]**2,
    #           self.points[0][0]*self.points[0][1],self.points[0][1]**2],
    #           [1, self.points[1][0], self.points[1][1],self.points[1][0]**2,
    #            self.points[1][0]*self.points[1][1],self.points[1][1]**2],
    #           [1, self.points[2][0], self.points[2][1],self.points[2][0]**2,
    #            self.points[2][0]*self.points[2][1],self.points[2][1]**2],
    #           [1, self.points[3][0], self.points[3][1],self.points[3][0]**2,
    #            self.points[3][0]*self.points[3][1],self.points[3][1]**2],
    #           [1, self.points[4][0], self.points[4][1],self.points[4][0]**2,
    #            self.points[4][0]*self.points[4][1],self.points[4][1]**2],
    #           [1, self.points[5][0], self.points[5][1],self.points[5][0]**2,
    #            self.points[5][0]*self.points[5][1],self.points[5][1]**2]]

    #     try :
    #         inv_a = a.inverse()
    #     except ValueError:
    #         self.plot()
    #         print(self._area())
    #         raise FlatElementError('form function bug')
    #     x1 = inv_a.dot([1,0,0,0,0,0])
    #     x2 = inv_a.dot([1,0,0,0,0,0])
    #     x3 = inv_a.dot([1,0,0,0,0,0])
    #     x4=inv_a.dot([1,0,0,0,0,0])

    #     return x1, x2, x3

    def _area(self):
        u_vect = self.points[1] - self.points[0]
        v_vect = self.points[2] - self.points[0]
        return abs(u_vect.cross(v_vect)) / 2

    # def point_belongs(self, point):
    #     polygon = volmdlr.wires.ClosedPolygon2D(self.points)
    #     point_belongs = polygon.point_belongs(point)
    #     return point_belongs

    # def rotation(self, center, angle, copy=True):
    #     if copy:
    #         return TriangularElement([pt.rotation(center, angle, copy=True)
    #                                   for pt in self.points])
    #     else:
    #         for pt in self.points:
    #             pt.Rotation(center, angle, copy=False)

    # def translation(self, offset, copy=True):
    #     if copy:
    #         return TriangularElement([pt.translation(offset, copy=True)
    #                                   for pt in self.points])
    #     else:
    #         for pt in self.points:
    #             pt.translation(offset, copy=False)

    def axial_symmetry(self, line):
        new_points = []
        for point in self.points:
            new_points.append(point.axial_symmetry(line))
        return self.__class__(new_points)

    # def plot(self, ax=None, color='k', width=None,
    #           plot_points=False, fill=False):
    #     if ax is None:
    #         fig, ax = plt.subplots()
    #         ax.set_aspect('equal')

    #     if fill:
    #         x = [p[0] for p in self.points]
    #         y = [p[1] for p in self.points]
    #         plt.fill(x, y, facecolor=color, edgecolor="k")
    #         return ax

    #     for p1, p2 in zip(self.points, self.points[1:]+[self.points[0]]):
    #         if width is None:
    #             width = 1
    #         if plot_points:
    #             ax.plot([p1.x, p2.x], [p1.y, p2.y], color=color,
    #                     marker='o', linewidth=width)
    #         else:
    #             ax.plot([p1.x, p2.x], [p1.y, p2.y], color=color,
    #                     linewidth=width)
    #     return ax

    # def triangle_to_polygon(self):
    #     points = self.points
    #     return volmdlr.wires.ClosedPolygon2D(points)

    def plot(self, ax=None, color='k', width=None,
             plot_points=False, fill=False):
        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect('equal')

        if fill:
            x = [p[0] for p in self.points]
            y = [p[1] for p in self.points]
            plt.fill(x, y, facecolor=color, edgecolor="k")
            return ax

        for p1, p2 in zip(self.points, self.points[1:] + [self.points[0]]):
            if width is None:
                width = 1
            if plot_points:
                ax.plot([p1.x, p2.x], [p1.y, p2.y], color=color,
                        marker='o', linewidth=width)
            else:
                ax.plot([p1.x, p2.x], [p1.y, p2.y], color=color,
                        linewidth=width)
        return ax


class QuadrilateralElement2D(vmw.ClosedPolygon2D):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points):
        self.points = points
        # self.linear_elements = self._to_linear_elements()
        # self.form_functions = self._form_functions()
        # self.line_segments = self.line_segments
        self.center = self.center_of_mass()

        self.area = self.area()
        self._line_segments = None
        vmw.ClosedPolygon2D.__init__(self, points)


class TriangularElement3D(TriangularElement, vmw.ClosedPolygon3D):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points):
        self.points = points
        # self.linear_elements = self._to_linear_elements()
        # self.form_functions = self._form_functions()
        # self.line_segments = self.line_segments
        self.center = (self.points[0] + self.points[1] + self.points[2]) / 3

        # self.area = self._area()
        self._line_segments = None
        TriangularElement.__init__(self, points)

    # def _to_linear_elements(self):
    #     vec1 = vm.Vector2D(self.points[1].x - self.points[0].x,
    #                        self.points[1].y - self.points[0].y)
    #     vec2 = vm.Vector2D(self.points[2].x - self.points[1].x,
    #                        self.points[2].y - self.points[1].y)
    #     vec3 = vm.Vector2D(self.points[0].x - self.points[2].x,
    #                        self.points[0].y - self.points[2].y)
    #     normal1 = vm.Vector2D(-vec1.y, vec1.x)
    #     normal2 = vm.Vector2D(-vec2.y, vec2.x)
    #     normal3 = vm.Vector2D(-vec3.y, vec3.x)
    #     normal1.normalize()
    #     normal2.normalize()
    #     normal3.normalize()
    #     if normal1.dot(vec2) < 0:
    #         normal1 = - normal1
    #     if normal2.dot(vec3) < 0:
    #         normal2 = - normal2
    #     if normal3.dot(vec1) < 0:
    #         normal3 = - normal3
    #     linear_element_1 = LinearElement(self.points[0], self.points[1],
    #                                      normal1)
    #     linear_element_2 = LinearElement(self.points[1], self.points[2],
    #                                      normal2)
    #     linear_element_3 = LinearElement(self.points[2], self.points[0],
    #                                      normal3)
    #     return [linear_element_1, linear_element_2, linear_element_3]

    # def _form_functions(self):
    #     a = vm.Matrix33(1, self.points[0].x, self.points[0].y,
    #                     1, self.points[1].x, self.points[1].y,
    #                     1, self.points[2].x, self.points[2].y)
    #     try:
    #         inv_a = a.inverse()
    #     except ValueError as expt:
    #         self.plot()
    #         print('buggy element area', self.area)
    #         raise FlatElementError('form function bug') from expt
    #     x1 = inv_a.vector_multiplication(vm.X3D)
    #     x2 = inv_a.vector_multiplication(vm.Y3D)
    #     x3 = inv_a.vector_multiplication(vm.Z3D)
    #     return x1, x2, x3

    # def _quadratic_form_functions(self):
    #     a = [[1, self.points[0][0], self.points[0][1],self.points[0][0]**2,
    #           self.points[0][0]*self.points[0][1],self.points[0][1]**2],
    #           [1, self.points[1][0], self.points[1][1],self.points[1][0]**2,
    #            self.points[1][0]*self.points[1][1],self.points[1][1]**2],
    #           [1, self.points[2][0], self.points[2][1],self.points[2][0]**2,
    #            self.points[2][0]*self.points[2][1],self.points[2][1]**2],
    #           [1, self.points[3][0], self.points[3][1],self.points[3][0]**2,
    #            self.points[3][0]*self.points[3][1],self.points[3][1]**2],
    #           [1, self.points[4][0], self.points[4][1],self.points[4][0]**2,
    #            self.points[4][0]*self.points[4][1],self.points[4][1]**2],
    #           [1, self.points[5][0], self.points[5][1],self.points[5][0]**2,
    #            self.points[5][0]*self.points[5][1],self.points[5][1]**2]]

    #     try :
    #         inv_a = a.inverse()
    #     except ValueError:
    #         self.plot()
    #         print(self._area())
    #         raise FlatElementError('form function bug')
    #     x1 = inv_a.dot([1,0,0,0,0,0])
    #     x2 = inv_a.dot([1,0,0,0,0,0])
    #     x3 = inv_a.dot([1,0,0,0,0,0])
    #     x4=inv_a.dot([1,0,0,0,0,0])

    #     return x1, x2, x3

    # def _area(self):
    #     u = self.points[1] - self.points[0]
    #     v = self.points[2] - self.points[0]
    #     return abs(u.cross(v)) / 2

    # def point_belongs(self, point):
    #     polygon = volmdlr.wires.ClosedPolygon2D(self.points)
    #     point_belongs = polygon.point_belongs(point)
    #     return point_belongs

    # def rotation(self, center, angle, copy=True):
    #     if copy:
    #         return TriangularElement([pt.rotation(center, angle, copy=True)
    #                                   for pt in self.points])
    #     else:
    #         for pt in self.points:
    #             pt.Rotation(center, angle, copy=False)

    # def translation(self, offset, copy=True):
    #     if copy:
    #         return TriangularElement([pt.translation(offset, copy=True)
    #                                   for pt in self.points])
    #     else:
    #         for pt in self.points:
    #             pt.translation(offset, copy=False)

    def axial_symmetry(self, line):
        new_points = []
        for point in self.points:
            new_points.append(point.axial_symmetry(line))
        return self.__class__(new_points)

    # def plot(self, ax=None, color='k', width=None,
    #           plot_points=False, fill=False):
    #     if ax is None:
    #         fig, ax = plt.subplots()
    #         ax.set_aspect('equal')

    #     if fill:
    #         x = [p[0] for p in self.points]
    #         y = [p[1] for p in self.points]
    #         plt.fill(x, y, facecolor=color, edgecolor="k")
    #         return ax

    #     for p1, p2 in zip(self.points, self.points[1:]+[self.points[0]]):
    #         if width is None:
    #             width = 1
    #         if plot_points:
    #             ax.plot([p1.x, p2.x], [p1.y, p2.y], color=color,
    #                     marker='o', linewidth=width)
    #         else:
    #             ax.plot([p1.x, p2.x], [p1.y, p2.y], color=color,
    #                     linewidth=width)
    #     return ax

    # def triangle_to_polygon(self):
    #     points = self.points
    #     return volmdlr.wires.ClosedPolygon2D(points)


class TetrahedralElement(TriangularElement, vmw.ClosedPolygon3D):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points, name: str = ''):
        self.points = points
        self.name = name
        # self.linear_elements = self._to_linear_elements()
        self.form_functions = self._form_functions()
        # self.line_segments = self._line_segments()
        self.center = (self.points[0] + self.points[1] + self.points[2] + self.points[3]) / 4
        self.triangular_elements = self._triangular_elements()

        self.volume = self._volume()
        DessiaObject.__init__(self, name=name)

    def _triangular_elements(self):

        indices_combinations = list(combinations(list(range(len(self.points))), r=3))
        # indices_combinations = [x for x in combinations(list(range(len(self.points))), r=3)]
        triangular_elements = []

        for indices in indices_combinations:
            triangular_elements.append(TriangularElement3D([self.points[indices[0]],
                                                            self.points[indices[1]],
                                                            self.points[indices[2]]]))

        return triangular_elements

    def plot(self, ax=None, color='k'):
        if ax is None:
            ax = plt.figure().add_subplot(projection='3d')
        for point in self.points:
            point.plot(ax=ax)
        for triangle in self.triangular_elements:
            triangle.plot(ax=ax)
        return ax

    def _volume(self):

        data = []
        for i in range(3):
            data.extend([*self.points[i + 1] - self.points[0]])

        return abs(1 / 6 * (npy.linalg.det(npy.array(data).reshape(3, 3))))

    def _form_functions(self):
        # coeff = [1, -1, 1, 1]
        # alpha = []
        # for i in range(4):
        #     data = []
        #     for c in range(4):
        #         if c != i:
        #             data.extend([self.points[c].x, self.points[c].y, self.points[c].z])
        #     alpha.append(coeff[i] * (npy.linalg.det(npy.array(data).reshape(3,3))))
        # gamma = []
        # for i in range(4):
        #     data = []
        #     for c in range(4):
        #         if c != i:
        #             data.extend([1, self.points[c].x, self.points[c].z])
        #     gamma.append(coeff[i] * (npy.linalg.det(npy.array(data).reshape(3,3))))

        # coeff = [-1, 1, -1, 1]
        # betha = []
        # for i in range(4):
        #     data = []
        #     for c in range(4):
        #         if c != i:
        #             data.extend([1, self.points[c].y, self.points[c].z])
        #     betha.append(coeff[i] * (npy.linalg.det(npy.array(data).reshape(3,3))))
        # delta = []
        # for i in range(4):
        #     data = []
        #     for c in range(4):
        #         if c != i:
        #             data.extend([1, self.points[c].x, self.points[c].y])
        #     delta.append(coeff[i] * (npy.linalg.det(npy.array(data).reshape(3,3))))

        coeff = [1, -1, 1, -1]
        form_funct = []
        for i in range(4):
            data_alpha, data_gamma, data_betha, data_delta = [], [], [], []
            for c_coef in range(4):
                if c_coef != i:
                    data_alpha.extend([self.points[c_coef].x, self.points[c_coef].y, self.points[c_coef].z])
                    data_gamma.extend([1, self.points[c_coef].x, self.points[c_coef].z])
                    data_betha.extend([1, self.points[c_coef].y, self.points[c_coef].z])
                    data_delta.extend([1, self.points[c_coef].x, self.points[c_coef].y])

            form_funct.append([(coeff[i] * (npy.linalg.det(npy.array(data_alpha).reshape(3, 3)))),
                               ((-1) * coeff[i] * (npy.linalg.det(npy.array(data_betha).reshape(3, 3)))),
                               (coeff[i] * (npy.linalg.det(npy.array(data_gamma).reshape(3, 3)))),
                               ((-1) * coeff[i] * (npy.linalg.det(npy.array(data_delta).reshape(3, 3))))])

        return form_funct[0], form_funct[1], form_funct[2], form_funct[3]


class ElementsGroup(DessiaObject):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, elements, name: str):
        self.elements = elements
        self.nodes = self._nodes()
        self.name = name

        self._elements_per_node = None

        DessiaObject.__init__(self, name=name)

    def _nodes(self):
        nodes = set()
        for element in self.elements:
            for point in element.points:
                nodes.add(point)
        return nodes

    def point_to_element(self, point):
        for element in self.elements:
            if element.point_belongs(point):
                return element
        return None

    @property
    def elements_per_node(self):
        if self._elements_per_node is not None:
            return self._elements_per_node

        dict_node_element = {}
        for element in self.elements:
            for point in element.points:
                try:
                    dict_node_element[point].append(element)
                except KeyError:
                    dict_node_element[point] = [element]
        self._elements_per_node = dict_node_element
        return dict_node_element

    # def rotation(self, center, angle, copy=True):
    #     if copy:
    #         return Mesh([elem.rotation(center, angle, copy=True)
    #                       for elem in self.elements])
    #     else:
    #         for elem in self.elements:
    #             elem.rotation(center, angle, copy=False)

    # def translation(self, offset, copy=True):
    #     if copy:
    #         return Mesh([elem.translation(offset, copy=True)
    #                       for elem in self.elements])
    #     else:
    #         for elem in self.elements:
    #             elem.translation(offset, copy=False)

    def plot(self, ax=None, color='k'):  # , fill=False):
        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect('equal')
        for element in self.elements:
            element.plot(ax=ax, color=color)  # fill=fill
        return ax


class Mesh(DessiaObject):
    _standalone_in_db = True
    _non_serializable_attributes = ['node_to_index']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, elements_groups: List[ElementsGroup]):
        self.elements_groups = elements_groups
        self.nodes = self._set_nodes_number()
        self.node_to_index = {self.nodes[i]: i for i in range(len(self.nodes))}

        DessiaObject.__init__(self, name='')

    # def __add__(self, other_mesh):
    #     new_nodes = self.nodes[:]
    #     new_nodes_index = {p: i for i, p in enumerate(self.points)}
    #     ip = len(new_nodes)
    #     for point in other_mesh.nodes:
    #         if point not in new_nodes_index:
    #             new_nodes_index[point] = ip
    #             ip += 1
    #             new_nodes.append(point)

    #     new_elements_groups = self.elements_groups[:]
    #     for i1, i2, i3 in other_mesh.elements_groups:
    #         p1 = other_mesh.nodes[i1]
    #         p2 = other_mesh.nodes[i2]
    #         p3 = other_mesh.nodes[i3]
    #         new_elements_groups.append((new_nodes_index[p1],
    #                                     new_nodes_index[p2],
    #                                     new_nodes_index[p3]))
    #     return self.__class__(new_elements_groups)

    def _set_nodes_number(self):
        nodes = set()
        for elements_group in self.elements_groups:
            for element in elements_group.elements:
                for point in element.points:
                    nodes.add(point)
                # nodes.add(element.points[0])
                # nodes.add(element.points[1])
                # nodes.add(element.points[2])
        return tuple(nodes)

    def point_to_element(self, point):
        for element_group in self.elements_groups:
            element = element_group.point_to_element(point)
            if element is not None:
                return element
        return None

    # def set_node_displacement_index(self):
    #     indexes = {}
    #     for node in self.nodes:
    #         indexes[node] = [2*self.node_to_index[node],
    #                           2*self.node_to_index[node]+1]
    #     return indexes

    # def boundary_dict(self):
    #     boundary_dict = {}
    #     for elements_group1, elements_group2 in combinations(
    #             self.elements_groups, 2):
    #         linear_elements1 = []
    #         linear_elements2 = []
    #         for element in elements_group1.elements:
    #             linear_elements1.extend(element.linear_elements)
    #         for element in elements_group2.elements:
    #             linear_elements2.extend(element.linear_elements)
    #         duplicate_linear_elements = find_duplicate_linear_element(
    #             linear_elements1, linear_elements2)
    #         if duplicate_linear_elements:
    #             boundary_dict[(elements_group1,
    #                             elements_group2)] = duplicate_linear_elements
    #     return boundary_dict

    def plot(self, ax=None):
        if ax is None:
            if self.elements_groups[0].elements[0].__class__.__name__[-2::] == '2D':
                _, ax = plt.subplots()
                ax.set_aspect('equal')
            else:
                ax = plt.figure().add_subplot(projection='3d')
        for elements_group in self.elements_groups:
            elements_group.plot(ax=ax)
        return ax

    # def plot_data(self, pos=0, quote=True, constructor=True, direction=1):
    #     plot_datas = []
    #     for element_group in self.elements_groups:
    #         for element in element_group.elements:
    #             c1 = volmdlr.wires.Contour2D([volmdlr.edges.LineSegment2D(
    #                 element.points[0], element.points[1])])
    #             c2 = volmdlr.wires.Contour2D([volmdlr.edges.LineSegment2D(
    #                 element.points[1], element.points[2])])
    #             c3 = volmdlr.wires.Contour2D([volmdlr.edges.LineSegment2D(
    #                 element.points[2], element.points[0])])
    #             plot_datas.append(c1.plot_data())
    #             plot_datas.append(c2.plot_data())
    #             plot_datas.append(c3.plot_data())
    #             # plot_datas.extend([c1, c2, c3])
    #     return plot_datas

    # def plot_displaced_mesh(self,
    #                         node_displacement: Dict[volmdlr.Point2D,
    #                                                 List[float]],
    #                         ax=None, amplification=0.5):

    #     deformed_mesh = self.copy()
    #     nodes = deformed_mesh.nodes

    #     for node in nodes:
    #         for displaced_node in node_displacement:
    #             if node == displaced_node:
    #                 node.x += amplification*node_displacement[
    #                     displaced_node][0]
    #                 node.y += amplification*node_displacement[
    #                     displaced_node][1]

    #     ax = deformed_mesh.plot(ax=ax)
    #     ax.set_aspect('equal')

    #     return ax

    def bounding_rectangle(self):
        nodes = self.nodes
        x, y = [], []
        for n in nodes:
            x.append(n.x)
            y.append(n.y)
        if len([*nodes[0]]) == 2:
            return min(x), max(x), min(y), max(y)

        z = [n.z for n in nodes]
        return min(x), max(x), min(y), max(y), min(z), max(z)

    def delete_duplicated_nodes(self, tol=1e-4):
        mesh = self.__class__(self.elements_groups[:])
        nodes_list = list(mesh.nodes[:])
        nodes_index = []

        for i, node in enumerate(nodes_list):
            for j in range(i + 1, len(nodes_list)):
                dist = node.point_distance(nodes_list[j])
                if dist < tol:
                    nodes_index.append((j, i))

        if nodes_index:
            nodes_index = sorted(nodes_index, key=lambda item: item[0], reverse=True)
            for _, index in enumerate(nodes_index):
                nodes_list.pop(index[0])
                for group in mesh.elements_groups:
                    if mesh.nodes[index[0]] in group.nodes:
                        dict_node_element = group.elements_per_node
                        for element in dict_node_element[mesh.nodes[index[0]]]:
                            element.points[element.points.index(mesh.nodes[index[0]])] = mesh.nodes[index[1]]

            mesh.nodes = nodes_list
            mesh.node_to_index = {mesh.nodes[i]: i for i in range(len(mesh.nodes))}

        return mesh
