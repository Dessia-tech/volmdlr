#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=invalid-unary-operand-type
"""
Module containing mesh and relative objects.
"""

import math
from itertools import combinations
from typing import List

import matplotlib.pyplot as plt
import numpy as npy

from dessia_common.core import DessiaObject  # isort: skip

import volmdlr as vm
import volmdlr.edges as vme
import volmdlr.gmsh_vm
import volmdlr.wires as vmw
from volmdlr.core import EdgeStyle


class FlatElementError(Exception):
    """An error in case an element is flat."""


# def find_duplicate_linear_element(linear_elements1, linear_elements2):
#     duplicates = []
#     for linear_element in linear_elements1:
#         if linear_element in linear_elements2 and linear_element not in duplicates:
#             duplicates.append(linear_element)
#     return duplicates


class Node2D(vm.Point2D):
    """
    A node is a Point2D with some hash capabilities for performance used for Mesh.

    """

    def __init__(self, x: float, y: float, name: str = ""):
        self.x = x
        self.y = y
        vm.Point2D.__init__(self, x, y, name)

    def __hash__(self):
        return int(1e6 * (self.x + self.y))

    def __eq__(self, other_node: 'Node2D'):
        if other_node.__class__.__name__ not in ['Vector2D', 'Point2D',
                                                 'Node2D']:
            return False
        return math.isclose(self.x, other_node.x, abs_tol=1e-12) \
            and math.isclose(self.y, other_node.y, abs_tol=1e-12)

    @classmethod
    def from_point(cls, point2d, name: str = ''):
        """
        Defines a node2d from a point2d.

        :param point2d: A point2d
        :type point2d: vm.Point2D.
        :param name: object's name.
        :return: A node2d
        :rtype: Node2D
        """

        return cls(point2d.x, point2d.y, name=name)


class Node3D(vm.Point3D):
    """
    A node is a Point3D with some hash capabilities for performance used for Mesh.
    """

    def __init__(self, x: float, y: float, z: float, name: str = ""):
        self.x = x
        self.y = y
        self.z = z
        vm.Point3D.__init__(self, x, y, z, name)

    def __hash__(self):
        return int(1e6 * (self.x + self.y + self.z))

    def __eq__(self, other_node: 'Node3D'):
        if other_node.__class__.__name__ not in ['Vector3D', 'Point3D',
                                                 'Node3D']:
            return False
        return math.isclose(self.x, other_node.x, abs_tol=1e-12) \
            and math.isclose(self.y, other_node.y, abs_tol=1e-12)

    @classmethod
    def from_point(cls, point3d, name: str = ''):
        """
        Defines a node3d from a point3d.

        :param point3d: A point3d
        :type point3d: vm.Point3D.
        :param name: object's name.
        :return: A node3d
        :rtype: Node3D
        """

        return cls(point3d.x, point3d.y, point3d.z, name=name)


class LinearElement(vme.LineSegment2D):
    """ A class that defines a linear element. """
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
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
    """
    A mesh element defined with 3 nodes.
    """
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points: List[volmdlr.Point2D]):
        # super().__init__(*points)
        self.points = points
        # self.linear_elements = self._to_linear_elements()
        # self.form_functions = self._form_functions()
        # self.line_segments = self.line_segments()
        self.center = (self.points[0] + self.points[1] + self.points[2]) / 3

        self._line_segments = None

        # self.area = self._area()

        vmw.Triangle.__init__(self, *points)

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
        normal1 = normal1.unit_vector()
        normal2 = normal2.unit_vector()
        normal3 = normal3.unit_vector()
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
    """ Class to define a 2D triangular element. """
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points, name: str = ''):
        super().__init__(points)
        # self.points = points
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
        normal1 = normal1.unit_vector()
        normal2 = normal2.unit_vector()
        normal3 = normal3.unit_vector()
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

    def _area(self):
        u_vect = self.points[1] - self.points[0]
        v_vect = self.points[2] - self.points[0]
        return abs(u_vect.cross(v_vect)) / 2

    def axial_symmetry(self, line):
        new_points = []
        for point in self.points:
            new_points.append(point.axial_symmetry(line))
        return self.__class__(new_points)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle(), point_numbering=False,
             fill=False, fill_color='w'):
        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect('equal')

        if fill:
            x = [p[0] for p in self.points]
            y = [p[1] for p in self.points]
            plt.fill(x, y, facecolor=fill_color, edgecolor="k")
            return ax
        if point_numbering:
            for index_point, point in enumerate(self.points):
                ax.text(*point, f'point {index_point + 1}', ha='center', va='top')
        for p1, p2 in zip(self.points, self.points[1:] + [self.points[0]]):
            if edge_style.width is None:
                edge_style.width = 1
            if edge_style.plot_points:
                ax.plot([p1.x, p2.x], [p1.y, p2.y], color=edge_style.color,
                        marker='o', linewidth=edge_style.linewidth)
            else:
                ax.plot([p1.x, p2.x], [p1.y, p2.y], color=edge_style.color,
                        linewidth=edge_style.linewidth)
        return ax


class QuadrilateralElement2D(vmw.ClosedPolygon2D):
    """ Class to define a 2D quadrilateral element. """

    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points: List[volmdlr.Point2D]):
        self.points = points
        # self.linear_elements = self._to_linear_elements()
        # self.form_functions = self._form_functions()
        # self.line_segments = self.line_segments
        self.center = self.center_of_mass()

        self.area = self.area()
        self._line_segments = None
        vmw.ClosedPolygon2D.__init__(self, points)


class TriangularElement3D(TriangularElement, vmw.ClosedPolygon3D):
    """ Class to define a 3D triangular element. """

    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points):
        self.points = points
        # self.linear_elements = self._to_linear_elements()
        # self.form_functions = self._form_functions()
        # self.line_segments = self.line_segments
        self.center = (self.points[0] + self.points[1] + self.points[2]) / 3

        # self.area = self._area()
        self._line_segments = None
        super().__init__(points)

    def axial_symmetry(self, line):
        """ Returns a symmetric new element with respect to the given line. """
        new_points = []
        for point in self.points:
            new_points.append(point.axial_symmetry(line))
        return self.__class__(new_points)


class TetrahedralElement(DessiaObject):
    """ Class to define a 3D tetrahedral element. """

    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
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

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            ax = plt.figure().add_subplot(projection='3d')
        for point in self.points:
            point.plot(ax=ax)
        for triangle in self.triangular_elements:
            triangle.plot(ax=ax, edge_style=edge_style)
        return ax

    def _volume(self):

        data = []
        for i in range(3):
            data.extend([*self.points[i + 1] - self.points[0]])

        return abs(1 / 6 * (npy.linalg.det(npy.array(data).reshape(3, 3))))

    def _form_functions(self):
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


class TetrahedralElementQuadratic(DessiaObject):
    """ Class to define a 3D quadratic tetrahedral element. """

    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points, name: str = ''):
        self.points = points
        self.name = name
        DessiaObject.__init__(self, name=name)


class ElementsGroup(DessiaObject):
    """Defines a group of elements."""

    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
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
        """
        Checks if a point belongs to an element.

        :param point: A point2d/3d
        :type point: Point2D/Point3D

        :return: The input Element OR None
        """

        for element in self.elements:
            if element.point_belongs(point):
                return element
        return None

    @property
    def elements_per_node(self):
        """Maps the elements of each node."""
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
        """
        Plot an ElementGroup.
        """

        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect('equal')
        for element in self.elements:
            element.plot(ax=ax, edge_style=EdgeStyle(color=color))  # fill=fill
        return ax


class Mesh(DessiaObject):
    """Defines a mesh."""

    _standalone_in_db = True
    _non_serializable_attributes = ['node_to_index']
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, elements_groups: List[ElementsGroup]):
        self.elements_groups = elements_groups
        self.nodes = self._set_nodes_number()
        self.node_to_index = {self.nodes[i]: i for i in range(len(self.nodes))}
        self._nodes_correction = {}
        self._gmsh = None
        DessiaObject.__init__(self, name='')

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

    # def delete_duplicated_nodes(self, tol=1e-4):
    #     mesh = self.__class__(self.elements_groups[:])
    #     nodes_list = list(mesh.nodes[:])
    #     nodes_index = []

    #     for i, node in enumerate(nodes_list):
    #         for j in range(i + 1, len(nodes_list)):
    #             dist = node.point_distance(nodes_list[j])
    #             if dist < tol:
    #                 nodes_index.append((j, i))

    #     if nodes_index:
    #         nodes_index = sorted(nodes_index, key=lambda item: item[0], reverse=True)
    #         for _, index in enumerate(nodes_index):
    #             nodes_list.pop(index[0])
    #             for group in mesh.elements_groups:
    #                 if mesh.nodes[index[0]] in group.nodes:
    #                     dict_node_element = group.elements_per_node
    #                     for element in dict_node_element[mesh.nodes[index[0]]]:
    #                         element.points[element.points.index(mesh.nodes[index[0]])] = mesh.nodes[index[1]]

    #         mesh.nodes = nodes_list
    #         mesh.node_to_index = {mesh.nodes[i]: i for i in range(len(mesh.nodes))}

    #     return mesh

    def nodes_correction(self, reference_index, tol=1e-4):
        if not self._nodes_correction:
            nodes_reference = self.elements_groups[reference_index].nodes
            groups = self.elements_groups[:]
            groups.pop(reference_index)
            nodes_correction = {}

            for group in groups:
                for node in group.nodes:
                    for node_ref in nodes_reference:
                        d = node.point_distance(node_ref)
                        if 1e-8 < d < 1e-4:
                            nodes_correction[node] = node_ref

            self._nodes_correction = nodes_correction

        return self._nodes_correction

    def delete_duplicated_nodes(self, reference_index, tol=1e-4):
        """Delete duplicated nodes."""
        groups = self.elements_groups[:]
        groups.pop(reference_index)

        nodes_correction = self.nodes_correction(reference_index, tol)

        count = 0
        old_elements, new_elements = set(), set()
        for index_g, group in enumerate(groups):
            elements = []
            for element in group.elements:
                new = False
                points = []

                for point in element.points:
                    correc_point = nodes_correction.get(point)
                    if correc_point is not None:
                        points.append(correc_point)
                        count += 1
                        old_elements.add(element)
                        new = True
                    else:
                        points.append(point)

                elements.append(element.__class__(points))
                if new:
                    new_elements.add(element.__class__(points))

            groups[index_g] = group.__class__(elements, name='')

        groups.insert(reference_index, self.elements_groups[reference_index])

        return self._helper_create_mesh(groups)

    def _helper_create_mesh(self, groups):
        mesh = self.__class__(groups)
        mesh.gmsh = self.gmsh
        mesh.set_nodes_correction(self.get_nodes_correction())
        return mesh

    def get_nodes_correction(self):
        """
        A getter method for nodes_correction private variable.

        :return: A dict of nodes_correction
        :rtype: dict
        """

        return self._nodes_correction

    def set_nodes_correction(self, nodes_correction):
        """
        A setter method for nodes_correction private variable.

        :param nodes_correction: A dict of nodes_correction
        :type nodes_correction: dict
        """

        if not isinstance(nodes_correction, dict):
            raise ValueError("It must be volmdlr.GmshParser class")
        self._nodes_correction = nodes_correction

    @property
    def gmsh(self):
        """
        A property to get gmsh (a private variable).

        :return: A gmsh_parser data
        :rtype: GmshParser
        """

        return self._gmsh

    @gmsh.setter
    def gmsh(self, gmsh_parser):
        """
        A setter for gmsh (a private variable).

        :param gmsh_parser: A gmsh_parser data
        :type gmsh_parser: GmshParser
        """

        if not isinstance(gmsh_parser, volmdlr.gmsh_vm.GmshParser):
            raise ValueError("It must be volmdlr.GmshParser class")
        self._gmsh = gmsh_parser

    def copy(self):
        mesh = self.__class__(elements_groups=self.elements_groups[:])
        mesh.set_nodes_correction(self.get_nodes_correction())
        mesh.gmsh = self.gmsh

        return mesh
