#!/usr/bin/env python3
# -*- coding: utf-8 -*-
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


class Node2D(vm.Point2D):
    """
    A node is a Point2D with some hash capabilities for performance used for Mesh.

    """

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
        """
        Defines a node2d from a point2d.

        :param point2d: A point2d
        :type point2d: vm.Point2D
        :return: A node2d
        :rtype: Node2D
        """

        return cls(point2d.x, point2d.y)


class Node3D(vm.Point3D):
    """
    A node is a Point3D with some hash capabilities for performance used for Mesh.
    """

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
        """
        Defines a node3d from a point3d.

        :param point3d: A point3d
        :type point3d: vm.Point3D
        :return: A node3d
        :rtype: Node3D
        """

        return cls(point3d.x, point3d.y, point3d.z)


class LinearElement(vme.LineSegment2D):
    """ A class that defines a linear element. """
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


class TriangularElement(vmw.Triangle):
    """
    A mesh element defined with 3 nodes.
    """
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points: List[volmdlr.Point2D]):
        super().__init__(*points)
        self.points = points
        self.center = (self.points[0] + self.points[1] + self.points[2]) / 3

        self._line_segments = None

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

    def _area(self):
        u_vect = self.points[1] - self.points[0]
        v_vect = self.points[2] - self.points[0]
        return abs(u_vect.cross(v_vect)) / 2

    def axial_symmetry(self, line):
        new_points = []
        for point in self.points:
            new_points.append(point.axial_symmetry(line))
        return self.__class__(new_points)


class TriangularElement2D(TriangularElement, vmw.ClosedPolygon2D):
    """ Class to define a 2D triangular element. """
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points, name: str = ''):
        super().__init__(points)
        # self.points = points
        self.name = name
        self.linear_elements = self._to_linear_elements()
        self.form_functions = self._form_functions()
        self.center = (self.points[0] + self.points[1] + self.points[2]) / 3

        self.area = self._area()


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
            for ip, point in enumerate(self.points):
                ax.text(*point, f'point {ip + 1}', ha='center', va='top')
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
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points: List[volmdlr.Point2D]):
        self.points = points
        self.center = self.center_of_mass()

        self.area = self.area()
        self._line_segments = None
        vmw.ClosedPolygon2D.__init__(self, points)


class TriangularElement3D(TriangularElement, vmw.ClosedPolygon3D):
    """ Class to define a 3D triangular element. """

    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points):
        self.points = points
        self.center = (self.points[0] + self.points[1] + self.points[2]) / 3

        self._line_segments = None
        TriangularElement.__init__(self, points)

    def axial_symmetry(self, line):
        """ Returns a symmetric new element with respect to the given line. """
        new_points = []
        for point in self.points:
            new_points.append(point.axial_symmetry(line))
        return self.__class__(new_points)


class TetrahedralElement(TriangularElement, vmw.ClosedPolygon3D):
    """ Class to define a 3D tetrahedral element. """

    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self, points, name: str = ''):
        self.points = points
        self.name = name
        self.form_functions = self._form_functions()
        self.center = (self.points[0] + self.points[1] + self.points[2] + self.points[3]) / 4
        self.triangular_elements = self._triangular_elements()

        self.volume = self._volume()
        DessiaObject.__init__(self, name=name)

    def _triangular_elements(self):

        indices_combinations = list(combinations(list(range(len(self.points))), r=3))
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


class ElementsGroup(DessiaObject):
    """Defines a group of elements."""

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

    def plot(self, ax=None, color='k'):
        """
        Plot an ElementGroup.
        """

        if ax is None:
            _, ax = plt.subplots()
            ax.set_aspect('equal')
        for element in self.elements:
            element.plot(ax=ax, edge_style=EdgeStyle(color=color))
        return ax


class Mesh(DessiaObject):
    """Defines a mesh."""

    _standalone_in_db = True
    _non_serializable_attributes = ['node_to_index']
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
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

        groups = self.elements_groups[:]
        groups.pop(reference_index)

        nodes_correction = self.nodes_correction(reference_index, tol)

        count = 0
        old_elements, new_elements = set(), set()
        for g, group in enumerate(groups):
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

            groups[g] = group.__class__(elements, name='')

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
        m = self.__class__(elements_groups=self.elements_groups[:])
        m.set_nodes_correction(self.get_nodes_correction())
        m.gmsh = self.gmsh

        return m
