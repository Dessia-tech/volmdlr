#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
from typing import List, Tuple
import math
import matplotlib.pyplot as plt
import dessia_common as dc
import volmdlr.edges

class Node2D(volmdlr.Point2D):
    def __hash__(self):
        return int(1e6*(self.x+self.y))

    def __eq__(self, other_node:'Node2D'):
        if other_node.__class__.__name__ not in ['Vector2D', 'Point2D',
                                                 'Node2D']:
            return False
        return math.isclose(self.x, other_node.x, abs_tol=1e-06) \
        and math.isclose(self.y, other_node.y, abs_tol=1e-06)


class Node3D(volmdlr.Point3D):
    def __hash__(self):
        return int(1e6*(self.x+self.y+self.z))

    def __eq__(self, other_node:'Node3D'):
        if other_node.__class__.__name__ not in ['Vector3D', 'Point3D',
                                                 'Node3D']:
            return False
        return math.isclose(self.x, other_node.x, abs_tol=1e-06) \
        and math.isclose(self.y, other_node.y, abs_tol=1e-06) \
        and math.isclose(self.z, other_node.z, abs_tol=1e-06)

class DisplayMesh(dc.DessiaObject):
    def __init__(self, points, triangles, edges=None, name=''):

        self.points = points
        self.triangles = triangles
        if edges is None:
            edges = []
        self.edges = edges
        self.name = name
        
        
    def check(self):
        npoints = len(self.points)
        for triangle in self.triangles:
            if max(triangle) >= npoints:
                return False
        return True
            

    def __add__(self, other_mesh):
        new_points = self.points[:]
        # print(self.points[0].__class__.__name__)
        new_point_index = {p: i for i, p in enumerate(self.points)}
        ip = len(new_points)
        for point in other_mesh.points:
            if not point in new_point_index:
                new_point_index[point] = ip
                ip += 1
                new_points.append(point)

        new_triangles = self.triangles[:]
        for i1, i2, i3 in other_mesh.triangles:
            p1 = other_mesh.points[i1]
            p2 = other_mesh.points[i2]
            p3 = other_mesh.points[i3]
            new_triangles.append((new_point_index[p1],
                                  new_point_index[p2],
                                  new_point_index[p3]))

        return self.__class__(new_points, new_triangles)

    def plot(self, ax=None, numbering=False):
        for ip, point in enumerate(self.points):
            ax = point.plot(ax=ax)
            if numbering:
                ax.text(*point, 'node {}'.format(ip+1),
                        ha='center', va='center')

        for i1, i2, i3 in self.triangles:
            self._linesegment_class(self.points[i1], self.points[i2]).plot(
                ax=ax)
            self._linesegment_class(self.points[i2], self.points[i3]).plot(
                ax=ax)
            self._linesegment_class(self.points[i1], self.points[i3]).plot(
                ax=ax)

        # for i, (i1, i2) in enumerate(self.edges):
        #     self._linesegment_class(self.points[i1], self.points[i2]).plot(
        #         ax=ax)
        #     if numbering:
        #         ax.text(*0.5*(self.points[i1]+self.points[i2]), 'edge {}'.format(i+1),
        #                 ha='center', va='center')

        return ax

        

class DisplayMesh2D(DisplayMesh):
    _linesegment_class = volmdlr.edges.LineSegment2D
    _point_class = volmdlr.Point2D

    def __init__(self, points: List[volmdlr.Point2D],
                 triangles: List[Tuple[int, int, int]],
                 edges: List[Tuple[int, int]]=None,
                 name: str=''):
        DisplayMesh.__init__(self, points, triangles, edges, name=name)



class DisplayMesh3D(DisplayMesh):
    _linesegment_class = volmdlr.edges.LineSegment3D
    _point_class = volmdlr.Point3D

    def __init__(self, points: List[volmdlr.Point3D],
                 triangles: List[Tuple[int, int, int]], name=''):
        DisplayMesh.__init__(self, points, triangles)

    def to_babylon(self):
        """
        return mesh in babylon format: https://doc.babylonjs.com/how_to/custom
        """
        positions = []
        for p in self.points:
            positions.extend([k for k in round(p, 6)])

        flatten_indices = []
        for i in self.triangles:
            flatten_indices.extend(i)
        return positions, flatten_indices