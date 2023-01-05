#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base classes.
"""

import os
import subprocess
import tempfile
import webbrowser
from datetime import datetime
from typing import List

import matplotlib.pyplot as plt
import numpy as npy

import dessia_common.core as dc
import dessia_common.files as dcf
import volmdlr
import volmdlr.templates

npy.seterr(divide='raise')

DEFAULT_COLOR = (0.8, 0.8, 0.8)

# TODO: put voldmlr metadata in this freecad header
STEP_HEADER = '''ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('{name}'),'2;1');
FILE_NAME('{filename}','{timestamp}',('Author'),(''),'Volmdlr v{version}','','Unknown');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN {{ 1 0 10303 214 1 1 1 1 }}'));
ENDSEC;
DATA;
#1 = APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#2);
#2 = APPLICATION_CONTEXT('core data for automotive mechanical design processes');
#3 = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );
#4 = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );
#5 = ( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() );
#6 = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(1.E-07),#3,'distance_accuracy_value','confusion accuracy');
#7 = ( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#6)) GLOBAL_UNIT_ASSIGNED_CONTEXT ((#3,#4,#5)) REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY') );
'''

STEP_FOOTER = '''ENDSEC;
END-ISO-10303-21;
'''


def determinant(vec1, vec2, vec3):
    """
    Calculates the determinant for a three vector matrix.

    """
    a = npy.array((vec1.vector, vec2.vector, vec3.vector))
    return npy.linalg.det(a)


def delete_double_point(list_point):
    """
    Delete duplicate points from a list of points.

    :param list_point: The initial list of points
    :type list_point: Union[List[:class:`volmdlr.Point2D`],
        List[:class:`volmdlr.Point3D`]]
    :return: The final list of points containing no duplicates
    :rtype: Union[List[:class:`volmdlr.Point2D`],
        List[:class:`volmdlr.Point3D`]]
    """
    # TODO : this method would be faster using sets
    points = []
    for pt in list_point:
        if pt not in points:
            points.append(pt)
        else:
            continue
    return points


def step_ids_to_str(ids):
    """
    Returns a string with a '#' in front of each ID and a comma separating
    eachone.

    :param ids: A list of step primitives IDs
    :type ids: List[int]
    :return: A string containing all the IDs
    :rtype: str
    """
    return ','.join(['#{}'.format(i) for i in ids])


class CompositePrimitive(dc.PhysicalObject):
    """
    A collection of simple primitives.

    :param name: The name of the collection of primitives.
    :type name: str
    """

    def __init__(self, primitives, name=''):
        self.primitives = primitives
        self.name = name
        self._primitives_to_index = None
        self._utd_primitives_to_index = False
        self.basis_primitives = []

        dc.PhysicalObject.__init__(self, name=name)

    def primitive_to_index(self, primitive):
        if not self._utd_primitives_to_index:
            self._primitives_to_index = {p: ip for ip, p in enumerate(self.primitives)}
            self._utd_primitives_to_index = True
        return self._primitives_to_index[primitive]

    def update_basis_primitives(self):
        basis_primitives = []
        for primitive in self.primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.basis_primitives)
            else:
                basis_primitives.append(primitive)

        self.basis_primitives = basis_primitives


class Primitive2D(dc.PhysicalObject):
    """
    Abstract class for 2D primitives.

    :param name: The name of the 2D primitive.
    :type name: str
    """

    def __init__(self, name=''):
        self.name = name

        dc.PhysicalObject.__init__(self, name=name)


class CompositePrimitive2D(CompositePrimitive):
    """
    A collection of simple 2D primitives.

    :param name: The name of the collection of 2D primitives.
    :type name: str
    """
    _non_serializable_attributes = ['name', '_utd_primitives_to_index',
                                    '_primitives_to_index']
    _non_data_hash_attributes = ['name', '_utd_primitives_to_index',
                                 '_primitives_to_index']

    def __init__(self, primitives, name=''):
        CompositePrimitive.__init__(self, primitives, name=name)
        # self.primitives = primitives
        self.update_basis_primitives()

        self._utd_primitives_to_index = False

    # def primitive_to_index(self, primitive):
    #     if not self._utd_primitives_to_index:
    #         self._primitives_to_index = {p: ip for ip, p in enumerate(self.primitives)}
    #         self._utd_primitives_to_index = True
    #     return self._primitives_to_index[primitive]

    # def update_basis_primitives(self):
    #     basis_primitives = []
    #     for primitive in self.primitives:
    #         if hasattr(primitive, 'basis_primitives'):
    #             basis_primitives.extend(primitive.basis_primitives)
    #         else:
    #             basis_primitives.append(primitive)

    #     self.basis_primitives = basis_primitives

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Rotates the CompositePrimitive2D.

        :param center: rotation center
        :param angle: angle rotation
        :return: a new rotated CompositePrimitive2D
        """
        return self.__class__([point.rotation(center, angle)
                               for point in self.primitives])

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        Rotates the CompositePrimitive2D. Object is updated inplace.

        :param center: rotation center
        :param angle: rotation angle
        """
        primitives = []
        for primitive in self.primitives:
            primitives.append(primitive.rotation(center, angle))
        self.primitives = primitives
        self.update_basis_primitives()

    def translation(self, offset: volmdlr.Vector2D):
        """
        Translates the CompositePrimitive2D.

        :param offset: translation vector
        :return: A new translated CompositePrimitive2D
        """
        return self.__class__([primitive.translation(offset)
                               for primitive in self.primitives])

    def translation_inplace(self, offset: volmdlr.Vector2D):
        """
        Translates the CompositePrimitive2D. Object is updated inplace.

        :param offset: translation vector
        """
        primitives = []
        for primitive in self.primitives:
            primitives.append(primitive.translation(offset))
        self.primitives = primitives
        self.update_basis_primitives()

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and return a new CompositePrimitive2D
        side = 'old' or 'new'
        """
        return self.__class__([primitive.frame_mapping(frame, side)
                               for primitive in self.primitives])

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and the object is updated inplace
        side = 'old' or 'new'
        """
        primitives = []
        for primitive in self.primitives:
            primitives.append(primitive.frame_mapping(frame, side))
        self.primitives = primitives
        self.update_basis_primitives()

    def plot(self, ax=None, color='k', alpha=1,
             plot_points=False, equal_aspect=True):

        if ax is None:
            _, ax = plt.subplots()

        if equal_aspect:
            ax.set_aspect('equal')

        for element in self.primitives:
            element.plot(ax=ax, color=color, alpha=alpha)  # , plot_points=plot_points)

        ax.margins(0.1)
        plt.show()

        return ax

    def plot_data(self, name, fill=None, color='black', stroke_width=0.2,
                  opacity=1):
        plot_data = {}
        plot_data['fill'] = fill
        plot_data['name'] = name
        plot_data['type'] = 'wire'
        plot_data['plot_data'] = []
        for item in self.primitives:
            plot_data['plot_data'].append(item.plot_data(color=color,
                                                         stroke_width=stroke_width,
                                                         opacity=opacity))
        return plot_data


class Primitive3D(dc.PhysicalObject):
    """
    Defines a Primitive3D.
    """

    def __init__(self, color=None, alpha=1, name=''):
        self.color = color
        self.alpha = alpha

        dc.PhysicalObject.__init__(self, name=name)

    def volmdlr_primitives(self):
        return [self]

    def babylon_param(self):
        babylon_param = {'alpha': self.alpha,
                         'name': self.name,
                         }
        if self.color is None:
            babylon_param['color'] = [0.8, 0.8, 0.8]
        else:
            babylon_param['color'] = list(self.color)

        return babylon_param

    def triangulation(self):
        raise NotImplementedError(
            'triangulation method should be implemented on class {}'.format(
                self.__class__.__name__))

    def babylon_meshes(self):
        mesh = self.triangulation()
        if mesh is None:
            return []
        positions, indices = mesh.to_babylon()

        babylon_mesh = {'positions': positions,
                        'indices': indices
                        }
        babylon_mesh.update(self.babylon_param())
        return [babylon_mesh]


class CompositePrimitive3D(CompositePrimitive, Primitive3D):
    """
    A collection of simple primitives3D
    """
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = ['basis_primitives']
    _non_data_eq_attributes = ['name', 'basis_primitives']
    _non_data_hash_attributes = []

    def __init__(self, primitives: List[Primitive3D], color=None, alpha=1, name: str = ''):
        CompositePrimitive.__init__(self, primitives=primitives, name=name)
        Primitive3D.__init__(self, color=color, alpha=alpha, name=name)
        self._utd_primitives_to_index = False

    def plot(self, ax=None, color='k', alpha=1, edge_details=False):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        for primitive in self.primitives:
            primitive.plot(ax=ax, color=color, alpha=alpha)
        return ax

    def babylon_points(self):
        points = []
        if hasattr(self, 'primitives') and \
                hasattr(self.primitives[0], 'start'):
            points = [[self.primitives[0].start.x,
                       self.primitives[0].start.y,
                       self.primitives[0].start.z]]
            for primitive in self.primitives:
                if hasattr(primitive, 'curve'):
                    points.extend(primitive.curve.evalpts)
                else:
                    points.append([primitive.end.x, primitive.end.y,
                                   primitive.end.z])
        elif hasattr(self, 'curve'):
            points = self.curve.evalpts
        elif hasattr(self, 'polygon_points'):
            points = self.polygon_points(50)
        return points

    def babylon_lines(self, points=None):
        if points is None:
            points = self.babylon_points()
        babylon_lines = {'points': points}
        babylon_lines.update(self.babylon_param())
        return [babylon_lines]

    def babylon_curves(self):
        points = self.babylon_points()
        babylon_curves = self.babylon_lines(points)[0]
        return babylon_curves


class BoundingRectangle(dc.DessiaObject):
    """
    Bounding rectangle.

    :param xmin: minimal x coordinate
    :type xmin: float
    :param xmax: maximal x coordinate
    :type xmax: float
    :param ymin: minimal y coordinate
    :type ymin: float
    :param ymax: maximal y coordinate
    :type ymax: float
    """

    def __init__(self, xmin: float, xmax: float, ymin: float, ymax: float, name: str = ''):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        dc.DessiaObject.__init__(self, name=name)

    def __getitem__(self, key):
        if key == 0:
            return self.xmin
        if key == 1:
            return self.xmax
        if key == 2:
            return self.ymin
        if key == 3:
            return self.ymax
        raise IndexError

    def bounds(self):
        """
        Return the bounds of the BoundingRectangle.
        """
        return self.xmin, self.xmax, self.ymin, self.ymax

    def plot(self, ax=None, color='k', linestyle='dotted'):
        """
        Plot of the bounding rectangle and its vertex.
        """

        if not ax:
            _, ax = plt.subplots()
        x = [self.xmin, self.xmax, self.xmax, self.xmin, self.xmin]
        y = [self.ymin, self.ymin, self.ymax, self.ymax, self.ymin]

        ax.plot(x, y, color=color, linestyle=linestyle)
        ax.scatter(x, y, color=color)
        return ax

    def area(self):
        """
        Calculates the area of the bounding rectangle.
        """
        return (self.xmax - self.xmin) * (self.ymax - self.ymin)

    def center(self):
        """
        Calculates the bounding rectangle center.
        """
        return volmdlr.Point2D(0.5 * (self.xmin + self.xmax), 0.5 * (self.ymin + self.ymax))

    def b_rectangle_intersection(self, b_rectangle2):
        """
        Returns True if there is an intersection with another specified bounding rectangle or False otherwise.

        :param b_rectangle2: bounding rectangle to verify intersection
        :type b_rectangle2: :class:`BoundingRectangle`
        """
        return self.xmin < b_rectangle2.xmax and self.xmax > b_rectangle2.xmin \
            and self.ymin < b_rectangle2.ymax and self.ymax > b_rectangle2.ymin

    def is_inside_b_rectangle(self, b_rectangle2, tol: float = 1e-6):
        """
        Returns True if the bounding rectangle is totally inside specified bounding rectangle and False otherwise.

        :param b_rectangle2: A bounding rectangle
        :type b_rectangle2: :class:`BoundingRectangle`
        :param tol: A tolerance for considering inside
        :type tol: float
        """
        return (self.xmin >= b_rectangle2.xmin - tol) and (self.xmax <= b_rectangle2.xmax + tol) \
            and (self.ymin >= b_rectangle2.ymin - tol) and (self.ymax <= b_rectangle2.ymax + tol)

    def point_belongs(self, point: volmdlr.Point2D):
        """
        Returns True if a specified point is inside the bounding rectangle and False otherwise.

        :param point: A 2 dimensional point
        :type point: :class:`volmdlr.Point2D`
        """
        return self.xmin < point.x < self.xmax and self.ymin < point.y < self.ymax

    def intersection_area(self, b_rectangle2):
        """
        Calculates the intersection area between two bounding rectangle.

        :param b_rectangle2: A bounding rectangle
        :type b_rectangle2: :class:`BoundingRectangle`
        """
        if not self.b_rectangle_intersection(b_rectangle2):
            return 0
        if self.is_inside_b_rectangle(b_rectangle2) or b_rectangle2.is_inside_b_rectangle(self):
            return min(self.area(), b_rectangle2.area())

        lx = min(self.xmax, b_rectangle2.xmax) - max(self.xmin, b_rectangle2.xmin)
        ly = min(self.ymax, b_rectangle2.ymax) - max(self.ymin, b_rectangle2.ymin)

        return lx * ly

    def distance_to_b_rectangle(self, b_rectangle2):
        """
        Calculates the minimal distance between two bounding rectangles.

        :param b_rectangle2: A bounding rectangle
        :type b_rectangle2: :class:`BoundingRectangle`
        """
        if self.b_rectangle_intersection(b_rectangle2):
            return 0

        permute_b_rec1 = self
        permute_b_rec2 = b_rectangle2

        if permute_b_rec2.xmin < permute_b_rec1.xmin:
            permute_b_rec1, permute_b_rec2 = permute_b_rec2, permute_b_rec1
        dx = max(permute_b_rec2.xmin - permute_b_rec1.xmax, 0)

        if permute_b_rec2.ymin < permute_b_rec1.ymin:
            permute_b_rec1, permute_b_rec2 = permute_b_rec2, permute_b_rec1
        dy = max(permute_b_rec2.ymin - permute_b_rec1.ymax, 0)

        return (dx ** 2 + dy ** 2) ** 0.5

    def distance_to_point(self, point: volmdlr.Point2D):
        """
        Calculate the minimal distance between the bounding rectangle and
        a specified point.

        :param point: A 2 dimensional point
        :type point: :class:`volmdlr.Point2D`
        """
        if self.point_belongs(point):
            return min([self.xmax - point.x, point.y - self.xmin,
                        self.ymax - point.y, point.y - self.ymin])

        if point.x < self.xmin:
            dx = self.xmin - point.x
        elif self.xmax < point.x:
            dx = point.x - self.xmax
        else:
            dx = 0

        if point.y < self.ymin:
            dy = self.ymin - point.y
        elif self.ymax < point.y:
            dy = point.y - self.ymax
        else:
            dy = 0

        return (dx ** 2 + dy ** 2) ** 0.5


class BoundingBox(dc.DessiaObject):
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

        dc.DessiaObject.__init__(self, name=name)

        self.center = volmdlr.Point3D(0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0.5 * (zmin + zmax))

    def __hash__(self):
        return sum(hash(point) for point in self.points)

    def __add__(self, other_bbox):
        return BoundingBox(min(self.xmin, other_bbox.xmin),
                           max(self.xmax, other_bbox.xmax),
                           min(self.ymin, other_bbox.ymin),
                           max(self.ymax, other_bbox.ymax),
                           min(self.zmin, other_bbox.zmin),
                           max(self.zmax, other_bbox.zmax))

    def to_dict(self, use_pointers: bool = True, memo=None, path: str = '#'):
        """
        Converts the bounding box to a dict.

        :param use_pointers: DESCRIPTION, defaults to True
        :type use_pointers: bool, optional
        :param memo: DESCRIPTION, defaults to None
        :type memo: TYPE, optional
        :param path: DESCRIPTION, defaults to '#'
        :type path: str, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """
        return {'object_class': 'volmdlr.core.BoundingBox',
                'name': self.name,
                'xmin': self.xmin,
                'xmax': self.xmax,
                'ymin': self.ymin,
                'ymax': self.ymax,
                'zmin': self.zmin,
                'zmax': self.zmax,
                }

    @property
    def points(self):
        return [volmdlr.Point3D(self.xmin, self.ymin, self.zmin),
                volmdlr.Point3D(self.xmax, self.ymin, self.zmin),
                volmdlr.Point3D(self.xmax, self.ymax, self.zmin),
                volmdlr.Point3D(self.xmin, self.ymax, self.zmin),
                volmdlr.Point3D(self.xmin, self.ymin, self.zmax),
                volmdlr.Point3D(self.xmax, self.ymin, self.zmax),
                volmdlr.Point3D(self.xmax, self.ymax, self.zmax),
                volmdlr.Point3D(self.xmin, self.ymax, self.zmax)]

    def plot(self, ax=None, color='gray'):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        bbox_edges = [[self.points[0], self.points[1]],
                      [self.points[0], self.points[3]],
                      [self.points[0], self.points[4]],
                      [self.points[1], self.points[2]],
                      [self.points[1], self.points[5]],
                      [self.points[2], self.points[3]],
                      [self.points[2], self.points[6]],
                      [self.points[3], self.points[7]],
                      [self.points[4], self.points[5]],
                      [self.points[5], self.points[6]],
                      [self.points[6], self.points[7]],
                      [self.points[7], self.points[4]]]

        # x = [p[0] for p in self.points]
        # y = [p[1] for p in self.points]
        # z = [p[2] for p in self.points]
        # ax.scatter(x, y, z, color)
        for edge in bbox_edges:
            ax.plot3D([edge[0][0], edge[1][0]],
                      [edge[0][1], edge[1][1]],
                      [edge[0][2], edge[1][2]],
                      color=color)
        ax.set_xlabel('X Label')
        ax.set_ylabel('Y Label')
        ax.set_zlabel('Z Label')
        return ax

    @classmethod
    def from_bounding_boxes(cls, bounding_boxes):
        xmin = min(bb.xmin for bb in bounding_boxes)
        xmax = max(bb.xmax for bb in bounding_boxes)
        ymin = min(bb.ymin for bb in bounding_boxes)
        ymax = max(bb.ymax for bb in bounding_boxes)
        zmin = min(bb.zmin for bb in bounding_boxes)
        zmax = max(bb.zmax for bb in bounding_boxes)
        return cls(xmin, xmax, ymin, ymax, zmin, zmax)

    @classmethod
    def from_points(cls, points):
        # if len(points) == 0:
        #     return (0, 0, 0, 0, 0, 0)
        xmin = min(pt.x for pt in points)
        xmax = max(pt.x for pt in points)
        ymin = min(pt.y for pt in points)
        ymax = max(pt.y for pt in points)
        zmin = min(pt.z for pt in points)
        zmax = max(pt.z for pt in points)
        return cls(xmin, xmax, ymin, ymax, zmin, zmax)

    def to_frame(self):
        x = volmdlr.Vector3D((self.xmax - self.xmin), 0, 0)
        y = volmdlr.Vector3D(0, (self.ymax - self.ymin), 0)
        z = volmdlr.Vector3D(0, 0, (self.zmax - self.zmin))
        return volmdlr.Frame3D(self.center, x, y, z)

    def volume(self):
        return (self.xmax - self.xmin) * (self.ymax - self.ymin) * (
                self.zmax - self.zmin)

    def bbox_intersection(self, bbox2):
        if self.xmin < bbox2.xmax and self.xmax > bbox2.xmin:
            if self.ymin < bbox2.ymax and self.ymax > bbox2.ymin\
                    and self.zmin < bbox2.zmax and self.zmax > bbox2.zmin:
                return True
        if self.xmin == bbox2.xmax and self.xmax == bbox2.xmin:
            if self.ymin < bbox2.ymax and self.ymax > bbox2.ymin \
                    and self.zmin < bbox2.zmax and self.zmax > bbox2.zmin:
                return True
        return False

    def is_inside_bbox(self, bbox2):
        return (self.xmin >= bbox2.xmin - 1e-6) and (self.xmax <= bbox2.xmax + 1e-6) \
               and (self.ymin >= bbox2.ymin - 1e-6) and (self.ymax <= bbox2.ymax + 1e-6) \
               and (self.zmin >= bbox2.zmin - 1e-6) and (self.zmax <= bbox2.zmax + 1e-6)

    def intersection_volume(self, bbox2):
        if not self.bbox_intersection(bbox2):
            return 0
        if self.is_inside_bbox(bbox2) or bbox2.is_inside_bbox(self):
            return min(self.volume(), bbox2.volume())

        lx = min(self.xmax, bbox2.xmax) - max(self.xmin, bbox2.xmin)
        ly = min(self.ymax, bbox2.ymax) - max(self.ymin, bbox2.ymin)
        lz = min(self.zmax, bbox2.zmax) - max(self.zmin, bbox2.zmin)

        return lx * ly * lz

    # def intersection_volume(self, bbox2):
    #     if not self.bbox_intersection(bbox2):
    #         return 0
    #
    #     permute_bbox1 = self
    #     permute_bbox2 = bbox2
    #
    #     if permute_bbox2.xmin < permute_bbox1.xmin:
    #         permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
    #     lx = permute_bbox1.xmax - permute_bbox2.xmin
    #
    #     if permute_bbox2.ymin < permute_bbox1.ymin:
    #         permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
    #     ly = permute_bbox1.ymax - permute_bbox2.ymin
    #
    #     if permute_bbox2.zmin < permute_bbox1.zmin:
    #         permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
    #     lz = permute_bbox1.zmax - permute_bbox2.zmin
    #
    #     return lx*ly*lz

    def distance_to_bbox(self, bbox2):
        if self.bbox_intersection(bbox2):
            return 0

        permute_bbox1 = self
        permute_bbox2 = bbox2

        if permute_bbox2.xmin < permute_bbox1.xmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dx = max(permute_bbox2.xmin - permute_bbox1.xmax, 0)

        if permute_bbox2.ymin < permute_bbox1.ymin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dy = max(permute_bbox2.ymin - permute_bbox1.ymax, 0)

        if permute_bbox2.zmin < permute_bbox1.zmin:
            permute_bbox1, permute_bbox2 = permute_bbox2, permute_bbox1
        dz = max(permute_bbox2.zmin - permute_bbox1.zmax, 0)

        return (dx ** 2 + dy ** 2 + dz ** 2) ** 0.5

    def point_belongs(self, point):
        return self.xmin < point[0] < self.xmax \
               and self.ymin < point[1] < self.ymax \
               and self.zmin < point[2] < self.zmax

    def distance_to_point(self, point):
        if self.point_belongs(point):
            return min([self.xmax - point[0], point[0] - self.xmin,
                        self.ymax - point[1], point[1] - self.ymin,
                        self.zmax - point[2], point[2] - self.zmin])

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
        return (dx ** 2 + dy ** 2 + dz ** 2) ** 0.5


class VolumeModel(dc.PhysicalObject):
    """
    A class containing one or several :class:`volmdlr.core.Primitive3D`.

    :param primitives: The vector's abscissa
    :type primitives: List[:class:`volmdlr.core.Primitive3D`]
    :param name: The VolumeModel's name
    :type name: str
    """
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = ['shells', 'bounding_box']
    _non_data_eq_attributes = ['name', 'shells', 'bounding_box', 'contours',
                               'faces']
    _non_data_hash_attributes = ['name', 'shells', 'bounding_box', 'contours',
                                 'faces']
    _dessia_methods = ['to_stl_model']

    def __init__(self, primitives: List[Primitive3D], name: str = ''):
        self.primitives = primitives
        self.name = name
        self.shells = []
        self._bbox = None
        dc.PhysicalObject.__init__(self, name=name)

    def __hash__(self):
        return sum(hash(point) for point in self.primitives)

    # def _extract_shells(self):
    #     shells = []
    #     for primitive in self.primitives:
    #         if isinstance(primitive, Shell3D):
    #             shells.append(primitive)
    #     return shells

    def __eq__(self, other):
        if self.__class__.__name__ != other.__class__.__name__:
            return False
        equ = True
        if len(self.primitives) != len(other.primitives):
            return False
        for p1, p2 in zip(self.primitives, other.primitives):
            equ = equ and p1 == p2
        return equ

    @property
    def bounding_box(self):
        """
        Returns the bounding box.
        """
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def _bounding_box(self) -> BoundingBox:
        """
        Computes the bounding box of the model.
        """
        return BoundingBox.from_bounding_boxes([p.bounding_box for p in self.primitives])

    def volume(self) -> float:
        """
        Return the sum of volumes of the primitives.

        It does not make any boolean operation in case of overlaping.
        """
        volume = 0
        for primitive in self.primitives:
            volume += primitive.volume()
        return volume

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Rotates the VolumeModel.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated VolumeModel
        """
        new_primitives = [
            primitive.rotation(center, axis, angle) for
            primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        Rotates the VolumeModel. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        for primitive in self.primitives:
            primitive.rotation_inplace(center, axis, angle)
        self.bounding_box = self._bounding_box()

    def translation(self, offset: volmdlr.Vector3D):
        """
        Translates the VolumeModel.

        :param offset: translation vector
        :return: A new translated VolumeModel
        """
        new_primitives = [primitive.translation(offset) for
                          primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Translates the VolumeModel. Object is updated inplace.

        :param offset: translation vector
        """
        for primitives in self.primitives:
            primitives.translation_inplace(offset)
        self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new VolumeModel
        side = 'old' or 'new'
        """
        new_primitives = [primitive.frame_mapping(frame, side)
                          for primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace
        side = 'old' or 'new'
        """
        for primitives in self.primitives:
            primitives.frame_mapping_inplace(frame, side)
        self.bounding_box = self._bounding_box()

    def copy(self, deep=True, memo=None):
        """
        Specific copy
        """
        new_primitives = [primitive.copy(deep=deep, memo=memo) for primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def plot2d(self, ax=None, color=None):
        fig = plt.figure()
        if ax is None:
            ax = fig.add_subplot(111, projection='3d')

        for i, shell in enumerate(self.shells):
            bbox = shell.bbox()
            bbox.plot(ax, color[i])

        return ax

    def plot(self, equal_aspect=True):
        """
        Matplotlib plot of model.
        To use for debug.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', adjustable='box')
        for primitive in self.primitives:
            primitive.plot(ax)
        if equal_aspect:
            # ax.set_aspect('equal')
            ax.set_aspect('auto')
        ax.margins(0.1)
        return ax

    def freecad_script(self, fcstd_filepath,
                       freecad_lib_path='/usr/lib/freecad/lib',
                       export_types=('fcstd',),
                       save_to='',
                       tolerance=0.0001):
        """
        Generates python a FreeCAD definition of model.

        :param fcstd_filename: a filename without extension to give the name at the fcstd part written in python code
        :type fcstd_filename:str
        """
        fcstd_filepath = os.path.abspath(fcstd_filepath)
        fcstd_filepath = fcstd_filepath.replace('\\', '\\\\')
        freecad_lib_path = freecad_lib_path.replace('\\', '\\\\')

        s = '# -*- coding: utf-8 -*-\n'
        if freecad_lib_path != '':
            s += "import sys\nsys.path.append('" + freecad_lib_path + "')\n"

        s += "import math\nimport FreeCAD as fc\nimport Part\n\ndoc=fc.newDocument('doc')\n\n"
        for ip, primitive in enumerate(self.primitives):
            if primitive.name == '':
                primitive_name = f'Primitive_{ip}'
            else:
                primitive_name = f'Primitive_{ip}_{primitive.name}'
            s += f"part = doc.addObject('App::Part','{primitive_name}')\n"
            if hasattr(primitive, 'FreeCADExport'):
                sp = primitive.FreeCADExport(ip)
                if sp != '':
                    #                        s += (sp+'\n')
                    s += (sp)
                    s += f'shapeobj = doc.addObject("Part::Feature","{primitive_name}")\n'
                    # if isinstance(primitive, BSplineCurve3D) \
                    #         or isinstance(primitive, BSplineSurface3D) \
                    #         or isinstance(primitive, Circle3D) \
                    #         or isinstance(primitive, LineSegment3D) \
                    #         or isinstance(primitive, Ellipse3D):
                    #     #                            print(primitive)
                    #     #                            s += 'S = Part.Shape([primitive{}])\n'.format(ip)
                    #     #                            s += 'shapeobj.Shape = S\n'
                    #     s += 'shapeobj.Shape = primitive{}.toShape()\n'.format(
                    #         ip)
                    # else:
                    s += f"shapeobj.Shape = primitive{ip}\n"
                    s += 'part.addObject(shapeobj)\n\n'
            # --------------------DEBUG-------------------
        #                else:
        #                    raise NotImplementedError
        # ---------------------------------------------

        s += 'doc.recompute()\n'
        if 'fcstd' in export_types:
            s += "doc.saveAs('" + fcstd_filepath + ".fcstd')\n\n"
        if 'stl' in export_types:
            s += f"import Mesh\nMesh.export(doc.Objects,'{fcstd_filepath}.stl', tolerance={tolerance})\n"
        if 'step' in export_types:
            s += f"Part.export(doc.Objects,'{fcstd_filepath}.step')\n"

        if save_to != '':
            with open(os.path.abspath(save_to), 'w') as file:
                file.write(s)
        return s

    def freecad_export(self, fcstd_filepath,
                       python_path='python3',
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
        fcstd_filepath = os.path.abspath(fcstd_filepath)
        s = self.freecad_script(fcstd_filepath,
                                freecad_lib_path=freecad_lib_path,
                                export_types=export_types,
                                tolerance=tolerance)
        with tempfile.NamedTemporaryFile(suffix=".py", delete=False) as f:
            f.write(bytes(s, 'utf8'))

        arg = f.name
        output = subprocess.call([python_path, arg])

        f.close()
        os.remove(f.name)
        return output

    def babylon_data(self):
        meshes = []
        lines = []
        for primitive in self.primitives:
            if hasattr(primitive, 'babylon_meshes'):
                meshes.extend(primitive.babylon_meshes())
            if hasattr(primitive, 'babylon_curves'):
                lines.append(primitive.babylon_curves())

        bbox = self.bounding_box
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])

        babylon_data = {'meshes': meshes,
                        'lines': lines,
                        'max_length': max_length,
                        'center': list(center)}

        return babylon_data

    @classmethod
    def babylonjs_script(cls, babylon_data, use_cdn=True,
                         debug=False):
        if use_cdn:
            script = volmdlr.templates.babylon_unpacker_cdn_header  # .substitute(name=page_name)
        else:
            script = volmdlr.templates.babylon_unpacker_embedded_header  # .substitute(name=page_name)

        script += volmdlr.templates.babylon_unpacker_body_template.substitute(
            babylon_data=babylon_data)
        return script

    def babylonjs(self, page_name=None, use_cdn=True, debug=False):
        babylon_data = self.babylon_data()
        script = self.babylonjs_script(babylon_data, use_cdn=use_cdn,
                                       debug=debug)
        if page_name is None:
            with tempfile.NamedTemporaryFile(suffix=".html",
                                             delete=False) as file:
                file.write(bytes(script, 'utf8'))
            page_name = file.name
        else:
            if not page_name.endswith('.html'):
                page_name += '.html'
            with open(page_name, 'w') as file:
                file.write(script)

        webbrowser.open('file://' + os.path.realpath(page_name))

        return page_name

    def save_babylonjs_to_file(self, filename: str = None,
                               use_cdn=True, debug=False):
        babylon_data = self.babylon_data()
        script = self.babylonjs_script(babylon_data, use_cdn=use_cdn,
                                       debug=debug)
        if filename is None:
            with tempfile.NamedTemporaryFile(suffix=".html",
                                             delete=False) as file:
                file.write(bytes(script, 'utf8'))
                return file.name

        if not filename.endswith('.html'):
            filename += '.html'

            with open(filename, 'w') as file:
                file.write(script)
            return filename

    def to_stl_model(self):
        mesh = self.primitives[0].triangulation()
        for primitive in self.primitives[1:]:
            mesh.merge_mesh(primitive.triangulation())
        stl = mesh.to_stl()
        return stl

    def to_stl(self, filepath: str):
        if not filepath.endswith('.stl'):
            filepath += '.stl'
        with open(filepath, 'wb') as file:
            self.to_stl_stream(file)

    def to_stl_stream(self, stream: dcf.BinaryFile):
        stl = self.to_stl_model()
        stl.save_to_stream(stream)
        return stream

    def to_step(self, filepath: str):
        if not (filepath.endswith('.step') or filepath.endswith('.stp')):
            filepath += '.step'
        with open(filepath, 'w') as file:
            self.to_step_stream(file)

    def to_step_stream(self, stream: dcf.StringFile):

        step_content = STEP_HEADER.format(name=self.name,
                                          filename='',
                                          timestamp=datetime.now().isoformat(),
                                          version=volmdlr.__version__)
        current_id = 8

        for primitive in self.primitives:
            if primitive.__class__.__name__ == 'OpenShell3D':
                primitive_content, primitive_id, face_ids = primitive.to_step_face_ids(current_id)
            else:
                primitive_content, primitive_id = primitive.to_step(current_id)

            step_content += primitive_content

            product_definition_context_id = primitive_id + 1
            step_content += f"#{product_definition_context_id} = PRODUCT_DEFINITION_CONTEXT('part definition',#2,'design');\n"

            product_context_id = product_definition_context_id + 1
            step_content += f"#{product_context_id} = PRODUCT_CONTEXT('',#2,'mechanical');\n"
            product_id = product_context_id + 1
            step_content += f"#{product_id} = PRODUCT('{primitive.name}'," \
                            f"'{primitive.name}','',(#{product_context_id}));\n"
            product_definition_formation_id = product_id + 1
            step_content += f"#{product_definition_formation_id} = PRODUCT_DEFINITION_FORMATION('','',#{product_id});\n"
            product_definition_id = product_definition_formation_id + 1
            step_content += f"#{product_definition_id} = PRODUCT_DEFINITION('design'," \
                            f"'',#{product_definition_formation_id},#{product_definition_context_id});\n"
            product_definition_shape_id = product_definition_id + 1
            step_content += f"#{product_definition_shape_id} = PRODUCT_DEFINITION_SHAPE(''," \
                            f"'',#{product_definition_id});\n"
            shape_definition_repr_id = product_definition_shape_id + 1
            step_content += "#{} = SHAPE_DEFINITION_REPRESENTATION(#{},#{});\n".format(shape_definition_repr_id,
                                                                                       product_definition_shape_id,
                                                                                       primitive_id
                                                                                       )
            product_related_category = shape_definition_repr_id + 1
            step_content += f"#{product_related_category} = PRODUCT_RELATED_PRODUCT_CATEGORY(" \
                            f"'part',$,(#{product_id}));\n"
            draughting_id = product_related_category + 1
            step_content += f"#{draughting_id} = DRAUGHTING_PRE_DEFINED_CURVE_FONT('continuous');\n"
            color_id = draughting_id + 1
            primitive_color = (1, 1, 1)
            if hasattr(primitive, 'color') and primitive.color is not None:
                primitive_color = primitive.color
            step_content += f"#{color_id} = COLOUR_RGB('',{round(float(primitive_color[0]), 4)}," \
                            f"{round(float(primitive_color[1]), 4)}, {round(float(primitive_color[2]), 4)});\n"

            curve_style_id = color_id + 1
            step_content += f"#{curve_style_id} = CURVE_STYLE('',#{draughting_id}," \
                            f"POSITIVE_LENGTH_MEASURE(0.1),#{color_id});\n"

            fill_area_color_id = curve_style_id + 1
            step_content += f"#{fill_area_color_id} = FILL_AREA_STYLE_COLOUR('',#{color_id});\n"

            fill_area_id = fill_area_color_id + 1
            step_content += f"#{fill_area_id} = FILL_AREA_STYLE('',#{fill_area_color_id});\n"

            suface_fill_area_id = fill_area_id + 1
            step_content += f"#{suface_fill_area_id} = SURFACE_STYLE_FILL_AREA(#{fill_area_id});\n"

            suface_side_style_id = suface_fill_area_id + 1
            step_content += f"#{suface_side_style_id} = SURFACE_SIDE_STYLE('',(#{suface_fill_area_id}));\n"

            suface_style_usage_id = suface_side_style_id + 1
            step_content += f"#{suface_style_usage_id} = SURFACE_STYLE_USAGE(.BOTH.,#{suface_side_style_id});\n"

            presentation_style_id = suface_style_usage_id + 1

            step_content += f"#{presentation_style_id} = PRESENTATION_STYLE_ASSIGNMENT((#{suface_style_usage_id}," \
                            f"#{curve_style_id}));\n"

            styled_item_id = presentation_style_id + 1
            if primitive.__class__.__name__ == 'OpenShell3D':
                for face_id in face_ids:
                    step_content += f"#{styled_item_id} = STYLED_ITEM('color',(#{presentation_style_id})," \
                                    f"#{face_id});\n"
                    styled_item_id += 1
                styled_item_id -= 1
            else:
                step_content += f"#{styled_item_id} = STYLED_ITEM('color',(#{presentation_style_id})," \
                                f"#{primitive_id});\n"

            current_id = styled_item_id + 1

        step_content += STEP_FOOTER

        stream.write(step_content)

    def volmdlr_volume_model(self):
        return self


class MovingVolumeModel(VolumeModel):
    """
    A volume model with possibility to declare time steps at which the primitives are positionned with frames.

    """

    def __init__(self, primitives, step_frames, name=''):
        VolumeModel.__init__(self, primitives=primitives, name=name)
        self.step_frames = step_frames

        if not self.is_consistent():
            raise dc.ConsistencyError

    def is_consistent(self):
        n_primitives = len(self.primitives)
        for frames in self.step_frames:
            if len(frames) != n_primitives:
                return False
        return True

    def step_volume_model(self, istep):
        primitives = []
        for primitive, frame in zip(self.primitives, self.step_frames[istep]):
            primitives.append(
                primitive.frame_mapping(frame, side='old'))
        return VolumeModel(primitives)

    def babylon_data(self):
        meshes = []
        primitives_to_meshes = []
        for ip, primitive in enumerate(self.primitives):
            if hasattr(primitive, 'babylon_meshes'):
                meshes.extend(primitive.babylon_meshes())
                primitives_to_meshes.append(ip)

        bbox = self._bounding_box()
        center = bbox.center
        max_length = max([bbox.xmax - bbox.xmin,
                          bbox.ymax - bbox.ymin,
                          bbox.zmax - bbox.zmin])

        steps = []
        for istep, frames in enumerate(self.step_frames):

            # step_positions = []
            # step_orientations = []
            step = {'time': istep}
            for iframe, frame in enumerate(frames):
                if iframe in primitives_to_meshes:
                    imesh = primitives_to_meshes.index(iframe)
                    step[imesh] = {}
                    step[imesh]['position'] = list(round(frame.origin, 6))
                    step[imesh]['orientations'] = [list(round(frame.u, 6)),
                                                   list(round(frame.v, 6)),
                                                   list(round(frame.w, 6))]

            steps.append(step)

        babylon_data = {'meshes': meshes,
                        'max_length': max_length,
                        'center': list(center),
                        'steps': steps}
        return babylon_data
