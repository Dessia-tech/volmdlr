#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base classes.
"""

import os
import tempfile
import warnings
import webbrowser
from dataclasses import dataclass
from datetime import datetime
from functools import lru_cache
from typing import List, Tuple

# import gmsh
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


def point_in_list(point, list_points, tol: float = 1e-6):
    """
    Verifies if a point is inside a list  of points, considering a certain tolerance.

    :param point: Point to be verified inside list.
    :param list_points: List of points to be used.
    :param tol: Tolerance to consider if two points are the same.
    :return: True if there is a point inside the list close to the point to given tolerance.
    """
    for point_i in list_points:
        if point.is_close(point_i, tol):
            return True
    return False


def get_point_index_in_list(point, list_points, tol: float = 1e-6):
    """
    Gets the index a point inside a list of points, considering a certain tolerance.

    :param point: Point to be verified inside list.
    :param list_points: List of points to be used.
    :param tol: Tolerance to consider if two points are the same.
    :return: The point index.
    """
    for i, point_i in enumerate(list_points):
        if point_i.is_close(point, tol):
            return i
    raise ValueError(f'{point} is not in list')


def determinant(vec1, vec2, vec3):
    """
    Calculates the determinant for a three vector matrix.

    """
    # TODO: to be removed
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
    for point in list_point:
        if point not in points:
            points.append(point)
        else:
            continue
    return points


def step_ids_to_str(ids):
    """
    Returns a string with a '#' in front of each ID and a comma separating each-one.

    :param ids: A list of step primitives IDs
    :type ids: List[int]
    :return: A string containing all the IDs
    :rtype: str
    """
    return ','.join([f"#{i}" for i in ids])


@dataclass
class EdgeStyle:
    """
    Data class for styling edges matplotlib plots.

    """
    color: str = 'k'
    alpha: float = 1
    edge_ends: bool = False
    edge_direction: bool = False
    width: float = None
    arrow: bool = False
    plot_points: bool = False
    dashed: bool = True
    linestyle: str = '-'
    linewidth: float = 1
    equal_aspect: bool = True


class CompositePrimitive(dc.PhysicalObject):
    """
    A collection of simple primitives.

    :param name: The name of the collection of primitives.
    :type name: str
    """

    def __init__(self, primitives, name: str = ''):
        self.primitives = primitives
        self.name = name
        self._primitives_to_index = None
        self._utd_primitives_to_index = False
        self.basis_primitives = []

        dc.PhysicalObject.__init__(self, name=name)

    def to_dict(self, *args, **kwargs):
        """Avoids storing points in memo that makes serialization slow."""
        return dc.PhysicalObject.to_dict(self, use_pointers=False)

    def primitive_to_index(self, primitive):
        """Constructs a dictionary associating primitive to its index."""
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

    def __init__(self, name: str = ''):
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

    def __init__(self, primitives: List[Primitive2D], name: str = ''):
        CompositePrimitive.__init__(self, primitives, name=name)
        self.update_basis_primitives()

        self._utd_primitives_to_index = False

    def rotation(self, center: volmdlr.Point2D, angle: float):
        """
        Rotates the CompositePrimitive2D.

        :param center: rotation center.
        :param angle: angle rotation.
        :return: a new rotated CompositePrimitive2D.
        """
        return self.__class__([point.rotation(center, angle)
                               for point in self.primitives])

    def rotation_inplace(self, center: volmdlr.Point2D, angle: float):
        """
        Rotates the CompositePrimitive2D. Object is updated in-place.

        :param center: rotation center.
        :param angle: rotation angle.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

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
        Translates the CompositePrimitive2D. Object is updated in-place.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        primitives = []
        for primitive in self.primitives:
            primitives.append(primitive.translation(offset))
        self.primitives = primitives
        self.update_basis_primitives()

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and return a new CompositePrimitive2D.

        side = 'old' or 'new'
        """
        return self.__class__([primitive.frame_mapping(frame, side)
                               for primitive in self.primitives])

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        primitives = []
        for primitive in self.primitives:
            primitives.append(primitive.frame_mapping(frame, side))
        self.primitives = primitives
        self.update_basis_primitives()

    def plot(self, ax=None, edge_style=EdgeStyle()):

        if ax is None:
            _, ax = plt.subplots()

        if edge_style.equal_aspect:
            ax.set_aspect('equal')

        for element in self.primitives:
            element.plot(ax=ax, edge_style=edge_style)

        ax.margins(0.1)
        plt.show()

        return ax

    def plot_data(self, name, fill=None, color='black', stroke_width=0.2,
                  opacity=1):
        # TODO: not working on 0.8.0
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

    def __init__(self, color: Tuple[float, float, float] = None, alpha: int = 1, name: str = ''):
        self.color = color
        self.alpha = alpha

        dc.PhysicalObject.__init__(self, name=name)

    def volmdlr_primitives(self):
        return [self]

    def babylon_param(self):
        """
        Returns babylonjs parameters.

        :return: babylonjs parameters (alpha, name, color)
        :rtype: dict

        """

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
            f"triangulation method should be implemented on class {self.__class__.__name__}")

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
    A collection of simple primitives3D.
    """
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = ['basis_primitives']
    _non_data_eq_attributes = ['name', 'basis_primitives']
    _non_data_hash_attributes = []

    def __init__(self, primitives: List[Primitive3D], color: Tuple[float, float, float] = None, alpha: float = 1,
                 name: str = ''):
        CompositePrimitive.__init__(self, primitives=primitives, name=name)
        Primitive3D.__init__(self, color=color, alpha=alpha, name=name)
        self._utd_primitives_to_index = False

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        for primitive in self.primitives:
            primitive.plot(ax=ax, edge_style=edge_style)
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
        # Disabling Dessia object init call for performance. Check when performance enhancement on dessia_common side
        # dc.DessiaObject.__init__(self, name=name)
        self.name = name

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

    def bounding_points(self):
        """
        Return the bounds of the BoundingRectangle.
        """
        return [volmdlr.Point2D(self.xmin, self.ymin), volmdlr.Point2D(self.xmax, self.ymin),
                volmdlr.Point2D(self.xmax, self.ymax), volmdlr.Point2D(self.xmin, self.ymax)]

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
        Calculate the minimal distance between the bounding rectangle and a specified point.

        :param point: A 2D point
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

    @classmethod
    def from_points(cls, points: List[volmdlr.Point2D]) -> "BoundingRectangle":
        """
        Initializes a bounding rectangle from a list of points.

        :param points: The list of points to create the bounding rectangle from.
        :type points: List[volmdlr.Point2D]
        :return: The bounding rectangle initialized from the list of points.
        :rtype: BoundingRectangle
        """
        xmin = min(pt.x for pt in points)
        xmax = max(pt.x for pt in points)
        ymin = min(pt.y for pt in points)
        ymax = max(pt.y for pt in points)
        return cls(xmin, xmax, ymin, ymax)


class BoundingBox(dc.DessiaObject):
    """
    An axis aligned boundary box.
    """

    def __init__(self, xmin: float, xmax: float, ymin: float, ymax: float, zmin: float, zmax: float, name: str = ""):
        """
        Initializes a bounding box.

        :param xmin: The x-coordinate of the lower-left corner.
        :type xmin: float
        :param xmax: The x-coordinate of the upper-right corner.
        :type xmax: float
        :param ymin: The y-coordinate of the lower-left corner.
        :type ymin: float
        :param ymax: The y-coordinate of the upper-right corner.
        :type ymax: float
        :param zmin: The z-coordinate of the lower-left corner.
        :type zmin: float
        :param zmax: The z-coordinate of the upper-right corner.
        :type zmax: float
        :param name: The name of the bounding box.
        :type name: str, optional
        """
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax
        self.zmin = zmin
        self.zmax = zmax

        # disabling super init call for efficiency, put back when dc disable kwargs
        # dc.DessiaObject.__init__(self, name=name)
        self.name = name

    @property
    @lru_cache
    def center(self):
        """
        Computes the center of the bounding box.

        TODO: change lru_cache to cached property when support for py3.7 is dropped.
        """
        return volmdlr.Point3D(0.5 * (self.xmin + self.xmax),
                               0.5 * (self.ymin + self.ymax),
                               0.5 * (self.zmin + self.zmax))

    def __hash__(self) -> int:
        return sum(hash(point) for point in self.points)

    def __add__(self, other_bbox) -> "BoundingBox":
        return BoundingBox(min(self.xmin, other_bbox.xmin),
                           max(self.xmax, other_bbox.xmax),
                           min(self.ymin, other_bbox.ymin),
                           max(self.ymax, other_bbox.ymax),
                           min(self.zmin, other_bbox.zmin),
                           max(self.zmax, other_bbox.zmax))

    def to_dict(self, *args, **kwargs) -> dict:
        """
        Converts the bounding box to a dictionary representation.

        :param use_pointers: DESCRIPTION, defaults to True
        :type use_pointers: bool, optional
        :param memo: DESCRIPTION, defaults to None
        :type memo: TYPE, optional
        :param path: A string representing the current position of the object in the serialized data structure.
        :type path: str, optional

        :return: The dictionary representation of the bounding box.
        :rtype: dict
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
    def points(self) -> List[volmdlr.Point3D]:
        """
        Returns the eight corner points of the bounding box.

        :return: A list of eight 3D points representing the corners of the bounding box.
        :rtype: list of volmdlr.Point3D
        """
        return [volmdlr.Point3D(self.xmin, self.ymin, self.zmin),
                volmdlr.Point3D(self.xmax, self.ymin, self.zmin),
                volmdlr.Point3D(self.xmax, self.ymax, self.zmin),
                volmdlr.Point3D(self.xmin, self.ymax, self.zmin),
                volmdlr.Point3D(self.xmin, self.ymin, self.zmax),
                volmdlr.Point3D(self.xmax, self.ymin, self.zmax),
                volmdlr.Point3D(self.xmax, self.ymax, self.zmax),
                volmdlr.Point3D(self.xmin, self.ymax, self.zmax)]

    def plot(self, ax=None, color='gray'):
        """
        Plot the bounding box on 3D axes.

        :param ax: The 3D axes to plot on. If not provided, a new figure will be created.
        :type ax: matplotlib.axes._subplots.Axes3DSubplot, optional
        :param color: The color of the lines used to plot the bounding box.
        :type color: str, optional
        :return: The 3D axes with the plotted bounding box.
        :rtype: matplotlib.axes._subplots.Axes3DSubplot
        """
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
    def from_bounding_boxes(cls, bounding_boxes: List["BoundingBox"]) -> "BoundingBox":
        """
        Creates a bounding box that contains multiple bounding boxes.

        :param bounding_boxes: A list of bounding boxes that need to be contained.
        :type bounding_boxes: list of BoundingBox
        :return: A new bounding box that contains all the input bounding boxes.
        :rtype: BoundingBox
        """
        xmin = min(bb.xmin for bb in bounding_boxes)
        xmax = max(bb.xmax for bb in bounding_boxes)
        ymin = min(bb.ymin for bb in bounding_boxes)
        ymax = max(bb.ymax for bb in bounding_boxes)
        zmin = min(bb.zmin for bb in bounding_boxes)
        zmax = max(bb.zmax for bb in bounding_boxes)
        return cls(xmin, xmax, ymin, ymax, zmin, zmax)

    @classmethod
    def from_points(cls, points: List[volmdlr.Point3D]) -> "BoundingBox":
        """
        Initializes a bounding box from a list of points.

        :param points: The list of points to create the bounding box from.
        :type points: List[volmdlr.Point3D]
        :return: The bounding box initialized from the list of points.
        :rtype: BoundingBox
        """
        # if len(points) == 0:
        #     return (0, 0, 0, 0, 0, 0)
        xmin = min(pt.x for pt in points)
        xmax = max(pt.x for pt in points)
        ymin = min(pt.y for pt in points)
        ymax = max(pt.y for pt in points)
        zmin = min(pt.z for pt in points)
        zmax = max(pt.z for pt in points)
        return cls(xmin, xmax, ymin, ymax, zmin, zmax)

    def to_frame(self) -> volmdlr.Frame3D:
        """
        Converts the bounding box to a 3D frame.

        :return: A 3D frame with origin at the center and axes aligned with the x, y, and z dimensions of
            the bounding box.
        :rtype: volmdlr.Frame3D
        """
        x = volmdlr.Vector3D((self.xmax - self.xmin), 0, 0)
        y = volmdlr.Vector3D(0, (self.ymax - self.ymin), 0)
        z = volmdlr.Vector3D(0, 0, (self.zmax - self.zmin))
        return volmdlr.Frame3D(self.center, x, y, z)

    def volume(self) -> float:
        """
        Calculates the volume of a bounding box.

        :return: The volume of the bounding box.
        :rtype: float
        """
        return (self.xmax - self.xmin) * (self.ymax - self.ymin) * (
                self.zmax - self.zmin)

    def bbox_intersection(self, bbox2: "BoundingBox") -> bool:
        """
        Calculates if there is an intersection between two bounding boxes.

        :param bbox2: The second bounding box to compare with the current bounding box (self).
        :type bbox2: BoundingBox
        :return: A boolean value indicating whether the two bounding boxes intersect (True) or not (False).
        :rtype: bool
        """
        if self.xmin < bbox2.xmax and self.xmax > bbox2.xmin:
            if self.ymin < bbox2.ymax and self.ymax > bbox2.ymin \
                    and self.zmin < bbox2.zmax and self.zmax > bbox2.zmin:
                return True
        if self.xmin == bbox2.xmax and self.xmax == bbox2.xmin:
            if self.ymin < bbox2.ymax and self.ymax > bbox2.ymin \
                    and self.zmin < bbox2.zmax and self.zmax > bbox2.zmin:
                return True
        return False

    def is_inside_bbox(self, bbox2: "BoundingBox") -> bool:
        """
        Checks if a bounding box is contained inside another bounding box.

        :param bbox2: The bounding box to check against.
        :type bbox2: BoundingBox
        :return: True if the bounding box is contained inside bbox2, False otherwise.
        :rtype: bool
        """
        return (self.xmin >= bbox2.xmin - 1e-6) and (self.xmax <= bbox2.xmax + 1e-6) \
            and (self.ymin >= bbox2.ymin - 1e-6) and (self.ymax <= bbox2.ymax + 1e-6) \
            and (self.zmin >= bbox2.zmin - 1e-6) and (self.zmax <= bbox2.zmax + 1e-6)

    def intersection_volume(self, bbox2: "BoundingBox") -> float:
        """
        Calculate the volume of the intersection of two bounding boxes.

        :param bbox2: The second bounding box to intersect with the first one.
        :type bbox2: BoundingBox
        :return: The volume of the intersection of two bounding boxes.
        :rtype: float
        """
        if not self.bbox_intersection(bbox2):
            return 0
        if self.is_inside_bbox(bbox2) or bbox2.is_inside_bbox(self):
            return min(self.volume(), bbox2.volume())

        lx = min(self.xmax, bbox2.xmax) - max(self.xmin, bbox2.xmin)
        ly = min(self.ymax, bbox2.ymax) - max(self.ymin, bbox2.ymin)
        lz = min(self.zmax, bbox2.zmax) - max(self.zmin, bbox2.zmin)

        return lx * ly * lz

    def distance_to_bbox(self, bbox2: "BoundingBox") -> float:
        """
        Calculates the distance between the bounding box and another bounding box.

        If the bounding boxes intersect, the distance is 0.
        Otherwise, the distance is the minimum Euclidean distance between their closest faces.

        :param bbox2: Another bounding box to compare with.
        :type bbox2: BoundingBox
        :return: The distance between the bounding boxes.
        :rtype: float
        """

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

    def point_belongs(self, point: volmdlr.Point3D) -> bool:
        """
        Determines if a point belongs to the bounding box.

        :param point: The point to check for inclusion.
        :type point: volmdlr.Point3D
        :return: True if the point belongs to the bounding box, False otherwise.
        :rtype: bool
        """
        return (
                self.xmin < point[0] < self.xmax
                and self.ymin < point[1] < self.ymax
                and self.zmin < point[2] < self.zmax
        )

    def distance_to_point(self, point: volmdlr.Point3D) -> float:
        """
        Calculates the minimum Euclidean distance between the bounding box and a point.

        :param point: The point to compare with.
        :type point: volmdlr.Point3D
        :return: The minimum distance between the point and the bounding box.
        :rtype: float
        """
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
            # TODO: if 2 volume models has exact same primitives but in a different order, they are different?
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

        It does not make any boolean operation in case of overlapping.

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
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

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
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for primitives in self.primitives:
            primitives.translation_inplace(offset)
        self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new VolumeModel.

        side = 'old' or 'new'
        """
        new_primitives = [primitive.frame_mapping(frame, side)
                          for primitive in self.primitives]
        return VolumeModel(new_primitives, self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        side = 'old' or 'new'.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for primitives in self.primitives:
            primitives.frame_mapping_inplace(frame, side)
        self.bounding_box = self._bounding_box()

    def copy(self, deep=True, memo=None):
        """
        Specific copy.
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

    def babylon_data(self):
        """
        Get babylonjs data.

        :return: Dictionary with babylon data.
        """

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
    def babylonjs_script(cls, babylon_data, use_cdn=True, debug=False):
        """
        Run babylonjs script.

        """
        if use_cdn:
            script = volmdlr.templates.BABYLON_UNPACKER_CDN_HEADER  # .substitute(name=page_name)
        else:
            script = volmdlr.templates.BABYLON_UNPACKER_EMBEDDED_HEADER  # .substitute(name=page_name)

        script += volmdlr.templates.BABYLON_UNPACKER_BODY_TEMPLATE.substitute(
            babylon_data=babylon_data)
        return script

    def babylonjs(self, page_name=None, use_cdn=True, debug=False):
        """
        Creates a html file using babylonjs to show a 3d model in the browser.

        """
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
            with open(page_name, 'w', encoding='utf-8') as file:
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

            with open(filename, 'w', encoding='utf-8') as file:
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
        with open(filepath, 'w', encoding='utf-8') as file:
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
            step_content += (f"#{product_definition_context_id} = "
                             + "PRODUCT_DEFINITION_CONTEXT('part definition',#2,'design');\n")

            product_context_id = product_definition_context_id + 1
            step_content += f"#{product_context_id} = PRODUCT_CONTEXT('',#2,'mechanical');\n"
            product_id = product_context_id + 1
            step_content += f"#{product_id} = PRODUCT('{primitive.name}'," \
                            f"'{primitive.name}','',(#{product_context_id}));\n"
            product_definition_formation_id = product_id + 1
            step_content += f"#{product_definition_formation_id} = " \
                            f"PRODUCT_DEFINITION_FORMATION('','',#{product_id});\n"
            product_definition_id = product_definition_formation_id + 1
            step_content += f"#{product_definition_id} = PRODUCT_DEFINITION('design'," \
                            f"'',#{product_definition_formation_id},#{product_definition_context_id});\n"
            product_definition_shape_id = product_definition_id + 1
            step_content += f"#{product_definition_shape_id} = PRODUCT_DEFINITION_SHAPE(''," \
                            f"'',#{product_definition_id});\n"
            shape_definition_repr_id = product_definition_shape_id + 1
            step_content += f"#{shape_definition_repr_id} = SHAPE_DEFINITION_REPRESENTATION(" \
                            f"#{product_definition_shape_id},#{primitive_id});\n"
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

    def get_geo_lines(self):
        """
        Gets the lines that define a VolumeModel geometry in a .geo file.

        :return: A list of lines that describe the geometry
        :rtype: List[str]

        """

        update_data = {'point_account': 0,
                       'line_account': 0,
                       'line_loop_account': 0,
                       'surface_account': 0,
                       'surface_loop_account': 0}

        lines = []
        volume = 0
        for primitive in self.primitives:
            if isinstance(primitive, volmdlr.faces.ClosedShell3D):
                volume += 1
                lines_primitives, update_data = primitive.get_geo_lines(update_data)
                lines.extend(lines_primitives)
                surface_loop = ((lines[-1].split('('))[1].split(')')[0])
                lines.append('Volume(' + str(volume) + ') = {' + surface_loop + '};')
            elif isinstance(primitive, volmdlr.faces.OpenShell3D):
                lines_primitives, update_data = primitive.get_geo_lines(update_data)
                lines.extend(lines_primitives)

        return lines

    def get_mesh_lines(self,
                       factor: float, **kwargs):
        # curvature_mesh_size: int = 0,
        # min_points: int = None,
        # initial_mesh_size: float = 5):
        """
        Gets the lines that define mesh parameters for a VolumeModel, to be added to a .geo file.

        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A list of lines that describe mesh parameters
        :rtype: List[str]
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        # try:
        #     curvature_mesh_size = kwargs['curvature_mesh_size']
        # except KeyError:
        #     curvature_mesh_size = 0
        # try:
        #     min_points = kwargs['min_points']
        # except KeyError:
        #     min_points = None
        # try:
        #     initial_mesh_size = kwargs['initial_mesh_size']
        # except KeyError:
        #     initial_mesh_size = 5

        # meshsizes_max = []
        field_num = 1
        field_nums = []
        lines = []

        lines.append('Mesh.CharacteristicLengthMin = 0;')
        lines.append('Mesh.CharacteristicLengthMax = 1e+22;')
        lines.append('Geometry.Tolerance = 1e-20;')
        lines.append('Mesh.AngleToleranceFacetOverlap = 0.01;')
        lines.append('General.Verbosity = 0;')

        for i, primitive in enumerate(self.primitives):
            if isinstance(primitive, volmdlr.faces.ClosedShell3D):
                bbx = primitive.bounding_box
                dim1, dim2, dim3 = (bbx.xmax - bbx.xmin), (bbx.ymax - bbx.ymin), (bbx.zmax - bbx.zmin)
                volume = dim1 * dim2 * dim3

                if factor == 0:
                    factor = 1e-3

                size = ((volume ** (1. / 3.)) / kwargs['initial_mesh_size']) * factor

                # meshsizes_max.append(size)

                if kwargs['min_points']:
                    lines.extend(primitive.get_mesh_lines_with_transfinite_curves(min_points=kwargs['min_points'],
                                                                                  size=size))

                    # primitives, primitives_length = [], []
                    # for face in primitive.faces:
                    #     for _, contour in enumerate(list(chain(*[[face.outer_contour3d], face.inner_contours3d]))):
                    #         if isinstance(contour, volmdlr.wires.Circle2D):
                    #             primitives.append(contour)
                    #             primitives.append(contour)
                    #             primitives_length.append(contour.length() / 2)
                    #             primitives_length.append(contour.length() / 2)
                    #         else:
                    #             for _, primitive_c in enumerate(contour.primitives):
                    #                 if ((primitive_c not in primitives)
                    #                         and (primitive_c.reverse() not in primitives)):
                    #                     primitives.append(primitive_c)
                    #                     primitives_length.append(primitive_c.length())

                    # for i, length in enumerate(primitives_length):
                    #     if length < kwargs['min_points'] * size:
                    #         lines.append('Transfinite Curve {' + str(i) + '} = ' +
                    #                      str(kwargs['min_points']) + ' Using Progression 1;')

                lines.append('Field[' + str(field_num) + '] = MathEval;')
                lines.append('Field[' + str(field_num) + '].F = "' + str(size) + '";')

                lines.append('Field[' + str(field_num + 1) + '] = Restrict;')
                lines.append('Field[' + str(field_num + 1) + '].InField = ' + str(field_num) + ';')
                lines.append('Field[' + str(field_num + 1) + '].VolumesList = {' + str(i + 1) + '};')
                field_nums.append(field_num + 1)
                field_num += 2

            elif isinstance(primitive, volmdlr.faces.OpenShell3D):
                continue

        # meshsize_max = max(meshsizes_max)
        # meshsize_min = meshsize_max/100

        # lines.append('Mesh.CharacteristicLengthMin = ' + str(meshsize_min) + ';')
        # lines.append('Mesh.CharacteristicLengthMax = ' + str(meshsize_max) + ';')

        lines.append('Field[' + str(field_num) + '] = MinAniso;')
        lines.append('Field[' + str(field_num) + '].FieldsList = {' + str(field_nums)[1:-1] + '};')
        lines.append('Background Field = ' + str(field_num) + ';')

        lines.append('Mesh.MeshSizeFromCurvature = ' + str(kwargs['curvature_mesh_size']) + ';')

        lines.append('Coherence;')

        return lines

    def to_geo_stream(self, stream: dcf.StringFile,
                      factor: float, **kwargs):
        """
        Gets the .geo file for the VolumeModel.

        :param file_name: The geo. file name
        :type file_name: str
        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A txt file
        :rtype: .txt
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        lines = self.get_geo_lines()
        lines.extend(self.get_mesh_lines(factor,
                                         curvature_mesh_size=kwargs['curvature_mesh_size'],
                                         min_points=kwargs['min_points'],
                                         initial_mesh_size=kwargs['initial_mesh_size']))

        content = ''
        for line in lines:
            content += line + '\n'

        stream.write(content)

    def to_geo(self, file_name: str = '',
               factor: float = 0.5, **kwargs):
        # curvature_mesh_size: int = 0,
        # min_points: int = None,
        # initial_mesh_size: float = 5):
        """
        Gets the .geo file for the VolumeModel.

        :param file_name: The geo. file name
        :type file_name: str
        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A txt file
        :rtype: .txt
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        if not (file_name.endswith('.geo') or file_name.endswith('.geo')):
            file_name += '.geo'

        with open(file_name, mode='w', encoding='utf-8') as file:

            self.to_geo_stream(file, factor,
                               curvature_mesh_size=kwargs['curvature_mesh_size'],
                               min_points=kwargs['min_points'],
                               initial_mesh_size=kwargs['initial_mesh_size'])

        # for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
        #     if element[0] not in kwargs:
        #         kwargs[element[0]] = element[1]

        # # try:
        # #     curvature_mesh_size = kwargs['curvature_mesh_size']
        # # except KeyError:
        # #     curvature_mesh_size = 0
        # # try:
        # #     min_points = kwargs['min_points']
        # # except KeyError:
        # #     min_points = None
        # # try:
        # #     initial_mesh_size = kwargs['initial_mesh_size']
        # # except KeyError:
        # #     initial_mesh_size = 5

        # lines = self.get_geo_lines()
        # lines.extend(self.get_mesh_lines(factor,
        #                                   curvature_mesh_size=kwargs['curvature_mesh_size'],
        #                                   min_points=kwargs['min_points'],
        #                                   initial_mesh_size=kwargs['initial_mesh_size']))
        # with open(file_name + '.geo', 'w', encoding="utf-8") as file:
        #     for line in lines:
        #         file.write(line)
        #         file.write('\n')
        # file.close()

    def to_geo_with_stl(self, file_name: str,
                        factor: float, **kwargs):
        # curvature_mesh_size: int = 0,
        # min_points: int = None,
        # initial_mesh_size: float = 5):
        """
        Gets the .geo file for the VolumeModel, with saving each closed shell in a stl file.

        :param file_name: The geo. file name
        :type file_name: str
        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A txt file
        :rtype: .txt
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        # try:
        #     curvature_mesh_size = kwargs['curvature_mesh_size']
        # except KeyError:
        #     curvature_mesh_size = 0
        # try:
        #     min_points = kwargs['min_points']
        # except KeyError:
        #     min_points = None
        # try:
        #     initial_mesh_size = kwargs['initial_mesh_size']
        # except KeyError:
        #     initial_mesh_size = 5

        lines = self.get_geo_lines()
        lines.extend(self.get_mesh_lines(factor,
                                         curvature_mesh_size=kwargs['curvature_mesh_size'],
                                         min_points=kwargs['min_points'],
                                         initial_mesh_size=kwargs['initial_mesh_size']))

        contours, faces_account = [], 0
        surfaces = []
        for i, primitive in enumerate(self.primitives):
            if i == 0:
                surfaces.append(list(range(1, 1 + len(primitive.faces))))
                face_contours = [face.outer_contour3d for face in primitive.faces]
                contours.append(face_contours)
                lines.append('Mesh 2;')
                lines.append('Physical Surface(' + str(i + 1) + ') = {' + str(surfaces[i])[1:-1] + '};')
                lines.append('Save "' + file_name + '.stl" ;')
                faces_account += len(primitive.faces) + 1
            else:
                surfaces.append(list(range(faces_account, faces_account + len(primitive.faces))))
                face_contours = [face.outer_contour3d for face in primitive.faces]
                surfaces = self.update_surfaces_list(face_contours, surfaces, contours, i)
                # for k, face_c in enumerate(face_contours):
                #     for l, contour_l in enumerate(contours):
                #         for c, contour in enumerate(contour_l):
                #             if face_c.is_superposing(contour):
                #                 surfaces[i][k] = surfaces[l][c]
                #                 continue
                lines.append('Mesh 2;')
                lines.append('Physical Surface(' + str(i + 1) + ') = {' + str(surfaces[i])[1:-1] + '};')
                lines.append('Save "' + file_name + '.stl" ;')
                faces_account += len(primitive.faces) + 1
                contours.append(face_contours)

        return lines

    @staticmethod
    def update_surfaces_list(face_contours, surfaces, contours, i):
        for k_f, face_c in enumerate(face_contours):
            for l_c, contour_l in enumerate(contours):
                for c_c, contour in enumerate(contour_l):
                    if face_c.is_superposing(contour):
                        surfaces[i][k_f] = surfaces[l_c][c_c]
                        continue
        return surfaces

    def to_msh(self, mesh_dimension: int,
               factor: float, file_name: str = '', **kwargs):
        # curvature_mesh_size: int = 0,
        # min_points: int = None,
        # initial_mesh_size: float = 5):
        """
        Gets .msh file for the VolumeModel generated by gmsh.

        :param file_name: The msh. file name
        :type file_name: str
        :param mesh_dimension: The mesh dimension (1: 1D-Edge, 2: 2D-Triangle, 3D-Tetrahedra)
        :type mesh_dimension: int
        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A txt file
        :rtype: .txt
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        # try:
        #     curvature_mesh_size = kwargs['curvature_mesh_size']
        # except KeyError:
        #     curvature_mesh_size = 0
        # try:
        #     min_points = kwargs['min_points']
        # except KeyError:
        #     min_points = None
        # try:
        #     initial_mesh_size = kwargs['initial_mesh_size']
        # except KeyError:
        #     initial_mesh_size = 5

        if file_name == '':
            with tempfile.NamedTemporaryFile(delete=False) as file:
                file_name = file.name

        self.to_geo(file_name=file_name,
                    factor=factor,
                    curvature_mesh_size=kwargs['curvature_mesh_size'],
                    min_points=kwargs['min_points'],
                    initial_mesh_size=kwargs['initial_mesh_size'])

        self.generate_msh_file(file_name, mesh_dimension)

        # gmsh.initialize()
        # gmsh.open(file_name + ".geo")

        # gmsh.model.geo.synchronize()
        # gmsh.model.mesh.generate(mesh_dimension)

        # gmsh.write(file_name + ".msh")

        # gmsh.finalize()

    # @staticmethod
    # def generate_msh_file(file_name, mesh_dimension):
    #     """
    #     Generates a mesh written in a .msh file using GMSH library.
    #
    #     :param file_name: DESCRIPTION
    #     :type file_name: TYPE
    #     :param mesh_dimension: DESCRIPTION
    #     :type mesh_dimension: TYPE
    #     :return: DESCRIPTION
    #     :rtype: TYPE
    #
    #     """
    #
    #     gmsh.initialize()
    #     gmsh.open(file_name + ".geo")
    #
    #     gmsh.model.geo.synchronize()
    #     gmsh.model.mesh.generate(mesh_dimension)
    #
    #     gmsh.write(file_name + ".msh")
    #
    #     gmsh.finalize()
    #
    # def to_msh_stream(self, mesh_dimension: int,
    #                   factor: float, stream: dcf.StringFile,
    #                   file_name: str = '', **kwargs):
    #     """
    #     Gets .msh file for the VolumeModel generated by gmsh.
    #
    #     :param file_name: The msh. file name
    #     :type file_name: str
    #     :param mesh_dimension: The mesh dimension (1: 1D-Edge, 2: 2D-Triangle, 3D-Tetrahedra)
    #     :type mesh_dimension: int
    #     :param factor: A float, between 0 and 1, that describes the mesh quality
    #     (1 for coarse mesh - 0 for fine mesh)
    #     :type factor: float
    #     :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
    #     (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
    #     :type curvature_mesh_size: int, optional
    #     :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
    #     on that edge), defaults to None
    #     :type min_points: int, optional
    #     :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
    #     :type initial_mesh_size: float, optional
    #
    #     :return: A txt file
    #     :rtype: .txt
    #     """
    #
    #     for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
    #         if element[0] not in kwargs:
    #             kwargs[element[0]] = element[1]
    #
    #     if file_name == '':
    #         with tempfile.NamedTemporaryFile(delete=False) as file:
    #             file_name = file.name
    #
    #     self.to_geo(file_name=file_name,
    #                 factor=factor,
    #                 curvature_mesh_size=kwargs['curvature_mesh_size'],
    #                 min_points=kwargs['min_points'],
    #                 initial_mesh_size=kwargs['initial_mesh_size'])
    #
    #     gmsh.initialize()
    #     gmsh.open(file_name + ".geo")
    #
    #     gmsh.model.geo.synchronize()
    #     gmsh.model.mesh.generate(mesh_dimension)
    #
    #     lines = []
    #     lines.append('$MeshFormat')
    #     lines.append('4.1 0 8')
    #     lines.append('$EndMeshFormat')
    #
    #     lines.extend(self.get_nodes_lines(gmsh))
    #     lines.extend(self.get_elements_lines(gmsh))
    #
    #     content = ''
    #     for line in lines:
    #         content += line + '\n'
    #
    #     stream.write(content)
    #
    #     # gmsh.finalize()

    def to_msh_file(self, mesh_dimension: int,
                    factor: float, file_name: str = '', **kwargs):
        """ Convert and write model to a .msh file. """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        if file_name == '':
            with tempfile.NamedTemporaryFile(delete=False) as file:
                file_name = file.name

        with open(file_name, mode='w', encoding='utf-8') as file:
            self.to_msh_stream(mesh_dimension,
                               factor, file,
                               curvature_mesh_size=kwargs['curvature_mesh_size'],
                               min_points=kwargs['min_points'],
                               initial_mesh_size=kwargs['initial_mesh_size'])

    @staticmethod
    def get_nodes_lines(gmsh_model):
        lines_nodes = []
        lines_nodes.append('$Nodes')

        tag = None
        entities = gmsh_model.model.getEntities()
        for dim, tag in entities:
            node_tags, node_coords, _ = gmsh_model.model.mesh.getNodes(dim, tag)

            lines_nodes.append(str(dim) + ' ' + str(tag) + ' ' + '0 ' + str(len(node_tags)))
            for tag in node_tags:
                lines_nodes.append(str(tag))
            for n in range(0, len(node_coords), 3):
                lines_nodes.append(str(node_coords[n:n + 3])[1:-1])

        lines_nodes.insert(1, str(len(entities)) + ' ' + str(tag) + ' 1 ' + str(tag))
        lines_nodes.append('$EndNodes')

        return lines_nodes

    @staticmethod
    def get_elements_lines(gmsh_model):
        lines_elements = []
        lines_elements.append('$Elements')

        entities = gmsh_model.model.getEntities()
        for dim, tag in entities:
            elem_types, elem_tags, elem_node_tags = gmsh_model.model.mesh.getElements(dim, tag)

            lines_elements.append(str(dim) + ' ' + str(tag) + ' ' + str(elem_types[0]) + ' ' + str(len(elem_tags[0])))
            range_list = int(len(elem_node_tags[0]) / len(elem_tags[0]))
            for n in range(0, len(elem_node_tags[0]), range_list):
                lines_elements.append(str(elem_tags[0][int(n / range_list)]) + ' ' +
                                      str(elem_node_tags[0][n:n + range_list])[1:-1])

        tag = str(elem_tags[0][int(n / range_list)])
        lines_elements.insert(1, str(len(entities)) + ' ' + tag + ' 1 ' + tag)
        lines_elements.append('$EndElements')

        return lines_elements


class MovingVolumeModel(VolumeModel):
    """
    A volume model with possibility to declare time steps at which the primitives are positioned with frames.

    """

    def __init__(self, primitives: List[Primitive3D], step_frames: List[List[volmdlr.Frame3D]], name: str = ''):
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

    def step_volume_model(self, istep: int):
        primitives = []
        for primitive, frame in zip(self.primitives, self.step_frames[istep]):
            primitives.append(
                primitive.frame_mapping(frame, side='old'))
        return VolumeModel(primitives)

    def babylon_data(self):
        """
        Get babylonjs data.

        :return: Dictionary with babylonjs data.
        """
        meshes = []
        primitives_to_meshes = []
        for i_prim, primitive in enumerate(self.primitives):
            if hasattr(primitive, 'babylon_meshes'):
                meshes.extend(primitive.babylon_meshes())
                primitives_to_meshes.append(i_prim)

        bbox = self._bounding_box()
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
                        'center': list(bbox.center),
                        'steps': steps}
        return babylon_data
