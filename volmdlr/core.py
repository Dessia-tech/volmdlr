#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Base classes.
"""
import warnings
from dataclasses import dataclass
from functools import lru_cache
from typing import List, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray

import dessia_common.core as dc

import volmdlr
import volmdlr.templates
from volmdlr.core_compiled import bbox_is_intersecting
from volmdlr.discrete_representation_compiled import triangle_intersects_voxel
from volmdlr.geometry import get_transfer_matrix_from_basis

np.seterr(divide='raise')

DEFAULT_COLOR = (0.8, 0.8, 0.8)


def element_in_list(element, list_elements, tol: float = 1e-6):
    """
    Verifies if a volmdlr element is inside a list  of elements, considering a certain tolerance.

    :param element: Element to be verified inside list.
    :param list_elements: List of elements to be used.
    :param tol: Tolerance to consider if two points are the same.
    :return: True if there is an element inside the list close to the element to given tolerance.
    """
    for element_i in list_elements:
        if element.is_close(element_i, tol):
            return True
    return False


def edge_in_list(edge, list_edges, tol: float = 1e-6):
    """
    Verifies if an edge is inside a list  of edges, considering a certain tolerance.

    :param edge: Edge to be verified inside list.
    :param list_edges: List of edges to be used.
    :param tol: Tolerance to consider if two points are the same.
    :return: True if there is an edge inside the list close to the edge to given tolerance.
    """

    return element_in_list(edge, list_edges, tol)


def get_element_index_in_list(element, list_elements, tol: float = 1e-6):
    """
    Gets the index an element inside a list of elements, considering a certain tolerance.

    :param element: Element to be verified inside list.
    :param list_elements: List of elements to be used.
    :param tol: Tolerance to consider if two elements are the same.
    :return: The element index.
    """
    for i, element_i in enumerate(list_elements):
        if element_i.is_close(element, tol):
            return i
    return None


def get_point_index_in_list(point, list_points, tol: float = 1e-6):
    """
    Gets the index a point inside a list of points, considering a certain tolerance.

    :param point: Point to be verified inside list.
    :param list_points: List of points to be used.
    :param tol: Tolerance to consider if two points are the same.
    :return: The point index.
    """

    return get_element_index_in_list(point, list_points, tol)


def get_edge_index_in_list(edge, list_edges, tol: float = 1e-6):
    """
    Gets the index a edge inside a list of edges, considering a certain tolerance.

    :param edge: Edge to be verified inside list.
    :param list_edges: List of edges to be used.
    :param tol: Tolerance to consider if two edges are the same.
    :return: The edge index.
    """

    return get_element_index_in_list(edge, list_edges, tol)


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
    points = []
    points_set = set()
    for point in list_point:
        if point not in points_set:
            points_set.add(point)
            points.append(point)
        else:
            continue
    return points


def map_primitive_with_initial_and_final_frames(primitive, initial_frame, final_frame):
    """
    Frame maps a primitive in an assembly to its good position.

    :param primitive: primitive to map
    :type primitive: Primitive3D
    :param initial_frame: Initial frame
    :type initial_frame: volmdlr.Frame3D
    :param final_frame: The frame resulted after applying a transformation to the initial frame
    :type final_frame: volmdlr.Frame3D
    :return: A new positioned primitive
    :rtype: Primitive3D

    """
    if initial_frame == final_frame:
        return primitive
    if initial_frame == primitive:
        return final_frame
    transfer_matrix = get_transfer_matrix_from_basis(initial_frame.basis(), final_frame.basis())
    u_vector = volmdlr.Vector3D(*transfer_matrix[0])
    v_vector = volmdlr.Vector3D(*transfer_matrix[1])
    w_vector = volmdlr.Vector3D(*transfer_matrix[2])
    new_frame = volmdlr.Frame3D(final_frame.origin, u_vector, v_vector, w_vector)
    if new_frame == volmdlr.OXYZ:
        return primitive
    new_primitive = primitive.frame_mapping(new_frame, 'old')
    return new_primitive


def helper_babylon_data(babylon_data, display_points):
    """Helper function to babylon_data."""
    # Compute max length in each direction
    all_positions = []
    all_points = []
    for mesh in babylon_data["meshes"]:
        all_positions += _extract_positions(mesh)

    for line in babylon_data["lines"]:
        points = line["points"]
        all_points.extend(points)
    if display_points:
        all_points.extend(display_points)

    # Convert to a NumPy array and reshape
    positions_array = np.array([])
    if all_points and all_positions:
        positions_array = np.concatenate((np.array(all_positions).reshape(-1, 3), np.array(all_points)))
    elif all_positions:
        positions_array = np.array(all_positions).reshape(-1, 3)
    elif all_points:
        positions_array = np.array(all_points)
    # Compute min and max for each dimension
    min_vals = positions_array.min(axis=0)
    max_vals = positions_array.max(axis=0)

    # Calculate max length of the bounding box
    max_length = np.max(max_vals - min_vals)

    # Calculate center point of the bounding box
    center = (0.5 * (min_vals + max_vals)).tolist()

    babylon_data['max_length'] = max_length
    babylon_data['center'] = center

    return babylon_data


def _extract_positions(mesh):
    """Helper function to extract positions from babylon_data."""
    all_positions = []

    for primitives_mesh in mesh.get("primitives_meshes", []):
        all_positions += _extract_positions(primitives_mesh)

    all_positions += mesh.get("positions", [])
    return all_positions


def get_babylon_data(shape, merge_meshes=True):
    """
    Get babylonjs data.

    :return: Dictionary with babylon data.
    """

    babylon_data_ = {'meshes': [],
                     'lines': []}
    display_points = []
    for primitive in shape.primitives:
        if hasattr(primitive, 'babylon_meshes'):
            babylon_data_['meshes'].extend(primitive.babylon_meshes(merge_meshes=merge_meshes))
        elif hasattr(primitive, 'babylon_curves'):
            curves = primitive.babylon_curves()
            if curves:
                babylon_data_['lines'].append(curves)
        elif hasattr(primitive, 'babylon_data'):
            data = primitive.babylon_data(merge_meshes=merge_meshes)
            babylon_data_['meshes'].extend(mesh for mesh in data.get("meshes"))
            babylon_data_['lines'].extend(line for line in data.get("lines"))
        elif isinstance(primitive, volmdlr.Point3D):
            display_points.append(primitive)
    return helper_babylon_data(babylon_data_, display_points)


@dataclass
class EdgeStyle:
    """
    Data class for styling edges Matplotlib plots.

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
    equal_aspect: bool = False


class Primitive3D(dc.PhysicalObject):
    """
    Defines a Primitive3D.
    """

    def __init__(self, color: Tuple[float, float, float] = None, alpha: float = 1.0,
                 reference_path: str = volmdlr.PATH_ROOT, name: str = ''):
        self.color = color
        self.alpha = alpha
        self.reference_path = reference_path

        dc.PhysicalObject.__init__(self, name=name)

    def volmdlr_primitives(self):
        """ Return a list of volmdlr primitives to build up volume model."""
        return [self]

    def babylon_param(self):
        """
        Returns babylonjs parameters.

        :return: babylonjs parameters (alpha, name, color)
        :rtype: dict
        """

        babylon_param = {
            'alpha': self.alpha,
            'name': self.name,
            'color': list(self.color) if self.color is not None else [0.8, 0.8, 0.8]
        }

        return babylon_param

    def triangulation(self, *args, **kwargs):
        """
        Get object triangulation.
        """
        raise NotImplementedError(
            f"triangulation method should be implemented on class {self.__class__.__name__}")

    def babylon_meshes(self, *args, **kwargs):
        """
        Returns the babylonjs mesh.
        """
        mesh = self.triangulation()
        if mesh is None:
            return []
        babylon_mesh = mesh.to_babylon()
        babylon_mesh.update(self.babylon_param())
        babylon_mesh["reference_path"] = self.reference_path
        return [babylon_mesh]


class CompositePrimitive3D(Primitive3D):
    """
    A collection of simple primitives3D.
    """
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = []

    def __init__(self, primitives: List[Primitive3D], color: Tuple[float, float, float] = None, alpha: float = 1,
                 reference_path: str = volmdlr.PATH_ROOT, name: str = ""):
        self.primitives = primitives
        Primitive3D.__init__(self, color=color, alpha=alpha, reference_path=reference_path, name=name)
        self._utd_primitives_to_index = False

    def to_dict(self, *args, **kwargs):
        """Avoids storing points in memo that makes serialization slow."""
        return dc.PhysicalObject.to_dict(self, use_pointers=False)

    def plot(self, ax=None, edge_style: EdgeStyle = EdgeStyle()):
        """
        Plot the 3D primitives onto the given Axes3D object.

        :param ax: optional
            The Axes3D object onto which to plot the primitives. If None, a new
            figure and Axes3D object will be created.
        :type ax: Matplotlib plot
        edge_style : optional
            The EdgeStyle to use when plotting the primitives.
        :type edge_style: vme.EdgeStyle
        :return: The Axes3D object onto which the primitives were plotted.
        :rtype: Matplotlib plot
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        for primitive in self.primitives:
            primitive.plot(ax=ax, edge_style=edge_style)
        return ax


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

    def is_intersecting(self, b_rectangle2):
        """
        Returns True if there is an intersection with another specified bounding rectangle or False otherwise.

        :param b_rectangle2: bounding rectangle to verify intersection
        :type b_rectangle2: :class:`BoundingRectangle`
        """
        return self.xmin < b_rectangle2.xmax and self.xmax > b_rectangle2.xmin \
            and self.ymin < b_rectangle2.ymax and self.ymax > b_rectangle2.ymin

    def b_rectangle_intersection(self, b_rectangle2):
        """
        Returns True if there is an intersection with another specified bounding rectangle or False otherwise.

        :param b_rectangle2: bounding rectangle to verify intersection
        :type b_rectangle2: :class:`BoundingRectangle`
        """
        warnings.warn('b_rectangle_intersection is deprecated, please use is_intersecting instead')
        return self.is_intersecting(b_rectangle2)

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

    def point_inside(self, point: volmdlr.Point2D):
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
        if not self.is_intersecting(b_rectangle2):
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
        if self.is_intersecting(b_rectangle2):
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
        if self.point_inside(point):
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
    def from_points(cls, points: List[volmdlr.Point2D], name: str = '') -> "BoundingRectangle":
        """
        Initializes a bounding rectangle from a list of points.

        :param points: The list of points to create the bounding rectangle from.
        :type points: List[volmdlr.Point2D].
        :param name: object's name.
        :return: The bounding rectangle initialized from the list of points.
        :rtype: BoundingRectangle
        """
        points_array = np.array(points)
        # Compute min and max for each dimension
        xmin, ymin = points_array.min(axis=0)
        xmax, ymax = points_array.max(axis=0)
        return cls(xmin, xmax, ymin, ymax, name=name)


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
        self._size = None
        self._octree = None
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
        :type ax: Matplotlib.axes._subplots.Axes3DSubplot, optional
        :param color: The color of the lines used to plot the bounding box.
        :type color: str, optional
        :return: The 3D axes with the plotted bounding box.
        :rtype: Matplotlib.axes._subplots.Axes3DSubplot
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
    def from_bounding_boxes(cls, bounding_boxes: List["BoundingBox"], name: str = "") -> "BoundingBox":
        """
        Create a bounding box that contains multiple bounding boxes.

        :param bounding_boxes: A list of bounding boxes that need to be contained.
        :type bounding_boxes: List[BoundingBox]
        :param name: A name for the bounding box, optional.
        :type name: str

        :return: A new bounding box that contains all the input bounding boxes.
        :rtype: BoundingBox
        """
        # Create a 2D NumPy array where each row corresponds to the coordinates of a bounding box
        # [xmin, xmax, ymin, ymax, zmin, zmax]
        coords = np.array([[bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax] for bb in bounding_boxes])

        # Find the global minimum and maximum for each axis
        mins = np.amin(coords, axis=0)
        maxs = np.amax(coords, axis=0)

        # Assign min and max for each axis
        xmin, xmax, ymin, ymax, zmin, zmax = mins[0], maxs[1], mins[2], maxs[3], mins[4], maxs[5]

        return cls(xmin, xmax, ymin, ymax, zmin, zmax, name=name)

    @classmethod
    def from_points(cls, points: Union[List[volmdlr.Point3D], NDArray], name: str = '') -> "BoundingBox":
        """
        Initializes a bounding box from a list of points.

        :param points: The list of points to create the bounding box from.
        :type points: List[volmdlr.Point3D].
        :param name: object's name.
        :return: The bounding box initialized from the list of points.
        :rtype: BoundingBox
        """
        points_array = np.asarray(points)
        # Compute min and max for each dimension
        xmin, ymin, zmin = points_array.min(axis=0)
        xmax, ymax, zmax = points_array.max(axis=0)

        return cls(xmin, xmax, ymin, ymax, zmin, zmax, name=name)

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

    def get_points_inside_bbox(self, points_x, points_y, points_z):
        """
        Gets points inside the BoudingBox.

        :param points_x: Number of points in x direction.
        :param points_y: Number of points in y direction.
        :param points_z: Number of points in z direction.
        :return: list of points inside bounding box.
        """
        _size = [self.size[0] / points_x, self.size[1] / points_y,
                 self.size[2] / points_z]
        initial_center = self.center.translation(
            -volmdlr.Vector3D(self.size[0] / 2 - _size[0] / 2,
                              self.size[1] / 2 - _size[1] / 2,
                              self.size[2] / 2 - _size[2] / 2))
        points = []
        for z_box in range(points_z):
            for y_box in range(points_y):
                for x_box in range(points_x):
                    translation_vector = volmdlr.Vector3D(x_box * _size[0], y_box * _size[1],
                                                          z_box * _size[2])
                    point = initial_center.translation(translation_vector)
                    points.append(point)
        return points

    @property
    def size(self):
        """Gets the Size of the Bounding Box."""

        if not self._size:
            self._size = [self.xmax - self.xmin, self.ymax - self.ymin, self.zmax - self.zmin]
        return self._size

    def volume(self) -> float:
        """
        Calculates the volume of a bounding box.

        :return: The volume of the bounding box.
        :rtype: float
        """
        return (self.xmax - self.xmin) * (self.ymax - self.ymin) * (self.zmax - self.zmin)

    def scale(self, factor: float) -> "BoundingBox":
        """
        Scales the bounding box by a given factor and returns a new BoundingBox.

        :param factor: The scaling factor.
        :type factor: float

        :return: A new scaled BoundingBox.
        :rtype: BoundingBox
        """
        x_center = (self.xmin + self.xmax) / 2
        y_center = (self.ymin + self.ymax) / 2
        z_center = (self.zmin + self.zmax) / 2
        x_size, y_size, z_size = self.size

        scaled_half_x_size = (x_size * factor) / 2
        scaled_half_y_size = (y_size * factor) / 2
        scaled_half_z_size = (z_size * factor) / 2

        # Calculate new min and max values
        new_xmin = x_center - scaled_half_x_size
        new_xmax = x_center + scaled_half_x_size
        new_ymin = y_center - scaled_half_y_size
        new_ymax = y_center + scaled_half_y_size
        new_zmin = z_center - scaled_half_z_size
        new_zmax = z_center + scaled_half_z_size

        # Return a new BoundingBox object
        return BoundingBox(new_xmin, new_xmax, new_ymin, new_ymax, new_zmin, new_zmax, self.name)

    def is_intersecting(self, bbox2, tol: float = 1e-6):
        """
        Checks if two bounding boxes are intersecting or touching.

        :param self: BoundingBox object representing the first bounding box.
        :param bbox2: BoundingBox object representing the second bounding box.
        :param tol: tolerance to be considered.

        :return: True if the bounding boxes are intersecting or touching, False otherwise.
        """
        return bbox_is_intersecting(self, bbox2, tol)

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
        if not self.is_intersecting(bbox2):
            return 0
        if self.is_inside_bbox(bbox2) or bbox2.is_inside_bbox(self):
            return min(self.volume(), bbox2.volume())

        lx = min(self.xmax, bbox2.xmax) - max(self.xmin, bbox2.xmin)
        ly = min(self.ymax, bbox2.ymax) - max(self.ymin, bbox2.ymin)
        lz = min(self.zmax, bbox2.zmax) - max(self.zmin, bbox2.zmin)

        return lx * ly * lz

    def is_intersecting_triangle(self, triangle: "Triangle3D") -> bool:
        """
        Check if the bounding box and a triangle are intersecting or touching.

        :param triangle: the triangle to check if there is an intersection with.
        :type triangle: Triangle3D

        :return: True if the bounding box and the triangle are intersecting or touching, False otherwise.
        :rtype: bool
        """
        _triangle = tuple((point.x, point.y, point.z) for point in triangle.points)
        _center = (self.center[0], self.center[1], self.center[2])
        _extents = tuple(size / 2 for size in self.size)

        return triangle_intersects_voxel(_triangle, _center, _extents)

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

        if self.is_intersecting(bbox2):
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

    def point_inside(self, point: volmdlr.Point3D, tol=1e-6) -> bool:
        """
        Determines if a point belongs to the bounding box.

        :param point: The point to check for inclusion.
        :type point: volmdlr.Point3D
        :param tol: tolerance.
        :return: True if the point belongs to the bounding box, False otherwise.
        :rtype: bool
        """
        return (
                self.xmin - tol <= point[0] <= self.xmax + tol
                and self.ymin - tol <= point[1] <= self.ymax + tol
                and self.zmin - tol <= point[2] <= self.zmax + tol
        )

    def distance_to_point(self, point: volmdlr.Point3D) -> float:
        """
        Calculates the minimum Euclidean distance between the bounding box and a point.

        :param point: The point to compare with.
        :type point: volmdlr.Point3D
        :return: The minimum distance between the point and the bounding box.
        :rtype: float
        """
        if self.point_inside(point):
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

    def is_close(self, other_bounding_box: "BoundingBox", tol: float = 1e-6) -> bool:
        """
        Check if two bounding boxes are close to each other considering the Euclidean distance of their corners points.

        The tolerance can be modified.

        :param other_bounding_box: the other bounding box.
        :type other_bounding_box: BoundingBox
        :param tol: The tolerance under which the Euclidean distance is considered equal to 0.
        :type tol: float

        :return: True if the bounding boxes are equal at the given tolerance, False otherwise.
        :rtype: bool
        """
        self_corner_min = volmdlr.Point3D(self.xmin, self.ymin, self.zmin)
        self_conrer_max = volmdlr.Point3D(self.xmax, self.ymax, self.zmax)
        other_corner_min = volmdlr.Point3D(other_bounding_box.xmin, other_bounding_box.ymin, other_bounding_box.zmin)
        other_corner_max = volmdlr.Point3D(other_bounding_box.xmax, other_bounding_box.ymax, other_bounding_box.zmax)

        return self_corner_min.is_close(other_corner_min, tol) and self_conrer_max.is_close(other_corner_max, tol)

    def octree(self):
        """Creates a simple octree structure for a bounding box."""
        if not self._octree:
            octants = []
            points_x, points_y, points_z = 2, 2, 2
            _size = [self.size[0] / points_x, self.size[1] / points_y,
                     self.size[2] / points_z]
            octants_center = self.get_points_inside_bbox(points_x, points_y, points_z)
            for octant_center in octants_center:
                mins_maxs = []
                for i, size_component in enumerate(_size):
                    mins_maxs.extend([octant_center[i] - size_component / 2, octant_center[i] + size_component / 2])
                octants.append(self.__class__(mins_maxs[0], mins_maxs[1], mins_maxs[2], mins_maxs[3],
                                              mins_maxs[4], mins_maxs[5]))
            self._octree = octants
        return self._octree


class VolumeModel:
    """
    VolumeModel is deprecated, please, use volmdlr.model.VolumeModel instead.
    """

    def __init__(self, primitives: List[Primitive3D], name: str = ''):
        warnings.warn("volmdlr.core.VolumeModel is deprecated and will be removed in future releases, please "
                      "use volmdlr.model.VolumeModel instead", UserWarning)
        # pylint: disable=cyclic-import, import-outside-toplevel
        from volmdlr import model
        self.instance = model.VolumeModel(primitives=primitives, name=name)

    def __getattr__(self, name):
        return getattr(self.instance, name)
