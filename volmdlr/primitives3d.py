#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common primitives 3D.
"""

import math
import warnings
from random import uniform
from typing import Dict, List, Tuple

import plot_data

import dessia_common.core as dc
import matplotlib.pyplot as plt
import numpy as npy
from scipy.optimize import Bounds, NonlinearConstraint, minimize
from scipy.stats import qmc

import volmdlr
import volmdlr.core
import volmdlr.edges
import volmdlr.faces
import volmdlr.primitives
import volmdlr.wires
from volmdlr import shells, surfaces, curves


npy.seterr(divide='raise')


class RoundedLineSegments3D(volmdlr.primitives.RoundedLineSegments):
    """
    A class representing a series of rounded line segments in 3D.

    This class inherits from the `RoundedLineSegments` class,
    and provides methods to work with rounded line segments in 3D.

    :param points: The list of points defining the line segments.
    :type points: List[volmdlr.Point3D]
    :param radius: The dictionary mapping segment indices to their respective radii.
    :type radius:Dict[int, float]
    :param adapt_radius: Flag indicating whether to adapt the radius based on segment length.
    Defaults to False.
    :type adapt_radius: bool, optional
    :param name: The name of the rounded line segments. Defaults to ''.
    :type name: str, optional
    """
    line_class = volmdlr.edges.LineSegment3D
    arc_class = volmdlr.edges.Arc3D

    def __init__(self, points: List[volmdlr.Point3D], radius: Dict[str, float],
                 adapt_radius: bool = False, name: str = ''):
        volmdlr.primitives.RoundedLineSegments.__init__(self, points, radius, adapt_radius=adapt_radius, name=name)

    def arc_features(self, point_index: int):
        radius = self.radius[point_index]
        point_1, point_i, point_2 = self.get_points(point_index)
        dist1 = (point_1 - point_i).norm()
        dist2 = (point_2 - point_i).norm()
        dist3 = (point_1 - point_2).norm()
        alpha = math.acos(-(dist3**2 - dist1**2 - dist2**2) / (2 * dist1 * dist2)) / 2.
        dist = radius / math.tan(alpha)

        u1 = (point_1 - point_i) / dist1
        u2 = (point_2 - point_i) / dist2

        p3 = point_i + u1 * dist
        p4 = point_i + u2 * dist

        n = u1.cross(u2)
        n /= n.norm()
        v1 = u1.cross(n)
        v2 = u2.cross(n)

        line1 = curves.Line3D(p3, p3 + v1)
        line2 = curves.Line3D(p4, p4 + v2)

        w = u1 + u2  # mean of v1 and v2
        w /= w.norm()

        interior = line1.minimum_distance_points(line2)[0] - w * radius
        return p3, interior, p4, dist, alpha

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        OpenRoundedLineSegments3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated OpenRoundedLineSegments3D
        """
        return self.__class__([point.rotation(center, axis, angle)
                               for point in self.points],
                              self.radius, self.closed, self.name)

    def translation(self, offset: volmdlr.Vector3D):
        """
        OpenRoundedLineSegments3D translation.

        :param offset: translation vector
        :return: A new translated OpenRoundedLineSegments3D
        """
        return self.__class__([point.translation(offset)
                               for point in self.points],
                              self.radius, self.closed, self.name)


class OpenRoundedLineSegments3D(volmdlr.wires.Wire3D, RoundedLineSegments3D):
    """
    Defines an open rounded line segments.

    :param points: Points used to draw the wire.
    :type points: List of Point3D.
    :param radius: Radius used to connect different parts of the wire.
    :type radius: {position1(n): float which is the radius linked the n-1 and.
    n+1 points, position2(n+1):...}
    """
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']

    def __init__(self, points: List[volmdlr.Point3D], radius: Dict[str, float],
                 adapt_radius: bool = False, name: str = ''):
        RoundedLineSegments3D.__init__(self, points, radius, adapt_radius=adapt_radius, name='')
        self.closed = False

        volmdlr.wires.Wire3D.__init__(self, self._primitives(), name)


class ClosedRoundedLineSegments3D(RoundedLineSegments3D, volmdlr.wires.Contour3D):
    """
    Defines a closed rounded line segment in 3D.

    :param points: Points used to draw the wire
    :type points: List of Point3D
    :param radius: Radius used to connect different parts of the wire
    :type radius: {position1(n): float which is the radius linked the n-1 and
    n+1 points, position2(n+1):...}
    """
    _non_serializable_attributes = []
    _non_data_eq_attributes = ['name']
    _non_data_hash_attributes = ['name']

    def __init__(self, points: List[volmdlr.Point3D], radius: float, adapt_radius: bool = False, name: str = ''):
        RoundedLineSegments3D.__init__(self, points, radius, adapt_radius, name)
        self.closed = True
        volmdlr.wires.Contour3D.__init__(self, primitives=self._primitives(), name=name)


class Block(shells.ClosedShell3D):
    """
    Creates a block.

    :param frame: a frame 3D. The origin of the frame is the center of the block,
     the 3 vectors are defining the edges. The frame has not to be orthogonal
    """
    _standalone_in_db = True
    _non_serializable_attributes = ['size', 'bounding_box', 'faces', 'contours', 'plane', 'points', 'polygon2D']
    _non_data_eq_attributes = ['name', 'color', 'alpha', 'size', 'bounding_box',
                               'faces', 'contours', 'plane', 'points', 'polygon2D']
    _non_data_hash_attributes = []

    def __init__(self, frame: volmdlr.Frame3D, *,
                 color: Tuple[float, float, float] = None, alpha: float = 1.,
                 name: str = ''):
        self.frame = frame
        self.size = (self.frame.u.norm(),
                     self.frame.v.norm(),
                     self.frame.w.norm())
        self._octree = None
        self._quadtree = None
        faces = self.shell_faces()
        for face in faces:
            face.alpha = alpha
            face.color = color
        shells.ClosedShell3D.__init__(self, faces, color=color, alpha=alpha, name=name)

    def __eq__(self, other):
        if self.__class__.__name__ != other.__class__.__name__:
            return False
        if self.frame != other.frame:
            return False
        return True

    def __hash__(self):
        return hash((self.__class__.__name__, self.frame))

    def to_dict(self, *args, **kwargs):
        """
        Custom to_dict for performance.

        """
        dict_ = dc.DessiaObject.base_dict(self)
        dict_.update({'color': self.color,
                      'alpha': self.alpha,
                      'frame': self.frame.to_dict()})

        return dict_

    def volume(self):
        """Returns the volume of the block."""
        return self.size[0] * self.size[1] * self.size[2]

    @classmethod
    def from_bounding_box(cls, bounding_box, name: str = ''):
        """
        Transform a bounding box into a block.
        """
        origin = bounding_box.center
        bbox_size = bounding_box.size
        frame = volmdlr.Frame3D(origin, bbox_size[0] * volmdlr.Vector3D(1, 0, 0),
                                bbox_size[1] * volmdlr.Vector3D(0, 1, 0),
                                bbox_size[2] * volmdlr.Vector3D(0, 0, 1))
        return cls(frame=frame, name=name)

    def vertices(self):
        """Computes the vertices of the block."""
        return [self.frame.origin - 0.5 * self.frame.u - 0.5 * self.frame.v - 0.5 * self.frame.w,
                self.frame.origin - 0.5 * self.frame.u + 0.5 * self.frame.v - 0.5 * self.frame.w,
                self.frame.origin + 0.5 * self.frame.u + 0.5 * self.frame.v - 0.5 * self.frame.w,
                self.frame.origin + 0.5 * self.frame.u - 0.5 * self.frame.v - 0.5 * self.frame.w,
                self.frame.origin - 0.5 * self.frame.u - 0.5 * self.frame.v + 0.5 * self.frame.w,
                self.frame.origin - 0.5 * self.frame.u + 0.5 * self.frame.v + 0.5 * self.frame.w,
                self.frame.origin + 0.5 * self.frame.u + 0.5 * self.frame.v + 0.5 * self.frame.w,
                self.frame.origin + 0.5 * self.frame.u - 0.5 * self.frame.v + 0.5 * self.frame.w]

    def edges(self):
        """Computes the edges of the block."""
        point1, point2, point3, point4, point5, point6, point7, point8 = self.vertices()
        return [volmdlr.edges.LineSegment3D(point1.copy(), point2.copy()),
                volmdlr.edges.LineSegment3D(point2.copy(), point3.copy()),
                volmdlr.edges.LineSegment3D(point3.copy(), point4.copy()),
                volmdlr.edges.LineSegment3D(point4.copy(), point1.copy()),
                volmdlr.edges.LineSegment3D(point5.copy(), point6.copy()),
                volmdlr.edges.LineSegment3D(point6.copy(), point7.copy()),
                volmdlr.edges.LineSegment3D(point7.copy(), point8.copy()),
                volmdlr.edges.LineSegment3D(point8.copy(), point5.copy()),
                volmdlr.edges.LineSegment3D(point1.copy(), point5.copy()),
                volmdlr.edges.LineSegment3D(point2.copy(), point6.copy()),
                volmdlr.edges.LineSegment3D(point3.copy(), point7.copy()),
                volmdlr.edges.LineSegment3D(point4.copy(), point8.copy())]

    def face_contours3d(self):
        edge1, edge2, edge3, edge4, edge5, edge6, edge7, edge8, edge9, edge10, edge11, edge12 = self.edges()
        e5_switch = edge5.reverse()
        e6_switch = edge6.reverse()
        e7_switch = edge7.reverse()
        e8_switch = edge8.reverse()
        e9_switch = edge9.reverse()
        e10_switch = edge10.reverse()
        e11_switch = edge11.reverse()
        e12_switch = edge12.reverse()
        return [volmdlr.wires.Contour3D([edge1.copy(), edge2.copy(), edge3.copy(), edge4.copy()]),
                volmdlr.wires.Contour3D([edge5.copy(), edge6.copy(), edge7.copy(), edge8.copy()]),
                # volmdlr.Contour3D([e1.copy(), e9.copy(), e5.copy(), e10.copy()]),
                volmdlr.wires.Contour3D([edge1.copy(), edge10.copy(), e5_switch.copy(), e9_switch.copy()]),
                # volmdlr.Contour3D([e2.copy(), e10.copy(), e6.copy(), e11.copy()]),
                volmdlr.wires.Contour3D([edge2.copy(), edge11.copy(), e6_switch.copy(), e10_switch.copy()]),
                # volmdlr.Contour3D([e3.copy(), e11.copy(), e7.copy(), e12.copy()]),
                volmdlr.wires.Contour3D([edge3.copy(), edge12.copy(), e7_switch.copy(), e11_switch.copy()]),
                # volmdlr.Contour3D([e4.copy(), e12.copy(), e8.copy(), e9.copy()])]
                volmdlr.wires.Contour3D([edge4.copy(), edge9.copy(), e8_switch.copy(), e12_switch.copy()])]

    def shell_faces(self):
        """Computes the faces of the block."""
        hlx = 0.5 * self.frame.u.norm()
        hly = 0.5 * self.frame.v.norm()
        hlz = 0.5 * self.frame.w.norm()
        frame = self.frame.copy()
        frame = frame.normalize()
        xm_frame = volmdlr.Frame3D(frame.origin - 0.5 * self.frame.u,
                                   frame.v, frame.w, frame.u)
        xp_frame = volmdlr.Frame3D(frame.origin + 0.5 * self.frame.u,
                                   frame.v, frame.w, frame.u)
        ym_frame = volmdlr.Frame3D(frame.origin - 0.5 * self.frame.v,
                                   frame.w, frame.u, frame.v)
        yp_frame = volmdlr.Frame3D(frame.origin + 0.5 * self.frame.v,
                                   frame.w, frame.u, frame.v)
        zm_frame = volmdlr.Frame3D(frame.origin - 0.5 * self.frame.w,
                                   frame.u, frame.v, frame.w)
        zp_frame = volmdlr.Frame3D(frame.origin + 0.5 * self.frame.w,
                                   frame.u, frame.v, frame.w)

        xm_face = volmdlr.faces.PlaneFace3D.from_surface_rectangular_cut(
            surfaces.Plane3D(xm_frame), -hly, hly, -hlz, hlz)
        xp_face = volmdlr.faces.PlaneFace3D.from_surface_rectangular_cut(
            surfaces.Plane3D(xp_frame), -hly, hly, -hlz, hlz)
        ym_face = volmdlr.faces.PlaneFace3D.from_surface_rectangular_cut(
            surfaces.Plane3D(ym_frame), -hlz, hlz, -hlx, hlx)
        yp_face = volmdlr.faces.PlaneFace3D.from_surface_rectangular_cut(
            surfaces.Plane3D(yp_frame), -hlz, hlz, -hlx, hlx)
        zm_face = volmdlr.faces.PlaneFace3D.from_surface_rectangular_cut(
            surfaces.Plane3D(zm_frame), -hlx, hlx, -hly, hly)
        zp_face = volmdlr.faces.PlaneFace3D.from_surface_rectangular_cut(
            surfaces.Plane3D(zp_frame), -hlx, hlx, -hly, hly)

        return [xm_face, xp_face, ym_face, yp_face, zm_face, zp_face]

    def faces_center(self):
        """Computes the faces center of the block."""
        return [self.frame.origin - 0.5 * self.frame.u,
                self.frame.origin + 0.5 * self.frame.u,
                self.frame.origin - 0.5 * self.frame.v,
                self.frame.origin + 0.5 * self.frame.v,
                self.frame.origin - 0.5 * self.frame.w,
                self.frame.origin + 0.5 * self.frame.w]

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Block rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Block
        """
        new_frame = self.frame.rotation(center, axis, angle)
        return Block(new_frame, color=self.color,
                     alpha=self.alpha, name=self.name)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Returns a new translated block.

        :param offset: translation vector.
        :return: A new translated Block.
        """
        new_frame = self.frame.translation(offset)
        return Block(new_frame, color=self.color,
                     alpha=self.alpha, name=self.name)

    def cut_by_orthogonal_plane(self, plane_3d: surfaces.Plane3D):
        bouding_box = self.bounding_box
        if plane_3d.frame.w.dot(volmdlr.Vector3D(1, 0, 0)) == 0:
            pass
        elif plane_3d.frame.w.dot(volmdlr.Vector3D(0, 1, 0)) == 0:
            pass
        elif plane_3d.frame.w.dot(volmdlr.Vector3D(0, 0, 1)) == 0:
            pass
        else:
            raise KeyError('plane is not orthogonal either with x, y or z')
        point_min = volmdlr.Point3D(bouding_box.xmin, bouding_box.ymin,
                                    bouding_box.zmin)
        point_max = volmdlr.Point3D(bouding_box.xmax, bouding_box.ymax,
                                    bouding_box.zmax)
        point_min_2d = plane_3d.point3d_to_2d(point_min)
        point_max_2d = plane_3d.point3d_to_2d(point_max)
        points = [point_min_2d, volmdlr.Point2D(point_max_2d.x, point_min_2d.y),
                  point_max_2d, volmdlr.Point2D(point_min_2d.x, point_max_2d.y)]
        contour_2d = surfaces.Surface2D(
            volmdlr.wires.ClosedPolygon2D(points), [])

        return volmdlr.faces.PlaneFace3D(plane_3d, contour_2d)

    def frame_mapping_parametres(self, frame: volmdlr.Frame3D,
                                 side: str):
        basis = frame.basis()
        if side == 'new':
            new_origin = frame.global_to_local_coordinates(self.frame.origin)
            new_u = basis.global_to_local_coordinates(self.frame.u)
            new_v = basis.global_to_local_coordinates(self.frame.v)
            new_w = basis.global_to_local_coordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
        elif side == 'old':
            new_origin = frame.local_to_global_coordinates(self.frame.origin)
            new_u = basis.local_to_global_coordinates(self.frame.u)
            new_v = basis.local_to_global_coordinates(self.frame.v)
            new_w = basis.local_to_global_coordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
        else:
            raise ValueError('side value not valid, please specify'
                             'a correct value: \'old\' or \'new\'')
        return new_frame

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Frame3D.

        :param side: 'old' or 'new'
        """
        new_frame = self.frame_mapping_parametres(frame, side)
        return Block(new_frame, color=self.color,
                     alpha=self.alpha, name=self.name)

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of a Block.

        """
        new_origin = self.frame.origin.copy()
        new_u = self.frame.u.copy()
        new_v = self.frame.v.copy()
        new_w = self.frame.w.copy()
        new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
        return Block(new_frame, color=self.color,
                     alpha=self.alpha, name=self.name)

    def plot_data(self, x3d, y3d, edge_style=plot_data.EdgeStyle):
        """Plot the 2D projections of a block."""
        lines = []
        for edge3d in self.edges():
            lines.append(edge3d.plot_data(x3d, y3d, edge_style))

        return lines

    def plot2d(self, x3d, y3d, ax=None):
        """
        Plot 2d with Matplotlib.

        """
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = None

        for edge3d in self.edges():
            edge3d.plot2d(x3d, y3d, ax=ax)

        return fig, ax

    def octree(self):
        """Subdivide block into eight other blocks."""
        if not self._octree:
            self._octree = self.subdivide_block(2, 2, 2)
        return self._octree

    def quadtree(self):
        """Subdivide block into four other blocks."""
        if not self._quadtree:
            self._quadtree = self.subdivide_block(2, 2, 1)
        return self._quadtree

    def subdivide_block(self, number_blocks_x, number_blocks_y, number_blocks_z):
        """Divide block into sub blocks."""
        filling_boxes_size = [self.size[0] / number_blocks_x, self.size[1] / number_blocks_y,
                              self.size[2] / number_blocks_z]
        initial_frame_center = self.frame.origin.copy(deep=True).translation(
            -volmdlr.Vector3D(self.size[0] / 2 - filling_boxes_size[0] / 2,
                              self.size[1] / 2 - filling_boxes_size[1] / 2,
                              self.size[2] / 2 - filling_boxes_size[2] / 2))
        xyz = [volmdlr.Vector3D(filling_boxes_size[0], 0, 0), volmdlr.Vector3D(0, filling_boxes_size[1], 0),
               volmdlr.Vector3D(0, 0, filling_boxes_size[2])]

        dividing_blocks = []
        for z_box in range(number_blocks_z):
            for y_box in range(number_blocks_y):
                for x_box in range(number_blocks_x):
                    translation_vector = volmdlr.Vector3D(x_box * filling_boxes_size[0], y_box * filling_boxes_size[1],
                                                          z_box * filling_boxes_size[2])
                    origin_point = initial_frame_center.translation(translation_vector)
                    block = Block(frame=volmdlr.Frame3D(origin_point, *xyz))
                    dividing_blocks.append(block)
        return dividing_blocks


class ExtrudedProfile(shells.ClosedShell3D):
    """
    Extrude a profile given by outer and inner contours.

    TODO: In the future change to a frame and a surface2D and an extrusion vector.
    """
    _non_serializable_attributes = ['faces', 'inner_contours3d',
                                    'outer_contour3d']

    def __init__(self, frame: volmdlr.Frame3D,
                 outer_contour2d: volmdlr.wires.Contour2D,
                 inner_contours2d: List[volmdlr.wires.Contour2D],
                 extrusion_length: float,
                 color: Tuple[float, float, float] = None, alpha: float = 1.,
                 name: str = ''):
        self.frame = frame

        self.outer_contour2d = outer_contour2d
        self.outer_contour3d = outer_contour2d.to_3d(self.frame.origin, self.frame.u, self.frame.v)
        self.extrusion_length = extrusion_length
        self.inner_contours2d = inner_contours2d
        self.extrusion_vector = self.frame.w * extrusion_length
        self.inner_contours3d = []
        self.color = color

        bool_areas = []
        for contour in inner_contours2d:
            self.inner_contours3d.append(contour.to_3d(self.frame.origin, self.frame.u, self.frame.v))
            if contour.area() > outer_contour2d.area():
                bool_areas.append(True)
            else:
                bool_areas.append(False)
        if any(bool_areas):
            raise ValueError('At least one inner contour is not contained in outer_contour.')

        faces = self.shell_faces()

        shells.ClosedShell3D.__init__(self, faces, color=color,
                                      alpha=alpha, name=name)

    def to_dict(self, *args, **kwargs):
        """
        Serialize the ExtrudedProfile.

        """
        dict_ = dc.DessiaObject.base_dict(self)
        dict_.update({'color': self.color,
                      'alpha': self.alpha,
                      'frame': self.frame.to_dict(),
                      'outer_contour2d': self.outer_contour2d.to_dict(),
                      'inner_contours2d': [c.to_dict() for c in self.inner_contours2d],
                      'extrusion_length': self.extrusion_length,
                      })

        return dict_

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of Extruded Profile.

        """
        return self.__class__(
            frame=self.frame.copy(),
            outer_contour2d=self.outer_contour2d.copy(),
            inner_contours2d=[c.copy() for c in self.inner_contours2d],
            extrusion_length=self.extrusion_length,
            color=self.color,
            alpha=self.alpha,
            name=self.name)

    def shell_faces(self):
        """
        Computes the shell faces from init data.

        """
        lower_plane = surfaces.Plane3D.from_plane_vectors(
            self.frame.origin, self.frame.u, self.frame.v)
        lower_face = volmdlr.faces.PlaneFace3D(
            lower_plane, surfaces.Surface2D(self.outer_contour2d,
                                            self.inner_contours2d))

        upper_face = lower_face.translation(self.extrusion_vector)
        lateral_faces = []
        for primitive in self.outer_contour3d.primitives:
            lateral_faces.extend(primitive.extrusion(self.extrusion_vector))

        for inner_contour in self.inner_contours3d:
            for primitive in inner_contour.primitives:
                lateral_faces.extend(primitive.extrusion(self.extrusion_vector))

        return [lower_face, upper_face] + lateral_faces

    def area(self):
        """Returns the area of the extruded 2D surface."""
        areas = self.outer_contour2d.area()
        areas -= sum(contour.area() for contour in self.inner_contours2d)
        return areas

    def volume(self):
        z = self.frame.w
        return self.area() * self.extrusion_vector.dot(z)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new ExtrudeProfile.

        :param side: = 'old' or 'new'.
        """
        return ExtrudedProfile(
            self.frame.frame_mapping(frame, side),
            self.outer_contour2d, self.inner_contours2d,
            self.extrusion_length)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Extruded Profile rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated ExtrudedProfile.
        """
        return self.__class__(
            frame=self.frame.rotation(center, axis, angle),
            outer_contour2d=self.outer_contour2d,
            inner_contours2d=self.inner_contours2d,
            extrusion_length=self.extrusion_length,
            color=self.color, alpha=self.alpha)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Extruded Profile translation.

        :param offset: translation vector
        :return: A new translated ExtrudedProfile
        """
        return self.__class__(
            frame=self.frame.translation(offset),
            outer_contour2d=self.outer_contour2d,
            inner_contours2d=self.inner_contours2d,
            extrusion_length=self.extrusion_length,
            color=self.color, alpha=self.alpha)


class RevolvedProfile(shells.ClosedShell3D):
    """
    Revolve a 2D profile along an axis around a certain angle.

    """
    _non_serializable_attributes = ['faces', 'contour3D']

    def __init__(self, frame: volmdlr.Frame3D,
                 contour2d: volmdlr.wires.Contour2D,
                 axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float = 2 * math.pi, *,
                 color: Tuple[float, float, float] = None, alpha: float = 1,
                 name: str = ''):
        if frame.w.cross(axis).is_close(volmdlr.Vector3D(0.0, 0.0, 0.0)):
            raise ValueError(f"The normal vector of the Revolution's contour frame should not be parallel \n"
                             f"to revolution axis. frame.w: {frame.w}; revolution_axis: {axis}")
        self.contour2d = contour2d
        self.axis_point = axis_point
        self.axis = axis
        self.angle = angle
        self.frame = frame
        self.contour3d = self.contour2d.to_3d(frame.origin, frame.u, frame.v)

        faces = self.shell_faces()
        shells.ClosedShell3D.__init__(self, faces, color=color,
                                      alpha=alpha, name=name)

    def __eq__(self, other):
        if not self.__class__.__name__ == other.__class__.__name__:
            return False
        for self_param, other_param in zip([self.frame,
                                            self.contour2d, self.axis_point, self.axis, self.angle],
                                           [other.frame,
                                            other.contour2d, other.axis_point, other.axis, other.angle]
                                           ):
            if not self_param == other_param:
                return False
        return True

    def to_dict(self, *args, **kwargs):
        """
        Custom to dict for performance.
        """
        dict_ = dc.DessiaObject.base_dict(self)
        dict_.update({'color': self.color,
                      'alpha': self.alpha,
                      'frame': self.frame,
                      'contour2d': self.contour2d.to_dict(),
                      'axis_point': self.axis_point.to_dict(),
                      'angle': self.angle,
                      'axis': self.axis.to_dict()
                      })

        return dict_

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of Revolved-profile.

        """
        return self.__class__(frame=self.frame.copy(),
                              contour2d=self.contour2d.copy(deep=deep, memo=memo),
                              axis=self.axis.copy(), angle=self.angle,
                              axis_point=self.axis_point.copy(),
                              color=self.color, alpha=self.alpha,
                              name=self.name)

    def shell_faces(self):
        """
        Computes the shell faces from init data.

        """
        faces = []

        for edge in self.contour3d.primitives:
            faces.extend(edge.revolution(self.axis_point, self.axis, self.angle))

        if not math.isclose(self.angle, volmdlr.TWO_PI, abs_tol=1e-9):
            # Adding contours face to close
            plane1 = surfaces.Plane3D(self.frame)
            face1 = volmdlr.faces.PlaneFace3D(
                plane1, surfaces.Surface2D(self.contour2d, []))
            face2 = face1.rotation(self.axis_point, self.axis, self.angle)
            faces.append(face1)
            faces.append(face2)

        return faces

    def volume(self):
        """
        Volume from Guldin formulae.

        """
        point1 = self.axis_point.plane_projection3d(self.frame.origin, self.frame.u, self.frame.v)
        p1_2d = point1.to_2d(self.axis_point, self.frame.u, self.frame.v)
        p2_3d = self.axis_point + volmdlr.Point3D(self.axis.vector)
        p2_2d = p2_3d.to_2d(self.frame.origin, self.frame.u, self.frame.v)
        axis_2d = curves.Line2D(p1_2d, p2_2d)
        com = self.contour2d.center_of_mass()
        if com is not False:
            dist = axis_2d.point_distance(com)
            return self.angle * dist * self.contour2d.area()
        return 0.

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Revolved profile rotation.

        :param center: rotation center
        :type center: volmdlr.Point3D
        :param axis: rotation axis
        :type axis: volmdlr.Vector3D
        :param angle: angle rotation
        :type angle: float
        :return: a new rotated Revolved Profile
        :rtype: Revolved Profile
        """
        return self.__class__(
            frame=self.frame.rotation(center, axis, angle),
            contour2d=self.contour2d,
            axis_point=self.axis_point.rotation(center, axis, angle),
            axis=self.axis.rotation(center=volmdlr.O3D, axis=axis,
                                    angle=angle),
            angle=self.angle,
            color=self.color, alpha=self.alpha)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Revolved Profile translation.

        :param offset: translation vector.
        :return: A new translated Revolved Profile.
        """
        return self.__class__(
            frame=self.frame.translation(offset),
            contour2d=self.contour2d,
            axis_point=self.axis_point.translation(offset),
            axis=self.axis,
            angle=self.angle,
            color=self.color, alpha=self.alpha)

    def frame_mapping_parameters(self, frame: volmdlr.Frame3D, side: str):
        """Apply transformation to object's parameters."""
        basis = frame.Basis()
        if side == 'old':
            axis = basis.local_to_global_coordinates(self.axis)
        elif side == 'new':
            axis = basis.global_to_local_coordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')

        return axis

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Revolved Profile.

        side = 'old' or 'new'
        """
        axis = self.frame_mapping_parameters(frame, side)
        return RevolvedProfile(
            self.frame.frame_mapping(frame, side),
            self.contour2d,
            self.axis_point.frame_mapping(frame, side),
            axis=axis, angle=self.angle)


class Cylinder(shells.ClosedShell3D):
    """
    Represents a 3D cylinder defined by its frame, radius, and length.
    """
    # pylint: disable=too-many-arguments

    def __init__(
        self,
        frame: volmdlr.Frame3D,
        radius: float,
        length: float,
        color: Tuple[float, float, float] = None,
        alpha: float = 1.0,
        name: str = "",
    ):
        """
        Initializes the Cylinder instance.

        The `Cylinder` class creates a cylinder with the specified radius and length, positioned using the given frame.
        The axis of revolution of the cylinder corresponds to the local z-axis (w-axis) of the provided frame.

        :param frame: The reference frame defining the position and orientation of the cylinder.
            The w-axis of the frame corresponds to the axis of revolution of the cylinder.
        :type frame: volmdlr.Frame3D
        :param radius: The radius of the cylinder.
        :type radius: float
        :param length: The length of the cylinder.
        :type length: float
        :param color: The color of the cylinder as an RGB tuple. Default is None.
        :type color: Tuple[float, float, float], optional
        :param alpha: The opacity of the cylinder (0.0 to 1.0). Default is 1.0.
        :type alpha: float, optional
        :param name: The name of the cylinder. Default is an empty string.
        :type name: str, optional
        """
        self.frame = frame
        self.position = frame.origin
        self.axis = frame.w
        self.radius = radius
        self.length = length

        faces = self.shell_faces()

        shells.ClosedShell3D.__init__(self, faces=faces, color=color, alpha=alpha, name=name)

    def shell_faces(self):
        """
        Computes the shell faces from init data.
        """
        surface3d = surfaces.CylindricalSurface3D(
            self.frame.translation(-self.frame.w * (self.length * 0.5)), self.radius
        )
        cylindrical_face = volmdlr.faces.CylindricalFace3D.from_surface_rectangular_cut(
            surface3d, 0, 2 * math.pi, 0, self.length
        )
        lower_plane = surfaces.Plane3D.from_plane_vectors(
            self.frame.origin.translation(-self.frame.w * (self.length * 0.5)), self.frame.u, self.frame.v
        )
        circle = volmdlr.curves.Circle2D(
            self.position.to_2d(self.frame.origin, self.frame.u, self.frame.v), self.radius
        )
        lower_face = volmdlr.faces.PlaneFace3D(
            lower_plane, surfaces.Surface2D(volmdlr.wires.Contour2D([volmdlr.edges.FullArc2D.from_curve(circle)]), [])
        )
        upper_face = lower_face.translation(self.frame.w * self.length)

        return [lower_face, cylindrical_face, upper_face]

    def get_bounding_box(self) -> volmdlr.core.BoundingBox:
        """
        Computes the bounding box of a cylinder.

        :return: The BoundingBox of the Cylinder.
        :rtype: :class:`volmdlr.core.BoundingBox`
        """
        # This was copied for HollowCylinder. Inheritance removed to avoid problems
        radius = self.radius

        point_a = self.position - self.length / 2 * self.axis
        point_b = self.position + self.length / 2 * self.axis

        dx2 = (point_a[0] - point_b[0])**2
        dy2 = (point_a[1] - point_b[1])**2
        dz2 = (point_a[2] - point_b[2])**2

        if point_a[0] > point_b[0]:
            point_a, point_b = point_b, point_a
        xmin = point_a[0] - (((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5) * radius
        xmax = point_b[0] + (((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5) * radius

        if point_a[1] > point_b[1]:
            point_a, point_b = point_b, point_a
        ymin = point_a[1] - (((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5) * radius
        ymax = point_b[1] + (((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5) * radius

        if point_a[2] > point_b[2]:
            point_a, point_b = point_b, point_a
        zmin = point_a[2] - (((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5) * radius
        zmax = point_b[2] + (((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5) * radius

        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def volume(self) -> float:
        """
        Compute the volume of the cylinder.

        :return: The computed volume of the Cylinder.
        :rtype: float
        """
        return self.length * math.pi * self.radius**2

    @classmethod
    def from_end_points(
        cls,
        point1: volmdlr.Point3D,
        point2: volmdlr.Point3D,
        radius: float,
        color: Tuple[float, float, float] = None,
        alpha: float = 1,
        name: str = "",
    ):
        """
        Create a cylinder from two end points.

        :param point1: The first end point defining the base of the cylinder.
        :type point1: volmdlr.Point3D
        :param point2: The second end point defining the top of the cylinder.
        :type point2: volmdlr.Point3D
        :param radius: The radius of the cylinder.
        :type radius: float
        :param color: The color of the cylinder as an RGB tuple. Default is None.
        :type color: Tuple[float, float, float], optional
        :param alpha: The opacity of the cylinder (0.0 to 1.0). Default is 1.0.
        :type alpha: float, optional
        :param name: The name of the cylinder. Default is an empty string.
        :type name: str, optional

        :return: A Cylinder instance created from the specified end points.
        :rtype: Cylinder
        """
        position = 0.5 * (point1 + point2)
        length = point1.point_distance(point2)
        axis = (point2 - point1).to_vector()
        axis = axis.unit_vector()
        u_vector = axis.deterministic_unit_normal_vector()
        v_vector = axis.cross(u_vector)

        frame = volmdlr.Frame3D(position, u_vector, v_vector, axis)

        return cls(frame=frame, radius=radius, length=length, color=color, alpha=alpha, name=name)

    @classmethod
    def from_extremal_points(
        cls,
        point1: volmdlr.Point3D,
        point2: volmdlr.Point3D,
        radius: float,
        color: Tuple[float, float, float] = None,
        alpha: float = 1,
        name: str = "",
    ):
        """Deprecated class method. Use 'from_end_points' instead."""
        warnings.warn("Deprecated classmethod. Use 'from_end_points' instead.", DeprecationWarning)

        return cls.from_end_points(point1, point2, radius, color, alpha, name)

    @classmethod
    def from_center_point_and_axis(
        cls,
        center_point: volmdlr.Point3D,
        axis: volmdlr.Vector3D,
        radius: float,
        length: float,
        color: Tuple[float, float, float] = None,
        alpha: float = 1,
        name: str = "",
    ) -> 'Cylinder':
        """
        Create a cylinder from a center point, an axis, radius, and length.

        :param center_point: The center point of the cylinder (i.e. the middle point of the axis of the cylinder).
        :type center_point: volmdlr.Point3D
        :param axis: The axis of revolution for the cylinder.
        :type axis: volmdlr.Vector3D
        :param radius: The radius of the cylinder.
        :type radius: float
        :param length: The length of the cylinder.
        :type length: float
        :param color: The color of the cylinder as an RGB tuple. Default is None.
        :type color: Tuple[float, float, float], optional
        :param alpha: The opacity of the cylinder (0.0 to 1.0). Default is 1.0.
        :type alpha: float, optional
        :param name: The name of the cylinder. Default is an empty string.
        :type name: str, optional

        :return: A Cylinder instance created from the specified center point, axis, radius, and length.
        :rtype: Cylinder
        """
        u_vector = axis.deterministic_unit_normal_vector()
        v_vector = axis.cross(u_vector)
        frame = volmdlr.Frame3D(center_point, u_vector, v_vector, axis)

        return cls(frame=frame, radius=radius, length=length, color=color, alpha=alpha, name=name)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float) -> 'Cylinder':
        """
        Cylinder rotation.

        :param center: The rotation center.
        :type center: volmdlr.Point3D
        :param axis: The rotation axis.
        :type axis: volmdlr.Vector3D
        :param angle: The angle of rotation.
        :type angle: float

        :return: A new rotated Cylinder.
        :rtype: Cylinder
        """
        return self.__class__(
            frame=self.frame.rotation(center, axis, angle),
            length=self.length,
            radius=self.radius,
            color=self.color,
            alpha=self.alpha,
            name=self.name,
        )

    def translation(self, offset: volmdlr.Vector3D) -> 'Cylinder':
        """
        Cylinder translation.

        :param offset: The translation vector.
        :type offset: volmdlr.Vector3D

        :return: A new translated Cylinder.
        :rtype: Cylinder
        """
        return self.__class__(
            frame=self.frame.translation(offset),
            length=self.length,
            radius=self.radius,
            color=self.color,
            alpha=self.alpha,
            name=self.name,
        )

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str) -> 'Cylinder':
        """
        Changes frame_mapping and return a new Frame3D.

        side = 'old' or 'new'
        """
        return Cylinder(
            frame=self.frame.frame_mapping(frame, side),
            radius=self.radius,
            length=self.length,
            color=self.color,
            alpha=self.alpha,
            name=self.name,
        )

    def copy(self, deep=True, memo=None) -> 'Cylinder':
        """
        Creates a copy of Cylinder.

        :return: A copy of a current Cylinder.
        :rtype: Cylinder
        """
        return Cylinder(
            frame=self.frame.copy(),
            radius=self.radius,
            length=self.length,
            color=self.color,
            alpha=self.alpha,
            name=self.name,
        )

    def min_distance_to_other_cylinder(self, other_cylinder: 'Cylinder') -> float:
        """
        Compute the minimal distance between two volmdlr cylinders.

        :param other_cylinder: The other cylinder to compute the distance with.
        :type other_cylinder: Cylinder

        :return: The minimal distance between the two 3D cylinders.
        :rtype: float
        """
        # Basic check
        if self.point_belongs(other_cylinder.position) or other_cylinder.point_belongs(self.position):
            return 0.

        # Local frames of cylinders
        frame0 = volmdlr.Frame3D.from_point_and_vector(
            point=self.position, vector=self.axis, main_axis=volmdlr.Z3D
        )
        frame1 = volmdlr.Frame3D.from_point_and_vector(
            point=other_cylinder.position,
            vector=other_cylinder.axis,
            main_axis=volmdlr.Z3D,
        )

        matrix0 = frame0.transfer_matrix()
        x0, y0, z0 = frame0.origin.x, frame0.origin.y, frame0.origin.z
        matrix1 = frame1.transfer_matrix()
        x1, y1, z1 = frame1.origin.x, frame1.origin.y, frame1.origin.z

        # Euclidean distance
        def dist(point0, point1):
            return math.sqrt(
                (point0[0] - point1[0]) ** 2 + (point0[1] - point1[1]) ** 2 + (point0[2] - point1[2]) ** 2
            )

        # Local coordinates to global coordinates
        def to_global_point(point, matrix, origin):
            return [
                matrix.M11 * point[0] + matrix.M12 * point[1] + matrix.M13 * point[2] + origin[0],
                matrix.M21 * point[0] + matrix.M22 * point[1] + matrix.M23 * point[2] + origin[1],
                matrix.M31 * point[0] + matrix.M32 * point[1] + matrix.M33 * point[2] + origin[2],
            ]

        # Objective function
        def objective(x_param):
            point0_ = to_global_point(x_param[:3], matrix0, [x0, y0, z0])
            point1_ = to_global_point(x_param[3:], matrix1, [x1, y1, z1])

            return dist(point0_, point1_)

        # Gradient of objective function
        def gradient_objective(x_param):
            point0_ = to_global_point(x_param[:3], matrix0, [x0, y0, z0])
            point1_ = to_global_point(x_param[3:], matrix1, [x1, y1, z1])

            distance = dist(point0_, point1_)

            return [
                (point0_[0] - point1_[0]) / distance * matrix0.M11
                + (point0_[1] - point1_[1]) / distance * matrix0.M21
                + (point0_[2] - point1_[2]) / distance * matrix0.M31,
                (point0_[0] - point1_[0]) / distance * matrix0.M12
                + (point0_[1] - point1_[1]) / distance * matrix0.M22
                + (point0_[2] - point1_[2]) / distance * matrix0.M32,
                (point0_[0] - point1_[0]) / distance * matrix0.M13
                + (point0_[1] - point1_[1]) / distance * matrix0.M23
                + (point0_[2] - point1_[2]) / distance * matrix0.M33,
                (point1_[0] - point0_[0]) / distance * matrix1.M11
                + (point1_[1] - point0_[1]) / distance * matrix1.M21
                + (point1_[2] - point0_[2]) / distance * matrix1.M31,
                (point1_[0] - point0_[0]) / distance * matrix1.M12
                + (point1_[1] - point0_[1]) / distance * matrix1.M22
                + (point1_[2] - point0_[2]) / distance * matrix1.M32,
                (point1_[0] - point0_[0]) / distance * matrix1.M13
                + (point1_[1] - point0_[1]) / distance * matrix1.M23
                + (point1_[2] - point0_[2]) / distance * matrix1.M33,
            ]

        # Initial vector
        initial_guess = npy.zeros(6)

        # Constraints
        def constraint_radius_0(x):
            # radius of cylinder 0
            return x[0] ** 2 + x[1] ** 2

        def constraint_radius_1(x):
            # radius of cylinder 1
            return x[3] ** 2 + x[4] ** 2

        def gradient_constraint_radius_0(x):
            # gradient of constraint_radius_0
            return [2 * x[0], 2 * x[1], 0, 0, 0, 0]

        def gradient_constraint_radius_1(x):
            # gradient of constraint_radius_1
            return [0, 0, 0, 2 * x[3], 2 * x[4], 0]

        constraints = [
            NonlinearConstraint(
                fun=constraint_radius_0,
                lb=0,
                ub=self.radius**2,
                jac=gradient_constraint_radius_0,
            ),
            NonlinearConstraint(
                fun=constraint_radius_1,
                lb=0,
                ub=other_cylinder.radius**2,
                jac=gradient_constraint_radius_1,
            ),
        ]

        # Bounds
        bounds = Bounds(
            lb=[
                -self.radius,
                -self.radius,
                -self.length / 2,
                -other_cylinder.radius,
                -other_cylinder.radius,
                -other_cylinder.length / 2,
            ],
            ub=[
                self.radius,
                self.radius,
                self.length / 2,
                other_cylinder.radius,
                other_cylinder.radius,
                other_cylinder.length / 2,
            ],
        )

        return minimize(
            fun=objective,
            x0=initial_guess,
            bounds=bounds,
            tol=1e-6,
            constraints=constraints,
            jac=gradient_objective,
        ).fun

    def is_intersecting_other_cylinder(self, other_cylinder: 'Cylinder') -> bool:
        """
        Verifies if two cylinders are intersecting or not.

        :param other_cylinder: The other cylinder to compute if there is an intersection with.
        :type other_cylinder: Cylinder

        :return: True if cylinders are intersecting, False otherwise
        :rtype: bool
        """
        dist = self.min_distance_to_other_cylinder(other_cylinder)

        return dist < 1e-5

    def random_point_inside(self) -> volmdlr.Point3D:
        """
        Return a random point inside a cylinder.

        :return: a random point inside the Cylinder
        :rtype: volmdlr.Point3D
        """
        theta = uniform(0, 2 * math.pi)
        radius = math.sqrt(uniform(0, 1)) * self.radius

        x_local = radius * math.cos(theta)
        y_local = radius * math.sin(theta)
        z_local = uniform(-self.length / 2, self.length / 2)

        local_frame = volmdlr.Frame3D.from_point_and_vector(
            point=self.position, vector=self.axis, main_axis=volmdlr.Z3D
        )

        return local_frame.local_to_global_coordinates(volmdlr.Point3D(x_local, y_local, z_local)).to_point()

    def lhs_points_inside(self, n_points: int) -> List[volmdlr.Point3D]:
        """
        Returns some points inside the cylinder from a LHS samplings.

        :param n_points: The number of points to generate.
        :type n_points: int

        :return: The Latin Hypercube Sampling points inside the cylinder.
        :rtype: list[volmdlr.Point3D]
        """
        local_frame = volmdlr.Frame3D.from_point_and_vector(
            point=self.position, vector=self.axis, main_axis=volmdlr.Z3D
        )

        # sampling point in Cartesian local coordinates
        sampler = qmc.LatinHypercube(d=3, seed=0)
        sample = qmc.scale(
            sampler.random(n=n_points),
            [0, 0, -self.length / 2],
            [1, 2 * math.pi, self.length / 2],
        )

        # converting sampled point in global coordinates volmdlr.Point3D points
        points = []
        for point in sample:
            radius = math.sqrt(point[0]) * self.radius
            theta = point[1]

            x_local = radius * math.cos(theta)
            y_local = radius * math.sin(theta)
            z_local = point[2]

            points.append(
                local_frame.local_to_global_coordinates(volmdlr.Point3D(x_local, y_local, z_local)).to_point()
            )

        return points

    def point_belongs(self, point3d: volmdlr.Point3D, **kwargs) -> bool:
        """
        Check if the point belongs to the cylinder.

        :param point3d: The point to check if it's the cylinder.
        :type point3d: volmdlr.Point3D

        :return: True if the given point is inside the cylinder, False otherwise.
        :rtype: bool
        """
        local_frame = volmdlr.Frame3D.from_point_and_vector(
            point=self.position, vector=self.axis, main_axis=volmdlr.Z3D
        )

        local_point = local_frame.global_to_local_coordinates(point3d)

        return (math.sqrt(local_point.x ** 2 + local_point.y ** 2) <= self.radius) and (
                -self.length / 2 <= local_point.z <= self.length / 2
        )

    def interference_volume_with_other_cylinder(self, other_cylinder: "Cylinder", n_points: int = 1000) -> float:
        """
        Estimation of the interpenetration volume using LHS sampling (inspired by Monte-Carlo method).

        :param other_cylinder: The other cylinder to compute the interference volume with.
        :type other_cylinder: Cylinder
        :param n_points: Optional parameter used for the number of random point used to discretize the cylinder
        :type n_points: int

        :return: An estimation of the interference volume.
        :rtype: float
        """

        # doing the discretization on the smallest cylinder to have better precision
        if self.volume() < other_cylinder.volume():
            smallest_cylinder = self
        else:
            smallest_cylinder = other_cylinder
            other_cylinder = self

        return (
                len(
                    [
                        point
                        for point in smallest_cylinder.lhs_points_inside(n_points)
                        if other_cylinder.point_belongs(point)
                    ]
                )
                / n_points
        ) * smallest_cylinder.volume()


class Cone(shells.ClosedShell3D):
    """
    Represents a 3D cone defined by its frame, radius, and length.
    """
    # pylint: disable=too-many-arguments

    def __init__(
        self,
        frame: volmdlr.Frame3D,
        radius: float,
        length: float,
        color: Tuple[float, float, float] = None,
        alpha: float = 1.0,
        name: str = "",
    ):
        """
        Initializes the Cone instance.

        The `Cone` class creates a cone with the specified radius and length, positioned using the given frame.
        The axis of revolution of the cone corresponds to the local z-axis (w-axis) of the provided frame.
        The top of the cone is oriented according to the local z-axis (w-axis) of the provided frame.
        The center of the cone axis is positioned at the provided frame origin point.

        :param frame: The reference frame defining the position and orientation of the cone.
            The w-axis of the frame corresponds to the axis of revolution of the cone.
        :type frame: volmdlr.Frame3D
        :param radius: The radius of the cone.
        :type radius: float
        :param length: The length of the cone.
        :type length: float
        :param color: The color of the cone as an RGB tuple. Default is None.
        :type color: Tuple[float, float, float], optional
        :param alpha: The opacity of the cone (0.0 to 1.0). Default is 1.0.
        :type alpha: float, optional
        :param name: The name of the cone. Default is an empty string.
        :type name: str, optional
        """
        self.frame = frame
        self.position = frame.origin
        self.axis = frame.w
        self.radius = radius
        self.length = length

        faces = self.shell_faces()

        shells.ClosedShell3D.__init__(self, faces=faces, color=color, alpha=alpha, name=name)

    def shell_faces(self):
        """
        Computes the shell faces from init data.
        """
        conical_surface_frame = self.frame.translation(self.frame.w * (self.length * 0.5))
        conical_surface_frame.w = -conical_surface_frame.w
        surface3d = surfaces.ConicalSurface3D(conical_surface_frame, math.atan(self.radius / self.length))

        conical_face = volmdlr.faces.ConicalFace3D.from_surface_rectangular_cut(
            surface3d, 0, 2 * math.pi, 0, self.length
        )
        lower_plane = surfaces.Plane3D.from_plane_vectors(
            self.frame.origin.translation(-self.frame.w * (self.length * 0.5)), self.frame.u, self.frame.v
        )
        circle = volmdlr.curves.Circle2D(
            self.position.to_2d(self.frame.origin, self.frame.u, self.frame.v), self.radius
        )
        lower_face = volmdlr.faces.PlaneFace3D(
            lower_plane, surfaces.Surface2D(volmdlr.wires.Contour2D([volmdlr.edges.FullArc2D.from_curve(circle)]), [])
        )

        return [lower_face, conical_face]

    def get_bounding_box(self) -> volmdlr.core.BoundingBox:
        """
        Compute the bounding box of the cone.

        A is the point at the basis.
        B is the top.

        :return: The BoundingBox of the Cone.
        :rtype: :class:`volmdlr.core.BoundingBox`
        """
        point_a = self.position - self.length / 2 * self.axis
        point_b = self.position + self.length / 2 * self.axis

        dx2 = (point_a[0] - point_b[0])**2
        dy2 = (point_a[1] - point_b[1])**2
        dz2 = (point_a[2] - point_b[2])**2

        x_bound = (point_a[0] - (((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5) * self.radius,
                   point_a[0] + (((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5) * self.radius, point_b[0])
        xmin = min(x_bound)
        xmax = max(x_bound)

        y_bound = (point_a[1] - (((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5) * self.radius,
                   point_a[1] + (((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5) * self.radius, point_b[1])
        ymin = min(y_bound)
        ymax = max(y_bound)

        z_bound = (point_a[2] - (((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5) * self.radius,
                   point_a[2] + (((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5) * self.radius, point_b[2])
        zmin = min(z_bound)
        zmax = max(z_bound)

        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float) -> 'Cone':
        """
        Cone rotation.

        :param center: The rotation center.
        :type center: volmdlr.Point3D
        :param axis: The rotation axis.
        :type axis: volmdlr.Vector3D
        :param angle: The angle of rotation.
        :type angle: float

        :return: A new rotated Cone.
        :rtype: Cone
        """
        return self.__class__(
            frame=self.frame.rotation(center, axis, angle),
            radius=self.radius,
            length=self.length,
            color=self.color,
            alpha=self.alpha,
        )

    def translation(self, offset: volmdlr.Vector3D) -> 'Cone':
        """
        Cone translation.

        :param offset: The translation vector.
        :type offset: volmdlr.Vector3D

        :return: A new translated Cone.
        :rtype: Cone
        """
        return self.__class__(
            frame=self.frame.translation(offset),
            radius=self.radius,
            length=self.length,
            color=self.color,
            alpha=self.alpha,
        )

    def volume(self) -> float:
        """
        Compute the volume of the cone.

        :return: The computed volume of the cone.
        :rtype: float
        """
        return self.length * math.pi * self.radius**2 / 3

    @classmethod
    def from_center_point_and_axis(
        cls,
        center_point: volmdlr.Point3D,
        axis: volmdlr.Vector3D,
        radius: float,
        length: float,
        color: Tuple[float, float, float] = None,
        alpha: float = 1,
        name: str = "",
    ) -> 'Cone':
        """
        Create a cone from a center point, an axis, radius, and length.

        :param center_point: The center point of the cone (i.e. the middle point of the axis of the cone).
        :type center_point: volmdlr.Point3D
        :param axis: The axis of revolution for the cone.
        :type axis: volmdlr.Vector3D
        :param radius: The radius of the cone.
        :type radius: float
        :param length: The length of the cone.
        :type length: float
        :param color: The color of the cone as an RGB tuple. Default is None.
        :type color: Tuple[float, float, float], optional
        :param alpha: The opacity of the cone (0.0 to 1.0). Default is 1.0.
        :type alpha: float, optional
        :param name: The name of the cone. Default is an empty string.
        :type name: str, optional

        :return: A Cone instance created from the specified center point, axis, radius, and length.
        :rtype: Cone
        """
        u_vector = axis.deterministic_unit_normal_vector()
        v_vector = axis.cross(u_vector)
        frame = volmdlr.Frame3D(center_point, u_vector, v_vector, axis)

        return cls(frame=frame, radius=radius, length=length, color=color, alpha=alpha, name=name)


class HollowCylinder(shells.ClosedShell3D):
    """
    Represents a 3D hollow cylinder defined by its frame, radii, and length.
    """
    # pylint: disable=too-many-arguments

    def __init__(
        self,
        frame: volmdlr.Frame3D,
        inner_radius: float,
        outer_radius: float,
        length: float,
        color: Tuple[float, float, float] = None,
        alpha: float = 1,
        name: str = "",
    ):
        """
        Initializes the HollowCylinder instance.

        The `HollowCylinder` class creates a hollow cylinder with the specified radii and length, positioned using the
        given frame.
        The axis of revolution of the hollow cylinder corresponds to the local z-axis (w-axis) of the provided frame.

        :param frame: The reference frame defining the position and orientation of the hollow cylinder.
            The w-axis of the frame corresponds to the axis of revolution of the hollow cylinder.
        :type frame: volmdlr.Frame3D
        :param inner_radius: The inner radius of the hollow cylinder.
        :type inner_radius: float
        :param outer_radius: The outer radius of the hollow cylinder.
        :type outer_radius: float
        :param length: The length of the hollow cylinder.
        :type length: float
        :param color: The color of the hollow cylinder as an RGB tuple. Default is None.
        :type color: Tuple[float, float, float], optional
        :param alpha: The opacity of the hollow cylinder (0.0 to 1.0). Default is 1.0.
        :type alpha: float, optional
        :param name: The name of the hollow cylinder. Default is an empty string.
        :type name: str, optional
        """
        self.frame = frame
        self.position = frame.origin
        self.axis = frame.w
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.length = length

        faces = self.shell_faces()
        shells.ClosedShell3D.__init__(self, faces=faces, color=color, alpha=alpha, name=name)

    def shell_faces(self):
        """
        Computes the shell faces from init data.
        """
        surface3d_1 = surfaces.CylindricalSurface3D(
            self.frame.translation(-self.frame.w * (self.length * 0.5)), self.outer_radius
        )
        surface3d_2 = surfaces.CylindricalSurface3D(
            self.frame.translation(-self.frame.w * (self.length * 0.5)), self.inner_radius
        )

        cylindrical_face1 = volmdlr.faces.CylindricalFace3D.from_surface_rectangular_cut(
            surface3d_1, 0, 2 * math.pi, 0, self.length
        )
        cylindrical_face2 = volmdlr.faces.CylindricalFace3D.from_surface_rectangular_cut(
            surface3d_2, 0, 2 * math.pi, 0, self.length
        )

        lower_plane = surfaces.Plane3D.from_plane_vectors(
            self.frame.origin.translation(-self.frame.w * (self.length * 0.5)), self.frame.u, self.frame.v
        )

        position_2d = self.position.to_2d(self.frame.origin, self.frame.u, self.frame.v)
        outer_circle = volmdlr.curves.Circle2D(position_2d, self.outer_radius)
        inner_circle = volmdlr.curves.Circle2D(position_2d, self.inner_radius)

        lower_face = volmdlr.faces.PlaneFace3D(
            lower_plane,
            surfaces.Surface2D(
                volmdlr.wires.Contour2D([volmdlr.edges.FullArc2D.from_curve(outer_circle)]),
                [volmdlr.wires.Contour2D([volmdlr.edges.FullArc2D.from_curve(inner_circle)])],
            ),
        )
        upper_face = lower_face.translation(self.frame.w * self.length)

        return [lower_face, cylindrical_face1, cylindrical_face2, upper_face]

    def get_bounding_box(self) -> volmdlr.core.BoundingBox:
        """
        Computes the bounding box of a hollow cylinder.

        :return: The BoundingBox of the HollowCylinder.
        :rtype: :class:`volmdlr.core.BoundingBox`
        """
        radius = self.outer_radius

        point_a = self.position - self.length / 2 * self.axis
        point_b = self.position + self.length / 2 * self.axis

        dx2 = (point_a[0] - point_b[0])**2
        dy2 = (point_a[1] - point_b[1])**2
        dz2 = (point_a[2] - point_b[2])**2

        if point_a[0] > point_b[0]:
            point_a, point_b = point_b, point_a
        xmin = point_a[0] - (((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5) * radius
        xmax = point_b[0] + (((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5) * radius

        if point_a[1] > point_b[1]:
            point_a, point_b = point_b, point_a
        ymin = point_a[1] - (((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5) * radius
        ymax = point_b[1] + (((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5) * radius

        if point_a[2] > point_b[2]:
            point_a, point_b = point_b, point_a
        zmin = point_a[2] - (((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5) * radius
        zmax = point_b[2] + (((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5) * radius

        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def volume(self) -> float:
        """
        Compute the volume of the hollow cylinder.

        :return: The computed volume of the Cylinder.
        :rtype: float
        """
        return self.length * math.pi * (self.outer_radius**2 - self.inner_radius**2)

    @classmethod
    def from_end_points(
        cls,
        point1: volmdlr.Point3D,
        point2: volmdlr.Point3D,
        inner_radius: float,
        outer_radius: float,
        color: Tuple[float, float, float] = None,
        alpha: float = 1,
        name: str = "",
    ):
        """
        Create a hollow cylinder from two end points.

        :param point1: The first end point defining the base of the hollow cylinder.
        :type point1: volmdlr.Point3D
        :param point2: The second end point defining the top of the hollow cylinder.
        :type point2: volmdlr.Point3D
        :param inner_radius: The inner radius of the hollow cylinder.
        :type inner_radius: float
        :param outer_radius: The outer radius of the hollow cylinder.
        :type outer_radius: float
        :param color: The color of the hollow cylinder as an RGB tuple. Default is None.
        :type color: Tuple[float, float, float], optional
        :param alpha: The opacity of the hollow cylinder (0.0 to 1.0). Default is 1.0.
        :type alpha: float, optional
        :param name: The name of the hollow cylinder. Default is an empty string.
        :type name: str, optional

        :return: A HollowCylinder instance created from the specified end points.
        :rtype: HollowCylinder
        """
        position = 0.5 * (point1 + point2)
        length = point1.point_distance(point2)
        axis = (point2 - point1).to_vector()
        axis = axis.unit_vector()
        u_vector = axis.deterministic_unit_normal_vector()
        v_vector = axis.cross(u_vector)

        frame = volmdlr.Frame3D(position, u_vector, v_vector, axis)

        return cls(
            frame=frame,
            inner_radius=inner_radius,
            outer_radius=outer_radius,
            length=length,
            color=color,
            alpha=alpha,
            name=name,
        )

    @classmethod
    def from_extremal_points(
        cls,
        point1: volmdlr.Point3D,
        point2: volmdlr.Point3D,
        inner_radius: float,
        outer_radius: float,
        color: Tuple[float, float, float] = None,
        alpha: float = 1,
        name: str = "",
    ):
        """Deprecated class method. Use 'from_end_points' instead."""
        warnings.warn("Deprecated classmethod. Use 'from_end_points' instead.", DeprecationWarning)

        return cls.from_end_points(point1, point2, inner_radius, outer_radius, color, alpha, name)

    @classmethod
    def from_center_point_and_axis(
        cls,
        center_point: volmdlr.Point3D,
        axis: volmdlr.Vector3D,
        inner_radius: float,
        outer_radius: float,
        length: float,
        color: Tuple[float, float, float] = None,
        alpha: float = 1,
        name: str = "",
    ) -> 'HollowCylinder':
        """
        Create a hollow cylinder from a center point, an axis, radius, and length.

        :param center_point: The center point of the hollow cylinder
            (i.e. the middle point of the axis of the hollow cylinder).
        :type center_point: volmdlr.Point3D
        :param axis: The axis of revolution for the hollow cylinder.
        :type axis: volmdlr.Vector3D
        :param inner_radius: The inner radius of the hollow cylinder.
        :type inner_radius: float
        :param outer_radius: The outer radius of the hollow cylinder.
        :type outer_radius: float
        :param length: The length of the hollow cylinder.
        :type length: float
        :param color: The color of the hollow cylinder as an RGB tuple. Default is None.
        :type color: Tuple[float, float, float], optional
        :param alpha: The opacity of the hollow cylinder (0.0 to 1.0). Default is 1.0.
        :type alpha: float, optional
        :param name: The name of the hollow cylinder. Default is an empty string.
        :type name: str, optional

        :return: A HollowCylinder instance created from the specified center point, axis, radii, and length.
        :rtype: HollowCylinder
        """
        u_vector = axis.deterministic_unit_normal_vector()
        v_vector = axis.cross(u_vector)
        frame = volmdlr.Frame3D(center_point, u_vector, v_vector, axis)

        return cls(
            frame=frame,
            inner_radius=inner_radius,
            outer_radius=outer_radius,
            length=length,
            color=color,
            alpha=alpha,
            name=name,
        )

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float) -> 'HollowCylinder':
        """
        Hollow cylinder rotation.

        :param center: The rotation center.
        :type center: volmdlr.Point3D
        :param axis: The rotation axis.
        :type axis: volmdlr.Vector3D
        :param angle: The angle of rotation.
        :type angle: float

        :return: A new rotated HollowCylinder.
        :rtype: HollowCylinder
        """
        return self.__class__(
            frame=self.frame.rotation(center, axis, angle),
            length=self.length,
            inner_radius=self.inner_radius,
            outer_radius=self.outer_radius,
            color=self.color,
            alpha=self.alpha,
            name=self.name,
        )

    def translation(self, offset: volmdlr.Vector3D) -> 'HollowCylinder':
        """
        Hollow cylinder translation.

        :param offset: The translation vector.
        :type offset: volmdlr.Vector3D

        :return: A new translated HollowCylinder.
        :rtype: HollowCylinder
        """
        return self.__class__(
            frame=self.frame.translation(offset),
            length=self.length,
            inner_radius=self.inner_radius,
            outer_radius=self.outer_radius,
            color=self.color,
            alpha=self.alpha,
            name=self.name,
        )

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str) -> 'HollowCylinder':
        """
        Changes frame_mapping and return a new HollowCylinder.

        side = 'old' or 'new'.
        """
        return HollowCylinder(
            frame=self.frame.frame_mapping(frame, side),
            inner_radius=self.inner_radius,
            outer_radius=self.outer_radius,
            length=self.length,
            color=self.color,
            alpha=self.alpha,
            name=self.name,
        )

    def copy(self, *args, **kwargs) -> 'HollowCylinder':
        """
        Creates a copy of HollowCylinder.

        :return: A copy of a current HollowCylinder.
        :rtype: HollowCylinder
        """
        return HollowCylinder(
            frame=self.frame.copy(),
            inner_radius=self.inner_radius,
            outer_radius=self.outer_radius,
            length=self.length,
            color=self.color,
            alpha=self.alpha,
            name=self.name,
        )


class Sweep(shells.ClosedShell3D):
    """
    Sweep a 2D contour along a Wire3D.

    """

    def __init__(self, contour2d: volmdlr.wires.Contour2D,
                 wire3d: volmdlr.wires.Wire3D, *,
                 color: Tuple[float, float, float] = None, alpha: float = 1,
                 name: str = ''):
        self.contour2d = contour2d
        self.wire3d = wire3d
        self.frames = []
        faces = self.shell_faces()
        shells.ClosedShell3D.__init__(self, faces, color=color,
                                      alpha=alpha, name=name)

    def to_dict(self, *args, **kwargs):
        """Custom serialization for performance."""
        dict_ = dc.DessiaObject.base_dict(self)
        dict_.update({'color': self.color,
                      'alpha': self.alpha,
                      'wire3d': self.wire3d.to_dict(),
                      'contour2d': self.contour2d.to_dict()
                      })

        return dict_

    def shell_faces(self):
        """
        Generates the shell faces.

        For now, it does not take into account rotation of sections.
        """
        def get_sweep_profile_section(profile_section_face, touching_faces):
            primitives_on_section_face = []
            for face in touching_faces:
                for prim in face.outer_contour3d.primitives:
                    if profile_section_face.edge3d_inside(prim):
                        primitives_on_section_face.append(prim)
            contour3d_ = volmdlr.wires.Contour3D.contours_from_edges(primitives_on_section_face)[0]
            return contour3d_

        # End  planar faces
        w = self.wire3d.primitives[0].unit_direction_vector(0.)
        u = self.wire3d.primitives[0].unit_normal_vector(0.)
        if not u:
            u = w.deterministic_unit_normal_vector()
        v = w.cross(u)

        start_plane = surfaces.Plane3D(
            volmdlr.Frame3D(self.wire3d.point_at_abscissa(0.), u, v, w))

        length_last_primitive = self.wire3d.primitives[-1].length()
        w = self.wire3d.primitives[-1].unit_direction_vector(length_last_primitive)
        u = self.wire3d.primitives[-1].unit_normal_vector(length_last_primitive)
        if not u:
            u = w.deterministic_unit_normal_vector()
        v = w.cross(u)

        end_plane = surfaces.Plane3D(
            volmdlr.Frame3D(self.wire3d.primitives[-1].point_at_abscissa(length_last_primitive), u, v, w))

        faces = [volmdlr.faces.PlaneFace3D(start_plane, surfaces.Surface2D(self.contour2d, []))]

        max_brectangle = max(self.contour2d.bounding_rectangle.bounds()) * 4
        cutting_face_contour = volmdlr.wires.Contour2D.from_points(
                    [volmdlr.Point2D(-max_brectangle, -max_brectangle),
                     volmdlr.Point2D(max_brectangle, -max_brectangle),
                     volmdlr.Point2D(max_brectangle, max_brectangle),
                     volmdlr.Point2D(-max_brectangle, max_brectangle)])
        new_faces = []
        last_end_tangent = self.wire3d.primitives[0].unit_direction_vector(0.)
        for i, wire_primitive in enumerate(self.wire3d.primitives):
            start_tangent = wire_primitive.unit_direction_vector(0.)
            normal = wire_primitive.unit_normal_vector(0.)
            if normal is None:
                normal = start_tangent.deterministic_unit_normal_vector()
            tangent_normal_orthonormal = start_tangent.cross(normal)

            if i == 0 or not start_tangent.is_close(last_end_tangent):
                new_faces = wire_primitive.sweep(self.contour2d, None)
            else:
                plane = surfaces.Plane3D(volmdlr.Frame3D(wire_primitive.start, normal,
                                                         tangent_normal_orthonormal, start_tangent))
                cutting_face = volmdlr.faces.PlaneFace3D(plane, surfaces.Surface2D(cutting_face_contour, []))
                section_contour = get_sweep_profile_section(cutting_face, new_faces)
                new_faces = wire_primitive.sweep(self.contour2d, section_contour)

            faces.extend(new_faces)
            last_end_tangent = wire_primitive.unit_direction_vector(wire_primitive.length())
        cutting_face = volmdlr.faces.PlaneFace3D(end_plane, surfaces.Surface2D(cutting_face_contour, []))
        contour3d = get_sweep_profile_section(cutting_face, new_faces)
        contour2d = cutting_face.surface3d.contour3d_to_2d(contour3d)
        end_face = volmdlr.faces.PlaneFace3D(end_plane, surfaces.Surface2D(contour2d, []))
        faces.append(end_face)
        return faces

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Sweep.

        :param side: 'old' or 'new'
        """
        new_wire = self.wire3d.frame_mapping(frame, side)
        return Sweep(self.contour2d, new_wire, color=self.color,
                     alpha=self.alpha, name=self.name)

    def copy(self, deep=True, memo=None):
        """Creates a copy of the Sweep."""
        new_contour2d = self.contour2d.copy()
        new_wire3d = self.wire3d.copy()
        return Sweep(new_contour2d, new_wire3d, color=self.color,
                     alpha=self.alpha, name=self.name)


class Sphere(shells.ClosedShell3D):
    """
    Defines a sphere at a given position & radius.
    """

    def __init__(self, center: volmdlr.Point3D, radius: float,
                 color: Tuple[float, float, float] = None, alpha: float = 1.,
                 name: str = ''):
        self.center = center
        self.radius = radius
        self.position = center

        self.frame = volmdlr.Frame3D(center, volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        spherical_surface = surfaces.SphericalSurface3D(self.frame, self.radius)
        spherical_face = volmdlr.faces.SphericalFace3D.from_surface_rectangular_cut(spherical_surface)
        shells.ClosedShell3D.__init__(self, faces=[spherical_face], color=color, alpha=alpha, name=name)

    def volume(self):
        """
        Computes the volume of the sphere.

        :return: sphere's volume (m³)
        """
        return 4 / 3 * math.pi * self.radius**3

    def point_belongs(self, point3d: volmdlr.Point3D, **kwargs) -> bool:
        """
        Returns if the point belongs to the sphere.

        :param point3d: volmdlr Point3D
        :return: True if the given point is inside the sphere, False otherwise
        """
        return self.center.point_distance(point3d) <= self.radius

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Sphere.

        :param side: 'old' or 'new'
        """
        return Sphere(self.center.frame_mapping(frame, side), self.radius)

    def skin_points(self, resolution: float = 1e-3):
        """Gives points on the skin with respect to a resolution."""
        if resolution > 2 * self.radius:
            return []

        theta = 2 * math.asin(resolution / (2 * self.radius))

        nb_floor = int(math.pi / theta) + 1
        rota_theta = [n * theta for n in range(nb_floor)]

        point1 = self.center + volmdlr.X3D * self.radius

        skin_points = []

        for theta_ in rota_theta:
            pt_floor_init = point1.rotation(self.center, volmdlr.Y3D, theta_)

            if math.isclose(theta_, 0, abs_tol=1e-6) or math.isclose(theta_, math.pi, abs_tol=1e-6):
                skin_points.append(pt_floor_init)

            else:
                center_floor = volmdlr.Point3D(volmdlr.X3D.dot(pt_floor_init),
                                               self.center.y,
                                               self.center.z)

                r_floor = center_floor.point_distance(pt_floor_init)
                theta_floor = resolution / r_floor

                rota_theta_floor = [n * theta_floor for n in range(int(2 * math.pi / theta_floor) + 1)]

                if (2 * math.pi - rota_theta_floor[-1]) / theta_floor <= 0.1:
                    rota_theta_floor.pop()

                for tetha_f in rota_theta_floor:
                    pt_floor = pt_floor_init.rotation(center_floor, volmdlr.X3D, tetha_f)
                    skin_points.append(pt_floor)

        return skin_points

    def inner_points(self, resolution: float = 1e-3):
        """Gives points inside the sphere with a sub-sphere strategy."""
        in_points = [self.center]
        nb_spheres = int(self.radius / resolution)
        if nb_spheres == 0:
            return in_points

        spheres_radius = [n * resolution for n in range(1, nb_spheres + 1)]

        if (self.radius - spheres_radius[-1]) / resolution <= 0.1:
            spheres_radius.pop()

        for srad in spheres_radius:
            in_sphere = Sphere(self.center, srad)
            in_points.extend(in_sphere.skin_points(resolution=resolution))

        return in_points


class Measure3D:
    """
    Used to create a measure between two points in 3D.
    """

    def __init__(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D,
                 color: Tuple[float, float, float] = (1., 0, 0)):
        self.point1, self.point2 = point1, point2
        self.color = color
        self.distance = (point1 - point2).norm()
        # self.bounding_box = self._bounding_box()

    # !!! no eq defined!
    def __hash__(self):
        return hash(self.point1) + hash(self.point2)


class BSplineExtrusion(volmdlr.core.Primitive3D):
    """
    Defines the extrusion of a BSpline.

    :param vectorextru: extrusion vector.
    """

    def __init__(self, obj, vectorextru: volmdlr.Vector3D, name: str = ''):
        self.obj = obj
        vectorextru = vectorextru.unit_vector()
        self.vectorextru = vectorextru
        if obj.__class__ is curves.Ellipse3D:
            self.points = obj.tessel_points
        else:
            self.points = obj.points

        volmdlr.core.Primitive3D.__init__(self, name=name)

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a BSplineExtrusion.

        :param arguments: The arguments of the step primitive. The last element represents the unit_conversion_factor.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding BSplineExtrusion object.
        :rtype: :class:`volmdlr.primitives3d.BSplineExtrusion`
        """
        name = arguments[0][1:-1]
        if object_dict[arguments[1]].__class__ is curves.Ellipse3D:
            ell = object_dict[arguments[1]]
            vectextru = -object_dict[arguments[2]]
            return cls(ell, vectextru, name)

        if object_dict[arguments[1]].__class__ is volmdlr.edges.BSplineCurve3D:
            bsplinecurve = object_dict[arguments[1]]
            vectextru = object_dict[arguments[2]]
            return cls(bsplinecurve, vectextru, name)
        raise NotImplementedError  # to be adapted to bspline
