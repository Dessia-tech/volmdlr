#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common primitives 3D.
"""

import math
import warnings
from random import uniform
from typing import Dict, List, Tuple

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

# import dessia_common.typings as dct

npy.seterr(divide='raise')


class OpenRoundedLineSegments3D(volmdlr.wires.Wire3D,
                                volmdlr.primitives.RoundedLineSegments):
    """
    Defines an open rounded line segments.

    :param points: Points used to draw the wire.
    :type points: List of Point3D.
    :param radius: Radius used to connect different parts of the wire.
    :type radius: {position1(n): float which is the radius linked the n-1 and.
    n+1 points, position2(n+1):...}
    """
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']

    line_class = volmdlr.edges.LineSegment3D
    arc_class = volmdlr.edges.Arc3D

    def __init__(self, points: List[volmdlr.Point3D], radius: Dict[str, float],
                 adapt_radius: bool = False, name: str = ''):
        volmdlr.primitives.RoundedLineSegments.__init__(
            self, points, radius, closed=False, adapt_radius=adapt_radius,
            name='')

        volmdlr.wires.Wire3D.__init__(self, self._primitives(), name)

    def arc_features(self, point_index: int):
        # raise NotImplementedError
        radius = self.radius[point_index]
        pt1, pti, pt2 = self.get_points(point_index)
        dist1 = (pt1 - pti).norm()
        dist2 = (pt2 - pti).norm()
        dist3 = (pt1 - pt2).norm()
        alpha = math.acos(-(dist3**2 - dist1**2 - dist2**2) / (2 * dist1 * dist2)) / 2.
        dist = radius / math.tan(alpha)

        u1 = (pt1 - pti) / dist1
        u2 = (pt2 - pti) / dist2

        p3 = pti + u1 * dist
        p4 = pti + u2 * dist

        n = u1.cross(u2)
        n /= n.norm()
        v1 = u1.cross(n)
        v2 = u2.cross(n)

        line1 = volmdlr.edges.Line3D(p3, p3 + v1)
        line2 = volmdlr.edges.Line3D(p4, p4 + v2)

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

    def rotation_inplace(self, center: volmdlr.Point3D,
                         axis: volmdlr.Vector3D,
                         angle: float):
        """
        OpenRoundedLineSegments3D rotation. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in self.points:
            point.rotation_inplace(center, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        OpenRoundedLineSegments3D translation.

        :param offset: translation vector
        :return: A new translated OpenRoundedLineSegments3D
        """
        return self.__class__([point.translation(offset)
                               for point in self.points],
                              self.radius, self.closed, self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        OpenRoundedLineSegments3D translation. Object is updated inplace.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for point in self.points:
            point.translation_inplace(offset)


class ClosedRoundedLineSegments3D(volmdlr.wires.Contour3D,
                                  OpenRoundedLineSegments3D):
    """
    Defines a closed rounded line segment in 3D.

    :param points: Points used to draw the wire
    :type points: List of Point3D
    :param radius: Radius used to connect different parts of the wire
    :type radius: {position1(n): float which is the radius linked the n-1 and
    n+1 points, position2(n+1):...}
    """
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']

    def __init__(self, points: List[volmdlr.Point3D], radius: float, adapt_radius: bool = False, name: str = ''):
        volmdlr.primitives.RoundedLineSegments.__init__(
                self, points, radius, 'volmdlr.edges.LineSegment3D',
                'volmdlr.edges.Arc3D', closed=True, adapt_radius=adapt_radius,
                name='')

        volmdlr.wires.Contour3D.__init__(self, primitives=self._primitives(), name=name)


class Block(volmdlr.faces.ClosedShell3D):
    """
    Creates a block.

    :param frame: a frame 3D. The origin of the frame is the center of the block,
     the 3 vectors are defining the edges. The frame has not to be orthogonal
    """
    _standalone_in_db = True
    _non_serializable_attributes = ['size', 'bounding_box', 'faces', 'contours', 'plane', 'points', 'polygon2D']
    _non_eq_attributes = ['name', 'color', 'alpha', 'size', 'bounding_box',
                          'faces', 'contours', 'plane', 'points', 'polygon2D']
    _non_hash_attributes = []

    def __init__(self, frame: volmdlr.Frame3D, *,
                 color: Tuple[float, float, float] = None, alpha: float = 1.,
                 name: str = ''):
        self.frame = frame
        self.size = (self.frame.u.norm(),
                     self.frame.v.norm(),
                     self.frame.w.norm())

        faces = self.shell_faces()
        volmdlr.faces.ClosedShell3D.__init__(self, faces, color=color, alpha=alpha, name=name)

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
        return self.size[0] * self.size[1] * self.size[2]

    @classmethod
    def from_bounding_box(cls, bounding_box):
        """
        Transform a bounding box into a block.
        """
        origin = bounding_box.center
        sx = bounding_box.xmax - bounding_box.xmin
        sy = bounding_box.ymax - bounding_box.ymin
        sz = bounding_box.zmax - bounding_box.zmin
        frame = volmdlr.Frame3D(origin, sx * volmdlr.Vector3D(1, 0, 0),
                                sy * volmdlr.Vector3D(0, 1, 0),
                                sz * volmdlr.Vector3D(0, 0, 1))
        return cls(frame=frame)

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
        frame.normalize()
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

        xm_face = volmdlr.faces.Plane3D(xm_frame)\
            .rectangular_cut(-hly, hly, -hlz, hlz)
        xp_face = volmdlr.faces.Plane3D(xp_frame)\
            .rectangular_cut(-hly, hly, -hlz, hlz)
        ym_face = volmdlr.faces.Plane3D(ym_frame)\
            .rectangular_cut(-hlz, hlz, -hlx, hlx)
        yp_face = volmdlr.faces.Plane3D(yp_frame)\
            .rectangular_cut(-hlz, hlz, -hlx, hlx)
        zm_face = volmdlr.faces.Plane3D(zm_frame)\
            .rectangular_cut(-hlx, hlx, -hly, hly)
        zp_face = volmdlr.faces.Plane3D(zp_frame)\
            .rectangular_cut(-hlx, hlx, -hly, hly)

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

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        Block rotation. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.frame.rotation_inplace(center, axis, angle)
        self.faces = self.shell_faces()

    def translation(self, offset: volmdlr.Vector3D):
        """
        Block translation.

        :param offset: translation vector.
        :return: A new translated Block.
        """
        new_frame = self.frame.translation(offset)
        return Block(new_frame, color=self.color,
                     alpha=self.alpha, name=self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Block translation. Object is updated inplace.

        :param offset: translation vector.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.frame.translation_inplace(offset)
        self.faces = self.shell_faces()

    def cut_by_orthogonal_plane(self, plane_3d: volmdlr.faces.Plane3D):
        bouding_box = self.bounding_box
        if plane_3d.frame.w.dot(volmdlr.Vector3D(1, 0, 0)) == 0:
            pass
        elif plane_3d.frame.w.dot(volmdlr.Vector3D(0, 1, 0)) == 0:
            pass
        elif plane_3d.frame.w.dot(volmdlr.Vector3D(0, 0, 1)) == 0:
            pass
        else:
            raise KeyError('plane is not orthogonal either with x, y or z')
        dir1 = plane_3d.frame.u
        dir2 = plane_3d.frame.v
        point_min = volmdlr.Point3D(bouding_box.xmin, bouding_box.ymin,
                                    bouding_box.zmin)
        point_max = volmdlr.Point3D(bouding_box.xmax, bouding_box.ymax,
                                    bouding_box.zmax)
        points = [volmdlr.Point2D(point_min.dot(dir1), point_min.dot(dir2)),
                  volmdlr.Point2D(point_min.dot(dir1), point_max.dot(dir2)),
                  volmdlr.Point2D(point_max.dot(dir1), point_max.dot(dir2)),
                  volmdlr.Point2D(point_max.dot(dir1), point_min.dot(dir2))]
        contour_2d = volmdlr.faces.Surface2D(
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

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_frame = self.frame_mapping_parametres(frame, side)
        self.frame = new_frame
        self.faces = self.shell_faces()

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

    def plot_data(self, x3d, y3d, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        lines = []
        for edge3d in self.edges():
            lines.append(edge3d.plot_data(x3d, y3d, marker, color,
                                          stroke_width, dash, opacity, arrow))

        return lines

    def plot2d(self, x3d, y3d, ax=None):
        """
        Plot 2d with matplotlib.

        """
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = None

        for edge3D in self.edges():
            # edge2D = edge3D.PlaneProjection2D()
            edge3D.plot2d(x3d, y3d, ax=ax)

        return fig, ax


class ExtrudedProfile(volmdlr.faces.ClosedShell3D):
    """
    Extrude a profile given by outer and inner contours.

    TODO: In the future change to a frame and a surface2D and an extrusion vector.
    """
    _non_serializable_attributes = ['faces', 'inner_contours3d',
                                    'outer_contour3d']

    def __init__(self, plane_origin: volmdlr.Point3D,
                 x: volmdlr.Vector3D, y: volmdlr.Vector3D,
                 outer_contour2d: volmdlr.wires.Contour2D,
                 inner_contours2d: List[volmdlr.wires.Contour2D],
                 extrusion_vector: volmdlr.Vector3D,
                 color: Tuple[float, float, float] = None, alpha: float = 1.,
                 name: str = ''):
        self.plane_origin = plane_origin

        self.outer_contour2d = outer_contour2d
        self.outer_contour3d = outer_contour2d.to_3d(plane_origin, x, y)

        self.inner_contours2d = inner_contours2d
        self.extrusion_vector = extrusion_vector
        self.inner_contours3d = []
        self.x = x
        self.y = y
        self.color = color

        bool_areas = []
        for contour in inner_contours2d:
            self.inner_contours3d.append(contour.to_3d(plane_origin, x, y))
            if contour.area() > outer_contour2d.area():
                bool_areas.append(True)
            else:
                bool_areas.append(False)
        if any(bool_areas):
            raise ValueError('At least one inner contour is not contained in outer_contour.')

        faces = self.shell_faces()

        volmdlr.faces.ClosedShell3D.__init__(self, faces, color=color,
                                             alpha=alpha, name=name)

    def to_dict(self, *args, **kwargs):
        """
        Serialize the ExtrudedProfile.

        """
        dict_ = dc.DessiaObject.base_dict(self)
        dict_.update({'color': self.color,
                      'alpha': self.alpha,
                      'plane_origin': self.plane_origin.to_dict(),
                      'outer_contour2d': self.outer_contour2d.to_dict(),
                      'inner_contours2d': [c.to_dict() for c in self.inner_contours2d],
                      'extrusion_vector': self.extrusion_vector.to_dict(),
                      'x': self.x.to_dict(),
                      'y': self.y.to_dict(),
                      })

        return dict_

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of Extruded Profile.

        """
        return self.__class__(plane_origin=self.plane_origin.copy(),
                              x=self.x.copy(),
                              y=self.y.copy(),
                              outer_contour2d=self.outer_contour2d.copy(),
                              inner_contours2d=[c.copy() for c in self.inner_contours2d],
                              extrusion_vector=self.extrusion_vector.copy(),
                              color=self.color,
                              alpha=self.alpha,
                              name=self.name)

    def shell_faces(self):
        """
        Computes the shell faces from init data.

        """
        lower_plane = volmdlr.faces.Plane3D.from_plane_vectors(
            self.plane_origin, self.x, self.y)
        lower_face = volmdlr.faces.PlaneFace3D(
            lower_plane, volmdlr.faces.Surface2D(self.outer_contour2d,
                                                 self.inner_contours2d))

        upper_face = lower_face.translation(self.extrusion_vector)
        lateral_faces = []
        for p in self.outer_contour3d.primitives:
            lateral_faces.extend(p.extrusion(self.extrusion_vector))

        for inner_contour in self.inner_contours3d:
            for p in inner_contour.primitives:
                lateral_faces.extend(p.extrusion(self.extrusion_vector))

        return [lower_face] + [upper_face] + lateral_faces

    # def plot(self, ax=None, color:str='k', alpha:float=1):
    #     if ax is None:
    #         fig, ax = plt.subplots()
    #         ax.set_aspect('equal')
    #     for contour in [self.outer_contour2d]+self.inner_contours2d:
    #         for primitive in contour.primitives:
    #             primitive.plot(ax)
    #     ax.margins(0.1)
    #     return ax

    def area(self):
        areas = self.outer_contour2d.area()
        areas -= sum([contour.area() for contour in self.inner_contours2d])

        # sic=list(npy.argsort(areas))[::-1]# sorted indices of contours
        # area=areas[sic[0]]

        # for i in sic[1:]:
        #     area-=self.contours2D[i].Area()
        return areas

    def volume(self):
        z = self.x.cross(self.y)
        z.normalize()
        return self.area() * self.extrusion_vector.dot(z)

    def frame_mapping_parameters(self, frame: volmdlr.Frame3D,
                                 side: str):
        basis = frame.basis()
        if side == 'old':
            extrusion_vector = basis.local_to_global_coordinates(self.extrusion_vector)
            x = basis.local_to_global_coordinates(self.x)
            y = basis.local_to_global_coordinates(self.y)
        elif side == 'new':
            extrusion_vector = basis.global_to_local_coordinates(self.extrusion_vector)
            x = basis.global_to_local_coordinates(self.x)
            y = basis.global_to_local_coordinates(self.y)
        else:
            raise ValueError('side must be either old or new')
        return extrusion_vector, x, y

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new ExtrudeProfile.

        :param side: = 'old' or 'new'.
        """
        extrusion_vector, x, y = self.frame_mapping_parameters(frame,
                                                               side)
        return ExtrudedProfile(
            self.plane_origin.frame_mapping(frame, side),
            x, y, self.outer_contour2d, self.inner_contours2d,
            extrusion_vector)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        :param side: = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.extrusion_vector, self.x, self.y =\
            self.frame_mapping_parameters(frame, side)
        self.plane_origin.frame_mapping_inplace(frame, side)

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
            plane_origin=self.plane_origin.rotation(center, axis, angle),
            x=self.x.rotation(volmdlr.O3D, axis, angle),
            y=self.y.rotation(volmdlr.O3D, axis, angle),
            outer_contour2d=self.outer_contour2d,
            inner_contours2d=self.inner_contours2d,
            extrusion_vector=self.extrusion_vector.rotation(volmdlr.O3D,
                                                            axis, angle),
            color=self.color, alpha=self.alpha)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        Extruded Profile rotation. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.plane_origin.rotation_inplace(center, axis, angle)
        self.x.rotation_inplace(volmdlr.O3D, axis, angle)
        self.y.rotation_inplace(volmdlr.O3D, axis, angle)
        self.extrusion_vector.rotation_inplace(volmdlr.O3D, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Extruded Profile translation.

        :param offset: translation vector
        :return: A new translated ExtrudedProfile
        """
        return self.__class__(
            plane_origin=self.plane_origin.translation(offset),
            x=self.x, y=self.y,
            outer_contour2d=self.outer_contour2d,
            inner_contours2d=self.inner_contours2d,
            extrusion_vector=self.extrusion_vector,
            color=self.color, alpha=self.alpha)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Extruded profile translation. Object is updated inplace.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.plane_origin.translation_inplace(offset)


class RevolvedProfile(volmdlr.faces.ClosedShell3D):
    """
    Revolve a 2D profile along an axis around a certain angle.

    """
    _non_serializable_attributes = ['faces', 'contour3D']

    def __init__(self, plane_origin: volmdlr.Point3D,
                 x: volmdlr.Vector3D, y: volmdlr.Vector3D,
                 contour2d: volmdlr.wires.Contour2D,
                 axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float = 2 * math.pi, *,
                 color: Tuple[float, float, float] = None, alpha: float = 1,
                 name: str = ''):

        self.contour2d = contour2d
        self.axis_point = axis_point
        self.axis = axis
        self.angle = angle
        self.plane_origin = plane_origin
        self.x = x
        self.y = y
        self.contour3d = self.contour2d.to_3d(plane_origin, x, y)

        faces = self.shell_faces()
        volmdlr.faces.ClosedShell3D.__init__(self, faces, color=color,
                                             alpha=alpha, name=name)

    def to_dict(self, *args, **kwargs):
        """
        Custom to dict for performance.
        """
        dict_ = dc.DessiaObject.base_dict(self)
        dict_.update({'color': self.color,
                      'alpha': self.alpha,
                      'plane_origin': self.plane_origin.to_dict(),
                      'contour2d': self.contour2d.to_dict(),
                      'axis_point': self.axis_point.to_dict(),
                      'x': self.x.to_dict(),
                      'y': self.y.to_dict(),
                      'angle': self.angle,
                      'axis': self.axis.to_dict()
                      })

        return dict_

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of Revolvedprofile.

        """
        return self.__class__(plane_origin=self.plane_origin.copy(),
                              x=self.x.copy(), y=self.y.copy(),
                              contour2d=self.contour2d.copy(deep=deep, memo=memo),
                              axis=self.axis.copy(), angle=self.angle,
                              axis_point=self.axis_point.copy(),
                              color=self.color, alpha=self.alpha,
                              name=self.name)

    def shell_faces(self):
        faces = []

        for edge in self.contour3d.primitives:
            faces.extend(edge.revolution(self.axis_point,
                                         self.axis, self.angle))

        if not math.isclose(self.angle, volmdlr.TWO_PI, abs_tol=1e-9):
            # Adding contours face to close
            w = self.x.cross(self.y)
            plane1 = volmdlr.faces.Plane3D(volmdlr.Frame3D(self.plane_origin,
                                                           self.x,
                                                           self.y,
                                                           w))
            face1 = volmdlr.faces.PlaneFace3D(
                plane1, volmdlr.faces.Surface2D(self.contour2d, []))
            face2 = face1.rotation(self.axis_point, self.axis, self.angle)
            faces.append(face1)
            faces.append(face2)

        return faces

    def volume(self):
        """
        Volume from Guldin formulae.
        """
        point1 = self.axis_point.PlaneProjection3D(self.plane_origin,
                                                   self.x, self.y)
        p1_2d = point1.To2D(self.axis_point, self.x, self.y)
        p2_3d = self.axis_point + volmdlr.Point3D(self.axis.vector)
        p2_2d = p2_3d.To2D(self.plane_origin, self.x, self.y)
        axis_2d = volmdlr.edges.Line2D(p1_2d, p2_2d)
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
            plane_origin=self.plane_origin.rotation(center, axis, angle),
            x=self.x.rotation(center=volmdlr.O3D, axis=axis, angle=angle),
            y=self.y.rotation(center=volmdlr.O3D, axis=axis, angle=angle),
            contour2d=self.contour2d,
            axis_point=self.axis_point.rotation(center, axis, angle),
            axis=self.axis.rotation(center=volmdlr.O3D, axis=axis,
                                    angle=angle),
            angle=self.angle,
            color=self.color, alpha=self.alpha)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        Revolved profile rotation. Object is updated inplace.

        :param center: rotation center.
        :type center: `volmdlr.Point3D`.
        :param axis: rotation axis.
        :type axis: `volmdlr.Vector3D`.
        :param angle: rotation angle.
        :type angle: float.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.plane_origin.rotation_inplace(center, axis, angle)
        self.x.rotation_inplace(center=volmdlr.O3D, axis=axis, angle=angle)
        self.y.rotation_inplace(center=volmdlr.O3D, axis=axis, angle=angle)
        self.axis_point.rotation_inplace(center, axis, angle)
        self.axis.rotation_inplace(volmdlr.O3D, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Revolved Profile translation.

        :param offset: translation vector.
        :return: A new translated Revolved Profile.
        """
        return self.__class__(
            plane_origin=self.plane_origin.translation(offset),
            x=self.x, y=self.y, contour2d=self.contour2d,
            axis_point=self.axis_point.translation(offset),
            axis=self.axis,
            angle=self.angle,
            color=self.color, alpha=self.alpha)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Revolved Profile translation. Object is updated inplace.

        :param offset: translation vector.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.plane_origin.translation_inplace(offset)
        self.axis_point.translation_inplace(offset)

    def frame_mapping_parameters(self, frame: volmdlr.Frame3D,
                                 side: str):
        basis = frame.Basis()
        if side == 'old':
            axis = basis.local_to_global_coordinates(self.axis)
            x = basis.local_to_global_coordinates(self.x)
            y = basis.local_to_global_coordinates(self.y)
        elif side == 'new':
            axis = basis.global_to_local_coordinates(self.axis)
            x = basis.global_to_local_coordinates(self.x)
            y = basis.global_to_local_coordinates(self.y)
        else:
            raise ValueError('side must be either old or new')

        return axis, x, y

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Revolved Profile.

        side = 'old' or 'new'
        """
        axis, x, y = self.frame_mapping_parameters(frame, side)
        return RevolvedProfile(
            self.plane_origin.frame_mapping(frame, side),
            x, y, self.contour2d,
            self.axis_point.frame_mapping(frame, side),
            axis, self.angle)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.axis, self.x, self.y = self.frame_mapping_parameters(frame, side)
        self.plane_origin.frame_mapping_inplace(frame, side)
        self.axis_point.frame_mapping_inplace(frame, side)


class Cylinder(RevolvedProfile):
    """
    Creates a full cylinder with the position, the axis of revolution the radius and the length.
    """

    def __init__(self, position: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 radius: float, length: float,
                 color: Tuple[float, float, float] = None, alpha: float = 1.,
                 name: str = ''):

        self.position = position
        axis.normalize()
        self.axis = axis
        self.radius = radius
        self.length = length
        self.bounding_box = self._bounding_box()

        # Revolved Profile
        point1 = volmdlr.Point2D(-0.5 * self.length, 0.)
        point2 = volmdlr.Point2D(0.5 * self.length, 0.)
        point3 = volmdlr.Point2D(0.5 * self.length, self.radius)
        point4 = volmdlr.Point2D(-0.5 * self.length, self.radius)
        line_seg1 = volmdlr.edges.LineSegment2D(point1, point2)
        line_seg2 = volmdlr.edges.LineSegment2D(point2, point3)
        line_seg3 = volmdlr.edges.LineSegment2D(point3, point4)
        line_seg4 = volmdlr.edges.LineSegment2D(point4, point1)
        contour = volmdlr.wires.Contour2D([line_seg1, line_seg2, line_seg3, line_seg4])
        y = axis.random_unit_normal_vector()
        RevolvedProfile.__init__(self, position, axis, y, contour, position,
                                 axis, color=color, alpha=alpha, name=name)

    def _bounding_box(self):
        """
        Computes the bounding box of a cylinder.

        :return: The BoundingBox
        :rtype: :class:`volmdlr.core.BoundingBox`
        """
        # This was copied for HollowCylinder. Inheritance removed to avoid problems
        radius = self.radius

        point_a = self.position - self.length / 2 * self.axis
        point_b = self.position + self.length / 2 * self.axis

        dx2 = (point_a[0] - point_b[0])**2
        dy2 = (point_a[1] - point_b[1])**2
        dz2 = (point_a[2] - point_b[2])**2

        # kx = ((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5
        # ky = ((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5
        # kz = ((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5

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

    def volume(self):
        """Computes the volume of the cylinder."""
        return self.length * math.pi * self.radius**2

    @classmethod
    def from_extremal_points(cls, point1: volmdlr.Point3D, point2: volmdlr.Point3D,
                             radius: float,
                             color: Tuple[float, float, float] = None, alpha: float = 1,
                             name: str = ''):
        position = 0.5 * (point1 + point2)
        length = point1.point_distance(point2)
        axis = point2 - point1
        axis.normalize()
        return cls(position, axis, radius, length=length,
                   color=color, alpha=alpha, name=name)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Cylinder rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Cylinder
        """
        return self.__class__(
            position=self.position.rotation(center, axis, angle),
            axis=self.axis.rotation(volmdlr.O3D, axis, angle),
            length=self.length, radius=self.radius)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        Cylinder rotation. Object is updated inplace.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: rotation angle.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.position.rotation_inplace(center, axis, angle)
        self.axis.rotation_inplace(volmdlr.O3D, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Cylinder translation.

        :param offset: translation vector
        :return: A new translated Cylinder
        """
        return self.__class__(
            position=self.position.translation(offset),
            axis=self.axis, length=self.length, radius=self.radius)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Cylinder translation. Object is updated inplace.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.position.translation_inplace(offset)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Frame3D.

        side = 'old' or 'new'
        """
        basis = frame.basis()
        if side == 'old':
            axis = basis.local_to_global_coordinates(self.axis)
        elif side == 'new':
            axis = basis.global_to_local_coordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')
        return Cylinder(self.position.frame_mapping(frame, side),
                        axis, self.radius, self.length,
                        color=self.color, alpha=self.alpha)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        basis = frame.basis()
        if side == 'old':
            axis = basis.local_to_global_coordinates(self.axis)
        elif side == 'new':
            axis = basis.global_to_local_coordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')
        self.position.frame_mapping_inplace(frame, side)
        self.axis = axis

    def to_dict(self, use_pointers: bool = False, memo: bool = None, path: str = '#'):
        """
        Call to DessiaObject.to_dict to avoid calling the to_dict of the inherited class Revolved Profile.
        """
        return dc.DessiaObject.to_dict(self, use_pointers, memo, path)

    def copy(self, deep=True, memo=None):
        """
        Creates a copy of Cylinder.
        """
        new_position = self.position.copy()
        new_axis = self.axis.copy()
        return Cylinder(new_position, new_axis, self.radius, self.length,
                        color=self.color, alpha=self.alpha, name=self.name)

    def min_distance_to_other_cylinder(self, other_cylinder: 'Cylinder') -> float:
        """
        Compute the minimal distance between two volmdlr cylinders.

        :param other_cylinder: volmdlr Cylinder
        :return: minimal distance between two 3D cylinders
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

        :param other_cylinder: volmdlr Cylinder
        :return: boolean, True if cylinders are intersecting, False otherwise
        """
        dist = self.min_distance_to_other_cylinder(other_cylinder)

        return dist < 1e-5

    def random_point_inside(self) -> volmdlr.Point3D:
        """
        Gets a random point inside a cylinder.

        :return: a random point inside the Cylinder
        """
        theta = uniform(0, 2 * math.pi)
        radius = math.sqrt(uniform(0, 1)) * self.radius

        x_local = radius * math.cos(theta)
        y_local = radius * math.sin(theta)
        z_local = uniform(-self.length / 2, self.length / 2)

        local_frame = volmdlr.Frame3D.from_point_and_vector(
            point=self.position, vector=self.axis, main_axis=volmdlr.Z3D
        )

        return local_frame.local_to_global_coordinates(volmdlr.Point3D(x_local, y_local, z_local))

    def lhs_points_inside(self, n_points: int) -> List[volmdlr.Point3D]:
        """
        Returns some points inside the cylinder from a LHS samplings.

        :param n_points: number of points
        :return: Latin hypercube sampling points inside the cylinder
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
                local_frame.local_to_global_coordinates(volmdlr.Point3D(x_local, y_local, z_local))
            )

        return points

    def point_belongs(self, point3d: volmdlr.Point3D, **kwargs) -> bool:
        """
        Returns if the point belongs to the cylinder.

        :param point3d: volmdlr Point3D
        :return: True if the given point is inside the cylinder, False otherwise
        """
        local_frame = volmdlr.Frame3D.from_point_and_vector(
            point=self.position, vector=self.axis, main_axis=volmdlr.Z3D
        )

        local_point = local_frame.global_to_local_coordinates(point3d)

        return (math.sqrt(local_point.x ** 2 + local_point.y ** 2) <= self.radius) and (
                -self.length / 2 <= local_point.z <= self.length / 2
        )

    def interference_volume_with_other_cylinder(
            self, other_cylinder: "Cylinder", n_points: int = 1000
    ) -> float:
        """
        Estimation of the interpenetration volume using LHS sampling (inspired by Monte-Carlo method).

        :param other_cylinder: volmdlr Cylinder
        :param n_points: optional parameter used for the number of random point used to discretize the cylinder
        :return: an estimation of the interference volume
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
                        # for point in (smallest_cylinder.random_point_inside() for _ in range(n_points))
                        for point in smallest_cylinder.lhs_points_inside(n_points)
                        if other_cylinder.point_belongs(point)
                    ]
                )
                / n_points
        ) * smallest_cylinder.volume()


class Cone(RevolvedProfile):
    """
    Defines a cone at a given position & axis.
    """

    def __init__(self, position: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 radius: float, length: float,
                 color: Tuple[float, float, float] = None, alpha: float = 1.,
                 name: str = ''):

        self.position = position
        axis.normalize()
        self.axis = axis
        self.radius = radius
        self.length = length
        self.bounding_box = self._bounding_box()

        # Revolved Profile
        point1 = volmdlr.Point2D(0., 0.)
        point2 = volmdlr.Point2D(0., self.radius)
        point3 = volmdlr.Point2D(self.length, 0.)

        contour = volmdlr.wires.Contour2D([volmdlr.edges.LineSegment2D(point1, point2),
                                           volmdlr.edges.LineSegment2D(point2, point3),
                                           volmdlr.edges.LineSegment2D(point3, point1)])
        y = axis.random_unit_normal_vector()
        RevolvedProfile.__init__(self, position, axis, y, contour, position,
                                 axis, color=color, alpha=alpha, name=name)

    def _bounding_box(self):
        """
        A is the point at the basis.

        B is the top.
        """
        point_a = self.position - self.length / 2 * self.axis
        point_b = self.position + self.length / 2 * self.axis

        dx2 = (point_a[0] - point_b[0])**2
        dy2 = (point_a[1] - point_b[1])**2
        dz2 = (point_a[2] - point_b[2])**2

        # kx = ((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5
        # ky = ((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5
        # kz = ((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5

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

    def translation(self, offset: volmdlr.Vector3D):
        """
        Cone translation.

        :param offset: translation vector
        :return: A new translated Cone
        """
        return self.__class__(position=self.position.translation(offset),
                              axis=self.axis,
                              radius=self.radius,
                              length=self.length,
                              color=self.color,
                              alpha=self.alpha)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Plane3D translation. Object is updated inplace.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.position.translation_inplace(offset)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Cone rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Cone
        """
        return self.__class__(position=self.position.rotation(
            center, axis, angle), axis=self.axis.rotation(center, axis, angle),
            radius=self.radius, length=self.length, color=self.color,
            alpha=self.alpha)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        Cone rotation. Object is updated inplace.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: rotation angle.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.position.rotation_inplace(center, axis, angle)
        self.axis.rotation_inplace(center, axis, angle)

    def volume(self):
        """
        Returns the volume of the cone.

        """
        return self.length * math.pi * self.radius**2 / 3


class HollowCylinder(RevolvedProfile):
    """
    Creates a hollow cylinder with the position, the axis of revolution the inner and outer radius and the length.

    """

    def __init__(self, position: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 inner_radius: float, outer_radius: float, length: float,
                 color: Tuple[float, float, float] = None, alpha: float = 1,
                 name: str = ''):
        self.position = position
        axis.normalize()
        self.axis = axis
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.length = length

        # Revolved Profile
        point1 = volmdlr.Point2D(-0.5 * self.length, self.inner_radius)
        point2 = volmdlr.Point2D(0.5 * self.length, self.inner_radius)
        point3 = volmdlr.Point2D(0.5 * self.length, self.outer_radius)
        point4 = volmdlr.Point2D(-0.5 * self.length, self.outer_radius)
        line_seg1 = volmdlr.edges.LineSegment2D(point1, point2)
        line_seg2 = volmdlr.edges.LineSegment2D(point2, point3)
        line_seg3 = volmdlr.edges.LineSegment2D(point3, point4)
        line_seg4 = volmdlr.edges.LineSegment2D(point4, point1)
        contour = volmdlr.wires.Contour2D([line_seg1, line_seg2, line_seg3, line_seg4])
        y = axis.random_unit_normal_vector()
        # contour.plot()
        RevolvedProfile.__init__(self, position, axis, y, contour, position,
                                 axis, color=color, alpha=alpha, name=name)

    def _bounding_box(self):

        radius = self.outer_radius

        point_a = self.position - self.length / 2 * self.axis
        point_b = self.position + self.length / 2 * self.axis

        dx2 = (point_a[0] - point_b[0])**2
        dy2 = (point_a[1] - point_b[1])**2
        dz2 = (point_a[2] - point_b[2])**2

        # kx = ((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5
        # ky = ((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5
        # kz = ((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5

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

    def volume(self):
        return self.length * math.pi * (self.outer_radius**2
                                        - self.inner_radius**2)

    def copy(self, *args, **kwargs):
        """
        Creates a copy of HollowCylinder.

        """
        new_position = self.position.copy()
        new_axis = self.axis.copy()
        return HollowCylinder(new_position, new_axis, self.inner_radius, self.outer_radius, self.length,
                              color=self.color, alpha=self.alpha, name=self.name)

    @classmethod
    def from_extremal_points(cls, point1: volmdlr.Point3D, point2: volmdlr.Point3D,
                             inner_radius: float, outer_radius: float,
                             color: Tuple[float, float, float] = None, alpha: float = 1,
                             name: str = ''):
        position = 0.5 * (point1 + point2)
        length = point1.point_distance(point2)
        axis = point2 - point1
        axis.normalize()
        return cls(position, axis, inner_radius=inner_radius, outer_radius=outer_radius, length=length,
                   color=color, alpha=alpha, name=name)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Hollow cylinder rotation.

        :param center: rotation center.
        :type center: volmdlr.Point3D
        :param axis: rotation axis.
        :type axis: volmdlr.Vector3D
        :param angle: angle rotation.
        :type angle: float
        :return: a new rotated HollowCylinder.
        :rtype: HollowCylinder
        """
        return self.__class__(
            position=self.position.rotation(center, axis, angle),
            axis=self.axis.rotation(volmdlr.O3D, axis, angle),
            length=self.length, inner_radius=self.inner_radius,
            outer_radius=self.outer_radius)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        Hollow cylinder rotation. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.position.rotation_inplace(center, axis, angle)
        self.axis.rotation_inplace(volmdlr.O3D, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Hollow cylinder translation.

        :param offset: translation vector.
        :return: A new translated HollowCylinder.
        """
        return self.__class__(
            position=self.position.translation(offset), axis=self.axis,
            length=self.length, inner_radius=self.inner_radius,
            outer_radius=self.outer_radius)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Hollow cylinder translation. Object is updated in-place.

        :param offset: translation vector.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.position.translation_inplace(offset)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new HollowCylinder.

        side = 'old' or 'new'.
        """
        basis = frame.basis()
        if side == 'old':
            axis = basis.local_to_global_coordinates(self.axis)
        elif side == 'new':
            axis = basis.global_to_local_coordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')

        return HollowCylinder(
            position=self.position.frame_mapping(frame, side),
            axis=axis, inner_radius=self.inner_radius,
            outer_radius=self.outer_radius, length=self.length)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        side = 'old' or 'new'.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        basis = frame.basis()
        if side == 'old':
            axis = basis.local_to_global_coordinates(self.axis)
        elif side == 'new':
            axis = basis.global_to_local_coordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')
        self.position.frame_mapping_inplace(frame, side)
        self.axis = axis


class Sweep(volmdlr.faces.ClosedShell3D):
    """
    Sweep a 2D contour along a Wire3D.

    """

    def __init__(self, contour2d: List[volmdlr.wires.Contour2D],
                 wire3d: volmdlr.wires.Wire3D, *,
                 color: Tuple[float, float, float] = None, alpha: float = 1,
                 name: str = ''):
        self.contour2d = contour2d
        self.wire3d = wire3d
        self.frames = []

        faces = self.shell_faces()
        volmdlr.faces.ClosedShell3D.__init__(self, faces, color=color,
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

        For now it does not take into account rotation of sections.
        """

        # End  planar faces
        w = self.wire3d.primitives[0].unit_direction_vector(0.)
        u = self.wire3d.primitives[0].unit_normal_vector(0.)
        if not u:
            u = w.deterministic_unit_normal_vector()
        v = w.cross(u)

        start_plane = volmdlr.faces.Plane3D(
            volmdlr.Frame3D(self.wire3d.point_at_abscissa(0.), u, v, w)
        )

        l_last_primitive = self.wire3d.primitives[-1].length()
        w = self.wire3d.primitives[-1].unit_direction_vector(l_last_primitive)
        u = self.wire3d.primitives[-1].unit_normal_vector(l_last_primitive)
        if not u:
            u = w.deterministic_unit_normal_vector()
        v = w.cross(u)

        end_plane = volmdlr.faces.Plane3D(
            volmdlr.Frame3D(self.wire3d.primitives[-1].point_at_abscissa(
                l_last_primitive),
                            u, v, w))

        faces = [volmdlr.faces.PlaneFace3D(
            start_plane,
            volmdlr.faces.Surface2D(self.contour2d, [])),
                 volmdlr.faces.PlaneFace3D(
            end_plane,
            volmdlr.faces.Surface2D(self.contour2d, []))]

        for wire_primitive in self.wire3d.primitives:
            # tangent, normal = wire_primitive.frenet(0.)
            tangent = wire_primitive.unit_direction_vector(0.)
            normal = wire_primitive.unit_normal_vector(0.)

            if normal is None:
                normal = tangent.deterministic_unit_normal_vector()
            n2 = tangent.cross(normal)
            contour3d = self.contour2d.to_3d(wire_primitive.start, normal, n2)

            if wire_primitive.__class__ is volmdlr.edges.LineSegment3D:
                for contour_primitive in contour3d.primitives:
                    faces.extend(contour_primitive.extrusion(
                        wire_primitive.length()
                        * wire_primitive.unit_direction_vector()))
            elif wire_primitive.__class__ is volmdlr.edges.Arc3D:
                for contour_primitive in contour3d.primitives:
                    faces.extend(contour_primitive.revolution(
                        wire_primitive.center,
                        wire_primitive.normal,
                        wire_primitive.angle))
            elif wire_primitive.__class__ is volmdlr.wires.Circle3D:
                for contour_primitive in contour3d.primitives:
                    faces.extend(contour_primitive.revolution(
                        wire_primitive.center,
                        wire_primitive.normal,
                        volmdlr.TWO_PI))

            elif wire_primitive.__class__ is volmdlr.edges.BSplineCurve3D or \
                    wire_primitive.__class__ is volmdlr.edges.BezierCurve3D:

                tangents = []
                for k, _ in enumerate(wire_primitive.points):
                    position = k / (len(wire_primitive.points) - 1)
                    tangents.append(wire_primitive.tangent(position))

                circles = []
                for pt, tan in zip(wire_primitive.points, tangents):
                    # TODO: replace circle by real contour!
                    circles.append(volmdlr.wires.Circle3D.from_center_normal(center=pt,
                                                                             normal=tan,
                                                                             radius=self.contour2d.radius))

                polys = [volmdlr.wires.ClosedPolygon3D(c.discretization_points()) for c in circles]

                size_v, size_u = len(polys[0].points) + 1, len(polys)
                degree_u, degree_v = 3, 3

                points_3d = []
                for poly in polys:
                    points_3d.extend(poly.points)
                    points_3d.append(poly.points[0])

                bezier_surface3d = volmdlr.faces.BezierSurface3D(degree_u,
                                                                 degree_v,
                                                                 points_3d,
                                                                 size_u,
                                                                 size_v)

                outer_contour = volmdlr.wires.Contour2D([volmdlr.edges.LineSegment2D(volmdlr.O2D, volmdlr.X2D),
                                                         volmdlr.edges.LineSegment2D(
                                                             volmdlr.X2D, volmdlr.X2D + volmdlr.Y2D),
                                                         volmdlr.edges.LineSegment2D(
                                                             volmdlr.X2D + volmdlr.Y2D, volmdlr.Y2D),
                                                         volmdlr.edges.LineSegment2D(volmdlr.Y2D, volmdlr.O2D)])
                surf2d = volmdlr.faces.Surface2D(outer_contour, [])

                bsface3d = volmdlr.faces.BSplineFace3D(bezier_surface3d, surf2d)
                faces.append(bsface3d)

            else:
                raise NotImplementedError(f'Unimplemented primitive for sweep: {wire_primitive.__class__.__name__}')

        return faces

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Sweep.

        :param side: 'old' or 'new'
        """
        new_wire = self.wire3d.frame_mapping(frame, side)
        return Sweep(self.contour2d, new_wire, color=self.color,
                     alpha=self.alpha, name=self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        :param side: 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.wire3d.frame_mapping_inplace(frame, side)
        for face in self.faces:
            face.frame_mapping_inplace(frame, side)

    def copy(self, deep=True, memo=None):
        """Creates a copy of the Sweep."""
        new_contour2d = self.contour2d.copy()
        new_wire3d = self.wire3d.copy()
        return Sweep(new_contour2d, new_wire3d, color=self.color,
                     alpha=self.alpha, name=self.name)


# class Sphere(volmdlr.Primitive3D):
class Sphere(RevolvedProfile):
    """
    Defines a sphere at a given position & radius.

    """

    def __init__(self, center: volmdlr.Point3D, radius: float,
                 color: Tuple[float, float, float] = None, alpha: float = 1.,
                 name: str = ''):
        self.center = center
        self.radius = radius
        self.position = center

        # Revolved Profile for complete sphere
        s = volmdlr.Point2D(-self.radius, 0.01 * self.radius)
        i = volmdlr.Point2D(0, 1.01 * self.radius)
        e = volmdlr.Point2D(self.radius, 0.01 * self.radius)  # Not coherent but it works at first, to change !!

        contour = volmdlr.wires.Contour2D([
            volmdlr.edges.Arc2D(s, i, e), volmdlr.edges.LineSegment2D(s, e)])

        axis = volmdlr.X3D
        y = axis.random_unit_normal_vector()
        RevolvedProfile.__init__(self, center, axis, y, contour, center, axis,
                                 color=color, alpha=alpha, name=name)

    def volume(self):
        return 4 / 3 * math.pi * self.radius**3

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Sphere.

        :param side: 'old' or 'new'
        """
        return Sphere(self.center.frame_mapping(frame, side), self.radius)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.center.frame_mapping_inplace(frame, side)

    def to_point_skin(self, resolution: float = 1e-3):
        if resolution > 2 * self.radius:
            return []

        theta = 2 * math.asin(resolution / (2 * self.radius))

        nb_floor = int(math.pi / theta) + 1
        rota_theta = [n * theta for n in range(nb_floor)]

        point1 = self.center + volmdlr.X3D * self.radius
        rota_axis = volmdlr.Y3D

        skin_points = []

        for theta_ in rota_theta:
            pt_floor_init = point1.rotation(self.center, rota_axis, theta_)

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

    def to_point_in(self, resolution: float = 1e-3):
        in_points = [self.center]
        nb_spheres = int(self.radius / resolution)
        if nb_spheres == 0:
            return in_points

        spheres_radius = [n * resolution for n in range(1, nb_spheres + 1)]

        if (self.radius - spheres_radius[-1]) / resolution <= 0.1:
            spheres_radius.pop()

        for srad in spheres_radius:
            in_sphere = Sphere(self.center, srad)
            in_points.extend(in_sphere.to_point_skin(resolution=resolution))

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
        self.bounding_box = self._bounding_box()

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
        vectorextru.normalize()
        self.vectorextru = vectorextru
        if obj.__class__ is volmdlr.wires.Ellipse3D:
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
        if object_dict[arguments[1]].__class__ is volmdlr.wires.Ellipse3D:
            ell = object_dict[arguments[1]]
            vectextru = -object_dict[arguments[2]]
            return cls(ell, vectextru, name)

        if object_dict[arguments[1]].__class__ is volmdlr.edges.BSplineCurve3D:
            bsplinecurve = object_dict[arguments[1]]
            vectextru = object_dict[arguments[2]]
            return cls(bsplinecurve, vectextru, name)
        raise NotImplementedError  # a adapter pour les bpsline
