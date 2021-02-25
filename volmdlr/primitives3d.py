#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common primitives 3D
"""

import math

import numpy as npy
npy.seterr(divide='raise')

import volmdlr
import volmdlr.core
import volmdlr.primitives
import volmdlr.faces
from typing import Tuple, List, Dict

# import dessia_common.typings as dct
import matplotlib.pyplot as plt


class OpenRoundedLineSegments3D(volmdlr.wires.Wire3D, volmdlr.primitives.RoundedLineSegments):
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    
    def __init__(self, points:List[volmdlr.Point3D], radius:Dict[str, float],
                 adapt_radius:bool=False, name:str=''):
        primitives = volmdlr.primitives.RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.edges.LineSegment3D,
                                                  volmdlr.edges.Arc3D,
                                                  closed=False,
                                                  adapt_radius=adapt_radius,
                                                  name='')

        volmdlr.wires.Wire3D.__init__(self, primitives, name)

    def arc_features(self, ipoint):
        radius = self.radius[ipoint]
        if self.closed:
            if ipoint == 0:
                pt1 = self.points[-1]
            else:
                pt1 = self.points[(ipoint) -1]
            pti = self.points[ipoint]
            if ipoint < self.npoints-1:
                pt2 = self.points[ipoint+1]
            else:
                pt2 = self.points[0]
        else:
            pt1 = self.points[ipoint - 1]
            pti = self.points[ipoint]
            pt2 = self.points[ipoint + 1]

        dist1 = (pt1 - pti).norm()
        dist2 = (pt2 - pti).norm()
        dist3 = (pt1 - pt2).norm()
        alpha = math.acos(-(dist3**2-dist1**2-dist2**2)/(2*dist1*dist2))/2.
        dist = radius/math.tan(alpha)

        u1 = (pt1 - pti) / dist1
        u2 = (pt2 - pti) / dist2

        p3 = pti + u1*dist
        p4 = pti + u2*dist

        n = u1.cross(u2)
        n /= n.norm()
        v1 = u1.cross(n)
        v2 = u2.cross(n)

        l1 = volmdlr.edges.Line3D(p3, p3+v1)
        l2 = volmdlr.edges.Line3D(p4, p4+v2)
        c, _ = l1.minimum_distance_points(l2)

        u3 = u1 + u2# mean of v1 and v2
        u3 /= u3.norm()

        interior = c - u3 * radius
        return p3, interior, p4, dist, alpha


    def rotation(self, center, angle, copy=True):
        if copy:
            return self.__class__([p.rotation(center, angle, copy=True)\
                                          for p in self.points],
                                         self.radius, self.closed, self.name)
        else:
            self.__init__([p.rotation(center, angle, copy=True)\
                           for p in self.points],
                          self.radius, self.closed, self.name)

    def translation(self, offset, copy=True):
        if copy:
            return self.__class__([p.translation(offset, copy=True)\
                                          for p in self.points],
                                         self.radius, self.closed, self.name)
        else:
            self.__init__([p.translation(offset, copy=True)\
                           for p in self.points],
                          self.radius, self.closed, self.name)


class ClosedRoundedLineSegments3D(volmdlr.wires.Contour3D,
                                  OpenRoundedLineSegments3D):
    """
    :param points: Points used to draw the wire 
    :type points: List of Point3D
    :param radius: Radius used to connect different parts of the wire
    :type radius: {position1(n): float which is the radius linked the n-1 and n+1 points, position2(n+1):...}
    """
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    
    def __init__(self, points, radius, adapt_radius=False, name=''):
        primitives = volmdlr.primitives.RoundedLineSegments.__init__(self,
                                                                     points,
                                                                     radius,
                                                                     volmdlr.edges.LineSegment3D,
                                                                     volmdlr.edges.Arc3D,
                                                                     closed=True,
                                                                     adapt_radius=adapt_radius,
                                                                     name='')

        volmdlr.wires.Wire3D.__init__(self, primitives, name)



class Block(volmdlr.faces.ClosedShell3D):
    _standalone_in_db = True
    _generic_eq = True
    _non_serializable_attributes  = ['size']
    _non_eq_attributes = ['name', 'color', 'alpha', 'size', 'bounding_box', 'faces', 'contours',
                          'plane', 'points', 'polygon2D']
    _non_hash_attributes = []

    """
    Creates a block
    :param frame: a frame 3D. The origin of the frame is the center of the block,
     the 3 vectors are defining the edges. The frame has not to be orthogonal
    """
    def __init__(self, frame:volmdlr.Frame3D, *,
                 color:Tuple[float, float, float]=None, alpha:float=1.,
                 name:str=''):
        self.frame = frame
        self.size = (self.frame.u.norm(), self.frame.v.norm(), self.frame.w.norm())

        faces = self.shell_faces()
        volmdlr.faces.OpenShell3D.__init__(self, faces,  color=color, alpha=alpha, name=name)

    # def __hash__(self):
    #     return hash(self.frame)

    def volume(self):
        return self.size[0]*self.size[1]*self.size[2]

    @classmethod
    def from_bounding_box(cls, bounding_box):
        bb = bounding_box
        xmin, xmax, ymin, ymax, zmin, zmax = bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax
        origin = bb.center
        sx, sy, sz = xmax - xmin, ymax - ymin, zmax - zmin
        frame = volmdlr.Frame3D(origin, sx * volmdlr.Vector3D(1, 0, 0),
                           sy * volmdlr.Vector3D(0, 1, 0), sz * volmdlr.Vector3D(0, 0, 1))
        return cls(frame = frame)

    def vertices(self):
        return [self.frame.origin - 0.5*self.frame.u - 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin - 0.5*self.frame.u + 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u + 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u - 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin - 0.5*self.frame.u - 0.5*self.frame.v + 0.5*self.frame.w,
                self.frame.origin - 0.5*self.frame.u + 0.5*self.frame.v + 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u + 0.5*self.frame.v + 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u - 0.5*self.frame.v + 0.5*self.frame.w]

    def edges(self):
        p1, p2, p3, p4, p5, p6, p7, p8 = self.vertices()
        return [volmdlr.edges.LineSegment3D(p1.copy(), p2.copy()),
                volmdlr.edges.LineSegment3D(p2.copy(), p3.copy()),
                volmdlr.edges.LineSegment3D(p3.copy(), p4.copy()),
                volmdlr.edges.LineSegment3D(p4.copy(), p1.copy()),
                volmdlr.edges.LineSegment3D(p5.copy(), p6.copy()),
                volmdlr.edges.LineSegment3D(p6.copy(), p7.copy()),
                volmdlr.edges.LineSegment3D(p7.copy(), p8.copy()),
                volmdlr.edges.LineSegment3D(p8.copy(), p5.copy()),
                volmdlr.edges.LineSegment3D(p1.copy(), p5.copy()),
                volmdlr.edges.LineSegment3D(p2.copy(), p6.copy()),
                volmdlr.edges.LineSegment3D(p3.copy(), p7.copy()),
                volmdlr.edges.LineSegment3D(p4.copy(), p8.copy())]

    def face_contours3d(self):
        e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12 = self.edges()
        e5_switch = e5.reverse()
        e6_switch = e6.reverse()
        e7_switch = e7.reverse()
        e8_switch = e8.reverse()
        e9_switch = e9.reverse()
        e10_switch = e10.reverse()
        e11_switch = e11.reverse()
        e12_switch = e12.reverse()
        return [volmdlr.wires.Contour3D([e1.copy(), e2.copy(), e3.copy(), e4.copy()]),
                volmdlr.wires.Contour3D([e5.copy(), e6.copy(), e7.copy(), e8.copy()]),
                # volmdlr.Contour3D([e1.copy(), e9.copy(), e5.copy(), e10.copy()]),
                volmdlr.wires.Contour3D([e1.copy(), e10.copy(), e5_switch.copy(), e9_switch.copy()]),
                # volmdlr.Contour3D([e2.copy(), e10.copy(), e6.copy(), e11.copy()]),
                volmdlr.wires.Contour3D([e2.copy(), e11.copy(), e6_switch.copy(), e10_switch.copy()]),
                # volmdlr.Contour3D([e3.copy(), e11.copy(), e7.copy(), e12.copy()]),
                volmdlr.wires.Contour3D([e3.copy(), e12.copy(), e7_switch.copy(), e11_switch.copy()]),
                # volmdlr.Contour3D([e4.copy(), e12.copy(), e8.copy(), e9.copy()])]
                volmdlr.wires.Contour3D([e4.copy(), e9.copy(), e8_switch.copy(), e12_switch.copy()])]

    def shell_faces(self):
        hlx = 0.5*self.frame.u.norm()
        hly = 0.5*self.frame.v.norm()
        hlz = 0.5*self.frame.w.norm()
        frame = self.frame.copy()
        frame.normalize()
        xm_frame = volmdlr.Frame3D(frame.origin-0.5*self.frame.u,
                                   frame.v, frame.w, frame.u)
        xp_frame = volmdlr.Frame3D(frame.origin+0.5*self.frame.u,
                                   frame.v, frame.w, frame.u)
        ym_frame = volmdlr.Frame3D(frame.origin-0.5*self.frame.v,
                                   frame.w, frame.u, frame.v)
        yp_frame = volmdlr.Frame3D(frame.origin+0.5*self.frame.v,
                                   frame.w, frame.u, frame.v)
        zm_frame = volmdlr.Frame3D(frame.origin-0.5*self.frame.w,
                                   frame.u, frame.v, frame.w)
        zp_frame = volmdlr.Frame3D(frame.origin+0.5*self.frame.w,
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

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            new_frame = self.frame.rotation(center, axis, angle, copy=True)
            return Block(new_frame, color=self.color, alpha=self.alpha, name=self.name)
        else:
            self.frame.rotation(center, axis, angle, copy=False)
            self.faces = self.shell_faces()

    def translation(self, offset, copy=True):
        if copy:
            new_frame = self.frame.translation(offset, copy=True)
            return Block(new_frame, color=self.color, alpha=self.alpha, name=self.name)
        else:
            self.frame.translation(offset, copy=False)
            self.faces = self.shell_faces()

    def cut_by_orthogonal_plane(self, plane_3d: volmdlr.faces.Plane3D):
        bb = self.bounding_box
        if plane_3d.frame.w.dot(volmdlr.Vector3D(1, 0, 0)) == 1:
            pass
        elif plane_3d.frame.w.dot(volmdlr.Vector3D(0, 1, 0)) == 1:
            pass
        elif plane_3d.frame.w.dot(volmdlr.Vector3D(0, 0, 1)) == 1:
            pass
        else:
            raise KeyError('plane is not orthogonal either with x, y or z')

        dir1 = plane_3d.frame.u
        dir2 = plane_3d.frame.v
        center2d = volmdlr.Point2D(bb.center.dot(dir1), bb.center.dot(dir2))
        point_min = volmdlr.Point3D(bb.xmin, bb.ymin, bb.zmin)
        point_max = volmdlr.Point3D(bb.xmax, bb.ymax, bb.zmax)
        points = [-center2d + volmdlr.Point2D(point_min.dot(dir1), point_min.dot(dir2)),
                  -center2d + volmdlr.Point2D(point_min.dot(dir1), point_max.dot(dir2)),
                  -center2d + volmdlr.Point2D(point_max.dot(dir1), point_max.dot(dir2)),
                  -center2d + volmdlr.Point2D(point_max.dot(dir1), point_min.dot(dir2))]
        contour_2d = volmdlr.faces.Surface2D(volmdlr.wires.ClosedPolygon2D(points), [])

        return volmdlr.faces.PlaneFace3D(plane_3d, contour_2d)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.basis()
        if side == 'new':
            new_origin = frame.new_coordinates(self.frame.origin)
            new_u = basis.new_coordinates(self.frame.u)
            new_v = basis.new_coordinates(self.frame.v)
            new_w = basis.new_coordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return Block(new_frame, color=self.color, alpha=self.alpha, name=self.name)
            else:
                self.frame = new_frame
                self.faces = self.shell_faces()

        if side == 'old':
            new_origin = frame.old_coordinates(self.frame.origin)
            new_u = basis.old_coordinates(self.frame.u)
            new_v = basis.old_coordinates(self.frame.v)
            new_w = basis.old_coordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return Block(new_frame, color=self.color, alpha=self.alpha, name=self.name)
            else:
                self.frame = new_frame
                self.faces = self.shell_faces()

    def copy(self):
        new_origin = self.frame.origin.copy()
        new_u = self.frame.u.copy()
        new_v = self.frame.v.copy()
        new_w = self.frame.w.copy()
        new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
        return Block(new_frame, color=self.color, alpha=self.alpha, name=self.name)

    def plot_data(self, x3D, y3D, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        lines = []
        for edge3D in self.Edges():
            lines.append(edge3D.plot_data(x3D, y3D, marker, color, stroke_width,
                         dash, opacity, arrow))

        return lines

    def plot2D(self, x3D, y3D, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = None

        for edge3D in self.Edges():
#            edge2D = edge3D.PlaneProjection2D()
            edge3D.plot2D(x3D, y3D, ax)

        return fig, ax


class Cone(volmdlr.core.Primitive3D):
    def __init__(self, position, axis, radius, length, name=''):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.position = position
        axis.normalize()
        self.axis = axis
        self.radius = radius
        self.length = length
        self.bounding_box = self._bounding_box()

    def _bounding_box(self):
        """
        A is the point at the basis
        B is the top
        """
        pointA = self.position - self.length/2 * self.axis
        pointB = self.position + self.length/2 * self.axis

        dx2 = (pointA[0]-pointB[0])**2
        dy2 = (pointA[1]-pointB[1])**2
        dz2 = (pointA[2]-pointB[2])**2

        kx = ((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5
        ky = ((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5
        kz = ((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5

        x_bound = (pointA[0] - kx * self.radius, pointA[0] + kx * self.radius, pointB[0])
        xmin = min(x_bound)
        xmax = max(x_bound)

        y_bound = (pointA[1] - ky * self.radius, pointA[1] + ky * self.radius, pointB[1])
        ymin = min(y_bound)
        ymax = max(y_bound)

        z_bound = (pointA[2] - kz * self.radius, pointA[2] + kz * self.radius, pointB[2])
        zmin = min(z_bound)
        zmax = max(z_bound)

        return volmdlr.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def volume(self):
        return self.length * math.pi * self.radius**2 / 3

    def babylon_script(self):
        new_axis = volmdlr.Vector3D((self.axis[0], self.axis[1], self.axis[2]))
        normal_vector1 = new_axis.RandomUnitnormalVector()
        normal_vector2 = new_axis.cross(normal_vector1)
        x, y, z = self.position
        s = 'var cone = BABYLON.MeshBuilder.CreateCylinder("cone", {{diameterTop:0, diameterBottom:{}, height: {}, tessellation: 100}}, scene);\n'.format(2*self.radius, self.length)
        s += 'cone.position = new BABYLON.Vector3({},{},{});\n;'.format(x,y,z)
        s += 'var axis1 = new BABYLON.Vector3({},{},{});\n'.format(new_axis[0], new_axis[1], new_axis[2])
        s += 'var axis2 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector1[0], normal_vector1[1], normal_vector1[2])
        s += 'var axis3 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector2[0], normal_vector2[1], normal_vector2[2])
        s += 'cone.rotation = BABYLON.Vector3.rotationFromAxis(axis3, axis1, axis2);\n'
        return s


class ExtrudedProfile(volmdlr.faces.ClosedShell3D):
    """
    
    """
    _non_serializable_attributes  = ['faces', 'inner_contours3d', 'outer_contour3d']
    def __init__(self, plane_origin:volmdlr.Point3D,
                 x:volmdlr.Vector3D, y:volmdlr.Vector3D,
                 outer_contour2d:volmdlr.wires.Contour2D,
                 inner_contours2d:List[volmdlr.wires.Contour2D],
                 extrusion_vector:volmdlr.Vector3D,
                 color:Tuple[float, float, float]=None, alpha=1., name=''):
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



    def shell_faces(self):
        lower_plane = volmdlr.faces.Plane3D.from_plane_vectors(self.plane_origin,
                                                               self.x,
                                                               self.y)
        lower_face = volmdlr.faces.PlaneFace3D(lower_plane,
                                               volmdlr.faces.Surface2D(
                                                    self.outer_contour2d,
                                                    self.inner_contours2d)
                                               )

        upper_face = lower_face.translation(self.extrusion_vector)
        lateral_faces = []
        for p in self.outer_contour3d.primitives:
            lateral_faces.extend(p.extrusion(self.extrusion_vector))
         
        for inner_contour in self.inner_contours3d:
            for p in inner_contour.primitives:
                lateral_faces.extend(p.extrusion(self.extrusion_vector))

        return [lower_face]+[upper_face]+lateral_faces

    def plot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        for contour in [self.outer_contour2d]+self.inner_contours2d:
            for primitive in contour.primitives:
                primitive.plot(ax)
        ax.margins(0.1)
        return ax

    def FreeCADExport(self, ip):
        name='primitive'+str(ip)
        s = 'Wo = []\n'
        s += 'Eo = []\n'
        for ip, primitive in enumerate(self.outer_contour3d.primitives):
            s += primitive.FreeCADExport('L{}'.format(ip))
            s += 'Eo.append(Part.Edge(L{}))\n'.format(ip)
        s += 'Wo.append(Part.Wire(Eo[:]))\n'
        s += 'Fo = Part.Face(Wo)\n'

        s += 'Fi = []\n'
        s += 'W = []\n'
        for ic,contour in enumerate(self.inner_contours3d):
            s+='E = []\n'
            for ip, primitive in enumerate(contour.primitives):
                s += primitive.FreeCADExport('L{}_{}'.format(ic, ip))
                s += 'E.append(Part.Edge(L{}_{}))\n'.format(ic, ip)
            s += 'Wi = Part.Wire(E[:])\n'
            s += 'Fi.append(Part.Face(Wi))\n'

        if len(self.inner_contours3d) != 0:
            s += 'Fo = Fo.cut(Fi)\n'
        e1, e2, e3 = round(1000*self.extrusion_vector, 6)

        s+='{} = Fo.extrude(fc.Vector({}, {}, {}))\n'.format(name,e1,e2,e3)
        return s

    def area(self):
        areas = self.outer_contour2d.area()
        areas -= sum([c.area() for c in self.inner_contours2d])
        # sic=list(npy.argsort(areas))[::-1]# sorted indices of contours
        # area=areas[sic[0]]

        # for i in sic[1:]:
        #     area-=self.contours2D[i].Area()
        return areas

    def volume(self):
        z = self.x.cross(self.y)
        z.normalize()
        coeff = npy.dot(self.extrusion_vector, z)
        return self.Area()*coeff


    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.basis()
        if side == 'old':
            extrusion_vector = basis.old_coordinates(self.extrusion_vector)
            x = basis.old_coordinates(self.x)
            y = basis.old_coordinates(self.y)
        elif side == 'new':
            extrusion_vector = basis.new_coordinates(self.extrusion_vector)
            x = basis.new_coordinates(self.x)
            y = basis.new_coordinates(self.y)
        else:
            raise ValueError('side must be either old or new')

        if copy:
            return ExtrudedProfile(self.plane_origin.frame_mapping(frame, side, copy),
                                   x,
                                   y,
                                   self.outer_contour2d,
                                   self.inner_contours2d,
                                   extrusion_vector)
        else:
            self.plane_origin.frame_mapping(frame, side, copy)
            self.__init__(self.plane_origin,
                          x,
                          y,
                          self.outer_contour2d,
                          self.inner_contours2d,
                          extrusion_vector)

    def translation(self, offset: volmdlr.Vector3D, copy=True):
        if copy:

            return self.__class__(plane_origin=self.plane_origin.translation(offset, copy=True),
                                  x=self.x, y=self.y,
                                  outer_contour2d=self.outer_contour2d,
                                  inner_contours2d=self.inner_contours2d,
                                  extrusion_vector=self.extrusion_vector,
                                  color=self.color, alpha=self.alpha)
        else:
            self.plane_origin.translation(offset, copy=False)
            self.axis_point.translation(offset, copy=False)

    def rotation(self, center, axis, angle, copy=True):
        if copy:

            return self.__class__(plane_origin=self.plane_origin.rotation(center, axis, angle, copy=True),
                                  x=self.x.rotation(volmdlr.O3D, axis, angle, copy=True),
                                  y=self.y.rotation(volmdlr.O3D, axis, angle, copy=False),
                                  outer_contour2d=self.outer_contour2d,
                                  inner_contours2d=self.inner_contours2d,
                                  extrusion_vector=self.extrusion_vector.rotation(volmdlr.O3D, axis, angle, copy=True),
                                  color=self.color, alpha=self.alpha)
        else:
            self.plane_origin.rotation(center, axis, angle, copy=False)
            self.x.rotation(volmdlr.O3D, axis, angle, copy=False)
            self.y.rotation(volmdlr.O3D, axis, angle, copy=False)
            self.extrusion_vector.rotation(volmdlr.O3D, axis, angle, copy=False)

class RevolvedProfile(volmdlr.faces.ClosedShell3D):
    """

    """
    _non_serializable_attributes  = ['faces', 'contour3D']

    def __init__(self, plane_origin:volmdlr.Point3D,
                 x:volmdlr.Vector3D, y:volmdlr.Vector3D,
                 contour2d:volmdlr.wires.Contour2D,
                 axis_point:volmdlr.Point3D, axis:volmdlr.Vector3D,
                 angle:float=2*math.pi, *,
                 color:Tuple[float, float, float]=None, alpha:float=1,
                 name:str=''):

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
            face1 = volmdlr.faces.PlaneFace3D(plane1,
                                              volmdlr.faces.Surface2D(self.contour2d,
                                                                      []))
            face2 = face1.rotation(self.axis_point, self.axis, self.angle)
            faces.append(face1)
            faces.append(face2)
            

        return faces

    def plot(self, ax=None):
        # if ax is None:
        #     fig, ax = plt.subplots()
        # for contour in self.contours3d:
        ax = self.contour3d.plot(ax)

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive'+str(ip)
        s = 'W=[]\n'
#        for ic, contour in enumerate(self.contours3D):
        s += 'L=[]\n'
        for ibp, basis_primitive in enumerate(self.contour3D.edges):
            s += basis_primitive.FreeCADExport('L{}_{}'.format(1, ibp), 8)
            s += 'L.append(L{}_{})\n'.format(1,ibp)
        s += 'S = Part.Shape(L)\n'
        s += 'W.append(Part.Wire(S.Edges))\n'
        s += 'F=Part.Face(W)\n'
        a1, a2, a3 = self.axis.vector
        ap1, ap2, ap3 = self.axis_point.vector
        ap1 = round(ap1*1000, ndigits)
        ap2 = round(ap2*1000, ndigits)
        ap3 = round(ap3*1000, ndigits)
        angle = self.angle/math.pi*180
        s += '{} = F.revolve(fc.Vector({},{},{}), fc.Vector({},{},{}),{})\n'.format(name, ap1,ap2,ap3,a1,a2,a3,angle)
        return s

    def volume(self):
        p1 = self.axis_point.PlaneProjection3D(self.plane_origin,self.x,self.y)
        p1_2D = p1.To2D(self.axis_point,self.x,self.y)
        p2_3D = self.axis_point+volmdlr.Point3D(self.axis.vector)
        p2_2D = p2_3D.To2D(self.plane_origin,self.x,self.y)
        axis_2D = volmdlr.Line2D(p1_2D,p2_2D)
        com = self.contour2D.CenterOfMass()
        if com is not False:
            rg = axis_2D.point_distance(com)
            return self.angle*rg*self.contour2D.Area()
        else:
            return 0

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'old':
            axis = basis.old_coordinates(self.axis)
            x = basis.old_coordinates(self.x)
            y = basis.old_coordinates(self.y)
        elif side == 'new':
            axis = basis.new_coordinates(self.axis)
            x = basis.new_coordinates(self.x)
            y = basis.new_coordinates(self.y)
        else:
            raise ValueError('side must be either old or new')

        if copy:

            return RevolvedProfile(self.plane_origin.frame_mapping(frame, side, copy),
                                   x,
                                   y,
                                   self.contour2D,
                                   self.axis_point.frame_mapping(frame, side, copy),
                                   axis,
                                   self.angle)
        else:
            self.__init__(self.plane_origin.frame_mapping(frame, side, copy),
                          x,
                          y,
                          self.contour2D,
                          self.axis_point.frame_mapping(frame, side, copy),
                          axis,
                          self.angle)

    def translation(self, offset: volmdlr.Vector3D, copy=True):
        if copy:

            return self.__class__(plane_origin=self.plane_origin.translation(offset, copy=True),
                                  x=self.x, y=self.y, contour2d=self.contour2d,
                                  axis_point=self.axis_point.translation(offset, copy=True),
                                  axis=self.axis,
                                  angle=self.angle,
                                  color=self.color, alpha=self.alpha)
        else:
            self.plane_origin.translation(offset, copy=False)
            self.axis_point.translation(offset, copy=False)

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            return self.__class__(plane_origin=self.plane_origin.rotation(center, axis, angle, copy=True),
                                  x=self.x.rotation(center=volmdlr.O3D, axis=axis, angle=angle, copy=True),
                                  y=self.y.rotation(center=volmdlr.O3D, axis=axis, angle=angle, copy=True),
                                  contour2d=self.contour2d,
                                  axis_point=self.axis_point.rotation(center, axis, angle, copy=True),
                                  axis=self.axis.rotation(center=volmdlr.O3D,
                                                          axis=axis,
                                                          angle=angle, copy=True),
                                  angle=self.angle,
                                  color=self.color, alpha=self.alpha)
        else:
            self.plane_origin.rotation(center, axis, angle, copy=False)
            self.x.rotation(center=volmdlr.O3D, axis=axis, angle=angle,
                            copy=False)
            self.y.rotation(center=volmdlr.O3D, axis=axis, angle=angle,
                            copy=False)
            self.axis_point.rotation(center, axis, angle, copy=False)
            self.axis.rotation(volmdlr.O3D, axis, angle, copy=False)

class Cylinder(RevolvedProfile):
    """
    Creates a full cylinder with the position, the axis of revolution, the radius and the length.
    """
    def __init__(self, position:volmdlr.Point3D, axis:volmdlr.Vector3D,
                 radius:float, length:float,
                 color:Tuple[float, float, float]=None, alpha:float=1.,
                 name:str=''):

        self.position = position
        axis.normalize()
        self.axis = axis
        self.radius = radius
        self.length = length
        self.bounding_box = self._bounding_box()

        # Revolved Profile
        p1 = volmdlr.Point2D(-0.5*self.length, 0.)
        p2 = volmdlr.Point2D(0.5*self.length, 0.)
        p3 = volmdlr.Point2D(0.5*self.length, self.radius)
        p4 = volmdlr.Point2D(-0.5*self.length, self.radius)
        l1 = volmdlr.edges.LineSegment2D(p1, p2)
        l2 = volmdlr.edges.LineSegment2D(p2, p3)
        l3 = volmdlr.edges.LineSegment2D(p3, p4)
        l4 = volmdlr.edges.LineSegment2D(p4, p1)
        contour = volmdlr.wires.Contour2D([l1, l2, l3, l4])
        y = axis.random_unit_normal_vector()
        RevolvedProfile.__init__(self, position, axis, y, contour, position, axis,
                                 color=color, alpha=alpha, name=name)


    def _bounding_box(self):
        
        if hasattr(self, 'radius'):
            radius = self.radius
        elif hasattr(self, 'outer_radius'):
            radius = self.outer_radius

            
        pointA = self.position - self.length/2 * self.axis
        pointB = self.position + self.length/2 * self.axis

        dx2 = (pointA[0]-pointB[0])**2
        dy2 = (pointA[1]-pointB[1])**2
        dz2 = (pointA[2]-pointB[2])**2

        kx = ((dy2 + dz2) / (dx2 + dy2 + dz2))**0.5
        ky = ((dx2 + dz2) / (dx2 + dy2 + dz2))**0.5
        kz = ((dx2 + dy2) / (dx2 + dy2 + dz2))**0.5

        if pointA[0] > pointB[0]:
            pointA, pointB = pointB, pointA
        xmin = pointA[0] - kx * radius
        xmax = pointB[0] + kx * radius

        if pointA[1] > pointB[1]:
            pointA, pointB = pointB, pointA
        ymin = pointA[1] - ky * radius
        ymax = pointB[1] + ky * radius

        if pointA[2] > pointB[2]:
            pointA, pointB = pointB, pointA
        zmin = pointA[2] - kz * radius
        zmax = pointB[2] + kz * radius

        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def volume(self):
        return self.length * math.pi * self.radius**2

    def FreeCADExport(self, ip):
        if self.radius > 0:
            name = 'primitive'+str(ip)
            e = str(1000*self.length)
            r = str(1000*self.radius)
            position = 1000*(self.position - self.axis*self.length/2.)
            x, y, z = position
            x = str(x)
            y = str(y)
            z = str(z)

            ax, ay, az = self.axis
            ax = str(ax)
            ay = str(ay)
            az = str(az)
            return name+'=Part.makeCylinder('+r+','+e+',fc.Vector('+x+','+y+','+z+'),fc.Vector('+ax+','+ay+','+az+'),360)\n'
        else:
            return ''

    def babylon_script(self, name='primitive_mesh'):
        normal_vector1 = self.axis.RandomUnitnormalVector()
#        normal_vector2 = new_axis.cross(normal_vector1)
#        x, y, z = self.position
#        s='var {} = BABYLON.Mesh.CreateCylinder("{}", {}, {}, {}, 30, 1, scene,false, BABYLON.Mesh.DEFAULTSIDE);'.format(name, self.name,self.length,2*self.radius,2*self.radius)
#        s += '{}.position = new BABYLON.Vector3({},{},{});\n;'.format(name, x,y,z)
#        s += 'var axis1 = new BABYLON.Vector3({},{},{});\n'.format(new_axis[0], new_axis[1], new_axis[2])
#        s += 'var axis2 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector1[0], normal_vector1[1], normal_vector1[2])
#        s += 'var axis3 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector2[0], normal_vector2[1], normal_vector2[2])
#        s += '{}.rotation = BABYLON.Vector3.rotationFromAxis(axis3, axis1, axis2);\n'.format(name)
        p1 = volmdlr.Point2D((-0.5*self.length, self.radius))
        p2 = volmdlr.Point2D((0.5*self.length, self.radius))
        p3 = volmdlr.Point2D((0.5*self.length, 0.))
        p4 = volmdlr.Point2D((-0.5*self.length, 0.))
        l1 = volmdlr.LineSegment2D(p1, p2)
        l2 = volmdlr.LineSegment2D(p2, p3)
        l3 = volmdlr.LineSegment2D(p3, p4)
        l4 = volmdlr.LineSegment2D(p4, p1)
        extruded_profile = RevolvedProfile(self.position, self.axis, normal_vector1,
                                           volmdlr.Contour2D([l1, l2, l3, l4]), self.position, self.axis, name=self.name)
        return extruded_profile.babylon_script(name=name)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'old':
            axis = basis.old_coordinates(self.axis)
        elif side == 'new':
            axis = basis.new_coordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')

        if copy:
            return Cylinder(self.position.frame_mapping(frame, side, copy),
                            axis,
                            self.radius, self.length, color=self.color,
                            alpha=self.alpha)
        else:
            self.position.frame_mapping(frame, side, copy)
            self.axis = axis
            Cylinder.__init__(self, self.position, self.axis, self.radius, 
                              self.length, color=self.color, alpha=self.alpha)


    def translation(self, offset: volmdlr.Vector3D, copy=True):
        if copy:

            return self.__class__(position=self.position.translation(offset, copy=True),
                                  axis=self.axis,
                                  length=self.length, radius=self.radius)
        else:
            self.position.translation(offset, copy=False)

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            return self.__class__(position=self.position.rotation(center, axis,
                                                                  angle, copy=True),
                                  axis=self.axis.rotation(volmdlr.O3D, axis,
                                                          angle,copy=True),
                                  length=self.length, radius=self.radius)
        else:
            self.position.rotation(center, axis, angle, copy=False)
            self.axis.rotation(volmdlr.O3D, axis, angle, copy=False)

    def copy(self):
        new_position = self.position.copy()
        new_axis = self.axis.copy()
        return Cylinder(new_position, new_axis, self.radius, self.length,
                        color=self.color, alpha=self.alpha, name=self.name)
        

class HollowCylinder(Cylinder):
    def __init__(self, position:volmdlr.Point3D, axis:volmdlr.Vector3D,
                 inner_radius:float, outer_radius:float, length:float,
                 color:Tuple[float, float, float]=None, alpha:float=1,
                 name:str=''):
        volmdlr.core.Primitive3D.__init__(self, name=name)
        self.position = position
        axis.normalize()
        self.axis = axis
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.length = length
        
        # Revolved Profile
        p1 = volmdlr.Point2D(-0.5*self.length, self.inner_radius)
        p2 = volmdlr.Point2D(0.5*self.length, self.inner_radius)
        p3 = volmdlr.Point2D(0.5*self.length, self.outer_radius)
        p4 = volmdlr.Point2D(-0.5*self.length, self.outer_radius)
        l1 = volmdlr.edges.LineSegment2D(p1, p2)
        l2 = volmdlr.edges.LineSegment2D(p2, p3)
        l3 = volmdlr.edges.LineSegment2D(p3, p4)
        l4 = volmdlr.edges.LineSegment2D(p4, p1)
        contour = volmdlr.wires.Contour2D([l1, l2, l3, l4])
        y = axis.random_unit_normal_vector()
        # contour.plot()
        RevolvedProfile.__init__(self, position, axis, y, contour, position, axis,
                                 color=color, alpha=alpha, name=name)

        

    def volume(self):
        return self.length * math.pi* (self.outer_radius**2 - self.inner_radius**2)


    def FreeCADExport(self, ip):
        if self.outer_radius > 0.:
            name = 'primitive'+str(ip)
            re = round(1000*self.outer_radius, 6)
            ri = round(1000*self.inner_radius, 6)
            x, y, z = round((1000*(self.position - self.axis*self.length/2)), 6)
            ax, ay, az = npy.round(self.axis.vector, 6)

            s='C2 = Part.makeCircle({}, fc.Vector({}, {}, {}),fc.Vector({}, {}, {}))\n'.format(re, x, y, z, ax, ay, az)
            s+='W2 = Part.Wire(C2.Edges)\n'
            s+='F2 = Part.Face(W2)\n'

            if self.inner_radius!=0.:
                s+='C1 = Part.makeCircle({}, fc.Vector({}, {}, {}),fc.Vector({}, {}, {}))\n'.format(ri, x, y, z, ax, ay, az)
                s+='W1 = Part.Wire(C1.Edges)\n'
                s+='F1 = Part.Face(W1)\n'
                s+='F2 = F2.cut(F1)\n'

            vx, vy, vz = round(self.axis*self.length*1000, 6)

            s += '{} = F2.extrude(fc.Vector({}, {}, {}))\n'.format(name, vx, vy, vz)
            return s

        else:
            return ''


    def babylon_script(self, name='primitive_mesh'):
        normal_vector1 = self.axis.RandomUnitnormalVector()
        p1 = volmdlr.Point2D((-0.5*self.length, self.outer_radius))
        p2 = volmdlr.Point2D((0.5*self.length, self.outer_radius))
        p3 = volmdlr.Point2D((0.5*self.length, self.inner_radius))
        p4 = volmdlr.Point2D((-0.5*self.length, self.inner_radius))
        l1 = volmdlr.LineSegment2D(p1, p2)
        l2 = volmdlr.LineSegment2D(p2, p3)
        l3 = volmdlr.LineSegment2D(p3, p4)
        l4 = volmdlr.LineSegment2D(p4, p1)
        extruded_profile = RevolvedProfile(self.position, self.axis, normal_vector1,
                                           volmdlr.Contour2D([l1, l2, l3, l4]),
                                           self.position, self.axis, name=self.name)
        return extruded_profile.babylon_script(name=name)


    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'old':
            axis = basis.old_coordinates(self.axis)
        elif side == 'new':
            axis = basis.new_coordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')

        if copy:
            return HollowCylinder(position=self.position.frame_mapping(frame, side, copy),
                                  axis=axis,
                                  inner_radius=self.inner_radius, 
                                  outer_radius=self.outer_radius,
                                  length=self.length)
        else:
            self.position.frame_mapping(frame, side, copy)
            self.axis = axis

    def translation(self, offset: volmdlr.Vector3D, copy=True):
        if copy:

            return self.__class__(position=self.position.translation(offset, copy=True),
                                  axis=self.axis,
                                  length=self.length,
                                  inner_radius=self.inner_radius,
                                  outer_radius=self.outer_radius)
        else:
            self.position.translation(offset, copy=False)

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            return self.__class__(position=self.position.rotation(center, axis, angle,
                                                                  copy=True),
                                  axis=self.axis.rotation(volmdlr.O3D, axis,
                                                          angle,copy=True),
                                  length=self.length,
                                  inner_radius=self.inner_radius,
                                  outer_radius=self.outer_radius)
        else:
            self.position.rotation(center, axis, angle, copy=False)
            self.axis.rotation(volmdlr.O3D, axis, angle, copy=False)

class Sweep(volmdlr.faces.ClosedShell3D):
    """
    Sweep a 2D contour along a Wire3D
    """

    def __init__(self, contour2d:List[volmdlr.wires.Contour2D],
                 wire3d:volmdlr.wires.Wire3D, *,
                 color:Tuple[float, float, float]=None, alpha:float=1,
                 name:str=''):
        self.contour2d = contour2d
        self.wire3d = wire3d
        self.frames = []
        
        faces = self.shell_faces()
        volmdlr.faces.ClosedShell3D.__init__(self, faces, color=color, alpha=alpha, name=name)


    def shell_faces(self):
        """
        For now it does not take into account rotation of sections
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
            volmdlr.Frame3D(self.wire3d.primitives[-1].point_at_abscissa(l_last_primitive),
                            u, v, w))

        faces = [volmdlr.faces.PlaneFace3D(start_plane,
                                           volmdlr.faces.Surface2D(self.contour2d, [])),
                 volmdlr.faces.PlaneFace3D(end_plane,
                                           volmdlr.faces.Surface2D(self.contour2d,
                                                             [])),
                 ]

        for wire_primitive in self.wire3d.primitives :
            # tangent, normal = wire_primitive.frenet(0.)
            tangent = wire_primitive.unit_direction_vector(0.)
            normal = wire_primitive.unit_normal_vector(0.)

            if normal is None:
                normal = tangent.deterministic_unit_normal_vector()
            n2 = tangent.cross(normal)
            contour3d = self.contour2d.to_3d(wire_primitive.start,
                                            normal,
                                            n2)

            if wire_primitive.__class__ is volmdlr.edges.LineSegment3D:
                for contour_primitive in contour3d.primitives:
                    faces.extend(contour_primitive.extrusion(
                        wire_primitive.direction_vector()))
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
            else:
                raise NotImplementedError(
                    'Unimplemented primitive for sweep: {}'\
                        .format(wire_primitive.__class__.__name__)
                )

        return faces


    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_wire = self.wire3d.frame_mapping(frame, side, copy)
            return Sweep(self.contour2d, new_wire, color=self.color, alpha=self.alpha, name=self.name)
        else:
            self.wire3d.frame_mapping(frame, side, copy=False)
            for face in self.faces :
                face.frame_mapping(frame, side, copy=False)
            

    def copy(self):
        new_contour2d = self.contour2d.copy()
        new_wire3d = self.wire3d.copy()
        return Sweep(new_contour2d, new_wire3d, color=self.color, alpha=self.alpha, name=self.name)
    
    
# class Sphere(volmdlr.Primitive3D):
class Sphere(RevolvedProfile):
    def __init__(self, center, radius,
                 color:Tuple[float, float, float]=None, alpha=1., name=''):
        volmdlr.core.Primitive3D.__init__(self, name=name)
        self.center = center
        self.radius = radius
        self.position = center
        
        # Revolved Profile for complete sphere
        s = volmdlr.Point2D(-self.radius, 0.01*self.radius)
        i = volmdlr.Point2D(0, 1.01*self.radius)
        e = volmdlr.Point2D(self.radius, 0.01*self.radius) #Not coherent but it works at first, to change !!
        
        # s = volmdlr.Point2D((-self.radius, 0))
        # i = volmdlr.Point2D(((math.sqrt(2)/2)*self.radius,(math.sqrt(2)/2)*self.radius))
        # e = volmdlr.Point2D(((-math.sqrt(2)/2)*self.radius,(-math.sqrt(2)/2)*self.radius)) 
        
        contour = volmdlr.wires.Contour2D([
            volmdlr.edges.Arc2D(s, i , e), volmdlr.edges.LineSegment2D(s, e)])
        # fig, ax = plt.subplots()
        # c.plot(ax=ax)
        
        # contour = volmdlr.Contour2D([c])
        axis = volmdlr.X3D
        y = axis.random_unit_normal_vector()
        RevolvedProfile.__init__(self, center, axis, y, contour, center, axis,
                                 color=color, alpha=alpha, name=name)

    def volume(self):
        return 4/3*math.pi*self.radius**3

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive'+str(ip)
        r = 1000*self.radius
        x, y, z = round(1000*self.center, ndigits)
        return '{} = Part.makeSphere({}, fc.Vector({}, {}, {}))\n'.format(name,r,x,y,z)

    def babylon_script(self, name='primitive_mesh'):
        p1 = volmdlr.Point2D((-self.radius, 0))
        p2 = volmdlr.Point2D((0, self.radius))
        p3 = volmdlr.Point2D((self.radius, 0))
        line = volmdlr.LineSegment2D(p1, p3)
        arc = volmdlr.Arc2D(p1, p2, p3)
        extruded_profile = RevolvedProfile(self.position, volmdlr.X3D, volmdlr.Y3D,
                                           volmdlr.Contour2D([line, arc]), self.position, volmdlr.X3D, name=self.name)
        return extruded_profile.babylon_script(name=name)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            return Sphere(self.center.frame_mapping(frame, side, copy),
                          self.radius)
        else:
            self.center.frame_mapping(frame, side, copy)




class Measure3D(volmdlr.edges.Line3D):
    def __init__(self, point1, point2, color=(1, 0, 0)):
        self.point1, self.point2 = point1, point2
        self.color = color
        self.distance = (point1 - point2).norm()
        self.bounding_box = self._bounding_box()

    # !!! no eq defined!
    def __hash__(self):
        return sum([hash(p) for p in self.points])

    def babylon_script(self):
        s = 'var myPoints = [];\n'
        s += 'var point1 = new BABYLON.Vector3({},{},{});\n'.format(
            self.points[0][0], self.points[0][1], self.points[0][2])
        s += 'myPoints.push(point1);\n'
        s += 'var point2 = new BABYLON.Vector3({},{},{});\n'.format(
            self.points[1][0], self.points[1][1], self.points[1][2])
        s += 'myPoints.push(point2);\n'
        s += 'var line = BABYLON.MeshBuilder.CreateLines("lines", {points: myPoints}, scene);\n'
        s += 'line.color = new BABYLON.Color3({}, {}, {});\n'.format(
            self.color[0], self.color[1], self.color[2])
        return s


class BSplineExtrusion(volmdlr.core.Primitive3D):

    def __init__(self, obj, vectorextru, name=''):
        self.obj = obj
        vectorextru.normalize()
        self.vectorextru = vectorextru
        if obj.__class__ is volmdlr.edges.Ellipse3D:
            self.points = obj.tessel_points
        else:
            self.points = obj.points

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]
        if object_dict[arguments[1]].__class__ is volmdlr.edges.Ellipse3D:
            ell = object_dict[arguments[1]]
            vectextru = -object_dict[arguments[2]]
            return cls(ell, vectextru, name)

        elif object_dict[arguments[1]].__class__ is volmdlr.edges.BSplineCurve3D:
            bsplinecurve = object_dict[arguments[1]]
            vectextru = object_dict[arguments[2]]
            return cls(bsplinecurve, vectextru, name)
        else:
            raise NotImplementedError  ## a adapter pour les bpsline
