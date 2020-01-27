#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common primitives 3D
"""

import math

import numpy as npy
npy.seterr(divide='raise')

import volmdlr
from volmdlr.primitives import RoundedLineSegments
from typing import Tuple


import matplotlib.pyplot as plt

class OpenedRoundedLineSegments3D(volmdlr.Wire3D, RoundedLineSegments):
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    
    def __init__(self, points, radius, adapt_radius=False, name=''):
        primitives = RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.LineSegment3D,
                                                  volmdlr.Arc3D,
                                                  closed=False,
                                                  adapt_radius=adapt_radius,
                                                  name='')

        volmdlr.Wire3D.__init__(self, primitives, name)

    def ArcFeatures(self, ipoint):
        radius = self.radius[ipoint]
        if self.closed:
            if ipoint == 0:
                pt1 = self.points[-1]
            else:
                pt1 = self.points[ipoint -1]
            pti = self.points[ipoint]
            if ipoint < self.npoints-1:
                pt2 = self.points[ipoint+1]
            else:
                pt2 = self.points[0]
        else:
            pt1 = self.points[ipoint - 1]
            pti = self.points[ipoint]
            pt2 = self.points[ipoint + 1]

        dist1 = (pt1 - pti).Norm()
        dist2 = (pt2 - pti).Norm()
        dist3 = (pt1 - pt2).Norm()
        alpha = math.acos(-(dist3**2-dist1**2-dist2**2)/(2*dist1*dist2))/2.
        dist = radius/math.tan(alpha)

        u1 = (pt1 - pti) / dist1
        u2 = (pt2 - pti) / dist2

        p3 = pti + u1*dist
        p4 = pti + u2*dist

        n = u1.Cross(u2)
        n /= n.Norm()
        v1 = u1.Cross(n)
        v2 = u2.Cross(n)

        l1 = volmdlr.Line3D(p3, p3+v1)
        l2 = volmdlr.Line3D(p4, p4+v2)
        c, _ = l1.MinimumDistancePoints(l2)

        u3 = u1 + u2# mean of v1 and v2
        u3 /= u3.Norm()

        interior = c - u3 * radius
        return p3, interior, p4, dist, alpha


    def Rotation(self, center, angle, copy=True):
        if copy:
            return self.__class__([p.Rotation(center, angle, copy=True)\
                                          for p in self.points],
                                         self.radius, self.closed, self.name)
        else:
            self.__init__([p.Rotation(center, angle, copy=True)\
                           for p in self.points],
                          self.radius, self.closed, self.name)

    def Translation(self, offset, copy=True):
        if copy:
            return self.__class__([p.Translation(offset, copy=True)\
                                          for p in self.points],
                                         self.radius, self.closed, self.name)
        else:
            self.__init__([p.Translation(offset, copy=True)\
                           for p in self.points],
                          self.radius, self.closed, self.name)


class ClosedRoundedLineSegments3D(volmdlr.Contour3D, OpenedRoundedLineSegments3D):
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True
    
    def __init__(self, points, radius, adapt_radius=False, name=''):
        primitives = RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.LineSegment3D,
                                                  volmdlr.Arc3D,
                                                  closed=True,
                                                  adapt_radius=adapt_radius,
                                                  name='')

        volmdlr.Contour3D.__init__(self, primitives, name=name)

class Sphere(volmdlr.Primitive3D):
    def __init__(self,center, radius, name=''):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.center = center
        self.radius = radius
        self.position = center

    def Volume(self):
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

class Block(volmdlr.Shell3D):
    _standalone_in_db = True
    _generic_eq = True
    _non_serializable_attributes  = ['size']
    _non_eq_attributes = ['name', 'color', 'size', 'bounding_box', 'faces', 'contours',
                          'plane', 'points', 'polygon2D']
    _non_hash_attributes = []

    """
    Creates a block
    :param frame: a frame 3D. The origin of the frame is the center of the block,
     the 3 vectors are defining the edges. The frame has not to be orthogonal
    """
    def __init__(self, frame:volmdlr.Frame3D, name:str='', color:Tuple[float, float, float]=None):
        self.frame = frame
        self.size = (self.frame.u.Norm(), self.frame.v.Norm(), self.frame.w.Norm())

        faces = self.shell_faces()
        volmdlr.Shell3D.__init__(self, faces, name, color)

    # def __hash__(self):
    #     return hash(self.frame)

    def Vertices(self):
        return [self.frame.origin - 0.5*self.frame.u - 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin - 0.5*self.frame.u + 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u + 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u - 0.5*self.frame.v - 0.5*self.frame.w,
                self.frame.origin - 0.5*self.frame.u - 0.5*self.frame.v + 0.5*self.frame.w,
                self.frame.origin - 0.5*self.frame.u + 0.5*self.frame.v + 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u + 0.5*self.frame.v + 0.5*self.frame.w,
                self.frame.origin + 0.5*self.frame.u - 0.5*self.frame.v + 0.5*self.frame.w]

    def Edges(self):
        p1, p2, p3, p4, p5, p6, p7, p8 = self.Vertices()
        return [volmdlr.LineSegment3D(p1.copy(), p2.copy()),
                volmdlr.LineSegment3D(p2.copy(), p3.copy()),
                volmdlr.LineSegment3D(p3.copy(), p4.copy()),
                volmdlr.LineSegment3D(p4.copy(), p1.copy()),
                volmdlr.LineSegment3D(p5.copy(), p6.copy()),
                volmdlr.LineSegment3D(p6.copy(), p7.copy()),
                volmdlr.LineSegment3D(p7.copy(), p8.copy()),
                volmdlr.LineSegment3D(p8.copy(), p5.copy()),
                volmdlr.LineSegment3D(p1.copy(), p5.copy()),
                volmdlr.LineSegment3D(p2.copy(), p6.copy()),
                volmdlr.LineSegment3D(p3.copy(), p7.copy()),
                volmdlr.LineSegment3D(p4.copy(), p8.copy())]

    def face_contours(self):
        e1, e2, e3, e4, e5, e6, e7, e8, e9, e10, e11, e12 = self.Edges()
        return [volmdlr.Contour3D([e1.copy(), e2.copy(), e3.copy(), e4.copy()]),
                volmdlr.Contour3D([e5.copy(), e6.copy(), e7.copy(), e8.copy()]),
                volmdlr.Contour3D([e1.copy(), e9.copy(), e5.copy(), e10.copy()]),
                volmdlr.Contour3D([e2.copy(), e10.copy(), e6.copy(), e11.copy()]),
                volmdlr.Contour3D([e3.copy(), e11.copy(), e7.copy(), e12.copy()]),
                volmdlr.Contour3D([e4.copy(), e12.copy(), e8.copy(), e9.copy()])]

    def shell_faces(self):
        c1, c2, c3, c4, c5, c6 = self.face_contours()
        return [volmdlr.Face3D([c1]),
                volmdlr.Face3D([c2]),
                volmdlr.Face3D([c3]),
                volmdlr.Face3D([c4]),
                volmdlr.Face3D([c5]),
                volmdlr.Face3D([c6])]

    def Rotation(self, center, axis, angle, copy=True):
        if copy:
            new_frame = self.frame.Rotation(center, axis, angle, copy=True)
            return Block(new_frame, self.name, self.color)
        else:
            self.frame.Rotation(center, axis, angle, copy=False)
            volmdlr.Shell3D.Rotation(self, center, axis, angle, copy=False)

    def Translation(self, offset, copy=True):
        if copy:
            new_frame = self.frame.Translation(offset, copy=True)
            return Block(new_frame, self.name, self.color)
        else:
            self.frame.Translation(offset, copy=False)
            volmdlr.Shell3D.Translation(self, offset, copy=False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'new':
            new_origin = frame.NewCoordinates(self.frame.origin)
            new_u = basis.NewCoordinates(self.frame.u)
            new_v = basis.NewCoordinates(self.frame.v)
            new_w = basis.NewCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return Block(new_frame, self.name, self.color)
            else:
                self.frame = new_frame
                volmdlr.Shell3D.frame_mapping(self, frame, side, copy=False)

        if side == 'old':
            new_origin = frame.OldCoordinates(self.frame.origin)
            new_u = basis.OldCoordinates(self.frame.u)
            new_v = basis.OldCoordinates(self.frame.v)
            new_w = basis.OldCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return Block(new_frame, self.name, self.color)
            else:
                self.frame = new_frame
                volmdlr.Shell3D.frame_mapping(self, frame, side, copy=False)

    def copy(self):
        new_origin = self.frame.origin.Copy()
        new_u = self.frame.u.Copy()
        new_v = self.frame.v.Copy()
        new_w = self.frame.w.Copy()
        new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
        return Block(new_frame, self.name, self.color)

    def plot_data(self, x3D, y3D, marker=None, color='black', stroke_width=1,
                  dash=False, opacity=1, arrow=False):
        lines = []
        for edge3D in self.Edges():
            lines.append(edge3D.plot_data(x3D, y3D, marker, color, stroke_width,
                         dash, opacity, arrow))

        return lines

    def MPLPlot2D(self, x3D, y3D, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = None

        for edge3D in self.Edges():
#            edge2D = edge3D.PlaneProjection2D()
            edge3D.MPLPlot2D(x3D, y3D, ax)

        return fig, ax



class Cone(volmdlr.Primitive3D):
    def __init__(self, position, axis, radius, length, name=''):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.position = position
        axis.Normalize()
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

    def Volume(self):
        return self.length * math.pi * self.radius**2 / 3

    def babylon_script(self):
        new_axis = volmdlr.Vector3D((self.axis[0], self.axis[1], self.axis[2]))
        normal_vector1 = new_axis.RandomUnitNormalVector()
        normal_vector2 = new_axis.Cross(normal_vector1)
        x, y, z = self.position
        s = 'var cone = BABYLON.MeshBuilder.CreateCylinder("cone", {{diameterTop:0, diameterBottom:{}, height: {}, tessellation: 100}}, scene);\n'.format(2*self.radius, self.length)
        s += 'cone.position = new BABYLON.Vector3({},{},{});\n;'.format(x,y,z)
        s += 'var axis1 = new BABYLON.Vector3({},{},{});\n'.format(new_axis[0], new_axis[1], new_axis[2])
        s += 'var axis2 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector1[0], normal_vector1[1], normal_vector1[2])
        s += 'var axis3 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector2[0], normal_vector2[1], normal_vector2[2])
        s += 'cone.rotation = BABYLON.Vector3.RotationFromAxis(axis3, axis1, axis2);\n'
        return s


class ExtrudedProfile(volmdlr.Shell3D):
    """

    """
    _non_serializable_attributes  = ['faces', 'inner_contours3d', 'outer_contour3d']
    def __init__(self, plane_origin, x, y, outer_contour2d, inner_contours2d,
                 extrusion_vector, name='', color=None):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.plane_origin = plane_origin
        self.outer_contour2d = outer_contour2d
        self.outer_contour3d = outer_contour2d.To3D(plane_origin, x, y)
        self.inner_contours2d = inner_contours2d
        self.extrusion_vector = extrusion_vector
        self.inner_contours3d = []
        self.x = x
        self.y = y
        self.color = color
        self.bounding_box = self._bounding_box()

        bool_areas = []
        for contour in inner_contours2d:
            self.inner_contours3d.append(contour.To3D(plane_origin, x, y))
            if contour.Area() > outer_contour2d.Area():
                bool_areas.append(True)
            else:
                bool_areas.append(False)
        if any(bool_areas):
            raise ValueError('At least one inner contour is not contained in outer_contour.')

        faces = self.shell_faces()
        volmdlr.Shell3D.__init__(self, faces, name)

    def shell_faces(self):

        lower_contours = [self.outer_contour3d]+self.inner_contours3d
        lower_face = volmdlr.Face3D(lower_contours)

        upper_contours = [contour.Translation(self.extrusion_vector, True) for contour in lower_contours]
        upper_face = volmdlr.Face3D(upper_contours)

        lateral_faces = []
        for i in range(len(self.inner_contours3d)+1):
            lower_points = lower_contours[i].points + [lower_contours[i].points[0]]
            upper_points = upper_contours[i].points + [upper_contours[i].points[0]]
            for j in range(len(lower_points[:-1])):
                lower_vertice1 = lower_points[j]
                lower_vertice2 = lower_points[j+1]
                upper_vertice1 = upper_points[j]
                upper_vertice2 = upper_points[j+1]
                edge1 = volmdlr.LineSegment3D(lower_vertice1, lower_vertice2)
                edge2 = volmdlr.LineSegment3D(lower_vertice2, upper_vertice2)
                edge3 = volmdlr.LineSegment3D(upper_vertice2, upper_vertice1)
                edge4 = volmdlr.LineSegment3D(upper_vertice1, lower_vertice1)
                contour = volmdlr.Contour3D([edge1, edge2, edge3, edge4])
                face = volmdlr.Face3D([contour])
                lateral_faces.append(face)
        return [lower_face]+[upper_face]+lateral_faces

    def _bounding_box(self):
        return volmdlr.BoundingBox.from_points(self.outer_contour3d.points)

    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        for contour in [self.outer_contour2d]+self.inner_contours2d:
            for primitive in contour.primitives:
                primitive.MPLPlot(ax)
        ax.margins(0.1)
        return ax

    def FreeCADExport(self, ip):
        name='primitive'+str(ip)
        s = 'Wo = []\n'
        s += 'Eo = []\n'
        for ip, primitive in enumerate(self.outer_contour3d.edges):
            s += primitive.FreeCADExport('L{}'.format(ip))
            s += 'Eo.append(Part.Edge(L{}))\n'.format(ip)
        s += 'Wo.append(Part.Wire(Eo[:]))\n'
        s += 'Fo = Part.Face(Wo)\n'

        s += 'Fi = []\n'
        s += 'W = []\n'
        for ic,contour in enumerate(self.inner_contours3d):
            s+='E = []\n'
            for ip, primitive in enumerate(contour.edges):
                s += primitive.FreeCADExport('L{}_{}'.format(ic, ip))
                s += 'E.append(Part.Edge(L{}_{}))\n'.format(ic, ip)
            s += 'Wi = Part.Wire(E[:])\n'
            s += 'Fi.append(Part.Face(Wi))\n'

        if len(self.inner_contours3d) != 0:
            s += 'Fo = Fo.cut(Fi)\n'
        e1, e2, e3 = round(1000*self.extrusion_vector, 6)

        s+='{} = Fo.extrude(fc.Vector({}, {}, {}))\n'.format(name,e1,e2,e3)
        return s

    def Area(self):
        areas = [c.Area() for c in self.contours2D]
        sic=list(npy.argsort(areas))[::-1]# sorted indices of contours
        area=areas[sic[0]]

        for i in sic[1:]:
            area-=self.contours2D[i].Area()
        return area

    def Volume(self):
#        e=extrusion_vector/norm(extrusion_vector)
        z = npy.cross(self.x.vector,self.y.vector)
        z.Normalize()
        coeff = npy.dot(self.extrusion_vector, z)

        return self.Area()*coeff


    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'old':
            extrusion_vector = basis.OldCoordinates(self.extrusion_vector)
            x = basis.OldCoordinates(self.x)
            y = basis.OldCoordinates(self.y)
        elif side == 'new':
            extrusion_vector = basis.NewCoordinates(self.extrusion_vector)
            x = basis.NewCoordinates(self.x)
            y = basis.NewCoordinates(self.y)
        else:
            raise ValueError('side must be either old or new')

        if copy:
            print('frame', frame)
            print('nxy', self.x, x, self.y, y)
            return ExtrudedProfile(self.plane_origin.frame_mapping(frame, side, copy),
                                   x,
                                   y,
                                   self.outer_contour2d,
                                   self.inner_contours2d,
                                   extrusion_vector)
        else:
            self.__init__(self.plane_origin.frame_mapping(frame, side, copy),
                          x,
                          y,
                          self.outer_contour2d,
                          self.inner_contours2d,
                          extrusion_vector)


class RevolvedProfile(volmdlr.Shell3D):
    """

    """
    _non_serializable_attributes  = ['faces', 'contour3D']
    def __init__(self, plane_origin, x, y, contour2D, axis_point,
                 axis, angle=2*math.pi, name='', color=None):
#        volmdlr.Primitive3D.__init__(self, name=name)
        self.contour2D = contour2D
        self.axis_point = axis_point
        self.axis = axis
        self.angle = angle
        self.plane_origin = plane_origin
        self.x = x
        self.y = y
        self.contour3D = self.contour2D.To3D(plane_origin, x, y)

        faces = self.shell_faces()
        volmdlr.Shell3D.__init__(self, faces, name, color)


    def shell_faces(self):
        faces = []
        number_points_for_circle = 30
        number_points_tesselation = math.ceil(number_points_for_circle*self.angle/2/math.pi)
        delta_angle = self.angle/number_points_tesselation

        for nb in range(number_points_tesselation):
            if nb == 0:
                points = self.contour3D.points
            else:
                points = [p.Rotation(self.axis_point, self.axis, nb*delta_angle, copy=True) for p in self.contour3D.points]

            rotated_points = [p.Rotation(self.axis_point, self.axis, delta_angle, copy=True) for p in points]

            points_pair = list(zip(points, rotated_points))
            for i, (pt1, pt2) in enumerate(points_pair):
                
                if i == len(points_pair)-1:
                    pt1_next = points_pair[0][0]
                    pt2_next = points_pair[0][1]
                else:
                    pt1_next = points_pair[i+1][0]
                    pt2_next = points_pair[i+1][1]

                if pt1 == pt2 and pt1_next == pt2_next:
                    continue
                
                elif pt1 == pt2:
                    edges = [volmdlr.LineSegment3D(pt1.copy(), pt2_next.copy()),
                             volmdlr.LineSegment3D(pt2_next.copy(), pt1_next.copy()),
                             volmdlr.LineSegment3D(pt1_next.copy(), pt1.copy())]
                
                elif pt1_next == pt2_next:
                    edges = [volmdlr.LineSegment3D(pt1.copy(), pt2.copy()),
                             volmdlr.LineSegment3D(pt2.copy(), pt1_next.copy()),
                             volmdlr.LineSegment3D(pt1_next.copy(), pt1.copy())]
                    
                else:
                    edges = [volmdlr.LineSegment3D(pt1.copy(), pt2.copy()),
                             volmdlr.LineSegment3D(pt2.copy(), pt2_next.copy()),
                             volmdlr.LineSegment3D(pt2_next.copy(), pt1_next.copy()),
                             volmdlr.LineSegment3D(pt1_next.copy(), pt1.copy())]

                contour = volmdlr.Contour3D(edges)
                faces.append(volmdlr.Face3D([contour]))

        return faces



    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
        for contour in self.contours3D:
#            for primitive in contour:
            contour.MPLPlot(ax)

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

    def Volume(self):
        p1 = self.axis_point.PlaneProjection3D(self.plane_origin,self.x,self.y)
        p1_2D = p1.To2D(self.axis_point,self.x,self.y)
        p2_3D = self.axis_point+volmdlr.Point3D(self.axis.vector)
        p2_2D = p2_3D.To2D(self.plane_origin,self.x,self.y)
        axis_2D = volmdlr.Line2D(p1_2D,p2_2D)
        com = self.contour2D.CenterOfMass()
        rg = axis_2D.PointDistance(com)

        return self.angle*rg*self.contour2D.Area()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        basis = frame.Basis()
        if side == 'old':
            axis = basis.OldCoordinates(self.axis)
            x = basis.OldCoordinates(self.x)
            y = basis.OldCoordinates(self.y)
        elif side == 'new':
            axis = basis.NewCoordinates(self.axis)
            x = basis.NewCoordinates(self.x)
            y = basis.NewCoordinates(self.y)
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

class Cylinder(RevolvedProfile):
    def __init__(self, position, axis, radius, length, name=''):
        self.position = position
        axis.Normalize()
        self.axis = axis
        self.radius = radius
        self.length = length
        self.bounding_box = self._bounding_box()

        # Revolved Profile
        p1 = volmdlr.Point2D((-0.5*self.length, 0))
        p2 = volmdlr.Point2D((0.5*self.length, 0))
        p3 = volmdlr.Point2D((0.5*self.length, self.radius))
        p4 = volmdlr.Point2D((-0.5*self.length, self.radius))
        l1 = volmdlr.LineSegment2D(p1, p2)
        l2 = volmdlr.LineSegment2D(p2, p3)
        l3 = volmdlr.LineSegment2D(p3, p4)
        l4 = volmdlr.LineSegment2D(p4, p1)
        contour = volmdlr.Contour2D([l1, l2, l3, l4])
        y = axis.RandomUnitNormalVector()
        RevolvedProfile.__init__(self, position, axis, y, contour, position, axis, name=name)


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

        return volmdlr.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def Volume(self):
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
        normal_vector1 = self.axis.RandomUnitNormalVector()
#        normal_vector2 = new_axis.Cross(normal_vector1)
#        x, y, z = self.position
#        s='var {} = BABYLON.Mesh.CreateCylinder("{}", {}, {}, {}, 30, 1, scene,false, BABYLON.Mesh.DEFAULTSIDE);'.format(name, self.name,self.length,2*self.radius,2*self.radius)
#        s += '{}.position = new BABYLON.Vector3({},{},{});\n;'.format(name, x,y,z)
#        s += 'var axis1 = new BABYLON.Vector3({},{},{});\n'.format(new_axis[0], new_axis[1], new_axis[2])
#        s += 'var axis2 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector1[0], normal_vector1[1], normal_vector1[2])
#        s += 'var axis3 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector2[0], normal_vector2[1], normal_vector2[2])
#        s += '{}.rotation = BABYLON.Vector3.RotationFromAxis(axis3, axis1, axis2);\n'.format(name)
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
            axis = basis.OldCoordinates(self.axis)
        elif side == 'new':
            axis = basis.NewCoordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')

        if copy:
            return Cylinder(self.position.frame_mapping(frame, side, copy),
                            axis,
                            self.radius, self.length)
        else:
            self.position.frame_mapping(frame, side, copy)
            self.axis = axis

class HollowCylinder(Cylinder):
    def __init__(self, position, axis, inner_radius, outer_radius, length, name=''):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.position = position
        axis.Normalize()
        self.axis = axis
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.length = length
        
        # Revolved Profile
        p1 = volmdlr.Point2D((-0.5*self.length, self.inner_radius))
        p2 = volmdlr.Point2D((0.5*self.length, self.inner_radius))
        p3 = volmdlr.Point2D((0.5*self.length, self.outer_radius))
        p4 = volmdlr.Point2D((-0.5*self.length, self.outer_radius))
        l1 = volmdlr.LineSegment2D(p1, p2)
        l2 = volmdlr.LineSegment2D(p2, p3)
        l3 = volmdlr.LineSegment2D(p3, p4)
        l4 = volmdlr.LineSegment2D(p4, p1)
        contour = volmdlr.Contour2D([l1, l2, l3, l4])
        y = axis.RandomUnitNormalVector()
        RevolvedProfile.__init__(self, position, axis, y, contour, position, axis, name=name)

        

    def Volume(self):
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
        normal_vector1 = self.axis.RandomUnitNormalVector()
        p1 = volmdlr.Point2D((-0.5*self.length, self.outer_radius))
        p2 = volmdlr.Point2D((0.5*self.length, self.outer_radius))
        p3 = volmdlr.Point2D((0.5*self.length, self.inner_radius))
        p4 = volmdlr.Point2D((-0.5*self.length, self.inner_radius))
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
            axis = basis.OldCoordinates(self.axis)
        elif side == 'new':
            axis = basis.NewCoordinates(self.axis)
        else:
            raise ValueError('side must be either old or new')

        if copy:
            return HollowCylinder(self.position.frame_mapping(frame, side, copy),
                                  axis,
                                  self.inner_radius, self.outer_radius, self.length)
        else:
            self.position.frame_mapping(frame, side, copy)
            self.axis = axis


class HelicalExtrudedProfile(volmdlr.Primitive3D):
    """

    """
    def __init__(self, plane_origin, x, y, axis_point,axis, pitch,
                 outer_contour2d, inner_contours2d=None, name=''):
        volmdlr.Primitive3D.__init__(self, name)
        if inner_contours2d is not None:
            self.inner_contours2d = inner_contours2d
        else:
            self.inner_contours2d = []
        self.outer_contour2d = outer_contour2d
        self.axis_point=axis_point
        self.axis=axis
        self.pitch=pitch

        self.inner_contours3d=[c.To3D(plane_origin,x,y) for c in inner_contours2d]
        self.outer_contour3d=outer_contour2d.To3D(plane_origin,x,y)


    def FreeCADExport(self,ip,ndigits=3):
        name='primitive{}'.format(ip)
        s="E = []\n"
        for icontour, contour in enumerate(self.outer_contour3d.basis_primitives):
            s += contour.FreeCADExport('L_{}'.format(icontour))
            s += 'E.append(Part.Edge(L_{}))\n'.format(icontour)
        s += 'W = Part.Wire(E[:])\n'

#        a1,a2,a3=self.axis
        ap1, ap2, ap3 = self.axis_point
        ap1 = round(ap1*1000, ndigits)
        ap2 = round(ap2*1000, ndigits)
        ap3 = round(ap3*1000, ndigits)

        width = self.axis.Norm()*1000
        direction = bool(self.pitch < 0)
        pitch = round(abs(self.pitch)*1000, ndigits)
        s += "helix = Part.makeHelix({}, {}, 50., 0, {})\n".format(pitch, width, direction)
        s += "helix.translate(fc.Vector({},{},{}))\n".format(ap1, ap2, ap3)

        s += '{} = helix.makePipeShell([W],True,True)\n'.format(name)
        for ic, contour in enumerate(self.inner_contours3d):
            s += "Ei=[]\n"
            s += "helix2 = Part.makeHelix({}, {}, 50., 0, {})\n".format(pitch,1.01*width+pitch,direction)
            s += "helix2.translate(fc.Vector({},{},{}))\n".format(ap1,ap2,ap3-pitch)
            for ip, primitive in enumerate(contour.basis_primitives):
                s += primitive.FreeCADExport('L_{}_{}'.format(ic,ip))
                s += 'Ei.append(Part.Edge(L_{}_{}))\n'.format(ic,ip)
            s+= 'Wi = Part.Wire(Ei[:])\n'
            s+= "{} = {}.cut(helix2.makePipeShell([Wi],True,True))\n".format(name,name)

        return s


class Sweep(volmdlr.Primitive3D):
    """
    Sweep a 3D contour along a Wire3D
    """
    def __init__(self, contour3d, wire3d, name=''):
        volmdlr.Primitive3D.__init__(self,name)
        self.contour3d = contour3d
        self.wire3d = wire3d

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        s = "E = []\n"
        for icontour, contour in enumerate(self.contour3d.edges):
            s += contour.FreeCADExport('L_{}'.format(icontour))
            s += 'E.append(Part.Edge(L_{}))\n'.format(icontour)
        s += 'contour = Part.Wire(E[:])\n'

        s += "E=[]\n"
        for iwire, wire in enumerate(self.wire3d.edges):
            s += wire.FreeCADExport('L_{}'.format(iwire))
            s += 'E.append(Part.Edge(L_{}))\n'.format(iwire)
        s += 'wire = Part.Wire(E[:])\n'

        s += '{} = wire.makePipeShell([contour],True, True)\n'.format(name)


        return s

class Cut(volmdlr.Primitive3D):
    """
    Cut primitive 1 by primitive 2
    """
    def __init__(self,primitive,cut_primitives,name=''):
        volmdlr.Primitive3D.__init__(self,name)
        self.primitive=primitive
        self.cut_primitives = cut_primitives


    def FreeCADExport(self,ip):
        name = 'primitive{}'.format(ip)

        s = self.primitive.FreeCADExport('{}'.format(ip))
        for icp, cut_primitive in enumerate(self.cut_primitives):
            s += cut_primitive.FreeCADExport('{}_{}'.format(ip, icp))
            s += "{} = {}.cut({}_{})\n".format(name, name, name, icp)

        return s

class Fuse(volmdlr.Primitive3D):
    """
    Fuse primitives
    """
    def __init__(self, primitives, name=''):
        volmdlr.Primitive3D.__init__(self, name)
        self.primitives = primitives


    def FreeCADExport(self,ip):
        name = 'primitive{}'.format(ip)


        s = self.primitives[0].FreeCADExport(ip)
        for primitive in self.primitives[1:]:
            s += primitive.FreeCADExport('{}_0'.format(ip))
            s += "{} = {}.fuse({}_0)\n".format(name,name,name)

        return s