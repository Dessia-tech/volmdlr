#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:08:23 2017

@author: steven
"""

import math

import numpy as npy
npy.seterr(divide='raise')

import volmdlr
from volmdlr.primitives import RoundedLineSegments


import matplotlib.pyplot as plt

class RoundedLineSegments3D(volmdlr.Wire3D, RoundedLineSegments):
    def __init__(self, points, radius, closed=False, adapt_radius=False, name=''):
        primitives = RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.LineSegment3D,
                                                  volmdlr.Arc3D,
                                                  closed, adapt_radius, name='')
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
            return RoundedLineSegments3D([p.Rotation(center, angle, copy=True)\
                                          for p in self.points],
                                         self.radius, self.closed, self.name)
        else:
            self.__init__([p.Rotation(center, angle, copy=True)\
                           for p in self.points],
                          self.radius, self.closed, self.name)

    def Translation(self, offset, copy=True):
        if copy:
            return RoundedLineSegments3D([p.Translation(offset, copy=True)\
                                          for p in self.points],
                                         self.radius, self.closed, self.name)
        else:
            self.__init__([p.Translation(offset, copy=True)\
                           for p in self.points],
                          self.radius, self.closed, self.name)

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
        x, y, z = npy.round(1000*self.center.vector, ndigits)
        return '{} = Part.makeSphere({}, fc.Vector({}, {}, {}))\n'.format(name,r,x,y,z)

class Block(volmdlr.Shell3D):
    """
    Creates a block
    :param frame: a frame 3D. The origin of the frame is the center of the block,
     the 3 vectors are defining the edges. The frame has not to be orthogonal
    """
    def __init__(self, frame, name='', color=None):
        self.frame = frame
        self.size = (self.frame.u.Norm(), self.frame.v.Norm(), self.frame.w.Norm())
        
        faces = self.shell_faces()
        volmdlr.Shell3D.__init__(self, faces, name, color)

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
        if side == 'new':
            new_origin = frame.NewCoordinates(self.frame.origin)
            new_u = frame.NewCoordinates(self.frame.u)
            new_v = frame.NewCoordinates(self.frame.v)
            new_w = frame.NewCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return Block(new_frame, self.name, self.color)
            else:
                self.frame = new_frame
                volmdlr.Shell3D.frame_mapping(self, frame, side, copy=False)
        
        if side == 'old':
            new_origin = frame.OldCoordinates(self.frame.origin)
            new_u = frame.OldCoordinates(self.frame.u)
            new_v = frame.OldCoordinates(self.frame.v)
            new_w = frame.OldCoordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return Block(new_frame, self.name, self.color)
            else:
                self.frame = new_frame
                volmdlr.Shell3D.frame_mapping(self, frame, side, copy=False)
            
    def copy(self):
        new_origin = self.frame.origin.copy()
        new_u = self.frame.u.copy()
        new_v = self.frame.v.copy()
        new_w = self.frame.w.copy()
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

class Cylinder(volmdlr.Primitive3D):
    def __init__(self, position, axis, radius, length, name=''):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.position = position
        axis.Normalize()
        self.axis = axis
        self.radius = radius
        self.length = length
        self.bounding_box = self._bounding_box()
        
    def _bounding_box(self):
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
        xmin = pointA[0] - kx * self.radius
        xmax = pointB[0] + kx * self.radius
        
        if pointA[1] > pointB[1]:
            pointA, pointB = pointB, pointA
        ymin = pointA[1] - ky * self.radius
        ymax = pointB[1] + ky * self.radius
        
        if pointA[2] > pointB[2]:
            pointA, pointB = pointB, pointA
        zmin = pointA[2] - kz * self.radius
        zmax = pointB[2] + kz * self.radius
        
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

    def Babylon(self):
        new_axis = volmdlr.Vector3D((self.axis[0], self.axis[1], self.axis[2]))
        normal_vector1 = new_axis.RandomUnitNormalVector()
        normal_vector2 = new_axis.Cross(normal_vector1)
        x, y, z = self.position
        s='var cylinder = BABYLON.Mesh.CreateCylinder("{}", {}, {}, {}, 30, 1, scene,false, BABYLON.Mesh.DEFAULTSIDE);'.format(self.name,self.length,2*self.radius,2*self.radius)
        s += 'cylinder.position = new BABYLON.Vector3({},{},{});\n;'.format(x,y,z)
        s += 'var axis1 = new BABYLON.Vector3({},{},{});\n'.format(new_axis[0], new_axis[1], new_axis[2])
        s += 'var axis2 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector1[0], normal_vector1[1], normal_vector1[2])
        s += 'var axis3 = new BABYLON.Vector3({},{},{});\n'.format(normal_vector2[0], normal_vector2[1], normal_vector2[2])
        s += 'cylinder.rotation = BABYLON.Vector3.RotationFromAxis(axis3, axis1, axis2);\n'
        return s

class HollowCylinder(volmdlr.Primitive3D):
    def __init__(self, position, axis, inner_radius, outer_radius, width, name=''):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.position = position
        axis.Normalize()
        self.axis = axis
        self.inner_radius = inner_radius
        self.outer_radius = outer_radius
        self.width = width

    def Volume(self):
        return self.width * math.pi* (self.outer_radius**2 - self.inner_radius**2)


    def FreeCADExport(self, ip):
        if self.outer_radius > 0.:
            name = 'primitive'+str(ip)
            re = round(1000*self.outer_radius, 6)
            ri = round(1000*self.inner_radius, 6)
            x, y, z = round((1000*(self.position - self.axis*self.width/2)), 6)
            ax, ay, az = npy.round(self.axis.vector, 6)

            s='C2 = Part.makeCircle({}, fc.Vector({}, {}, {}),fc.Vector({}, {}, {}))\n'.format(re, x, y, z, ax, ay, az)
            s+='W2 = Part.Wire(C2.Edges)\n'
            s+='F2 = Part.Face(W2)\n'

            if self.inner_radius!=0.:
                s+='C1 = Part.makeCircle({}, fc.Vector({}, {}, {}),fc.Vector({}, {}, {}))\n'.format(ri, x, y, z, ax, ay, az)
                s+='W1 = Part.Wire(C1.Edges)\n'
                s+='F1 = Part.Face(W1)\n'
                s+='F2 = F2.cut(F1)\n'

            vx, vy, vz = round(self.axis*self.width*1000, 6)

            s += '{} = F2.extrude(fc.Vector({}, {}, {}))\n'.format(name, vx, vy, vz)
            return s

        else:
            return ''

    def Babylon(self):
        ya,xa,za=self.axis# to counter y definition in babylon
        theta=math.acos(za/self.width)
        phi=math.atan(ya/xa)
        x,z,y=self.position
        s='var cylinder = BABYLON.Mesh.CreateCylinder("{}", {}, {}, {}, 30, 1, scene,false, BABYLON.Mesh.DEFAULTSIDE);'.format(self.name,self.width,2*self.outer_radius,2*self.outer_radius)
        s+='cylinder.position = new BABYLON.Vector3({},{},{});\n;'.format(x,y,z)
        s+='cylinder.rotation.x={}\n;'.format(-theta*math.sin(phi))
        s+='cylinder.rotation.y={}\n;'.format(theta*math.cos(phi))
        s+='cylinder.rotation.z={}\n;'.format(phi)
        return s
    

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

    def Babylon(self):
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
    def __init__(self, plane_origin, x, y, outer_contour2d, inner_contours2d,
                 extrusion_vector, name='', color=None):
        volmdlr.Primitive3D.__init__(self, name=name)
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
            
            
    def _bounding_box(self):
        return volmdlr.BoundingBox.from_points(self.outer_contour3d.points)

    def MPLPlot(self, ax=None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        for contour in [self.outer_contour2d]+self.inner_contours2d:
            for primitive in contour.primitives:
                primitive.MPLPlot(ax)
        return ax

    def FreeCADExport(self, ip):
        name='primitive'+str(ip)
        s = 'Wo = []\n'
        s += 'Eo = []\n'
        for ip, primitive in enumerate(self.outer_contour3d.basis_primitives):
            s += primitive.FreeCADExport('L{}'.format(ip))
            s += 'Eo.append(Part.Edge(L{}))\n'.format(ip)
        s += 'Wo.append(Part.Wire(Eo[:]))\n'
        s += 'Fo = Part.Face(Wo)\n'

        s += 'Fi = []\n'
        s += 'W = []\n'
        for ic,contour in enumerate(self.inner_contours3d):
            s+='E = []\n'
            for ip, primitive in enumerate(contour.basis_primitives):
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
    
    def Babylon(self):
        s = 'var mat = new BABYLON.StandardMaterial("mat", scene);\n'
        s += 'mat.backFaceCulling = false;\n'
        if self.color is not None:
            s += 'mat.diffuseColor = new BABYLON.Color3({},{},{});\n'.format(self.color[0], self.color[1], self.color[2])
            
            
        lower_outer_ribbon_points = self.outer_contour3d.points
        upper_outer_ribbon_points = [p.Translation(self.extrusion_vector) for p in lower_outer_ribbon_points]
        
        s += 'var LowerOuterPathRibbon = [];\n'
        for point in lower_outer_ribbon_points:
            s += 'LowerOuterPathRibbon.push(new BABYLON.Vector3({},{},{}));\n'.format(round(point[0],6), round(point[1],6), round(point[2],6))
        s += 'var UpperOuterPathRibbon = [];\n'
        for point in upper_outer_ribbon_points:
            s += 'UpperOuterPathRibbon.push(new BABYLON.Vector3({},{},{}));\n'.format(round(point[0],6), round(point[1],6), round(point[2],6))
        s += 'var outerRibbon = BABYLON.MeshBuilder.CreateRibbon("ribbon", {pathArray: [LowerOuterPathRibbon, UpperOuterPathRibbon]}, scene);\n'
        s += 'outerRibbon.material = mat;\n'
        
        
        lower_inner_ribbon_contours_points = []
        upper_inner_ribbon_contours_points = []
        for inner_contour3d in self.inner_contours3d:
            lower_inner_ribbon_points = inner_contour3d.points
            upper_inner_ribbon_points = [p.Translation(self.extrusion_vector) for p in lower_inner_ribbon_points]
            lower_inner_ribbon_contours_points.append(lower_inner_ribbon_points)
            upper_inner_ribbon_contours_points.append(upper_inner_ribbon_points)
                
            s += 'var LowerInnerPathRibbon = [];\n'
            for point in lower_inner_ribbon_points:
                s += 'LowerInnerPathRibbon.push(new BABYLON.Vector3({},{},{}));\n'.format(round(point[0],6), round(point[1],6), round(point[2],6))
            s += 'var UpperInnerPathRibbon = [];\n'
            for point in upper_inner_ribbon_points:
                s += 'UpperInnerPathRibbon.push(new BABYLON.Vector3({},{},{}));\n'.format(round(point[0],6), round(point[1],6), round(point[2],6))
            s += 'var innerRibbon = BABYLON.MeshBuilder.CreateRibbon("ribbon", {pathArray: [LowerInnerPathRibbon, UpperInnerPathRibbon]}, scene);\n'
            s += 'innerRibbon.material = mat;\n'
            
            
        lower_plane = volmdlr.Plane3D.from_points(lower_outer_ribbon_points[:-1])
        lower_outer_polygon_points = [p.To2D(lower_plane.origin, lower_plane.vectors[0], lower_plane.vectors[1]) for p in lower_outer_ribbon_points[:-1]]
        s += 'LowerPolygonPoints = [];\n'
        for point in lower_outer_polygon_points:
            s += 'LowerPolygonPoints.push(new BABYLON.Vector2({},{}));\n'.format(round(point[0],6), round(point[1],6))
        s += 'var lowerPoly_tri = new BABYLON.PolygonMeshBuilder("lowerPolygon", LowerPolygonPoints, scene);\n'
        for inner_contour in lower_inner_ribbon_contours_points:
            lower_inner_polygon_points = [p.To2D(lower_plane.origin, lower_plane.vectors[0], lower_plane.vectors[1]) for p in inner_contour]
            s += 'var lower_hole = [];\n'
            for point in lower_inner_polygon_points[:0:-1]:
                s += 'lower_hole.push(new BABYLON.Vector2({},{}));\n'.format(round(point[0],6), round(point[1],6))
            s += 'lowerPoly_tri.addHole(lower_hole);\n'
        s += 'var lowerPolygon = lowerPoly_tri.build(true, 0);\n'
        s += 'lowerPolygon.material = mat;\n'
        
        x = lower_outer_ribbon_points[0][0]
        y = lower_outer_ribbon_points[0][1]
        z = lower_outer_ribbon_points[0][2]
        axis1 = lower_plane.vectors[0]
        axis2 = lower_plane.vectors[1]
        axis3 = -lower_plane.normal
        s += 'lowerPolygon.position = new BABYLON.Vector3({},{},{});\n'.format(x,y,z)
        s += 'var axis1 = new BABYLON.Vector3({},{},{});\n'.format(*axis3)
        s += 'var axis2 = new BABYLON.Vector3({},{},{});\n'.format(*axis2)
        s += 'var axis3 = new BABYLON.Vector3({},{},{});\n'.format(*axis1)
        s += 'lowerPolygon.rotation = BABYLON.Vector3.RotationFromAxis(axis3, axis1, axis2);\n'
        
        
        upper_plane = volmdlr.Plane3D.from_points(upper_outer_ribbon_points[:-1])
        upper_outer_polygon_points = [p.To2D(upper_plane.origin, upper_plane.vectors[0], upper_plane.vectors[1]) for p in upper_outer_ribbon_points[:-1]]
        s += 'UpperPolygonPoints = [];\n'
        for point in upper_outer_polygon_points:
            s += 'UpperPolygonPoints.push(new BABYLON.Vector2({},{}));\n'.format(round(point[0],6), round(point[1],6))
        s += 'var upperPoly_tri = new BABYLON.PolygonMeshBuilder("upperPolygon", UpperPolygonPoints, scene);\n'
        for inner_contour in upper_inner_ribbon_contours_points:
            upper_inner_polygon_points = [p.To2D(upper_plane.origin, upper_plane.vectors[0], upper_plane.vectors[1]) for p in inner_contour]
            s += 'var upper_hole = [];\n'
            for point in upper_inner_polygon_points[:0:-1]:
                s += 'upper_hole.push(new BABYLON.Vector2({},{}));\n'.format(round(point[0],6), round(point[1],6))
            s += 'upperPoly_tri.addHole(upper_hole);\n'
        s += 'var upperPolygon = upperPoly_tri.build(true, 0);\n'
        s += 'upperPolygon.material = mat;\n'
        
        x = upper_outer_ribbon_points[0][0]
        y = upper_outer_ribbon_points[0][1]
        z = upper_outer_ribbon_points[0][2]
        axis1 = upper_plane.vectors[0]
        axis2 = upper_plane.vectors[1]
        axis3 = -upper_plane.normal
        s += 'upperPolygon.position = new BABYLON.Vector3({},{},{});\n'.format(x,y,z)
        s += 'var axis1 = new BABYLON.Vector3({},{},{});\n'.format(*axis3)
        s += 'var axis2 = new BABYLON.Vector3({},{},{});\n'.format(*axis2)
        s += 'var axis3 = new BABYLON.Vector3({},{},{});\n'.format(*axis1)
        s += 'upperPolygon.rotation = BABYLON.Vector3.RotationFromAxis(axis3, axis1, axis2);\n'

        return s

class RevolvedProfile(volmdlr.Shell3D):
    """

    """
    def __init__(self, plane_origin, x, y, contour2D, axis_point,
                 axis, angle=2*math.pi, name='', color=None):
        volmdlr.Primitive3D.__init__(self, name=name)
        self.contour2D = contour2D
        self.axis_point = axis_point
        self.axis = axis
        self.angle = angle
        self.plane_origin = plane_origin
        self.x = x
        self.y = y

#        self.contours3D = []
#        for contour in contour2D:
#            self.contours3D.append(contour.To3D(plane_origin, x, y))
        self.contour3D = self.contour2D.To3D(plane_origin, x, y)
            
        faces = self.shell_faces()
        volmdlr.Shell3D.__init__(self, faces, name, color)

    
    def shell_faces(self):
        faces = []
        TESSELLATION = 40
        delta_angle = self.angle/TESSELLATION
        
        for nb in range(TESSELLATION):
            if nb == 0:
                points = self.contour3D.points
            else:
                points = [p.Rotation(self.axis_point, self.axis, nb*delta_angle, copy=True) for p in self.contour3D.points]
                
            rotated_points = [p.Rotation(self.axis_point, self.axis, delta_angle, copy=True) for p in points]
            
            points_pair = list(zip(points, rotated_points))
            for i, (pt1, pt2) in enumerate(points_pair[:-1]):
                pt1_next = points_pair[i+1][0]
                pt2_next = points_pair[i+1][1]
                
                if pt1 == pt2 and pt1_next == pt2_next:
                    continue
                
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
        for ibp, basis_primitive in enumerate(self.contour2D.basis_primitives):
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

#            myObject.Shape = Sweep
        return s

    def Volume(self):
#        areas=[c.Area() for c in self.contours2D]
#        # Maximum area is main surface, others cut into it
#        sic = list(npy.argsort(areas))[::-1]# sorted indices of contours
        p1=self.axis_point.PlaneProjection3D(self.plane_origin,self.x,self.y)
#        if self.axis_point.PointDistance(p1)!=0:
#            raise NotImplementedError
        p1_2D=p1.To2D(self.axis_point,self.x,self.y)
        p2_3D=self.axis_point+volmdlr.Point3D(self.axis.vector)
#        p2=p2_3D.PlaneProjection3D(self.plane_origin,self.x,self.y)
#        if p2_3D.PointDistance(p2)!=0:
#            raise NotImplementedError
        p2_2D=p2_3D.To2D(self.plane_origin,self.x,self.y)
        axis_2D=volmdlr.Line2D(p1_2D,p2_2D)
        com = self.contour2D.CenterOfMass()
        rg = axis_2D.PointDistance(com)
#        volume=areas[sic[0]]*rg
#
#        for i in sic[1:]:
#            com=self.contours2D[i].CenterOfMass()
#            rg=axis_2D.PointDistance(com)
#            volume-=areas[i]*rg

#        return self.angle*volume
        return self.angle*rg*self.contour2D.Area()



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
        for icontour, contour in enumerate(self.contour3d.basis_primitives):
            s += contour.FreeCADExport('L_{}'.format(icontour))
            s += 'E.append(Part.Edge(L_{}))\n'.format(icontour)
        s += 'contour = Part.Wire(E[:])\n'

        s += "E=[]\n"
        for iwire, wire in enumerate(self.wire3d.basis_primitives):
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