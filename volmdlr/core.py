#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:07:37 2017

@author: steven
"""

import math
import numpy as npy
import matplotlib.pyplot as plt
from matplotlib.patches import Arc
from mpl_toolkits.mplot3d import Axes3D 

from .vmcy import PolygonPointBelongs

from scipy.linalg import solve,LinAlgError, inv

import volmdlr.geometry as geometry

from jinja2 import Environment, PackageLoader, select_autoescape

import webbrowser
import os

import tempfile
import subprocess

class Vector:
    """
    Abstract class of vector
    """
    def __setitem__(self, key, item):
        self.vector[key] = item

    def __getitem__(self, key):
        return self.vector[key]
    
    def __repr__(self):
        return '{}: {}'.format(self.__class__.__name__, self.vector)
    
    def __add__(self, other_vector):
        return self.__class__(self.vector + other_vector.vector)
    
    def __radd__(self, other_vector):
        return self + other_vector

    def __sub__(self, other_vector):
        return self.__class__(self.vector - other_vector.vector)
    
    def __rsub__(self, other_vector):
        return self - other_vector    
    
    def __neg__(self):
        return self.__class__(-self.vector)
    
    def __mul__(self, value):
        return self.__class__(self.vector * value)
    
    def __rmul__(self, value):
        return self * value

    def __truediv__(self, value):
        return self.__class__(self.vector / value)
    
    def __rtruediv__(self, value):
        return self / value
    
    def __lt__(self, other_vector):
        return self.Norm() < other_vector.Norm()

    def __le__(self, other_vector):
        return self.Norm() <= other_vector.Norm()
    
    def __ne__(self, other_vector):
        return not npy.allclose(self.vector, other_vector.vector)

    def __eq__(self, other_vector):
        try:
            return npy.allclose(self.vector, other_vector.vector)
        except AttributeError:
            return False

    def __hash__(self):
        return int(1000*npy.sum(npy.round(self.vector,3)))
    
    def Normalize(self):
        """
        Normalize the vector modifying it's coordinate
        """
        self.vector /= self.Norm()


class Vector2D(Vector):
    def __init__(self, vector):
        self.vector = npy.zeros(2)
        self.vector[0] = vector[0]        
        self.vector[1] = vector[1]
        
    def __round__(self, ndigits):
        return self.__class__((round(self.vector[0], ndigits),
                               round(self.vector[1], ndigits)))
        
    def Norm(self):
        """
        :returns: norm of vector
        """
        x, y = self.vector
        return (x**2 + y**2)**0.5

    def Dot(self, other_vector):
        u1, u2 = self.vector
        v1, v2 = other_vector.vector
        return u1*v1 + u2*v2
    
    def Rotation(self,center, angle, copy=True):
        vector2 = (npy.dot(npy.array([[math.cos(angle),-math.sin(angle)], 
                                       [math.sin(angle),math.cos(angle)]]),
                          (self.vector-center.vector))
                   + center.vector)
        if copy:
            return self.__class__(vector2)
        else:
            self.vector=vector2

    def Translation(self, offset, copy=True):
        """
        :param offset: an other Vector2D
        """
        vector2 = self.vector + offset.vector
        if copy:
            return self.__class__(vector2)
        else:
            self.vector = vector2
        
    def To3D(self,plane_origin, x1, x2):
        x, y = self.vector
        return Vector3D(plane_origin.vector + x1.vector*x + x2.vector*y)
    
    def NormalVector(self, unit=False):
        n = Vector2D((-self.vector[1], self.vector[0]))
        if unit:
            n.Normalize()
        return n
    
o2D = Vector2D((0, 0))
x2D = Vector2D((1, 0))
y2D = Vector2D((0, 1))
    
class Point2D(Vector2D):
    def __init__(self, vector, name=''):
        Vector2D.__init__(self, vector)
        self.name = name
        
    def __add__(self,point2d):
        return Point2D(self.vector + point2d.vector)
    
    def __sub__(self,point2d):
        return Point2D(self.vector - point2d.vector)
        
    def __mul__(self,value):
        return Point2D(self.vector * value)
    
    def __truediv__(self,value):
        return Point2D(self.vector / value)
    
    def To3D(self,plane_origin, x1, x2):
        x, y = self.vector
        return Point3D(plane_origin.vector + x1.vector*x + x2.vector*y)

    def MPLPlot(self, ax, style='ob'):
        x1=self.vector
        ax.plot([x1[0]], [x1[1]], style)        
        return []
    
    def PointDistance(self, point2):
        return (self-point2).Norm()

    #def Distance(self,point2):
    #    return norm(self.vector-point2.vector)

    @classmethod
    def LinesIntersection(cls, line1,line2, curvilinear_abscissa=False):
        p11=line1.points[0].vector
        p12=line1.points[1].vector
        p21=line2.points[0].vector
        p22=line2.points[1].vector
        A=npy.array([[p12[0]-p11[0],p21[0]-p22[0]],[p12[1]-p11[1],p21[1]-p22[1]]])
        x=npy.array([p21[0]-p11[0],p21[1]-p11[1]])
        try:
            t=solve(A,x)
            if not curvilinear_abscissa:
                return cls(p11+t[0]*(p12-p11))
            else:
                return (cls(p11+t[0]*(p12-p11)),t[0],t[1])
        except LinAlgError:
            return None
        
    @classmethod
    def MiddlePoint(cls,point1,point2):
        p1=point1.vector
        p2=point2.vector
        return cls((p1+p2)*0.5)

    @classmethod
    def LineProjection(cls,point,line):
        p1,p2=line.points
        d=(p2-p1)/p2.PointDistance(p1)
        n=d.Rotation(Point2D((0,0)),math.pi/2).vector
        pp1=point.vector-p1.vector
        p=pp1-npy.dot(pp1,n)*n+p1.vector
        return Point2D(p)


class Basis2D:
    """
    Defines a 2D basis
    :param u: first vector of the basis
    :param v: second vector of the basis
    """
    def __init__(self, u, v):
        self.u = u
        self.v = v
        
    def __repr__(self):
        return '{}: U={}, V={}'.format(self.__class__.__name__, self.u, self.v)
    
   
    def TransfertMatrix(self):
        return npy.array([[self.u[0], self.v[0]],
                          [self.u[1], self.v[1]]])
    
    def InverseTransfertMatrix(self):
        # Todo: cache for performance
        return inv(self.TransfertMatrix())
        
    def NewCoordinates(self, vector):
        return Vector2D(npy.dot(self.InverseTransfertMatrix(), vector.vector))

    def OldCoordinates(self, vector):
        return Vector2D(npy.dot(self.TransfertMatrix(), vector.vector))
    
xy = Basis2D(x2D, y2D)
     
class Frame2D(Basis2D):
    """
    Defines a 2D basis
    :param origin: origin of the basis
    :param u: first vector of the basis
    :param v: second vector of the basis
    """
    def __init__(self, origin, u, v):
        self.origin = origin
        Basis2D.__init__(self, u, v)

    def __repr__(self):
        return '{}: O={} U={}, V={}'.format(self.__class__.__name__, self.origin, self.u, self.v)

    def Basis(self):
        return Basis3D(self.u, self.v)
        
    def NewCoordinates(self, vector):
        return Basis2D.NewCoordinates(self, vector - self.origin)

    def OldCoordinates(self, vector):
        return Basis2D.OldCoordinates(self, vector) + self.origin

        
oxy = Frame2D(o2D, x2D, y2D)    
    
class Primitive2D:
    def __init__(self, name=''):
        self.name=name
        
class CompositePrimitive2D(Primitive2D):
    """
    A collection of simple primitives
    """
    def __init__(self,primitives, name=''):
        Primitive2D.__init__(self, name)        
        self.primitives=primitives
        self.UpdateBasisPrimitives()
        
    def UpdateBasisPrimitives(self):
        basis_primitives=[]
        for primitive in self.primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.basis_primitives)
            else:
                basis_primitives.append(primitive)
                
        self.basis_primitives = basis_primitives
        
        
    def Rotation(self,center,angle,copy=False):
        if copy:
            return self.__class__([p.Rotation(center,angle,copy=True) for p in self.primitives])
        else:
            for p in self.basis_primitives:
                p.Rotation(center,angle,copy=False)
            self.UpdateBasisPrimitives()
            
    def Translation(self, offset, copy=False):
        if copy:
            return self.__class__([p.Translation(offset, copy=True) for p in self.primitives])
        else:
            for p in self.basis_primitives:
                p.Translation(offset, copy=False)
            self.UpdateBasisPrimitives()
    
    def To3D(self, plane_origin, x, y, name = None):
        if name is None:
            name = '3D of {}'.format(self.name)
        primitives3D = [p.To3D(plane_origin, x, y) for p in self.primitives]
        return CompositePrimitive3D(primitives3D, name)
        
    def MPLPlot(self, ax = None):
        if ax is None:
            fig, ax = plt.subplots()
            ax.set_aspect('equal')
        else:
            fig = None
            
        for element in self.basis_primitives:
            element.MPLPlot(ax)
        
        ax.margins(0.1)
        plt.show() 
        
        return fig, ax
        
        
class Wire2D(CompositePrimitive2D):
    """
    A collection of simple primitives, following each other making a wire
    """
    def __init__(self, primitives, name=''):
        CompositePrimitive2D.__init__(self, primitives, name)        
    
    # TODO: method to check if it is a wire
        
    def Length(self):
        length = 0.
        for primitive in self.basis_primitives:
            length += primitive.Length()
        return length
    
    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        length = 0.
        for primitive in self.basis_primitives:
            primitive_length = primitive.Length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(curvilinear_abscissa - length)
            length += primitive_length
        return ValueError
    
    
class Contour2D(Wire2D):
    """
    A collection of 3D primitives forming a closed wire3D
    """
    def __init__(self, primitives, name=''):
        Wire2D.__init__(self, primitives, name)
        
    def To3D(self, plane_origin, x, y, name = None):
        if name is None:
            name = '3D of {}'.format(self.name)
        primitives3D = [p.To3D(plane_origin, x, y) for p in self.primitives]
        return Contour3D(primitives3D, name)
    
    def Area(self):
        if len(self.primitives)==1:
            return self.primitives[0].Area()
        arcs=[]
        points_polygon=[]
        for primitive in self.primitives:
            if primitive.__class__.__name__=='LineSegment2D':
                points_polygon.extend(primitive.points)
            elif primitive.__class__.__name__=='Arc2D':
                points_polygon.append(primitive.center)
                arcs.append(primitive)
        polygon=Polygon2D(points_polygon)
        A=polygon.Area()

        for arc in arcs:
            if polygon.PointBelongs(arc.interior):
                A-=arc.Area()
            else:
                A+=arc.Area()
                
        return A
    
    def CenterOfMass(self):
        if len(self.primitives)==1:
            return self.primitives[0].CenterOfMass()

        arcs=[]
        points_polygon=[]
        for primitive in self.primitives:
            if primitive.__class__.__name__=='LineSegment2D':
                points_polygon.extend(primitive.points)
            elif primitive.__class__.__name__=='Arc2D':
                points_polygon.append(primitive.center)
                arcs.append(primitive)
        polygon=Polygon2D(points_polygon)
        
        area=polygon.Area()
        c=area*polygon.CenterOfMass()
        
        for arc in arcs:
            arc_area=arc.Area()
            if polygon.PointBelongs(arc.middle):
                c-=arc_area*arc.CenterOfMass()
                area-=arc_area
            else:
                c+=arc_area*arc.CenterOfMass()
                area+=arc_area
        return c/area



    def SecondMomentArea(self,point):
        if len(self.primitives)==1:
            return self.primitives[0].SecondMomentArea(point)

        arcs=[]
        points_polygon=[]
        for primitive in self.primitives:
            if primitive.__class__.__name__=='Line2D':
                points_polygon.extend(primitive.points)
            elif primitive.__class__.__name__=='Arc2D':
                points_polygon.append(primitive.center)
                arcs.append(primitive)
        polygon=Polygon2D(points_polygon)
        A=polygon.SecondMomentArea(point)
        for arc in arcs:
            if polygon.PointBelongs(arc.middle):
                A-=arc.SecondMomentArea(point)
            else:
                A+=arc.SecondMomentArea(point)
        return A
    
    


class Mesh2D:
    def __init__(self, contours, points_densities, default_density):
        self.contours = contours
        self.points_densities = points_densities
        self.default_density = default_density
        
    def GeoScript(self, filepath=''):
        s=''
        ipt=1# point index
        ipr=1# primitive index 
        points_index={}
        #assigning an index to point    
        for contour in self.contours:
            for primitive in contour.primitives:
                for point in primitive.geo_points:
                    try:
                        points_index[point]
                    except KeyError:
                        points_index[point]=ipt
                        try:
                            d=self.points_densities[point]
                        except KeyError:
                            d=self.default_density
                        s+='Point({})={{{},{},0.,{}}};\n'.format(ipt,*point.vector,d)
                        ipt+=1
        contours_indices=[]
        for contour in self.contours:
            contour_iprs=[]
            for primitive in contour.primitives:
                spr,ipr2=primitive.GeoScript(ipr,[points_index[p] for p in primitive.geo_points])
                s+=spr
                contour_iprs.extend(range(ipr,ipr2))
                ipr=ipr2
            s+='Line Loop({}) = {{{}}};\n'.format(ipr,str(contour_iprs)[1:-1])
            contours_indices.append(ipr)
            ipr+=1
        s+='Plane Surface({}) = {{{}}};\n'.format(ipr,str(contours_indices)[1:-1])
        # Saving to file if required
        if filepath!='':
            with open(filepath,'w') as file:
                file.write(s)
        return s

class Line:
    def DirectionVector(self, unit=False):
        u = self.points[1] - self.points[0]
        if unit:
            u.Normalize()
        return u
    
    def NormalVector(self, unit=False):
        return self.DirectionVector(unit).NormalVector()

class Line2D(Primitive2D, Line):
    """
    Define an infinte line given by two points.
    """
    def __init__(self, point1, point2, name=''):
        Primitive2D.__init__(self, name)        
        self.points=[point1, point2]

    def To3D(self, plane_origin, x1, x2):
        p3D = [p.To3D(plane_origin, x1, x2) for p in self.points]
        return Line2D(*p3D, self.name)
    
    def Rotation(self, center, angle, copy=False):
        if copy:
            return Line2D(*[p.Rotation(center, angle, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center, angle, copy=False)
            
    def Translation(self, offset, copy=False):
        if copy:
            return Line2D(*[p.Translation(offset, copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset, copy=False)
    
    def PointDistance(self, point):
        """
        Computes the distance of a point to line 
        """
        p1, p2=self.points
        t = (point-p1).Dot(p2-p1)/ (p2-p1).Norm()**2
        projection = p1 + t * (p2-p1)# Projection falls on the segment
        return (point-projection).Norm()
    
    def PointProjection(self, point, curvilinear_abscissa=False):
        p1, p2 = self.points
        t = (point - p1).Dot(p2 - p1) / (p2-p1).Norm()**2
        projection = p1 + t * (p2-p1)
        if curvilinear_abscissa:
            return projection,t
        return projection
    
    def MPLPlot(self, ax):
        p1, p2 = self.points
        u = p2 - p1
        plt.plot([p1[0], p2[0]], [p1[1], p2[1]], 'ok')        
        p3 = p1 - 3* u
        p4 = p2 + 4*u
        ax.plot([p3[0], p4[0]], [p3[1], p4[1]], '-k', linestyle = '-.')        
        return []

class LineSegment2D(Line2D):
    """
    Define a line segment limited by two points
    """
    def __init__(self,point1, point2,name=''):
        Line2D.__init__(self, point1, point2, name = name)
        
    def _get_geo_points(self):
        return self.points

    geo_points=property(_get_geo_points)      
    
    def Length(self):
        return self.points[1].PointDistance(self.points[0])
    
    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.points[0] + self.DirectionVector(unit=True) * curvilinear_abscissa
    
    def PointDistance(self, point):
        """
        Computes the distance of a point to segment of line 
        """
        p1,p2=self.points
        v=p1.vector
        w=p2.vector
        p=point.vector
        nwv2=(w[1]-v[1])**2+(w[0]-v[0])**2# replace norm(w-v)**2
        t = max(0, min(1, npy.dot(p - v, w - v) / nwv2))

        projection = v + t * (w - v)# Projection falls on the segment
        return ((p[1]-projection[1])**2+(p[0]-projection[0])**2)**0.5

    def PointProjection(self, point, curvilinear_abscissa=False):
        point, curv_abs = Line2D.PointProjection(self, point, True)
        if curv_abs <= 0.:
            point = self.points[0]
            curv_abs = 0.
        elif curv_abs >= 1.:
            point = self.points[1]
            curv_abs = 1.
            
        if curvilinear_abscissa:
            return point, curv_abs
        else:
            return point
        
    def MPLPlot(self, ax, style='-k'):
        p1, p2 = self.points
        ax.plot([p1[0], p2[0]], [p1[1], p2[1]], style)        
        return []
    
    def To3D(self, plane_origin, x1, x2):
        p3D=[p.To3D(plane_origin,x1,x2) for p in self.points]
        return LineSegment3D(*p3D,self.name)
    
    def Rotation(self, center, angle, copy=False):
        if copy:
            return LineSegment2D(*[p.Rotation(center,angle,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center,angle,copy=False)
            
    def Translation(self, offset, copy=False):
        if copy:
            return LineSegment2D(*[p.Translation(offset,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset,copy=False)
                
    def GeoScript(self, primitive_index, points_indices):
        s='Line({}) = {{{}, {}}};\n'.format(primitive_index,*points_indices)
        return s,primitive_index+1
                

                
class Arc2D(Primitive2D):
    def __init__(self, start, interior, end, name=''):        
        Primitive2D.__init__(self, name)        
        self.interior = interior
        self.start = start
        self.end = end
        xi, yi = interior.vector
        xe, ye = end.vector
        xs, ys = start.vector
        A = npy.array([[2*(xs-xi), 2*(ys-yi)],
                       [2*(xs-xe), 2*(ys-ye)]])
        b = - npy.array([xi**2 + yi**2 - xs**2 - ys**2,
                         xe**2 + ye**2 - xs**2 - ys**2])
        self.center = Point2D(solve(A,b))
        r1 = self.start - self.center
        r2 = self.end - self.center
        ri = self.interior - self.center
        
        self.radius = r1.Norm()
        angle1 = npy.arctan2(r1.vector[1], r1.vector[0])
        angle2 = npy.arctan2(r2.vector[1], r2.vector[0])
        
        anglei = npy.arctan2(ri.vector[1], ri.vector[0])
        order = [y for x, y in sorted(zip([angle1, anglei, angle2], [0, 1, 2]))]
        order = order*2
        i = order.index(0)
        if order[i+1] == 1:
            # Trigo wise angle should be plus
            self.angle1 = angle1
            self.angle2 = angle2
            if angle1 > angle2:
                self.angle = angle2 - angle1 + 2 * math.pi
            else:
                self.angle = angle2 - angle1
#            self.angle = abs(self.angle2 - self.angle1)
        else:
            # Clock wise
            self.angle1 = angle2
            self.angle2 = angle1
#            self.angle = -abs(self.angle2 - self.angle1)
            if angle1 < angle2:
                self.angle = angle2 - angle1 + 2 * math.pi
            else:
                self.angle = angle2 - angle1 

#        if self.angle1 < 0:
#            a1 = self.angle1 + 2* math.pi
#        else:
#            a1 = self.angle1
#        if self.angle2 < 0:
#            a2 = self.angle2 + 2* math.pi
#        else:
#            a2 = self.angle2
        
        
    def _get_geo_points(self):
        return [self.start,self.center,self.end]

    geo_points=property(_get_geo_points)        
    
    def Length(self):
        return self.radius * abs(self.angle)
        
    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        if self.angle>0:
            return self.start.Rotation(self.center, curvilinear_abscissa/self.radius)
        else:
            return self.start.Rotation(self.center, -curvilinear_abscissa/self.radius)
            
    def GeoScript(self, primitive_index, points_indices):
        s='Circle({}) = {{{}, {}, {}}};\n'.format(primitive_index,*points_indices)
        return s,primitive_index+1
            
    def Area(self):
        if self.angle2<self.angle1:
            angle=self.angle2+2*math.pi-self.angle1
        else:
            angle=self.angle2-self.angle1
        return self.radius**2*angle/2
            
    def CenterOfMass(self):
        u=self.middle.vector-self.center.vector
        u.Normalize()
        alpha=abs(self.angle1-self.angle2)
        return Point2D(self.center.vector+4/(3*alpha)*self.radius*math.sin(alpha*0.5)*u)
        
    def MPLPlot(self, ax, style='-k'):
        pc = self.center.vector
#        ax.plot([pc[0]], [pc[1]], 'or')
#        ax.plot([self.interior[0]], [self.interior[1]], 'ob')
        ax.add_patch(Arc(pc, 2*self.radius, 2*self.radius, angle=0, 
                    theta1=self.angle1*0.5/math.pi*360,
                    theta2=self.angle2*0.5/math.pi*360,
                    color='k'))

    def To3D(self,plane_origin, x, y):
        ps = self.start.To3D(plane_origin, x, y)
        pi = self.interior.To3D(plane_origin, x, y)
        pe = self.end.To3D(plane_origin, x, y)
        
#        pe = Point2D(self.center.vector + self.radius*npy.array((math.cos(self.angle1),math.sin(self.angle1))))
#        pe3 = pe.To3D(plane_origin, x, y)
#        ps = Point2D(self.center.vector+self.radius*npy.array((math.cos(self.angle2),math.sin(self.angle2))))
#        ps3 = ps.To3D(plane_origin, x, y)
        return Arc3D(ps, pi, pe, self.name)
    
    def Rotation(self, center, angle, copy=False):
        if copy:
            return Arc2D(*[p.Rotation(center,angle,copy=True) for p in [self.start,self.interior,self.end]])
        else:
            self.__init__(*[p.Rotation(center,angle,copy=True) for p in [self.start,self.interior,self.end]])
            
    def Translation(self, offset, copy=False):
        if copy:
            return Arc2D(*[p.Translation(offset,copy=True) for p in [self.start,self.interior,self.end]])
        else:
            self.__init__(*[p.Translation(offset,copy=True) for p in [self.start,self.interior,self.end]])

    def SecondMomentArea(self, point):
        """ 
        Second moment area of part of disk
        """
        if self.angle2<self.angle1:
            angle2=self.angle2+2*math.pi
            
        else:
            angle2=self.angle2
        angle1=self.angle1

        Ix=self.radius**4/8*(angle2-angle1+0.5*(math.sin(2*angle1)-math.sin(2*angle2)))
        Iy=self.radius**4/8*(angle2-angle1+0.5*(math.sin(2*angle2)-math.sin(2*angle1)))
        Ixy=self.radius**4/8*(math.cos(angle1)**2-math.cos(angle2)**2)
        Ic=npy.array([[Ix,Ixy],[Ixy,Iy]])
        return geometry.Huygens2D(Ic, self.Area(), self.center, point)



class Circle2D(Primitive2D):
    def __init__(self,center,radius,name=''):        
        Primitive2D.__init__(self,name)        
        self.center=center
        self.radius=radius
        self.utd_geo_points=False
        
    def _get_geo_points(self):
        if not self.utd_geo_points:
            self._geo_start=self.center+self.radius*Point2D((1,0))
            self.utd_geo_points=True
        return [self._geo_start,self.center,self._geo_start]

    geo_points = property(_get_geo_points)        
        
    def Length(self):
        return 2* math.pi * self.radius
    
    def GeoScript(self, primitive_index, points_indices):
        s = 'Circle({}) = {{{}, {}, {}}};\n'.format(primitive_index,*points_indices)
        return s, primitive_index+1
    
    def MPLPlot(self, ax, style='-k'):
        pc = self.center.vector
        return [Arc(pc,2*self.radius,2*self.radius,angle=0,theta1=0,theta2=360,color='k')]

    def To3D(self, plane_origin, x, y):
        normal = Vector3D(npy.cross(x.vector, y.vector))
        pc=self.center.To3D(plane_origin, x, y)
        return Circle3D(pc,self.radius,normal, self.name)

    def Rotation(self, center, angle, copy=False):
        if copy:
            return Circle2D(self.center.Rotation(center,angle,copy=True),self.radius)
        else:
            self.center.Rotation(center,angle,copy=False) 
            self.utd_geo_points=False
            
    def Translation(self,offset,copy=False):
        if copy:
            return Circle2D(self.center.Translation(offset,copy=True),self.radius)
        else:
            self.center.Translation(offset,copy=False)
            self.utd_geo_points=False

    def Area(self):
        return math.pi*self.radius**2

    def SecondMomentArea(self, point):
        """ 
        Second moment area of part of disk
        """
        I = math.pi*self.radius**4/4
        Ic = npy.array([[I,0],[0,I]])
        return geometry.Huygens2D(Ic,self.Area(),self.center,point)
    
    def CenterOfMass(self):
        return self.center

class Polygon2D(CompositePrimitive2D):
    # TODO: inherit from contour?
    def __init__(self,points,name=''):     
        self.points=points
        primitives=[]
        for p1,p2 in zip(points,points[1:]+[points[0]]):
            primitives.append(LineSegment2D(p1,p2))
            
        self.line_segments = self._LineSegments()

        CompositePrimitive2D.__init__(self, primitives, name)
        
        
    def Area(self):
        
        x=[point.vector[0]for point in self.points]
        y=[point.vector[1]for point in self.points]

        return 0.5*npy.abs(npy.dot(x,npy.roll(y,1))-npy.dot(y,npy.roll(x,1)))
    
    def CenterOfMass(self):
        
        x = [point.vector[0] for point in self.points]
        y = [point.vector[1] for point in self.points]

        
        xi_xi1 = x+npy.roll(x,-1)
        yi_yi1 = y+npy.roll(y,-1)
        xi_yi1 = npy.multiply(x,npy.roll(y,-1))
        xi1_yi = npy.multiply(npy.roll(x,-1),y)
        
        a=0.5*npy.sum(xi_yi1-xi1_yi)# signed area!
#        a=self.Area()
        
        cx=npy.sum(npy.multiply(xi_xi1,(xi_yi1-xi1_yi)))/6./a  
        cy=npy.sum(npy.multiply(yi_yi1,(xi_yi1-xi1_yi)))/6./a  

        return Point2D((cx, cy)) 
    
    def PointBelongs(self, point):
        """
        Ray casting algorithm copied from internet...
        """
        return PolygonPointBelongs(point.vector, [p.vector for p in self.points])

    def SecondMomentArea(self, point):
        Ix, Iy, Ixy = 0, 0, 0
        for pi, pj in zip(self.points,self.points[1:]+[self.points[0]]):
            xi, yi = pi.vector-point.vector
            xj, yj = pj.vector-point.vector
            Ix += (yi**2 + yi*yj + yj**2)*(xi*yj - xj*yi)
            Iy += (xi**2 + xi*xj + xj**2)*(xi*yj - xj*yi)
            Ixy += (xi*yj + 2*xi*yi + 2*xj*yj + xj*yi)*(xi*yj - xj*yi)
        if Ix < 0:
            Ix =- Ix
            Iy =- Iy
            Ixy =- Ixy
        return npy.array([[Ix/12., Ixy/24.], [Ixy/24., Iy/12.]])

    def _LineSegments(self):
        lines=[]
        for p1,p2 in zip(self.points,self.points[1:]+[self.points[0]]):
            lines.append(LineSegment2D(p1,p2))
        return lines
    
    def Rotation(self, center, angle, copy=False):
        if copy:
            return Polygon2D([p.Rotation(center,angle,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center,angle,copy=False)
            
    def Translation(self, offset, copy=False):
        if copy:
            return Polygon2D([p.Translation(offset,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset,copy=False)

    def PointBorderDistance(self, point):
        """
        Compute the distance to the border distance of polygon
        Output is always positive, even if the point belongs to the polygon
        """
        d_min = self.line_segments[0].PointDistance(point)
        for line in self.line_segments[1:]:
            d=line.PointDistance(point)
            if d<d_min:
                d_min=d
        return d_min
                

class Primitive3D:
    def __init__(self, name=''):
        self.name = name        

        
class Vector3D(Vector):
    def __init__(self, vector):
        self.vector=npy.zeros(3)
        self.vector[0] = vector[0]
        self.vector[1] = vector[1]
        self.vector[2] = vector[2]
    
    def __round__(self, ndigits):
        return self.__class__((round(self.vector[0], ndigits),
                               round(self.vector[1], ndigits),
                               round(self.vector[2], ndigits)))
    
    def Dot(self, other_vector):
        u1, u2, u3 = self.vector
        v1, v2, v3 = other_vector.vector
        return u1*v1 + u2*v2 + u3*v3
    
    def Cross(self, other_vector):
        u1, u2, u3 = self.vector
        v1, v2, v3 = other_vector.vector
        return Vector3D((u2*v3 - u3*v2, u3*v1 - u1*v3, u1*v2 - u2*v1))
    
    def Norm(self):
        x,y,z = self.vector
        return (x**2 + y**2 + z**2)**0.5
    
    def Rotation(self, center, axis, angle, copy=True):
        u = axis.vector
        ux = npy.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
        R = math.cos(angle)*npy.eye(3)+math.sin(angle)*ux+(1-math.cos(angle))*npy.tensordot(u,u,axes=0)
        vector2 = npy.dot(R,(self.vector-center.vector))+center.vector
        if copy:
            return Point3D(vector2)
        else:
            self.vector = vector2

    def Translation(self, offset, copy=True):
        vector2 = self.vector+offset
        if copy:
            return Point3D(vector2)
        else:
            self.vector = vector2    
            
    def RandomUnitNormalVector(self):
        """
        Returns a random normal vector
        """
        v = npy.random.random(3)
        
        v = Vector3D(v-npy.dot(v,self.vector)*self.vector/(self.Norm()**2))
        v.vector = v.vector/v.Norm()
        return v

o3D = Vector3D((0, 0, 0))
x3D = Vector3D((1, 0, 0))
y3D = Vector3D((0, 1, 0))
z3D = Vector3D((0, 0, 1))

        
class Point3D(Vector3D):
    def __init__(self, vector, name=''):
        Vector3D.__init__(self, vector)
        self.name=name
        
    def MPLPlot(self, ax):
        ax.scatter(*self.vector)
        
    def PlaneProjection3D(self, plane_origin, x, y):
        z = npy.cross(x.vector,y.vector)
        z = x.Cross(y)
        z /= z.Norm()
        return Point3D(self.vector-npy.dot(self.vector-plane_origin.vector,z)*z)

    def PlaneProjection2D(self, x, y):
        z = npy.cross(x.vector,y.vector)
        z = x.Cross(y)
        z.Normalize()
        p3d = self - self.Dot(z)*z
        u1 = p3d.Dot(x)
        u2 = p3d.Dot(y)
        return Point2D((u1, u2))

        
    def To2D(self, plane_origin, x, y):
        if npy.dot(x.vector, y.vector) != 0:
            raise NotImplementedError
        x2d = npy.dot(self.vector, x.vector) - npy.dot(plane_origin.vector, x.vector)
        y2d = npy.dot(self.vector, y.vector) - npy.dot(plane_origin.vector, y.vector)
        return Point2D((x2d,y2d))
    
    def PointDistance(self, point2):
        return (self-point2).Norm()
        
class Basis3D:
    """
    Defines a 3D basis
    :param origin: origin of the basis
    :param u: first vector of the basis
    :param v: second vector of the basis
    :param w: third vector of the basis
    """
    # TODO: create a Basis and Frame class to mutualize between 2D and 2D
    def __init__(self, u, v, w):
        self.u = u
        self.v = v
        self.w = w
        
    def __repr__(self):
        return '{}: U={}, V={}, W={}'.format(self.__class__.__name__, self.u, self.v, self.w)

        
    def TransfertMatrix(self):
        return npy.array([[self.u[0], self.v[0], self.w[0]],
                          [self.u[1], self.v[1], self.w[1]],
                          [self.u[2], self.v[2], self.w[2]]])
    
    def InverseTransfertMatrix(self):
        # Todo: cache for performance
        return inv(self.TransfertMatrix())

    def NewCoordinates(self, vector):
        return Vector3D(npy.dot(self.InverseTransfertMatrix(), vector.vector))

    def OldCoordinates(self, vector):
        return Vector3D(npy.dot(self.TransfertMatrix(), vector.vector))

        
xyz = Basis3D(x3D, y3D, z3D)

class Frame3D(Basis3D):
    """
    Defines a 3D frame
    :param origin: origin of the basis
    :param u: first vector of the basis
    :param v: second vector of the basis
    :param w: third vector of the basis
    """
    def __init__(self, origin, u, v, w):
        self.origin = origin
        Basis3D.__init__(self, u, v, w)
        
    def __repr__(self):
        return '{}: O= {} U={}, V={}, W={}'.format(self.__class__.__name__, self.origin, self.u, self.v, self.w)
        
    def Basis(self):
        return Basis3D(self.u, self.v, self.w)
    
    def NewCoordinates(self, vector):
        return Basis3D.NewCoordinates(self, vector - self.origin)

    def OldCoordinates(self, vector):
        return Basis3D.OldCoordinates(self, vector) + self.origin

        
oxyz = Frame3D(o3D, x3D, y3D, z3D)    
  
class Line3D(Primitive3D, Line):
    """
    Define an infinite line passing through the 2 points
    """
    def __init__(self, point1, point2, name=''):
        Primitive3D.__init__(self, name)        
        self.points = [point1, point2]    
        
    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        return self.points[0] + (self.points[1]-self.points[0]) * curvilinear_abscissa
        
    def MPLPlot(self, ax):
        # Line segment
        x = [p.vector[0] for p in self.points]
        y = [p.vector[1] for p in self.points]
        z = [p.vector[2] for p in self.points]
        ax.plot(x,y,z, 'ok')
        
        # Drawing 3 times length of segment on each side
        u = self.points[1] - self.points[0]
        x1, y1, z1 = (self.points[0] - 3*u).vector
        x2, y2, z2 = (self.points[1] + 3*u).vector
        ax.plot([x1, x2], [y1, y2], [z1, z2], '-k')
        
    def MinimumDistancePoints(self, other_line):
        """
        Returns the points on this line and the other line that are the closest
        of lines
        """
        u = self.points[1] - self.points[0]
        v = other_line.points[1] - other_line.points[0]
        w = self.points[0] - other_line.points[0]
        a = u.Dot(u)
        b = u.Dot(v)
        c = v.Dot(v)
        d = u.Dot(w)
        e = v.Dot(w)
    
        s = (b*e -c*d) / (a*c - b**2)
        t = (a*e -b*d) / (a*c - b**2)
        p1 = self.points[0] + s*u
        p2 = other_line.points[0] + t*v
        return p1, p2
    
class LineSegment3D(Line3D):
    """
    Define a line segment limited by two points
    """
    def __init__(self, point1, point2, name=''):
        Line3D.__init__(self, point1, point2, name)
        
    def Length(self):
        return self.points[1].PointDistance(self.points[0])

    def PlaneProjection2D(self, x, y):
        return LineSegment2D(self.points[0].PlaneProjection2D(x, y),
                             self.points[1].PlaneProjection2D(x, y))
        
    def MPLPlot(self, ax):
        x=[p.vector[0] for p in self.points]
        y=[p.vector[1] for p in self.points]
        z=[p.vector[2] for p in self.points]
        ax.plot(x,y,z, 'o-k')
        
    def MPLPlot2D(self, x3D, y3D, ax):
        edge2D =  self.PlaneProjection2D(x3D, y3D)
        edge2D.MPLPlot(ax)
        
        
    def FreeCADExport(self, name, ndigits=6):
        x1, y1, z1 = npy.round(1000*self.points[0].vector, ndigits)
        x2, y2, z2 = npy.round(1000*self.points[1].vector, ndigits)
        return '{} = Part.LineSegment(fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,x1,y1,z1,x2,y2,z2)


class Circle3D(Primitive3D):
    def __init__(self, center, radius, normal, name=''):        
        Primitive2D.__init__(self, name)        
        self.center = center
        self.radius = radius
        self.normal = normal
        
    def Length(self):
        return 2* math.pi * self.radius
        

    def FreeCADExport(self,name,ndigits=3):
        xc,yc,zc = npy.round(1000*self.center.vector,ndigits)
        xn,yn,zn = npy.round(self.normal.vector,ndigits)
        return '{} = Part.Circle(fc.Vector({},{},{}),fc.Vector({},{},{}),{})\n'.format(name,xc,yc,zc,xn,yn,zn,1000*self.radius)
        

class Arc3D(Primitive3D):
    """
    An arc is defined by a starting point, an end point and an interior point
    """
    def __init__(self, start, interior, end, name=''):        
        Primitive2D.__init__(self, name)        
        self.start = start
        self.interior = interior
        self.end = end
       
        u1 = (self.interior - self.start)
        u2 = (self.interior - self.end)
        n = u1.Cross(u2)
        v1 = u1.Cross(n)
        v2 = u2.Cross(n)

        p11 = 0.5 * (start + interior)# Mid point of segment s,m
        p12 = p11 + v1
        p21 = 0.5 * (end + interior)# Mid point of segment s,m
        p22 = p21 + v2
        
        l1 = Line3D(p11, p12)
        l2 = Line3D(p21, p22)
        
        c1, c2 = l1.MinimumDistancePoints(l2)
        self.center = c1
        self.radius = (self.center - self.start).Norm()
        self.angle = math.atan((self.center - self.start).Norm() / self.radius)
        
    def Length(self):
        return self.radius * self.angle

    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        u = self.start - self.center
        v = self.end - self.center
        w = u.Cross(v)
        w.Normalize()
        return self.start.Rotation(self.center, w, curvilinear_abscissa/self.radius)
        
    def MPLPlot(self, ax):
        # TODO: there's a bug. Maybe in rotation?
        ax.scatter(*self.center.vector,c='b')
        ax.scatter(*self.start.vector,c='r')
        ax.scatter(*self.end.vector,c='r')
        ax.scatter(*self.interior.vector,c='g')
        x = []
        y = []
        z = []
        l = self.Length()
        for i in range(31):
            p = self.PointAtCurvilinearAbscissa(l*(i)/30)
            x.append(p[0])
            y.append(p[1])
            z.append(p[2])
            
        ax.plot(x, y, z, 'k')
    
    def FreeCADExport(self, name, ndigits=6):
        xs, ys, zs = npy.round(1000*self.start.vector, ndigits)
        xm, ym, zm = npy.round(1000*self.interior.vector, ndigits)
        xe, ye, ze = npy.round(1000*self.end.vector, ndigits)
        return '{} = Part.Arc(fc.Vector({},{},{}),fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,xs,ys,zs,xm,ym,zm,xe,ye,ze)


class CompositePrimitive3D(Primitive3D):
    """
    A collection of simple primitives3D
    """
    def __init__(self, primitives, name=''):
        Primitive3D.__init__(self, name)        
        
        basis_primitives=[]
        for primitive in primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.basis_primitives)
            else:
                basis_primitives.append(primitive)
                
        self.basis_primitives = basis_primitives
        
    def UpdateBasisPrimitives(self):
        # TODO: This is a copy/paste from CompositePrimitive2D, in the future make a Common abstract class
        basis_primitives=[]
        for primitive in self.primitives:
            if hasattr(primitive, 'basis_primitives'):
                basis_primitives.extend(primitive.basis_primitives)
            else:
                basis_primitives.append(primitive)
                
        self.basis_primitives = basis_primitives
            
    def MPLPlot(self, ax = None):
        if ax is None:
            fig = plt.figure()
#            ax = fig.add_subplot(111, projection='3d', adjustable='box')
            ax = Axes3D(fig)
        else:
            fig = None

        for primitive in self.basis_primitives:
            primitive.MPLPlot(ax)
        
        ax.set_aspect('equal')
        
        return fig, ax
    
class Wire3D(CompositePrimitive3D):
    """
    A collection of simple primitives, following each other making a wire
    """
    def __init__(self, primitives, name=''):
        CompositePrimitive2D.__init__(self, primitives, name)        
    
    def Length(self):
        length = 0.
        for primitive in self.basis_primitives:
            length += primitive.Length()
        return length
    
    def PointAtCurvilinearAbscissa(self, curvilinear_abscissa):
        length = 0.
        for primitive in self.basis_primitives:
            primitive_length = primitive.Length()
            if length + primitive_length > curvilinear_abscissa:
                return primitive.PointAtCurvilinearAbscissa(curvilinear_abscissa - length)
            length += primitive_length
        # Outside of length
        raise ValueError

    # TODO: method to check if it is a wire
    
    
    def FreeCADExport(self, ip):
        name='primitive'+str(ip)
        
        s = 'E = []\n'
        for ip, primitive in enumerate(self.basis_primitives):
            s += primitive.FreeCADExport('L{}'.format(ip))
            s += 'E.append(Part.Edge(L{}))\n'.format(ip)
        s += '{} = Part.Wire(E[:])\n'.format(name)

        return s
        
class Contour3D(Wire3D):
    """
    A collection of 3D primitives forming a closed wire3D
    """
    def __init__(self, primitives, name=''):    
        primitives2=[]
        for primitive in primitives:
            try:
                primitives2.extend(primitive.primitives)
            except AttributeError:
                primitives2.append(primitive)

        CompositePrimitive3D.__init__(self,primitives2, name)
        

class VolumeModel:
    """
    :param groups: A list of two element tuple. The first element is a string naming the group and the second element is a list of primitives of the group
    """
    def __init__(self, groups, name=''):
        self.groups = groups
        self.name=name
        
    def Volume(self):
        volume=0
        for group_name, primitives_group in self.groups:
            for primitive in primitives_group:
                volume+=primitive.Volume()
        return volume
    
    def MPLPlot(self):
        """
        Matplotlib plot of model.
        To use for debug.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d', adjustable='box')
#        ax.set_aspect('equal')
        for name, primitive_group in self.groups:
            for primitive in primitive_group:
                primitive.MPLPlot(ax)
        ax.set_aspect('equal')
        return fig, ax
    
    def FreeCADScript(self, fcstd_filepath,
                      path_lib_freecad='/usr/lib/freecad/lib',
                      export_types=['fcstd'],
                      save_to = ''):
        """
        Generate python a FreeCAD definition of model 
        :param fcstd_filename: a filename without extension to give the name at the fcstd part written in python code
        :type fcstd_filename:str
        """
        fcstd_filepath = os.path.abspath(fcstd_filepath)
        fcstd_filepath = fcstd_filepath.replace('\\','\\\\')
        path_lib_freecad = path_lib_freecad.replace('\\','\\\\')
        
        s=''
        if path_lib_freecad != '':
            s+="import sys\nsys.path.append('"+path_lib_freecad+"')\n"

        s+="import math\nimport FreeCAD as fc\nimport Part\n\ndoc=fc.newDocument('doc')\n\n"
        
        for ig, (group_name, primitives_group) in enumerate(self.groups):
            if group_name == '':
                group_name = 'Group_{}'.format(ig)
            else:
                group_name = 'Group_{}_{}'.format(ig, group_name)
            s += "part = doc.addObject('App::Part','{}')\n".format(group_name)
            for ip, primitive in enumerate(primitives_group):
                sp = primitive.FreeCADExport(ip)
                if sp != '':
                    s += (sp+'\n')
                    if primitive.name != '':
                        primitive_name = 'primitive_{}_{}_{}'.format(ig, ip, primitive.name)
                    else:
                        primitive_name = 'primitive_{}_{}'.format(ig, ip)
                    s += 'shapeobj = doc.addObject("Part::Feature","{}")\n'.format(primitive_name)
                    s += "shapeobj.Shape = primitive{}\n".format(ip)
                    s += 'part.addObject(shapeobj)\n'.format(ip, primitive.name)
#                    
        s+='doc.recompute()\n'
        if 'fcstd' in export_types:
            s+="doc.saveAs('"+fcstd_filepath+".fcstd')\n\n"
        if 'stl' in export_types:
            s+="import Mesh\nMesh.export(doc.Objects,'{}.stl')\n".format(fcstd_filepath)
        if 'step' in export_types:
            s+="Part.export(doc.Objects,'{}.step')\n".format(fcstd_filepath)
                

        if save_to != '':
            with open(os.path.abspath(save_to),'w') as file:
                file.write(s)
        return s
            
    
    
    def FreeCADExport(self,fcstd_filepath,
                      python_path='python',
                      path_lib_freecad='/usr/lib/freecad/lib', 
                      export_types=['fcstd']):
        """
        Export model to .fcstd FreeCAD standard
        
        :param python_path: path of python binded to freecad
        
            * on windows: something like C:\\\\Program Files\\\\FreeCAD X.XX\\\\bin\\\\python
            * on linux: python if installed by a dstribution package
        :param filepath: path of fcstd file (without extension)
        :param path_lib_freecad: FreeCAD.so lib path (/usr/lib/freecad/lib in general)

        """
        fcstd_filepath=os.path.abspath(fcstd_filepath)
        s=self.FreeCADScript(fcstd_filepath,
                             path_lib_freecad = path_lib_freecad,
                             export_types = export_types)
        with tempfile.NamedTemporaryFile(suffix=".py",delete=False) as f:
            f.write(bytes(s,'utf8'))

        arg=f.name
        output=subprocess.call([python_path,arg])

        f.close()
        os.remove(f.name)
        return output
        

    
    
    def BabylonScript(self):

        env = Environment(loader=PackageLoader('volmdlr', 'templates'),
                          autoescape=select_autoescape(['html', 'xml']))
        
        template = env.get_template('babylon.html')
        
        center,max_length=self.ModelCaracteristicLengths()
        
        primitives_strings=[]
        for primitive in self.primitives:
            try:
                primitives_strings.append(primitive.Babylon())
            except AttributeError:
                pass
        return template.render(name=self.name,center=tuple(center),length=2*max_length,
                               primitives_strings=primitives_strings)
    
    def BabylonShow(self,page='vm_babylonjs'):
        page+='.html'
        with open(page,'w') as file:
            file.write(self.BabylonScript())
        
        webbrowser.open('file://' + os.path.realpath(page))
        
    def ModelCaracteristicLengths(self):
        min_vect = self.primitives[0].position
        max_vect = self.primitives[0].position
        center = self.primitives[0].position
        n=1
        for primitive in self.primitives[1:]:
            try:
                for i,(xmin,xmax,xi) in enumerate(zip(min_vect, max_vect, primitive.position)):

                    if xi<xmin:
                        min_vect[i]=xi

                    if xi>xmax:
                        max_vect[i]=xi
                center += primitive.position
                n+=1
            except AttributeError:
                pass
                
        center=center/n

        max_length = (min_vect-max_vect).Norm()

        return center,max_length

