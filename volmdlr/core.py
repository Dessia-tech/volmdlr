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
#import vmcy
from .vmcy import PolygonPointBelongs

from scipy.linalg import norm,solve,LinAlgError

import volmdlr.geometry as geometry

from jinja2 import Environment, PackageLoader, select_autoescape

import webbrowser
import os

import tempfile
import subprocess


class Vector2D:
    def __init__(self,vector):
#        print(vector,npy.array(vector))
#        self.vector=npy.array([vector[0],vector[1]])
        self.vector=npy.zeros(2)
        self.vector[0]=vector[0]        
        self.vector[1]=vector[1]
        
    def To3D(self,plane_origin,x1,x2):
        x,y=self.vector
        return Point3D(plane_origin.vector+x1.vector*x+x2.vector*y)
    
    def __add__(self,point2d):
        return Vector2D(self.vector+point2d.vector)
    
    def __radd__(self,point2d):
        return self+point2d

    def __sub__(self,point2d):
        return Vector2D(self.vector-point2d.vector)
    
    def __rsub__(self,point2d):
        return self-point2d    
    
    def __mul__(self,value):
        return Vector2D(self.vector*value)
    
    def __rmul__(self,value):
        return self*value

    def __truediv__(self,value):
        return Vector2D(self.vector/value)
    
    def __rtruediv__(self,value):
        return self/value
    
    def Rotation(self,center,angle,copy=True):
        vector2=npy.dot(npy.array([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]]),(self.vector-center.vector))+center.vector
        if copy:
            return Point2D(vector2)
        else:
            self.vector=vector2

    def Translation(self,offset,copy=True):
        vector2=self.vector+offset
        if copy:
            return Point2D(vector2)
        else:
            self.vector=vector2
        

class Point2D(Vector2D):
    def __init__(self,vector,name=''):
        Vector2D.__init__(self,vector)
        self.name=name
        
    def __add__(self,point2d):
        return Point2D(self.vector+point2d.vector)
    
    def __sub__(self,point2d):
        return Point2D(self.vector-point2d.vector)
        
    def __mul__(self,value):
        return Point2D(self.vector*value)
    
    def __truediv__(self,value):
        return Point2D(self.vector/value)
    

    def MPLPlot(self,style='ob'):
        x1=self.vector
        plt.plot([x1[0]],[x1[1]],style)        
        return []
    
    def PointDistance(self,point2):
        return norm(self.vector-point2.vector)

    #def Distance(self,point2):
    #    return norm(self.vector-point2.vector)

    @classmethod
    def LinesIntersection(cls,line1,line2,curvilinear_abscissa=False):
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
        
class Primitive2D:
    def __init__(self,name=''):
        self.name=name
        
class CompositePrimitive2D(Primitive2D):
    """
    A collection of simple primitives
    """
    def __init__(self,primitives,name=''):
        Primitive2D.__init__(self,name)        
        self.primitives=primitives
        
    def Rotation(self,center,angle,copy=False):
        if copy:
            return self.__class__([p.Rotation(center,angle,copy=True) for p in self.primitives])
        else:
            for p in self.primitives:
                p.Rotation(center,angle,copy=False)
            
    def Translation(self,offset,copy=False):
        if copy:
            return self.__class__([p.Translation(offset,copy=True) for p in self.primitives])
        else:
            for p in self.primitives:
                p.Translation(offset,copy=False)
    
        
    def MPLPlot(self):
#        lines,arcs=self.Export2D(point,x1,x2)
        fig, ax = plt.subplots()
        ax.set_aspect('equal')
        ps=[]
        for element in self.primitives:
            ps.extend(element.MPLPlot())
#        p = PatchCollection(ps)
#        print(arcs)
        for p in ps:
            ax.add_patch(p)
        
        ax.margins(0.1)
        plt.show() 
        
class Contour2D(CompositePrimitive2D):
    """
    A collection of primitives intended to be closed
    """
    def __init__(self,primitives,name=''):
        
        primitives2=[]
        for primitive in primitives:
            try:
                primitives2.extend(primitive.primitives)
            except AttributeError:
                primitives2.append(primitive)

        CompositePrimitive2D.__init__(self,primitives2,name)
        
    def To3D(self,plane_origin,x,y):
        return [primitive.To3D(plane_origin,x,y) for primitive in self.primitives]
    
    def Area(self):
        if len(self.primitives)==1:
            return self.primitives[0].Area()
        arcs=[]
        points_polygon=[]
        for primitive in self.primitives:
#            print(primitive.__class__.__name__)
            if primitive.__class__.__name__=='Line2D':
                points_polygon.extend(primitive.points)
            elif primitive.__class__.__name__=='Arc2D':
                points_polygon.append(primitive.center)
                arcs.append(primitive)
        polygon=Polygon2D(points_polygon)
        A=polygon.Area()
#        print('A',A)
#        print('arcs',arcs)
        for arc in arcs:
#            print(arc,arc.Area())
            if polygon.PointBelongs(arc.middle):
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
            if primitive.__class__.__name__=='Line2D':
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
            if polygon.BelongsTo(arc.middle):
                A-=arc.SecondMomentArea(point)
            else:
                A+=arc.SecondMomentArea(point)
        return A

class Mesh2D:
    def __init__(self,contours,points_densities,default_density):
        self.contours=contours
        self.points_densities=points_densities
        self.default_density=default_density
        
    def GeoScript(self,filepath=''):
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

class Line2D(Primitive2D):
    def __init__(self,point1,point2,name=''):
        Primitive2D.__init__(self,name)        
        self.points=[point1,point2]
        
    def _get_geo_points(self):
        return self.points

    geo_points=property(_get_geo_points)      
    
    def PointProjection(self,point,curvilinear_abscissa=False):
        p1,p2=self.points
        v=p1.vector
        w=p2.vector
        p=point.vector
        t = npy.dot(p - v, w - v) / norm(w-v)**2
#        print('t', t)
        projection = Point2D(v + t * (w - v))# Projection falls on the segment
        if curvilinear_abscissa:
            return projection,t
        return projection
    
    def PointSegmentDistance(self,point):
        """
        Computes the distance of a point to segment of line 
        """
        p1,p2=self.points
        v=p1.vector
        w=p2.vector
        p=point.vector
#        print('point vector',point.vector)
        t = max(0, min(1, npy.dot(p - v, w - v) / norm(w-v)**2))
#        print('t', t)
        projection = v + t * (w - v)# Projection falls on the segment
#        print(p,projection)
        return norm(p-projection);
    
    def PointDistance(self,point):
        """
        Computes the distance of a point to line 
        """
        p1,p2=self.points
        v=p1.vector
        w=p2.vector
        p=point.vector
#        print('point vector',point.vector)
        t = npy.dot(p - v, w - v) / norm(w-v)**2
#        print('t', t)
        projection = v + t * (w - v)# Projection falls on the segment
#        print(p,projection)
        return norm(p-projection);


    def MPLPlot(self):
        p1,p2=self.points
        x1=p1.vector
        x2=p2.vector
        plt.plot([x1[0],x2[0]],[x1[1],x2[1]],'black')        
        return []
    
    def To3D(self,plane_origin,x1,x2):
        p3D=[p.To3D(plane_origin,x1,x2) for p in self.points]
        return Line3D(*p3D,self.name)
    
    def Rotation(self,center,angle,copy=False):
        if copy:
            return Line2D(*[p.Rotation(center,angle,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Rotation(center,angle,copy=False)
            
    def Translation(self,offset,copy=False):
        if copy:
            return Line2D(*[p.Translation(offset,copy=True) for p in self.points])
        else:
            for p in self.points:
                p.Translation(offset,copy=False)
                
    def GeoScript(self,primitive_index,points_indices):
        s='Line({}) = {{{}, {}}};\n'.format(primitive_index,*points_indices)
        return s,primitive_index+1
                
                
class Arc2D(Primitive2D):
    def __init__(self,start,middle,end,name=''):        
        Primitive2D.__init__(self,name)        
        self.middle=middle
        self.start=start
        self.end=end
        xm,ym=middle.vector
        xe,ye=end.vector
        xs,ys=start.vector
        A=npy.array([[2*(xs-xm),2*(ys-ym)],
                   [2*(xs-xe),2*(ys-ye)]])
        b=-npy.array([xm**2+ym**2-xs**2-ys**2,xe**2+ye**2-xs**2-ys**2])
        self.center=Point2D(solve(A,b))
        r1=self.start.vector-self.center.vector
        r2=self.end.vector-self.center.vector
        rc=self.middle.vector-self.center.vector
        
        self.radius=norm(r1)
        angle1=npy.arctan2(r1[1],r1[0])
        angle2=npy.arctan2(r2[1],r2[0])
#        angle_min=min(angle1,angle2)
#        angle_max=max(angle1,angle2)
        anglem=npy.arctan2(rc[1],rc[0])
        order=[y for x,y in sorted(zip([angle1,anglem,angle2],[0,1,2]))]
        order=order*2
        i=order.index(0)
        if order[i+1]==1:
            self.angle1=angle1
            self.angle2=angle2
        else:
            self.angle1=angle2
            self.angle2=angle1
            
        
    def _get_geo_points(self):
        return [self.start,self.center,self.end]

    geo_points=property(_get_geo_points)        
    
    def GeoScript(self,primitive_index,points_indices):
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
        u/=norm(u)
        alpha=abs(self.angle1-self.angle2)
        return Point2D(self.center.vector+4/(3*alpha)*self.radius*math.sin(alpha*0.5)*u)
        
    def MPLPlot(self):
        pc=self.center.vector
        return [Arc(pc,2*self.radius,2*self.radius,angle=0,theta1=self.angle1*0.5/math.pi*360,theta2=self.angle2*0.5/math.pi*360,color='black')]

    def To3D(self,plane_origin,x,y):
        pc=self.center.To3D(plane_origin,x,y)
        pe=Point2D(self.center.vector+self.radius*npy.array((math.cos(self.angle1),math.sin(self.angle1))))
        pe3=pe.To3D(plane_origin,x,y)
        ps=Point2D(self.center.vector+self.radius*npy.array((math.cos(self.angle2),math.sin(self.angle2))))
        ps3=ps.To3D(plane_origin,x,y)
        return Arc3D(pe3,pc,ps3,self.name)
    
    def Rotation(self,center,angle,copy=False):
        if copy:
            return Arc2D(*[p.Rotation(center,angle,copy=True) for p in [self.start,self.middle,self.end]])
        else:
            self.__init__(*[p.Rotation(center,angle,copy=True) for p in [self.start,self.middle,self.end]])
            
    def Translation(self,offset,copy=False):
        if copy:
            return Arc2D(*[p.Translation(offset,copy=True) for p in [self.start,self.middle,self.end]])
        else:
            self.__init__(*[p.Translation(offset,copy=True) for p in [self.start,self.middle,self.end]])

    def SecondMomentArea(self,point):
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
        return geometry.Huygens2D(Ic,self.Area(),self.center,point)



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

    geo_points=property(_get_geo_points)        
            
    def GeoScript(self,primitive_index,points_indices):
        s='Circle({}) = {{{}, {}, {}}};\n'.format(primitive_index,*points_indices)
        return s,primitive_index+1
    
    def MPLPlot(self):
        pc=self.center.vector
        return [Arc(pc,2*self.radius,2*self.radius,angle=0,theta1=0,theta2=360,color='black')]

    def To3D(self,plane_origin,x,y):
        normal=Vector3D(npy.cross(x.vector,y.vector))
        pc=self.center.To3D(plane_origin,x,y)
        return Circle3D(pc,self.radius,normal,self.name)

    def Rotation(self,center,angle,copy=False):
        if copy:
            return Circle2D(self.center.vector.Rotation(center,angle,copy=True),self.radius)
        else:
            self.center.Rotation(center,angle,copy=False) 
            self.utd_geo_points=False
            
    def Translation(self,offset,copy=False):
        if copy:
            return Circle2D(self.center.vector.Translation(offset,copy=True),self.radius)
        else:
            self.center.Translation(offset,copy=False)
            self.utd_geo_points=False

    def Area(self):
        return math.pi*self.radius**2

    def SecondMomentArea(self,point):
        """ 
        Second moment area of part of disk
        """
        I=math.pi*self.radius**4/4
        Ic=npy.array([[I,0],[0,I]])
        return geometry.Huygens2D(Ic,self.Area(),self.center,point)
    
    def CenterOfMass(self):
        return self.center

class Polygon2D(CompositePrimitive2D):
    def __init__(self,points,name=''):     
        self.points=points
        primitives=[]
        for p1,p2 in zip(points,points[1:]+[points[0]]):
            primitives.append(Line2D(p1,p2))
            
            
        CompositePrimitive2D.__init__(self,primitives,name)
        
        
    def Area(self):
        
        x=[point.vector[0]for point in self.points]
        y=[point.vector[1]for point in self.points]

        return 0.5*npy.abs(npy.dot(x,npy.roll(y,1))-npy.dot(y,npy.roll(x,1)))
    
    def CenterOfMass(self):
        
        x=[point.vector[0] for point in self.points]
        y=[point.vector[1] for point in self.points]

        
        xi_xi1=x+npy.roll(x,-1)
        yi_yi1=y+npy.roll(y,-1)
        xi_yi1=npy.multiply(x,npy.roll(y,-1))
        xi1_yi=npy.multiply(npy.roll(x,-1),y)
        
        a=0.5*npy.sum(xi_yi1-xi1_yi)# signed area!
#        a=self.Area()
        
        cx=npy.sum(npy.multiply(xi_xi1,(xi_yi1-xi1_yi)))/6./a  
        cy=npy.sum(npy.multiply(yi_yi1,(xi_yi1-xi1_yi)))/6./a  

        return Point2D((cx,cy)) 
    
    def PointBelongs(self,point):
        """
        Ray casting algorithm copied from internet...
        """
        return PolygonPointBelongs(point.vector,[p.vector for p in self.points])

    def SecondMomentArea(self,point):
        Ix,Iy,Ixy=0,0,0
        for pi,pj in zip(self.points,self.points[1:]+[self.points[0]]):
            xi,yi=pi.vector-point.vector
            xj,yj=pj.vector-point.vector
            Ix+=(yi**2+yi*yj+yj**2)*(xi*yj-xj*yi)
            Iy+=(xi**2+xi*xj+xj**2)*(xi*yj-xj*yi)
            Ixy+=(xi*yj+2*xi*yi+2*xj*yj+xj*yi)*(xi*yj-xj*yi)
        if Ix<0:
            Ix=-Ix
            Iy=-Iy
            Ixy=-Ixy
        return npy.array([[Ix,Ixy],[Ixy,Iy]])

    def Lines(self):
        lines=[]
        for p1,p2 in zip(self.points,self.points[1:]+[self.points[0]]):
            lines.append(Line2D(p1,p2))
        return lines
    
    def MPLPlot(self):
        x=[]
        y=[]
        for p in self.points+[self.points[0]]:
            x.append(p.vector[0])
            y.append(p.vector[1])
        plt.plot(x,y,label=self.name)
            
#        for line in self.Lines():
#            line.MPLPlot()
        return []

    def PointDistance(self,point):
        lines=self.Lines()
        d_min=lines[0].PointSegmentDistance(point)
#        print('d: ',d_min,0)
        for line in lines[1:]:
            d=line.PointSegmentDistance(point)
#            print('d: ',d)
            if d<d_min:
                d_min=d
        return d_min
                

class Primitive3D:
    def __init__(self,name=''):
        self.name=name        

        
class Vector3D:
    def __init__(self,vector):
#        self.vector=npy.array(vector)
        self.vector=npy.zeros(3)
        self.vector[0]=vector[0]        
        self.vector[1]=vector[1]
        self.vector[2]=vector[2]
    
    def __add__(self,point2d):
        return self.__class__(self.vector+point2d.vector)
    
    def __radd__(self,point2d):
        return self+point2d

    def __sub__(self,point2d):
        return Vector3D(self.vector-point2d.vector)
    
    def __rsub__(self,point2d):
        return self-point2d    
    
    def __mul__(self,value):
        return Vector3D(self.vector*value)
    
    def __rmul__(self,value):
        return self*value

    def __div__(self,value):
        return Vector3D(self.vector/value)
    
    def __rdiv__(self,value):
        return self/value
        
    def Rotation(self,center,axis,angle,copy=True):
        u=axis.vector
        ux=npy.array([[0,-u[2],u[1]],[u[2],0,-u[0]],[-u[1],u[0],0]])
        R=math.cos(angle)*npy.eye(3)+math.sin(angle)*ux+(1-math.cos(angle))*npy.tensordot(u,u,axes=0)
        vector2=npy.dot(R,(self.vector-center.vector))+center.vector
        if copy:
            return Point3D(vector2)
        else:
            self.vector=vector2

    def Translation(self,offset,copy=True):
        vector2=self.vector+offset
        if copy:
            return Point3D(vector2)
        else:
            self.vector=vector2    

        
class Point3D(Vector3D):
    def __init__(self,vector,name=''):
        Vector3D.__init__(self,vector)
        self.name=name
        
    def MPLPlot(self,ax):
        ax.scatter(*self.vector)
        
    def PlaneProjection(self,plane_origin,x,y):
        z=npy.cross(x.vector,y.vector)
        z/=norm(z)
        return Point3D(self.vector-npy.dot(self.vector-plane_origin.vector,z)*z)
        
    def To2D(self,plane_origin,x,y):
        if npy.dot(x.vector,y.vector)!=0:
            raise NotImplemented
        x2d=npy.dot(self.vector,x.vector)-npy.dot(plane_origin.vector,x.vector)
        y2d=npy.dot(self.vector,y.vector)-npy.dot(plane_origin.vector,y.vector)
        return Point2D((x2d,y2d))
    
    def PointDistance(self,point2):
        return norm(self.vector-point2.vector)
        
class Line3D(Primitive3D):
    def __init__(self,point1,point2,name=''):
        Primitive3D.__init__(self,name)        
        self.points=[point1,point2]
        
    def MPLPlot(self,ax):
        x=[p.vector[0] for p in self.points]
        y=[p.vector[1] for p in self.points]
        z=[p.vector[2] for p in self.points]
        ax.plot(x,y,z)
        
        
    def FreeCADExport(self,name,ndigits=3):
        x1,y1,z1=npy.round(1000*self.points[0].vector,ndigits)
        x2,y2,z2=npy.round(1000*self.points[1].vector,ndigits)
        return '{}=Part.Line(fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,x1,y1,z1,x2,y2,z2)


class Circle3D(Primitive3D):
    def __init__(self,center,radius,normal,name):        
        Primitive2D.__init__(self,name)        
        self.center=center
        self.radius=radius
        self.normal=normal

    def FreeCADExport(self,name,ndigits=3):
        xc,yc,zc=npy.round(1000*self.center.vector,ndigits)
        xn,yn,zn=npy.round(self.normal.vector,ndigits)
        return '{}=Part.Circle(fc.Vector({},{},{}),fc.Vector({},{},{}),{})\n'.format(name,xc,yc,zc,xn,yn,zn,1000*self.radius)
        

class Arc3D(Primitive3D):
    def __init__(self,start,center,end,name):        
        Primitive2D.__init__(self,name)        
        self.center=center
        self.start=start
        self.end=end
        r=norm(start.vector-center.vector)
        n=(start.vector-2*center.vector+end.vector)
        n=n/norm(n)
        self.middle=Point3D(r*n+center.vector)

    def MPLPlot(self,ax):
        ax.scatter(*self.center.vector,c='b')
        ax.scatter(*self.start.vector,c='r')
        ax.scatter(*self.end.vector,c='r')
        ax.scatter(*self.middle.vector,c='g')
        
        
    
    def FreeCADExport(self,name,ndigits=3):
        xs,ys,zs=npy.round(1000*self.start.vector,ndigits)
        xm,ym,zm=npy.round(1000*self.middle.vector,ndigits)
        xe,ye,ze=npy.round(1000*self.end.vector,ndigits)
        return '{}=Part.Arc(fc.Vector({},{},{}),fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,xs,ys,zs,xm,ym,zm,xe,ye,ze)

        

class VolumeModel:
    def __init__(self,primitives,name=''):
        self.primitives=primitives
        self.name=name
        
    def Volume(self):
        volume=0
        for primitive in self.primitives:
            volume+=primitive.Volume()
        return volume
    
    def MPLPlot(self):
        """
        Matplotlib plot of model.
        To use for debug.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d',adjustable='box')
        ax.set_aspect('equal')
        for primitive in self.primitives:
            primitive.MPLPlot(ax)
#        ax.set_aspect('equal')
    
    def FreeCADScript(self,fcstd_filepath,path_lib_freecad='',py_filepath='',export_types=['fcstd']):
        """
        Generate python a FreeCAD definition of model 
        :param fcstd_filename: a filename without extension to give the name at the fcstd part written in python code
        :param py_filepath: python file to be exported. If none give, only return in string
        :type fcstd_filename:str
        :type py_filepath:str
        """
        fcstd_filepath=os.path.abspath(fcstd_filepath)
        s=''
        if path_lib_freecad!='':
            s+="import sys\nsys.path.append('"+path_lib_freecad+"')\n"

        s+="import math\nimport FreeCAD as fc\nimport Part\n\ndoc=fc.newDocument('doc')\n\n"
        
        for ip,primitive in enumerate(self.primitives):
            s+=(primitive.FreeCADExport(ip)+'\n')
            s+='shapeobj = doc.addObject("Part::Feature","p{} {}")\n'.format(ip,primitive.name)
            s+="shapeobj.Shape = primitive{}\n".format(ip)
        s+='doc.recompute()\n'
        if 'fcstd' in export_types:
            s+="doc.saveAs('"+fcstd_filepath+".fcstd')\n\n"
        if 'stl' in export_types:
#            p_list=str(['primitive{}'.format(i) for i in range(len(self.primitives))]).replace("'","")
            s+="import Mesh\nMesh.export(doc.Objects,'{}.stl')\n".format(fcstd_filepath)
                

        if py_filepath!='':
            with open(py_filepath,'w') as fichier:
                fichier.write(s)
        return s
            
    
    
    def FreeCADExport(self,python_path,fcstd_filepath,path_lib_freecad='',export_types=['fcstd']):
        """
        Export model to .fcstd FreeCAD standard
        
        :param python_path: path of python binded to freecad
             - on windows: 'something like C:\Program Files\FreeCAD X.XX\bin\python'
             - on linux: python (in general)
        :param filepath: path of fcstd file (without extension)
        :param path_lib_freecad: FreeCAD.so lib path (/usr/lib/freecad/lib in general)

        """
        fcstd_filepath=os.path.abspath(fcstd_filepath)
        s=self.FreeCADScript(fcstd_filepath,path_lib_freecad,'',export_types)
#        print(s)
        with tempfile.NamedTemporaryFile(suffix=".py",delete=False) as f:
#        with open(f,'w') as fw:
            f.write(bytes(s,'utf8'))
#            print(f.name)
#            import time
#            time.sleep(20)


#            f.read()
#            print(python_path,arg)
        arg=f.name
        output=subprocess.call([python_path,arg])
#        print(f.name)
        f.close()
#        p = subprocess.Popen('{} {}'.format(python_path,arg), shell=True, stdin=subprocess.PIPE, 
#                  stdout=subprocess.PIPE, stderr=subprocess.PIPE, close_fds=True)
#        output = p.stdout.read()
#        output += p.stderr.read()
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
#        print(self.BabylonScript())
        with open(page,'w') as file:
            file.write(self.BabylonScript())
        
        webbrowser.open('file://' + os.path.realpath(page))
        
    def ModelCaracteristicLengths(self):
        min_vect=npy.array(self.primitives[0].position)
        max_vect=npy.array(self.primitives[0].position)
        center=npy.array(self.primitives[0].position)
        n=1
        for primitive in self.primitives[1:]:
            try:
                for i,(xmin,xmax,xi) in enumerate(zip(min_vect,max_vect,primitive.position)):
    #                print(i,xmin,xmax,xi)
                    if xi<xmin:
                        min_vect[i]=xi
    #                    print('min',min_vect)
                    if xi>xmax:
                        max_vect[i]=xi
    #                    print('max',max_vect)
    #            print(linkage,linkage.position)
                center+=primitive.position
                n+=1
            except AttributeError:
                pass
                
        center=center/n
#        print(min_vect,max_vect)
        max_length=norm(min_vect-max_vect)
#        print(center,max_length)
        return center,max_length

