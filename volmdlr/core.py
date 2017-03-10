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

from scipy.linalg import norm,solve,LinAlgError


class Vector2D:
    def __init__(self,vector):
        self.vector=npy.array(vector)
        
    def To3D(self,plane_origin,x1,x2):
        x,y=self.vector
        return Point3D(plane_origin.vector+x1.vector*x+x2.vector*y)
    
    def __add__(self,point2d):
        return Vector2D(self.vector+point2d)
    
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

    def __div__(self,value):
        return Vector2D(self.vector/value)
    
    def __rdiv__(self,value):
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

    def MPLPlot(self):
        x1=self.vector
        plt.plot([x1[0]],[x1[1]],'ob')        
        return []
    
    def Distance(self,point2):
        return norm(self.vector-point2.vector)

    @classmethod
    def LinesIntersection(cls,line1,line2):
        p11=line1.points[0].vector
        p12=line1.points[1].vector
        p21=line2.points[0].vector
        p22=line2.points[1].vector
        A=npy.array([[p11[0]-p12[0],p22[0]-p21[0]],[p11[1]-p12[1],p22[1]-p21[1]]])
        x=npy.array([p21[0]-p11[0],p21[1]-p11[1]])
        try:
            t=solve(A,x)
            return cls(p11-t[0]*(p12-p11))
        except LinAlgError:
            return None
        
    @classmethod
    def MiddlePoint(cls,point1,point2):
        p1=point1.vector
        p2=point2.vector
        return cls((p1+p2)*0.5)
        
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
        
    def MPLPlot(self):
#        lines,arcs=self.Export2D(point,x1,x2)
        fig, ax = plt.subplots()
        ps=[]
        for element in self.primitives:
            ps.extend(element.MPLPlot())
#        p = PatchCollection(ps)
#        print(arcs)
        for p in ps:
            ax.add_patch(p)
        ax.set_aspect('equal')
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
        arcs=[]
        points_polygon=[]
        for primitive in self.primitives:
            if primitive.__class__.__name__=='Line2D':
                points_polygon.extend(primitive.points)
            elif primitive.__class__.__name__=='Arc2D':
                points_polygon.append(primitive.center)
#                arcs.append(primitive)
        polygon=Polygon2D(points_polygon)
        A=polygon.Area()
        for arc in arcs:
            if polygon.BelongsTo(arc.center):
                A-=arc.Area()
            else:
                A+=arc.Area()
                
        return A 


class Line2D(Primitive2D):
    def __init__(self,point1,point2,name=''):
        Primitive2D.__init__(self,name)        
        self.points=[point1,point2]

    def MPLPlot(self):
        p1,p2=self.points
        x1=p1.vector
        x2=p2.vector
        plt.plot([x1[0],x2[0]],[x1[1],x2[1]],'black')        
        return []
    
    def To3D(self,plane_origin,x1,x2):
        p3D=[p.To3D(plane_origin,x1,x2) for p in self.points]
        return Line3D(*p3D,self.name)
    
    def Rotation(self,center,angle,copy=True):
        if copy:
            return Line2D(*[p.vector.Rotation(center,angle,copy=True) for p in self.points])
        else:
            self.points=[p.vector.Rotation(center,angle,copy=True) for p in self.points]
            
    def Translation(self,offset,copy=True):
        if copy:
            return Line2D(*[p.vector.Translation(offset,copy=True) for p in self.points])
        else:
            self.points=[p.vector.Translation(offset,copy=True) for p in self.points]
    
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
            
    def Area(self):
        if self.angle2<self.angle1:
            angle=self.angle2+2*math.pi-self.angle1
        else:
            angle=self.angle2-self.angle1
        return self.radius**2*angle/2
            
        
    def MPLPlot(self):
        pc=self.center.vector
        return [Arc(pc,2*self.radius,2*self.radius,angle=0,theta1=self.angle1*0.5/math.pi*360,theta2=self.angle2*0.5/math.pi*360,color='black')]

    def To3D(self,plane_origin,x,y):
        pc=self.center.To3D(plane_origin,x,y)
        pe=self.center+self.radius*npy.array((math.cos(self.angle1),math.sin(self.angle1)))
        pe3=pe.To3D(plane_origin,x,y)
        ps=self.center+self.radius*npy.array((math.cos(self.angle2),math.sin(self.angle2)))
        ps3=ps.To3D(plane_origin,x,y)
        return Arc3D(pe3,pc,ps3,self.name)
    
    def Rotation(self,center,angle,copy=True):
        if copy:
            return Arc2D(*[p.Rotation(center,angle,copy=True) for p in [self.start,self.middle,self.end]])
        else:
            self.__init__(*[p.Rotation(center,angle,copy=True) for p in [self.start,self.middle,self.end]])
            
    def Translation(self,offset,copy=True):
        if copy:
            return Arc2D(*[p.Translation(offset,copy=True) for p in [self.start,self.middle,self.end]])
        else:
            self.__init__(*[p.Translation(offset,copy=True) for p in [self.start,self.middle,self.end]])

class Circle2D(Primitive2D):
    def __init__(self,center,radius,name=''):        
        Primitive2D.__init__(self,name)        
        self.center=center
        self.radius=radius
        
    def MPLPlot(self):
        pc=self.center.vector
        return [Arc(pc,2*self.radius,2*self.radius,angle=0,theta1=0,theta2=360,color='black')]

    def To3D(self,plane_origin,x,y):
        normal=Vector3D(npy.cross(x.vector,y.vector))
        pc=self.center.To3D(plane_origin,x,y)
        print(normal,pc)
        return Circle3D(pc,self.radius,normal,self.name)

    def Rotation(self,center,angle,copy=True):
        if copy:
            return Circle2D(self.center.vector.Rotation(center,angle,copy=True),self.radius)
        else:
            self.center.Rotation(center,angle,copy=False) 
            
    def Translation(self,offset,copy=True):
        if copy:
            return Circle2D(self.center.vector.Translation(offset,copy=True),self.radius)
        else:
            self.center.Translation(offset,copy=False)


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

        return  (0.5*npy.abs(npy.dot(x,npy.roll(y,1))-npy.dot(y,npy.roll(x,1))))
        
    
    def PointBelongs(self,point):
        """
        Ray casting algorithm copied from internet...
        """
        n = len(self.points)
        inside = False
        x,y=point.vector
    
        p1x,p1y = self.points[0].vector
        
        for i in range(n+1):
            p2x,p2y = self.points[i % n].vector
            if y > min(p1y,p2y):
                if y <= max(p1y,p2y):
                    if x <= max(p1x,p2x):
                        if p1y != p2y:
                            xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                        if p1x == p2x or x <= xints:
                            inside = not inside
            p1x,p1y = p2x,p2y
    
        return inside

class Primitive3D:
    def __init__(self,name=''):
        self.name=name        

        
class Vector3D:
    def __init__(self,vector):
        self.vector=npy.array(vector)

    
    def __add__(self,point2d):
        return Vector3D(self.vector+point2d.vector)
    
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
        

class Point3D(Vector3D):
    def __init__(self,vector,name=''):
        Vector3D.__init__(self,vector)
        self.name=name
        
class Line3D(Primitive3D):
    def __init__(self,point1,point2,name=''):
        Primitive3D.__init__(self,name)        
        self.points=[point1,point2]
        
    def MPLPlot(self,ax):
        x=[p.vector[0] for p in self.points]
        y=[p.vector[1] for p in self.points]
        z=[p.vector[2] for p in self.points]
        ax.plot(x,y,z)
        
        
    def FreeCADExport(self,name):
        x1,y1,z1=self.points[0].vector
        x2,y2,z2=self.points[1].vector
        return '{}=Part.Line(fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,x1,y1,z1,x2,y2,z2)


class Circle3D(Primitive3D):
    def __init__(self,center,radius,normal,name):        
        Primitive2D.__init__(self,name)        
        self.center=center
        self.radius=radius
        self.normal=normal

    def FreeCADExport(self,name):
        xc,yc,zc=self.center.vector
        xn,yn,zn=self.normal.vector
        return '{}=Part.Circle(fc.Vector({},{},{}),fc.Vector({},{},{}),{})\n'.format(name,xc,yc,zc,xn,yn,zn,self.radius)
        

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
        
        
    
    def FreeCADExport(self,name):
        xs,ys,zs=self.start.vector
        xm,ym,zm=self.middle.vector
        xe,ye,ze=self.end.vector
        return '{}=Part.Arc(fc.Vector({},{},{}),fc.Vector({},{},{}),fc.Vector({},{},{}))\n'.format(name,xs,ys,zs,xm,ym,zm,xe,ye,ze)

        

class VolumeModel:
    def __init__(self,primitives):
        self.primitives=primitives
    
    def MPLPlot(self):
        """
        Matplotlib plot of model.
        To use for debug.
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
#        ax.set_aspect('equal')
        for primitive in self.primitives:
            primitive.MPLPlot(ax)
    
    def FreeCADScript(self,fcstd_filepath,path_lib_freecad='',py_filepath=''):
        """
        Generate python a FreeCAD definition of model 
        :param fcstd_filename: a filename without extension to give the name at the fcstd part written in python code
        :param py_filepath: python file to be exported. If none give, only return in string
        :type fcstd_filename:str
        :type py_filepath:str
        """
        s=''
        if path_lib_freecad!='':
            s+="import sys\nsys.path.append('"+path_lib_freecad+"')\n"

        s+="import FreeCAD as fc\nimport Part\n\ndoc=fc.newDocument('doc')\n\n"
        
        for ip,primitive in enumerate(self.primitives):
            print(ip)
            s+=(primitive.FreeCADExport(ip)+'\n')
            s+=('shapeobj = doc.addObject("Part::Feature","primitive'+str(ip)+' '+primitive.name+'")\n')
            s+=("shapeobj.Shape = primitive"+str(ip)+'\n\n')
        s+='doc.recompute()\n'
        s+="doc.saveAs('"+fcstd_filepath+".fcstd')"

        if py_filepath!='':
            with open(py_filepath,'w') as fichier:
                fichier.write(s)
        return s
            
    
    
    def FreeCADExport(self,python_path,fcstd_filepath,path_lib_freecad=''):
        """
        Export model to .fcstd FreeCAD standard
        
        :param python_path: path of python binded to freecad
             - on windows: 'something like C:\Program Files\FreeCAD X.XX\bin\python'
             - on linux: python (in general)
        :param filepath: path of fcstd file (without extension)
        :param path_lib_freecad: FreeCAD.so lib path (/usr/lib/freecad/lib in general)

        """
        s=self.FreeCADScript(fcstd_filepath,path_lib_freecad)
        print(s)
        import subprocess
        arg='-c\n'+s
        rep=subprocess.call([python_path,arg])
        print(rep)