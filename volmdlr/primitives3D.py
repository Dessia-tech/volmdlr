#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:08:23 2017

@author: steven
"""
import numpy as npy
from scipy.linalg import norm
import volmdlr.geometry as geometry
import math

import matplotlib.pyplot as plt
#import matplotlib.patches as patches
#from matplotlib.collections import PatchCollection

from volmdlr.primitives2D import Line2D,Arc2D,Point2D

        
class Vector3D:
    def __init__(self,vector):
        self.vector=npy.array(vector)

    
    def __add__(self,point2d):
        return Vector3D(self.point+point2d.point)
    
    def __radd__(self,point2d):
        return self+point2d

    def __sub__(self,point2d):
        return Vector3D(self.point-point2d.point)
    
    def __rsub__(self,point2d):
        return self-point2d    
    
    def __mul__(self,value):
        return Vector3D(self.point*value)
    
    def __rmul__(self,value):
        return self*value

    def __div__(self,value):
        return Vector3D(self.point/value)
    
    def __rdiv__(self,value):
        return self/value
        

class Point3D(Vector3D):
    def __init__(self,vector,name=''):
        Vector3D.__init__(self,vector)
        self.name=name
        

class Primitive3D:
    def __init__(self,name=''):
        self.name=name        
        
class Cylinder(Primitive3D):
    def __init__(self,position,axis,radius,width,name=''):
        Primitive3D.__init__(self,name)
        self.position=position
        self.axis=axis/norm(axis)
        self.radius=radius
        self.width=width
        
    def FreeCADExport(self,ip):
        name='primitive'+str(ip)
        e=str(self.width)
        r=str(self.radius)
        position=(self.position-self.axis*self.width/2)
        x,y,z=position
        ax,ay,az=self.axis
        x=str(x)
        y=str(y)
        z=str(z)
        ax=str(ax)
        ay=str(ay)
        az=str(az)
        return name+'=Part.makeCylinder('+r+','+e+',fc.Vector('+x+','+y+','+z+'),fc.Vector('+ax+','+ay+','+az+'),360)'
  
class HollowCylinder(Primitive3D):
    def __init__(self,position,axis,inner_radius,outer_radius,width,name=''):
        Primitive3D.__init__(self,name)
        self.position=position
        self.axis=axis/norm(axis)
        self.inner_radius=inner_radius
        self.outer_radius=outer_radius
        self.width=width

    def FreeCADExport(self,ip):
        name='primitive'+str(ip)
        re=str(self.outer_radius)
        ri=str(self.inner_radius)        
        position=self.position-self.axis*self.width/2
        x,y,z=position
        ax,ay,az=self.axis
        x=str(x)
        y=str(y)
        z=str(z)
        ax=str(ax)
        ay=str(ay)
        az=str(az)
#        return 'Part.makeCylinder('+r+','+e+',fc.Vector('+x+','+y+','+z+'),fc.Vector('+ox+','+oy+','+oz+'),360)'

        s='C2= Part.makeCircle('+re+',fc.Vector('+x+','+y+','+z+'),fc.Vector('+ax+','+ay+','+az+'))\n'
        s+='W2=Part.Wire(C2.Edges)\n'
        s+='F2=Part.Face(W2)\n'
        
        if self.inner_radius!=0.:
            s+='C1= Part.makeCircle('+ri+',fc.Vector('+x+','+y+','+z+'),fc.Vector('+ax+','+ay+','+az+'))\n'        
            s+='W1=Part.Wire(C1.Edges)\n'
            s+='F1=Part.Face(W1)\n'        
            s+='F2=F2.cut(F1)\n'        

        vx,vy,vz=self.axis*self.width
        vx=str(vx)
        vy=str(vy)
        vz=str(vz)
        
        s+=name+'=F2.extrude(fc.Vector('+vx+','+vy+','+vz+'))\n'
        return s
    
    
class ExtrudedProfile(Primitive3D):
    """
    :param points: a list of 3D numpy arrays
    :param radius: a dict containing link between index of rounded angles (keys) and radius (values)
    """
    def __init__(self,plane_origin,x1,x2,primitives,radius,extrusion_vector,name=''):
        Primitive3D.__init__(self,name)
        self.plane_origin=Point3D(plane_origin)
        self.x1=Vector3D(x1)
        self.x2=Vector3D(x2)
        self.primitives=primitives
        self.extrusion_vector=extrusion_vector
        
#    def Export2D(self,point,x1,x2):
#        """
#        export in plane defined by point,x1,x2
#        """
#        point=Point3D(point)
#        x1=Point3D(x1)
#        x2=Point3D(x2)
#
#        points2D=[geometry.PointLocalProjectionPlane(p,point,x1,x2) for p in self.points3D]
#        p1s=[points2D[-1]]+points2D[:-1]
#        pis=points2D
#        p2s=points2D[1:]+[points2D[0]]
#        points_l=points2D[:]
#        arcs=[]
#        for i in range(len(points2D)):
#            try:
#                r=self.radius[i]
#                pt1=p1s[i].vector
#                pti=pis[i].vector
#                pt2=p2s[i].vector
#
#                dist1=norm(pt1-pti)
#                dist2=norm(pt2-pti)
#                dist3=norm(pt1-pt2)
#                alpha=math.acos(-(dist3**2-dist1**2-dist2**2)/(2*dist1*dist2))/2
#                dist=r/math.tan(alpha)
#
#                vec1=(pt1-pti)/dist1
#                vec2=(pt2-pti)/dist2
#
#                indice=-vec1[1]*vec2[0]+vec1[0]*vec2[1]
#
#                indice=indice/abs(indice)
#
#                p3=pti+vec1*dist
#                p4=pti+vec2*dist
#
#                ptcx=p3[0]-indice*vec1[1]*r
#                ptcy=p3[1]+indice*vec1[0]*r
#                pc=npy.array((ptcx,ptcy))
#
#                p3c=(p3-pc)
#                p4c=(p4-pc)
#                theta1=npy.arctan2(p3c[1],p3c[0])
#                theta2=npy.arctan2(p4c[1],p4c[0])
##                theta1,theta2=sorted([theta1,theta2])
#                points_l[i]=(Point2D(p3),Point2D(p4))                           
#                arcs.append(Arc2D(Point2D(pc),r,theta1,theta2))
#            except KeyError:
#                pass
#
#        lines=[]
#        try:
#            last_point=points_l[0][1]
#        except IndexError:
#            last_point=points_l[0]
#        for p in points_l[1:]+[points_l[0]]:
#            if type(p)==tuple:
#                lines.append(Line2D(last_point,p[0]))
#                last_point=p[1]
#            else:
#                lines.append(Line2D(last_point,p))
#                last_point=p
#                
#        return lines,arcs
    
    def MPLPlot(self,point,x1,x2):
        lines,arcs=self.Export2D(point,x1,x2)
        fig, ax = plt.subplots()
        ps=[]
        for element in lines+arcs:
            ps.extend(element.MPLPlot())
#        p = PatchCollection(ps)
#        print(arcs)
        for p in ps:
            ax.add_patch(p)
        ax.set_aspect('equal')
        ax.margins(0.1)
        plt.show() 
        
    def FreeCADExport(self,ip):
        name='primitive'+str(ip)
        s='L=[]\n'
        for p1,p2 in zip(self.points3D,self.points3D[1:]+[self.points3D[0]]):
            s+='L.append(Part.makeLine({p1},{p2}))\n'.format(p1=tuple(p1.vector),p2=tuple(p2.vector))
        s+='C=Part.Wire(L)\n'
        s+='F=Part.Face(C)\n'
        e1,e2,e3=self.extrusion_vector
        s+=name+'=F.extrude(fc.Vector({},{},{}))\n'.format(e1,e2,e3)

        return s
