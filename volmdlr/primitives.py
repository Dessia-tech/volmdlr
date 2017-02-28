#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 14:08:23 2017

@author: steven
"""
import numpy as npy
from scipy.linalg import norm

class Primitive:
    def __init__(self,name=''):
        self.name=name

        
        
class Cylinder:
    def __init__(self,position,axis,radius,width,name=''):
        Primitive.__init__(self,name)
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
  
class HollowCylinder:
    def __init__(self,position,axis,inner_radius,outer_radius,width,name=''):
        Primitive.__init__(self,name)
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
    
    
class ExtrudedProfile:
    def __init__(self,points,rounded_vertices,radius,extrusion_vector,name=''):
        Primitive.__init__(self,name)
        self.points=points
        self.rounded_vertices=rounded_vertices
        self.radius=radius
        self.extrusion_vector=extrusion_vector
        
    def FreeCADExport(self,ip):
        name='primitive'+str(ip)
        s='L=[]\n'
        for p1,p2 in zip(self.points,self.points[1:]+[self.points[0]]):
            s+='L.append(Part.makeLine({p1},{p2}))\n'.format(p1=p1,p2=p2)
        s+='C=Part.Wire(L)\n'
        s+='F=Part.Face(C)\n'
        e1,e2,e3=self.extrusion_vector
        s+=name+'=F.extrude(fc.Vector({},{},{}))\n'.format(e1,e2,e3)

        return s
