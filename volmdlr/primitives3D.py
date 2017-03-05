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

#from volmdlr.primitives2D import Line2D,Arc2D,Point2D
import volmdlr
        
class Cylinder(volmdlr.Primitive3D):
    def __init__(self,position,axis,radius,width,name=''):
        volmdlr.Primitive3D.__init__(self,name)
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
  
class HollowCylinder(volmdlr.Primitive3D):
    def __init__(self,position,axis,inner_radius,outer_radius,width,name=''):
        volmdlr.Primitive3D.__init__(self,name)
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
    
    
class ExtrudedProfile(volmdlr.Primitive3D):
    """
    :param points: a list of 3D numpy arrays
    :param radius: a dict containing link between index of rounded angles (keys) and radius (values)
    """
    def __init__(self,plane_origin,x,y,contours2D,extrusion_vector,name=''):
        volmdlr.Primitive3D.__init__(self,name)
        self.contours2D=contours2D
        self.extrusion_vector=extrusion_vector
        self.contours3D=[]
        for contour in contours2D:
#            print(contour)
            self.contours3D.append(contour.To3D(plane_origin,x,y))
        
    def MPLPlot(self,ax):
        for contour in self.contours3D:
            for primitive in contour:
                primitive.MPLPlot(ax)
        
    def FreeCADExport(self,ip):
        name='primitive'+str(ip)
        for ic,contour in enumerate(self.contours3D):
            s='L=[]\n'
            for ip,primitive in enumerate(contour):
                s+=primitive.FreeCADExport('L'+str(ip))
                s+='L.append(L{})\n'.format(ip)
            s+='S = Part.Shape(L)\n' 
            s+='C=Part.Wire(S.Edges)\n'
            s+='F=Part.Face(C)\n'
        e1,e2,e3=self.extrusion_vector
        s+=name+'=F.extrude(fc.Vector({},{},{}))\n'.format(e1,e2,e3)

        return s
