#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 14:10:59 2017

@author: steven
"""
import numpy as npy
import math
import matplotlib.pyplot as plt
from scipy.linalg import norm

from matplotlib.patches import Arc


class Vector2D:
    def __init__(self,vector):
        self.vector=npy.array(vector)
    
    def __add__(self,point2d):
        return Vector2D(self.point+point2d.point)
    
    def __radd__(self,point2d):
        return self+point2d

    def __sub__(self,point2d):
        return Vector2D(self.point-point2d.point)
    
    def __rsub__(self,point2d):
        return self-point2d    
    
    def __mul__(self,value):
        return Vector2D(self.point*value)
    
    def __rmul__(self,value):
        return self*value

    def __div__(self,value):
        return Vector2D(self.point/value)
    
    def __rdiv__(self,value):
        return self/value
        

class Point2D(Vector2D):
    def __init__(self,vector,name=''):
        Vector2D.__init__(self,vector)
        self.name=name


class Primitive2D:
    def __init__(self,name=''):
        self.name=name
        
class Profile2D(Primitive2D):
    def __init__(self,primitives,name=''):
        Primitive2D.__init__(self,name)        
        self.primitives=primitives

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
    
class Arc2D(Primitive2D):
    def __init__(self,center,radius,angle1,angle2,name=''):        
        Primitive2D.__init__(self,name)        
        self.center=center
        self.radius=radius
        self.angle1=angle1
        self.angle2=angle2
        
    def MPLPlot(self):
        pc=self.center.vector
        return [Arc(pc,2*self.radius,2*self.radius,angle=0,theta1=self.angle1*0.5/math.pi*360,theta2=self.angle2*0.5/math.pi*360,color='black')]


class RoundedLines2D(Profile2D):
    def __init__(self,points,radius,name=''):        
        
        # Construncting Arcs and lines of profile
        p1s=[points[-1]]+points[:-1]
        pis=points
        p2s=points[1:]+[points[0]]
        points_l=points[:]
#        primitives=[]
        arcs=[]
        for i in range(len(points)):
            try:
                r=self.radius[i]
                pt1=p1s[i].vector
                pti=pis[i].vector
                pt2=p2s[i].vector

                dist1=norm(pt1-pti)
                dist2=norm(pt2-pti)
                dist3=norm(pt1-pt2)
                alpha=math.acos(-(dist3**2-dist1**2-dist2**2)/(2*dist1*dist2))/2
                dist=r/math.tan(alpha)

                vec1=(pt1-pti)/dist1
                vec2=(pt2-pti)/dist2

                indice=-vec1[1]*vec2[0]+vec1[0]*vec2[1]

                indice=indice/abs(indice)

                p3=pti+vec1*dist
                p4=pti+vec2*dist

                ptcx=p3[0]-indice*vec1[1]*r
                ptcy=p3[1]+indice*vec1[0]*r
                pc=npy.array((ptcx,ptcy))

                p3c=(p3-pc)
                p4c=(p4-pc)
                theta1=npy.arctan2(p3c[1],p3c[0])
                theta2=npy.arctan2(p4c[1],p4c[0])
#                theta1,theta2=sorted([theta1,theta2])
                points_l[i]=(Point2D(p3),Point2D(p4))                           
                arcs.append(Arc2D(Point2D(pc),r,theta1,theta2))
            except KeyError:
                pass

        lines=[]
        try:
            last_point=points_l[0][1]
        except IndexError:
            last_point=points_l[0]
        for p in points_l[1:]+[points_l[0]]:
            if type(p)==tuple:
                lines.append(Line2D(last_point,p[0]))
                last_point=p[1]
            else:
                lines.append(Line2D(last_point,p))
                last_point=p
        primitives=lines+arcs
        
        Profile2D.__init__(self,primitives,name)        
