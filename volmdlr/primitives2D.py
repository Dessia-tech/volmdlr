#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 14:10:59 2017

@author: steven
"""
import numpy as npy
import matplotlib.pyplot as plt

from matplotlib.patches import Arc


class Point2D:
    def __init__(self,point,name=''):
        self.point=npy.array(point)
        self.name=name     
    
    def __add__(self,point2d):
        return Point2D(self.point+point2d.point)
    
    def __radd__(self,point2d):
        return self+point2d

    def __sub__(self,point2d):
        return Point2D(self.point-point2d.point)
    
    def __rsub__(self,point2d):
        return self-point2d    
    
    def __mul__(self,value):
        return Point2D(self.point*value)
    
    def __rmul__(self,value):
        return self*value

    def __div__(self,value):
        return Point2D(self.point/value)
    
    def __rdiv__(self,value):
        return self/value


class Primitive2D:
    def __init__(self,name=''):
        self.name=name

class Line2D(Primitive2D):
    def __init__(self,point1,point2,name=''):
        Primitive2D.__init__(self,name)        
        self.points=[point1,point2]

    def PltPlot(self):
        p1,p2=self.points
        x1=p1.point
        x2=p2.point
        plt.plot([x1[0],x2[0]],[x1[1],x2[1]],'black')        
        return []
    
class Arc2D(Primitive2D):
    def __init__(self,center,radius,angle1,angle2,name=''):        
        Primitive2D.__init__(self,name)        
        self.center=center
        self.radius=radius
        self.angle1=angle1
        self.angle2=angle2
        
    def PltPlot(self):
        pc=self.center.point
        print(pc,self.radius,self.angle1,self.angle2)
        return [Arc(pc,self.radius,self.radius,angle=0,theta1=self.angle1,theta2=self.angle2,color='r')]
#        plt.plot(pc[0],pc[1],'ob')
#        return[Arc((0,0),0.1,0.1,angle=0,theta1=90,theta2=235,color='red')]