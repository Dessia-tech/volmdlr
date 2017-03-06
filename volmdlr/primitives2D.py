#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 14:10:59 2017

@author: steven
"""
import numpy as npy
import math

from scipy.linalg import norm


import volmdlr



class RoundedLines2D(volmdlr.CompositePrimitive2D):
    def __init__(self,points,radius,closed=False,name=''):        
        self.points=points
        self.radius=radius
        # Construncting Arcs and lines of profile
        if closed:
            p1s=[points[-1]]+points[:-1]
            pis=points
            p2s=points[1:]+[points[0]]
        else:
            p1s=[points[-1]]+points[:-2]
            pis=points[:-1]
            p2s=points[1:]
    
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
                
                d=p3-2*pc+p4
                d=d/norm(d)
                pm=pc+r*d

#                p3c=(p3-pc)
#                p4c=(p4-pc)
#                theta1=npy.arctan2(p3c[1],p3c[0])
#                theta2=npy.arctan2(p4c[1],p4c[0])
#                theta1,theta2=sorted([theta1,theta2])
                points_l[i]=(volmdlr.Point2D(p3),volmdlr.Point2D(p4))                           
                arcs.append(volmdlr.Arc2D(volmdlr.Point2D(p3),volmdlr.Point2D(pm),volmdlr.Point2D(p4)))
            except KeyError:
                pass

        lines=[]
#        if not closed:
#            del points_l[-1]
        try:
            last_point=points_l[0][1]
#            print(points_l)
        except TypeError:
            last_point=points_l[0]
        for p in points_l[1:]+[points_l[0]]:
            if type(p)==tuple:
                lines.append(volmdlr.Line2D(last_point,p[0]))
                last_point=p[1]
            else:
                lines.append(volmdlr.Line2D(last_point,p))
                last_point=p
        if not closed:
            del lines[-1]
        primitives=lines+arcs
        
        
        volmdlr.CompositePrimitive2D.__init__(self,primitives,name)        
