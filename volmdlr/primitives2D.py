#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 14:10:59 2017

@author: steven
"""
import numpy as npy
import math

#from scipy.linalg import norm
import volmdlr


class RoundedLineSegments2D(volmdlr.CompositePrimitive2D):
    def __init__(self,points,radius,closed=False,name=''):        
        self.points=points
        self.radius=radius
        self.closed=closed
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
        arcs=[]
        for i in sorted(self.radius.keys()):
            r=self.radius[i]
            pt1=p1s[i]
            pti=pis[i]
            pt2=p2s[i]

            dist1 = (pt1-pti).Norm()
            dist2 = (pt2-pti).Norm()
            dist3 = (pt1-pt2).Norm()
            alpha = math.acos(-(dist3**2 - dist1**2 - dist2**2) / (2*dist1*dist2)) / 2.
            dist = r / math.tan(alpha)

            u1 = (pt1-pti) / dist1
            u2 = (pt2-pti) / dist2

            p3 = pti + u1 * dist
            p4 = pti + u2 * dist

            w = (u1+u2)
            w.Normalize()

            v1 = u1.NormalVector()
            if v1.Dot(w) < 0:
                v1 = -v1

            pc = p3 + v1 * r
            pm = pc - r*w
            
            points_l[i] = (volmdlr.Point2D(p3) ,volmdlr.Point2D(p4))                           
            arcs.append(volmdlr.Arc2D(volmdlr.Point2D(p3),
                                      volmdlr.Point2D(pm),
                                      volmdlr.Point2D(p4)))
            
        lines = []
        if type(points_l[0]) == tuple:
            last_point = points_l[0][1]# Getting p4 of first point
        else:
            last_point = points_l[0]
            
        for p in points_l[1:]+[points_l[0]]:
            if type(p)==tuple:
                lines.append(volmdlr.LineSegment2D(last_point, p[0]))
                last_point=p[1]
            else:
                lines.append(volmdlr.LineSegment2D(last_point, p))
                last_point=p
        if not closed:
            del lines[-1]
            
        primitives=lines[:]
        nlines=len(lines)
        narcs=len(arcs)
        for ii,i in enumerate(sorted(radius.keys(),reverse=True)):
            if i>nlines:
                primitives.append(arcs[narcs-ii-1])
            else:
                primitives.insert(i,arcs[narcs-ii-1])
        
        volmdlr.CompositePrimitive2D.__init__(self,primitives,name)        
        
        
    def Rotation(self, center, angle, copy=True):
        if copy:
            return RoundedLineSegments2D([p.Rotation(center,angle,copy=True) for p in self.points],
                                          self.radius, self.closed, self.name)
        else:
            self.__init__([p.Rotation(center,angle,copy=True) for p in self.points],
                           self.radius, self.closed, self.name)
            
    def Translation(self, offset, copy=True):
        if copy:
            return RoundedLineSegments2D([p.Translation(offset, copy=True) for p in self.points],
                                          self.radius, self.closed, self.name)
        else:
            self.__init__([p.Translation(offset,copy=True) for p in self.points],
                           self.radius, self.closed, self.name)


