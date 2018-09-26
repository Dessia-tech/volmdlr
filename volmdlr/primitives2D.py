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


class RoundedLineSegments2D(volmdlr.Wire2D):
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

#    def Offset(self, offset):
#        '''
#        Offset a RoundedLineSegments2D profil of the value offset (with the sign of offset)
#        '''        
#        node_list = []
#        dict_radius = {}
#        for ind_pt, pt in enumerate(self.points[:-1]):
#            ind_ptp = ind_pt + 1
#            if ind_pt == len(self.points)-1:
#                ind_ptp = 0
#            li = volmdlr.Line2D(pt,self.points[ind_ptp])
#            vect_dir = li.points[1].vector - li.points[0].vector
#            vect_normal = npy.cross(1/npy.linalg.norm(vect_dir)*vect_dir,(0,0,1))
#            li_trans = li.Translation(offset*vect_normal[0:2],True)
#            if ind_pt > 0:
#                r1 = lm.points[1] - lm.points[0]
#                angle1 = npy.arctan2(r1.vector[1], r1.vector[0])
#                r2 = li_trans.points[1] - li_trans.points[0]
#                angle2 = npy.arctan2(r2.vector[1], r2.vector[0])
#                if angle1 != angle2:
#                    ptp = volmdlr.Point2D.LinesIntersection(lm, li_trans)
#                else:
#                    ptp = li_trans.points[0]
#                node_list.append(ptp)
#            else:
#                node_list.append(li_trans.points[0])
#            lm = li_trans
#        node_list.append(node_list[0])
#        
#        dict_radius = {}
#        for num_node, radius in self.radius.items():
#            dict_radius[num_node] = radius - offset
#        
#        return RoundedLineSegments2D(node_list, dict_radius, False)
            
    def Offset(self, offset):
        offset_lines = []
        if self.closed:
            lines = zip(self.points, self.points[1:]+[self.points[0]])
        else:
            lines = zip(self.points[:-1], self.points[1:])

        # Offseting lines
        for point1, point2 in lines:
            n = (point2 - point1).NormalVector()
            n.Normalize()
            point1_offset = point1 + offset*n
            point2_offset = point2 + offset*n
            offset_lines.append(volmdlr.Line2D(point1_offset, point2_offset))
            
        # Computing
        if self.closed:
            zip_offset_lines = zip([offset_lines[-1]]+offset_lines[:-1], offset_lines)
        else:
            zip_offset_lines = zip(offset_lines[:-1], offset_lines[1:])
            
        offset_points = []
        new_radii = {}
        for ipoint, (line1, line2) in enumerate(zip_offset_lines):
            offset_points.append(volmdlr.Point2D.LinesIntersection(line1, line2))
            if ipoint in self.radius:
                n1 = line1.NormalVector()
                u2 = line2.DirectionVector()
                if n1.Dot(u2)*offset > 0.:
                    new_radius = self.radius[ipoint] - abs(offset)
                    if new_radius > 0:
                        new_radii[ipoint] = new_radius
                else:
                    new_radii[ipoint] = self.radius[ipoint] + abs(offset)
        
        return RoundedLineSegments2D(offset_points, new_radii, self.closed)