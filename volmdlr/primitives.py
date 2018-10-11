#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common abstract primitives
"""

#from volmdlr.core import CompositePrimitive3D
from scipy.optimize import linprog
import math
from numpy import zeros

class RoundedLineSegments:
    def __init__(self, points, radius, line_class, arc_class,
                 closed=False, adapt_radius=False, name=''):

        self.points=points
        self.radius=radius
        self.closed=closed
        self.adapt_radius = adapt_radius
        self.npoints = len(points)
        primitives = self.Primitives(line_class, arc_class)
        
        return primitives
    
    def Primitives(self, line_class, arc_class):
        alpha = {}
        dist = {}
        lines_length = {}
        # Computing optimal radii
        rounded_points_indices = sorted(self.radius.keys())
        groups = []
        arcs = {}
        primitives = []

        
        if self.radius != {}:
            group = [rounded_points_indices[0]]
            _, _, _, dist0, alpha0 = self.ArcFeatures(rounded_points_indices[0])
            dist[rounded_points_indices[0]] = dist0
            alpha[rounded_points_indices[0]] = alpha0
            
            for i in rounded_points_indices[1:]:
                # Computing the arc
                ps2, pi2, pe2, dist2, alpha2 = self.ArcFeatures(i)
                dist[i] = dist2
                alpha[i] = alpha2
                if i-1 in self.radius:
                    p1 = self.points[i-1]
                    p2 = self.points[i]
                    l = (p2 - p1).Norm()
                    lines_length[i-1] = l
                    dist1 = dist[i-1]

                    if dist1 + dist2 <= l:
                        groups.append(group)
                        group = [i]
                    else:
                        if not self.adapt_radius:
                            raise ValueError
                        group.append(i)                
                else:
                    if group != []:
                        groups.append(group)
                    group = [i]
            if group != []:
                groups.append(group)
            if self.adapt_radius:
                if self.closed:
                    if 0 in groups[0]:
                        if self.npoints in groups[-1]:
                            new_group = groups[0] + groups[-1]
                            groups[0] = new_group
                            del groups[-1]
        
                groups2 = []
                ndof = 0
                dof = {}
                neq_ub = 0
                bounds = []
                for group in groups:
                    lg = len(group)
                    if lg == 1:
                        # Single point, reducing its radius by simple computation if needed
                        ipoint = group[0]
                        if self.closed:
                            if ipoint == 0:
                                p1 = self.points[-1]
                                p2 = self.points[0]
                                p3 = self.points[1]
                            elif ipoint == self.npoints-1:
                                p1 = self.points[-2]
                                p2 = self.points[-1]
                                p3 = self.points[0]
                            else:
                                p1 = self.points[ipoint - 1]
                                p2 = self.points[ipoint]
                                p3 = self.points[ipoint +1]
        
                        else:
                            p1 = self.points[ipoint - 1]
                            p2 = self.points[ipoint]
                            p3 = self.points[ipoint +1]
                                
        
                        d1 = p1.PointDistance(p2)
                        d2 = p2.PointDistance(p3)

                        if dist[ipoint] > (min(d1, d2)):
                            self.radius[ipoint] = min(self.radius[ipoint], min(d1, d2) * math.tan(alpha[ipoint]))
        
                    else:
                        # Adding to dof
                        bounds.extend([(0, self.radius[j]/math.tan(alpha[j])) for j in group])
                        dof.update({j: ndof+i for i, j in enumerate(group)})
                        ndof += lg
                        groups2.append(group)
                        neq_ub += lg-1
                
                # Constructing simplex problem
                # Concstructing C
                C = zeros(ndof)
                for j, i in dof.items():
                    C[i] = -math.tan(alpha[j])
                    
                A_ub = zeros((neq_ub, ndof))
                b_ub = zeros(neq_ub)
                ieq_ub = 0
        
                for group in groups2:                    
                    for ip1, ip2 in zip(group[:-1], group[1:]):
                        A_ub[ieq_ub, dof[ip1]] = 1
                        A_ub[ieq_ub, dof[ip2]] = 1
                        b_ub[ieq_ub] = lines_length[ip1]
                        ieq_ub += 1
        
                d = linprog(C, A_ub, b_ub, bounds = bounds)
        
                for ipoint, dof_point in dof.items():
                    r = d.x[dof_point]*math.tan(alpha[ipoint])
                    if r > 1e-10:
                        self.radius[ipoint] = r
                    else:
                        del self.radius[ipoint]
    
            # Creating geometry
            # Creating arcs
            for ipoint, r in self.radius.items():
                ps, pi, pe, _, _ = self.ArcFeatures(ipoint)
                arcs[ipoint] = arc_class(ps, pi, pe)
        
        # Creating lines
        for iline in range(self.npoints-1):
            if iline in self.radius:
                arc1 = arcs[iline]
                primitives.append(arc1)
                if iline+1 in self.radius:
                    arc2 = arcs[iline+1]
                    if arc1.end != arc2.start:
                        primitives.append(line_class(arc1.end, arc2.start))
                else:
                    if arc1.end != self.points[iline+1]:
                        primitives.append(line_class(arc1.end, self.points[iline+1]))
            else:
                p1 = self.points[iline]
                if iline+1 in self.radius:
                    arc2 = arcs[iline+1]
                    if p1 != arc2.start:
                        primitives.append(line_class(p1, arc2.start))
                else:
                    primitives.append(line_class(p1, self.points[iline+1]))

        if self.closed:
            if self.npoints-1 in self.radius:
                arc1 = arcs[self.npoints-1]
                primitives.append(arc1)
                if 0 in self.radius:
                    arc2 = arcs[0]
                    if arc1.end != arc2.start:
                        primitives.append(line_class(arc1.end, arc2.start))
                else:
                    primitives.append(line_class(arc1.end, self.points[iline+1]))
            else:
                p1 = self.points[self.npoints-1]
                if 0 in self.radius:
                    arc2 = arcs[0]
                    if p1 != arc2.start:
                        primitives.append(line_class(p1, arc2.start))
                else:
                    primitives.append(line_class(p1, self.points[0]))
                    
        return primitives