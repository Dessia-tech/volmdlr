#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common abstract primitives
"""

from scipy.optimize import linprog
import math
from numpy import zeros
import dessia_common as dc

class Edge(dc.DessiaObject):
    def __init__(self, start, end, name=''):
        self.start = start
        self.end = end
        dc.DessiaObject.__init__(self, name=name)

    @classmethod
    def from_step(cls, arguments, object_dict):
        if object_dict[arguments[3]].__class__ is volmdlr.primitives3D.Line3D:
            return LineSegment3D(object_dict[arguments[1]],
                                 object_dict[arguments[2]], arguments[0][1:-1])

        elif object_dict[arguments[3]].__class__ is volmdlr.primitives3D.Circle3D:
            # We supposed that STEP file is reading on trigo way
            center = object_dict[arguments[3]].center
            normal = object_dict[arguments[3]].normal
            normal.Normalize()
            radius = object_dict[arguments[3]].radius
            p1 = object_dict[arguments[1]]
            p2 = object_dict[arguments[2]]
            other_vec = object_dict[arguments[3]].other_vec
            if other_vec is None:
                other_vec = p1 - center
            other_vec.Normalize()
            frame = Frame3D(center, other_vec, normal.Cross(other_vec), normal)
            if p1 == p2:
                angle = math.pi
            else:
                # p1_new, p2_new = frame.NewCoordinates(p1), frame.NewCoordinates(p2)
                # #Angle for p1
                # u1, u2 = p1_new.vector[0]/radius, p1_new.vector[1]/radius
                # theta1 = sin_cos_angle(u1, u2)
                # #Angle for p2
                # u3, u4 = p2_new.vector[0]/radius, p2_new.vector[1]/radius
                # theta2 = sin_cos_angle(u3, u4)
                theta1, theta2 = posangle_arc(p1, p2, radius, frame)
                if theta1 > theta2:  # sens trigo
                    angle = math.pi + (theta1 + theta2) / 2
                else:
                    angle = (theta1 + theta2) / 2
            p_3 = Point3D(
                (radius * math.cos(angle), radius * math.sin(angle), 0))
            p3 = frame.OldCoordinates(p_3)
            if p1 == p3 or p2 == p3:
                p_3 = Point3D((radius * math.cos(0), radius * math.sin(0), 0))
                p3 = frame.OldCoordinates(p_3)
            arc = volmdlr.primitives3D.Arc3D(p1, p3, p2, normal, arguments[0][1:-1], other_vec)
            if math.isclose(arc.radius, 0, abs_tol=1e-9):
                if p1 == p2:
                    p_3 = Point3D(
                        (radius * math.cos(0), radius * math.sin(0), 0))
                    p3 = frame.OldCoordinates(p_3)
                    arc = volmdlr.primitives3D.Arc3D(p1, p3, p2, normal, arguments[0][1:-1],
                                other_vec)
            return arc

        elif object_dict[arguments[3]].__class__ is volmdlr.primitives3D.Ellipse3D:
            majorax = object_dict[arguments[3]].major_axis
            minorax = object_dict[arguments[3]].minor_axis
            center = object_dict[arguments[3]].center
            normal = object_dict[arguments[3]].normal
            normal.Normalize()
            majordir = object_dict[arguments[3]].major_dir
            majordir.Normalize()
            minordir = normal.Cross(majordir)
            minordir.Normalize()
            frame = Frame3D(center, majordir, minordir, normal)
            p1 = object_dict[
                arguments[1]]  # on part du principe que p1 suivant majordir
            p2 = object_dict[arguments[2]]
            if p1 == p2:
                angle = 5 * math.pi / 4
                xtra = Point3D((majorax * math.cos(math.pi / 2),
                                minorax * math.sin(math.pi / 2), 0))
                extra = frame.OldCoordinates(xtra)
            else:
                extra = None
                ## Positionnement des points dans leur frame
                p1_new, p2_new = frame.NewCoordinates(
                    p1), frame.NewCoordinates(p2)
                # Angle pour le p1
                u1, u2 = p1_new.vector[0] / majorax, p1_new.vector[1] / minorax
                theta1 = sin_cos_angle(u1, u2)
                # Angle pour le p2
                u3, u4 = p2_new.vector[0] / majorax, p2_new.vector[1] / minorax
                theta2 = sin_cos_angle(u3, u4)

                if theta1 > theta2:  # sens trigo
                    angle = math.pi + (theta1 + theta2) / 2
                else:
                    angle = (theta1 + theta2) / 2

            p_3 = Point3D(
                (majorax * math.cos(angle), minorax * math.sin(angle), 0))
            p3 = frame.OldCoordinates(p_3)

            arcellipse = ArcEllipse3D(p1, p3, p2, center, majordir, normal,
                                      arguments[0][1:-1], extra)

            return arcellipse

        elif object_dict[arguments[3]].__class__ is volmdlr.primitives3D.BSplineCurve3D:
            # print(object_dict[arguments[1]], object_dict[arguments[2]])
            # BSplineCurve3D à couper à gauche et à droite avec les points ci dessus ?
            return object_dict[arguments[3]]

        else:
            print(object_dict[arguments[3]])
            raise NotImplementedError

class Line(dc.DessiaObject):
    """
    Abstract class
    """
    def __init__(self, point1, point2, name=''):
        self.point1 = point1
        self.point2 = point2

    def unit_direction_vector(self):
        u = self.direction_vector()
        u.Normalize()
        return u

    def direction_vector(self):
        return self.point2 - self.point1

    def normal_vector(self):
        return self.unit_direction_vector().normal_vector()

    def point_projection(self, point):

        u = self.point2 - self.point1
        norm_u = u.Norm()
        t = (point - self.point1).Dot(u) / norm_u ** 2
        projection = self.point1 + t * u

        return projection, t * norm_u

    def split(self, split_point):
        return [self.__class__(self.point1, split_point),
                self.__class__(split_point, self.point2)]

class LineSegment(Edge):
    """
    Abstract class
    """

    def unit_direction_vector(self):
        u = self.direction_vector()
        u.Normalize()
        return u

    def direction_vector(self):
        return self.end - self.start

    def normal_vector(self):
        return self.unit_direction_vector().normal_vector()

    def point_projection(self, point):
        p1, p2 = self.points
        u = p2 - p1
        norm_u = u.Norm()
        t = (point - p1).Dot(u) / norm_u ** 2
        projection = p1 + t * u

        return projection, t * norm_u

    def split(self, split_point):
        return [self.__class__(self.start, split_point),
                self.__class__(split_point, self.end)]


class RoundedLineSegments:
    def __init__(self, points, radius, line_class, arc_class,
                 closed=False, adapt_radius=False, name=''):

        self.points=points
        self.radius=radius
        self.closed=closed
        self.adapt_radius = adapt_radius
        self.name = name
        self.npoints = len(points)
        primitives = self._primitives(line_class, arc_class)
        
        return primitives
    
    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            return self.__class__([p.frame_mapping(frame, side, copy=True)\
                                   for p in self.points], radius=self.radius,
                                    adapt_radius=self.adapt_radius,
                                    name=self.name)
        else:
            for p in self.points:
                p.frame_mapping(frame, side, copy=False)
                    
    def _primitives(self, line_class, arc_class):
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
                                
        
                        d1 = p1.point_distance(p2)
                        d2 = p2.point_distance(p3)

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
                if ndof > 0:
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