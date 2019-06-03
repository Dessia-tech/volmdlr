#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import math
import volmdlr
from volmdlr.primitives import RoundedLineSegments


class RoundedLineSegments2D(volmdlr.Wire2D, RoundedLineSegments):
    def __init__(self, points, radius, closed=False, adapt_radius=False, name=''):
        primitives = RoundedLineSegments.__init__(self, points, radius,
                                                  volmdlr.LineSegment2D, volmdlr.Arc2D,
                                                  closed, adapt_radius, name='')

        volmdlr.Wire2D.__init__(self, primitives, name)


    def ArcFeatures(self, ipoint):
        radius = self.radius[ipoint]
        if self.closed:
            if ipoint == 0:
                pt1 = self.points[-1]
            else:
                pt1 = self.points[ipoint -1]
            pti = self.points[ipoint]
            if ipoint < self.npoints-1:
                pt2 = self.points[ipoint+1]
            else:
                pt2 = self.points[0]
        else:
            pt1 = self.points[ipoint - 1]
            pti = self.points[ipoint]
            pt2 = self.points[ipoint + 1]

        # TODO: change to PointDistance
        dist1 = (pt1-pti).Norm()
        dist2 = (pt2-pti).Norm()
        dist3 = (pt1-pt2).Norm()
        alpha = math.acos(-(dist3**2 - dist1**2 - dist2**2) / (2*dist1*dist2)) / 2.
        dist = radius / math.tan(alpha)

        u1 = (pt1-pti) / dist1
        u2 = (pt2-pti) / dist2

        p3 = pti + u1 * dist
        p4 = pti + u2 * dist

        w = (u1+u2)
        if w != volmdlr.Vector2D((0, 0)):
            w.Normalize()

        v1 = u1.NormalVector()
        if v1.Dot(w) < 0:
            v1 = -v1

        pc = p3 + v1 * radius
        pm = pc - radius*w

        return p3, pm, p4, dist, alpha


    def Rotation(self, center, angle, copy=True):
        if copy:
            return RoundedLineSegments2D([p.Rotation(center, angle, copy=True)\
                                          for p in self.points],
                                          self.radius, self.closed,
                                          adapt_radius =self.adapt_radius,
                                          name = self.name)
        else:
            self.__init__([p.Rotation(center, angle, copy=True) for p in self.points],
                           self.radius, self.closed,
                           adapt_radius=self.adapt_radius, name = self.name)

    def Translation(self, offset, copy=True):
        if copy:
            return RoundedLineSegments2D([p.Translation(offset, copy=True) for p in self.points],
                                          self.radius, self.closed, adapt_radius=self.adapt_radius, name = self.name)
        else:
            self.__init__([p.Translation(offset,copy=True) for p in self.points],
                           self.radius, self.closed, adapt_radius =self.adapt_radius, name = self.name)

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

#    def Offset(self, offset):
#        offset_lines = []
#        if self.closed:
#            lines = zip(self.points, self.points[1:]+[self.points[0]])
#        else:
#            lines = zip(self.points[:-1], self.points[1:])
#
#        # Offseting lines
#        for point1, point2 in lines:
#            n = (point2 - point1).NormalVector()
#            n.Normalize()
#            point1_offset = point1 + offset*n
#            point2_offset = point2 + offset*n
#            offset_lines.append(volmdlr.Line2D(point1_offset, point2_offset))
#
#        # Computing
#        if self.closed:
#            zip_offset_lines = zip([offset_lines[-1]]+offset_lines[:-1], offset_lines)
#        else:
#            zip_offset_lines = zip(offset_lines[:-1], offset_lines[1:])
#
#        offset_points = []
#        new_radii = {}
#        for ipoint, (line1, line2) in enumerate(zip_offset_lines):
#            if volmdlr.Point2D.LinesIntersection(line1, line2) is not None:
#                offset_points.append(volmdlr.Point2D.LinesIntersection(line1, line2))
#            #######################
#            # VERIFIER CES LIGNES #
#            #######################
#            ipoint2 = ipoint+(not self.closed)
#            if ipoint2 in self.radius:
#                n1 = line1.NormalVector()
#                u2 = line2.DirectionVector()
#
#                if n1.Dot(u2)*offset < 0.:
#                    new_radius = self.radius[ipoint2] + abs(offset)
#                else:
#                    new_radius = self.radius[ipoint+(not self.closed)] - abs(offset)
#
#                if new_radius > 0:
#                    new_radii[ipoint2] = new_radius
#                else:
#                    if self.adapt_radius:
#                        new_radii[ipoint2] = 1e-6
#
#        # If not closed, end points should be offset also
#        if not self.closed:
#            vs = (self.points[1]-self.points[0]).NormalVector(unit=True)
#            ps = self.points[0] + vs*offset
#
#            ve = (self.points[-1]-self.points[-2]).NormalVector(unit=True)
#            pe = self.points[-1] + ve*offset
#            offset_points.insert(0, ps)
#            offset_points.append(pe)
#
#        return RoundedLineSegments2D(offset_points, new_radii, self.closed,
#                                     adapt_radius=self.adapt_radius)

    def Offset(self, offset):
        nb = len(self.points)
        vectors = []
        for i in range(nb-1):
            v1 = volmdlr.Vector2D((self.points[i+1] - self.points[i]))
            v2 = volmdlr.Vector2D((self.points[i] - self.points[i+1]))
            v1.Normalize()
            v2.Normalize()
            vectors.append(v1)
            vectors.append(v2)

        if self.closed:
            v1 = volmdlr.Vector2D((self.points[0] - self.points[-1]))
            v2 = volmdlr.Vector2D((self.points[-1] - self.points[0]))
            v1.Normalize()
            v2.Normalize()
            vectors.append(v1)
            vectors.append(v2)


        offset_vectors = []
        new_radii = {}
        offset_points = []

        for i in range((not self.closed),nb-(not self.closed)):

            check = False
            ni = vectors[2*i-1] + vectors[2*i]
            if ni == volmdlr.Vector2D((0,0)):
                ni = vectors[2*i]
                ni = ni.NormalVector()
                offset_vectors.append(ni)
            else:
                ni.Normalize()
                if ni.Dot(vectors[2*i-1].NormalVector()) > 0:
                    ni = - ni
                    check = True
                offset_vectors.append(ni)

            if i in self.radius:
                if (check and offset > 0) or (not check and offset < 0):
                    new_radius = self.radius[i] + abs(offset)
                else:
                    new_radius = self.radius[i] - abs(offset)
                if new_radius > 0:
                    new_radii[i] = new_radius
                else:
                    if self.adapt_radius:
                        new_radii[i] = 1e-6

            normal_vector1 = - vectors[2*i-1].NormalVector()
            normal_vector2 =   vectors[ 2*i ].NormalVector()
            normal_vector1.Normalize()
            normal_vector2.Normalize()
            alpha = math.acos(normal_vector1.Dot(normal_vector2))

            offset_point = self.points[i] + offset/math.cos(alpha/2)*offset_vectors[i-(not self.closed)]
            offset_points.append(offset_point)


        if not self.closed:
            n1 = vectors[0].NormalVector(unit=True)
            offset_vectors.insert(0, n1)
            offset_points.insert(0, self.points[0] + offset*offset_vectors[0])

            n_last = vectors[-1].NormalVector(unit=True)
            n_last = - n_last
            offset_vectors.append(n_last)
            offset_points.append(self.points[-1] + offset*offset_vectors[-1])


        return RoundedLineSegments2D(offset_points, new_radii, self.closed,
                                         adapt_radius=self.adapt_radius)




























