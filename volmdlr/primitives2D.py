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
        
        
    def OffsetSingleLine(self, line_index, offset):
        """
        line_index = 0 being the 1st line 
        """
        new_linesegment2D_points = []
        dont_add_last_point = False

        for i, point in enumerate(self.points[:-1]+(self.closed)*[self.points[-1]]):

            if i == line_index:
                # Not closed RLS2D and the offset line is the last one
                if i == len(self.points)-2:
                    dir_vec_1 = volmdlr.Vector2D(point-self.points[i-1])
                    dir_vec_1.Normalize()
                    dir_vec_2 = dir_vec_1
                    dont_add_last_point = True
                # The offset line is the first one
                elif i == 0:
                    dir_vec_2 = volmdlr.Vector2D(self.points[i+1]-self.points[i+2])
                    dir_vec_2.Normalize()
                    if not self.closed:
                        dir_vec_1 = dir_vec_2
                    else:
                        dir_vec_1 = volmdlr.Vector2D(point-self.points[i-1])
                        dir_vec_1.Normalize()
                # Closed RLS2D and the offset line is the last one
                elif i == len(self.points)-1:
                    dir_vec_1 = volmdlr.Vector2D(point-self.points[i-1])
                    dir_vec_1.Normalize()
                    dir_vec_2 = volmdlr.Vector2D(self.points[0]-self.points[1])
                    dir_vec_2.Normalize()
                    dont_add_last_point = True
                else:
                    dir_vec_1 = volmdlr.Vector2D(point-self.points[i-1])
                    dir_vec_1.Normalize()
                    dir_vec_2 = volmdlr.Vector2D(self.points[i+1]-self.points[i+2])
                    dir_vec_2.Normalize()
                
                
                if self.closed and line_index == len(self.points)-1:
                    normal_vector = volmdlr.Vector2D(self.points[0]-point).NormalVector(unit=True)
                else:
                    normal_vector = volmdlr.Vector2D(self.points[i+1]-point).NormalVector(unit=True)
                    
                alpha1 = math.acos(dir_vec_1.Dot(normal_vector))
                alpha2 = math.acos(dir_vec_2.Dot(normal_vector))
                
                # If 3 segments are aligned and the middle one have to be offset
                if math.isclose(math.cos(alpha1), 0, abs_tol=1e-9) or math.isclose(math.cos(alpha2), 0, abs_tol=1e-9):
                    print('ca sort direct')
                    print('direction vector', dir_vec_1, dir_vec_2)
                    print('normal vector', normal_vector)
                    print('alpha', alpha1*180/math.pi, alpha2*180/math.pi)
                    print(point, self.points[i+1])
                    return self
#                    distance_dir1 = offset
#                    distance_dir2 = offset
                    
                distance_dir1 = offset / math.cos(alpha1)
                distance_dir2 = offset / math.cos(alpha2)

                new_point1 = point + distance_dir1 * dir_vec_1
                if self.closed and line_index == len(self.points)-1:
                    new_point2 = self.points[0] + distance_dir2 * dir_vec_2
                else:
                    new_point2 = self.points[i+1] + distance_dir2 * dir_vec_2
                
                new_linesegment2D_points.append(new_point1)
                new_linesegment2D_points.append(new_point2)
                
            elif i == line_index+1:
                pass
            
            elif line_index == len(self.points)-1 and i == 0:
                pass
            else:
                new_linesegment2D_points.append(point)
            
            
        if not dont_add_last_point and not self.closed:
            new_linesegment2D_points.append(self.points[-1])
        
        rls2D = RoundedLineSegments2D(new_linesegment2D_points, self.radius,
                                                 self.closed, adapt_radius=self.adapt_radius)
        
        return rls2D
        
    
    def OffsetLines(self, line_indexes, offset):
        """
        line_indexes is a list of consecutive line indexes
        These line should all be aligned
        line_indexes = 0 being the 1st line 
        
        if self.close last line_index can be len(self.points)-1
        if not, last line_index can be len(self.points)-2
        """  
        new_linesegment2D_points = []
        
# =============================================================================
# COMPUTES THE DIRECTIVE VECTORS BETWEEN WHICH THE OFFSET WILL BE DRAWN 
# =============================================================================
        dir_vec_1 = None
        dir_vec_2 = None
            
        if line_indexes[0] == 0 and not self.closed:
            pass
        else: 
            dir_vec_1 = volmdlr.Vector2D((self.points[line_indexes[0]] - self.points[line_indexes[0]-1]))
        
        if line_indexes[-1] == len(self.points)-(2-self.closed):
            if not self.closed:
                pass
            else:
                dir_vec_2 = volmdlr.Vector2D((self.points[0] - self.points[1]))
        elif self.closed and line_indexes[-1] == len(self.points)-2:
            dir_vec_2 = volmdlr.Vector2D((self.points[line_indexes[-1]+1] - self.points[0]))
        else:
            dir_vec_2 = volmdlr.Vector2D((self.points[line_indexes[-1]+1] - self.points[line_indexes[-1]+2]))

        if dir_vec_1 is None:
            dir_vec_1 = dir_vec_2
        if dir_vec_2 is None:
            dir_vec_2 = dir_vec_1
                
        dir_vec_1.Normalize()
        dir_vec_2.Normalize()
            
    
# =============================================================================
# COMPUTES THE ANGLE BETWEEN THE NORMAL VECTOR OF THE SURFACE TO OFFSET AND 
# THE DIRETIVE VECTOR IN ORDER TO SET THE NEW POINT AT THE RIGHT DISTANCE
# =============================================================================
        normal_vectors = []
        for index in line_indexes:
            if index == len(self.points)-1:
                normal_vectors.append(volmdlr.Vector2D(self.points[0]-self.points[index]).NormalVector(unit=True))
            else:
                normal_vectors.append(volmdlr.Vector2D(self.points[index+1]-self.points[index]).NormalVector(unit=True))

        dot1 = dir_vec_1.Dot(normal_vectors[0])
        dot2 = dir_vec_2.Dot(normal_vectors[-1])
    
        if math.isclose(dot1, 0, abs_tol=1e-9):
            # call function considering the line before, because the latter and 
            # the first offset segment are parallel
            return self.OffsetLines([line_indexes[0]-1]+line_indexes, offset)
        if math.isclose(dot2, 0, abs_tol=1e-9):
            # call function considering the line after, because the latter and 
            # the last offset segment are parallel
            return self.OffsetLines(line_indexes+[line_indexes[-1]+1], offset)
        
        distance_dir1 = offset / dot1
        distance_dir2 = offset / dot2
        
        if len(line_indexes) > 1:
            intersection = volmdlr.Point2D.LinesIntersection(volmdlr.Line2D(self.points[line_indexes[0]], self.points[line_indexes[0]]+dir_vec_1), volmdlr.Line2D(self.points[line_indexes[-1]+1], self.points[line_indexes[-1]+1]+dir_vec_2))
            vec1 = intersection.PointDistance(self.points[line_indexes[0]]) * dir_vec_1
            vec2 = intersection.PointDistance(self.points[line_indexes[-1]+1]) * dir_vec_2
        
# =============================================================================
# COMPUTES THE NEW POINTS AFTER THE OFFSET
# =============================================================================
        new_points = {}
        
        new_points[line_indexes[0]] = self.points[line_indexes[0]] + distance_dir1 * dir_vec_1
        
        for nb, index in enumerate(line_indexes[1:]):
            coeff_vec_2 = volmdlr.Point2D.PointDistance(self.points[line_indexes[0]], self.points[index]) / volmdlr.Point2D.PointDistance(self.points[line_indexes[0]], self.points[line_indexes[-1]+1])
            coeff_vec_1 = 1 - coeff_vec_2
            if dir_vec_1.Dot(normal_vectors[nb+1]) < 0:
                coeff_vec_1 = - coeff_vec_1
            if dir_vec_2.Dot(normal_vectors[nb+1]) < 0:
                coeff_vec_2 = - coeff_vec_2
            index_dir_vector = coeff_vec_1 * vec1 + coeff_vec_2 * vec2
            index_dot = index_dir_vector.Dot(normal_vectors[nb+1])
            index_distance_dir = offset / index_dot
            new_points[index] = self.points[index] + index_distance_dir * index_dir_vector
            
        if self.closed and line_indexes[-1] == len(self.points)-1:
            new_points[0] = self.points[0] + distance_dir2 * dir_vec_2
        else:
            new_points[line_indexes[-1]+1] = self.points[line_indexes[-1]+1] + distance_dir2 * dir_vec_2
        
# =============================================================================
# CREATE THE NEW POINTS' LIST 
# ============================================================================= 
        for i in range(len(self.points)):
            if i in new_points.keys():
                new_linesegment2D_points.append(new_points[i])
            else:
                new_linesegment2D_points.append(self.points[i])
                
        
        rls2D = RoundedLineSegments2D(new_linesegment2D_points, self.radius,
                                                 self.closed, adapt_radius=self.adapt_radius)
        

        return rls2D