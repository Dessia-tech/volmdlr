#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common abstract primitives.

"""

import math
from typing import Dict, List

from numpy import zeros
from scipy.optimize import linprog

# import dessia_common as dc
import volmdlr.edges


class RoundedLineSegments:
    """
    Rounded Line Segments class.

    """
    _non_serializable_attributes = ['line_class', 'arc_class', 'basis_primitives', 'primitives']

    line_class = volmdlr.edges.LineSegment
    arc_class = volmdlr.edges.ArcMixin

    def __init__(self, points: List[volmdlr.Point3D], radius: Dict[str, float],
                 closed: bool = False, adapt_radius: bool = False, name: str = ''):

        self.points = points
        self.radius = {int(k): v for k, v in radius.items()}
        self.closed = closed
        self.adapt_radius = adapt_radius
        self.name = name
        self.npoints = len(points)

    def get_points(self, point_index):
        """
        Gets points to calculate the arc features.
        """
        if self.closed:
            if point_index == 0:
                pt1 = self.points[-1]
            else:
                pt1 = self.points[point_index - 1]
            pti = self.points[point_index]
            if point_index < self.npoints - 1:
                pt2 = self.points[point_index + 1]
            else:
                pt2 = self.points[0]
        else:
            pt1 = self.points[point_index - 1]
            pti = self.points[point_index]
            pt2 = self.points[point_index + 1]
        return pt1, pti, pt2

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new RoundedLineSegments.

        side = 'old' or 'new'
        """
        return self.__class__([point.frame_mapping(frame, side)
                               for point in self.points], radius=self.radius,
                              adapt_radius=self.adapt_radius,
                              name=self.name)

    def arc_features(self, point_index: int):
        raise NotImplementedError('The method arc_features should be overloaded.')

    def _primitives(self):
        alpha = {}
        dist = {}
        lines_length = {}
        # Computing optimal radii
        rounded_points_indices = [int(i) for i in sorted(self.radius.keys())]
        groups = []
        arcs = {}

        if self.radius != {}:
            group = [rounded_points_indices[0]]
            _, _, _, dist0, alpha0 = self.arc_features(rounded_points_indices[0])
            dist[rounded_points_indices[0]] = dist0
            alpha[rounded_points_indices[0]] = alpha0

            for i in rounded_points_indices[1:]:
                # Computing the arc
                _, _, _, dist2, alpha2 = self.arc_features(i)
                dist[i] = dist2
                alpha[i] = alpha2
                if i - 1 in self.radius:
                    point1 = self.points[i - 1]
                    point2 = self.points[i]
                    length = (point2 - point1).norm()
                    lines_length[i - 1] = length
                    dist1 = dist[i - 1]

                    if dist1 + dist2 <= length:
                        groups.append(group)
                        group = [i]
                    else:
                        if not self.adapt_radius:
                            raise ValueError
                        group.append(i)
                else:
                    if group:
                        groups.append(group)
                    group = [i]
            if group:
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
                    len_group = len(group)
                    if len_group == 1:
                        # Single point, reducing its radius by simple computation if needed
                        ipoint = group[0]
                        if self.closed:
                            if ipoint == 0:
                                point1 = self.points[-1]
                                point2 = self.points[0]
                                point3 = self.points[1]
                            elif ipoint == self.npoints - 1:
                                point1 = self.points[-2]
                                point2 = self.points[-1]
                                point3 = self.points[0]
                            else:
                                point1 = self.points[ipoint - 1]
                                point2 = self.points[ipoint]
                                point3 = self.points[ipoint + 1]

                        else:
                            point1 = self.points[ipoint - 1]
                            point2 = self.points[ipoint]
                            point3 = self.points[ipoint + 1]

                        distance_1 = point1.point_distance(point2)
                        distance_2 = point2.point_distance(point3)

                        if dist[ipoint] > (min(distance_1, distance_2)):
                            self.radius[ipoint] = min(self.radius[ipoint],
                                                      min(distance_1, distance_2) * math.tan(alpha[ipoint]))

                    else:
                        # Adding to dof
                        bounds.extend([(0, self.radius[j] / math.tan(alpha[j])) for j in group])
                        dof.update({j: ndof + i for i, j in enumerate(group)})
                        ndof += len_group
                        groups2.append(group)
                        neq_ub += len_group - 1

                # Constructing simplex problem
                # C matrix:
                if ndof > 0:
                    c = zeros(ndof)
                    for j, i in dof.items():
                        c[i] = -math.tan(alpha[j])

                    a_ub = zeros((neq_ub, ndof))
                    b_ub = zeros(neq_ub)
                    ieq_ub = 0

                    for group in groups2:
                        for ip1, ip2 in zip(group[:-1], group[1:]):
                            a_ub[ieq_ub, dof[ip1]] = 1
                            a_ub[ieq_ub, dof[ip2]] = 1
                            b_ub[ieq_ub] = lines_length[ip1]
                            ieq_ub += 1
                    optimized_radius_solution = linprog(c, a_ub, b_ub, bounds=bounds)

                    for ipoint, dof_point in dof.items():
                        radius = optimized_radius_solution .x[dof_point] * math.tan(alpha[ipoint])
                        if radius > 1e-10:
                            self.radius[ipoint] = radius
                        else:
                            del self.radius[ipoint]

            # Creating geometry
            # Creating arcs
            for ipoint, radius in self.radius.items():
                p_start, p_iterior, p_end, _, _ = self.arc_features(ipoint)
                arcs[ipoint] = self.arc_class.from_3_points(p_start, p_iterior, p_end)

        return self.primitives_from_arcs(arcs)

    def primitives_from_arcs(self, arcs):
        primitives = []
        # Creating lines
        for index_line in range(self.npoints - 1):
            if index_line in self.radius:
                arc1 = arcs[index_line]
                primitives.append(arc1)
                if index_line + 1 in self.radius:
                    arc2 = arcs[index_line + 1]
                    if not arc1.end.is_close(arc2.start):
                        primitives.append(self.line_class(arc1.end, arc2.start))
                else:
                    if not arc1.end.is_close(self.points[index_line + 1]):
                        primitives.append(self.line_class(arc1.end, self.points[index_line + 1]))
            else:
                p1 = self.points[index_line]
                if index_line + 1 in self.radius:
                    arc2 = arcs[index_line + 1]
                    if not p1.is_close(arc2.start):
                        primitives.append(self.line_class(p1, arc2.start))
                else:
                    primitives.append(self.line_class(p1, self.points[index_line + 1]))

        if self.closed:
            if self.npoints - 1 in self.radius:
                arc1 = arcs[self.npoints - 1]
                primitives.append(arc1)
                if 0 in self.radius:
                    arc2 = arcs[0]
                    if not arc1.end.is_close(arc2.start):
                        primitives.append(self.line_class(arc1.end, arc2.start))
                else:
                    primitives.append(self.line_class(arc1.end, self.points[0]))
            else:
                p1 = self.points[self.npoints - 1]
                if 0 in self.radius:
                    arc2 = arcs[0]
                    if not p1.is_close(arc2.start):
                        primitives.append(self.line_class(p1, arc2.start))
                else:
                    primitives.append(self.line_class(p1, self.points[0]))

        return primitives
