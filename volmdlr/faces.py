#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import triangle
from typing import List, Tuple
import math
import numpy as npy
import scipy as scp
import matplotlib.pyplot as plt
import dessia_common as dc
from geomdl import BSpline
from geomdl import utilities
import volmdlr.core
import volmdlr.core_compiled
import volmdlr.edges as vme
import volmdlr.wires
import volmdlr.display
import networkx as nx


class Surface2D(volmdlr.core.Primitive2D):
    """
    A surface bounded by an outer contour
    """

    def __init__(self, outer_contour: volmdlr.wires.Contour2D,
                 inner_contours: List[volmdlr.wires.Contour2D],
                 name: str = 'name'):
        self.outer_contour = outer_contour
        self.inner_contours = inner_contours

        volmdlr.core.Primitive2D.__init__(self, name=name)

    def area(self):
        return self.outer_contour.area() - sum(
            [c.area() for c in self.inner_contours])

    def point_belongs(self, point2d: volmdlr.Point2D):
        if not self.outer_contour.point_belongs(point2d):
            return False

        for inner_contour in self.inner_contours:
            if inner_contour.point_belongs(point2d):
                return False

        return True

    def triangulation(self, min_x_density=None, min_y_density=None):

        if self.area() == 0.:
            return volmdlr.display.DisplayMesh2D([], triangles=[])

        outer_polygon = self.outer_contour.to_polygon(angle_resolution=10)

        points = [volmdlr.display.Node2D(*p) for p in outer_polygon.points]
        vertices = [(p.x, p.y) for p in points]
        n = len(outer_polygon.points)
        segments = [(i, i + 1) for i in range(n - 1)]
        segments.append((n - 1, 0))
        point_index = {p: i for i, p in enumerate(points)}
        holes = []

        for inner_contour in self.inner_contours:
            inner_polygon = inner_contour.to_polygon(angle_resolution=10)

            for point in inner_polygon.points:
                if point not in point_index:
                    points.append(point)
                    vertices.append((point.x, point.y))
                    point_index[point] = n
                    n += 1
            for point1, point2 in zip(inner_polygon.points[:-1],
                                      inner_polygon.points[1:]):
                segments.append((point_index[point1],
                                 point_index[point2]))
            segments.append((point_index[inner_polygon.points[-1]],
                             point_index[inner_polygon.points[0]]))
            rpi = inner_contour.random_point_inside()
            holes.append((rpi.x, rpi.y))

        tri = {'vertices': npy.array(vertices).reshape((-1, 2)),
               'segments': npy.array(segments).reshape((-1, 2)),
               }
        if holes:
            tri['holes'] = npy.array(holes).reshape((-1, 2))
        t = triangle.triangulate(tri, 'p')
        triangles = t['triangles'].tolist()
        np = t['vertices'].shape[0]
        points = [volmdlr.display.Node2D(*t['vertices'][i, :]) for i in
                  range(np)]

        return volmdlr.display.DisplayMesh2D(points, triangles=triangles,
                                             edges=None)
        return volmdlr.display.DisplayMesh2D([], [])

    def split_by_lines(self, lines):
        cutted_surfaces = []
        iteration_surfaces = self.cut_by_line(lines[0])

        for line in lines[1:]:
            iteration_surfaces2 = []
            for surface in iteration_surfaces:
                line_cutted_surfaces = surface.cut_by_line(line)
                llcs = len(line_cutted_surfaces)

                if llcs == 1:
                    cutted_surfaces.append(line_cutted_surfaces[0])
                else:
                    iteration_surfaces2.extend(line_cutted_surfaces)

            iteration_surfaces = iteration_surfaces2[:]

        cutted_surfaces.extend(iteration_surfaces)
        return cutted_surfaces

    def split_regularly(self, n):
        """
        Split in n slices
        """
        xmin, xmax, ymin, ymax = self.outer_contour.bounding_rectangle()
        lines = []
        for i in range(n - 1):
            xi = xmin + (i + 1) * (xmax - xmin) / n
            lines.append(vme.Line2D(volmdlr.Point2D(xi, 0),
                                    volmdlr.Point2D(xi, 1)))
        return self.split_by_lines(lines)

    def cut_by_line(self, line: vme.Line2D):
        """
        This method makes inner contour disappear for now
        """
        # try:
        splitted_outer_contours = self.outer_contour.cut_by_line(line)
        # except IndexError:
        #     ax = self.outer_contour.plot()
        #     line.plot(ax=ax, color='r')
        return [Surface2D(oc, []) for oc in splitted_outer_contours]

    def split_at_centers(self):
        """
        Split in n slices
        """
        xmin, xmax, ymin, ymax = self.outer_contour.bounding_rectangle()

        cutted_contours = []
        iteration_contours = []
        c1 = self.inner_contours[0].center_of_mass()
        c2 = self.inner_contours[1].center_of_mass()
        cut_line = vme.Line2D(c1, c2)

        iteration_contours2 = []

        sc = self.cut_by_line2(cut_line)

        iteration_contours2.extend(sc)

        iteration_contours = iteration_contours2[:]
        cutted_contours.extend(iteration_contours)

        return cutted_contours

    def cut_by_line2(self, line):
        all_contours = []
        inner_1 = self.inner_contours[0]
        inner_2 = self.inner_contours[1]

        inner_intersections_1 = inner_1.line_intersections(line)
        inner_intersections_2 = inner_2.line_intersections(line)

        Arc1, Arc2 = inner_1.split(inner_intersections_1[1],
                                   inner_intersections_1[0])
        Arc3, Arc4 = inner_2.split(inner_intersections_2[1],
                                   inner_intersections_2[0])
        new_inner_1 = volmdlr.wires.Contour2D([Arc1, Arc2])
        new_inner_2 = volmdlr.wires.Contour2D([Arc3, Arc4])

        intersections = []
        intersections.append((inner_intersections_1[0], Arc1))
        intersections.append((inner_intersections_1[1], Arc2))
        intersections += self.outer_contour.line_intersections(line)
        intersections.append((inner_intersections_2[0], Arc3))
        intersections.append((inner_intersections_2[1], Arc4))
        intersections += self.outer_contour.line_intersections(line)

        if not intersections:
            all_contours.extend([self])
        if len(intersections) < 4:
            return [self]
        elif len(intersections) >= 4:
            if isinstance(intersections[0][0], volmdlr.Point2D) and \
                    isinstance(intersections[1][0], volmdlr.Point2D):
                ip1, ip2 = sorted(
                    [new_inner_1.primitives.index(intersections[0][1]),
                     new_inner_1.primitives.index(intersections[1][1])])
                ip5, ip6 = sorted(
                    [new_inner_2.primitives.index(intersections[4][1]),
                     new_inner_2.primitives.index(intersections[5][1])])
                ip3, ip4 = sorted(
                    [self.outer_contour.primitives.index(intersections[2][1]),
                     self.outer_contour.primitives.index(intersections[3][1])])

                sp11, sp12 = intersections[2][1].split(intersections[2][0])
                sp21, sp22 = intersections[3][1].split(intersections[3][0])
                sp33, sp34 = intersections[6][1].split(intersections[6][0])
                sp44, sp43 = intersections[7][1].split(intersections[7][0])

                primitives1 = []
                primitives1.append(
                    volmdlr.edges.LineSegment2D(intersections[6][0],
                                                intersections[1][0]))
                primitives1.append(new_inner_1.primitives[ip1])
                primitives1.append(
                    volmdlr.edges.LineSegment2D(intersections[0][0],
                                                intersections[5][0]))
                primitives1.append(new_inner_2.primitives[ip5])
                primitives1.append(
                    volmdlr.edges.LineSegment2D(intersections[4][0],
                                                intersections[7][0]))
                primitives1.append(sp44)
                primitives1.extend(self.outer_contour.primitives[ip3 + 1:ip4])
                primitives1.append(sp34)

                primitives2 = []
                primitives2.append(
                    volmdlr.edges.LineSegment2D(intersections[7][0],
                                                intersections[4][0]))
                primitives2.append(new_inner_2.primitives[ip6])
                primitives2.append(
                    volmdlr.edges.LineSegment2D(intersections[5][0],
                                                intersections[0][0]))
                primitives2.append(new_inner_1.primitives[ip2])
                primitives2.append(
                    volmdlr.edges.LineSegment2D(intersections[1][0],
                                                intersections[6][0]))
                primitives2.append(sp33)
                a = self.outer_contour.primitives[:ip3]
                a.reverse()
                primitives2.extend(a)
                primitives2.append(sp43)

                all_contours.extend([volmdlr.wires.Contour2D(primitives1),
                                     volmdlr.wires.Contour2D(primitives2)])


            else:
                raise NotImplementedError(
                    'Non convex contour not supported yet')

                raise NotImplementedError(
                    '{} intersections not supported yet'.format(
                        len(intersections)))

        return all_contours

    def cut_by_line3(self, line):
        # ax=self.outer_contour.plot()
        all_contours = []
        inner = self.inner_contours[0]
        inner_2 = self.inner_contours[1]
        inner_3 = self.inner_contours[2]

        c = inner.center_of_mass()
        c_2 = inner_2.center_of_mass()
        c_3 = inner_3.center_of_mass()
        direction_vector = line.normal_vector()
        direction_line = volmdlr.edges.Line2D(c, volmdlr.Point2D(
            (direction_vector.y * c.x - direction_vector.x * c.y) / (
                direction_vector.y), 0))
        direction_line_2 = volmdlr.edges.Line2D(c_2, volmdlr.Point2D(
            (direction_vector.y * c_2.x - direction_vector.x * c_2.y) / (
                direction_vector.y), 0))

        direction_line_3 = volmdlr.edges.Line2D(c_3, volmdlr.Point2D(
            (direction_vector.y * c_3.x - direction_vector.x * c_3.y) / (
                direction_vector.y), 0))
        inner_intersections = inner.line_intersections(direction_line)
        inner_intersections_2 = inner_2.line_intersections(direction_line_2)
        inner_intersections_3 = inner_3.line_intersections(direction_line_3)
        Arc1, Arc2 = inner.split(inner_intersections[1],
                                 inner_intersections[0])
        Arc3, Arc4 = inner_2.split(inner_intersections_2[1],
                                   inner_intersections_2[0])
        Arc5, Arc6 = inner_3.split(inner_intersections_3[1],
                                   inner_intersections_3[0])
        new_inner = volmdlr.wires.Contour2D([Arc1, Arc2])
        new_inner_2 = volmdlr.wires.Contour2D([Arc3, Arc4])
        new_inner_3 = volmdlr.wires.Contour2D([Arc5, Arc6])
        intersections = []

        intersections.append((inner_intersections[0], Arc1))
        intersections.append((inner_intersections[1], Arc2))
        if len(self.outer_contour.line_intersections(direction_line)) > 2:

            intersections.append(
                self.outer_contour.line_intersections(direction_line)[0])
            intersections.append(
                self.outer_contour.line_intersections(direction_line)[2])
        else:
            intersections.append(
                self.outer_contour.line_intersections(direction_line)[0])
            intersections.append(
                self.outer_contour.line_intersections(direction_line)[1])
        intersections.append((inner_intersections_2[0], Arc3))
        intersections.append((inner_intersections_2[1], Arc4))
        if len(self.outer_contour.line_intersections(direction_line_2)) > 2:
            intersections.append(
                self.outer_contour.line_intersections(direction_line_2)[0])
            intersections.append(
                self.outer_contour.line_intersections(direction_line_2)[2])
        else:
            intersections.append(
                self.outer_contour.line_intersections(direction_line_2)[0])
            intersections.append(
                self.outer_contour.line_intersections(direction_line_2)[1])
        intersections.append((inner_intersections_3[0], Arc5))
        intersections.append((inner_intersections_3[1], Arc6))
        if len(self.outer_contour.line_intersections(direction_line_3)) > 2:

            intersections.append(
                self.outer_contour.line_intersections(direction_line_3)[0])
            intersections.append(
                self.outer_contour.line_intersections(direction_line_3)[2])
        else:
            intersections.append(
                self.outer_contour.line_intersections(direction_line_3)[0])
            intersections.append(
                self.outer_contour.line_intersections(direction_line_3)[1])

        if isinstance(intersections[0][0], volmdlr.Point2D) and \
                isinstance(intersections[1][0], volmdlr.Point2D):
            ip1, ip2 = sorted([new_inner.primitives.index(intersections[0][1]),
                               new_inner.primitives.index(
                                   intersections[1][1])])
            ip5, ip6 = sorted(
                [new_inner_2.primitives.index(intersections[4][1]),
                 new_inner_2.primitives.index(intersections[5][1])])
            ip7, ip8 = sorted(
                [new_inner_3.primitives.index(intersections[8][1]),
                 new_inner_3.primitives.index(intersections[9][1])])
            ip3, ip4 = sorted(
                [self.outer_contour.primitives.index(intersections[2][1]),
                 self.outer_contour.primitives.index(intersections[3][1])])

            sp11, sp12 = intersections[2][1].split(intersections[2][0])
            sp21, sp22 = intersections[3][1].split(intersections[3][0])
            sp33, sp34 = intersections[6][1].split(intersections[6][0])
            sp44, sp43 = intersections[7][1].split(intersections[7][0])
            sp55, sp56 = intersections[10][1].split(intersections[10][0])
            sp66, sp65 = intersections[11][1].split(intersections[11][0])

            primitives1 = []
            primitives1.append(volmdlr.edges.LineSegment2D(intersections[7][0],
                                                           intersections[5][
                                                               0]))
            primitives1.append(new_inner_2.primitives[ip5])
            primitives1.append(volmdlr.edges.LineSegment2D(intersections[6][0],
                                                           intersections[4][
                                                               0]))
            primitives1.append(sp33)
            primitives1.append(sp43)

            primitives2 = []
            primitives2.append(volmdlr.edges.LineSegment2D(intersections[6][0],
                                                           intersections[4][
                                                               0]))
            primitives2.append(new_inner_2.primitives[ip6])
            primitives2.append(volmdlr.edges.LineSegment2D(intersections[5][0],
                                                           intersections[7][
                                                               0]))
            primitives2.append(volmdlr.edges.LineSegment2D(intersections[7][0],
                                                           intersections[11][
                                                               0]))
            primitives2.append(
                volmdlr.edges.LineSegment2D(intersections[11][0],
                                            intersections[9][0]))
            primitives2.append(new_inner_3.primitives[ip7])
            primitives2.append(volmdlr.edges.LineSegment2D(intersections[8][0],
                                                           intersections[10][
                                                               0]))
            primitives2.append(sp34)

            primitives3 = []
            primitives3.append(
                volmdlr.edges.LineSegment2D(intersections[10][0],
                                            intersections[8][0]))
            primitives3.append(new_inner_3.primitives[ip8])
            primitives3.append(volmdlr.edges.LineSegment2D(intersections[9][0],
                                                           intersections[11][
                                                               0]))
            primitives3.append(sp22)
            primitives3.append(volmdlr.edges.LineSegment2D(intersections[3][0],
                                                           intersections[1][
                                                               0]))
            primitives3.append(new_inner.primitives[ip1])
            primitives3.append(volmdlr.edges.LineSegment2D(intersections[0][0],
                                                           intersections[2][
                                                               0]))
            primitives3.append(volmdlr.edges.LineSegment2D(intersections[2][0],
                                                           intersections[10][
                                                               0]))

            primitives4 = []
            primitives4.append(volmdlr.edges.LineSegment2D(intersections[3][0],
                                                           intersections[1][
                                                               0]))
            a = volmdlr.edges.Arc2D(new_inner.primitives[ip2].end,
                                    new_inner.primitives[ip2].interior,
                                    new_inner.primitives[ip2].start)
            primitives4.append(a)
            primitives4.append(volmdlr.edges.LineSegment2D(intersections[0][0],
                                                           intersections[2][
                                                               0]))
            primitives4.append(sp12)
            primitives4.append(sp21)

            # Contour2D(primitives1),Contour2D(primitives2),
            #                      Contour2D(primitives3),
            all_contours.extend([volmdlr.wires.Contour2D(primitives4)])

        else:
            raise NotImplementedError(
                '{} intersections not supported yet'.format(
                    len(intersections)))

        return all_contours

    def bounding_rectangle(self):
        return self.outer_contour.bounding_rectangle()

    def plot(self, ax=None, color='k', alpha=1, equal_aspect=False):

        if ax is None:
            fig, ax = plt.subplots()
        self.outer_contour.plot(ax=ax, color=color, alpha=alpha,
                                equal_aspect=equal_aspect)
        for inner_contour in self.inner_contours:
            inner_contour.plot(ax=ax, color=color, alpha=alpha,
                               equal_aspect=equal_aspect)

        if equal_aspect:
            ax.set_aspect('equal')

        ax.margins(0.1)
        return ax


class Surface3D(dc.DessiaObject):
    x_periodicity = None
    y_periodicity = None
    """
    Abstract class
    """

    def face_from_contours3d(self,
                             contours3d: List[volmdlr.wires.Contour3D],
                             name: str = ''):
        """
        """

        area = -1
        inner_contours2d = []
        for contour3d in contours3d:
            contour2d = self.contour3d_to_2d(contour3d)
            inner_contours2d.append(contour2d)
            contour_area = contour2d.area()
            if contour_area > area:
                area = contour_area
                outer_contour2d = contour2d

        inner_contours2d.remove(outer_contour2d)

        if isinstance(self.face_class, str):
            class_ = globals()[self.face_class]
        else:
            class_ = self.face_class

        surface2d = Surface2D(outer_contour=outer_contour2d,
                              inner_contours=inner_contours2d)
        return class_(self,
                      surface2d=surface2d,
                      name=name)

    def contour3d_to_2d(self, contour3d):
        primitives2d = []
        last_primitive = None

        should_study_periodicity = self.x_periodicity or self.y_periodicity
        for primitive3d in contour3d.primitives:
            method_name = '{}_to_2d'.format(
                primitive3d.__class__.__name__.lower())
            if hasattr(self, method_name):
                primitives = getattr(self, method_name)(primitive3d)
                if last_primitive:
                    delta_x1 = abs(
                        primitives[0].start.x - last_primitive.end.x)
                    delta_x2 = abs(primitives[-1].end.x - last_primitive.end.x)
                    delta_y1 = abs(
                        primitives[0].start.y - last_primitive.end.y)
                    delta_y2 = abs(primitives[0].end.y - last_primitive.end.y)
                    if self.x_periodicity:
                        delta_x1 = delta_x1 % self.x_periodicity
                        delta_x2 = delta_x2 % self.x_periodicity
                        if math.isclose(delta_x1, self.x_periodicity,
                                        abs_tol=1e-4):
                            delta_x1 = 0.
                        if math.isclose(delta_x2, self.x_periodicity,
                                        abs_tol=1e-4):
                            delta_x2 = 0.

                    if self.y_periodicity:
                        delta_y1 = delta_y1 % self.y_periodicity
                        delta_y2 = delta_y2 % self.y_periodicity
                        if math.isclose(delta_y1, self.y_periodicity,
                                        abs_tol=1e-4):
                            delta_y1 = 0.
                        if math.isclose(delta_y2, self.y_periodicity,
                                        abs_tol=1e-4):
                            delta_y2 = 0.

                    end_match = False
                    if (math.isclose(delta_x1, 0., abs_tol=5e-5)
                            and math.isclose(delta_y1, 0., abs_tol=5e-5)):
                        end_match = True
                    elif (math.isclose(delta_x2, 0., abs_tol=5e-5) and
                          math.isclose(delta_y2, 0., abs_tol=5e-5)):
                        end_match = True
                        primitives = [p.reverse() for p in primitives[::-1]]
                    else:
                        ax2 = contour3d.plot()
                        primitive3d.plot(ax=ax2, color='r')
                        last_primitive3d.plot(ax=ax2, color='b')
                        ax = last_primitive.plot(color='b', plot_points=True)
                        # primitives[0].plot(ax=ax ,color='r', plot_points=True)
                        for p in primitives:
                            p.plot(ax=ax, color='r', plot_points=True)
                        raise ValueError(
                            'Primitives not following each other in contour: delta={}, {}, {}, {}'.format(
                                delta_x1, delta_x2, delta_y1, delta_y2))

                    # if not end_match and should_study_periodicity:
                    #     # Study if translating does the trick
                    #     if self.x_periodicity:
                    #         math.isclose(abs(delta_x1), self.x_periodicity,
                    #                      abs_tol=1e-4)
                    #
                    # if not :
                    #     # TODO: lower abs tol, but need to have more precise points?
                    #     if math.isclose(abs(delta_x1), self.x_periodicity, abs_tol=1e-4):
                    #         # primitives = [p.translation(-delta_x*volmdlr.X2D)\
                    #         #               for p in primitives[:]]
                    #         primitives[0].start.translation(
                    #             -delta_x * volmdlr.X2D, copy=False)
                    #     elif math.isclose(abs(delta_x2), self.x_periodicity, abs_tol=1e-4):
                    #     else:
                    #         print('sn', self.__class__.__name__)
                    #         print('lp', len(primitives))
                    #         contour3d.plot(edge_details=True)
                    #         ax = last_primitive.plot(color='b')
                    #         primitives[0].plot(ax=ax ,color='r')
                    #         for p in primitives[1:]:
                    #             print(p)
                    #             p.plot(ax=ax, color='r', ends=True)
                    #         raise ValueError('Primitives not following each other in contour: deltax={}'.format(delta_x))
                    #
                    # delta_y = primitives[0].start.y - last_primitive.end.y
                    # if not math.isclose(delta_y, 0., abs_tol=1e-4):
                    #     if abs(delta_y) == self.y_periodicity:
                    #         # primitives = [p.translation(-delta_y*volmdlr.Y2D)\
                    #         #               for p in primitives[:]]
                    #         primitives[0].start.translation(
                    #             -delta_y*volmdlr.Y2D, copy=False)
                    #     else:
                    #         contour3d.plot()
                    #         raise ValueError('Primitives not following each other in contour: deltay={}'.format(delta_y))

                if primitives:
                    last_primitive = primitives[-1]
                    last_primitive3d = primitive3d
                    primitives2d.extend(primitives)
            else:
                raise NotImplementedError(
                    'Class {} does not implement {}'.format(
                        self.__class__.__name__,
                        method_name))
        return volmdlr.wires.Contour2D(primitives2d)

    def contour2d_to_3d(self, contour2d):
        primitives3d = []
        for primitive2d in contour2d.primitives:
            method_name = '{}_to_3d'.format(
                primitive2d.__class__.__name__.lower())
            if hasattr(self, method_name):
                primitives3d.extend(getattr(self, method_name)(primitive2d))
            else:
                raise NotImplementedError(
                    'Class {} does not implement {}'.format(
                        self.__class__.__name__,
                        method_name))

        return volmdlr.wires.Contour3D(primitives3d)

    def linesegment3d_to_2d(self, linesegment3d):
        """
        a line segment on a surface will be in any case a line in 2D?
        """
        return [vme.LineSegment2D(self.point3d_to_2d(linesegment3d.start),
                                  self.point3d_to_2d(linesegment3d.end))]

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        """
        Is this right?
        """
        control_points = [self.point3d_to_2d(p) \
                          for p in bspline_curve3d.control_points]
        return [vme.BSplineCurve2D(
                    bspline_curve3d.degree,
                    control_points=control_points,
                    knot_multiplicities=bspline_curve3d.knot_multiplicities,
                    knots=bspline_curve3d.knots,
                    weights=bspline_curve3d.weights,
                    periodic=bspline_curve3d.periodic)]

class Plane3D(Surface3D):
    face_class = 'PlaneFace3D'

    def __init__(self, frame: volmdlr.Frame3D, name: str = ''):
        """
        :param frame: u and v of frame describe the plane, w is the normal
        """
        self.frame = frame
        self.name = name

    def __hash__(self):
        return hash(self.frame)

    def __eq__(self, other_plane):
        if other_plane.__class__.__name__ != self.__class__.__name__:
            return False
        return (self.frame.origin == other_plane.frame.origin and
                self.frame.w.is_colinear_to(other_plane.frame.w))

    def to_dict(self):
        # improve the object structure ?
        dict_ = dc.DessiaObject.base_dict(self)
        dict_['frame'] = self.frame.to_dict()
        dict_['name'] = self.name
        dict_['object_class'] = 'volmdlr.core.Plane3D'
        return dict_

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        frame3d.normalize()
        frame = volmdlr.Frame3D(frame3d.origin,
                                frame3d.v, frame3d.w, frame3d.u)
        return cls(frame, arguments[0][1:-1])

    def to_step(self, current_id):
        frame = volmdlr.Frame3D(self.frame.origin, self.frame.w, self.frame.u,
                                self.frame.v)
        content, frame_id = frame.to_step(current_id)
        plane_id = frame_id + 1
        content += "#{} = PLANE('{}',#{});\n".format(plane_id, self.name,
                                                     frame_id)
        return content, plane_id

    @classmethod
    def from_3_points(cls, point1, point2, point3):
        """
        Point 1 is used as origin of the plane
        """
        vector1 = point2 - point1
        vector2 = point3 - point1
        vector1.normalize()
        vector2.normalize()
        normal = vector1.cross(vector2)
        normal.normalize()
        vector = normal.cross(vector1)
        frame = volmdlr.Frame3D(point1, vector1, normal.cross(vector1), normal)
        return cls(frame)

    @classmethod
    def from_normal(cls, point, normal):
        v1 = normal.deterministic_unit_normal_vector()
        v2 = v1.cross(normal)
        return cls(volmdlr.Frame3D(point, v1, v2, normal))

    @classmethod
    def from_plane_vectors(cls, plane_origin: volmdlr.Point3D,
                           plane_x: volmdlr.Vector3D,
                           plane_y: volmdlr.Vector3D):
        normal = plane_x.cross(plane_y)
        return cls(volmdlr.Frame3D(plane_origin, plane_x, plane_y, normal))

    @classmethod
    def from_points(cls, points):
        if len(points) < 3:
            raise ValueError
        elif len(points) == 3:
            return cls.from_3_points(volmdlr.Point3D(points[0].vector),
                                     volmdlr.Vector3D(points[1].vector),
                                     volmdlr.Vector3D(points[2].vector))
        else:
            points = [p.copy() for p in points]
            indexes_to_del = []
            for i, point in enumerate(points[1:]):
                if point == points[0]:
                    indexes_to_del.append(i)
            for index in indexes_to_del[::-1]:
                del points[index + 1]

            origin = volmdlr.Point3D(points[0].vector)
            vector1 = volmdlr.Vector3D(points[1] - origin)
            vector1.normalize()
            vector2_min = volmdlr.Vector3D(points[2] - origin)
            vector2_min.normalize()
            dot_min = abs(vector1.dot(vector2_min))
            for point in points[3:]:
                vector2 = volmdlr.Vector3D(point - origin)
                vector2.normalize()
                dot = abs(vector1.dot(vector2))
                if dot < dot_min:
                    vector2_min = vector2
                    dot_min = dot
            return cls.from_3_points(origin, vector1 + origin,
                                     vector2_min + origin)

    def point_on_plane(self, point):
        if math.isclose(self.frame.w.dot(point - self.frame.origin), 0,
                        abs_tol=1e-6):
            return True
        return False

    def line_intersections(self, line):
        u = line.points[1] - line.points[0]
        w = line.points[0] - self.frame.origin
        if math.isclose(self.frame.w.dot(u), 0, abs_tol=1e-08):
            return []
        intersection_abscissea = - self.frame.w.dot(w) / self.frame.w.dot(u)
        return [line.points[0] + intersection_abscissea * u]

    def linesegment_intersections(self, linesegment: vme.LineSegment3D) \
            -> List[volmdlr.Point3D]:
        u = linesegment.end - linesegment.start
        w = linesegment.start - self.frame.origin
        normaldotu = self.frame.w.dot(u)
        if math.isclose(normaldotu, 0, abs_tol=1e-08):
            return []
        intersection_abscissea = - self.frame.w.dot(w) / normaldotu
        if intersection_abscissea < 0 or intersection_abscissea > 1:
            return []
        return [linesegment.start + intersection_abscissea * u]

    def equation_coefficients(self):
        """
        returns the a,b,c,d coefficient from equation ax+by+cz+d = 0
        """
        a, b, c = self.frame.w
        d = -self.frame.origin.dot(self.frame.w)
        return (a, b, c, d)

    def plane_intersection(self, other_plane):
        line_direction = self.frame.w.cross(other_plane.frame.w)

        if line_direction.norm() < 1e-6:
            return None

        a1, b1, c1, d1 = self.equation_coefficients()
        a2, b2, c2, d2 = other_plane.equation_coefficients()

        if a1 * b2 - a2 * b1 != 0.:
            x0 = (b1 * d2 - b2 * d1) / (a1 * b2 - a2 * b1)
            y0 = (a2 * d1 - a1 * d2) / (a1 * b2 - a2 * b1)
            point1 = volmdlr.Point3D((x0, y0, 0))
        else:
            y0 = (b2 * d2 - c2 * d1) / (b1 * c2 - c1 * b2)
            z0 = (c1 * d1 - b1 * d2) / (b1 * c2 - c1 * b2)
            point1 = volmdlr.Point3D((0, y0, z0))

        point2 = point1 + line_direction
        return volmdlr.Line3D(point1, point2)

    def rotation(self, center, axis, angle, copy=True):
        # center_frame = self.frame.origin.copy()
        # center_frame.rotation(center, axis, angle, copy=False)
        if copy:
            new_frame = self.frame.rotation(center=center, axis=axis,
                                            angle=angle, copy=True)
            # new_frame.origin = center_frame
            return Plane3D(new_frame)
        else:
            self.frame.rotation(center=center, axis=axis, angle=angle, copy=False)
            # self.frame.origin = center_frame

    def translation(self, offset, copy=True):
        if copy:
            new_frame = self.frame.translation(offset, True)
            return Plane3D(new_frame)
        else:
            self.frame.translation(offset, False)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if side == 'old':
            new_origin = frame.old_coordinates(self.frame.origin)
            new_vector1 = frame.basis().old_coordinates(self.frame.u)
            new_vector2 = frame.basis().old_coordinates(self.frame.v)
            new_vector3 = frame.basis().old_coordinates(self.frame.w)
            if copy:
                return Plane3D(
                    volmdlr.Frame3D(new_origin, new_vector1, new_vector2,
                                    new_vector3), self.name)
            else:
                # self.origin = new_origin
                # self.vectors = [new_vector1, new_vector2]
                # self.normal = frame.Basis().old_coordinates(self.normal)
                # self.normal.normalize()
                self.frame.origin = new_origin
                self.frame.u = new_vector1
                self.frame.v = new_vector2
                self.frame.w = new_vector3

        if side == 'new':
            new_origin = frame.new_coordinates(self.frame.origin)
            new_vector1 = frame.basis().new_coordinates(self.frame.u)
            new_vector2 = frame.basis().new_coordinates(self.frame.v)
            new_vector3 = frame.basis().new_coordinates(self.frame.w)
            if copy:
                return Plane3D(
                    volmdlr.Frame3D(new_origin, new_vector1, new_vector2,
                                    new_vector3), self.name)
            else:
                self.frame.origin = new_origin
                self.frame.u = new_vector1
                self.frame.v = new_vector2
                self.frame.w = new_vector3

    def copy(self):
        new_frame = self.frame.copy()
        return Plane3D(new_frame, self.name)

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        self.origin.plot(ax)
        self.vectors[0].plot(ax, starting_point=self.origin, color='r')
        self.vectors[1].plot(ax, starting_point=self.origin, color='g')
        return ax

    def babylon_script(self):
        s = 'var myPlane = BABYLON.MeshBuilder.CreatePlane("myPlane", {width: 0.5, height: 0.5, sideOrientation: BABYLON.Mesh.DOUBLESIDE}, scene);\n'
        s += 'myPlane.setPositionWithLocalVector(new BABYLON.Vector3({},{},{}));\n'.format(
            self.origin[0], self.origin[1], self.origin[2])

        s += 'var axis1 = new BABYLON.Vector3({}, {}, {});\n'.format(
            self.vectors[0][0], self.vectors[0][1], self.vectors[0][2])
        s += 'var axis2 = new BABYLON.Vector3({}, {}, {});\n'.format(
            self.vectors[1][0], self.vectors[1][1], self.vectors[1][2])
        s += 'var axis3 = new BABYLON.Vector3({}, {}, {});\n'.format(
            self.normal[0], self.normal[1], self.normal[2])
        s += 'var orientation = BABYLON.Vector3.rotationFromAxis(axis1, axis2, axis3);\n'
        s += 'myPlane.rotation = orientation;\n'

        s += 'var planemat = new BABYLON.StandardMaterial("planemat", scene);\n'
        s += 'planemat.alpha = 0.4;\n'
        s += 'myPlane.material = planemat;\n'

        return s

    def point2d_to_3d(self, point2d):
        return point2d.to_3d(self.frame.origin, self.frame.u, self.frame.v)

    def point3d_to_2d(self, point3d):
        return point3d.to_2d(self.frame.origin, self.frame.u, self.frame.v)

    def contour2d_to_3d(self, contour2d):
        return contour2d.to_3d(self.frame.origin, self.frame.u, self.frame.v)

    def contour3d_to_2d(self, contour3d):
        return contour3d.to_2d(self.frame.origin, self.frame.u, self.frame.v)

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        control_points = [self.point3d_to_2d(p) \
                          for p in bspline_curve3d.control_points]
        return [vme.BSplineCurve2D(
            bspline_curve3d.degree,
            control_points=control_points,
            knot_multiplicities=bspline_curve3d.knot_multiplicities,
            knots=bspline_curve3d.knots,
            weights=bspline_curve3d.weights,
            periodic=bspline_curve3d.periodic)]

    def bsplinecurve2d_to_3d(self, bspline_curve2d):
        control_points = [self.point2d_to_3d(p) \
                          for p in bspline_curve2d.control_points]
        return [vme.BSplineCurve3D(
            bspline_curve2d.degree,
            control_points=control_points,
            knot_multiplicities=bspline_curve2d.knot_multiplicities,
            knots=bspline_curve2d.knots,
            weights=bspline_curve2d.weights,
            periodic=bspline_curve2d.periodic)]

    def rectangular_cut(self, x1: float, x2: float,
                        y1: float, y2: float, name: str = ''):

        p1 = volmdlr.Point2D(x1, y1)
        p2 = volmdlr.Point2D(x2, y1)
        p3 = volmdlr.Point2D(x2, y2)
        p4 = volmdlr.Point2D(x1, y2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([p1, p2, p3, p4])
        surface = Surface2D(outer_contour, [])
        return PlaneFace3D(self, surface, name)


PLANE3D_OXY = Plane3D(volmdlr.OXYZ)
PLANE3D_OYZ = Plane3D(volmdlr.OYZX)
PLANE3D_OZX = Plane3D(volmdlr.OZXY)


class CylindricalSurface3D(Surface3D):
    face_class = 'CylindricalFace3D'
    x_periodicity = volmdlr.TWO_PI
    """
    The local plane is defined by (theta, z)
    :param frame: frame.w is axis, frame.u is theta=0 frame.v theta=pi/2
    :param radius: Cylinder's radius
    """

    def __init__(self, frame, radius, name=''):
        self.frame = frame
        self.radius = radius
        self.name = name

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        p = volmdlr.Point3D(self.radius * math.cos(point2d.x),
                            self.radius * math.sin(point2d.x),
                            point2d.y)
        return self.frame.old_coordinates(p)

    def point3d_to_2d(self, point3d):
        x, y, z = self.frame.new_coordinates(point3d)
        u1 = x / self.radius
        u2 = y / self.radius
        # theta = volmdlr.core.sin_cos_angle(u1, u2)
        theta = math.atan2(u2, u1)
        return volmdlr.Point2D(theta, z)

    def arc3d_to_2d(self, arc3d):
        start = self.point3d_to_2d(arc3d.start)
        end = self.point3d_to_2d(arc3d.end)
        # angle = abs(start.x-end.x)
        # if arc3d.is_trigo:
        # end = start + volmdlr.Point2D(arc3d.angle, 0)
        # else:
        #     end = start + volmdlr.Point2D(-arc3d.angle, 0)
        # interior = self.point3d_to_2d(arc3d.interior)
        # if start.x < interior.x:
        #     end = start + volmdlr.Point2D(arc3d.angle, 0)
        # else:
        #     end = start - volmdlr.Point2D(arc3d.angle, 0)
        return [vme.LineSegment2D(start, end)]

    def linesegment2d_to_3d(self, linesegment2d):
        theta1, z1 = linesegment2d.start
        theta2, z2 = linesegment2d.end
        if theta1 == theta2:
            return [vme.LineSegment3D(
                self.point2d_to_3d(linesegment2d.start),
                self.point2d_to_3d(linesegment2d.end),
            )]
        elif z1 == z2:
            if abs(theta1 - theta2) == volmdlr.TWO_PI:
                return [vme.FullArc3D(center=self.frame.origin + z1 * self.frame.w,
                                      start_end=self.point2d_to_3d(linesegment2d.start),
                                      normal=self.frame.w)]
            else:
                interior = self.point2d_to_3d(linesegment2d.point_at_abscissa(linesegment2d.length() * 0.5))
                return [vme.Arc3D(
                    self.point2d_to_3d(linesegment2d.start),
                    self.point2d_to_3d(
                        volmdlr.Point2D(0.5 * (theta1 + theta2), z1)),
                    self.point2d_to_3d(linesegment2d.end),
                )]
        else:
            raise NotImplementedError('Ellipse?')

    def fullarc3d_to_2d(self, fullarc3d):
        if self.frame.w.is_colinear_to(fullarc3d.normal):
            p1 = self.point3d_to_2d(fullarc3d.start)
            return [vme.LineSegment2D(p1, p1 + volmdlr.TWO_PI * volmdlr.X2D)]
        else:
            print(fullarc3d.normal, self.frame.w)
            raise ValueError('Impossible!')

    def circle3d_to_2d(self, circle3d):
        return []

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        # TODO: enhance this, this is a non exact method!
        l = bspline_curve3d.length()
        points = [self.point3d_to_2d(bspline_curve3d.point_at_abscissa(i / 10 * l)) \
                  for i in range(11)]
        return [vme.LineSegment2D(p1, p2) \
                for p1, p2 in zip(points[:-1], points[1:])]

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, -frame3d.u
        U.normalize()
        W.normalize()
        V = W.cross(U)
        frame_direct = volmdlr.Frame3D(frame3d.origin, U, V, W)
        radius = float(arguments[2]) / 1000
        return cls(frame_direct, radius, arguments[0][1:-1])

    def to_step(self, current_id):
        frame = volmdlr.Frame3D(self.frame.origin, self.frame.w, self.frame.u,
                                self.frame.v)
        content, frame_id = frame.to_step(current_id)
        current_id = frame_id + 1
        content += "#{} = CYLINDRICAL_SURFACE('{}',#{},{});\n" \
            .format(current_id, self.name, frame_id,
                    round(1000 * self.radius, 3))
        return content, current_id

    def frame_mapping(self, frame, side, copy=True):
        basis = frame.basis()
        if side == 'new':
            new_origin = frame.new_coordinates(self.frame.origin)
            new_u = basis.new_coordinates(self.frame.u)
            new_v = basis.new_coordinates(self.frame.v)
            new_w = basis.new_coordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return CylindricalSurface3D(new_frame, self.radius,
                                            name=self.name)
            else:
                self.frame = new_frame

        if side == 'old':
            new_origin = frame.old_coordinates(self.frame.origin)
            new_u = basis.old_coordinates(self.frame.u)
            new_v = basis.old_coordinates(self.frame.v)
            new_w = basis.old_coordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return CylindricalSurface3D(new_frame, self.radius,
                                            name=self.name)
            else:
                self.frame = new_frame

    def rectangular_cut(self, theta1: float, theta2: float,
                        z1: float, z2: float, name: str = ''):

        if theta1 == theta2:
            theta2 += volmdlr.TWO_PI

        p1 = volmdlr.Point2D(theta1, z1)
        p2 = volmdlr.Point2D(theta2, z1)
        p3 = volmdlr.Point2D(theta2, z2)
        p4 = volmdlr.Point2D(theta1, z2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([p1, p2, p3, p4])
        surface2d = Surface2D(outer_contour, [])
        return volmdlr.faces.CylindricalFace3D(self, surface2d, name)

    def translation(self, offset: volmdlr.Vector3D, copy=True):
        if copy:
            return self.__class__(self.frame.translation(offset, copy=True),
                                  self.radius)
        else:
            self.frame.translation(offset, copy=False)

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            new_frame = self.frame.rotation(center=center, axis=axis,
                                            angle=angle, copy=True)
            return self.__class__(new_frame, self.radius)
        else:
            self.frame.rotation(center, axis, angle, copy=False)

class ToroidalSurface3D(Surface3D):
    face_class = 'ToroidalFace3D'
    x_periodicity = volmdlr.TWO_PI
    y_periodicity = volmdlr.TWO_PI
    """
    The local plane is defined by (theta, phi)
    theta is the angle around the big (R) circle and phi around the small(r)

    :param frame: Tore's frame: origin is the center, u is pointing at
                    theta=0
    :param R: Tore's radius
    :param r: Circle to revolute radius
    Definitions of R and r according to https://en.wikipedia.org/wiki/Torus
    """

    def __init__(self, frame: volmdlr.Frame3D,
                 R: float, r: float, name: str = ''):
        self.frame = frame
        self.R = R
        self.r = r
        self.name = name

    def _bounding_box(self):
        d = self.R + self.r
        p1 = self.frame.origin + self.frame.u * d + self.frame.v * d + self.frame.w * self.r
        p2 = self.frame.origin + self.frame.u * d + self.frame.v * d - self.frame.w * self.r
        p3 = self.frame.origin + self.frame.u * d - self.frame.v * d + self.frame.w * self.r
        p4 = self.frame.origin + self.frame.u * d - self.frame.v * d - self.frame.w * self.r
        p5 = self.frame.origin - self.frame.u * d + self.frame.v * d + self.frame.w * self.r
        p6 = self.frame.origin - self.frame.u * d + self.frame.v * d - self.frame.w * self.r
        p7 = self.frame.origin - self.frame.u * d - self.frame.v * d + self.frame.w * self.r
        p8 = self.frame.origin - self.frame.u * d - self.frame.v * d - self.frame.w * self.r

        return volmdlr.core.BoundingBox.from_points(
            [p1, p2, p3, p4, p5, p6, p7, p8])

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        theta, phi = point2d
        x = (self.R + self.r * math.cos(phi)) * math.cos(theta)
        y = (self.R + self.r * math.cos(phi)) * math.sin(theta)
        z = self.r * math.sin(phi)
        return self.frame.old_coordinates(volmdlr.Point3D(x, y, z))

    def point3d_to_2d(self, point3d):
        # points_2D = []
        x, y, z = self.frame.new_coordinates(point3d)
        if z < -self.r:
            z = -self.r
        elif z > self.r:
            z = self.r

        zr = z / self.r
        phi = math.asin(zr)

        u = self.R + math.sqrt((self.r ** 2) - (z ** 2))
        u1, u2 = round(x / u, 5), round(y / u, 5)
        theta = volmdlr.core.sin_cos_angle(u1, u2)

        return volmdlr.Point2D(theta, phi)

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, -frame3d.u
        U.normalize()
        W.normalize()
        V = W.cross(U)
        frame_direct = volmdlr.Frame3D(frame3d.origin, U, V, W)
        rcenter = float(arguments[2]) / 1000
        rcircle = float(arguments[3]) / 1000
        return cls(frame_direct, rcenter, rcircle, arguments[0][1:-1])

    def to_step(self, current_id):
        frame = volmdlr.Frame3D(self.frame.origin, self.frame.w, self.frame.u,
                                self.frame.v)
        content, frame_id = frame.to_step(current_id)
        current_id = frame_id + 1
        content += "#{} = TOROIDAL_SURFACE('{}',#{},{},{});\n" \
            .format(current_id, self.name, frame_id,
                    round(1000 * self.R, 3),
                    round(1000 * self.r, 3))
        return content, current_id

    def frame_mapping(self, frame, side, copy=True):
        basis = frame.Basis()
        if side == 'new':
            new_origin = frame.new_coordinates(self.frame.origin)
            new_u = basis.new_coordinates(self.frame.u)
            new_v = basis.new_coordinates(self.frame.v)
            new_w = basis.new_coordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ToroidalSurface3D(new_frame,
                                         self.R, self.r,
                                         name=self.name)
            else:
                self.frame = new_frame

        if side == 'old':
            new_origin = frame.old_coordinates(self.frame.origin)
            new_u = basis.old_coordinates(self.frame.u)
            new_v = basis.old_coordinates(self.frame.v)
            new_w = basis.old_coordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ToroidalSurface3D(new_frame,
                                         self.R, self.r,
                                         name=self.name)
            else:
                self.frame = new_frame

    def rectangular_cut(self, theta1, theta2, phi1, phi2, name=''):
        if phi1 == phi2:
            phi2 += volmdlr.TWO_PI
        elif phi2 < phi1:
            phi2 += volmdlr.TWO_PI
        if theta1 == theta2:
            theta2 += volmdlr.TWO_PI
        elif theta2 < theta1:
            theta2 += volmdlr.TWO_PI

        p1 = volmdlr.Point2D(theta1, phi1)
        p2 = volmdlr.Point2D(theta1, phi2)
        p3 = volmdlr.Point2D(theta2, phi2)
        p4 = volmdlr.Point2D(theta2, phi1)
        outer_contour = volmdlr.wires.ClosedPolygon2D([p1, p2, p3, p4])
        return ToroidalFace3D(self,
                              Surface2D(outer_contour, []),
                              name)

    def linesegment2d_to_3d(self, linesegment2d):
        theta1, phi1 = linesegment2d.start
        theta2, phi2 = linesegment2d.end
        if theta1 == theta2:
            if math.isclose(phi1 - phi2, volmdlr.TWO_PI, abs_tol=1e-9):
                u = self.frame.u.rotation(self.frame.origin, self.frame.w,
                                          angle=theta1)
                v = self.frame.u.rotation(self.frame.origin, self.frame.w,
                                          angle=theta1)
                center = self.frame.origin+self.R*u
                return [vme.FullArc3D(center=center,
                                      start_end=center + self.r * u,
                                      normal=v)]
            else:
                return [vme.Arc3D(
                    self.point2d_to_3d(linesegment2d.start),
                    self.point2d_to_3d(volmdlr.Point2D(theta1, 0.5 * (phi1 + phi2))),
                    self.point2d_to_3d(linesegment2d.end),
                )]
        elif math.isclose(phi1, phi2, abs_tol=1e-9):
            if abs(theta1 - theta2) == volmdlr.TWO_PI:
                center = self.frame.origin + self.r * math.sin(phi1) * self.frame.w
                start_end = center + self.frame.u * (self.r + self.R)
                return [vme.FullArc3D(center=center,
                                      start_end=start_end,
                                      normal=self.frame.w)]
            else:
                return [vme.Arc3D(
                    self.point2d_to_3d(linesegment2d.start),
                    self.point2d_to_3d(volmdlr.Point2D(0.5 * (theta1 + theta2), phi1)),
                    self.point2d_to_3d(linesegment2d.end),
                )]
        else:
            raise NotImplementedError('Ellipse?')

    def fullarc3d_to_2d(self, fullarc3d):
        if self.frame.w.is_colinear_to(fullarc3d.normal):
            p1 = self.point3d_to_2d(fullarc3d.start)
            return [vme.LineSegment2D(p1, p1 + volmdlr.TWO_PI * volmdlr.X2D)]
        elif fullarc3d.normal.dot(self.frame.w):
            p1 = self.point3d_to_2d(fullarc3d.start)
            return [vme.LineSegment2D(p1, p1 + volmdlr.TWO_PI * volmdlr.Y2D)]
        else:
            raise ValueError('Impossible!')

    def circle3d_to_2d(self, circle3d):
        return []

    def triangulation(self):
        face = self.rectangular_cut(0, volmdlr.TWO_PI, 0, volmdlr.TWO_PI)
        return face.triangulation()

    def translation(self, offset: volmdlr.Vector3D, copy=True):
        if copy:
            return self.__class__(self.frame.translation(offset, copy=True),
                                  self.R,
                                  self.r)
        else:
            self.frame.translation(offset, copy=False)

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            new_frame = self.frame.rotation(center=center, axis=axis,
                                            angle=angle, copy=True)
            return self.__class__(new_frame, self.R, self.r)
        else:
            self.frame.rotation(center, axis, angle, copy=False)

class ConicalSurface3D(Surface3D):
    face_class = 'ConicalFace3D'
    x_periodicity = volmdlr.TWO_PI
    """
    The local plane is defined by (theta, z)
    :param frame: Cone's frame to position it: frame.w is axis of cone
                    frame.origin is at the angle of the cone
    :param semi_angle: Cone's semi-angle
    """

    def __init__(self, frame: volmdlr.Frame3D, semi_angle: float,
                 name: str = ''):
        self.frame = frame
        self.semi_angle = semi_angle
        self.name = name

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, frame3d.u
        U.normalize()
        W.normalize()
        V = W.cross(U)
        radius = float(arguments[2]) / 1000
        semi_angle = float(arguments[3])
        origin = frame3d.origin - radius / math.tan(semi_angle) * W
        frame_direct = volmdlr.Frame3D(origin, U, V, W)
        return cls(frame_direct, semi_angle, arguments[0][1:-1])

    def to_step(self, current_id):
        frame = volmdlr.Frame3D(self.frame.origin, self.frame.w, self.frame.u,
                                self.frame.v)
        content, frame_id = frame.to_step(current_id)
        current_id = frame_id + 1
        content += "#{} = CONICAL_SURFACE('{}',#{},{},{});\n" \
            .format(current_id, self.name, frame_id,
                    0.,
                    round(self.semi_angle, 3))
        return content, current_id

    def frame_mapping(self, frame, side, copy=True):
        basis = frame.Basis()
        if side == 'new':
            new_origin = frame.new_coordinates(self.frame.origin)
            new_u = basis.new_coordinates(self.frame.u)
            new_v = basis.new_coordinates(self.frame.v)
            new_w = basis.new_coordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ConicalSurface3D(new_frame, self.radius, name=self.name)
            else:
                self.frame = new_frame

        if side == 'old':
            new_origin = frame.old_coordinates(self.frame.origin)
            new_u = basis.old_coordinates(self.frame.u)
            new_v = basis.old_coordinates(self.frame.v)
            new_w = basis.old_coordinates(self.frame.w)
            new_frame = volmdlr.Frame3D(new_origin, new_u, new_v, new_w)
            if copy:
                return ConicalSurface3D(new_frame, self.radius, name=self.name)
            else:
                self.frame = new_frame

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        theta, z = point2d
        r = math.tan(self.semi_angle) * z
        new_point = volmdlr.Point3D(r * math.cos(theta),
                                    r * math.sin(theta),
                                    z)
        return self.frame.old_coordinates(new_point)

    # def point3d_to_2d(self, point3d: volmdlr.Point3D):
    #     z = self.frame.w.dot(point3d)
    #     x, y = point3d.plane_projection2d(self.frame.origin, self.frame.u,
    #                                       self.frame.v)
    #     theta = math.atan2(y, x)
    #     return volmdlr.Point2D(theta, z+0.003)

    def point3d_to_2d(self, point3d: volmdlr.Point3D):
        x, y, z = self.frame.new_coordinates(point3d)
        # x, y = point3d.plane_projection2d(self.frame.origin, self.frame.u,
        #                                   self.frame.v)
        theta = math.atan2(y, x)
        return volmdlr.Point2D(theta, z)

    def rectangular_cut(self, theta1: float, theta2: float,
                        z1: float, z2: float, name: str = ''):
        # theta1 = angle_principal_measure(theta1)
        # theta2 = angle_principal_measure(theta2)
        if theta1 == theta2:
            theta2 += volmdlr.TWO_PI

        p1 = volmdlr.Point2D(theta1, z1)
        p2 = volmdlr.Point2D(theta2, z1)
        p3 = volmdlr.Point2D(theta2, z2)
        p4 = volmdlr.Point2D(theta1, z2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([p1, p2, p3, p4])
        return ConicalFace3D(self, Surface2D(outer_contour, []), name)

    def fullarc3d_to_2d(self, fullarc3d):
        if self.frame.w.is_colinear_to(fullarc3d.normal):
            p1 = self.point3d_to_2d(fullarc3d.start)
            return [vme.LineSegment2D(p1, p1 + volmdlr.TWO_PI * volmdlr.X2D)]
        else:
            raise ValueError('Impossible!')

    def circle3d_to_2d(self, circle3d):
        return []

    def linesegment2d_to_3d(self, linesegment2d):
        theta1, z1 = linesegment2d.start
        theta2, z2 = linesegment2d.end
        if math.isclose(z1, z2, abs_tol=1e-9) and math.isclose(z1, 0.,
                                                               abs_tol=1e-9):
            return []
        elif math.isclose(abs(theta1 - theta2) % volmdlr.TWO_PI, 0., abs_tol=1e-9):
            return [vme.LineSegment3D(
                self.point2d_to_3d(linesegment2d.start),
                self.point2d_to_3d(linesegment2d.end),
            )]
        elif math.isclose(z1, z2, abs_tol=1e-9):

            if abs(theta1 - theta2) % volmdlr.TWO_PI == 0.:
                return [vme.FullArc3D(center=self.frame.origin + z1 * self.frame.w,
                                      start_end=self.point2d_to_3d(linesegment2d.start),
                                      normal=self.frame.w)]
            else:
                return [vme.Arc3D(
                    self.point2d_to_3d(linesegment2d.start),
                    self.point2d_to_3d(
                        volmdlr.Point2D(0.5 * (theta1 + theta2), z1)),
                    self.point2d_to_3d(linesegment2d.end))
                ]
        else:
            raise NotImplementedError('Ellipse?')

    def translation(self, offset: volmdlr.Vector3D, copy=True):
        if copy:
            return self.__class__(self.frame.translation(offset, copy=True),
                                  self.semi_angle)
        else:
            self.frame.translation(offset, copy=False)

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            new_frame = self.frame.rotation(center=center, axis=axis, angle=angle, copy=True)
            return self.__class__(new_frame, self.semi_angle)
        else:
            self.frame.rotation(center, axis, angle, copy=False)

class SphericalSurface3D(Surface3D):
    face_class = 'SphericalFace3D'
    """
    :param frame: Sphere's frame to position it
    :type frame: volmdlr.Frame3D
    :param radius: Sphere's radius
    :type radius: float
    """

    def __init__(self, frame, radius, name=''):
        self.frame = frame
        self.radius = radius
        self.name = name
        # V = frame.v
        # V.normalize()
        # W = frame.w
        # W.normalize()
        # self.plane = Plane3D(frame.origin, V, W)

    def _bounding_box(self):
        points = [self.frame.origin + volmdlr.Point3D(-self.radius,
                                                      -self.radius,
                                                      -self.radius),
                  self.frame.origin + volmdlr.Point3D(self.radius,
                                                      self.radius,
                                                      self.radius),

                  ]
        return volmdlr.core.BoundingBox.from_points(points)

    @classmethod
    def from_step(cls, arguments, object_dict):
        frame3d = object_dict[arguments[1]]
        U, W = frame3d.v, frame3d.u
        U.normalize()
        W.normalize()
        V = W.cross(U)
        frame_direct = volmdlr.Frame3D(frame3d.origin, U, V, W)
        radius = float(arguments[2]) / 1000
        return cls(frame_direct, radius, arguments[0][1:-1])

    def point2d_to_3d(self, point2d):
        # source mathcurve.com/surfaces/sphere
        # -pi<theta<pi, -pi/2<phi<pi/2
        theta, phi = point2d
        x = self.radius * math.cos(phi) * math.cos(theta)
        y = self.radius * math.cos(phi) * math.sin(theta)
        z = self.radius * math.sin(phi)
        return self.frame.old_coordinates(volmdlr.Point3D(x, y, z))

    def point3d_to_2d(self, point3d):
        x, y, z = point3d
        if z < -self.radius:
            z = -self.radius
        elif z > self.radius:
            z = self.radius

        zr = z / self.radius
        phi = math.asin(zr)

        u = math.sqrt((self.radius ** 2) - (z ** 2))
        if u == 0:
            u1, u2 = x, y
        else:
            u1, u2 = round(x / u, 5), round(y / u, 5)
        theta = volmdlr.sin_cos_angle(u1, u2)
        return volmdlr.Point2D(theta, phi)

    def linesegment2d_to_3d(self, linesegment2d):
        start = self.point2d_to_3d(linesegment2d.start)
        interior = self.point2d_to_3d(0.5 * (linesegment2d.start + linesegment2d.end))
        end = self.point2d_to_3d(linesegment2d.end)
        if start == end:
            u = start - self.frame.origin
            u.normalize()
            v = interior - self.frame.origin
            v.normalize()
            normal = u.cross(v)
            return [vme.FullArc3D(self.frame.origin, start, normal)]
        return [vme.Arc3D(start, interior, end)]

    def plot(self, ax=None, color='grey', alpha=0.5):
        points = []
        for i in range(20):
            theta = i / 20. * volmdlr.TWO_PI
            t_points = []
            for j in range(20):
                phi = j / 20. * volmdlr.TWO_PI
                t_points.append(self.point2d_to_3d(volmdlr.Point2D(theta, phi)))
            ax = volmdlr.wires.ClosedPolygon3D(t_points).plot(ax=ax, color=color, alpha=alpha)

        return ax

    def rectangular_cut(self, theta1, theta2, phi1, phi2, name=''):
        if phi1 == phi2:
            phi2 += volmdlr.TWO_PI
        elif phi2 < phi1:
            phi2 += volmdlr.TWO_PI
        if theta1 == theta2:
            theta2 += volmdlr.TWO_PI
        elif theta2 < theta1:
            theta2 += volmdlr.TWO_PI

        p1 = volmdlr.Point2D(theta1, phi1)
        p2 = volmdlr.Point2D(theta1, phi2)
        p3 = volmdlr.Point2D(theta2, phi2)
        p4 = volmdlr.Point2D(theta2, phi1)
        outer_contour = volmdlr.wires.ClosedPolygon2D([p1, p2, p3, p4])
        return SphericalFace3D(self,
                               Surface2D(outer_contour, []),
                               name=name)


class RuledSurface3D(Surface3D):
    face_class = 'RuledFace3D'
    """
    :param frame: frame.w is axis, frame.u is theta=0 frame.v theta=pi/2
    :type frame: volmdlr.Frame3D
    :param radius: Cylinder's radius
    :type radius: float
    """

    def __init__(self,
                 wire1: volmdlr.wires.Wire3D,
                 wire2: volmdlr.wires.Wire3D,
                 name: str = ''):
        self.wire1 = wire1
        self.wire2 = wire2
        self.length1 = wire1.length()
        self.length2 = wire2.length()
        self.name = name

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        x, y = point2d
        point1 = self.wire1.point_at_abscissa(x * self.length1)
        point2 = self.wire2.point_at_abscissa(x * self.length2)
        joining_line = vme.LineSegment3D(point1, point2)
        point = joining_line.point_at_abscissa(y * joining_line.length())
        return point

    def point3d_to_2d(self, point3d):
        raise NotImplementedError

    def rectangular_cut(self, x1: float, x2: float,
                        y1: float, y2: float, name: str = ''):
        p1 = volmdlr.Point2D(x1, y1)
        p2 = volmdlr.Point2D(x2, y1)
        p3 = volmdlr.Point2D(x2, y2)
        p4 = volmdlr.Point2D(x1, y2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([p1, p2, p3, p4])
        surface2d = Surface2D(outer_contour, [])
        return volmdlr.faces.RuledFace3D(self, surface2d, name)


class BSplineSurface3D(Surface3D):
    face_class = 'BSplineFace3D'

    def __init__(self, degree_u, degree_v, control_points, nb_u, nb_v,
                 u_multiplicities, v_multiplicities, u_knots, v_knots,
                 weights=None, name=''):
        self.control_points = control_points
        self.degree_u = degree_u
        self.degree_v = degree_v
        self.nb_u = nb_u
        self.nb_v = nb_v

        u_knots = vme.standardize_knot_vector(u_knots)
        v_knots = vme.standardize_knot_vector(v_knots)
        self.u_knots = u_knots
        self.v_knots = v_knots
        self.u_multiplicities = u_multiplicities
        self.v_multiplicities = v_multiplicities
        self.weights = weights

        self.control_points_table = []
        points_row = []
        i = 1
        for pt in control_points:
            points_row.append(pt)
            if i == nb_v:
                self.control_points_table.append(points_row)
                points_row = []
                i = 1
            else:
                i += 1
        surface = BSpline.Surface()
        surface.degree_u = degree_u
        surface.degree_v = degree_v
        if weights is None:
            P = [(control_points[i][0], control_points[i][1],
                  control_points[i][2]) for i in range(len(control_points))]
            surface.set_ctrlpts(P, nb_u, nb_v)
        else:
            Pw = [(control_points[i][0] * weights[i],
                   control_points[i][1] * weights[i],
                   control_points[i][2] * weights[i],
                   weights[i]) for i in range(len(control_points))]
            surface.set_ctrlpts(Pw, nb_u, nb_v)
        knot_vector_u = []
        for i, u_knot in enumerate(u_knots):
            knot_vector_u.extend([u_knot] * u_multiplicities[i])
        knot_vector_v = []
        for i, v_knot in enumerate(v_knots):
            knot_vector_v.extend([v_knot] * v_multiplicities[i])
        surface.knotvector_u = knot_vector_u
        surface.knotvector_v = knot_vector_v
        surface.delta = 0.05
        surface_points = surface.evalpts

        self.surface = surface
        # self.points = [volmdlr.Point3D(*p) for p in surface_points]
        volmdlr.core.Primitive3D.__init__(self, name=name)

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        x, y = point2d
        return volmdlr.Point3D(*self.surface.evaluate_single((x, y)))

    def point3d_to_2d(self, point3d: volmdlr.Point3D):
        # x, y, z = point3d
        def f(x):
            return (point3d - self.point2d_to_3d(
                volmdlr.Point2D(x[0], x[1]))).norm()

        for x0 in [(0, 0), (0, 1), (1, 0), (1, 1), (0.5, 0.5)]:
            sol = scp.optimize.minimize(f, x0=x0,
                                        bounds=[(0, 1), (0, 1)],
                                        options={'eps': 1e-12})
            if sol.fun < 1e-3:
                return volmdlr.Point2D(*sol.x)

        raise RuntimeError(
            'No convergence in point3d to 2d of bspline surface')

    def linesegment2d_to_3d(self, linesegment2d):
        # TODO: this is a non exact method!
        l = linesegment2d.length()
        points = [self.point2d_to_3d(linesegment2d.point_at_abscissa(i / l / 10.)) for i in range(11)]

        return [vme.LineSegment3D(p1, p2) \
                for p1, p2 in zip(points[:-1], points[1:])]

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        # TODO: enhance this, it is a non exact  method!
        l = bspline_curve3d.length()
        points = [self.point3d_to_2d(bspline_curve3d.point_at_abscissa(i / 10 * l)) \
                  for i in range(11)]
        return [vme.LineSegment2D(p1, p2) \
                for p1, p2 in zip(points[:-1], points[1:])]

    def arc3d_to_2d(self, arc3d):
        number_points = math.ceil(arc3d.angle * 7) + 1  # 7 points per radian
        l = arc3d.length()
        points = [self.point3d_to_2d(arc3d.point_at_abscissa(i * l / (number_points - 1))) \
                  for i in range(number_points)]
        return [vme.LineSegment2D(p1, p2) \
                for p1, p2 in zip(points[:-1], points[1:])]

    def _bounding_box(self):
        return volmdlr.core.BoundingBox.from_points(self.control_points)

    def rectangular_cut(self, u1: float, u2: float,
                        v1: float, v2: float, name: str = ''):
        p1 = volmdlr.Point2D(u1, v1)
        p2 = volmdlr.Point2D(u2, v1)
        p3 = volmdlr.Point2D(u2, v2)
        p4 = volmdlr.Point2D(u1, v2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([p1, p2, p3, p4])
        surface = Surface2D(outer_contour, [])
        return PlaneFace3D(self, surface, name)

    def FreeCADExport(self, ip, ndigits=3):
        name = 'primitive{}'.format(ip)
        script = ""
        points = '['
        for i, pts_row in enumerate(self.control_points_table):
            pts = '['
            for j, pt in enumerate(pts_row):
                point = 'fc.Vector({},{},{}),'.format(pt[0], pt[1], pt[2])
                pts += point
            pts = pts[:-1] + '],'
            points += pts
        points = points[:-1] + ']'

        script += '{} = Part.BSplineSurface()\n'.format(name)
        if self.weights is None:
            script += '{}.buildFromPolesMultsKnots({},{},{},udegree={},vdegree={},uknots={},vknots={})\n'.format(
                name, points, self.u_multiplicities, self.v_multiplicities,
                self.degree_u, self.degree_v, self.u_knots, self.v_knots)
        else:
            script += '{}.buildFromPolesMultsKnots({},{},{},udegree={},vdegree={},uknots={},vknots={},weights={})\n'.format(
                name, points, self.u_multiplicities, self.v_multiplicities,
                self.degree_u, self.degree_v, self.u_knots, self.v_knots,
                self.weights)

        return script

    def rotation(self, center, axis, angle, copy=True):
        new_control_points = [p.rotation(center, axis, angle, copy=True) for p in
                              self.control_points]
        new_bsplinesurface3d = BSplineSurface3D(self.degree_u, self.degree_v,
                                                new_control_points, self.nb_u,
                                                self.nb_v,
                                                self.u_multiplicities,
                                                self.v_multiplicities,
                                                self.u_knots, self.v_knots,
                                                self.weights, self.name)
        if copy:
            return new_bsplinesurface3d
        else:
            self.control_points = new_control_points
            self.surface = new_bsplinesurface3d.surface
            # self.points = new_BSplineSurface3D.points

    def translation(self, offset, copy=True):
        new_control_points = [p.translation(offset, True) for p in
                              self.control_points]
        new_bsplinesurface3d = BSplineSurface3D(self.degree_u, self.degree_v,
                                                new_control_points, self.nb_u,
                                                self.nb_v,
                                                self.u_multiplicities,
                                                self.v_multiplicities,
                                                self.u_knots, self.v_knots,
                                                self.weights, self.name)
        if copy:
            return new_bsplinesurface3d
        else:
            self.control_points = new_control_points
            self.surface = new_bsplinesurface3d.surface
            # self.points = new_BSplineSurface3D.points

    def frame_mapping(self, frame, side, copy=True):
        new_control_points = [p.frame_mapping(frame, side, True) for p in
                              self.control_points]
        new_bsplinesurface3d = BSplineSurface3D(self.degree_u, self.degree_v,
                                                new_control_points, self.nb_u,
                                                self.nb_v,
                                                self.u_multiplicities,
                                                self.v_multiplicities,
                                                self.u_knots, self.v_knots,
                                                self.weights, self.name)
        if copy:
            return new_bsplinesurface3d
        else:
            self.control_points = new_control_points
            self.surface = new_bsplinesurface3d.surface
            # self.points = new_BSplineSurface3D.points

    def plot(self, ax=None):
        for p in self.control_points:
            ax = p.plot(ax=ax)
        return ax

    @classmethod
    def from_step(cls, arguments, object_dict):
        name = arguments[0][1:-1]

        degree_u = int(arguments[1])
        degree_v = int(arguments[2])
        points_sets = arguments[3][1:-1].split("),")
        points_sets = [elem + ")" for elem in points_sets[:-1]] + [
            points_sets[-1]]
        control_points = []
        for points_set in points_sets:
            points = [object_dict[int(i[1:])] for i in
                      points_set[1:-1].split(",")]
            nb_v = len(points)
            control_points.extend(points)
        nb_u = int(len(control_points) / nb_v)
        surface_form = arguments[4]
        if arguments[5] == '.F.':
            u_closed = False
        elif arguments[5] == '.T.':
            u_closed = True
        else:
            raise ValueError
        if arguments[6] == '.F.':
            v_closed = False
        elif arguments[6] == '.T.':
            v_closed = True
        else:
            raise ValueError
        self_intersect = arguments[7]
        u_multiplicities = [int(i) for i in arguments[8][1:-1].split(",")]
        v_multiplicities = [int(i) for i in arguments[9][1:-1].split(",")]
        u_knots = [float(i) for i in arguments[10][1:-1].split(",")]
        v_knots = [float(i) for i in arguments[11][1:-1].split(",")]
        knot_spec = arguments[12]

        if 13 in range(len(arguments)):
            weight_data = [float(i) for i in
                           arguments[13][1:-1].replace("(", "").replace(")",
                                                                        "").split(
                               ",")]
        else:
            weight_data = None

        return cls(degree_u, degree_v, control_points, nb_u, nb_v,
                   u_multiplicities, v_multiplicities, u_knots, v_knots,
                   weight_data, name)


class BezierSurface3D(BSplineSurface3D):

    def __init__(self, degree_u: int, degree_v: int,
                 control_points: List[List[volmdlr.Point3D]],
                 nb_u: int, nb_v: int, name=''):

        u_knots = utilities.generate_knot_vector(degree_u, nb_u)
        v_knots = utilities.generate_knot_vector(degree_v, nb_v)

        u_multiplicities = [1] * len(u_knots)
        v_multiplicities = [1] * len(v_knots)

        BSplineSurface3D.__init__(self, degree_u, degree_v,
                                  control_points, nb_u, nb_v,
                                  u_multiplicities, v_multiplicities,
                                  u_knots, v_knots, None, name)


class Face3D(volmdlr.core.Primitive3D):
    min_x_density = 1
    min_y_density = 1

    def __init__(self, surface3d, surface2d: Surface2D,
                 name: str = ''):
        self.surface3d = surface3d
        self.surface2d = surface2d
        self.bounding_box = self._bounding_box()

        volmdlr.core.Primitive3D.__init__(self, name=name)

    def __hash__(self):
        return hash(self.surface3d) + hash(self.surface2d)

    def __eq__(self, other_):
        if other_.__class__.__name__ != self.__class__.__name__:
            return False
        equal = (self.surface3d == other_.surface3d
                 and self.surface2d == other_.surface2d)
        return equal

    def point_belongs(self, point3d: volmdlr.Point3D):
        """
        Tells you if a point is on the 3D face and inside its contour
        """
        point2d = self.surface3d.point3d_to_2d(point3d)
        check_point3d = self.surface3d.point2d_to_3d(point2d)
        if check_point3d.point_distance(point3d) > 1e-6:
            return False

        return self.surface2d.point_belongs(point2d)

    @property
    def outer_contour3d(self):
        """

        """
        return self.surface3d.contour2d_to_3d(self.surface2d.outer_contour)

    @property
    def inner_contours3d(self):
        """

        """
        return [self.surface3d.contour2d_to_3d(c) for c in
                self.surface2d.inner_contours]

    def _bounding_box(self):
        """
        this error is raised to enforce overloading of this method
        """
        raise NotImplementedError(
            '_bounding_box method must be overloaded by {}'.format(
                self.__class__.__name__))

    @classmethod
    def from_step(cls, arguments, object_dict):
        contours = [object_dict[int(arguments[1][0][1:])]]

        # Detecting inner and outer contours
        name = arguments[0][1:-1]
        surface = object_dict[int(arguments[2])]

        if hasattr(surface, 'face_from_contours3d'):
            if (len(contours) == 1) and isinstance(contours[0],
                                                   volmdlr.Point3D):
                return surface

            return surface.face_from_contours3d(contours, name)
        else:
            raise NotImplementedError(
                'Not implemented :face_from_contours3d in {}'.format(surface))

    def to_step(self, current_id):
        xmin, xmax, ymin, ymax = self.surface2d.bounding_rectangle()
        subsurfaces2d = [self.surface2d]
        line_x = None
        if self.surface3d.x_periodicity and (xmax - xmin) >= 0.45 * self.surface3d.x_periodicity:
            line_x = vme.Line2D(volmdlr.Point2D(0.5 * (xmin + xmax), 0),
                                volmdlr.Point2D(
                                    0.5 * (xmin + xmax), 1))
        line_y = None
        if self.surface3d.y_periodicity and (
                ymax - ymin) >= 0.45 * self.surface3d.y_periodicity:
            line_y = vme.Line2D(
                volmdlr.Point2D(0., 0.5 * (ymin + ymax)),
                volmdlr.Point2D(1, 0.5 * (ymin + ymax)))

        if line_x:
            subsurfaces2 = []
            for subsurface2d in subsurfaces2d:
                subsurfaces2.extend(subsurface2d.cut_by_line(line_x))
            subsurfaces2d = subsurfaces2

        if line_y:
            subsurfaces2 = []
            for subsurface2d in subsurfaces2d:
                subsurfaces2.extend(subsurface2d.cut_by_line(line_y))
            subsurfaces2d = subsurfaces2

        if len(subsurfaces2d) > 1:
            content = ''
            face_ids = []
            for subsurface2d in subsurfaces2d:
                face = self.__class__(self.surface3d, subsurface2d)
                try:
                    face_content, face_id = face.to_step_without_splitting(
                        current_id)
                    face_ids.append(face_id[0])
                    content += face_content
                    current_id = face_id[0] + 1
                except NotImplementedError:
                    print('Warning: a face of class {} has not been exported due to NotImplementedError'.format(
                        face.__class__.__name__))
                except AttributeError:
                    print(
                        'Warning: a face of class {} has not been exported due to AttributeError'.format(
                         face.__class__.__name__))
            return content, face_ids
        else:
            return self.to_step_without_splitting(current_id)

    def to_step_without_splitting(self, current_id):
        
        content, surface3d_id = self.surface3d.to_step(current_id)
        current_id = surface3d_id + 1

        outer_contour_content, outer_contour_id = self.outer_contour3d.to_step(
            current_id)
        # surface_id=surface3d_id)
        content += outer_contour_content
        content += "#{} = FACE_BOUND('{}',#{},.T.);\n".format(
            outer_contour_id + 1, self.name, outer_contour_id)
        contours_ids = [outer_contour_id + 1]
        current_id = outer_contour_id + 2
        for inner_contour3d in self.inner_contours3d:
            inner_contour_content, inner_contour_id = inner_contour3d.to_step(
                current_id)
            # surface_id=surface3d_id)
            content += inner_contour_content
            face_bound_id = inner_contour_id + 1
            content += "#{} = FACE_BOUND('',#{},.T.);\n".format(
                face_bound_id, inner_contour_id)
            contours_ids.append(face_bound_id)
            current_id = face_bound_id + 1

        content += "#{} = ADVANCED_FACE('{}',({}),#{},.T.);\n".format(
            current_id,
            self.name,
            volmdlr.core.step_ids_to_str(contours_ids),
            surface3d_id)
        return content, [current_id]

    def triangulation_lines(self):
        return [], []

    def triangulation(self):

        lines_x, lines_y = self.triangulation_lines()
        if lines_x and lines_y:
            surfaces = []
            for surface in self.surface2d.split_by_lines(lines_x):
                surfaces.extend(surface.split_by_lines(lines_y))

        elif lines_x:
            surfaces = self.surface2d.split_by_lines(lines_x)
        elif lines_y:
            surfaces = self.surface2d.split_by_lines(lines_y)
        else:
            surfaces = [self.surface2d]

        mesh2d = surfaces[0].triangulation()
        for subsurface in surfaces[1:]:
            mesh2d += subsurface.triangulation()

        return volmdlr.display.DisplayMesh3D(
            [volmdlr.display.Node3D(*self.surface3d.point2d_to_3d(p)) for p in
             mesh2d.points],
            mesh2d.triangles)

    def plot2d(self, ax=None, color='k', alpha=1):
        if ax is None:
            _, ax = plt.subplots()

        self.outer_contour.plot()

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            new_surface = self.surface3d.rotation(center=center, axis=axis,
                                                  angle=angle, copy=True)
            return self.__class__(new_surface, self.surface2d)
        else:
            self.surface3d.rotation(center=center, axis=axis,
                                    angle=angle, copy=False)
            self.bounding_box = self._bounding_box()

    def translation(self, offset, copy=True):
        if copy:
            new_surface3d = self.surface3d.translation(offset=offset,
                                                       copy=True)
            return self.__class__(new_surface3d, self.surface2d)
        else:
            self.surface3d.translation(offset=offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_surface = self.surface3d.frame_mapping(frame, side, copy=True)
            return self.__class__(new_surface, self.surface2d.copy(),
                                  self.name)
        else:
            self.surface3d.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()

    def copy(self):
        return Face3D(self.surface3d.copy(), self.surface2d.copy(), self.name)

    def linesegment_intersections(self,
                                  linesegment: vme.LineSegment3D,
                                  ) -> List[volmdlr.Point3D]:
        intersections = []
        for intersection in self.surface3d.linesegment_intersections(
                linesegment):
            if self.point_belongs(intersection):
                intersections.append(intersection)

        return intersections

    def plot(self, ax=None, color='k', alpha=1):
        if not ax:
            ax = plt.figure().add_subplot(111, projection='3d')
        self.outer_contour3d.plot(ax=ax, color=color, alpha=alpha)
        return ax


class PlaneFace3D(Face3D):
    """
    :param contours: The face's contour2D
    :type contours: volmdlr.Contour2D
    :param plane: Plane used to place your face
    :type plane: Plane3D
    """
    _standalone_in_db = False
    _generic_eq = True
    _non_serializable_attributes = ['bounding_box', 'polygon2D']
    _non_eq_attributes = ['name', 'bounding_box', 'outer_contour3d',
                          'inner_contours3d']
    _non_hash_attributes = []

    def __init__(self, surface3d: Plane3D, surface2d: Surface2D,
                 name: str = ''):
        # if not isinstance(outer_contour2d, volmdlr.Contour2D):
        #     raise ValueError('Not a contour2D: {}'.format(outer_contour2d))
        Face3D.__init__(self,
                        surface3d=surface3d,
                        surface2d=surface2d,
                        name=name)

    # @classmethod
    # def _repair_points_and_polygon2d(cls, points, plane):
    #     if points[0] == points[-1]:
    #         points = points[:-1]
    #     polygon_points = [
    #         p.to_2d(plane.origin, plane.vectors[0], plane.vectors[1]) for p in
    #         points]
    #     repaired_points = [p.copy() for p in points]
    #     polygon2D = volmdlr.ClosedPolygon2D(polygon_points)
    #     if polygon2D.SelfIntersect()[0]:
    #         repaired_points = [repaired_points[1]] + [
    #             repaired_points[0]] + repaired_points[2:]
    #         polygon_points = [polygon_points[1]] + [
    #             polygon_points[0]] + polygon_points[2:]
    #         if polygon_points[0] == polygon_points[-1]:
    #             repaired_points = repaired_points[:-1]
    #             polygon_points = polygon_points[:-1]
    #         polygon2D = volmdlr.ClosedPolygon2D(polygon_points)
    #     return repaired_points, polygon2D

    @classmethod
    def dict_to_object(cls, dict_):
        plane3d = Plane3D.dict_to_object(dict_['surface3d'])
        surface2d = Surface2D.dict_to_object(dict_['surface2d'])
        return cls(plane3d, surface2d, dict_['name'])

    def copy(self):
        return PlaneFace3D(self.surface3d.copy(), self.surface2d.copy(),
                           self.name)

    def _bounding_box(self):
        """
        """
        return self.outer_contour3d._bounding_box()

    # def average_center_point(self):
    #     """
    #     excluding holes
    #     """
    #     points = self.points
    #     nb = len(points)
    #     x = npy.sum([p[0] for p in points]) / nb
    #     y = npy.sum([p[1] for p in points]) / nb
    #     z = npy.sum([p[2] for p in points]) / nb
    #     return volmdlr.Point3D((x, y, z))

    def distance_to_point(self, point, return_other_point=False):
        # """
        # Only works if the surface is planar
        # TODO : this function does not take into account if Face has holes
        # """
        # On projette le point sur la surface plane
        # Si le point est à l'intérieur de la face,
        # on retourne la distance de projection
        # Si le point est à l'extérieur, on projette le point sur le plan
        # On calcule en 2D la distance entre la projection
        # et le polygone contour
        # On utilise le theroeme de Pythagore pour calculer
        # la distance minimale entre le point et le contour

        projected_pt = point.PlaneProjection3D(self.plane.origin,
                                               self.plane.vectors[0],
                                               self.plane.vectors[1])
        projection_distance = point.point_distance(projected_pt)

        if self.point_on_face(projected_pt):
            if return_other_point:
                return projection_distance, projected_pt
            return projection_distance

        point_2D = point.to_2d(self.plane.origin, self.plane.vectors[0],
                               self.plane.vectors[1])

        border_distance, other_point = self.polygon2D.PointBorderDistance(
            point_2D, return_other_point=True)

        other_point = other_point.to_3d(self.plane.origin,
                                        self.plane.vectors[0],
                                        self.plane.vectors[1])

        if return_other_point:
            return (projection_distance ** 2 + border_distance ** 2) ** 0.5, \
                   other_point
        return (projection_distance ** 2 + border_distance ** 2) ** 0.5

    def minimum_distance_points_plane(self, other_plane_face,
                                      return_points=False):
        ## """
        ## Only works if the surface is planar
        ## TODO : this function does not take into account if Face has holes
        ## TODO : TRAITER LE CAS OU LA DISTANCE LA PLUS COURTE N'EST PAS D'UN SOMMET
        ## """
        # On calcule la distance entre la face 1 et chaque point de la face 2
        # On calcule la distance entre la face 2 et chaque point de la face 1

        # if self.face_intersection(other_plane_face) is not None:
        #     return 0, None, None
        #
        # polygon1_points_3D = [volmdlr.Point3D(p.vector) for p in
        #                       self.contours3d[0].tessel_points]
        # polygon2_points_3D = [volmdlr.Point3D(p.vector) for p in
        #                       other_plane_face.contours3d[0].tessel_points]
        #
        # distances = []
        # if not return_points:
        #     d_min = other_plane_face.distance_to_point(polygon1_points_3D[0])
        #     for point1 in polygon1_points_3D[1:]:
        #         d = other_plane_face.distance_to_point(point1)
        #         if d < d_min:
        #             d_min = d
        #     for point2 in polygon2_points_3D:
        #         d = self.distance_to_point(point2)
        #         if d < d_min:
        #             d_min = d
        #     return d_min
        #
        # else:
        #     for point1 in polygon1_points_3D:
        #         d, other_point = other_plane_face.distance_to_point(
        #             point1,
        #             return_other_point=True)
        #         distances.append((d, point1, other_point))
        #     for point2 in polygon2_points_3D:
        #         d, other_point = self.distance_to_point(
        #             point2,
        #             return_other_point=True
        #         )
        #         distances.append((d, point2, other_point))
        #
        # d_min, point_min, other_point_min = distances[0]
        # for distance in distances[1:]:
        #     if distance[0] < d_min:
        #         d_min = distance[0]
        #         point_min = distance[1]
        #         other_point_min = distance[2]
        #
        # return point_min, other_point_min

        min_distance = math.inf
        for edge1 in self.outer_contour3d.primitives:
            for edge2 in other_plane_face.outer_contour3d.primitives:
                dist = edge1.minimum_distance(edge2,
                                              return_points=return_points)
                if return_points:
                    if dist[0] < min_distance:
                        min_distance = dist[0]
                        p1, p2 = dist[1], dist[2]
                else:
                    if dist < min_distance:
                        min_distance = dist
        if return_points:
            return min_distance, p1, p2
        else:
            return min_distance

    def edge_intersections(self, edge):
        intersections = []
        linesegment = vme.LineSegment3D(edge.start, edge.end)
        for surface3d_inter in self.surface3d.linesegment_intersections(linesegment):
            point2d = self.surface3d.point3d_to_2d(surface3d_inter)
            if self.surface2d.point_belongs(point2d):
                intersections.append(surface3d_inter)

        return intersections

    def face_intersections(self, face2):
        ## """
        ## Only works if the surface is planar
        ## TODO : this function does not take into account if Face has holes
        ## """
        bbox1 = self.bounding_box
        bbox2 = face2.bounding_box
        if not bbox1.bbox_intersection(bbox2):
            return []

        intersections = []

        for edge2 in face2.outer_contour3d.primitives:
            intersection_points = self.edge_intersections(edge2)
            if intersection_points:
                intersections.extend(intersection_points)

        for edge1 in self.outer_contour3d.primitives:
            intersection_points = face2.edge_intersections(edge1)
            if intersection_points:
                intersections.extend(intersection_points)

        return intersections

    def minimum_distance(self, other_face, return_points=False):
        if other_face.__class__ is CylindricalFace3D:
            p1, p2 = other_face.minimum_distance_points_cyl(self)
            if return_points:
                return p1.point_distance(p2), p1, p2
            else:
                return p1.point_distance(p2)

        if other_face.__class__ is PlaneFace3D:
            if return_points:
                dist, p1, p2 = self.minimum_distance_points_plane(other_face,
                                                                  return_points=return_points)
                return dist, p1, p2
            else:
                dist = self.minimum_distance_points_plane(other_face,
                                                          return_points=return_points)
                return dist

        if other_face.__class__ is ToroidalFace3D:
            p1, p2 = other_face.minimum_distance_points_plane(self)
            if return_points:
                return p1.point_distance(p2), p1, p2
            else:
                return p1.point_distance(p2)

        else:
            return NotImplementedError


class CylindricalFace3D(Face3D):
    """
    :param contours2d: The cylinder's contour2D
    :type contours2d: volmdlr.Contour2D
    :param cylindricalsurface3d: Information about the Cylinder
    :type cylindricalsurface3d: CylindricalSurface3D
    :param points: contours2d's point
    :type points: List of volmdlr.Point2D

    :Example:
        >>> contours2d is rectangular and will create a classic cylinder with x= 2*pi*radius, y=h
    """
    min_x_density = 5
    min_y_density = 1

    def __init__(self,
                 cylindricalsurface3d: CylindricalSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):

        self.radius = cylindricalsurface3d.radius
        self.center = cylindricalsurface3d.frame.origin
        self.normal = cylindricalsurface3d.frame.w
        Face3D.__init__(self, surface3d=cylindricalsurface3d,
                        surface2d=surface2d,
                        name=name)

    def _bounding_box(self):
        theta_min, theta_max, zmin, zmax = self.surface2d.outer_contour.bounding_rectangle()

        xp = (volmdlr.X3D.dot(self.surface3d.frame.u) * self.surface3d.frame.u
              + volmdlr.X3D.dot(
                    self.surface3d.frame.v) * self.surface3d.frame.v)
        xp_norm = xp.norm()
        if xp_norm != 0:
            xp = xp / xp_norm

        yp = (volmdlr.Y3D.dot(self.surface3d.frame.u) * self.surface3d.frame.u
              + volmdlr.Y3D.dot(
                    self.surface3d.frame.v) * self.surface3d.frame.v)
        yp_norm = yp.norm()
        if yp_norm != 0:
            yp = yp / yp_norm

        zp = (volmdlr.Z3D.dot(self.surface3d.frame.u) * self.surface3d.frame.u
              + volmdlr.Z3D.dot(
                    self.surface3d.frame.v) * self.surface3d.frame.v)
        zp_norm = zp.norm()
        if zp_norm != 0:
            zp = zp / zp_norm

        lower_center = self.surface3d.frame.origin + zmin * self.surface3d.frame.w
        upper_center = self.surface3d.frame.origin + zmax * self.surface3d.frame.w

        points = [lower_center - self.surface3d.radius * xp,
                  lower_center + self.surface3d.radius * xp,
                  lower_center - self.surface3d.radius * yp,
                  lower_center + self.surface3d.radius * yp,
                  lower_center - self.surface3d.radius * zp,
                  lower_center + self.surface3d.radius * zp,
                  upper_center - self.surface3d.radius * xp,
                  upper_center + self.surface3d.radius * xp,
                  upper_center - self.surface3d.radius * yp,
                  upper_center + self.surface3d.radius * yp,
                  upper_center - self.surface3d.radius * zp,
                  upper_center + self.surface3d.radius * zp,
                  ]

        return volmdlr.core.BoundingBox.from_points(points)

    def triangulation_lines(self, angle_resolution=5):
        theta_min, theta_max, zmin, zmax = self.surface2d.bounding_rectangle()
        delta_theta = theta_max - theta_min
        nlines = math.ceil(delta_theta * angle_resolution)
        lines = []
        for i in range(nlines):
            theta = theta_min + (i + 1) / (nlines + 1) * delta_theta
            lines.append(vme.Line2D(volmdlr.Point2D(theta, zmin),
                                    volmdlr.Point2D(theta, zmax)))
        return lines, []

    def range_closest(list_point, radius, frame):
        points_set = volmdlr.delete_double_point(list_point)
        points_set3D = CylindricalFace3D.points2d_to3d(None, [points_set],
                                                       radius, frame)

        points_3dint = [points_set3D[0]]
        points_2dint = [points_set[0]]
        s = 1
        for k in range(1, len(points_set)):
            closest = points_set3D[s]
            while closest is None:
                s += 1
                closest = points_set3D[s]
            dist_min = (points_3dint[-1] - closest).norm()
            pos = s
            for i in range(s + 1, len(points_set3D)):
                close_test = points_set3D[i]
                if close_test is None:
                    continue
                else:
                    dist_test = (points_3dint[-1] - close_test).norm()
                    if dist_test <= dist_min:
                        dist_min = dist_test
                        closest = close_test
                        pos = i
            points_2dint.append(points_set[pos])
            points_set3D[pos] = None

        return points_2dint

    # def frame_mapping(self, frame, side, copy=True):
    #     if copy:
    #         new_cylindricalsurface3d = CylindricalSurface3D.frame_mapping(
    #             frame, side, copy)
    #         return CylindricalFace3D(self.contours2d, new_cylindricalsurface3d,
    #                                  points=self.points, name=self.name)
    #     else:
    #         self.cylindricalsurface3d.frame_mapping(frame, side, copy=False)

    def minimum_maximum(self, contour2d, radius):
        points = contour2d.tessel_points

        min_h, min_theta = min([pt[1] for pt in points]), min(
            [pt[0] for pt in points])
        max_h, max_theta = max([pt[1] for pt in points]), max(
            [pt[0] for pt in points])
        return min_h, min_theta, max_h, max_theta

    def minimum_distance_points_cyl(self, other_cyl):
        r1, r2 = self.radius, other_cyl.radius
        min_h1, min_theta1, max_h1, max_theta1 = self.minimum_maximum(
            self.contours2d[0], r1)

        n1 = self.normal
        u1 = self.cylindricalsurface3d.frame.u
        v1 = self.cylindricalsurface3d.frame.v
        frame1 = volmdlr.Frame3D(self.center, u1, v1, n1)

        min_h2, min_theta2, max_h2, max_theta2 = self.minimum_maximum(
            other_cyl.contours2d[0], r2)

        n2 = other_cyl.normal
        u2 = other_cyl.cylindricalsurface3d.frame.u
        v2 = other_cyl.cylindricalsurface3d.frame.v
        frame2 = volmdlr.Frame3D(other_cyl.center, u2, v2, n2)
        # st2 = volmdlr.Point3D((r2*math.cos(min_theta2), r2*math.sin(min_theta2), min_h2))
        # start2 = frame2.old_coordinates(st2)

        w = other_cyl.center - self.center

        n1n1, n1u1, n1v1, n1n2, n1u2, n1v2 = n1.dot(n1), n1.dot(u1), n1.dot(
            v1), n1.dot(n2), n1.dot(u2), n1.dot(v2)
        u1u1, u1v1, u1n2, u1u2, u1v2 = u1.dot(u1), u1.dot(v1), u1.dot(
            n2), u1.dot(u2), u1.dot(v2)
        v1v1, v1n2, v1u2, v1v2 = v1.dot(v1), v1.dot(n2), v1.dot(u2), v1.dot(v2)
        n2n2, n2u2, n2v2 = n2.dot(n2), n2.dot(u2), n2.dot(v2)
        u2u2, u2v2, v2v2 = u2.dot(u2), u2.dot(v2), v2.dot(v2)

        w2, wn1, wu1, wv1, wn2, wu2, wv2 = w.dot(w), w.dot(n1), w.dot(
            u1), w.dot(v1), w.dot(n2), w.dot(u2), w.dot(v2)

        # x = (theta1, h1, theta2, h2)
        def distance_squared(x):
            return (n1n1 * (x[1] ** 2) + u1u1 * ((math.cos(x[0])) ** 2) * (
                    r1 ** 2) + v1v1 * ((math.sin(x[0])) ** 2) * (r1 ** 2)
                    + w2 + n2n2 * (x[3] ** 2) + u2u2 * (
                            (math.cos(x[2])) ** 2) * (r2 ** 2) + v2v2 * (
                            (math.sin(x[2])) ** 2) * (r2 ** 2)
                    + 2 * x[1] * r1 * math.cos(x[0]) * n1u1 + 2 * x[
                        1] * r1 * math.sin(x[0]) * n1v1 - 2 * x[1] * wn1
                    - 2 * x[1] * x[3] * n1n2 - 2 * x[1] * r2 * math.cos(
                        x[2]) * n1u2 - 2 * x[1] * r2 * math.sin(x[2]) * n1v2
                    + 2 * math.cos(x[0]) * math.sin(x[0]) * u1v1 * (
                            r1 ** 2) - 2 * r1 * math.cos(x[0]) * wu1
                    - 2 * r1 * x[3] * math.cos(
                        x[0]) * u1n2 - 2 * r1 * r2 * math.cos(x[0]) * math.cos(
                        x[2]) * u1u2
                    - 2 * r1 * r2 * math.cos(x[0]) * math.sin(
                        x[2]) * u1v2 - 2 * r1 * math.sin(x[0]) * wv1
                    - 2 * r1 * x[3] * math.sin(
                        x[0]) * v1n2 - 2 * r1 * r2 * math.sin(x[0]) * math.cos(
                        x[2]) * v1u2
                    - 2 * r1 * r2 * math.sin(x[0]) * math.sin(
                        x[2]) * v1v2 + 2 * x[3] * wn2 + 2 * r2 * math.cos(
                        x[2]) * wu2
                    + 2 * r2 * math.sin(x[2]) * wv2 + 2 * x[3] * r2 * math.cos(
                        x[2]) * n2u2 + 2 * x[3] * r2 * math.sin(x[2]) * n2v2
                    + 2 * math.cos(x[2]) * math.sin(x[2]) * u2v2 * (r2 ** 2))

        x01 = npy.array([(min_theta1 + max_theta1) / 2, (min_h1 + max_h1) / 2,
                         (min_theta2 + max_theta2) / 2, (min_h2 + max_h2) / 2])
        x02 = npy.array([min_theta1, (min_h1 + max_h1) / 2,
                         min_theta2, (min_h2 + max_h2) / 2])
        x03 = npy.array([max_theta1, (min_h1 + max_h1) / 2,
                         max_theta2, (min_h2 + max_h2) / 2])

        minimax = [(min_theta1, min_h1, min_theta2, min_h2),
                   (max_theta1, max_h1, max_theta2, max_h2)]

        res1 = scp.optimize.least_squares(distance_squared, x01,
                                          bounds=minimax)
        res2 = scp.optimize.least_squares(distance_squared, x02,
                                          bounds=minimax)
        res3 = scp.optimize.least_squares(distance_squared, x03,
                                          bounds=minimax)

        pt1 = volmdlr.Point3D(
            (r1 * math.cos(res1.x[0]), r1 * math.sin(res1.x[0]), res1.x[1]))
        p1 = frame1.old_coordinates(pt1)
        pt2 = volmdlr.Point3D(
            (r2 * math.cos(res1.x[2]), r2 * math.sin(res1.x[2]), res1.x[3]))
        p2 = frame2.old_coordinates(pt2)
        d = p1.point_distance(p2)
        result = res1

        res = [res2, res3]
        for couple in res:
            pttest1 = volmdlr.Point3D((r1 * math.cos(couple.x[0]),
                                       r1 * math.sin(couple.x[0]),
                                       couple.x[1]))
            pttest2 = volmdlr.Point3D((r2 * math.cos(couple.x[2]),
                                       r2 * math.sin(couple.x[2]),
                                       couple.x[3]))
            ptest1 = frame1.old_coordinates(pttest1)
            ptest2 = frame2.old_coordinates(pttest2)
            dtest = ptest1.point_distance(ptest2)
            if dtest < d:
                result = couple
                p1, p2 = ptest1, ptest2

        pt1_2d, pt2_2d = volmdlr.Point2D(
            (result.x[0], result.x[1])), volmdlr.Point2D(
            (result.x[2], result.x[3]))

        if not (self.contours2d[0].point_belongs(pt1_2d)):
            # Find the closest one
            points_contours1 = self.contours2d[0].tessel_points

            poly1 = volmdlr.ClosedPolygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d,
                                                       return_other_point=True)
            pt1 = volmdlr.Point3D((r1 * math.cos(new_pt1_2d.vector[0]),
                                   r1 * math.sin(new_pt1_2d.vector[0]),
                                   new_pt1_2d.vector[1]))
            p1 = frame1.old_coordinates(pt1)

        if not (other_cyl.contours2d[0].point_belongs(pt2_2d)):
            # Find the closest one
            points_contours2 = other_cyl.contours2d[0].tessel_points

            poly2 = volmdlr.ClosedPolygon2D(points_contours2)
            d2, new_pt2_2d = poly2.PointBorderDistance(pt2_2d,
                                                       return_other_point=True)
            pt2 = volmdlr.Point3D((r2 * math.cos(new_pt2_2d.vector[0]),
                                   r2 * math.sin(new_pt2_2d.vector[0]),
                                   new_pt2_2d.vector[1]))
            p2 = frame2.old_coordinates(pt2)

        return p1, p2

    def minimum_distance_points_plane(self,
                                      planeface):  # Planeface with contour2D
        #### ADD THE FACT THAT PLANEFACE.CONTOURS : [0] = contours totale, le reste = trous
        r = self.radius
        min_h1, min_theta1, max_h1, max_theta1 = self.minimum_maximum(
            self.contours2d[0], r)

        n1 = self.normal
        u1 = self.cylindricalsurface3d.frame.u
        v1 = self.cylindricalsurface3d.frame.v
        frame1 = volmdlr.Frame3D(self.center, u1, v1, n1)
        # st1 = volmdlr.Point3D((r*math.cos(min_theta1), r*math.sin(min_theta1), min_h1))
        # start1 = frame1.old_coordinates(st1)

        poly2d = planeface.polygon2D
        pfpoints = poly2d.points
        xmin, ymin = min([pt[0] for pt in pfpoints]), min(
            [pt[1] for pt in pfpoints])
        xmax, ymax = max([pt[0] for pt in pfpoints]), max(
            [pt[1] for pt in pfpoints])
        origin, vx, vy = planeface.plane.origin, planeface.plane.vectors[0], \
                         planeface.plane.vectors[1]
        pf1_2d, pf2_2d = volmdlr.Point2D((xmin, ymin)), volmdlr.Point2D(
            (xmin, ymax))
        pf3_2d, pf4_2d = volmdlr.Point2D((xmax, ymin)), volmdlr.Point2D(
            (xmax, ymax))
        pf1, pf2 = pf1_2d.to_3d(origin, vx, vy), pf2_2d.to_3d(origin, vx, vy)
        pf3, _ = pf3_2d.to_3d(origin, vx, vy), pf4_2d.to_3d(origin, vx, vy)

        u, v = (pf3 - pf1), (pf2 - pf1)
        u.normalize()
        v.normalize()

        w = pf1 - self.center

        n1n1, n1u1, n1v1, n1u, n1v = n1.dot(n1), n1.dot(u1), n1.dot(
            v1), n1.dot(u), n1.dot(v)
        u1u1, u1v1, u1u, u1v = u1.dot(u1), u1.dot(v1), u1.dot(u), u1.dot(v)
        v1v1, v1u, v1v = v1.dot(v1), v1.dot(u), v1.dot(v)
        uu, uv, vv = u.dot(u), u.dot(v), v.dot(v)

        w2, wn1, wu1, wv1, wu, wv = w.dot(w), w.dot(n1), w.dot(u1), w.dot(
            v1), w.dot(u), w.dot(v)

        # x = (h, theta, x, y)
        def distance_squared(x):
            return (n1n1 * (x[0] ** 2) + ((math.cos(x[1])) ** 2) * u1u1 * (
                    r ** 2) + ((math.sin(x[1])) ** 2) * v1v1 * (r ** 2)
                    + w2 + uu * (x[2] ** 2) + vv * (x[3] ** 2) + 2 * x[
                        0] * math.cos(x[1]) * r * n1u1
                    + 2 * x[0] * math.sin(x[1]) * r * n1v1 - 2 * x[
                        0] * wn1 - 2 * x[0] * x[2] * n1u
                    - 2 * x[0] * x[3] * n1v + 2 * math.sin(x[1]) * math.cos(
                        x[1]) * u1v1 * (r ** 2)
                    - 2 * r * math.cos(x[1]) * wu1 - 2 * r * x[2] * math.cos(
                        x[1]) * u1u
                    - 2 * r * x[3] * math.sin(x[1]) * u1v - 2 * r * math.sin(
                        x[1]) * wv1
                    - 2 * r * x[2] * math.sin(x[1]) * v1u - 2 * r * x[
                        3] * math.sin(x[1]) * v1v
                    + 2 * x[2] * wu + 2 * x[3] * wv + 2 * x[2] * x[3] * uv)

        x01 = npy.array([(min_h1 + max_h1) / 2, (min_theta1 + max_theta1) / 2,
                         (xmax - xmin) / 2, (ymax - ymin) / 2])

        minimax = [(min_h1, min_theta1, 0, 0),
                   (max_h1, max_theta1, xmax - xmin, ymax - ymin)]

        res1 = scp.optimize.least_squares(distance_squared, x01,
                                          bounds=minimax)

        pt1 = volmdlr.Point3D(
            (r * math.cos(res1.x[1]), r * math.sin(res1.x[1]), res1.x[0]))
        p1 = frame1.old_coordinates(pt1)
        p2 = pf1 + res1.x[2] * u + res1.x[3] * v
        pt1_2d = volmdlr.Point2D((res1.x[1], res1.x[0]))
        pt2_2d = p2.to_2d(pf1, u, v)

        if not (self.contours2d[0].point_belongs(pt1_2d)):
            # Find the closest one
            points_contours1 = self.contours2d[0].tessel_points

            poly1 = volmdlr.ClosedPolygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d,
                                                       return_other_point=True)
            pt1 = volmdlr.Point3D((r * math.cos(new_pt1_2d.vector[0]),
                                   r * math.sin(new_pt1_2d.vector[0]),
                                   new_pt1_2d.vector[1]))
            p1 = frame1.old_coordinates(pt1)

        if not (planeface.contours[0].point_belongs(pt2_2d)):
            # Find the closest one
            d2, new_pt2_2d = planeface.polygon2D.PointBorderDistance(pt2_2d,
                                                                     return_other_point=True)

            p2 = new_pt2_2d.to_3d(pf1, u, v)

        return p1, p2

    def minimum_distance(self, other_face, return_points=False):
        if other_face.__class__ is CylindricalFace3D:
            p1, p2 = self.minimum_distance_points_cyl(other_face)
            if return_points:
                return p1.point_distance(p2), p1, p2
            else:
                return p1.point_distance(p2)

        if other_face.__class__ is PlaneFace3D:
            p1, p2 = self.minimum_distance_points_plane(other_face)
            if return_points:
                return p1.point_distance(p2), p1, p2
            else:
                return p1.point_distance(p2)

        if other_face.__class__ is ToroidalFace3D:
            p1, p2 = other_face.minimum_distance_points_cyl(self)
            if return_points:
                return p1.point_distance(p2), p1, p2
            else:
                return p1.point_distance(p2)

        else:
            return NotImplementedError


class ToroidalFace3D(Face3D):
    """
    :param contours2d: The Tore's contour2D
    :type contours2d: volmdlr.Contour2D
    :param toroidalsurface3d: Information about the Tore
    :type toroidalsurface3d: ToroidalSurface3D
    :param theta: angle of cut in main circle direction
    :param phi: angle of cut in secondary circle direction
    :type points: List of float

    Example
        contours2d is rectangular and will create a classic tore with x:2*pi, y:2*pi
        x is for exterior, and y for the circle to revolute
        points = [pi, 2*pi] for an half tore
    """
    min_x_density = 5
    min_y_density = 1

    def __init__(self, toroidalsurface3d: ToroidalSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):

        # self.toroidalsurface3d = toroidalsurface3d

        self.center = toroidalsurface3d.frame.origin
        self.normal = toroidalsurface3d.frame.w

        theta_min, theta_max, phi_min, phi_max = surface2d.outer_contour.bounding_rectangle()

        self.theta_min = theta_min
        self.theta_max = theta_max
        self.phi_min = phi_min
        self.phi_max = phi_max

        # contours3d = [self.toroidalsurface3d.contour2d_to_3d(c)\
        #               for c in [outer_contour2d]+inners_contours2d]

        Face3D.__init__(self,
                        surface3d=toroidalsurface3d,
                        surface2d=surface2d,
                        name=name)

    def points_resolution(self, line, pos,
                          resolution):  # With a resolution wished
        points = []
        points.append(line.points[0])
        limit = line.points[1].vector[pos]
        start = line.points[0].vector[pos]
        vec = [0, 0]
        vec[pos] = start
        echelon = [line.points[0].vector[0] - vec[0],
                   line.points[0].vector[1] - vec[1]]
        flag = start + resolution
        while flag < limit:
            echelon[pos] = flag
            flag += resolution
            points.append(volmdlr.Point2D(echelon))
        points.append(line.points[1])
        return points

    def _bounding_box(self):
        return self.surface3d._bounding_box()

    def triangulation_lines(self, angle_resolution=5):
        theta_min, theta_max, phi_min, phi_max = self.surface2d.bounding_rectangle()

        delta_theta = theta_max - theta_min
        nlines_x = int(delta_theta * angle_resolution)
        lines_x = []
        for i in range(nlines_x):
            theta = theta_min + (i + 1) / (nlines_x + 1) * delta_theta
            lines_x.append(vme.Line2D(volmdlr.Point2D(theta, phi_min),
                                      volmdlr.Point2D(theta, phi_max)))
        delta_phi = phi_max - phi_min
        nlines_y = int(delta_phi * angle_resolution)
        lines_y = []
        for i in range(nlines_y):
            phi = phi_min + (i + 1) / (nlines_y + 1) * delta_phi
            lines_y.append(vme.Line2D(volmdlr.Point2D(theta_min, phi),
                                      volmdlr.Point2D(theta_max, phi)))
        return lines_x, lines_y

    def minimum_maximum_tore(self, contour2d):
        points = contour2d.tessel_points

        min_phi, min_theta = min([pt[1] for pt in points]), min(
            [pt[0] for pt in points])
        max_phi, max_theta = max([pt[1] for pt in points]), max(
            [pt[0] for pt in points])
        return min_phi, min_theta, max_phi, max_theta

    def minimum_distance_points_tore(self, other_tore):
        R1, r1, R2, r2 = self.rcenter, self.rcircle, other_tore.rcenter, other_tore.rcircle

        min_phi1, min_theta1, max_phi1, max_theta1 = self.minimum_maximum_tore(
            self.contours2d[0])

        # start1 = self.start
        n1 = self.normal
        u1 = self.toroidalsurface3d.frame.u
        v1 = self.toroidalsurface3d.frame.v
        frame1 = volmdlr.Frame3D(self.center, u1, v1, n1)
        # start1 = self.points2d_to3d([[min_theta1, min_phi1]], R1, r1, frame1)

        min_phi2, min_theta2, max_phi2, max_theta2 = self.minimum_maximum_tore(
            other_tore.contours2d[0])

        # start2 = other_tore.start
        n2 = other_tore.normal
        u2 = other_tore.toroidalsurface3d.frame.u
        v2 = other_tore.toroidalsurface3d.frame.v
        frame2 = volmdlr.Frame3D(other_tore.center, u2, v2, n2)
        # start2 = other_tore.points2d_to3d([[min_theta2, min_phi2]], R2, r2, frame2)

        w = other_tore.center - self.center

        n1n1, n1u1, n1v1, n1n2, n1u2, n1v2 = n1.dot(n1), n1.dot(u1), n1.dot(
            v1), n1.dot(n2), n1.dot(u2), n1.dot(v2)
        u1u1, u1v1, u1n2, u1u2, u1v2 = u1.dot(u1), u1.dot(v1), u1.dot(
            n2), u1.dot(u2), u1.dot(v2)
        v1v1, v1n2, v1u2, v1v2 = v1.dot(v1), v1.dot(n2), v1.dot(u2), v1.dot(v2)
        n2n2, n2u2, n2v2 = n2.dot(n2), n2.dot(u2), n2.dot(v2)
        u2u2, u2v2, v2v2 = u2.dot(u2), u2.dot(v2), v2.dot(v2)

        w2, wn1, wu1, wv1, wn2, wu2, wv2 = w.dot(w), w.dot(n1), w.dot(
            u1), w.dot(v1), w.dot(n2), w.dot(u2), w.dot(v2)

        # x = (phi1, theta1, phi2, theta2)
        def distance_squared(x):
            return (u1u1 * (((R1 + r1 * math.cos(x[0])) * math.cos(x[1])) ** 2)
                    + v1v1 * (((R1 + r1 * math.cos(x[0])) * math.sin(
                        x[1])) ** 2)
                    + n1n1 * ((math.sin(x[0])) ** 2) * (r1 ** 2) + w2
                    + u2u2 * (((R2 + r2 * math.cos(x[2])) * math.cos(
                        x[3])) ** 2)
                    + v2v2 * (((R2 + r2 * math.cos(x[2])) * math.sin(
                        x[3])) ** 2)
                    + n2n2 * ((math.sin(x[2])) ** 2) * (r2 ** 2)
                    + 2 * u1v1 * math.cos(x[1]) * math.sin(x[1]) * (
                            (R1 + r1 * math.cos(x[0])) ** 2)
                    + 2 * (R1 + r1 * math.cos(x[0])) * math.cos(
                        x[1]) * r1 * math.sin(x[0]) * n1u1
                    - 2 * (R1 + r1 * math.cos(x[0])) * math.cos(x[1]) * wu1
                    - 2 * (R1 + r1 * math.cos(x[0])) * (
                            R2 + r2 * math.cos(x[2])) * math.cos(
                        x[1]) * math.cos(x[3]) * u1u2
                    - 2 * (R1 + r1 * math.cos(x[0])) * (
                            R2 + r2 * math.cos(x[2])) * math.cos(
                        x[1]) * math.sin(x[3]) * u1v2
                    - 2 * (R1 + r1 * math.cos(x[0])) * math.cos(
                        x[1]) * r2 * math.sin(x[2]) * u1n2
                    + 2 * (R1 + r1 * math.cos(x[0])) * math.sin(
                        x[1]) * r1 * math.sin(x[0]) * n1v1
                    - 2 * (R1 + r1 * math.cos(x[0])) * math.sin(x[1]) * wv1
                    - 2 * (R1 + r1 * math.cos(x[0])) * (
                            R2 + r2 * math.cos(x[2])) * math.sin(
                        x[1]) * math.cos(x[3]) * v1u2
                    - 2 * (R1 + r1 * math.cos(x[0])) * (
                            R2 + r2 * math.cos(x[2])) * math.sin(
                        x[1]) * math.sin(x[3]) * v1v2
                    - 2 * (R1 + r1 * math.cos(x[0])) * math.sin(
                        x[1]) * r2 * math.sin(x[2]) * v1n2
                    - 2 * r1 * math.sin(x[0]) * wn1
                    - 2 * r1 * math.sin(x[0]) * (
                            R2 + r2 * math.cos(x[2])) * math.cos(
                        x[3]) * n1u2
                    - 2 * r1 * math.sin(x[0]) * (
                            R2 + r2 * math.cos(x[2])) * math.sin(
                        x[3]) * n1v2
                    - 2 * r1 * r2 * math.sin(x[0]) * math.sin(x[2]) * n1n2
                    + 2 * (R2 + r2 * math.cos(x[2])) * math.cos(x[3]) * wu2
                    + 2 * (R2 + r2 * math.cos(x[2])) * math.sin(x[3]) * wv2
                    + 2 * r2 * math.sin(x[2]) * wn2
                    + 2 * u2v2 * math.cos(x[3]) * math.sin(x[3]) * (
                            (R2 + r2 * math.cos(x[2])) ** 2)
                    + 2 * math.cos(x[3]) * (
                            R2 + r2 * math.cos(x[2])) * r2 * math.sin(
                        x[2]) * n2u2
                    + 2 * math.sin(x[3]) * (
                            R2 + r2 * math.cos(x[2])) * r2 * math.sin(
                        x[2]) * n2v2)

        x01 = npy.array(
            [(min_phi1 + max_phi1) / 2, (min_theta1 + max_theta1) / 2,
             (min_phi2 + max_phi2) / 2, (min_theta2 + max_theta2) / 2])
        x02 = npy.array([min_phi1, min_theta1,
                         min_phi2, min_theta2])
        x03 = npy.array([max_phi1, max_theta1,
                         max_phi2, max_theta2])

        minimax = [(min_phi1, min_theta1, min_phi2, min_theta2),
                   (max_phi1, max_theta1, max_phi2, max_theta2)]

        res1 = scp.optimize.least_squares(distance_squared, x01,
                                          bounds=minimax)
        res2 = scp.optimize.least_squares(distance_squared, x02,
                                          bounds=minimax)
        res3 = scp.optimize.least_squares(distance_squared, x03,
                                          bounds=minimax)

        # frame1, frame2 = volmdlr.Frame3D(self.center, u1, v1, n1), volmdlr.Frame3D(other_tore.center, u2, v2, n2)
        pt1 = self.points2d_to3d([[res1.x[1], res1.x[0]]], R1, r1, frame1)
        pt2 = self.points2d_to3d([[res1.x[3], res1.x[2]]], R2, r2, frame2)
        p1, p2 = pt1[0], pt2[0]
        d = p1.point_distance(p2)
        result = res1

        res = [res2, res3]
        for couple in res:
            ptest1 = self.points2d_to3d([[couple.x[1], couple.x[0]]], R1, r1,
                                        frame1)
            ptest2 = self.points2d_to3d([[couple.x[3], couple.x[2]]], R2, r2,
                                        frame2)
            dtest = ptest1[0].point_distance(ptest2[0])
            if dtest < d:
                result = couple
                p1, p2 = ptest1[0], ptest2[0]

        pt1_2d, pt2_2d = volmdlr.Point2D(
            (result.x[1], result.x[0])), volmdlr.Point2D(
            (result.x[3], result.x[2]))

        if not (self.contours2d[0].point_belongs(pt1_2d)):
            # Find the closest one
            points_contours1 = self.contours2d[0].tessel_points

            poly1 = volmdlr.ClosedPolygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d,
                                                       return_other_point=True)

            pt1 = self.points2d_to3d([new_pt1_2d], R1, r1, frame1)
            p1 = pt1[0]

        if not (other_tore.contours2d[0].point_belongs(pt2_2d)):
            # Find the closest one
            points_contours2 = other_tore.contours2d[0].tessel_points

            poly2 = volmdlr.ClosedPolygon2D(points_contours2)
            d2, new_pt2_2d = poly2.PointBorderDistance(pt2_2d,
                                                       return_other_point=True)

            pt2 = self.points2d_to3d([new_pt2_2d], R2, r2, frame2)
            p2 = pt2[0]

        return p1, p2

    def minimum_distance_points_cyl(self, cyl):
        R2, r2, r = self.rcenter, self.rcircle, cyl.radius

        min_h, min_theta, max_h, max_theta = cyl.minimum_maximum(
            cyl.contours2d[0], r)

        n1 = cyl.normal
        u1 = cyl.cylindricalsurface3d.frame.u
        v1 = cyl.cylindricalsurface3d.frame.v
        frame1 = volmdlr.Frame3D(cyl.center, u1, v1, n1)
        # st1 = volmdlr.Point3D((r*math.cos(min_theta), r*math.sin(min_theta), min_h))
        # start1 = frame1.old_coordinates(st1)

        min_phi2, min_theta2, max_phi2, max_theta2 = self.minimum_maximum_tore(
            self.contours2d[0])

        n2 = self.normal
        u2 = self.toroidalsurface3d.frame.u
        v2 = self.toroidalsurface3d.frame.v
        frame2 = volmdlr.Frame3D(self.center, u2, v2, n2)
        # start2 = self.points2d_to3d([[min_theta2, min_phi2]], R2, r2, frame2)

        w = self.center - cyl.center

        n1n1, n1u1, n1v1, n1n2, n1u2, n1v2 = n1.dot(n1), n1.dot(u1), n1.dot(
            v1), n1.dot(n2), n1.dot(u2), n1.dot(v2)
        u1u1, u1v1, u1n2, u1u2, u1v2 = u1.dot(u1), u1.dot(v1), u1.dot(
            n2), u1.dot(u2), u1.dot(v2)
        v1v1, v1n2, v1u2, v1v2 = v1.dot(v1), v1.dot(n2), v1.dot(u2), v1.dot(v2)
        n2n2, n2u2, n2v2 = n2.dot(n2), n2.dot(u2), n2.dot(v2)
        u2u2, u2v2, v2v2 = u2.dot(u2), u2.dot(v2), v2.dot(v2)

        w2, wn1, wu1, wv1, wn2, wu2, wv2 = w.dot(w), w.dot(n1), w.dot(
            u1), w.dot(v1), w.dot(n2), w.dot(u2), w.dot(v2)

        # x = (theta, h, phi2, theta2)
        def distance_squared(x):
            return (u1u1 * ((math.cos(x[0]) * r) ** 2) + v1v1 * (
                    (math.sin(x[0]) * r) ** 2)
                    + n1n1 * (x[1] ** 2) + w2
                    + u2u2 * (((R2 + r2 * math.cos(x[2])) * math.cos(
                        x[3])) ** 2)
                    + v2v2 * (((R2 + r2 * math.cos(x[2])) * math.sin(
                        x[3])) ** 2)
                    + n2n2 * ((math.sin(x[2])) ** 2) * (r2 ** 2)
                    + 2 * u1v1 * math.cos(x[0]) * math.sin(x[0]) * (r ** 2)
                    + 2 * r * math.cos(x[0]) * x[1] * n1u1 - 2 * r * math.cos(
                        x[0]) * wu1
                    - 2 * r * math.cos(x[0]) * (
                            R2 + r2 * math.cos(x[2])) * math.cos(
                        x[3]) * u1u2
                    - 2 * r * math.cos(x[0]) * (
                            R2 + r2 * math.cos(x[2])) * math.sin(
                        x[3]) * u1v2
                    - 2 * r * math.cos(x[0]) * r2 * math.sin(x[2]) * u1n2
                    + 2 * r * math.sin(x[0]) * x[1] * n1v1 - 2 * r * math.sin(
                        x[0]) * wv1
                    - 2 * r * math.sin(x[0]) * (
                            R2 + r2 * math.cos(x[2])) * math.cos(
                        x[3]) * v1u2
                    - 2 * r * math.sin(x[0]) * (
                            R2 + r2 * math.cos(x[2])) * math.sin(
                        x[3]) * v1v2
                    - 2 * r * math.sin(x[0]) * r2 * math.sin(x[2]) * v1n2 - 2 *
                    x[1] * wn1
                    - 2 * x[1] * (R2 + r2 * math.cos(x[2])) * math.cos(
                        x[3]) * n1u2
                    - 2 * x[1] * (R2 + r2 * math.cos(x[2])) * math.sin(
                        x[3]) * n1v2
                    - 2 * x[1] * r2 * math.sin(x[2]) * n1n2
                    + 2 * (R2 + r2 * math.cos(x[2])) * math.cos(x[3]) * wu2
                    + 2 * (R2 + r2 * math.cos(x[2])) * math.sin(x[3]) * wv2
                    + 2 * r2 * math.sin(x[2]) * wn2
                    + 2 * u2v2 * math.cos(x[3]) * math.sin(x[3]) * (
                            (R2 + r2 * math.cos(x[2])) ** 2)
                    + 2 * math.cos(x[3]) * (
                            R2 + r2 * math.cos(x[2])) * r2 * math.sin(
                        x[2]) * n2u2
                    + 2 * math.sin(x[3]) * (
                            R2 + r2 * math.cos(x[2])) * r2 * math.sin(
                        x[2]) * n2v2)

        x01 = npy.array([(min_theta + max_theta) / 2, (min_h + max_h) / 2,
                         (min_phi2 + max_phi2) / 2,
                         (min_theta2 + max_theta2) / 2])
        x02 = npy.array([min_theta, min_h,
                         min_phi2, min_theta2])
        x03 = npy.array([max_theta, max_h,
                         max_phi2, max_theta2])

        minimax = [(min_theta, min_h, min_phi2, min_theta2),
                   (max_theta, max_h, max_phi2, max_theta2)]

        res1 = scp.optimize.least_squares(distance_squared, x01,
                                          bounds=minimax)
        res2 = scp.optimize.least_squares(distance_squared, x02,
                                          bounds=minimax)
        res3 = scp.optimize.least_squares(distance_squared, x03,
                                          bounds=minimax)

        pt1 = volmdlr.Point3D(
            (r * math.cos(res1.x[0]), r * math.sin(res1.x[0]), res1.x[1]))
        p1 = frame1.old_coordinates(pt1)
        pt2 = self.points2d_to3d([[res1.x[3], res1.x[2]]], R2, r2, frame2)
        p2 = pt2[0]
        d = p1.point_distance(p2)
        result = res1

        res = [res2, res3]
        for couple in res:
            pttest1 = volmdlr.Point3D((r * math.cos(couple.x[0]),
                                       r * math.sin(couple.x[0]), couple.x[1]))
            ptest1 = frame1.old_coordinates(pttest1)
            ptest2 = self.points2d_to3d([[couple.x[3], couple.x[2]]], R2, r2,
                                        frame2)
            dtest = ptest1.point_distance(ptest2[0])
            if dtest < d:
                result = couple
                p1, p2 = ptest1, ptest2[0]

        pt1_2d, pt2_2d = volmdlr.Point2D(
            (result.x[0], result.x[1])), volmdlr.Point2D(
            (result.x[3], result.x[2]))

        if not (self.contours2d[0].point_belongs(pt2_2d)):
            # Find the closest one
            points_contours2 = self.contours2d[0].tessel_points

            poly2 = volmdlr.ClosedPolygon2D(points_contours2)
            d2, new_pt2_2d = poly2.PointBorderDistance(pt2_2d,
                                                       return_other_point=True)

            pt2 = self.points2d_to3d([new_pt2_2d], R2, r2, frame2)
            p2 = pt2[0]

        if not (cyl.contours2d[0].point_belongs(pt1_2d)):
            # Find the closest one
            points_contours1 = cyl.contours2d[0].tessel_points

            poly1 = volmdlr.ClosedPolygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d,
                                                       return_other_point=True)

            pt1 = volmdlr.Point3D((r * math.cos(new_pt1_2d.vector[0]),
                                   r * math.sin(new_pt1_2d.vector[0]),
                                   new_pt1_2d.vector[1]))
            p1 = frame1.old_coordinates(pt1)

        return p1, p2

    def minimum_distance_points_plane(self,
                                      planeface):  # Planeface with contour2D
        # TODO: check that it takes into account holes

        poly2d = planeface.polygon2D
        pfpoints = poly2d.points
        xmin, ymin = min([pt[0] for pt in pfpoints]), min(
            [pt[1] for pt in pfpoints])
        xmax, ymax = max([pt[0] for pt in pfpoints]), max(
            [pt[1] for pt in pfpoints])
        origin, vx, vy = planeface.plane.origin, planeface.plane.vectors[0], \
                         planeface.plane.vectors[1]
        pf1_2d, pf2_2d = volmdlr.Point2D((xmin, ymin)), volmdlr.Point2D(
            (xmin, ymax))
        pf3_2d, pf4_2d = volmdlr.Point2D((xmax, ymin)), volmdlr.Point2D(
            (xmax, ymax))
        pf1, pf2 = pf1_2d.to_3d(origin, vx, vy), pf2_2d.to_3d(origin, vx, vy)
        pf3, _ = pf3_2d.to_3d(origin, vx, vy), pf4_2d.to_3d(origin, vx, vy)

        u, v = (pf3 - pf1), (pf2 - pf1)
        u.normalize()
        v.normalize()

        R1, r1 = self.rcenter, self.rcircle
        min_phi1, min_theta1, max_phi1, max_theta1 = self.minimum_maximum_tore(
            self.contours2d[0])

        n1 = self.normal
        u1 = self.toroidalsurface3d.frame.u
        v1 = self.toroidalsurface3d.frame.v
        frame1 = volmdlr.Frame3D(self.center, u1, v1, n1)
        # start1 = self.points2d_to3d([[min_theta1, min_phi1]], R1, r1, frame1)

        w = self.center - pf1

        n1n1, n1u1, n1v1, n1u, n1v = n1.dot(n1), n1.dot(u1), n1.dot(
            v1), n1.dot(u), n1.dot(v)
        u1u1, u1v1, u1u, u1v = u1.dot(u1), u1.dot(v1), u1.dot(u), u1.dot(v)
        v1v1, v1u, v1v = v1.dot(v1), v1.dot(u), v1.dot(v)
        uu, uv, vv = u.dot(u), u.dot(v), v.dot(v)

        w2, wn1, wu1, wv1, wu, wv = w.dot(w), w.dot(n1), w.dot(u1), w.dot(
            v1), w.dot(u), w.dot(v)

        # x = (x, y, phi1, theta1)
        def distance_squared(x):
            return (uu * (x[0] ** 2) + vv * (x[1] ** 2) + w2
                    + u1u1 * (((R1 + r1 * math.cos(x[2])) * math.cos(
                        x[3])) ** 2)
                    + v1v1 * (((R1 + r1 * math.cos(x[2])) * math.sin(
                        x[3])) ** 2)
                    + n1n1 * ((math.sin(x[2])) ** 2) * (r1 ** 2)
                    + 2 * x[0] * x[1] * uv - 2 * x[0] * wu
                    - 2 * x[0] * (R1 + r1 * math.cos(x[2])) * math.cos(
                        x[3]) * u1u
                    - 2 * x[0] * (R1 + r1 * math.cos(x[2])) * math.sin(
                        x[3]) * v1u
                    - 2 * x[0] * math.sin(x[2]) * r1 * n1u - 2 * x[1] * wv
                    - 2 * x[1] * (R1 + r1 * math.cos(x[2])) * math.cos(
                        x[3]) * u1v
                    - 2 * x[1] * (R1 + r1 * math.cos(x[2])) * math.sin(
                        x[3]) * v1v
                    - 2 * x[1] * math.sin(x[2]) * r1 * n1v
                    + 2 * (R1 + r1 * math.cos(x[2])) * math.cos(x[3]) * wu1
                    + 2 * (R1 + r1 * math.cos(x[2])) * math.sin(x[3]) * wv1
                    + 2 * math.sin(x[2]) * r1 * wn1
                    + 2 * u1v1 * math.cos(x[3]) * math.sin(x[3]) * (
                            (R1 + r1 * math.cos(x[2])) ** 2)
                    + 2 * (R1 + r1 * math.cos(x[2])) * math.cos(
                        x[3]) * r1 * math.sin(x[2]) * n1u1
                    + 2 * (R1 + r1 * math.cos(x[2])) * math.sin(
                        x[3]) * r1 * math.sin(x[2]) * n1v1)

        x01 = npy.array([(xmax - xmin) / 2, (ymax - ymin) / 2,
                         (min_phi1 + max_phi1) / 2,
                         (min_theta1 + max_theta1) / 2])

        minimax = [(0, 0, min_phi1, min_theta1),
                   (xmax - xmin, ymax - ymin, max_phi1, max_theta1)]

        res1 = scp.optimize.least_squares(distance_squared, x01,
                                          bounds=minimax)

        # frame1 = volmdlr.Frame3D(self.center, u1, v1, n1)
        pt1 = self.points2d_to3d([[res1.x[3], res1.x[2]]], R1, r1, frame1)
        p1 = pt1[0]
        p2 = pf1 + res1.x[2] * u + res1.x[3] * v

        pt1_2d = volmdlr.Point2D((res1.x[3], res1.x[2]))
        pt2_2d = p2.to_2d(pf1, u, v)

        if not (self.contours2d[0].point_belongs(pt1_2d)):
            # Find the closest one
            points_contours1 = self.contours2d[0].tessel_points

            poly1 = volmdlr.ClosedPolygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d,
                                                       return_other_point=True)

            pt1 = self.points2d_to3d([new_pt1_2d], R1, r1, frame1)
            p1 = pt1[0]

        if not (planeface.contours[0].point_belongs(pt2_2d)):
            # Find the closest one
            d2, new_pt2_2d = planeface.polygon2D.PointBorderDistance(pt2_2d,
                                                                     return_other_point=True)

            p2 = new_pt2_2d.to_3d(pf1, u, v)

        return p1, p2

    def minimum_distance(self, other_face, return_points=False):
        if other_face.__class__ is ToroidalFace3D:
            p1, p2 = self.minimum_distance_points_tore(other_face)
            if return_points:
                return p1.point_distance(p2), p1, p2
            else:
                return p1.point_distance(p2)

        if other_face.__class__ is CylindricalFace3D:
            p1, p2 = self.minimum_distance_points_cyl(other_face)
            if return_points:
                return p1.point_distance(p2), p1, p2
            else:
                return p1.point_distance(p2)

        if other_face.__class__ is PlaneFace3D:
            p1, p2 = self.minimum_distance_points_plane(other_face)
            if return_points:
                return p1.point_distance(p2), p1, p2
            else:
                return p1.point_distance(p2)
        else:
            return NotImplementedError


class ConicalFace3D(Face3D):
    """
    :param contours2d: The Cone's contour2D
    :type contours2d: volmdlr.Contour2D
    :param conicalsurface3d: Information about the Cone
    :type conicalsurface3d: ConicalSurface3D
    :param points: Contour2d's parameter Cone
    :type points: List of float

    """
    min_x_density = 5
    min_y_density = 1

    def __init__(self, conicalsurface3d: ConicalSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):

        Face3D.__init__(self,
                        surface3d=conicalsurface3d,
                        surface2d=surface2d,
                        name=name)

    def _bounding_box(self):
        theta_min, theta_max, zmin, zmax = self.surface2d.outer_contour.bounding_rectangle()

        xp = (volmdlr.X3D.dot(self.surface3d.frame.u) * self.surface3d.frame.u
              + volmdlr.X3D.dot(
                    self.surface3d.frame.v) * self.surface3d.frame.v)
        try:
            xp.normalize()
        except ZeroDivisionError:
            pass
        yp = (volmdlr.Y3D.dot(self.surface3d.frame.u) * self.surface3d.frame.u
              + volmdlr.Y3D.dot(
                    self.surface3d.frame.v) * self.surface3d.frame.v)

        try:
            yp.normalize()
        except ZeroDivisionError:
            pass

        zp = (volmdlr.Z3D.dot(self.surface3d.frame.u) * self.surface3d.frame.u
              + volmdlr.Z3D.dot(
                    self.surface3d.frame.v) * self.surface3d.frame.v)
        try:
            zp.normalize()
        except ZeroDivisionError:
            pass

        lower_center = self.surface3d.frame.origin + zmin * self.surface3d.frame.w
        upper_center = self.surface3d.frame.origin + zmax * self.surface3d.frame.w
        lower_radius = math.tan(self.surface3d.semi_angle) * zmin
        upper_radius = math.tan(self.surface3d.semi_angle) * zmax

        points = [lower_center - lower_radius * xp,
                  lower_center + lower_radius * xp,
                  lower_center - lower_radius * yp,
                  lower_center + lower_radius * yp,
                  lower_center - lower_radius * zp,
                  lower_center + lower_radius * zp,
                  upper_center - upper_radius * xp,
                  upper_center + upper_radius * xp,
                  upper_center - upper_radius * yp,
                  upper_center + upper_radius * yp,
                  upper_center - upper_radius * zp,
                  upper_center + upper_radius * zp,
                  ]

        return volmdlr.core.BoundingBox.from_points(points)

    def triangulation_lines(self, angle_resolution=5):
        theta_min, theta_max, zmin, zmax = self.surface2d.bounding_rectangle()
        delta_theta = theta_max - theta_min
        nlines = int(delta_theta * angle_resolution)
        lines_x = []
        for i in range(nlines):
            theta = theta_min + (i + 1) / (nlines + 1) * delta_theta
            lines_x.append(vme.Line2D(volmdlr.Point2D(theta, zmin),
                                      volmdlr.Point2D(theta, zmax)))

        if zmin < 1e-9:
            delta_z = zmax - zmin
            lines_y = [vme.Line2D(volmdlr.Point2D(theta_min, zmin + 0.1 * delta_z),
                                  volmdlr.Point2D(theta_max, zmin + 0.1 * delta_z))]
        else:
            lines_y = []
        return lines_x, lines_y

    # def create_triangle(self, all_contours_points, part):
    #     Triangles, ts = [], []
    #     pts, h_list = [], []
    #     for listpt in all_contours_points:
    #         for pt in listpt:
    #             pts.append(pt)
    #             h_list.append(pt[1])
    #     if part == 'bot':
    #         h_concerned = min(h_list)
    #     else:
    #         h_concerned = max(h_list)
    #     peak_list, other = [], []
    #     for pt in pts:
    #         if pt[1] == h_concerned:
    #             peak_list.append(pt)
    #         else:
    #             other.append(pt)
    #     points = [peak_list[0]] + other
    #
    #     for i in range(1, len(points)):
    #         if i == len(points) - 1:
    #             vertices = [points[i].vector, points[0].vector,
    #                         points[1].vector]
    #             segments = [[0, 1], [1, 2], [2, 0]]
    #             listindice = [i, 0, 1]
    #         else:
    #             vertices = [points[i].vector, points[0].vector,
    #                         points[i + 1].vector]
    #             segments = [[0, 1], [1, 2], [2, 0]]
    #             listindice = [i, 0, i + 1]
    #         tri = {'vertices': vertices, 'segments': segments}
    #         t = triangle.triangulate(tri, 'p')
    #         if 'triangles' in t:
    #             triangles = t['triangles'].tolist()
    #             triangles[0] = listindice
    #             Triangles.append(triangles)
    #         else:
    #             Triangles.append(None)
    #         ts.append(t)
    #
    #     return points, Triangles


class SphericalFace3D(Face3D):
    """
    :param contours2d: The Sphere's contour2D
    :type contours2d: volmdlr.Contour2D
    :param sphericalsurface3d: Information about the Sphere
    :type sphericalsurface3d: SphericalSurface3D
    :param points: Angle's Sphere
    :type points: List of float

    """
    min_x_density = 5
    min_y_density = 5

    def __init__(self, spherical_surface3d: SphericalSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):
        Face3D.__init__(self,
                        surface3d=spherical_surface3d,
                        surface2d=surface2d,
                        name=name)

    def _bounding_box(self):
        # To be enhanced
        return self.surface3d._bounding_box()

    def triangulation_lines(self, angle_resolution=7):
        theta_min, theta_max, phi_min, phi_max = self.surface2d.bounding_rectangle()

        delta_theta = theta_max - theta_min
        nlines_x = int(delta_theta * angle_resolution)
        lines_x = []
        for i in range(nlines_x):
            theta = theta_min + (i + 1) / (nlines_x + 1) * delta_theta
            lines_x.append(vme.Line2D(volmdlr.Point2D(theta, phi_min),
                                      volmdlr.Point2D(theta, phi_max)))
        delta_phi = phi_max - phi_min
        nlines_y = int(delta_phi * angle_resolution)
        lines_y = []
        for i in range(nlines_y):
            phi = phi_min + (i + 1) / (nlines_y + 1) * delta_phi
            lines_y.append(vme.Line2D(volmdlr.Point2D(theta_min, phi),
                                      volmdlr.Point2D(theta_max, phi)))
        return lines_x, lines_y


class RuledFace3D(Face3D):
    """

    """
    min_x_density = 50
    min_y_density = 1

    def __init__(self,
                 ruledsurface3d: RuledSurface3D,
                 surface2d: Surface2D,
                 name: str = '',
                 color=None):
        Face3D.__init__(self, surface3d=ruledsurface3d,
                        surface2d=surface2d,
                        name=name)

    def _bounding_box(self):
        # To be enhance by restricting wires to cut
        # xmin, xmax, ymin, ymax = self.surface2d.outer_contour.bounding_rectangle()
        points = [self.surface3d.point2d_to_3d(volmdlr.Point2D(i / 30, 0.)) for
                  i in range(31)]
        points.extend(
            [self.surface3d.point2d_to_3d(volmdlr.Point2D(i / 30, 1.)) for i
             in range(31)])

        return volmdlr.core.BoundingBox.from_points(points)

    def triangulation_lines(self, angle_resolution=10):
        xmin, xmax, ymin, ymax = self.surface2d.bounding_rectangle()
        delta_x = xmax - xmin
        nlines = int(delta_x * angle_resolution)
        lines = []
        for i in range(nlines):
            x = xmin + (i + 1) / (nlines + 1) * delta_x
            lines.append(vme.Line2D(volmdlr.Point2D(x, ymin),
                                    volmdlr.Point2D(x, ymax)))
        return lines, []


class BSplineFace3D(Face3D):
    def __init__(self, bspline_surface: BSplineSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):
        Face3D.__init__(self,
                        surface3d=bspline_surface,
                        surface2d=surface2d,
                        name=name)

    def _bounding_box(self):
        return self.surface3d._bounding_box()

    def triangulation_lines(self, resolution=10):
        u_min, u_max, v_min, v_max = self.surface2d.bounding_rectangle()

        delta_u = u_max - u_min
        nlines_x = int(delta_u * resolution)
        lines_x = []
        for i in range(nlines_x):
            u = u_min + (i + 1) / (nlines_x + 1) * delta_u
            lines_x.append(vme.Line2D(volmdlr.Point2D(u, v_min),
                                      volmdlr.Point2D(u, v_max)))
        delta_v = v_max - v_min
        nlines_y = int(delta_v * resolution)
        lines_y = []
        for i in range(nlines_y):
            v = v_min + (i + 1) / (nlines_y + 1) * delta_v
            lines_y.append(vme.Line2D(volmdlr.Point2D(v_min, v),
                                      volmdlr.Point2D(v_max, v)))
        return lines_x, lines_y


class OpenShell3D(volmdlr.core.CompositePrimitive3D):
    _standalone_in_db = True
    _non_serializable_attributes = ['bounding_box']
    _non_eq_attributes = ['name', 'color', 'alpha' 'bounding_box']
    _non_hash_attributes = []
    STEP_FUNCTION = 'OPEN_SHELL'

    def __init__(self, faces: List[Face3D],
                 color: Tuple[float, float, float] = None,
                 alpha: float = 1., name: str = ''):
        self.faces = faces
        self.name = name
        if not color:
            self.color = (0.8, 0.8, 0.8)
        else:
            self.color = color
        self.alpha = alpha
        self.bounding_box = self._bounding_box()

    def __hash__(self):
        return sum([hash(f) for f in self.faces])

    def __eq__(self, other_):
        if self.__class__ != other_.__class__:
            return False
        equal = True
        for face, other_face in zip(self.faces, other_.faces):
            equal = (equal and face == other_face)
        return equal

    @classmethod
    def from_step(cls, arguments, object_dict):
        faces = []
        for face in arguments[1]:
            faces.append(object_dict[int(face[1:])])
        return cls(faces, name=arguments[0][1:-1])

    def to_step(self, current_id):
        step_content = ''
        face_ids = []
        for face in self.faces:
            face_content, face_sub_ids = face.to_step(current_id)
            step_content += face_content
            face_ids.extend(face_sub_ids)
            current_id = max(face_sub_ids) + 1

        shell_id = current_id
        step_content += "#{} = {}('{}',({}));\n".format(current_id,
                                                        self.STEP_FUNCTION,
                                                        self.name,
                                                        volmdlr.core.step_ids_to_str(
                                                            face_ids))
        manifold_id = shell_id + 1
        step_content += "#{} = MANIFOLD_SOLID_BREP('{}',#{});\n".format(
            manifold_id,
            self.name,
            shell_id)

        frame_content, frame_id = volmdlr.OXYZ.to_step(manifold_id + 1)
        step_content += frame_content
        brep_id = frame_id + 1
        step_content += "#{} = ADVANCED_BREP_SHAPE_REPRESENTATION('',(#{},#{}),#7);\n".format(
            brep_id, frame_id, manifold_id)

        return step_content, brep_id

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            new_faces = [face.rotation(center, axis, angle, copy=True) for face
                         in self.faces]
            return OpenShell3D(new_faces, color=self.color, alpha=self.alpha, name=self.name)
        else:
            for face in self.faces:
                face.rotation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def translation(self, offset, copy=True):
        if copy:
            new_faces = [face.translation(offset, copy=True) for face in
                         self.faces]
            return OpenShell3D(new_faces, color=self.color, alpha=self.alpha, name=self.name)
        else:
            for face in self.faces:
                face.translation(offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_faces = [face.frame_mapping(frame, side, copy=True) for face in
                         self.faces]
            return self.__class__(new_faces, name=self.name)
        else:
            for face in self.faces:
                face.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()

    def copy(self):
        new_faces = [face.copy() for face in self.faces]
        return self.__class__(new_faces, color=self.color, alpha=self.alpha,
                              name=self.name)

    def union(self, shell2):
        new_faces = [face for face in self.faces + shell2.faces]
        new_name = self.name + ' union ' + shell2.name
        new_color = self.color
        return self.__class__(new_faces, name=new_name, color=new_color)

    def volume(self):
        """
        Does not consider holes
        """
        volume = 0
        for i, face in enumerate(self.faces):
            points_3D, triangles_indexes = face.triangulation()
            for triangle_indexes in triangles_indexes[0]:
                point1 = points_3D[triangle_indexes[0]]
                point2 = points_3D[triangle_indexes[1]]
                point3 = points_3D[triangle_indexes[2]]

                v321 = point3[0] * point2[1] * point1[2]
                v231 = point2[0] * point3[1] * point1[2]
                v312 = point3[0] * point1[1] * point2[2]
                v132 = point1[0] * point3[1] * point2[2]
                v213 = point2[0] * point1[1] * point3[2]
                v123 = point1[0] * point2[1] * point3[2]
                volume_tetraedre = 1 / 6 * (
                        -v321 + v231 + v312 - v132 - v213 + v123)

                volume += volume_tetraedre

        return abs(volume)

    def _bounding_box(self):
        """
        Returns the boundary box
        """
        bbox = self.faces[0]._bounding_box()

        for face in self.faces[1:]:
            bbox += face._bounding_box()

        return bbox

    def cut_by_plane(self, plane_3d: Plane3D):
        graph = nx.Graph()
        intersections = []

        frame_block = self.bounding_box.to_frame()
        frame_block.u = 1.1 * frame_block.u
        frame_block.v = 1.1 * frame_block.v
        frame_block.w = 1.1 * frame_block.w

        for face in self.faces:
            block = volmdlr.primitives3d.Block(frame_block)
            face_3d = block.cut_by_orthogonal_plane(plane_3d)
            inters = face.face_intersections(face_3d)
            if inters:
                graph.add_edges_from([(inters[0], inters[1])])
                intersections.append(inters)
        pts = list(nx.dfs_edges(graph, intersections[0][0]))
        points = []
        u = plane_3d.frame.u
        v = plane_3d.frame.v
        for pt1, pt2 in pts:
            if pt1 not in points:
                points.append(pt1)
            if pt2 not in points:
                points.append(pt2)
        center_2d = volmdlr.Point2D(plane_3d.frame.origin.dot(u), plane_3d.frame.origin.dot(v))
        points_2d = [volmdlr.Point2D(p.dot(u), p.dot(v)) - center_2d for p in points]
        contour_2d = volmdlr.faces.Surface2D(volmdlr.wires.ClosedPolygon2D(points_2d), [])

        return volmdlr.faces.PlaneFace3D(plane_3d, contour_2d)


    def linesegment_intersections(self,
                                  linesegment3d: vme.LineSegment3D) \
            -> List[Tuple[Face3D, List[volmdlr.Point3D]]]:
        intersections = []
        for face in self.faces:
            face_intersections = face.linesegment_intersections(linesegment3d)
            if face_intersections:
                intersections.append((face, face_intersections))
        return intersections

    def minimum_distance_points(self, shell2, resolution):
        """
        Returns a Mesure object if the distance is not zero, otherwise returns None
        """
        shell2_inter = self.shell_intersection(shell2, resolution)
        if shell2_inter is not None and shell2_inter != 1:
            return None

        # distance_min, point1_min, point2_min = self.faces[0].distance_to_face(shell2.faces[0], return_points=True)
        distance_min, point1_min, point2_min = self.faces[0].minimum_distance(
            shell2.faces[0], return_points=True)
        for face1 in self.faces:
            bbox1 = face1.bounding_box
            for face2 in shell2.faces:
                bbox2 = face2.bounding_box
                bbox_distance = bbox1.distance_to_bbox(bbox2)

                if bbox_distance < distance_min:
                    # distance, point1, point2 = face1.distance_to_face(face2, return_points=True)
                    distance, point1, point2 = face1.minimum_distance(face2,
                                                                      return_points=True)
                    if distance == 0:
                        return None
                    elif distance < distance_min:
                        distance_min, point1_min, point2_min = distance, point1, point2

        return point1_min, point2_min

    def distance_to_shell(self, other_shell: 'OpenShell3D', resolution: float):
        min_dist = self.minimum_distance_points(other_shell, resolution)
        if min_dist is not None:
            p1, p2 = min_dist
            return p1.point_distance(p2)
        else:
            return None

    def minimum_distance_point(self,
                               point: volmdlr.Point3D) -> volmdlr.Point3D:
        """
        Computes the distance of a point to a Shell3D, whether it is inside or outside the Shell3D
        """
        distance_min, point1_min = self.faces[0].distance_to_point(point,
                                                                   return_other_point=True)
        for face in self.faces[1:]:
            bbox_distance = self.bounding_box.distance_to_point(point)
            if bbox_distance < distance_min:
                distance, point1 = face.distance_to_point(point,
                                                          return_other_point=True)
                if distance < distance_min:
                    distance_min, point1_min = distance, point1

        return point1_min

    def intersection_internal_aabb_volume(self, shell2: 'OpenShell3D',
                                          resolution: float):
        """
        aabb made of the intersection points and the points of self internal to shell2
        """
        intersections_points = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersections(face2)
                if intersection_points is not None:
                    intersections_points.extend(intersection_points)

        shell1_points_inside_shell2 = []
        for face in self.faces:
            for point in face.outer_contour3d.discretization_points(
                    resolution):
                if shell2.point_belongs(point):
                    shell1_points_inside_shell2.append(point)

        if len(intersections_points + shell1_points_inside_shell2) == 0:
            return 0
        bbox = volmdlr.core.BoundingBox.from_points(
            intersections_points + shell1_points_inside_shell2)
        return bbox.volume()

    def intersection_external_aabb_volume(self, shell2: 'OpenShell3D',
                                          resolution: float):
        """
        aabb made of the intersection points and the points of self external to shell2
        """
        intersections_points = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersections(face2)
                if intersection_points is not None:
                    intersections_points.extend(intersection_points)

        shell1_points_outside_shell2 = []
        for face in self.faces:
            for point in face.outer_contour3d.discretization_points(
                    resolution):
                if not shell2.point_belongs(point):
                    shell1_points_outside_shell2.append(point)

        if len(intersections_points + shell1_points_outside_shell2) == 0:
            return 0
        bbox = volmdlr.core.BoundingBox.from_points(
            intersections_points + shell1_points_outside_shell2)
        return bbox.volume()

    def primitive_inside_bbox(self, bounding_box: volmdlr.core.BoundingBox):
        for primitive in self.primitives:
            bbox = primitive.bounding_box

    def triangulation(self):
        mesh = volmdlr.display.DisplayMesh3D([], [])
        for i, face in enumerate(self.faces):
            try:
                face_mesh = face.triangulation()
                mesh += face_mesh
            except NotImplementedError:
                print('Warning: a face has been skipped in rendering')
        return mesh

    def babylon_script(self, name='primitive_mesh'):
        s = 'var {} = new BABYLON.Mesh("{}", scene);\n'.format(name, name)

        mesh = self.babylon_meshes()[0]

        s += 'var positions = {};\n'.format(mesh['positions'])
        s += 'var indices = {};\n'.format(mesh['indices'])
        s += 'var normals = [];\n'
        s += 'var vertexData = new BABYLON.VertexData();\n'
        s += 'BABYLON.VertexData.ComputeNormals(positions, indices, normals);\n'
        s += 'vertexData.positions = positions;\n'
        s += 'vertexData.indices = indices;\n'
        s += 'vertexData.normals = normals;\n'
        s += 'vertexData.applyToMesh({});\n'.format(name)
        s += '{}.enableEdgesRendering(0.9);\n'.format(name)
        s += '{}.edgesWidth = 0.1;\n'.format(name)
        s += '{}.edgesColor = new BABYLON.Color4(0, 0, 0, 0.6);\n'.format(name)
        s += 'var mat = new BABYLON.StandardMaterial("mat", scene);\n'
        #        s += 'mat.diffuseColor = BABYLON.Color3.Green();\n'
        #        s += 'mat.specularColor = new BABYLON.Color3(0.5, 0.6, 0.87);\n'
        #        s += 'mat.emissiveColor = new BABYLON.Color3(1, 1, 1);\n'
        #        s += 'mat.ambientColor = new BABYLON.Color3(0.23, 0.98, 0.53);\n'
        s += 'mat.backFaceCulling = false;\n'
        s += 'mat.alpha = {};\n'.format(self.alpha)
        s += '{}.material = mat;\n'.format(name)
        if self.color is not None:
            s += 'mat.diffuseColor = new BABYLON.Color3({}, {}, {});\n'.format(
                *self.color)
        return s

    def plot(self, ax=None, equal_aspect=True, color='k', alpha=1):
        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')

        for face in self.faces:
            face.plot(ax=ax, color=color, alpha=alpha)

        return ax


class ClosedShell3D(OpenShell3D):
    _standalone_in_db = True
    _non_serializable_attributes = ['bounding_box']
    _non_eq_attributes = ['name', 'color', 'alpha' 'bounding_box']
    STEP_FUNCTION = 'CLOSED_SHELL'

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            new_faces = [face.rotation(center, axis, angle, copy=True) for face
                         in self.faces]
            return ClosedShell3D(new_faces, color=self.color, alpha=self.alpha, name=self.name)
        else:
            for face in self.faces:
                face.rotation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def translation(self, offset, copy=True):
        if copy:
            new_faces = [face.translation(offset, copy=True) for face in
                         self.faces]
            return ClosedShell3D(new_faces, color=self.color, alpha=self.alpha, name=self.name)
        else:
            for face in self.faces:
                face.translation(offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_faces = [face.frame_mapping(frame, side, copy=True) for face in
                         self.faces]
            return ClosedShell3D(new_faces, name=self.name)
        else:
            for face in self.faces:
                face.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()

    def copy(self):
        new_faces = [face.copy() for face in self.faces]
        return ClosedShell3D(new_faces, color=self.color, alpha=self.alpha,
                             name=self.name)

    def shell_intersection(self, shell2: 'OpenShell3D', resolution: float):
        """
        Return None if disjointed
        Return (1, 0) or (0, 1) if one is inside the other
        Return (n1, n2) if intersection

        4 cases :
            (n1, n2) with face intersection             => (n1, n2)
            (0, 0) with face intersection               => (0, 0)
            (0, 0) with no face intersection            => None
            (1, 0) or (0, 1) with no face intersection  => 1
        """
        # Check if boundary boxes don't intersect
        bbox1 = self.bounding_box
        bbox2 = shell2.bounding_box
        if not bbox1.bbox_intersection(bbox2):
            # print("No intersection of shells' BBox")
            return None

        # Check if any point of the first shell is in the second shell
        points1 = []
        for face in self.faces:
            points1.extend(
                face.outer_contour3d.discretization_points(resolution))
        points2 = []
        for face in shell2.faces:
            points2.extend(
                face.outer_contour3d.discretization_points(resolution))

        nb_pts1 = len(points1)
        nb_pts2 = len(points2)
        compteur1 = 0
        compteur2 = 0
        for point1 in points1:
            if shell2.point_belongs(point1):
                compteur1 += 1
        for point2 in points2:
            if self.point_belongs(point2):
                compteur2 += 1

        inter1 = compteur1 / nb_pts1
        inter2 = compteur2 / nb_pts2

        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersections(face2)
                if intersection_points:
                    return inter1, inter2

        if inter1 == 0. and inter2 == 0.:
            return None
        return 1

    def point_belongs(self, point3d: volmdlr.Point3D, nb_rays: int = 1):
        """
        Ray Casting algorithm
        Returns True if the point is inside the Shell, False otherwise
        """

        bbox = self.bounding_box
        if not bbox.point_belongs(point3d):
            return False

        min_ray_length = 2 * max((bbox.xmax - bbox.xmin,
                                  bbox.ymax - bbox.ymin,
                                  bbox.zmax - bbox.zmin))
        two_min_ray_length = 2 * min_ray_length

        rays = []
        for k in range(0, nb_rays):
            rays.append(vme.LineSegment3D(
                point3d,
                point3d + volmdlr.Point3D.random(min_ray_length,
                                                 two_min_ray_length,
                                                 min_ray_length,
                                                 two_min_ray_length,
                                                 min_ray_length,
                                                 two_min_ray_length)))
        rays = sorted(rays, key=lambda ray: ray.length())

        rays_intersections = []
        tests = []

        # for ray in rays[:3]:
        for ray in rays[:nb_rays]:
            #
            count = 0
            ray_intersection = []
            is_inside = True
            for face, point_inters in self.linesegment_intersections(ray):
                count += len(point_inters)

            if count % 2 == 0:
                is_inside = False
            tests.append(is_inside)
            rays_intersections.append(ray_intersection)

        for test1, test2 in zip(tests[:-1], tests[1:]):
            if test1 != test2:
                raise ValueError
        return tests[0]

    def is_inside_shell(self, shell2, resolution: float):
        """
        Returns True if all the points of self are inside shell2 and no face \
        are intersecting
        This method is not exact
        """
        bbox1 = self.bounding_box
        bbox2 = shell2.bounding_box
        if not bbox1.is_inside_bbox(bbox2):
            return False

        points = []
        for face in self.faces:
            points.extend(
                face.outer_contour3d.discretization_points(resolution))
        for point in points:
            if not shell2.point_belongs(point):
                return False

        # Check if any faces are intersecting
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersections(face2)
                if intersection_points != []:
                    return False

        return True
