"""
Surfaces & faces
"""


from typing import List, Tuple
import math

from itertools import product

import triangle
import numpy as npy

import scipy as scp
import scipy.optimize as opt

import matplotlib.pyplot as plt
# import matplotlib.tri as plt_tri
# from pygeodesic import geodesic

import networkx as nx

from geomdl import BSpline
from geomdl import utilities
from geomdl.fitting import interpolate_surface, approximate_surface
from geomdl.operations import split_surface_u, split_surface_v

import dessia_common as dc
import volmdlr.core
import volmdlr.core_compiled
import volmdlr.edges as vme
import volmdlr.wires
import volmdlr.display as vmd
import volmdlr.geometry


def knots_vector_inv(knots_vector):
    '''
    compute knot elements and multiplicities based on the global knot vector
    '''

    knots = sorted(set(knots_vector))
    multiplicities = []
    for knot in knots:
        multiplicities.append(knots_vector.count(knot))

    return (knots, multiplicities)


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

    def second_moment_area(self, point: volmdlr.Point2D):
        Ix, Iy, Ixy = self.outer_contour.second_moment_area(point)
        for contour in self.inner_contours:
            Ixc, Iyc, Ixyc = contour.second_moment_area(point)
            Ix -= Ixc
            Iy -= Iyc
            Ixy -= Ixyc
        return Ix, Iy, Ixy

    def center_of_mass(self):
        center = self.outer_contour.area() * self.outer_contour.center_of_mass()
        for contour in self.inner_contours:
            center -= contour.area() * contour.center_of_mass()
        return center / self.area()

    def point_belongs(self, point2d: volmdlr.Point2D):
        if not self.outer_contour.point_belongs(point2d):
            return False

        for inner_contour in self.inner_contours:
            if inner_contour.point_belongs(point2d):
                return False

        return True

    def random_point_inside(self):
        '''
             returns a random point inside surface2d. Considers if it has holes
        '''
        valid_point = False
        point_inside_outer_contour = None
        while not valid_point:
            point_inside_outer_contour = self.outer_contour.random_point_inside()
            inside_inner_contour = False
            for inner_contour in self.inner_contours:
                if inner_contour.point_belongs(point_inside_outer_contour):
                    inside_inner_contour = True
            if not inside_inner_contour and\
                    point_inside_outer_contour is not None:
                valid_point = True

        return point_inside_outer_contour

    def triangulation(self, min_x_density=None, min_y_density=None):

        if self.area() == 0.:
            return vmd.DisplayMesh2D([], triangles=[])

        outer_polygon = self.outer_contour.to_polygon(angle_resolution=10)

        if not self.inner_contours:  # No holes
            return outer_polygon.triangulation()
        points = [vmd.Node2D(*p) for p in outer_polygon.points]
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
        points = [vmd.Node2D(*t['vertices'][i, :]) for i in
                  range(np)]

        return vmd.DisplayMesh2D(points, triangles=triangles, edges=None)

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
        bounding_rectangle = self.outer_contour.bounding_rectangle()
        lines = []
        for i in range(n - 1):
            xi = bounding_rectangle[0] + (i + 1) * (bounding_rectangle[1] - bounding_rectangle[0]) / n
            lines.append(vme.Line2D(volmdlr.Point2D(xi, 0),
                                    volmdlr.Point2D(xi, 1)))
        return self.split_by_lines(lines)

    def cut_by_line(self, line: vme.Line2D):
        """
        This method makes inner contour disappear for now
        """
        splitted_outer_contours = self.outer_contour.cut_by_line(line)

        return [Surface2D(oc, []) for oc in splitted_outer_contours]

    def split_at_centers(self):
        """
        Split in n slices
        """
        # xmin, xmax, ymin, ymax = self.outer_contour.bounding_rectangle()

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

                # sp11, sp12 = intersections[2][1].split(intersections[2][0])
                # sp21, sp22 = intersections[3][1].split(intersections[3][0])
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
                # raise NotImplementedError(
                #     '{} intersections not supported yet'.format(
                #         len(intersections)))

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

        lc3d = len(contours3d)

        if lc3d == 1:
            outer_contour2d = self.contour3d_to_2d(contours3d[0])
            inner_contours2d = []
        elif lc3d > 1:
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
        else:
            raise ValueError('Must have at least one contour')

        if isinstance(self.face_class, str):
            class_ = globals()[self.face_class]
        else:
            class_ = self.face_class

        surface2d = Surface2D(outer_contour=outer_contour2d,
                              inner_contours=inner_contours2d)
        return class_(self,
                      surface2d=surface2d,
                      name=name)

    def repair_primitives_periodicity(self, primitives, last_primitive):
        delta_x1 = abs(primitives[0].start.x
                       - last_primitive.end.x)
        delta_x2 = abs(primitives[-1].end.x
                       - last_primitive.end.x)
        delta_y1 = abs(primitives[0].start.y
                       - last_primitive.end.y)
        delta_y2 = abs(primitives[-1].end.y
                       - last_primitive.end.y)

        if self.x_periodicity \
                and not (math.isclose(delta_x1, 0,
                                      abs_tol=5e-5)
                         or math.isclose(delta_x2, 0,
                                         abs_tol=5e-5)):
            delta_x1 = delta_x1 % self.x_periodicity
            delta_x2 = delta_x2 % self.x_periodicity
            if math.isclose(delta_x1, self.x_periodicity,
                            abs_tol=1e-4):
                delta_x1 = 0.
            if math.isclose(delta_x2, self.x_periodicity,
                            abs_tol=1e-4):
                delta_x2 = 0.
            for prim in primitives:
                prim.start.x = abs(self.x_periodicity
                                   - prim.start.x)
                prim.end.x = abs(self.x_periodicity
                                 - prim.end.x)

        if self.y_periodicity \
                and not (math.isclose(delta_y1, 0,
                                      abs_tol=5e-5)
                         or math.isclose(delta_y2, 0,
                                         abs_tol=5e-5)):
            delta_y1 = delta_y1 % self.y_periodicity
            delta_y2 = delta_y2 % self.y_periodicity
            if math.isclose(delta_y1, self.y_periodicity,
                            abs_tol=1e-4):
                delta_y1 = 0.
            if math.isclose(delta_y2, self.y_periodicity,
                            abs_tol=1e-4):
                delta_y2 = 0.
            for prim in primitives:
                prim.start.y = abs(self.y_periodicity
                                   - prim.start.y)
                prim.end.y = abs(self.y_periodicity
                                 - prim.end.y)

        return primitives, delta_x1, delta_x2, delta_y1, delta_y2

    def contour3d_to_2d(self, contour3d):
        primitives2d = []
        last_primitive = None

        for primitive3d in contour3d.primitives:
            method_name = '{}_to_2d'.format(
                primitive3d.__class__.__name__.lower())
            if hasattr(self, method_name):
                primitives = getattr(self, method_name)(primitive3d)

                if primitives is None:
                    continue

                if last_primitive:
                    primitives, delta_x1, delta_x2, delta_y1, delta_y2 = \
                        self.repair_primitives_periodicity(primitives,
                                                           last_primitive)

                    dist1 = primitive3d.start.point_distance(
                        last_primitive3d.end)
                    dist2 = primitive3d.end.point_distance(
                        last_primitive3d.end)
                    if (math.isclose(delta_x1, 0., abs_tol=1e-3)
                            and math.isclose(delta_y1, 0., abs_tol=1e-3)
                            and math.isclose(dist1, 0, abs_tol=5e-5)):
                        pass
                    elif (math.isclose(delta_x2, 0., abs_tol=1e-3)
                            and math.isclose(delta_y2, 0., abs_tol=1e-3)
                            and math.isclose(dist2, 0, abs_tol=5e-5)):
                        primitives = [p.reverse() for p in primitives[::-1]]
                    else:
                        ax2 = contour3d.plot()
                        primitive3d.plot(ax=ax2, color='r')
                        last_primitive3d.plot(ax=ax2, color='b')
                        self.plot(ax=ax2)

                        ax = last_primitive.plot(color='b', plot_points=True)
                        # primitives[0].plot(ax=ax, color='r', plot_points=True)
                        # primitives[-1].plot(ax=ax, color='r', plot_points=True)
                        for p in primitives:
                            p.plot(ax=ax, color='r', plot_points=True)
                        if self.x_periodicity:
                            vme.Line2D(volmdlr.Point2D(self.x_periodicity, 0),
                                       volmdlr.Point2D(self.x_periodicity, 1))\
                                .plot(ax=ax)
                        print('Surface 3D:', self)
                        print('3D primitive in red:', primitive3d)
                        print('Previous 3D primitive:', last_primitive3d)
                        raise ValueError(
                            'Primitives not following each other in contour:',
                            'delta1={}, {}, {} ; '
                            'delta2={}, {}, {}'.format(
                                delta_x1, delta_y1, dist1,
                                delta_x2, delta_y2, dist2))

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
                try:
                    primitives3d.extend(getattr(self, method_name)(primitive2d))
                except NotImplementedError:
                    print('Error NotImplementedError')
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
        control_points = [self.point3d_to_2d(p)
                          for p in bspline_curve3d.control_points]
        return [vme.BSplineCurve2D(
                    bspline_curve3d.degree,
                    control_points=control_points,
                    knot_multiplicities=bspline_curve3d.knot_multiplicities,
                    knots=bspline_curve3d.knots,
                    weights=bspline_curve3d.weights,
                    periodic=bspline_curve3d.periodic)]

    def bsplinecurve2d_to_3d(self, bspline_curve2d):
        """
        Is this right?
        """
        control_points = [self.point2d_to_3d(p)
                          for p in bspline_curve2d.control_points]
        return [vme.BSplineCurve3D(
                    bspline_curve2d.degree,
                    control_points=control_points,
                    knot_multiplicities=bspline_curve2d.knot_multiplicities,
                    knots=bspline_curve2d.knots,
                    weights=bspline_curve2d.weights,
                    periodic=bspline_curve2d.periodic)]

    def normal_from_point2d(self, point2d):

        raise NotImplementedError('NotImplemented')

    def normal_from_point3d(self, point3d):
        """
        evaluates the normal vector of the bspline surface at this point3d.
        """

        return (self.normal_from_point2d(self.point3d_to_2d(point3d)))[1]

    def geodesic_distance_from_points2d(self, point1_2d: volmdlr.Point2D,
                                        point2_2d: volmdlr.Point2D, number_points: int = 50):
        """
        Approximation of geodesic distance via linesegments length sum in 3D
        """
        # points = [point1_2d]
        current_point3d = self.point2d_to_3d(point1_2d)
        distance = 0.
        for i in range(number_points):
            next_point3d = self.point2d_to_3d(point1_2d + (i + 1) / (number_points) * (point2_2d - point1_2d))
            distance += next_point3d.point_distance(current_point3d)
            current_point3d = next_point3d
        return distance

    def geodesic_distance(self, point1_3d: volmdlr.Point3D, point2_3d: volmdlr.Point3D):
        """
        Approximation of geodesic distance between 2 3D points supposed to be on the surface
        """
        point1_2d = self.point3d_to_2d(point1_3d)
        point2_2d = self.point3d_to_2d(point2_3d)
        return self.geodesic_distance_from_points2d(point1_2d, point2_2d)


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

    def to_dict(self, use_pointers: bool = True, memo=None, path: str = '#'):
        # improve the object structure ?
        dict_ = dc.DessiaObject.base_dict(self)
        dict_['frame'] = self.frame.to_dict(use_pointers=use_pointers, memo=memo, path=path + '/frame')
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
        return content, [plane_id]

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

            origin = points[0]
            vector1 = points[1] - origin
            vector1.normalize()
            vector2_min = points[2] - origin
            vector2_min.normalize()
            dot_min = abs(vector1.dot(vector2_min))
            for point in points[3:]:
                vector2 = point - origin
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
        u = line.point2 - line.point1
        w = line.point1 - self.frame.origin
        if math.isclose(self.frame.w.dot(u), 0, abs_tol=1e-08):
            return []
        intersection_abscissea = - self.frame.w.dot(w) / self.frame.w.dot(u)
        return [line.point1 + intersection_abscissea * u]

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

        # point2 = point1 + line_direction
        # return volmdlr.Line3D(point1, point2)
        return volmdlr.Line3D(point1, point1 + line_direction)

    def is_coincident(self, plane2):
        """
        Verifies if two planes are parallel and coincident
        """
        if self.frame.w.is_colinear_to(plane2.frame.w):
            if plane2.point_on_plane(self.frame.origin):
                return True
        return False

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

    def copy(self, deep=True, memo=None):
        new_frame = self.frame.copy()
        return Plane3D(new_frame, self.name)

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        else:
            fig = ax.figure

        self.frame.origin.plot(ax)
        self.frame.u.plot(ax, starting_point=self.frame.origin, color='r')
        self.frame.v.plot(ax, starting_point=self.frame.origin, color='g')
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
        control_points = [self.point3d_to_2d(p)
                          for p in bspline_curve3d.control_points]
        return [vme.BSplineCurve2D(
            bspline_curve3d.degree,
            control_points=control_points,
            knot_multiplicities=bspline_curve3d.knot_multiplicities,
            knots=bspline_curve3d.knots,
            weights=bspline_curve3d.weights,
            periodic=bspline_curve3d.periodic)]

    def bsplinecurve2d_to_3d(self, bspline_curve2d):
        control_points = [self.point2d_to_3d(p)
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
        if math.isclose(theta1, theta2, abs_tol=1e-9):
            return [vme.LineSegment3D(
                self.point2d_to_3d(linesegment2d.start),
                self.point2d_to_3d(linesegment2d.end),
            )]
        elif math.isclose(z1, z2, abs_tol=1e-9):
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
            # TODO: this is a non exact method!
            return [vme.LineSegment3D(self.point2d_to_3d(linesegment2d.start), self.point2d_to_3d(linesegment2d.end))]
            # raise NotImplementedError('Ellipse? delta_theta={} delta_z={}'.format(abs(theta2-theta1), abs(z1-z2)))

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
        points = [self.point3d_to_2d(bspline_curve3d.point_at_abscissa(i / 10 * l))
                  for i in range(11)]
        return [vme.LineSegment2D(p1, p2)
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
        return content, [current_id]

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

    def grid3d(self, points_x, points_y, xmin, xmax, ymin, ymax):
        '''
        generate 3d grid points of a Cylindrical surface, based on a 2d grid points parameters
        (xmin,xmax,points_x) limits and number of points in x,
        (ymin,ymax,points_y) limits and number of points in y
        '''

        points_2d = volmdlr.Point2D.grid2d(points_x, points_y, xmin, xmax, ymin, ymax)

        points_3d = []
        for j in range(0, len(points_2d)):
            points_3d.append(self.point2d_to_3d(points_2d[j]))

        return points_3d


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
        return content, [current_id]

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
        p2 = volmdlr.Point2D(theta2, phi1)
        p3 = volmdlr.Point2D(theta2, phi2)
        p4 = volmdlr.Point2D(theta1, phi2)
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
                center = self.frame.origin + self.R * u
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
        return content, [current_id]

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
        # points = []
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
        p2 = volmdlr.Point2D(theta2, phi1)
        p3 = volmdlr.Point2D(theta2, phi2)
        p4 = volmdlr.Point2D(theta1, phi2)
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
    _non_serializable_attributes = ['surface']

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
        # surface_points = surface.evalpts

        self.surface = surface
        # self.points = [volmdlr.Point3D(*p) for p in surface_points]
        volmdlr.core.Primitive3D.__init__(self, name=name)

        # Hidden Attributes
        self._displacements = ()
        self._grids2d = ()
        self._grids2d_deformed = ()

    @property
    def x_periodicity(self):
        p3d_x1 = self.point2d_to_3d(volmdlr.Point2D(1., 0.5))
        p2d_x0 = self.point3d_to_2d(p3d_x1, 0., 0.5)
        if self.point2d_to_3d(p2d_x0) == p3d_x1 and \
                not math.isclose(p2d_x0.x, 1, abs_tol=1e-3):
            return 1 - p2d_x0.x
        else:
            return None

    @property
    def y_periodicity(self):
        p3d_y1 = self.point2d_to_3d(volmdlr.Point2D(0.5, 1))
        p2d_y0 = self.point3d_to_2d(p3d_y1, 0., 0.5)
        if self.point2d_to_3d(p2d_y0) == p3d_y1 and \
                not math.isclose(p2d_y0.y, 1, abs_tol=1e-3):
            return 1 - p2d_y0.y
        else:
            return None

    def control_points_matrix(self, coordinates):
        '''
        define control points like a matrix, for each coordinate: x:0, y:1, z:2
        '''

        P = npy.empty((self.nb_u, self.nb_v))
        for i in range(0, self.nb_u):
            for j in range(0, self.nb_v):
                P[i][j] = self.control_points_table[i][j][coordinates]
        return P

    # Knots_vector
    def knots_vector_u(self):
        '''
        compute the global knot vector (u direction) based on knot elements and multiplicities
        '''

        knots = self.u_knots
        multiplicities = self.u_multiplicities

        knots_vec = []
        for i in range(0, len(knots)):
            for j in range(0, multiplicities[i]):
                knots_vec.append(knots[i])
        return knots_vec

    def knots_vector_v(self):
        '''
        compute the global knot vector (v direction) based on knot elements and multiplicities
        '''

        knots = self.v_knots
        multiplicities = self.v_multiplicities

        knots_vec = []
        for i in range(0, len(knots)):
            for j in range(0, multiplicities[i]):
                knots_vec.append(knots[i])
        return knots_vec

    def basis_functions_u(self, u, k, i):
        '''
        compute basis functions Bi in u direction for u=u and degree=k
        '''

        # k = self.degree_u
        t = self.knots_vector_u()

        if k == 0:
            return 1.0 if t[i] <= u < t[i + 1] else 0.0
        if t[i + k] == t[i]:
            c1 = 0.0
        else:
            c1 = (u - t[i]) / (t[i + k] - t[i]) * self.basis_functions_u(u, k - 1, i)
        if t[i + k + 1] == t[i + 1]:
            c2 = 0.0
        else:
            c2 = (t[i + k + 1] - u) / (t[i + k + 1] - t[i + 1]) * self.basis_functions_u(u, k - 1, i + 1)
        return c1 + c2

    def basis_functions_v(self, v, k, i):
        '''
        compute basis functions Bi in v direction for v=v and degree=k
        '''

        # k = self.degree_u
        t = self.knots_vector_v()

        if k == 0:
            return 1.0 if t[i] <= v < t[i + 1] else 0.0
        if t[i + k] == t[i]:
            c1 = 0.0
        else:
            c1 = (v - t[i]) / (t[i + k] - t[i]) * self.basis_functions_v(v, k - 1, i)
        if t[i + k + 1] == t[i + 1]:
            c2 = 0.0
        else:
            c2 = (t[i + k + 1] - v) / (t[i + k + 1] - t[i + 1]) * self.basis_functions_v(v, k - 1, i + 1)
        return c1 + c2

    def blending_vector_u(self, u):
        '''
        compute a vector of basis_functions in u direction for u=u
        '''

        blending_vect = npy.empty((1, self.nb_u))
        for j in range(0, self.nb_u):
            blending_vect[0][j] = self.basis_functions_u(u, self.degree_u, j)

        return blending_vect

    def blending_vector_v(self, v):
        '''
        compute a vector of basis_functions in v direction for v=v
        '''

        blending_vect = npy.empty((1, self.nb_v))
        for j in range(0, self.nb_v):
            blending_vect[0][j] = self.basis_functions_v(v, self.degree_v, j)

        return blending_vect

    def blending_matrix_u(self, u):
        '''
        compute a matrix of basis_functions in u direction for a vector u like [0,1]
        '''

        blending_mat = npy.empty((len(u), self.nb_u))
        for i in range(0, len(u)):
            for j in range(0, self.nb_u):
                blending_mat[i][j] = self.basis_functions_u(u[i], self.degree_u, j)
        return blending_mat

    def blending_matrix_v(self, v):
        '''
        compute a matrix of basis_functions in v direction for a vector v like [0,1]
        '''

        blending_mat = npy.empty((len(v), self.nb_v))
        for i in range(0, len(v)):
            for j in range(0, self.nb_v):
                blending_mat[i][j] = self.basis_functions_v(v[i], self.degree_v, j)
        return blending_mat

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        x, y = point2d
        if -1e-3 < x < 0:
            x = 0.
        elif 1 < x < 1 + 1e-3:
            x = 1
        if -1e-3 < y < 0:
            y = 0
        elif 1 < y < 1 + 1e-3:
            y = 1
        return volmdlr.Point3D(*self.surface.evaluate_single((x, y)))

    def point3d_to_2d(self, point3d: volmdlr.Point3D, min_bound_x: float = 0.,
                      max_bound_x: float = 1., min_bound_y: float = 0.,
                      max_bound_y: float = 1., tol=1e-9):
        def f(x):
            p3d = self.point2d_to_3d(volmdlr.Point2D(x[0], x[1]))
            return point3d.point_distance(p3d)

        results = []

        delta_bound_x = max_bound_x - min_bound_x
        delta_bound_y = max_bound_y - min_bound_y
        x0s = [((min_bound_x + max_bound_x) / 2, (min_bound_y + max_bound_y) / 2),
               (min_bound_x + delta_bound_x / 10, min_bound_y + delta_bound_y / 10),
               (min_bound_x + delta_bound_x / 10, max_bound_y - delta_bound_y / 10),
               (max_bound_x - delta_bound_x / 10, min_bound_y + delta_bound_y / 10),
               (max_bound_x - delta_bound_x / 10, max_bound_y - delta_bound_y / 10)]

        for x0 in x0s:
            z = scp.optimize.least_squares(f, x0=x0, bounds=([min_bound_x,
                                                              min_bound_y],
                                                             [max_bound_x,
                                                              max_bound_y]),
                                           ftol=tol / 10,
                                           xtol=tol / 10,
                                           # loss='soft_l1'
                                           )
            # z.cost represent the value of the cost function at the solution
            if z.fun < tol:
                return volmdlr.Point2D(*z.x)

            res = scp.optimize.minimize(f, x0=npy.array(x0),
                                        bounds=[(min_bound_x, max_bound_x),
                                                (min_bound_y, max_bound_y)],
                                        tol=tol)
            # res.fun represent the value of the objective function
            if res.fun < tol:
                return volmdlr.Point2D(*res.x)

            results.append((z.x, z.fun))
            results.append((res.x, res.fun))
        return (volmdlr.Point2D(*min(results, key=lambda r: r[1])[0]))

    def linesegment2d_to_3d(self, linesegment2d):
        # TODO: this is a non exact method!
        lth = linesegment2d.length()
        points = [self.point2d_to_3d(
            linesegment2d.point_at_abscissa(i * lth / 10.)) for i in range(11)]

        linesegment = vme.LineSegment3D(points[0], points[-1])
        flag = True
        for pt in points:
            if not linesegment.point_belongs(pt):
                flag = False
                break

        periodic = False
        if self.x_periodicity is not None and \
                math.isclose(lth, self.x_periodicity, abs_tol=1e-6) and \
                math.isclose(linesegment2d.start.y, linesegment2d.end.y,
                             abs_tol=1e-6):
            periodic = True
        elif self.y_periodicity is not None and \
                math.isclose(lth, self.y_periodicity, abs_tol=1e-6) and \
                math.isclose(linesegment2d.start.x, linesegment2d.end.x,
                             abs_tol=1e-6):
            periodic = True

        if flag:
            # All the points are on the same LineSegment3D
            linesegments = [linesegment]
        else:
            linesegments = [vme.BSplineCurve3D.from_points_interpolation(
                points, max(self.degree_u, self.degree_v), periodic=periodic)]
            # linesegments = [vme.LineSegment3D(p1, p2)
            #                 for p1, p2 in zip(points[:-1], points[1:])]
        return linesegments

    def linesegment3d_to_2d(self, linesegment3d):
        """
        a line segment on a BSplineSurface3D will be in any case a line in 2D?
        """
        x_perio = self.x_periodicity if self.x_periodicity is not None else 1.
        y_perio = self.y_periodicity if self.y_periodicity is not None else 1.
        return [vme.LineSegment2D(self.point3d_to_2d(linesegment3d.start,
                                                     max_bound_x=x_perio,
                                                     max_bound_y=y_perio),
                                  self.point3d_to_2d(linesegment3d.end,
                                                     max_bound_x=x_perio,
                                                     max_bound_y=y_perio))]

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        # TODO: enhance this, it is a non exact method!
        # TODO: bsplinecurve can be periodic but not around the bsplinesurface
        bsc_linesegment = vme.LineSegment3D(bspline_curve3d.points[0],
                                            bspline_curve3d.points[-1])
        flag = True
        for pt in bspline_curve3d.points:
            if not bsc_linesegment.point_belongs(pt):
                flag = False
                break

        if self.x_periodicity and not self.y_periodicity \
                and bspline_curve3d.periodic:
            p1 = self.point3d_to_2d(bspline_curve3d.points[0], min_bound_x=0.,
                                    max_bound_x=self.x_periodicity)
            p1_sup = self.point3d_to_2d(bspline_curve3d.points[0],
                                        min_bound_x=1 - self.x_periodicity)
            new_x = p1.x - p1_sup.x + self.x_periodicity
            new_x = new_x if 0 <= new_x else 0
            reverse = False
            if new_x < 0:
                new_x = 0
            elif math.isclose(new_x, self.x_periodicity, abs_tol=1e-5):
                new_x = 0
                reverse = True

            linesegments = [
                vme.LineSegment2D(
                    volmdlr.Point2D(new_x, p1.y),
                    volmdlr.Point2D(self.x_periodicity, p1.y))]
            if reverse:
                linesegments[0] = linesegments[0].reverse()

        elif self.y_periodicity and not self.x_periodicity \
                and bspline_curve3d.periodic:
            p1 = self.point3d_to_2d(bspline_curve3d.points[0], min_bound_y=0.,
                                    max_bound_y=self.y_periodicity)
            p1_sup = self.point3d_to_2d(bspline_curve3d.points[0],
                                        min_bound_y=1 - self.y_periodicity)
            new_y = p1.y - p1_sup.y + self.y_periodicity
            new_y = new_y if 0 <= new_y else 0
            reverse = False
            if new_y < 0:
                new_y = 0
            elif math.isclose(new_y, self.y_periodicity, abs_tol=1e-5):
                new_y = 0
                reverse = True

            linesegments = [
                vme.LineSegment2D(
                    volmdlr.Point2D(p1.x, new_y),
                    volmdlr.Point2D(p1.x, self.y_periodicity))]
            if reverse:
                linesegments[0] = linesegments[0].reverse()

        elif self.x_periodicity and self.y_periodicity \
                and bspline_curve3d.periodic:
            raise NotImplementedError

        elif flag:
            x_perio = self.x_periodicity if self.x_periodicity is not None\
                else 1.
            y_perio = self.y_periodicity if self.y_periodicity is not None\
                else 1.
            p1 = self.point3d_to_2d(bspline_curve3d.points[0],
                                    max_bound_x=x_perio,
                                    max_bound_y=y_perio)
            p2 = self.point3d_to_2d(bspline_curve3d.points[-1],
                                    max_bound_x=x_perio,
                                    max_bound_y=y_perio)

            if p1 == p2:
                print('BSplineCruve3D skipped because it is too small')
                linesegments = None
            else:
                p1_sup = self.point3d_to_2d(bspline_curve3d.points[0],
                                            min_bound_x=1 - x_perio,
                                            min_bound_y=1 - y_perio)
                p2_sup = self.point3d_to_2d(bspline_curve3d.points[-1],
                                            min_bound_x=1 - x_perio,
                                            min_bound_y=1 - y_perio)
                if self.x_periodicity and p1.point_distance(p1_sup) > 1e-5:
                    p1.x -= p1_sup.x - x_perio
                    p2.x -= p2_sup.x - x_perio
                if self.y_periodicity and p1.point_distance(p1_sup) > 1e-5:
                    p1.y -= p1_sup.y - y_perio
                    p2.y -= p2_sup.y - y_perio
                linesegments = [vme.LineSegment2D(p1, p2)]
            # How to check if end of surface overlaps start or the opposite ?
        else:
            lth = bspline_curve3d.length()
            if lth > 1e-5:
                points = [self.point3d_to_2d(
                        bspline_curve3d.point_at_abscissa(i / 10 * lth)
                        # max_bound_x=self.x_periodicity,
                        # max_bound_y=self.y_periodicity
                ) for i in range(11)]
                # linesegments = [vme.LineSegment2D(p1, p2)
                #                 for p1, p2 in zip(points[:-1], points[1:])]
                linesegments = [vme.BSplineCurve2D.from_points_interpolation(
                    points, max(self.degree_u, self.degree_v))]
            elif 1e-6 < lth <= 1e-5:
                linesegments = [vme.LineSegment2D(
                    self.point3d_to_2d(bspline_curve3d.start),
                    self.point3d_to_2d(bspline_curve3d.end))]
            else:
                print('BSplineCruve3D skipped because it is too small')
                linesegments = None

        # print(bspline_curve3d.start, bspline_curve3d.end)
        # print([(l.start, l.end) for l in linesegments])
        # print()
        return linesegments

    def arc3d_to_2d(self, arc3d):
        number_points = math.ceil(arc3d.angle * 10) + 1  # 10 points per radian
        l = arc3d.length()
        points = [self.point3d_to_2d(arc3d.point_at_abscissa(
            i * l / (number_points - 1))) for i in range(number_points)]
        # return [vme.LineSegment2D(p1, p2)
        #         for p1, p2 in zip(points[:-1], points[1:])]
        return [vme.BSplineCurve2D.from_points_interpolation(
                    points, max(self.degree_u, self.degree_v))]

    def arc2d_to_3d(self, arc2d):
        number_points = math.ceil(arc2d.angle * 7) + 1  # 7 points per radian
        l = arc2d.length()
        points = [self.point2d_to_3d(arc2d.point_at_abscissa(
            i * l / (number_points - 1))) for i in range(number_points)]
        return [vme.BSplineCurve3D.from_points_interpolation(
                    points, max(self.degree_u, self.degree_v))]

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
        return BSplineFace3D(self, surface, name)  # PlaneFace3D

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
            weight_data = [
                float(i) for i in
                arguments[13][1:-1].replace("(", "").replace(")", "").split(",")
            ]
        else:
            weight_data = None

        bsplinesurface = cls(degree_u, degree_v, control_points, nb_u, nb_v,
                             u_multiplicities, v_multiplicities, u_knots,
                             v_knots, weight_data, name)
        # if u_closed:
        #     bsplinesurface.x_periodicity = bsplinesurface.get_x_periodicity()
        # if v_closed:
        #     bsplinesurface.y_periodicity = bsplinesurface.get_y_periodicity()
        return bsplinesurface

    def to_step(self, current_id):
        content = ''
        point_matrix_ids = '('
        for points in self.control_points_table:
            point_ids = '('
            for point in points:
                point_content, point_id = point.to_step(current_id)
                content += point_content
                point_ids += '#{},'.format(point_id)
                current_id = point_id + 1
            point_ids = point_ids[:-1]
            point_ids += '),'
            point_matrix_ids += point_ids
        point_matrix_ids = point_matrix_ids[:-1]
        point_matrix_ids += ')'

        u_close = '.T.' if self.x_periodicity else '.F.'
        v_close = '.T.' if self.y_periodicity else '.F.'

        content += "#{} = B_SPLINE_SURFACE_WITH_KNOTS('{}',{},{},{},.UNSPECIFIED.,{},{},.F.,{},{},{},{},.UNSPECIFIED.);\n" \
            .format(current_id, self.name, self.degree_u, self.degree_v,
                    point_matrix_ids, u_close, v_close,
                    tuple(self.u_multiplicities), tuple(self.v_multiplicities),
                    tuple(self.u_knots), tuple(self.v_knots))
        return content, [current_id]

    def grid3d(self, points_x, points_y, xmin, xmax, ymin, ymax):
        '''
        generate 3d grid points of a Bspline surface, based on a 2d grid points parameters
        (xmin,xmax,points_x) limits and number of points in x,
        (ymin,ymax,points_y) limits and number of points in y
        '''

        points_2d = volmdlr.Point2D.grid2d(points_x, points_y, xmin, xmax, ymin, ymax)
        if not self._grids2d:
            self._grids2d = ([points_x, points_y, xmin, xmax, ymin, ymax], points_2d)

        points_3d = []
        for j in range(0, len(points_2d)):
            points_3d.append(self.point2d_to_3d(points_2d[j]))

        return points_3d

    def grid2d_deformed(self, points_x, points_y, xmin, xmax, ymin, ymax):
        '''
        dimension and deform a 2d grid points based on a Bspline surface
        (xmin,xmax,points_x) limits and number of points in x,
        (ymin,ymax,points_y) limits and number of points in y

        '''

        points_2d = volmdlr.Point2D.grid2d(points_x, points_y, xmin, xmax, ymin, ymax)
        points_3d = self.grid3d(points_x, points_y, xmin, xmax, ymin, ymax)

        # Parameters
        index_x = {}  # grid point position(i,j), x coordinates position in X(unknown variable)
        index_y = {}  # grid point position(i,j), y coordinates position in X(unknown variable)
        index_points = {}  # grid point position(j,i), point position in points_2d (or points_3d)
        k, p = 0, 0
        for i in range(0, points_x):
            for j in range(0, points_y):
                index_x.update({(j, i): k})
                index_y.update({(j, i): k + 1})
                index_points.update({(j, i): p})
                k = k + 2
                p = p + 1

        equation_points = []  # points combination to compute distances between 2D and 3D grid points
        # for i in range(0,points_y): #row from (0,i)
        #     for j in range(1,points_x):
        #         equation_points.append(((0,i),(j,i)))
        # for i in range(0,points_x): #column from (i,0)
        #     for j in range(1,points_y):
        #         equation_points.append(((i,0),(i,j)))
        for i in range(0, points_y):  # row
            for j in range(0, points_x - 1):
                equation_points.append(((j, i), (j + 1, i)))
        for i in range(0, points_x):  # column
            for j in range(0, points_x - 1):
                equation_points.append(((i, j), (i, j + 1)))
        for i in range(0, points_y - 1):  # diagonal
            for j in range(0, points_x - 1):
                equation_points.append(((j, i), (j + 1, i + 1)))

        for i in range(0, points_y):  # row 2segments (before.point.after)
            for j in range(1, points_x - 1):
                equation_points.append(((j - 1, i), (j + 1, i)))

        for i in range(0, points_x):  # column 2segments (before.point.after)
            for j in range(1, points_y - 1):
                equation_points.append(((i, j - 1), (i, j + 1)))

        # Euclidean distance
        # D=[] # distances between 3D grid points (based on points combination [equation_points])
        # for i in range(0, len(equation_points)):
        #     D.append((points_3d[index_points[equation_points[i][0]]].point_distance(points_3d[index_points[equation_points[i][1]]]))**2)

        # Geodesic distance
        # xx=[]
        # for p in points_2d:
        #     xx.append(p.x)
        # yy=[]
        # for p in points_2d:
        #     yy.append(p.y)

        # triang = plt_tri.Triangulation(xx, yy)
        # faces = triang.triangles
        # points = npy.empty([len(points_3d),3])
        # for i in range(0,len(points_3d)):
        #     points[i] = npy.array([points_3d[i].x,points_3d[i].y,points_3d[i].z])

        # geoalg = geodesic.PyGeodesicAlgorithmExact(points, faces)
        D = []  # geodesic distances between 3D grid points (based on points combination [equation_points])
        for i in range(0, len(equation_points)):
            D.append((self.geodesic_distance(
                points_3d[index_points[equation_points[i][0]]], points_3d[index_points[equation_points[i][1]]]))**2)

        # System of nonlinear equations
        def non_linear_equations(X):
            F = npy.empty(len(equation_points) + 2)
            for i in range(0, len(equation_points)):
                F[i] = abs((X[index_x[equation_points[i][0]]]**2 +
                            X[index_x[equation_points[i][1]]]**2 +
                            X[index_y[equation_points[i][0]]]**2 +
                            X[index_y[equation_points[i][1]]]**2 -
                            2 *
                            X[index_x[equation_points[i][0]]] *
                            X[index_x[equation_points[i][1]]] -
                            2 *
                            X[index_y[equation_points[i][0]]] *
                            X[index_y[equation_points[i][1]]] -
                            D[i]) /
                           D[i])

            F[i + 1] = X[0] * 1000
            F[i + 2] = X[1] * 1000
            # i=i+2
            # # F[i+3] = X[(len(points_2d)-points_x)*2]
            # l= 3
            # for f in range(1, points_x):
            #     F[i+f] = X[l]*1000
            #     l = l+2
            ## F[i+3] = X[3]*1000
            ## F[i+4] = X[5]*1000
            ## F[i+4] = X[points_x*2]*1000

            return F

        # Solution with "least_squares"
        x_init = []  # initial guess (2D grid points)
        for i in range(0, len(points_2d)):
            x_init.append(points_2d[i][0])
            x_init.append(points_2d[i][1])
        z = opt.least_squares(non_linear_equations, x_init)

        points_2d_deformed = []  # deformed 2d grid points
        for i in range(0, len(z.x), 2):
            points_2d_deformed.append(volmdlr.Point2D(z.x[i], z.x[i + 1]))

        self._grids2d_deformed = ([points_x, points_y, xmin, xmax, ymin, ymax], points_2d_deformed)

        return points_2d_deformed

    def grid2d_deformation(self, points_x, points_y, xmin, xmax, ymin, ymax):
        '''
        compute the deformation/displacement (dx/dy) of a 2d grid points based on a Bspline surface with:
        (xmin,xmax,points_x) limits and number of points in x,
        (ymin,ymax,points_y) limits and number of points in y

        '''

        points_2d = volmdlr.Point2D.grid2d(points_x, points_y, xmin, xmax, ymin, ymax)

        if [points_x, points_y, xmin, xmax, ymin, ymax] in self._grids2d_deformed:
            points_2d_deformed = self._grids2d_deformed[1]
        else:
            points_2d_deformed = self.grid2d_deformed(points_x, points_y, xmin, xmax, ymin, ymax)
            self._grids2d_deformed = ([points_x, points_y, xmin, xmax, ymin, ymax], points_2d_deformed)

        # Displacement,Deformation dx/dy
        displacement = npy.ones(shape=(len(points_2d), 2))  # 2D grid points displacement
        for i in range(0, len(displacement)):
            displacement[i][0] = points_2d_deformed[i][0] - points_2d[i][0]
            displacement[i][1] = points_2d_deformed[i][1] - points_2d[i][1]

        self._displacements = ([points_x, points_y, xmin, xmax, ymin, ymax], displacement)

        return displacement

    def point2d_parametric_to_dimension(self, point2d: volmdlr.Point3D, points_x, points_y, xmin, xmax, ymin, ymax):
        '''
        convert a point2d from the parametric to the dimensioned frame
        '''

        # Check if the 0<point2d.x<1 and 0<point2d.y<1
        if point2d.x < 0:
            point2d.x = 0
        elif point2d.x > 1:
            point2d.x = 1
        if point2d.y < 0:
            point2d.y = 0
        elif point2d.y > 1:
            point2d.y = 1

        if [points_x, points_y, xmin, xmax, ymin, ymax] in self._grids2d:
            points_2d = self._grids2d[1]
        else:
            points_2d = volmdlr.Point2D.grid2d(points_x, points_y, xmin, xmax, ymin, ymax)
            self._grids2d = ([points_x, points_y, xmin, xmax, ymin, ymax], points_2d)

        if [points_x, points_y, xmin, xmax, ymin, ymax] in self._displacements:
            displacement = self._displacements[1]
        else:
            displacement = self.grid2d_deformation(points_x, points_y, xmin, xmax, ymin, ymax)
            self._displacements = ([points_x, points_y, xmin, xmax, ymin, ymax], displacement)

        # Parameters
        index_points = {}  # grid point position(j,i), point position in points_2d (or points_3d)
        p = 0
        for i in range(0, points_x):
            for j in range(0, points_y):
                index_points.update({(j, i): p})
                p = p + 1

        # Form function "Finite Elements"
        def form_function(s, t):
            N = npy.empty(4)
            N[0] = (1 - s) * (1 - t) / 4
            N[1] = (1 + s) * (1 - t) / 4
            N[2] = (1 + s) * (1 + t) / 4
            N[3] = (1 - s) * (1 + t) / 4
            return N

        finite_elements_points = []  # 2D grid points index that define one element
        for j in range(0, points_y - 1):
            for i in range(0, points_x - 1):
                finite_elements_points.append(((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)))
        finite_elements = []  # finite elements defined with closed polygon
        for i in range(0, len(finite_elements_points)):
            finite_elements.append(volmdlr.wires.ClosedPolygon2D((points_2d[index_points[finite_elements_points[i][0]]],
                                                                  points_2d[index_points[finite_elements_points[i][1]]],
                                                                  points_2d[index_points[finite_elements_points[i][2]]],
                                                                  points_2d[index_points[finite_elements_points[i][3]]])))

        for k in range(0, len(finite_elements_points)):
            if (volmdlr.wires.Contour2D(finite_elements[k].primitives).point_belongs(point2d)  # finite_elements[k].point_belongs(point2d)
                or volmdlr.wires.Contour2D(finite_elements[k].primitives).point_over_contour(point2d)
                or ((points_2d[index_points[finite_elements_points[k][0]]][0] < point2d.x < points_2d[index_points[finite_elements_points[k][1]]][0])
                    and point2d.y == points_2d[index_points[finite_elements_points[k][0]]][1])
                or ((points_2d[index_points[finite_elements_points[k][1]]][1] < point2d.y < points_2d[index_points[finite_elements_points[k][2]]][1])
                    and point2d.x == points_2d[index_points[finite_elements_points[k][1]]][0])
                or ((points_2d[index_points[finite_elements_points[k][3]]][0] < point2d.x < points_2d[index_points[finite_elements_points[k][2]]][0])
                    and point2d.y == points_2d[index_points[finite_elements_points[k][1]]][1])
                or ((points_2d[index_points[finite_elements_points[k][0]]][1] < point2d.y < points_2d[index_points[finite_elements_points[k][3]]][1])
                    and point2d.x == points_2d[index_points[finite_elements_points[k][0]]][0])):
                break

        x0 = points_2d[index_points[finite_elements_points[k][0]]][0]
        y0 = points_2d[index_points[finite_elements_points[k][0]]][1]
        x1 = points_2d[index_points[finite_elements_points[k][1]]][0]
        y2 = points_2d[index_points[finite_elements_points[k][2]]][1]
        x = point2d.x
        y = point2d.y
        s = 2 * ((x - x0) / (x1 - x0)) - 1
        t = 2 * ((y - y0) / (y2 - y0)) - 1

        N = form_function(s, t)
        dx = npy.array([displacement[index_points[finite_elements_points[k][0]]][0],
                        displacement[index_points[finite_elements_points[k][1]]][0],
                        displacement[index_points[finite_elements_points[k][2]]][0],
                        displacement[index_points[finite_elements_points[k][3]]][0]])
        dy = npy.array([displacement[index_points[finite_elements_points[k][0]]][1],
                        displacement[index_points[finite_elements_points[k][1]]][1],
                        displacement[index_points[finite_elements_points[k][2]]][1],
                        displacement[index_points[finite_elements_points[k][3]]][1]])

        return volmdlr.Point2D(point2d.x + npy.transpose(N).dot(dx), point2d.y + npy.transpose(N).dot(dy))

    def point3d_to_2d_with_dimension(self, point3d: volmdlr.Point3D, points_x, points_y, xmin, xmax, ymin, ymax):
        '''
        compute the point2d of a point3d, on a Bspline surface, in the dimensioned frame
        '''

        point2d = self.point3d_to_2d(point3d)

        point2d_with_dimension = self.point2d_parametric_to_dimension(
            point2d, points_x, points_y, xmin, xmax, ymin, ymax)

        return point2d_with_dimension

    def point2d_with_dimension_to_parametric_frame(self, point2d, points_x, points_y, xmin, xmax, ymin, ymax):
        '''
        convert a point2d from the dimensioned to the parametric frame
        '''

        if [points_x, points_y, xmin, xmax, ymin, ymax] in self._grids2d:
            points_2d = self._grids2d[1]
        else:
            points_2d = volmdlr.Point2D.grid2d(points_x, points_y, xmin, xmax, ymin, ymax)
            self._grids2d = ([points_x, points_y, xmin, xmax, ymin, ymax], points_2d)

        if [points_x, points_y, xmin, xmax, ymin, ymax] in self._grids2d_deformed:
            points_2d_deformed = self._grids2d_deformed[1]
        else:
            points_2d_deformed = self.grid2d_deformed(points_x, points_y, xmin, xmax, ymin, ymax)
            self._grids2d_deformed = ([points_x, points_y, xmin, xmax, ymin, ymax], points_2d_deformed)

        # Parameters
        index_points = {}  # grid point position(j,i), point position in points_2d (or points_3d)
        p = 0
        for i in range(0, points_x):
            for j in range(0, points_y):
                index_points.update({(j, i): p})
                p = p + 1

        finite_elements_points = []  # 2D grid points index that define one element
        for j in range(0, points_y - 1):
            for i in range(0, points_x - 1):
                finite_elements_points.append(((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)))
        finite_elements = []  # finite elements defined with closed polygon  DEFORMED
        for i in range(0, len(finite_elements_points)):
            finite_elements.append(volmdlr.wires.ClosedPolygon2D((points_2d_deformed[index_points[finite_elements_points[i][0]]],
                                                                  points_2d_deformed[index_points[finite_elements_points[i][1]]],
                                                                  points_2d_deformed[index_points[finite_elements_points[i][2]]],
                                                                  points_2d_deformed[index_points[finite_elements_points[i][3]]])))

        finite_elements_initial = []  # finite elements defined with closed polygon  INITIAL
        for i in range(0, len(finite_elements_points)):
            finite_elements_initial.append(volmdlr.wires.ClosedPolygon2D((points_2d[index_points[finite_elements_points[i][0]]],
                                                                          points_2d[index_points[finite_elements_points[i][1]]],
                                                                          points_2d[index_points[finite_elements_points[i][2]]],
                                                                          points_2d[index_points[finite_elements_points[i][3]]])))

        for k in range(0, len(finite_elements_points)):
            if (finite_elements[k].point_belongs(point2d)
                or ((points_2d_deformed[index_points[finite_elements_points[k][0]]][0] < point2d.x < points_2d_deformed[index_points[finite_elements_points[k][1]]][0])
                    and point2d.y == points_2d_deformed[index_points[finite_elements_points[k][0]]][1])
                or ((points_2d_deformed[index_points[finite_elements_points[k][1]]][1] < point2d.y < points_2d_deformed[index_points[finite_elements_points[k][2]]][1])
                    and point2d.x == points_2d_deformed[index_points[finite_elements_points[k][1]]][0])
                or ((points_2d_deformed[index_points[finite_elements_points[k][3]]][0] < point2d.x < points_2d_deformed[index_points[finite_elements_points[k][2]]][0])
                    and point2d.y == points_2d_deformed[index_points[finite_elements_points[k][1]]][1])
                or ((points_2d_deformed[index_points[finite_elements_points[k][0]]][1] < point2d.y < points_2d_deformed[index_points[finite_elements_points[k][3]]][1])
                    and point2d.x == points_2d_deformed[index_points[finite_elements_points[k][0]]][0])
                or finite_elements[k].primitives[0].point_belongs(point2d) or finite_elements[k].primitives[1].point_belongs(point2d)
                    or finite_elements[k].primitives[2].point_belongs(point2d) or finite_elements[k].primitives[3].point_belongs(point2d)):

                break

        frame_deformed = volmdlr.Frame2D(finite_elements[k].center_of_mass(),
                                         volmdlr.Vector2D(finite_elements[k].primitives[1].middle_point()[0] - finite_elements[k].center_of_mass()[0],
                                                          finite_elements[k].primitives[1].middle_point()[1] - finite_elements[k].center_of_mass()[1]),
                                         volmdlr.Vector2D(finite_elements[k].primitives[0].middle_point()[0] - finite_elements[k].center_of_mass()[0],
                                                          finite_elements[k].primitives[0].middle_point()[1] - finite_elements[k].center_of_mass()[1]))

        point2d_frame_deformed = volmdlr.Point2D(point2d.frame_mapping(frame_deformed, 'new')[0],
                                                 point2d.frame_mapping(frame_deformed, 'new')[1])

        frame_inital = volmdlr.Frame2D(finite_elements_initial[k].center_of_mass(),
                                       volmdlr.Vector2D(finite_elements_initial[k].primitives[1].middle_point()[0] - finite_elements_initial[k].center_of_mass()[0],
                                                        finite_elements_initial[k].primitives[1].middle_point()[1] - finite_elements_initial[k].center_of_mass()[1]),
                                       volmdlr.Vector2D(finite_elements_initial[k].primitives[0].middle_point()[0] - finite_elements_initial[k].center_of_mass()[0],
                                                        finite_elements_initial[k].primitives[0].middle_point()[1] - finite_elements_initial[k].center_of_mass()[1]))

        X = point2d_frame_deformed.frame_mapping(frame_inital, 'old')[0]
        if X < 0:
            X = 0
        elif X > 1:
            X = 1
        Y = point2d_frame_deformed.frame_mapping(frame_inital, 'old')[1]
        if Y < 0:
            Y = 0
        elif Y > 1:
            Y = 1

        return volmdlr.Point2D(X, Y)

    def point2d_with_dimension_to_3d(self, point2d, points_x, points_y, xmin, xmax, ymin, ymax):
        '''
        compute the point3d, on a Bspline surface, of a point2d define in the dimensioned frame
        '''

        point2d_01 = self.point2d_with_dimension_to_parametric_frame(
            point2d, points_x, points_y, xmin, xmax, ymin, ymax)

        return self.point2d_to_3d(point2d_01)

    def linesegment2d_parametric_to_dimension(self, linesegment2d, points_x, points_y):
        '''
        convert a linesegment2d from the parametric to the dimensioned frame
        '''

        xmin, xmax, ymin, ymax = 0, 1, 0, 1

        points = linesegment2d.discretization_points(20)
        points_dim = [
            self.point2d_parametric_to_dimension(
                p,
                points_x,
                points_y,
                xmin,
                xmax,
                ymin,
                ymax) for p in points]

        return vme.BSplineCurve2D.from_points_interpolation(
                points_dim, max(self.degree_u, self.degree_v))

    def linesegment3d_to_2d_with_dimension(self, linesegment3d, points_x, points_y):
        '''
        compute the linesegment2d of a linesegment3d, on a Bspline surface, in the dimensioned frame
        '''

        linesegment2d = self.linesegment3d_to_2d(linesegment3d)
        bsplinecurve2d_with_dimension = self.linesegment2d_parametric_to_dimension(linesegment2d, points_x, points_y)

        return bsplinecurve2d_with_dimension

    def linesegment2d_with_dimension_to_parametric_frame(self, linesegment2d):
        '''
        convert a linesegment2d from the dimensioned to the parametric frame
        '''

        [points_x, points_y, xmin, xmax, ymin, ymax] = self._grids2d[0]

        try:
            linesegment2d = volmdlr.edges.LineSegment2D(
                self.point2d_with_dimension_to_parametric_frame(
                    linesegment2d.start, points_x, points_y, xmin, xmax, ymin, ymax), self.point2d_with_dimension_to_parametric_frame(
                    linesegment2d.end, points_x, points_y, xmin, xmax, ymin, ymax))
        except NotImplementedError:
            return None

        return linesegment2d

    def linesegment2d_with_dimension_to_3d(self, linesegment2d):
        '''
        compute the linesegment3d, on a Bspline surface, of a linesegment2d defined in the dimensioned frame
        '''

        linesegment2d_01 = self.linesegment2d_with_dimension_to_parametric_frame(linesegment2d)
        linesegment3d = self.linesegment2d_to_3d(linesegment2d_01)

        return linesegment3d

    def bsplinecurve2d_parametric_to_dimension(self, bsplinecurve2d, points_x, points_y):
        '''
        convert a bsplinecurve2d from the parametric to the dimensioned frame
        '''

        xmin, xmax, ymin, ymax = 0, 1, 0, 1
        # check if bsplinecurve2d is in a list
        if isinstance(bsplinecurve2d, list):
            bsplinecurve2d = bsplinecurve2d[0]
        points = bsplinecurve2d.control_points
        points_dim = []

        for p in points:
            points_dim.append(self.point2d_parametric_to_dimension(p, points_x, points_y, xmin, xmax, ymin, ymax))

        bsplinecurve2d_with_dimension = volmdlr.edges.BSplineCurve2D(bsplinecurve2d.degree, points_dim,
                                                                     bsplinecurve2d.knot_multiplicities,
                                                                     bsplinecurve2d.knots,
                                                                     bsplinecurve2d.weights,
                                                                     bsplinecurve2d.periodic)

        return bsplinecurve2d_with_dimension

    def bsplinecurve3d_to_2d_with_dimension(self, bsplinecurve3d, points_x, points_y):
        '''
        compute the bsplinecurve2d of a bsplinecurve3d, on a Bspline surface, in the dimensioned frame
        '''

        bsplinecurve2d_01 = self.bsplinecurve3d_to_2d(bsplinecurve3d)
        bsplinecurve2d_with_dimension = self.bsplinecurve2d_parametric_to_dimension(
            bsplinecurve2d_01, points_x, points_y)

        return bsplinecurve2d_with_dimension

    def bsplinecurve2d_with_dimension_to_parametric_frame(self, bsplinecurve2d):
        '''
        convert a bsplinecurve2d from the dimensioned to the parametric frame
        '''

        points_dim = bsplinecurve2d.control_points
        points = []
        [points_x, points_y, xmin, xmax, ymin, ymax] = self._grids2d[0]

        for p in points_dim:
            points.append(
                self.point2d_with_dimension_to_parametric_frame(
                    p, points_x, points_y, xmin, xmax, ymin, ymax))

        bsplinecurve2d = volmdlr.edges.BSplineCurve2D(bsplinecurve2d.degree, points,
                                                      bsplinecurve2d.knot_multiplicities,
                                                      bsplinecurve2d.knots,
                                                      bsplinecurve2d.weights,
                                                      bsplinecurve2d.periodic)
        return bsplinecurve2d

    def bsplinecurve2d_with_dimension_to_3d(self, bsplinecurve2d):
        '''
        compute the bsplinecurve3d, on a Bspline surface, of a bsplinecurve2d defined in the dimensioned frame
        '''

        bsplinecurve2d_01 = self.bsplinecurve2d_with_dimension_to_parametric_frame(bsplinecurve2d)
        bsplinecurve3d = self.bsplinecurve2d_to_3d(bsplinecurve2d_01)

        return bsplinecurve3d

    def arc2d_parametric_to_dimension(self, arc2d, points_x, points_y):
        '''
        convert a arc2d from the parametric to the dimensioned frame
        '''

        xmin, xmax, ymin, ymax = 0, 1, 0, 1

        number_points = math.ceil(arc2d.angle * 7) + 1
        l = arc2d.length()
        points = [self.point2d_parametric_to_dimension(arc2d.point_at_abscissa(
            i * l / (number_points - 1)), points_x, points_y, xmin, xmax, ymin, ymax) for i in range(number_points)]

        return vme.BSplineCurve2D.from_points_interpolation(
                    points, max(self.degree_u, self.degree_v))

    def arc3d_to_2d_with_dimension(self, arc3d, points_x, points_y):
        '''
        compute the arc2d of a arc3d, on a Bspline surface, in the dimensioned frame
        '''

        bsplinecurve2d = self.arc3d_to_2d(arc3d)[0]  # it's a bsplinecurve2d
        arc2d_with_dimension = self.bsplinecurve2d_parametric_to_dimension(bsplinecurve2d, points_x, points_y)

        return arc2d_with_dimension  # it's a bsplinecurve2d-dimension

    def arc2d_with_dimension_to_parametric_frame(self, arc2d):
        '''
        convert a arc2d from the dimensioned to the parametric frame
        '''

        [points_x, points_y, xmin, xmax, ymin, ymax] = self._grids2d[0]

        number_points = math.ceil(arc2d.angle * 7) + 1
        l = arc2d.length()

        points = [self.point2d_with_dimension_to_parametric_frame(arc2d.point_at_abscissa(
                i * l / (number_points - 1)), points_x, points_y, xmin, xmax, ymin, ymax) for i in range(number_points)]

        return vme.BSplineCurve2D.from_points_interpolation(
                    points, max(self.degree_u, self.degree_v))

    def arc2d_with_dimension_to_3d(self, arc2d):
        '''
        compute the  arc3d, on a Bspline surface, of a arc2d in the dimensioned frame
        '''

        arc2d_01 = self.arc2d_with_dimension_to_parametric_frame(arc2d)
        arc3d = self.arc2d_to_3d(arc2d_01)

        return arc3d  # it's a bsplinecurve3d

    def contour2d_parametric_to_dimension(self, contour2d: volmdlr.wires.Contour2D, points_x, points_y):
        '''
        convert a contour2d from the parametric to the dimensioned frame
        '''

        primitives2d_dim = []

        for primitive2d in contour2d.primitives:
            # method_name = '{}_parametric_to_dimension'.format(
            #     primitive2d.__class__.__name__.lower())
            method_name = f'{primitive2d.__class__.__name__.lower()}_parametric_to_dimension'

            if hasattr(self, method_name):
                primitives = getattr(self, method_name)(primitive2d, points_x, points_y)
                if primitives:
                    primitives2d_dim.append(primitives)

            else:
                raise NotImplementedError(
                    f'Class {self.__class__.__name__} does not implement {method_name}')

        return volmdlr.wires.Contour2D(primitives2d_dim)

    def contour3d_to_2d_with_dimension(self, contour3d: volmdlr.wires.Contour3D, points_x, points_y):
        '''
        compute the contou2d of a contour3d, on a Bspline surface, in the dimensioned frame
        '''

        contour2d_01 = self.contour3d_to_2d(contour3d)

        return self.contour2d_parametric_to_dimension(contour2d_01, points_x, points_y)

    def contour2d_with_dimension_to_parametric_frame(self, contour2d):
        '''
        convert a contour2d from the dimensioned to the parametric frame
        '''

        # TODO: check and avoid primitives with start=end
        primitives2d = []

        for primitive2d in contour2d.primitives:
            method_name = f'{primitive2d.__class__.__name__.lower()}_with_dimension_to_parametric_frame'

            if hasattr(self, method_name):
                primitives = getattr(self, method_name)(primitive2d)
                if primitives:
                    primitives2d.append(primitives)

            else:
                raise NotImplementedError(
                    # 'Class {} does not implement {}'.format(self.__class__.__name__,
                    #                                         method_name))
                    f'Class {self.__class__.__name__} does not implement {method_name}')

        # #Avoid to have primitives with start=end
        # start_points = []
        # for i in range(0, len(new_start_points)-1):
        #     if new_start_points[i] != new_start_points[i+1]:
        #         start_points.append(new_start_points[i])
        # if new_start_points[-1] != new_start_points[0]:
        #     start_points.append(new_start_points[-1])

        return volmdlr.wires.Contour2D(primitives2d)

    def contour2d_with_dimension_to_3d(self, contour2d):
        '''
        compute the contour3d, on a Bspline surface, of a contour2d define in the dimensioned frame
        '''

        contour01 = self.contour2d_with_dimension_to_parametric_frame(contour2d)

        return self.contour2d_to_3d(contour01)

    @classmethod
    def from_geomdl_surface(cls, surface):
        '''
        create a volmdlr's BSpline_Surface3D from a geomdl's one
        '''

        control_points = []
        for i in range(0, len(surface.ctrlpts)):
            control_points.append(volmdlr.Point3D(surface.ctrlpts[i][0], surface.ctrlpts[i][1], surface.ctrlpts[i][2]))

        (u_knots, u_multiplicities) = knots_vector_inv((surface.knotvector_u))
        (v_knots, v_multiplicities) = knots_vector_inv((surface.knotvector_v))

        bspline_surface = cls(degree_u=surface.degree_u,
                              degree_v=surface.degree_v,
                              control_points=control_points,
                              nb_u=surface.ctrlpts_size_u,
                              nb_v=surface.ctrlpts_size_v,
                              u_multiplicities=u_multiplicities,
                              v_multiplicities=v_multiplicities,
                              u_knots=u_knots,
                              v_knots=v_knots)

        return bspline_surface

    @classmethod
    def points_fitting_into_bspline_surface(cls, points_3d, size_u, size_v, degree_u, degree_v):
        '''
        Bspline Surface interpolation through 3d points

        Parameters
        ----------
        points_3d : volmdlr.Point3D
            data points
        size_u : int
            number of data points on the u-direction.
        size_v : int
            number of data points on the v-direction.
        degree_u : int
            degree of the output surface for the u-direction.
        degree_v : int
            degree of the output surface for the v-direction.

        Returns
        -------
        B-spline surface

        '''

        points = []
        for i in range(0, len(points_3d)):
            points.append((points_3d[i].x, points_3d[i].y, points_3d[i].z))

        surface = interpolate_surface(points, size_u, size_v, degree_u, degree_v)

        return cls.from_geomdl_surface(surface)

    @classmethod
    def points_approximate_into_bspline_surface(cls, points_3d, size_u, size_v, degree_u, degree_v, **kwargs):
        '''
        Bspline Surface approximate through 3d points

        Parameters
        ----------
        points_3d : volmdlr.Point3D
            data points
        size_u : int
            number of data points on the u-direction.
        size_v : int
            number of data points on the v-direction.
        degree_u : int
            degree of the output surface for the u-direction.
        degree_v : int
            degree of the output surface for the v-direction.

        Keyword Arguments:
            * ``ctrlpts_size_u``: number of control points on the u-direction. *Default: size_u - 1*
            * ``ctrlpts_size_v``: number of control points on the v-direction. *Default: size_v - 1*

        Returns
        -------
        B-spline surface: volmdlr.faces.BSplineSurface3D

        '''

        # Keyword arguments
        num_cpts_u = kwargs.get('ctrlpts_size_u', size_u - 1)  # number of datapts, r + 1 > number of ctrlpts, n + 1
        num_cpts_v = kwargs.get('ctrlpts_size_v', size_v - 1)  # number of datapts, s + 1 > number of ctrlpts, m + 1

        points = [tuple([*pt]) for pt in points_3d]

        surface = approximate_surface(points, size_u, size_v, degree_u, degree_v,
                                      ctrlpts_size_u=num_cpts_u, num_cpts_v=num_cpts_v)

        return cls.from_geomdl_surface(surface)

    @classmethod
    def from_cylindrical_faces(cls, cylindrical_faces, degree_u, degree_v,
                               points_x: int = 10, points_y: int = 10):
        '''
        define a bspline surface from a list of cylindrical faces

        Parameters
        ----------
        cylindrical_faces : List[volmdlr.faces.CylindricalFace3D]
            faces 3d
        degree_u : int
            degree of the output surface for the u-direction
        degree_v : int
            degree of the output surface for the v-direction
        points_x : int
            number of points in x-direction
        points_y : int
            number of points in y-direction

        Returns
        -------
        B-spline surface

        '''

        if len(cylindrical_faces) == 1:

            return cls.from_cylindrical_face(cylindrical_faces[0], degree_u, degree_v, 50, 50)

        if len(cylindrical_faces) > 1:
            bspline_surfaces = []
            direction = cylindrical_faces[0].adjacent_direction(cylindrical_faces[1])

            if direction == 'x':
                bounding_rectangle_0 = cylindrical_faces[0].surface2d.outer_contour.bounding_rectangle()
                ymin = bounding_rectangle_0[2]
                ymax = bounding_rectangle_0[3]
                for face in cylindrical_faces:
                    bounding_rectangle = face.surface2d.outer_contour.bounding_rectangle()
                    ymin = min(ymin, bounding_rectangle[2])
                    ymax = max(ymax, bounding_rectangle[3])
                for face in cylindrical_faces:
                    bounding_rectangle = face.surface2d.outer_contour.bounding_rectangle()

                    points_3d = face.surface3d.grid3d(points_x, points_y,
                                                      bounding_rectangle[0], bounding_rectangle[1],
                                                      ymin, ymax)
                    bspline_surfaces.append(
                        cls.points_fitting_into_bspline_surface(
                            points_3d, points_x, points_y, degree_u, degree_v))

            elif direction == 'y':
                bounding_rectangle_0 = cylindrical_faces[0].surface2d.outer_contour.bounding_rectangle()
                xmin = bounding_rectangle_0[0]
                xmax = bounding_rectangle_0[1]
                for face in cylindrical_faces:
                    bounding_rectangle = face.surface2d.outer_contour.bounding_rectangle()
                    xmin = min(xmin, bounding_rectangle[0])
                    xmax = max(xmax, bounding_rectangle[1])
                for face in cylindrical_faces:
                    bounding_rectangle = face.surface2d.outer_contour.bounding_rectangle()

                    points_3d = face.surface3d.grid3d(points_x, points_y, xmin, xmax,
                                                      bounding_rectangle[2], bounding_rectangle[3])
                    bspline_surfaces.append(
                        cls.points_fitting_into_bspline_surface(
                            points_3d, points_x, points_y, degree_u, degree_v))

            to_be_merged = bspline_surfaces[0]
            for i in range(0, len(bspline_surfaces) - 1):
                merged = to_be_merged.merge_with(bspline_surfaces[i + 1])
                to_be_merged = merged

            bspline_surface = to_be_merged

            return bspline_surface

    @classmethod
    # points_x: int = 50, points_y: int = 50):
    def from_cylindrical_face(cls, cylindrical_face, degree_u, degree_v, **kwargs):
        '''
        define a bspline surface from a cylindrical face

        Parameters
        ----------
        cylindrical_face : volmdlr.faces.CylindricalFace3D
            face 3d
        degree_u : int
            degree of the output surface for the u-direction.
        degree_v : int
            degree of the output surface for the v-direction.
        points_x : int
            number of points in x-direction
        points_y : int
            number of points in y-direction

        Returns
        -------
        B-spline surface

        '''

        points_x = kwargs['points_x']
        points_y = kwargs['points_y']
        bounding_rectangle = cylindrical_face.surface2d.outer_contour.bounding_rectangle()
        points_3d = cylindrical_face.surface3d.grid3d(points_x, points_y,
                                                      bounding_rectangle[0],
                                                      bounding_rectangle[1],
                                                      bounding_rectangle[2],
                                                      bounding_rectangle[3])

        return cls.points_fitting_into_bspline_surface(points_3d, points_x, points_x, degree_u, degree_v)

    def intersection_with(self, other_bspline_surface3d):
        '''
        compute intersection points between two Bspline surfaces
        return u,v parameters for intersection points for both surfaces
        '''

        def f(X):
            return (self.point2d_to_3d(volmdlr.Point2D(X[0], X[1])) -
                    other_bspline_surface3d.point2d_to_3d(volmdlr.Point2D(X[2], X[3]))).norm()

        x = npy.linspace(0, 1, 10)
        x_init = []
        for xi in x:
            for yi in x:
                x_init.append((xi, yi, xi, yi))

        u1, v1, u2, v2 = [], [], [], []
        solutions = []
        for x0 in x_init:
            z = scp.optimize.least_squares(f, x0=x0, bounds=([0, 1]))
            # print(z.cost)
            if z.fun < 1e-5:
                solution = z.x
                if solution not in solutions:
                    solutions.append(solution)
                    u1.append(solution[0])
                    v1.append(solution[1])
                    u2.append(solution[2])
                    v2.append(solution[3])

        # uv1 = [[min(u1),max(u1)],[min(v1),max(v1)]]
        # uv2 = [[min(u2),max(u2)],[min(v2),max(v2)]]

        return ((u1, v1), (u2, v2))  # (uv1, uv2)

    def plane_intersection(self, plane3d):
        '''
        compute intersection points between a Bspline surface and a plane3d
        '''

        def f(X):
            return ((self.surface.evaluate_single((X[0], X[1]))[0]) * plane3d.equation_coefficients()[0] +
                    (self.surface.evaluate_single((X[0], X[1]))[1]) * plane3d.equation_coefficients()[1] +
                    (self.surface.evaluate_single((X[0], X[1]))[2]) * plane3d.equation_coefficients()[2] +
                    plane3d.equation_coefficients()[3])

        x = npy.linspace(0, 1, 20)
        x_init = []
        for xi in x:
            for yi in x:
                x_init.append((xi, yi))

        # x_init = volmdlr.Point2D.grid2d(20, 20, 0, 1, 0, 1)

        intersection_points = []
        # solutions = []
        # u, v =[],  []

        for x0 in x_init:
            z = scp.optimize.least_squares(f, x0=x0, bounds=([0, 1]))
            if z.fun < 1e-20:
                #     cost.append(z.cost)
                # # print(z.cost)
                # if z.cost<1e-20:
                solution = z.x
                intersection_points.append(volmdlr.Point3D(self.surface.evaluate_single((solution[0], solution[1]))[0],
                                                           self.surface.evaluate_single((solution[0], solution[1]))[1],
                                                           self.surface.evaluate_single((solution[0], solution[1]))[2]))
        # intersection_points.sort()
                # u.append(solution[0])
                # v.append(solution[1])
                # solutions.append(solution)

        # return (u,v)
        return intersection_points

    def error_with_point3d(self, point3d):
        '''
        compute the error/distance between the Bspline surface and a point3d
        '''

        def f(x):
            return (point3d - self.point2d_to_3d(volmdlr.Point2D(x[0], x[1]))).norm()

        cost = []

        for x0 in [(0, 0), (0, 1), (1, 0), (1, 1), (0.5, 0.5)]:
            z = scp.optimize.least_squares(f, x0=x0, bounds=([0, 1]))
            cost.append(z.fun)

        return min(cost)

    def error_with_edge3d(self, edge3d):
        '''
        compute the error/distance between the Bspline surface and an edge3d
        it's the mean of the start and end points errors'
        '''

        return (self.error_with_point3d(edge3d.start) + self.error_with_point3d(edge3d.end)) / 2

    def nearest_edges3d(self, contour3d, threshold: float):
        '''
        compute the nearest edges of a contour3d to a Bspline_surface3d based on a threshold
        '''

        nearest = []
        for primitive in contour3d.primitives:
            if self.error_with_edge3d(primitive) <= threshold:
                nearest.append(primitive)
        nearest_primitives = volmdlr.wires.Wire3D(nearest)

        return nearest_primitives

    def edge3d_to_2d_with_dimension(self, edge3d, points_x, points_y):
        '''
        compute the edge2d of a edge3d, on a Bspline surface, in the dimensioned frame
        '''

        # method_name = '{}_to_2d_with_dimension'.format(edge3d.__class__.__name__.lower())
        method_name = f'{edge3d.__class__.__name__.lower()}_to_2d_with_dimension'

        if hasattr(self, method_name):
            edge2d_dim = getattr(self, method_name)(edge3d, points_x, points_y)
            if edge2d_dim:
                return edge2d_dim
            else:
                raise NotImplementedError
        else:
            raise NotImplementedError(
                # 'Class {} does not implement {}'.format(self.__class__.__name__,
                #                                         method_name))
                f'Class {self.__class__.__name__} does not implement {method_name}')

    def wire3d_to_2d(self, wire3d):
        '''
        compute the 2d of a wire3d, on a Bspline surface
        '''

        contour = self.contour3d_to_2d(wire3d)

        return volmdlr.wires.Wire2D(contour.primitives)

    def wire3d_to_2d_with_dimension(self, wire3d):
        '''
        compute the 2d of a wire3d, on a Bspline surface, in the dimensioned frame
        '''

        points_x = self._grids2d[0][0]
        points_y = self._grids2d[0][1]

        contour = self.contour3d_to_2d_with_dimension(wire3d, points_x, points_y)

        return volmdlr.wires.Wire2D(contour.primitives)

    def split_surface_u(self, u: float):
        '''
        split the surface at the input parametric coordinate on the u-direction

        Parameters
        ----------
        u : float
            Parametric coordinate u choosen between 0 and 1

        Returns
        -------
        surfaces : list
            Two splitted surfaces

        '''

        surfaces_geo = split_surface_u(self.surface, u)
        surfaces = []
        for s in surfaces_geo:
            surfaces.append(volmdlr.faces.BSplineSurface3D.from_geomdl_surface(s))

        return surfaces

    def split_surface_v(self, v: float):
        '''
        split the surface at the input parametric coordinate on the v-direction

        Parameters
        ----------
        v : float
            Parametric coordinate v choosen between 0 and 1

        Returns
        -------
        surfaces : list
            Two splitted surfaces

        '''

        surfaces_geo = split_surface_v(self.surface, v)
        surfaces = []
        for s in surfaces_geo:
            surfaces.append(volmdlr.faces.BSplineSurface3D.from_geomdl_surface(s))

        return surfaces

    def split_surface_with_bspline_curve(self, bspline_curve3d: volmdlr.edges.BSplineCurve3D):
        '''
        cuts the surface into two pieces with a bspline curve

        Parameters
        ----------
        bspline_curve3d : volmdlr.edges.BSplineCurve3D


        Returns
        -------
        surfaces : list
            Two splitted surfaces

        '''

        surfaces = []
        bspline_curve2d = self.bsplinecurve3d_to_2d(bspline_curve3d)[0]
        # if type(bspline_curve2d) == list:
        #     points = [bspline_curve2d[0].start]
        #     for edge in bspline_curve2d:
        #         points.append(edge.end)
        #     bspline_curve2d = vme.BSplineCurve2D.from_points_approximation(points, 2, ctrlpts_size = 5)
        contour = self.rectangular_cut(0, 1, 0, 1).surface2d.outer_contour
        contours = contour.cut_by_bspline_curve(bspline_curve2d)

        du, dv = bspline_curve2d.end - bspline_curve2d.start
        resolution = 8

        for contour in contours:
            u_min, u_max, v_min, v_max = contour.bounding_rectangle()
            if du > dv:
                delta_u = u_max - u_min
                nlines_x = int(delta_u * resolution)
                lines_x = [vme.Line2D(volmdlr.Point2D(u_min, v_min),
                                      volmdlr.Point2D(u_min, v_max))]
                for i in range(nlines_x):
                    u = u_min + (i + 1) / (nlines_x + 1) * delta_u
                    lines_x.append(vme.Line2D(volmdlr.Point2D(u, v_min),
                                              volmdlr.Point2D(u, v_max)))
                lines_x.append(vme.Line2D(volmdlr.Point2D(u_max, v_min),
                                          volmdlr.Point2D(u_max, v_max)))
                lines = lines_x

            else:
                delta_v = v_max - v_min
                nlines_y = int(delta_v * resolution)
                lines_y = [vme.Line2D(volmdlr.Point2D(v_min, v_min),
                                      volmdlr.Point2D(v_max, v_min))]
                for i in range(nlines_y):
                    v = v_min + (i + 1) / (nlines_y + 1) * delta_v
                    lines_y.append(vme.Line2D(volmdlr.Point2D(v_min, v),
                                              volmdlr.Point2D(v_max, v)))
                lines_y.append(vme.Line2D(volmdlr.Point2D(v_min, v_max),
                                          volmdlr.Point2D(v_max, v_max)))
                lines = lines_y

            pt0 = volmdlr.O2D
            points = []

            for l in lines:
                inter = contour.line_intersections(l)
                if inter:
                    pt = [inter[0][0], inter[1][0]]
                else:
                    raise NotImplementedError

                pt = sorted(pt, key=lambda p: pt0.point_distance(p))
                pt0 = pt[0]
                edge = volmdlr.edges.LineSegment2D(pt[0], pt[1])

                points.extend(edge.discretization_points(10))

            points3d = []
            for p in points:
                points3d.append(self.point2d_to_3d(p))

            size_u, size_v, degree_u, degree_v = 10, 10, self.degree_u, self.degree_v
            surfaces.append(
                volmdlr.faces.BSplineSurface3D.points_fitting_into_bspline_surface(
                    points3d, size_u, size_v, degree_u, degree_v))

        return surfaces

    def point_belongs(self, point3d):
        '''
        check if a point3d belongs to the bspline_surface or not
        '''

        def f(x):
            p3d = self.point2d_to_3d(volmdlr.Point2D(x[0], x[1]))
            return point3d.point_distance(p3d)

        x = npy.linspace(0, 1, 5)
        x_init = []
        for xi in x:
            for yi in x:
                x_init.append((xi, yi))

        for x0 in x_init:
            z = scp.optimize.least_squares(f, x0=x0, bounds=([0, 1]))
            if z.fun < 1e-10:
                return True
        return False

    def is_intersected_with(self, other_bspline_surface3d):
        '''
        check if the two surfaces are intersected or not
        return True, when there are more 50points on the intersection zone
        '''

        # intersection_results = self.intersection_with(other_bspline_surface3d)
        # if len(intersection_results[0][0]) >= 50:
        #     return True
        # else:
        #     return False

        def f(X):
            return (self.point2d_to_3d(volmdlr.Point2D(X[0], X[1])) -
                    other_bspline_surface3d.point2d_to_3d(volmdlr.Point2D(X[2], X[3]))).norm()

        x = npy.linspace(0, 1, 10)
        x_init = []
        for xi in x:
            for yi in x:
                x_init.append((xi, yi, xi, yi))

        i = 0
        for x0 in x_init:
            z = scp.optimize.least_squares(f, x0=x0, bounds=([0, 1]))
            if z.fun < 1e-5:
                i += 1
                if i >= 50:
                    return True
        return False

    def merge_with(self, other_bspline_surface3d):
        '''
        merge two adjacent surfaces based on their faces

        Parameters
        ----------
        other_bspline_face3d : volmdlr.faces.BSplineSurface3D

        Returns
        -------
        merged_surface : volmdlr.faces.BSplineSurface3D

        '''

        bspline_face3d = self.rectangular_cut(0, 1, 0, 1)
        other_bspline_face3d = other_bspline_surface3d.rectangular_cut(0, 1, 0, 1)

        bsplines = [self, other_bspline_surface3d]
        bsplines_new = bsplines

        center = [bspline_face3d.surface2d.outer_contour.center_of_mass(),
                  other_bspline_face3d.surface2d.outer_contour.center_of_mass()]
        grid2d_direction = (bspline_face3d.pair_with(other_bspline_face3d))[1]

        if bspline_face3d.outer_contour3d.is_sharing_primitives_with(other_bspline_face3d.outer_contour3d):

            xmin, xmax, ymin, ymax = self.xy_limits(other_bspline_surface3d)

        elif self.is_intersected_with(other_bspline_surface3d):
            # find pimitives to split with
            contour1 = bspline_face3d.outer_contour3d
            contour2 = other_bspline_face3d.outer_contour3d

            distances = []
            for p1 in contour1.primitives:
                dis = []
                for p2 in contour2.primitives:
                    point1 = (p1.start + p1.end) / 2
                    point2 = (p2.start + p2.end) / 2
                    dis.append(point1.point_distance(point2))
                distances.append(dis)

            i = distances.index((min(distances)))
            j = distances[i].index(min(distances[i]))

            curves = [contour2.primitives[j], contour1.primitives[i]]

            # split surface
            for i, bspline in enumerate(bsplines):
                surfaces = bspline.split_surface_with_bspline_curve(curves[i])

                errors = []
                for s in surfaces:
                    errors.append(s.error_with_point3d(bsplines[i].point2d_to_3d(center[i])))

                bsplines_new[i] = surfaces[errors.index(min(errors))]

            xmin, xmax, ymin, ymax = [0] * len(bsplines_new), [1] * len(bsplines_new), [0] * \
                len(bsplines_new), [1] * len(bsplines_new)

            grid2d_direction = (
                bsplines_new[0].rectangular_cut(
                    0, 1, 0, 1).pair_with(
                    bsplines_new[1].rectangular_cut(
                        0, 1, 0, 1)))[1]

        else:
            xmin, xmax, ymin, ymax = [0] * len(bsplines_new), [1] * len(bsplines_new), [0] * \
                                               len(bsplines_new), [1] * len(bsplines_new)

        # grid3d
        points3d = []
        for i, bspline in enumerate(bsplines_new):
            grid2d = volmdlr.Point2D.grid2d_with_direction(
                50, 50, xmin[i], xmax[i], ymin[i], ymax[i], grid2d_direction[i])[0]
            grid3d = []
            for p in grid2d:
                grid3d.append(bspline.point2d_to_3d(p))

            points3d.extend(grid3d)

            # fitting
        size_u, size_v, degree_u, degree_v = 100, 50, max(
            bsplines[0].degree_u, bsplines[1].degree_u), max(
            bsplines[0].degree_v, bsplines[1].degree_v)

        merged_surface = volmdlr.faces.BSplineSurface3D.points_fitting_into_bspline_surface(
            points3d, size_u, size_v, degree_u, degree_v)

        return merged_surface

    def xy_limits(self, other_bspline_surface3d):
        '''
        compute x, y limits to define grid2d
        '''

        grid2d_direction = (
            self.rectangular_cut(
                0, 1, 0, 1).pair_with(
                other_bspline_surface3d.rectangular_cut(
                    0, 1, 0, 1)))[1]

        xmin, xmax, ymin, ymax = [], [], [], []
        if grid2d_direction[0][1] == '+y':
            xmin.append(0)
            xmax.append(1)
            ymin.append(0)
            ymax.append(0.99)
        elif grid2d_direction[0][1] == '+x':
            xmin.append(0)
            xmax.append(0.99)
            ymin.append(0)
            ymax.append(1)
        elif grid2d_direction[0][1] == '-x':
            xmin.append(0.01)
            xmax.append(1)
            ymin.append(0)
            ymax.append(1)
        elif grid2d_direction[0][1] == '-y':
            xmin.append(0)
            xmax.append(1)
            ymin.append(0.01)
            ymax.append(1)

        xmin.append(0)
        xmax.append(1)
        ymin.append(0)
        ymax.append(1)

        return xmin, xmax, ymin, ymax


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
            for i, subsurface2d in enumerate(subsurfaces2d):
                face = self.__class__(self.surface3d, subsurface2d)
                face_content, face_id = face.to_step_without_splitting(
                    current_id)
                face_ids.append(face_id[0])
                content += face_content
                current_id = face_id[0] + 1
            return content, face_ids
        else:
            return self.to_step_without_splitting(current_id)

    def to_step_without_splitting(self, current_id):
        content, surface3d_ids = self.surface3d.to_step(current_id)
        current_id = max(surface3d_ids) + 1

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
            surface3d_ids[0])
        # TODO: create an ADVANCED_FACE for each surface3d_ids ?
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
            # try:
            surfaces = self.surface2d.split_by_lines(lines_x)
            # except:
            #     self.plot()
            #     raise NotImplementedError
        elif lines_y:
            surfaces = self.surface2d.split_by_lines(lines_y)
        else:
            surfaces = [self.surface2d]

        # mesh2d = surfaces[0].triangulation()
        # print('ls', len(surfaces))
        # for subsurface in surfaces[1:]:
        #     # mesh2d += subsurface.triangulation()
        #     mesh2d.merge_mesh(subsurface.triangulation())

        meshes = [s.triangulation() for s in surfaces]
        mesh2d = vmd.DisplayMesh2D.merge_meshes(meshes)
        return vmd.DisplayMesh3D(
            [vmd.Node3D(*self.surface3d.point2d_to_3d(p)) for p in
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

    def copy(self, deep=True, memo=None):
        return self.__class__(self.surface3d.copy(), self.surface2d.copy(),
                              self.name)

    def line_intersections(self,
                           line: vme.Line3D,
                           ) -> List[volmdlr.Point3D]:
        intersections = []
        for intersection in self.surface3d.line_intersections(line):
            if self.point_belongs(intersection):
                intersections.append(intersection)

        return intersections

    def linesegment_intersections(self,
                                  linesegment: vme.LineSegment3D,
                                  ) -> List[volmdlr.Point3D]:
        intersections = []
        for intersection in self.surface3d.linesegment_intersections(
                linesegment):
            if self.point_belongs(intersection):
                intersections.append(intersection)

        return intersections

    def plot(self, ax=None, color='k', alpha=1, edge_details=False):
        if not ax:
            ax = plt.figure().add_subplot(111, projection='3d')
        self.outer_contour3d.plot(ax=ax, color=color, alpha=alpha,
                                  edge_details=edge_details)
        [contour3d.plot(ax=ax, color=color, alpha=alpha,
                        edge_details=edge_details)
         for contour3d in self.inner_contours3d]
        return ax

    def random_point_inside(self):
        point_inside2d = self.surface2d.random_point_inside()
        return self.surface3d.point2d_to_3d(point_inside2d)


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
    def dict_to_object(cls, dict_, global_dict=None, pointers_memo=None):
        plane3d = Plane3D.dict_to_object(dict_['surface3d'],
                                         global_dict=global_dict,
                                         pointers_memo=pointers_memo)
        surface2d = Surface2D.dict_to_object(dict_['surface2d'],
                                             global_dict=global_dict,
                                             pointers_memo=pointers_memo)
        return cls(plane3d, surface2d, dict_['name'])

    def area(self):
        return self.surface2d.outer_contour.area()

    def copy(self, deep=True, memo=None):
        return PlaneFace3D(self.surface3d.copy(), self.surface2d.copy(),
                           self.name)

    def _bounding_box(self):
        """
        """
        return self.outer_contour3d._bounding_box()

    def face_inside(self, face2):
        """
        verifies if a face is inside another face.
        It returns True if face2 is inside or False if the opposite
        """

        if self.surface3d.is_coincident(face2.surface3d):
            self_contour2d = self.outer_contour3d.to_2d(
                self.surface3d.frame.origin, self.surface3d.frame.u, self.surface3d.frame.v)
            face2_contour2d = face2.outer_contour3d.to_2d(
                self.surface3d.frame.origin, self.surface3d.frame.u, self.surface3d.frame.v)
            if self_contour2d.is_inside(face2_contour2d):
                # ax=self_contour2d.plot()
                # face2_contour2d.plot(ax=ax, color='r')
                return True
        return False

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

        projected_pt = point.plane_projection3d(self.surface3d.frame.origin,
                                                self.surface3d.frame.u,
                                                self.surface3d.frame.v)
        projection_distance = point.point_distance(projected_pt)

        if self.point_belongs(projected_pt):
            if return_other_point:
                return projection_distance, projected_pt
            return projection_distance

        point_2D = point.to_2d(self.surface3d.frame.origin, self.surface3d.frame.u,
                               self.surface3d.frame.v)

        polygon2D = self.surface2d.outer_contour.to_polygon(angle_resolution=10)
        border_distance, other_point = polygon2D.point_border_distance(point_2D, return_other_point=True)

        other_point = self.surface3d.point2d_to_3d(volmdlr.Point2D(*other_point))

        if return_other_point:
            return (projection_distance ** 2 + border_distance ** 2) ** 0.5, \
                   other_point
        return (projection_distance ** 2 + border_distance ** 2) ** 0.5

    def minimum_distance_points_plane(self, other_plane_face,
                                      return_points=False):
        # """
        # Only works if the surface is planar
        # TODO : this function does not take into account if Face has holes
        # TODO : TRAITER LE CAS OU LA DISTANCE LA PLUS COURTE N'EST PAS D'UN SOMMET
        # """
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
            if self.surface2d.point_belongs(
                    point2d) or self.surface2d.outer_contour.point_over_contour(point2d, abs_tol=1e-7):
                if surface3d_inter not in intersections:
                    intersections.append(surface3d_inter)

        return intersections

    def face_intersections(self, face2, tol=1e-8) -> List[volmdlr.wires.Wire3D]:
        # """
        # Only works if the surface is planar
        # TODO : this function does not take into account if Face has holes
        # """

        bbox1 = self.bounding_box
        bbox2 = face2.bounding_box
        if not bbox1.bbox_intersection(bbox2) and bbox1.distance_to_bbox(bbox2) >= tol:
            return []

        intersections = []

        for edge2 in face2.outer_contour3d.primitives:
            intersection_points = self.edge_intersections(edge2)
            if intersection_points:
                intersections.extend(intersection_points)
        for edge1 in self.outer_contour3d.primitives:
            intersection_points = face2.edge_intersections(edge1)
            if intersection_points:
                for pt in intersection_points:
                    if pt not in intersections:
                        intersections.append(pt)
        if len(intersections) > 1:
            if intersections[0] == intersections[1]:
                return []
            primitive = volmdlr.edges.LineSegment3D(intersections[0], intersections[1])
            mid_point = primitive.middle_point()

            if self.outer_contour3d.point_over_contour(mid_point) and\
                    face2.outer_contour3d.point_over_contour(mid_point):
                return []

            intersections = [volmdlr.wires.Wire3D([primitive])]
            return intersections
        return []

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

    def get_face_cutting_contours(self, dict_intersecting_combinations):
        '''
            :param face: A face object
            :param dict_intersecting_combinations: dictionary containing as keys the combination of intersecting faces
             and as the values the resulting primitive from the intersection of these two faces

            return a list all contours cutting one particular face
        '''
        face_intersecting_primitives2d = []
        for intersecting_combination in dict_intersecting_combinations.keys():
            if self in intersecting_combination:
                primitive2 = dict_intersecting_combinations[intersecting_combination].primitives[0]
                primitive2_2d = volmdlr.edges.LineSegment2D(
                    self.surface3d.point3d_to_2d(
                        primitive2.start), self.surface3d.point3d_to_2d(
                        primitive2.end))

                if not self.surface2d.outer_contour.primitive_over_contour(primitive2_2d, tol=1e-7):
                    face_intersecting_primitives2d.append(primitive2_2d)

        if not face_intersecting_primitives2d:
            return []

        list_cutting_contours = volmdlr.wires.Contour2D.contours_from_edges(face_intersecting_primitives2d[:])

        return list_cutting_contours

    def divide_face(self, list_cutting_contours, inside, intersection_method=True):
        '''
            :param list_cutting_contours: list of contours cutting the face
            :param inside: when extracting a contour from another contour. It defines the extracted contour as being between the two points if True and outside these points if False
            return a list new faces resulting from face division
        '''
        list_faces = []
        list_open_cutting_contours = []
        list_closed_cutting_contours = []
        for cutting_contour in list_cutting_contours:
            if cutting_contour.primitives[0].start != cutting_contour.primitives[-1].end:
                list_open_cutting_contours.append(cutting_contour)
            else:
                list_closed_cutting_contours.append(cutting_contour)

        if list_open_cutting_contours:
            new_faces_contours = self.surface2d.outer_contour.divide(list_open_cutting_contours, inside)
            for contour in new_faces_contours:
                list_faces.append(
                    PlaneFace3D(self.surface3d, Surface2D(contour, [])))

        if list_closed_cutting_contours:
            new_contour = list_closed_cutting_contours[0]
            if len(new_contour.primitives) >= 3 and new_contour.primitives[0].start == new_contour.primitives[-1].end:
                surf3d = self.surface3d
                surf2d = Surface2D(self.surface2d.outer_contour, [new_contour])
                new_plane = PlaneFace3D(surf3d, surf2d)
                list_faces.append(new_plane)
                list_faces.append(PlaneFace3D(surf3d, Surface2D(new_contour, [])))
            else:
                surf3d = self.surface3d
                surf2d = Surface2D(self.surface2d.outer_contour, [])
                new_plane = PlaneFace3D(surf3d, surf2d)
                list_faces.append(new_plane)

        return list_faces

    def is_adjacent(self, face2: Face3D):
        contour1 = self.outer_contour3d.to_2d(
            self.surface3d.frame.origin,
            self.surface3d.frame.u,
            self.surface3d.frame.v)
        contour2 = face2.outer_contour3d.to_2d(
            self.surface3d.frame.origin,
            self.surface3d.frame.u,
            self.surface3d.frame.v)
        if contour1.shares_primitives(contour2):
            return True
        return False

    @staticmethod
    def merge_faces(list_coincident_faces: List[Face3D]):
        valid_coicident_faces = list_coincident_faces[:]
        list_new_faces = []
        list_inner_contours = []
        merge_finished = False
        face0 = valid_coicident_faces[0]
        merged_contour = face0.outer_contour3d.to_2d(
            face0.surface3d.frame.origin,
            face0.surface3d.frame.u,
            face0.surface3d.frame.v)
        valid_coicident_faces.remove(face0)
        while not merge_finished:
            adjacent_faces = False
            list_inner_contours = []
            for face in valid_coicident_faces:
                adjacent_faces = False
                contour = face.outer_contour3d.to_2d(
                    face0.surface3d.frame.origin,
                    face0.surface3d.frame.u,
                    face0.surface3d.frame.v)
                if contour.shares_primitives(merged_contour):
                    merged_contour_results = merged_contour.merge_with(
                        contour)
                    merged_contour = merged_contour_results[0]
                    merged_inner_contours = merged_contour_results[1:]
                    list_inner_contours.extend(merged_inner_contours)
                    list_inner_contours.extend(face.surface2d.inner_contours)
                    valid_coicident_faces.remove(face)
                    adjacent_faces = True
                    break
            if not adjacent_faces and valid_coicident_faces:
                list_new_faces.append(
                    PlaneFace3D(face0.surface3d,
                                Surface2D(merged_contour.copy(),
                                          face0.surface2d.inner_contours +
                                          list_inner_contours)))
                merged_contour = \
                    valid_coicident_faces[0].outer_contour3d.to_2d(
                        face0.surface3d.frame.origin,
                        face0.surface3d.frame.u,
                        face0.surface3d.frame.v)
                valid_coicident_faces.remove(valid_coicident_faces[0])

            if not valid_coicident_faces:
                merge_finished = True
        list_new_faces.append(
            PlaneFace3D(face0.surface3d,
                        Surface2D(merged_contour,
                                  face0.surface2d.inner_contours +
                                  list_inner_contours)))

        return list_new_faces

    def set_operations_new_faces(self, intersecting_combinations,
                                 contour_extract_inside):
        list_cutting_contours = self.get_face_cutting_contours(
            intersecting_combinations)
        if not list_cutting_contours:
            return [self]
        return self.divide_face(list_cutting_contours, contour_extract_inside)


class Triangle3D(PlaneFace3D):
    """
    :param point1: The first point
    :type point1: volmdlr.Point3D
    :param point2: The second point
    :type point2: volmdlr.Point3D
    :param point3: The third point
    :type point3: volmdlr.Point3D
    """
    _standalone_in_db = False
    # _generic_eq = True
    # _non_serializable_attributes = ['bounding_box', 'polygon2D']
    # _non_eq_attributes = ['name', 'bounding_box', 'outer_contour3d',
    #                       'inner_contours3d']
    # _non_hash_attributes = []

    def __init__(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D,
                 point3: volmdlr.Point3D, alpha=1, color=None, name: str = ''):
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
        self.points = [self.point1, self.point2, self.point3]
        self.color = color
        self.alpha = alpha
        self.name = name

        self._utd_surface3d = False
        self._utd_surface2d = False
        self.bounding_box = self._bounding_box()

        dc.DessiaObject.__init__(self, name=name)

        # Don't use inheritence for performance: class method fakes face3D behavior
        # Face3D.__init__(self,
        #                 surface3d=plane3d,
        #                 surface2d=surface2d,
        #                 name=name)

    def _bounding_box(self):
        return volmdlr.core.BoundingBox.from_points([self.point1, self.point2, self.point3])

    @property
    def surface3d(self):
        if not self._utd_surface3d:
            self._surface3d = Plane3D.from_3_points(self.point1, self.point2, self.point3)
            self._utd_surface3d = True
        return self._surface3d

    @property
    def surface2d(self):
        if not self._utd_surface2d:
            plane3d = self.surface3d
            contour3d = volmdlr.wires.Contour3D([vme.LineSegment3D(self.point1, self.point2),
                                                 vme.LineSegment3D(self.point2, self.point3),
                                                 vme.LineSegment3D(self.point3, self.point1)])

            contour2d = contour3d.to_2d(plane3d.frame.origin,
                                        plane3d.frame.u, plane3d.frame.v)

            self._surface2d = Surface2D(outer_contour=contour2d, inner_contours=[])

            self._utd_surface2d = True
        return self._surface2d

    def to_dict(self, use_pointers: bool = False, memo=None, path: str = '#'):
        dict_ = dc.DessiaObject.base_dict(self)
        dict_['point1'] = self.point1.to_dict()
        dict_['point2'] = self.point2.to_dict()
        dict_['point3'] = self.point3.to_dict()
        dict_['name'] = self.name

        return dict_

    @classmethod
    def dict_to_object(cls, dict_):
        point1 = volmdlr.Point3D.dict_to_object(dict_['point1'])
        point2 = volmdlr.Point3D.dict_to_object(dict_['point2'])
        point3 = volmdlr.Point3D.dict_to_object(dict_['point3'])
        return cls(point1, point2, point3, dict_['name'])

    def area(self) -> float:
        """

        :return: area triangle
        :rtype: float

        Formula explained here: https://www.triangle-calculator.com/?what=vc

        """
        a = self.point1.point_distance(self.point2)
        b = self.point2.point_distance(self.point3)
        c = self.point3.point_distance(self.point1)

        semi_perimeter = (a + b + c) / 2

        try:
            # Area with Heron's formula
            area = math.sqrt(semi_perimeter * (semi_perimeter - a) * (semi_perimeter - b) * (semi_perimeter - c))
        except ValueError:
            area = 0

        return area

    def height(self):
        # Formula explained here: https://www.triangle-calculator.com/?what=vc
        # Basis = vector point1 to point2d
        return 2 * self.area() / self.point1.point_distance(self.point2)

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            np1 = self.point1.frame_mapping(frame, side, copy=True)
            np2 = self.point2.frame_mapping(frame, side, copy=True)
            np3 = self.point3.frame_mapping(frame, side, copy=True)
            return self.__class__(np1, np2, np3, self.name)
        else:
            self.point1.frame_mapping(frame, side, copy=False)
            self.point2.frame_mapping(frame, side, copy=False)
            self.point3.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()

    def copy(self, deep=True, memo=None):
        return Triangle3D(self.point1.copy(), self.point2.copy(), self.point3.copy(),
                          self.name)

    def triangulation(self):
        return vmd.DisplayMesh3D([vmd.Node3D.from_point(self.point1),
                                  vmd.Node3D.from_point(self.point2),
                                  vmd.Node3D.from_point(self.point3)],
                                 [(0, 1, 2)])

    def translation(self, offset, copy=True):
        new_point1 = self.point1.translation(offset, True)
        new_point2 = self.point2.translation(offset, True)
        new_point3 = self.point3.translation(offset, True)

        new_triangle = Triangle3D(new_point1, new_point2, new_point3,
                                  self.alpha, self.color, self.name)
        if copy:
            return new_triangle
        else:
            self.point1 = new_point1
            self.point2 = new_point2
            self.point3 = new_point3

    def rotation(self, center, axis, angle, copy=True):
        new_point1 = self.point1.rotation(center, axis, angle, copy=True)
        new_point2 = self.point2.rotation(center, axis, angle, copy=True)
        new_point3 = self.point3.rotation(center, axis, angle, copy=True)

        new_triangle = Triangle3D(new_point1, new_point2, new_point3,
                                  self.alpha, self.color, self.name)
        if copy:
            return new_triangle
        else:
            self.point1 = new_point1
            self.point2 = new_point2
            self.point3 = new_point3

    def subdescription(self, resolution=0.01):
        frame = self.surface3d.frame
        pts2d = [pt.to_2d(frame.origin, frame.u, frame.v) for pt in self.points]

        t_poly2d = volmdlr.wires.ClosedPolygon2D(pts2d)

        xmin, xmax = min(pt.x for pt in pts2d), max(pt.x for pt in pts2d)
        ymin, ymax = min(pt.y for pt in pts2d), max(pt.y for pt in pts2d)

        nbx, nby = int(((xmax - xmin) / resolution) + 2), int(((ymax - ymin) / resolution) + 2)
        points_box = []
        for i in range(nbx):
            x = min(xmin + i * resolution, xmax)
            if x == xmin:
                x = xmin + 0.01 * resolution
            for j in range(nby):
                y = min(ymin + j * resolution, ymax)
                if y == ymin:
                    y = ymin + 0.01 * resolution
                points_box.append(volmdlr.Point2D(x, y))

        points = [pt.copy() for pt in self.points]
        for pt in points_box:
            if t_poly2d.point_belongs(pt):
                points.append(pt.to_3d(frame.origin, frame.u, frame.v))
            elif t_poly2d.point_over_contour(pt):
                points.append(pt.to_3d(frame.origin, frame.u, frame.v))

        return volmdlr.Vector3D.remove_duplicate(points)

    def subdescription_to_triangles(self, resolution=0.01):
        """
        This function will return a list of Triangle3D with resolution as max
        length of subtriangles side
        """

        frame = self.surface3d.frame
        pts2d = [pt.to_2d(frame.origin, frame.u, frame.v) for pt in self.points]

        t_poly2d = volmdlr.wires.ClosedPolygon2D(pts2d)

        sub_triangles2d = [t_poly2d]
        done = False
        while not done:
            triangles2d = []
            for t, subtri in enumerate(sub_triangles2d):
                ls_length = [ls.length() for ls in subtri.line_segments]
                ls_max = max(ls_length)

                if ls_max > resolution:
                    pos_ls_max = ls_length.index(ls_max)
                    taller = subtri.line_segments[pos_ls_max]
                    p1, p2 = taller.start, taller.end
                    p3 = list(set(subtri.points) - set([p1, p2]))[0]

                    pt_mid = (p1 + p2) / 2
                    new_triangles2d = [volmdlr.wires.ClosedPolygon2D([p1, pt_mid, p3]),
                                       volmdlr.wires.ClosedPolygon2D([p2, pt_mid, p3])]

                    triangles2d.extend(new_triangles2d)
                else:
                    triangles2d.append(subtri)

            if len(sub_triangles2d) == len(triangles2d):
                done = True
                break
            sub_triangles2d = triangles2d

        triangles3d = [Triangle3D(tri.points[0].to_3d(frame.origin, frame.u, frame.v),
                                  tri.points[1].to_3d(frame.origin, frame.u, frame.v),
                                  tri.points[2].to_3d(frame.origin, frame.u, frame.v)) for tri in sub_triangles2d]

        return triangles3d

    def middle(self):
        return (self.point1 + self.point2 + self.point3) / 3

    def normal(self):
        '''

        Returns
        -------
        normal to the face

        '''
        normal = self.surface3d.frame.w
        # vec12 = self.point2 - self.point1
        # vec13 = self.point3 - self.point1
        # normal  = vec12.cross(vec13)
        normal.normalize()
        return normal


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
                 surface3d: CylindricalSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):

        self.radius = surface3d.radius
        self.center = surface3d.frame.origin
        self.normal = surface3d.frame.w
        Face3D.__init__(self, surface3d=surface3d,
                        surface2d=surface2d,
                        name=name)

    def copy(self, deep=True, memo=None):
        return CylindricalFace3D(self.surface3d.copy(), self.surface2d.copy(),
                                 self.name)

    def _bounding_box(self):
        theta_min, theta_max, zmin, zmax = self.surface2d.outer_contour.bounding_rectangle()

        lower_center = self.surface3d.frame.origin + zmin * self.surface3d.frame.w
        upper_center = self.surface3d.frame.origin + zmax * self.surface3d.frame.w

        xmin, xmax = volmdlr.geometry.cos_image(theta_min, theta_max)
        ymin, ymax = volmdlr.geometry.sin_image(theta_min, theta_max)

        points = [(lower_center
                   + xmin * self.surface3d.radius * self.surface3d.frame.u
                   + ymin * self.surface3d.radius * self.surface3d.frame.v),
                  (lower_center
                   + xmax * self.surface3d.radius * self.surface3d.frame.u
                   + ymin * self.surface3d.radius * self.surface3d.frame.v),
                  (lower_center
                   + xmin * self.surface3d.radius * self.surface3d.frame.u
                   + ymax * self.surface3d.radius * self.surface3d.frame.v),
                  (lower_center
                   + xmax * self.surface3d.radius * self.surface3d.frame.u
                   + ymax * self.surface3d.radius * self.surface3d.frame.v),
                  (upper_center
                   + xmin * self.surface3d.radius * self.surface3d.frame.u
                   + ymin * self.surface3d.radius * self.surface3d.frame.v),
                  (upper_center
                   + xmax * self.surface3d.radius * self.surface3d.frame.u
                   + ymin * self.surface3d.radius * self.surface3d.frame.v),
                  (upper_center
                   + xmin * self.surface3d.radius * self.surface3d.frame.u
                   + ymax * self.surface3d.radius * self.surface3d.frame.v),
                  (upper_center
                   + xmax * self.surface3d.radius * self.surface3d.frame.u
                   + ymax * self.surface3d.radius * self.surface3d.frame.v)]

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

    def range_closest(self, list_points):
        """
        This method has be edited as it was really bad coded:
            * parameter removed, use of self data instead
        """
        points_set = volmdlr.delete_double_point(list_points)
        points_set3D = CylindricalFace3D.points2d_to3d(None, [points_set],
                                                       self.radius, self.surface3d.frame)

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

        if not self.contours2d[0].point_belongs(pt1_2d):
            # Find the closest one
            points_contours1 = self.contours2d[0].tessel_points

            poly1 = volmdlr.ClosedPolygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d,
                                                       return_other_point=True)
            pt1 = volmdlr.Point3D((r1 * math.cos(new_pt1_2d.vector[0]),
                                   r1 * math.sin(new_pt1_2d.vector[0]),
                                   new_pt1_2d.vector[1]))
            p1 = frame1.old_coordinates(pt1)

        if not other_cyl.contours2d[0].point_belongs(pt2_2d):
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
        # ADD THE FACT THAT PLANEFACE.CONTOURS : [0] = contours totale, le reste = trous
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

        if not self.contours2d[0].point_belongs(pt1_2d):
            # Find the closest one
            points_contours1 = self.contours2d[0].tessel_points

            poly1 = volmdlr.ClosedPolygon2D(points_contours1)
            d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d,
                                                       return_other_point=True)
            pt1 = volmdlr.Point3D((r * math.cos(new_pt1_2d.vector[0]),
                                   r * math.sin(new_pt1_2d.vector[0]),
                                   new_pt1_2d.vector[1]))
            p1 = frame1.old_coordinates(pt1)

        if not planeface.contours[0].point_belongs(pt2_2d):
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

    def adjacent_direction(self, other_face3d):
        '''
        find out in which direction the faces are adjacent
        Parameters
        ----------
        other_face3d : volmdlr.faces.CylindricalFace3D
        Returns
        -------
        adjacent_direction
        '''

        contour1 = self.outer_contour3d
        contour2 = other_face3d.outer_contour3d
        point1, point2 = contour1.shared_primitives_extremities(contour2)

        coord = point1 - point2
        coord = [abs(coord.x), abs(coord.y)]

        if coord.index(max(coord)) == 0:
            return 'x'
        else:
            return 'y'


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

    def __init__(self, surface3d: ToroidalSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):

        # self.toroidalsurface3d = toroidalsurface3d

        self.center = surface3d.frame.origin
        self.normal = surface3d.frame.w

        theta_min, theta_max, phi_min, phi_max = surface2d.outer_contour.bounding_rectangle()

        self.theta_min = theta_min
        self.theta_max = theta_max
        self.phi_min = phi_min
        self.phi_max = phi_max

        # contours3d = [self.toroidalsurface3d.contour2d_to_3d(c)\
        #               for c in [outer_contour2d]+inners_contours2d]

        Face3D.__init__(self,
                        surface3d=surface3d,
                        surface2d=surface2d,
                        name=name)

    def copy(self, deep=True, memo=None):
        return ToroidalFace3D(self.surface3d.copy(), self.surface2d.copy(),
                              self.name)

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

# =============================================================================
#  This code seems buggy...
# =============================================================================

    # def minimum_maximum_tore(self, contour2d):
    #     points = contour2d.tessel_points

    #     min_phi, min_theta = min([pt[1] for pt in points]), min(
    #         [pt[0] for pt in points])
    #     max_phi, max_theta = max([pt[1] for pt in points]), max(
    #         [pt[0] for pt in points])
    #     return min_phi, min_theta, max_phi, max_theta

    # def minimum_distance_points_tore(self, other_tore):
    #     raise NotImplementedError('This method seems unused, its code has been commented')
    #     R1, r1, R2, r2 = self.rcenter, self.rcircle, other_tore.rcenter, other_tore.rcircle

    #     min_phi1, min_theta1, max_phi1, max_theta1 = self.minimum_maximum_tore(
    #         self.contours2d[0])

    #     # start1 = self.start
    #     n1 = self.normal
    #     u1 = self.toroidalsurface3d.frame.u
    #     v1 = self.toroidalsurface3d.frame.v
    #     frame1 = volmdlr.Frame3D(self.center, u1, v1, n1)
    #     # start1 = self.points2d_to3d([[min_theta1, min_phi1]], R1, r1, frame1)

    #     min_phi2, min_theta2, max_phi2, max_theta2 = self.minimum_maximum_tore(
    #         other_tore.contours2d[0])

    #     # start2 = other_tore.start
    #     n2 = other_tore.normal
    #     u2 = other_tore.toroidalsurface3d.frame.u
    #     v2 = other_tore.toroidalsurface3d.frame.v
    #     frame2 = volmdlr.Frame3D(other_tore.center, u2, v2, n2)
    #     # start2 = other_tore.points2d_to3d([[min_theta2, min_phi2]], R2, r2, frame2)

    #     w = other_tore.center - self.center

    #     n1n1, n1u1, n1v1, n1n2, n1u2, n1v2 = n1.dot(n1), n1.dot(u1), n1.dot(
    #         v1), n1.dot(n2), n1.dot(u2), n1.dot(v2)
    #     u1u1, u1v1, u1n2, u1u2, u1v2 = u1.dot(u1), u1.dot(v1), u1.dot(
    #         n2), u1.dot(u2), u1.dot(v2)
    #     v1v1, v1n2, v1u2, v1v2 = v1.dot(v1), v1.dot(n2), v1.dot(u2), v1.dot(v2)
    #     n2n2, n2u2, n2v2 = n2.dot(n2), n2.dot(u2), n2.dot(v2)
    #     u2u2, u2v2, v2v2 = u2.dot(u2), u2.dot(v2), v2.dot(v2)

    #     w2, wn1, wu1, wv1, wn2, wu2, wv2 = w.dot(w), w.dot(n1), w.dot(
    #         u1), w.dot(v1), w.dot(n2), w.dot(u2), w.dot(v2)

    #     # x = (phi1, theta1, phi2, theta2)
    #     def distance_squared(x):
    #         return (u1u1 * (((R1 + r1 * math.cos(x[0])) * math.cos(x[1])) ** 2)
    #                 + v1v1 * (((R1 + r1 * math.cos(x[0])) * math.sin(
    #                     x[1])) ** 2)
    #                 + n1n1 * ((math.sin(x[0])) ** 2) * (r1 ** 2) + w2
    #                 + u2u2 * (((R2 + r2 * math.cos(x[2])) * math.cos(
    #                     x[3])) ** 2)
    #                 + v2v2 * (((R2 + r2 * math.cos(x[2])) * math.sin(
    #                     x[3])) ** 2)
    #                 + n2n2 * ((math.sin(x[2])) ** 2) * (r2 ** 2)
    #                 + 2 * u1v1 * math.cos(x[1]) * math.sin(x[1]) * (
    #                         (R1 + r1 * math.cos(x[0])) ** 2)
    #                 + 2 * (R1 + r1 * math.cos(x[0])) * math.cos(
    #                     x[1]) * r1 * math.sin(x[0]) * n1u1
    #                 - 2 * (R1 + r1 * math.cos(x[0])) * math.cos(x[1]) * wu1
    #                 - 2 * (R1 + r1 * math.cos(x[0])) * (
    #                         R2 + r2 * math.cos(x[2])) * math.cos(
    #                     x[1]) * math.cos(x[3]) * u1u2
    #                 - 2 * (R1 + r1 * math.cos(x[0])) * (
    #                         R2 + r2 * math.cos(x[2])) * math.cos(
    #                     x[1]) * math.sin(x[3]) * u1v2
    #                 - 2 * (R1 + r1 * math.cos(x[0])) * math.cos(
    #                     x[1]) * r2 * math.sin(x[2]) * u1n2
    #                 + 2 * (R1 + r1 * math.cos(x[0])) * math.sin(
    #                     x[1]) * r1 * math.sin(x[0]) * n1v1
    #                 - 2 * (R1 + r1 * math.cos(x[0])) * math.sin(x[1]) * wv1
    #                 - 2 * (R1 + r1 * math.cos(x[0])) * (
    #                         R2 + r2 * math.cos(x[2])) * math.sin(
    #                     x[1]) * math.cos(x[3]) * v1u2
    #                 - 2 * (R1 + r1 * math.cos(x[0])) * (
    #                         R2 + r2 * math.cos(x[2])) * math.sin(
    #                     x[1]) * math.sin(x[3]) * v1v2
    #                 - 2 * (R1 + r1 * math.cos(x[0])) * math.sin(
    #                     x[1]) * r2 * math.sin(x[2]) * v1n2
    #                 - 2 * r1 * math.sin(x[0]) * wn1
    #                 - 2 * r1 * math.sin(x[0]) * (
    #                         R2 + r2 * math.cos(x[2])) * math.cos(
    #                     x[3]) * n1u2
    #                 - 2 * r1 * math.sin(x[0]) * (
    #                         R2 + r2 * math.cos(x[2])) * math.sin(
    #                     x[3]) * n1v2
    #                 - 2 * r1 * r2 * math.sin(x[0]) * math.sin(x[2]) * n1n2
    #                 + 2 * (R2 + r2 * math.cos(x[2])) * math.cos(x[3]) * wu2
    #                 + 2 * (R2 + r2 * math.cos(x[2])) * math.sin(x[3]) * wv2
    #                 + 2 * r2 * math.sin(x[2]) * wn2
    #                 + 2 * u2v2 * math.cos(x[3]) * math.sin(x[3]) * (
    #                         (R2 + r2 * math.cos(x[2])) ** 2)
    #                 + 2 * math.cos(x[3]) * (
    #                         R2 + r2 * math.cos(x[2])) * r2 * math.sin(
    #                     x[2]) * n2u2
    #                 + 2 * math.sin(x[3]) * (
    #                         R2 + r2 * math.cos(x[2])) * r2 * math.sin(
    #                     x[2]) * n2v2)

    #     x01 = npy.array(
    #         [(min_phi1 + max_phi1) / 2, (min_theta1 + max_theta1) / 2,
    #          (min_phi2 + max_phi2) / 2, (min_theta2 + max_theta2) / 2])
    #     x02 = npy.array([min_phi1, min_theta1,
    #                      min_phi2, min_theta2])
    #     x03 = npy.array([max_phi1, max_theta1,
    #                      max_phi2, max_theta2])

    #     minimax = [(min_phi1, min_theta1, min_phi2, min_theta2),
    #                (max_phi1, max_theta1, max_phi2, max_theta2)]

    #     res1 = scp.optimize.least_squares(distance_squared, x01,
    #                                       bounds=minimax)
    #     res2 = scp.optimize.least_squares(distance_squared, x02,
    #                                       bounds=minimax)
    #     res3 = scp.optimize.least_squares(distance_squared, x03,
    #                                       bounds=minimax)

    #     # frame1, frame2 = volmdlr.Frame3D(self.center, u1, v1, n1), volmdlr.Frame3D(other_tore.center, u2, v2, n2)
    #     pt1 = self.points2d_to3d([[res1.x[1], res1.x[0]]], R1, r1, frame1)
    #     pt2 = self.points2d_to3d([[res1.x[3], res1.x[2]]], R2, r2, frame2)
    #     p1, p2 = pt1[0], pt2[0]
    #     d = p1.point_distance(p2)
    #     result = res1

    #     res = [res2, res3]
    #     for couple in res:
    #         ptest1 = self.points2d_to3d([[couple.x[1], couple.x[0]]], R1, r1,
    #                                     frame1)
    #         ptest2 = self.points2d_to3d([[couple.x[3], couple.x[2]]], R2, r2,
    #                                     frame2)
    #         dtest = ptest1[0].point_distance(ptest2[0])
    #         if dtest < d:
    #             result = couple
    #             p1, p2 = ptest1[0], ptest2[0]

    #     pt1_2d, pt2_2d = volmdlr.Point2D(
    #         (result.x[1], result.x[0])), volmdlr.Point2D(
    #         (result.x[3], result.x[2]))

    #     if not self.contours2d[0].point_belongs(pt1_2d):
    #         # Find the closest one
    #         points_contours1 = self.contours2d[0].tessel_points

    #         poly1 = volmdlr.ClosedPolygon2D(points_contours1)
    #         d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d,
    #                                                    return_other_point=True)

    #         pt1 = self.points2d_to3d([new_pt1_2d], R1, r1, frame1)
    #         p1 = pt1[0]

    #     if not other_tore.contours2d[0].point_belongs(pt2_2d):
    #         # Find the closest one
    #         points_contours2 = other_tore.contours2d[0].tessel_points

    #         poly2 = volmdlr.ClosedPolygon2D(points_contours2)
    #         d2, new_pt2_2d = poly2.PointBorderDistance(pt2_2d,
    #                                                    return_other_point=True)

    #         pt2 = self.points2d_to3d([new_pt2_2d], R2, r2, frame2)
    #         p2 = pt2[0]

    #     return p1, p2

    # def minimum_distance_points_cyl(self, cyl):
    #     R2, r2, r = self.rcenter, self.rcircle, cyl.radius

    #     min_h, min_theta, max_h, max_theta = cyl.minimum_maximum(
    #         cyl.contours2d[0], r)

    #     n1 = cyl.normal
    #     u1 = cyl.cylindricalsurface3d.frame.u
    #     v1 = cyl.cylindricalsurface3d.frame.v
    #     frame1 = volmdlr.Frame3D(cyl.center, u1, v1, n1)
    #     # st1 = volmdlr.Point3D((r*math.cos(min_theta), r*math.sin(min_theta), min_h))
    #     # start1 = frame1.old_coordinates(st1)

    #     min_phi2, min_theta2, max_phi2, max_theta2 = self.minimum_maximum_tore(
    #         self.contours2d[0])

    #     n2 = self.normal
    #     u2 = self.toroidalsurface3d.frame.u
    #     v2 = self.toroidalsurface3d.frame.v
    #     frame2 = volmdlr.Frame3D(self.center, u2, v2, n2)
    #     # start2 = self.points2d_to3d([[min_theta2, min_phi2]], R2, r2, frame2)

    #     w = self.center - cyl.center

    #     n1n1, n1u1, n1v1, n1n2, n1u2, n1v2 = n1.dot(n1), n1.dot(u1), n1.dot(
    #         v1), n1.dot(n2), n1.dot(u2), n1.dot(v2)
    #     u1u1, u1v1, u1n2, u1u2, u1v2 = u1.dot(u1), u1.dot(v1), u1.dot(
    #         n2), u1.dot(u2), u1.dot(v2)
    #     v1v1, v1n2, v1u2, v1v2 = v1.dot(v1), v1.dot(n2), v1.dot(u2), v1.dot(v2)
    #     n2n2, n2u2, n2v2 = n2.dot(n2), n2.dot(u2), n2.dot(v2)
    #     u2u2, u2v2, v2v2 = u2.dot(u2), u2.dot(v2), v2.dot(v2)

    #     w2, wn1, wu1, wv1, wn2, wu2, wv2 = w.dot(w), w.dot(n1), w.dot(
    #         u1), w.dot(v1), w.dot(n2), w.dot(u2), w.dot(v2)

    #     # x = (theta, h, phi2, theta2)
    #     def distance_squared(x):
    #         return (u1u1 * ((math.cos(x[0]) * r) ** 2) + v1v1 * (
    #                 (math.sin(x[0]) * r) ** 2)
    #                 + n1n1 * (x[1] ** 2) + w2
    #                 + u2u2 * (((R2 + r2 * math.cos(x[2])) * math.cos(
    #                     x[3])) ** 2)
    #                 + v2v2 * (((R2 + r2 * math.cos(x[2])) * math.sin(
    #                     x[3])) ** 2)
    #                 + n2n2 * ((math.sin(x[2])) ** 2) * (r2 ** 2)
    #                 + 2 * u1v1 * math.cos(x[0]) * math.sin(x[0]) * (r ** 2)
    #                 + 2 * r * math.cos(x[0]) * x[1] * n1u1 - 2 * r * math.cos(
    #                     x[0]) * wu1
    #                 - 2 * r * math.cos(x[0]) * (
    #                         R2 + r2 * math.cos(x[2])) * math.cos(
    #                     x[3]) * u1u2
    #                 - 2 * r * math.cos(x[0]) * (
    #                         R2 + r2 * math.cos(x[2])) * math.sin(
    #                     x[3]) * u1v2
    #                 - 2 * r * math.cos(x[0]) * r2 * math.sin(x[2]) * u1n2
    #                 + 2 * r * math.sin(x[0]) * x[1] * n1v1 - 2 * r * math.sin(
    #                     x[0]) * wv1
    #                 - 2 * r * math.sin(x[0]) * (
    #                         R2 + r2 * math.cos(x[2])) * math.cos(
    #                     x[3]) * v1u2
    #                 - 2 * r * math.sin(x[0]) * (
    #                         R2 + r2 * math.cos(x[2])) * math.sin(
    #                     x[3]) * v1v2
    #                 - 2 * r * math.sin(x[0]) * r2 * math.sin(x[2]) * v1n2 - 2 *
    #                 x[1] * wn1
    #                 - 2 * x[1] * (R2 + r2 * math.cos(x[2])) * math.cos(
    #                     x[3]) * n1u2
    #                 - 2 * x[1] * (R2 + r2 * math.cos(x[2])) * math.sin(
    #                     x[3]) * n1v2
    #                 - 2 * x[1] * r2 * math.sin(x[2]) * n1n2
    #                 + 2 * (R2 + r2 * math.cos(x[2])) * math.cos(x[3]) * wu2
    #                 + 2 * (R2 + r2 * math.cos(x[2])) * math.sin(x[3]) * wv2
    #                 + 2 * r2 * math.sin(x[2]) * wn2
    #                 + 2 * u2v2 * math.cos(x[3]) * math.sin(x[3]) * (
    #                         (R2 + r2 * math.cos(x[2])) ** 2)
    #                 + 2 * math.cos(x[3]) * (
    #                         R2 + r2 * math.cos(x[2])) * r2 * math.sin(
    #                     x[2]) * n2u2
    #                 + 2 * math.sin(x[3]) * (
    #                         R2 + r2 * math.cos(x[2])) * r2 * math.sin(
    #                     x[2]) * n2v2)

    #     x01 = npy.array([(min_theta + max_theta) / 2, (min_h + max_h) / 2,
    #                      (min_phi2 + max_phi2) / 2,
    #                      (min_theta2 + max_theta2) / 2])
    #     x02 = npy.array([min_theta, min_h,
    #                      min_phi2, min_theta2])
    #     x03 = npy.array([max_theta, max_h,
    #                      max_phi2, max_theta2])

    #     minimax = [(min_theta, min_h, min_phi2, min_theta2),
    #                (max_theta, max_h, max_phi2, max_theta2)]

    #     res1 = scp.optimize.least_squares(distance_squared, x01,
    #                                       bounds=minimax)
    #     res2 = scp.optimize.least_squares(distance_squared, x02,
    #                                       bounds=minimax)
    #     res3 = scp.optimize.least_squares(distance_squared, x03,
    #                                       bounds=minimax)

    #     pt1 = volmdlr.Point3D(
    #         (r * math.cos(res1.x[0]), r * math.sin(res1.x[0]), res1.x[1]))
    #     p1 = frame1.old_coordinates(pt1)
    #     pt2 = self.points2d_to3d([[res1.x[3], res1.x[2]]], R2, r2, frame2)
    #     p2 = pt2[0]
    #     d = p1.point_distance(p2)
    #     result = res1

    #     res = [res2, res3]
    #     for couple in res:
    #         pttest1 = volmdlr.Point3D((r * math.cos(couple.x[0]),
    #                                    r * math.sin(couple.x[0]), couple.x[1]))
    #         ptest1 = frame1.old_coordinates(pttest1)
    #         ptest2 = self.points2d_to3d([[couple.x[3], couple.x[2]]], R2, r2,
    #                                     frame2)
    #         dtest = ptest1.point_distance(ptest2[0])
    #         if dtest < d:
    #             result = couple
    #             p1, p2 = ptest1, ptest2[0]

    #     pt1_2d, pt2_2d = volmdlr.Point2D(
    #         (result.x[0], result.x[1])), volmdlr.Point2D(
    #         (result.x[3], result.x[2]))

    #     if not self.contours2d[0].point_belongs(pt2_2d):
    #         # Find the closest one
    #         points_contours2 = self.contours2d[0].tessel_points

    #         poly2 = volmdlr.ClosedPolygon2D(points_contours2)
    #         d2, new_pt2_2d = poly2.PointBorderDistance(pt2_2d,
    #                                                    return_other_point=True)

    #         pt2 = self.points2d_to3d([new_pt2_2d], R2, r2, frame2)
    #         p2 = pt2[0]

    #     if not cyl.contours2d[0].point_belongs(pt1_2d):
    #         # Find the closest one
    #         points_contours1 = cyl.contours2d[0].tessel_points

    #         poly1 = volmdlr.ClosedPolygon2D(points_contours1)
    #         d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d,
    #                                                    return_other_point=True)

    #         pt1 = volmdlr.Point3D((r * math.cos(new_pt1_2d.vector[0]),
    #                                r * math.sin(new_pt1_2d.vector[0]),
    #                                new_pt1_2d.vector[1]))
    #         p1 = frame1.old_coordinates(pt1)

    #     return p1, p2

    # def minimum_distance_points_plane(self,
    #                                   planeface):  # Planeface with contour2D
    #     # TODO: check that it takes into account holes

    #     poly2d = planeface.polygon2D
    #     pfpoints = poly2d.points
    #     xmin, ymin = min([pt[0] for pt in pfpoints]), min(
    #         [pt[1] for pt in pfpoints])
    #     xmax, ymax = max([pt[0] for pt in pfpoints]), max(
    #         [pt[1] for pt in pfpoints])
    #     origin, vx, vy = planeface.plane.origin, planeface.plane.vectors[0], \
    #                      planeface.plane.vectors[1]
    #     pf1_2d, pf2_2d = volmdlr.Point2D((xmin, ymin)), volmdlr.Point2D(
    #         (xmin, ymax))
    #     pf3_2d, pf4_2d = volmdlr.Point2D((xmax, ymin)), volmdlr.Point2D(
    #         (xmax, ymax))
    #     pf1, pf2 = pf1_2d.to_3d(origin, vx, vy), pf2_2d.to_3d(origin, vx, vy)
    #     pf3, _ = pf3_2d.to_3d(origin, vx, vy), pf4_2d.to_3d(origin, vx, vy)

    #     u, v = (pf3 - pf1), (pf2 - pf1)
    #     u.normalize()
    #     v.normalize()

    #     R1, r1 = self.rcenter, self.rcircle
    #     min_phi1, min_theta1, max_phi1, max_theta1 = self.minimum_maximum_tore(
    #         self.contours2d[0])

    #     n1 = self.normal
    #     u1 = self.toroidalsurface3d.frame.u
    #     v1 = self.toroidalsurface3d.frame.v
    #     frame1 = volmdlr.Frame3D(self.center, u1, v1, n1)
    #     # start1 = self.points2d_to3d([[min_theta1, min_phi1]], R1, r1, frame1)

    #     w = self.center - pf1

    #     n1n1, n1u1, n1v1, n1u, n1v = n1.dot(n1), n1.dot(u1), n1.dot(
    #         v1), n1.dot(u), n1.dot(v)
    #     u1u1, u1v1, u1u, u1v = u1.dot(u1), u1.dot(v1), u1.dot(u), u1.dot(v)
    #     v1v1, v1u, v1v = v1.dot(v1), v1.dot(u), v1.dot(v)
    #     uu, uv, vv = u.dot(u), u.dot(v), v.dot(v)

    #     w2, wn1, wu1, wv1, wu, wv = w.dot(w), w.dot(n1), w.dot(u1), w.dot(
    #         v1), w.dot(u), w.dot(v)

    #     # x = (x, y, phi1, theta1)
    #     def distance_squared(x):
    #         return (uu * (x[0] ** 2) + vv * (x[1] ** 2) + w2
    #                 + u1u1 * (((R1 + r1 * math.cos(x[2])) * math.cos(
    #                     x[3])) ** 2)
    #                 + v1v1 * (((R1 + r1 * math.cos(x[2])) * math.sin(
    #                     x[3])) ** 2)
    #                 + n1n1 * ((math.sin(x[2])) ** 2) * (r1 ** 2)
    #                 + 2 * x[0] * x[1] * uv - 2 * x[0] * wu
    #                 - 2 * x[0] * (R1 + r1 * math.cos(x[2])) * math.cos(
    #                     x[3]) * u1u
    #                 - 2 * x[0] * (R1 + r1 * math.cos(x[2])) * math.sin(
    #                     x[3]) * v1u
    #                 - 2 * x[0] * math.sin(x[2]) * r1 * n1u - 2 * x[1] * wv
    #                 - 2 * x[1] * (R1 + r1 * math.cos(x[2])) * math.cos(
    #                     x[3]) * u1v
    #                 - 2 * x[1] * (R1 + r1 * math.cos(x[2])) * math.sin(
    #                     x[3]) * v1v
    #                 - 2 * x[1] * math.sin(x[2]) * r1 * n1v
    #                 + 2 * (R1 + r1 * math.cos(x[2])) * math.cos(x[3]) * wu1
    #                 + 2 * (R1 + r1 * math.cos(x[2])) * math.sin(x[3]) * wv1
    #                 + 2 * math.sin(x[2]) * r1 * wn1
    #                 + 2 * u1v1 * math.cos(x[3]) * math.sin(x[3]) * (
    #                         (R1 + r1 * math.cos(x[2])) ** 2)
    #                 + 2 * (R1 + r1 * math.cos(x[2])) * math.cos(
    #                     x[3]) * r1 * math.sin(x[2]) * n1u1
    #                 + 2 * (R1 + r1 * math.cos(x[2])) * math.sin(
    #                     x[3]) * r1 * math.sin(x[2]) * n1v1)

    #     x01 = npy.array([(xmax - xmin) / 2, (ymax - ymin) / 2,
    #                      (min_phi1 + max_phi1) / 2,
    #                      (min_theta1 + max_theta1) / 2])

    #     minimax = [(0, 0, min_phi1, min_theta1),
    #                (xmax - xmin, ymax - ymin, max_phi1, max_theta1)]

    #     res1 = scp.optimize.least_squares(distance_squared, x01,
    #                                       bounds=minimax)

    #     # frame1 = volmdlr.Frame3D(self.center, u1, v1, n1)
    #     pt1 = self.points2d_to3d([[res1.x[3], res1.x[2]]], R1, r1, frame1)
    #     p1 = pt1[0]
    #     p2 = pf1 + res1.x[2] * u + res1.x[3] * v

    #     pt1_2d = volmdlr.Point2D((res1.x[3], res1.x[2]))
    #     pt2_2d = p2.to_2d(pf1, u, v)

    #     if not self.contours2d[0].point_belongs(pt1_2d):
    #         # Find the closest one
    #         points_contours1 = self.contours2d[0].tessel_points

    #         poly1 = volmdlr.ClosedPolygon2D(points_contours1)
    #         d1, new_pt1_2d = poly1.PointBorderDistance(pt1_2d,
    #                                                    return_other_point=True)

    #         pt1 = self.points2d_to3d([new_pt1_2d], R1, r1, frame1)
    #         p1 = pt1[0]

    #     if not planeface.contours[0].point_belongs(pt2_2d):
    #         # Find the closest one
    #         d2, new_pt2_2d = planeface.polygon2D.PointBorderDistance(pt2_2d,
    #                                                                  return_other_point=True)

    #         p2 = new_pt2_2d.to_3d(pf1, u, v)

    #     return p1, p2

    # def minimum_distance(self, other_face, return_points=False):
    #     if other_face.__class__ is ToroidalFace3D:
    #         p1, p2 = self.minimum_distance_points_tore(other_face)
    #         if return_points:
    #             return p1.point_distance(p2), p1, p2
    #         else:
    #             return p1.point_distance(p2)

    #     if other_face.__class__ is CylindricalFace3D:
    #         p1, p2 = self.minimum_distance_points_cyl(other_face)
    #         if return_points:
    #             return p1.point_distance(p2), p1, p2
    #         else:
    #             return p1.point_distance(p2)

    #     if other_face.__class__ is PlaneFace3D:
    #         p1, p2 = self.minimum_distance_points_plane(other_face)
    #         if return_points:
    #             return p1.point_distance(p2), p1, p2
    #         else:
    #             return p1.point_distance(p2)
    #     else:
    #         return NotImplementedError


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

    def __init__(self, surface3d: ConicalSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):

        Face3D.__init__(self,
                        surface3d=surface3d,
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

    def __init__(self, surface3d: SphericalSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):
        Face3D.__init__(self,
                        surface3d=surface3d,
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
                 surface3d: RuledSurface3D,
                 surface2d: Surface2D,
                 name: str = '',
                 color=None):
        Face3D.__init__(self, surface3d=surface3d,
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
    def __init__(self, surface3d: BSplineSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):
        Face3D.__init__(self,
                        surface3d=surface3d,
                        surface2d=surface2d,
                        name=name)

    def _bounding_box(self):
        return self.surface3d._bounding_box()

    def triangulation_lines(self, resolution=25):
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

    def pair_with(self, other_bspline_face3d):
        '''
        find out how the uv parametric frames are located compared to each other, and also how grid3d can be defined respected to these directions

        Parameters
        ----------
        other_bspline_face3d : volmdlr.faces.BSplineFace3D
        Returns
        -------
        corresponding_direction
        grid2d_direction
        '''

        adjacent_direction1, diff1, adjacent_direction2, diff2 = self.adjacent_direction(other_bspline_face3d)
        corresponding_directions = []
        if (diff1 > 0 and diff2 > 0) or (diff1 < 0 and diff2 < 0):
            corresponding_directions.append(('+' + adjacent_direction1, '+' + adjacent_direction2))
        else:
            corresponding_directions.append(('+' + adjacent_direction1, '-' + adjacent_direction2))

        if adjacent_direction1 == 'u' and adjacent_direction2 == 'u':
            corresponding_directions, grid2d_direction = self.adjacent_direction_uu(
                other_bspline_face3d, corresponding_directions)
        elif adjacent_direction1 == 'v' and adjacent_direction2 == 'v':
            corresponding_directions, grid2d_direction = self.adjacent_direction_vv(
                other_bspline_face3d, corresponding_directions)
        elif adjacent_direction1 == 'u' and adjacent_direction2 == 'v':
            corresponding_directions, grid2d_direction = self.adjacent_direction_uv(
                other_bspline_face3d, corresponding_directions)
        elif adjacent_direction1 == 'v' and adjacent_direction2 == 'u':
            corresponding_directions, grid2d_direction = self.adjacent_direction_vu(
                other_bspline_face3d, corresponding_directions)

        return corresponding_directions, grid2d_direction

    def adjacent_direction_uu(self, other_bspline_face3d, corresponding_directions):

        extremities = self.extremities(other_bspline_face3d)
        start1, start2 = extremities[0], extremities[2]
        borders_points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 0),
                          volmdlr.Point2D(1, 1), volmdlr.Point2D(0, 1)]

        # TODO: compute nearest_point in 'bounding_box points' instead of borders_points
        nearest_start1 = start1.nearest_point(borders_points)
        # nearest_end1 = end1.nearest_point(borders_points)
        nearest_start2 = start2.nearest_point(borders_points)
        # nearest_end2 = end2.nearest_point(borders_points)

        v1 = nearest_start1[1]
        v2 = nearest_start2[1]

        if (v1 == 0 and v2 == 0):
            corresponding_directions.append(('+v', '-v'))
            grid2d_direction = [['+x', '-y'], ['+x', '+y']]

        elif (v1 == 1 and v2 == 1):
            if corresponding_directions == [('+u', '-u')]:
                grid2d_direction = [['+x', '+y'], ['-x', '-y']]
            else:
                grid2d_direction = [['+x', '+y'], ['+x', '-y']]
            corresponding_directions.append(('+v', '-v'))

        elif (v1 == 1 and v2 == 0):
            corresponding_directions.append(('+v', '+v'))
            grid2d_direction = [['+x', '+y'], ['+x', '+y']]

        elif (v1 == 0 and v2 == 1):
            corresponding_directions.append(('+v', '+v'))
            grid2d_direction = [['+x', '-y'], ['+x', '-y']]

        return corresponding_directions, grid2d_direction

    def adjacent_direction_vv(self, other_bspline_face3d, corresponding_directions):

        extremities = self.extremities(other_bspline_face3d)
        start1, start2 = extremities[0], extremities[2]
        borders_points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 0),
                          volmdlr.Point2D(1, 1), volmdlr.Point2D(0, 1)]

        # TODO: compute nearest_point in 'bounding_box points' instead of borders_points
        nearest_start1 = start1.nearest_point(borders_points)
        # nearest_end1 = end1.nearest_point(borders_points)
        nearest_start2 = start2.nearest_point(borders_points)
        # nearest_end2 = end2.nearest_point(borders_points)

        u1 = nearest_start1[0]
        u2 = nearest_start2[0]

        if (u1 == 0 and u2 == 0):
            corresponding_directions.append(('+u', '-v'))
            grid2d_direction = [['-y', '-x'], ['-y', '+x']]

        elif (u1 == 1 and u2 == 1):
            corresponding_directions.append(('+u', '-v'))
            grid2d_direction = [['+y', '+x'], ['+y', '-x']]

        elif (u1 == 0 and u2 == 1):
            corresponding_directions.append(('+u', '+u'))
            grid2d_direction = [['+y', '-x'], ['+y', '-x']]

        elif (u1 == 1 and u2 == 0):
            corresponding_directions.append(('+u', '+u'))
            grid2d_direction = [['+y', '+x'], ['+y', '+x']]

        return corresponding_directions, grid2d_direction

    def adjacent_direction_uv(self, other_bspline_face3d, corresponding_directions):

        extremities = self.extremities(other_bspline_face3d)
        start1, start2 = extremities[0], extremities[2]
        borders_points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 0),
                          volmdlr.Point2D(1, 1), volmdlr.Point2D(0, 1)]

        # TODO: compute nearest_point in 'bounding_box points' instead of borders_points
        nearest_start1 = start1.nearest_point(borders_points)
        # nearest_end1 = end1.nearest_point(borders_points)
        nearest_start2 = start2.nearest_point(borders_points)
        # nearest_end2 = end2.nearest_point(borders_points)

        v1 = nearest_start1[1]
        u2 = nearest_start2[0]

        if (v1 == 1 and u2 == 0):
            corresponding_directions.append(('+v', '+u'))
            grid2d_direction = [['+x', '+y'], ['+y', '+x']]

        elif (v1 == 0 and u2 == 1):
            corresponding_directions.append(('+v', '+u'))
            grid2d_direction = [['-x', '-y'], ['-y', '-x']]

        elif (v1 == 1 and u2 == 1):
            corresponding_directions.append(('+v', '-u'))
            grid2d_direction = [['+x', '+y'], ['-y', '-x']]

        elif (v1 == 0 and u2 == 0):
            corresponding_directions.append(('+v', '-u'))
            grid2d_direction = [['-x', '-y'], ['-y', '+x']]

        return corresponding_directions, grid2d_direction

    def adjacent_direction_vu(self, other_bspline_face3d, corresponding_directions):

        extremities = self.extremities(other_bspline_face3d)
        start1, start2 = extremities[0], extremities[2]
        borders_points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 0),
                          volmdlr.Point2D(1, 1), volmdlr.Point2D(0, 1)]

        # TODO: compute nearest_point in 'bounding_box points' instead of borders_points
        nearest_start1 = start1.nearest_point(borders_points)
        # nearest_end1 = end1.nearest_point(borders_points)
        nearest_start2 = start2.nearest_point(borders_points)
        # nearest_end2 = end2.nearest_point(borders_points)

        u1 = nearest_start1[0]
        v2 = nearest_start2[1]

        if (u1 == 1 and v2 == 0):
            corresponding_directions.append(('+u', '+v'))
            grid2d_direction = [['+y', '+x'], ['+x', '+y']]

        elif (u1 == 0 and v2 == 1):
            corresponding_directions.append(('+u', '+v'))
            grid2d_direction = [['-y', '-x'], ['+x', '-y']]

        elif (u1 == 0 and v2 == 0):
            corresponding_directions.append(('+u', '-v'))
            grid2d_direction = [['+y', '-x'], ['+x', '+y']]

        elif (u1 == 1 and v2 == 1):
            if corresponding_directions == [('+v', '-u')]:
                grid2d_direction = [['+y', '+x'], ['-x', '-y']]
            else:
                grid2d_direction = [['+y', '+x'], ['+x', '-y']]
            corresponding_directions.append(('+u', '-v'))

        return corresponding_directions, grid2d_direction

    def extremities(self, other_bspline_face3d):
        '''
        find points extremities for nearest edges of two faces
        '''
        contour1 = self.outer_contour3d
        contour2 = other_bspline_face3d.outer_contour3d

        contour1_2d = self.surface2d.outer_contour
        contour2_2d = other_bspline_face3d.surface2d.outer_contour

        points1 = [p.start for p in contour1.primitives]
        points2 = [p.start for p in contour2.primitives]

        dis, ind = [], []
        for p in points1:
            pt = p.nearest_point(points2)
            ind.append(points2.index(pt))
            dis.append(p.point_distance(pt))

        dis_sorted = sorted(dis)

        shared = []
        for k, p1 in enumerate(contour1.primitives):
            if dis_sorted[0] == dis_sorted[1]:
                indices = npy.where(npy.array(dis) == dis_sorted[0])[0]
                index1 = indices[0]
                index2 = indices[1]
            else:
                index1 = dis.index(dis_sorted[0])
                index2 = dis.index(dis_sorted[1])
            if ((p1.start == points1[index1] and p1.end == points1[index2])
                or
                    (p1.end == points1[index1] and p1.start == points1[index2])):

                shared.append(p1)
                i = k

        for k, p2 in enumerate(contour2.primitives):
            if ((p2.start == points2[ind[index1]] and p2.end == points2[ind[index2]])
                or
                    (p2.end == points2[ind[index1]] and p2.start == points2[ind[index2]])):

                shared.append(p2)
                j = k

        points = [contour2.primitives[j].start, contour2.primitives[j].end]

        if points.index(contour1.primitives[i].start.nearest_point(points)) == 1:
            start1 = contour1_2d.primitives[i].start
            end1 = contour1_2d.primitives[i].end

            start2 = contour2_2d.primitives[j].end
            end2 = contour2_2d.primitives[j].start

        else:
            start1 = contour1_2d.primitives[i].start
            end1 = contour1_2d.primitives[i].end

            start2 = contour2_2d.primitives[j].start
            end2 = contour2_2d.primitives[j].end

        return start1, end1, start2, end2

    def adjacent_direction(self, other_bspline_face3d):
        '''
        find directions (u or v) between two faces, in the nearest edges between them
        '''

        start1, end1, start2, end2 = self.extremities(other_bspline_face3d)

        du1 = abs((end1 - start1)[0])
        dv1 = abs((end1 - start1)[1])

        if du1 < dv1:
            adjacent_direction1 = 'v'
            diff1 = (end1 - start1)[1]
        else:
            adjacent_direction1 = 'u'
            diff1 = (end1 - start1)[0]

        du2 = abs((end2 - start2)[0])
        dv2 = abs((end2 - start2)[1])

        if du2 < dv2:
            adjacent_direction2 = 'v'
            diff2 = (end2 - start2)[1]
        else:
            adjacent_direction2 = 'u'
            diff2 = (end2 - start2)[0]

        return adjacent_direction1, diff1, adjacent_direction2, diff2

    def adjacent_direction_xy(self, other_face3d):
        '''
        find out in which direction the faces are adjacent
        Parameters
        ----------
        other_face3d : volmdlr.faces.BSplineFace3D
        Returns
        -------
        adjacent_direction
        '''

        contour1 = self.outer_contour3d
        contour2 = other_face3d.outer_contour3d
        point1, point2 = contour1.shared_primitives_extremities(contour2)

        coord = point1 - point2
        coord = [abs(coord.x), abs(coord.y)]

        if coord.index(max(coord)) == 0:
            return 'x'
        else:
            return 'y'

    def merge_with(self, other_bspline_face3d):
        '''
        merge two adjacent faces
        Parameters
        ----------
        other_bspline_face3d : volmdlr.faces.BSplineFace3D
        Returns
        -------
        merged_face : volmdlr.faces.BSplineFace3D
        '''

        merged_surface = self.surface3d.merge_with(other_bspline_face3d.surface3d)
        contours = self.outer_contour3d.merge_with(other_bspline_face3d.outer_contour3d)
        contours.extend(self.inner_contours3d)
        contours.extend(other_bspline_face3d.inner_contours3d)
        merged_face = merged_surface.face_from_contours3d(contours)

        return merged_face


class OpenShell3D(volmdlr.core.CompositePrimitive3D):
    _standalone_in_db = True
    _non_serializable_attributes = ['bounding_box']
    _non_eq_attributes = ['name', 'color', 'alpha', 'bounding_box']
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

    # def __hash__(self):
    #     return sum([hash(f) for f in self.faces])

    # def __eq__(self, other_):
    #     if self.__class__ != other_.__class__:
    #         return False
    #     equal = True
    #     for face, other_face in zip(self.faces, other_.faces):
    #         equal = (equal and face == other_face)
    #     return equal

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

    def copy(self, deep=True, memo=None):
        new_faces = [face.copy() for face in self.faces]
        return self.__class__(new_faces, color=self.color, alpha=self.alpha,
                              name=self.name)

    def union(self, shell2):
        new_faces = self.faces + shell2.faces
        new_name = self.name + ' union ' + shell2.name
        new_color = self.color
        return self.__class__(new_faces, name=new_name, color=new_color)

    def volume(self):
        """
        Does not consider holes
        """
        volume = 0
        for i, face in enumerate(self.faces):
            display3d = face.triangulation()
            points_3D, triangles_indexes = display3d.points, display3d.triangles
            for triangle_index in triangles_indexes:
                point1 = points_3D[triangle_index[0]]
                point2 = points_3D[triangle_index[1]]
                point3 = points_3D[triangle_index[2]]

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
                graph.add_edges_from([(inters[0].primitives[0].start, inters[0].primitives[0].start)])
                intersections.append([inters[0].primitives[0].start, inters[0].primitives[0].start])
        pts = list(nx.dfs_edges(graph, intersections[0][0]))
        # print(pts)
        # print(intersections)
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

    def line_intersections(self,
                           line3d: vme.Line3D) \
            -> List[Tuple[Face3D, List[volmdlr.Point3D]]]:
        intersections = []
        for face in self.faces:
            face_intersections = face.line_intersections(line3d)
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
                if intersection_points:
                    intersection_points = [
                        intersection_points[0].primitives[0].start,
                        intersection_points[0].primitives[0].end]
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
                if intersection_points:
                    intersection_points = [
                        intersection_points[0].primitives[0].start,
                        intersection_points[0].primitives[0].end]
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
        # mesh = vmd.DisplayMesh3D([], [])
        meshes = []
        for i, face in enumerate(self.faces):
            try:
                face_mesh = face.triangulation()
                meshes.append(face_mesh)
                # mesh.merge_mesh(face_mesh)
            except NotImplementedError:
                print('Warning: a face has been skipped in rendering')
        return vmd.DisplayMesh3D.merge_meshes(meshes)

    def babylon_script(self, name='primitive_mesh'):
        s = f'var {name} = new BABYLON.Mesh("{name}", scene);\n'

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

    def plot(self, ax=None, color: str = 'k', alpha: float = 1):
        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')

        for face in self.faces:
            face.plot(ax=ax, color=color, alpha=alpha)

        return ax


class ClosedShell3D(OpenShell3D):
    _standalone_in_db = True
    _non_serializable_attributes = ['bounding_box']
    _non_eq_attributes = ['name', 'color', 'alpha', 'bounding_box']
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

    def copy(self, deep=True, memo=None):
        new_faces = [face.copy() for face in self.faces]
        return self.__class__(new_faces, color=self.color, alpha=self.alpha,
                              name=self.name)

    def face_on_shell(self, face):
        """
        Verifies if a face lies on the shell's surface
        """
        for fc in self.faces:
            if fc.face_inside(face):
                return True
        return False

    def is_face_inside(self, face: Face3D):
        for point in face.outer_contour3d.discretization_points(0.01):
            point_inside_shell = self.point_belongs(point)
            point_in_shells_faces = self.point_in_shell_face(point)
            if (not point_inside_shell) and (not point_in_shells_faces):
                return False
        return True

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

    def point_in_shell_face(self, point: volmdlr.Point3D):

        for face in self.faces:
            point2d = face.surface3d.point3d_to_2d(point)
            if face.point_belongs(point) or \
                    face.surface2d.outer_contour.point_over_contour(
                        point2d, abs_tol=1e-7):
                return True
        return False

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
        for face in self.faces:
            if not shell2.is_face_inside(face):
                return False
        return True

    def is_disjoint_from(self, shell2, tol=1e-8):
        '''
             verifies and rerturns a bool if two shells are disjointed or not.
        '''
        disjoint = True
        if self.bounding_box.bbox_intersection(shell2.bounding_box) or\
                self.bounding_box.distance_to_bbox(shell2.bounding_box) <= tol:
            return False
        return disjoint

    def intersecting_faces_combinations(self, shell2, tol=1e-8):
        '''
            :param shell2: ClosedShell3D
            for two closed shells, it calculates and return a list of face
            combinations (list = [(face_shell1, face_shell2),...])
            for intersecting faces. if two faces can not be intersected,
            there is no combination for those
            :param tol: Corresponde to the tolerance to consider two faces as intersecting faces
        '''
        list_coicident_faces = self.get_coincident_faces(shell2)
        face_combinations = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                if (volmdlr.faces.ClosedShell3D([face1]).bounding_box.bbox_intersection(volmdlr.faces.ClosedShell3D([face2]).bounding_box) or
                        volmdlr.faces.ClosedShell3D(
                            [face1]).bounding_box.distance_to_bbox(
                            volmdlr.faces.ClosedShell3D([face2]).bounding_box) <= tol) and \
                        (face1, face2) not in list_coicident_faces:
                    face_combinations.append((face1, face2))

        return face_combinations

    @staticmethod
    def dict_intersecting_combinations(intersecting_faces_combinations, tol=1e-8):
        '''
            :param intersecting_faces_combinations: list of face combinations (list = [(face_shell1, face_shell2),...]) for intersecting faces.
            :type intersecting_faces_combinations: list of face objects combinaitons

            returns a dictionary containing as keys the combination of intersecting faces
            and as the values the resulting primitive from the two intersecting faces.
            It is done so it is not needed to calculate the same intersecting primitive twice.
        '''
        intersecting_combinations = {}
        for k, combination in enumerate(intersecting_faces_combinations):
            face_intersection = combination[0].face_intersections(combination[1], tol)
            if face_intersection:
                intersecting_combinations[combination] = face_intersection[0]

        return intersecting_combinations

    @staticmethod
    def get_intersecting_faces(dict_intersecting_combinations):
        '''
            :param dict_intersecting_combinations: dictionary containing as keys the combination of intersecting faces
            and as the values the resulting primitive from the two intersecting faces

            returns two lists. One for the intersecting faces in shell1 and the other for the shell2
        '''
        intersecting_faces_shell1 = []
        intersecting_faces_shell2 = []
        for face in list(dict_intersecting_combinations.keys()):
            if face[0] not in intersecting_faces_shell1:
                intersecting_faces_shell1.append(face[0])
            if face[1] not in intersecting_faces_shell2:
                intersecting_faces_shell2.append(face[1])
        return intersecting_faces_shell1, intersecting_faces_shell2

    def get_non_intersecting_faces(self, shell2, intersecting_faces,
                                   intersection_method=False):
        '''
            :param shell2: ClosedShell3D
            :param intersecting_faces:
            returns a list of all the faces that never intersect any
            face of the other shell
        '''
        non_intersecting_faces = []

        for face in self.faces:
            if (face not in intersecting_faces) and (face not in non_intersecting_faces):
                if not intersection_method:
                    if not ClosedShell3D([face]).is_inside_shell(shell2,
                                                                 resolution=0.01):
                        coincident_plane = False
                        for face2 in shell2.faces:
                            if face.surface3d.is_coincident(face2.surface3d) and \
                                    ClosedShell3D([face]).bounding_box.is_inside_bbox(
                                    ClosedShell3D([face2]).bounding_box):
                                coincident_plane = True
                                break
                        if not coincident_plane:
                            non_intersecting_faces.append(face)
                else:
                    if ClosedShell3D([face]).is_inside_shell(shell2, resolution=0.01):
                        non_intersecting_faces.append(face)

        return non_intersecting_faces

    def get_coincident_faces(self, shell2):
        """
        Finds all pairs of faces that are coincidents faces, that is,
        faces lying on the same plane

        returns a List of tuples with the face pairs
        """
        list_coincident_faces = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                if face1.surface3d.is_coincident(face2.surface3d):
                    list_coincident_faces.append((face1, face2))

        return list_coincident_faces

    def two_shells_intersecting_contour(self, shell2, dict_intersecting_combinations=None):
        '''
            :param shell2: ClosedShell3D
            :param dict_intersecting_combinations: dictionary containing as keys the combination of intersecting faces
             and as the values the resulting primitive from the two intersecting faces

            :returns: intersecting contour for two intersecting shells
        '''
        if dict_intersecting_combinations is None:
            face_combinations = self.intersecting_faces_combinations(shell2)
            dict_intersecting_combinations = \
                self.dict_intersecting_combinations(face_combinations)
        intersecting_lines = list(dict_intersecting_combinations.values())
        intersecting_contour = \
            volmdlr.wires.Contour3D([wire.primitives[0] for
                                     wire in intersecting_lines])
        return intersecting_contour

    def reference_shell(self, shell2, face):
        if face in shell2.faces:
            contour_extract_inside = True
            reference_shell = self
        else:
            contour_extract_inside = False
            reference_shell = shell2
        return contour_extract_inside, reference_shell

    def set_operations_valid_exterior_faces(self, new_faces: List[Face3D],
                                            valid_faces: List[Face3D],
                                            shell2, reference_shell):
        for new_face in new_faces:
            inside_reference_shell = reference_shell.point_belongs(
                new_face.random_point_inside())
            if self.set_operations_exterior_face(new_face, valid_faces,
                                                 inside_reference_shell,
                                                 shell2):
                valid_faces.append(new_face)
        return valid_faces

    def union_faces(self, shell2, intersecting_faces,
                    intersecting_combinations):
        faces = []
        # list_coincident_faces = self.get_coincident_faces(shell2)
        for k, face in enumerate(intersecting_faces):
            contour_extract_inside, reference_shell = \
                self.reference_shell(shell2, face)
            new_faces = face.set_operations_new_faces(
                intersecting_combinations, contour_extract_inside)
            faces = self.set_operations_valid_exterior_faces(
                new_faces, faces, shell2, reference_shell)
        return faces

    def get_subtraction_valid_faces(self, new_faces, valid_faces,
                                    reference_shell,
                                    list_coincident_faces,
                                    shell2, keep_interior_faces):
        faces = []
        for new_face in new_faces:
            inside_reference_shell = reference_shell.point_belongs(
                new_face.random_point_inside())
            if keep_interior_faces:
                if self.set_operations_interior_face(new_face, valid_faces,
                                                     inside_reference_shell,
                                                     list_coincident_faces):
                    faces.append(new_face)
            elif self.set_operations_exterior_face(new_face, faces,
                                                   inside_reference_shell,
                                                   shell2):
                faces.append(new_face)
        return faces

    def subtraction_faces(self, shell2, intersecting_faces,
                          intersecting_combinations):
        faces = []
        list_coincident_faces = self.get_coincident_faces(shell2)
        for k, face in enumerate(intersecting_faces):
            keep_interior_faces = False
            if face in shell2.faces:
                keep_interior_faces = True
            contour_extract_inside, reference_shell = \
                self.reference_shell(shell2, face)
            new_faces = face.set_operations_new_faces(
                intersecting_combinations, contour_extract_inside)
            faces.extend(self.get_subtraction_valid_faces(
                new_faces, faces, reference_shell,
                list_coincident_faces, shell2, keep_interior_faces))

        return faces

    def valid_intersection_faces(self, new_faces, valid_faces,
                                 reference_shell, list_coincident_faces):
        faces = []
        for new_face in new_faces:
            inside_reference_shell = reference_shell.point_belongs(
                new_face.random_point_inside())
            if self.set_operations_interior_face(new_face, valid_faces,
                                                 inside_reference_shell,
                                                 list_coincident_faces):
                faces.append(new_face)

        return faces

    def intersection_faces(self, shell2, intersecting_faces,
                           intersecting_combinations):
        faces = []
        list_coincident_faces = self.get_coincident_faces(shell2)
        for face in intersecting_faces:
            contour_extract_inside, reference_shell = \
                self.reference_shell(shell2, face)
            new_faces = face.set_operations_new_faces(
                intersecting_combinations, contour_extract_inside)
            faces.extend(self.valid_intersection_faces(
                new_faces, faces, reference_shell, list_coincident_faces))

        valid_faces = []
        for i, fc1 in enumerate(faces):
            valid_face = True
            for j, fc2 in enumerate(faces):
                if i != j:
                    if fc2.face_inside(fc1):
                        valid_face = False
            if valid_face and fc1 not in valid_faces:
                valid_faces.append(fc1)
        return valid_faces

    @staticmethod
    def set_operations_interior_face(new_face, faces, inside_reference_shell,
                                     list_coincident_faces):
        if inside_reference_shell:
            if new_face not in faces:
                return True
        for coin_f1, coin_f2 in list_coincident_faces:
            if coin_f1.face_inside(new_face) and coin_f2.face_inside(
                    new_face):
                valid = True
                for fc in faces:
                    if fc.face_inside(new_face) or new_face.face_inside(
                            fc):
                        valid = False
                if valid:
                    return True
        return False

    def is_face_between_shells(self, shell2, face):
        centroide = face.surface2d.outer_contour.center_of_mass()
        normal1 = face.surface3d.point2d_to_3d(
            centroide) - 0.001 * face.surface3d.frame.w
        normal2 = face.surface3d.point2d_to_3d(
            centroide) + 0.001 * face.surface3d.frame.w
        if (self.point_belongs(normal1) and
            shell2.point_belongs(normal2)) or \
                (shell2.point_belongs(normal1) and
                 self.point_belongs(normal2)):
            return True
        return False

    def set_operations_exterior_face(self, new_face, valid_faces,
                                     inside_reference_shell, shell2):
        if new_face.area() < 1e-8:
            return False
        if new_face not in valid_faces and not inside_reference_shell:
            for fc in valid_faces:
                if self.is_face_between_shells(shell2, new_face) or\
                        (fc.face_inside(new_face) and
                         new_face.area() == fc.area()):
                    return False
            return True
        return False

    def validate_set_operation(self, shell2, tol):
        '''
        Verifies if two shells are valid for union or subtractions operations,
        that is, if they are disjointed or if one is totaly inside the other
        If it returns an empty list, it means the two shells are valid to continue the
        operation.
        '''
        if self.is_disjoint_from(shell2, tol):
            return [self, shell2]
        if self.is_inside_shell(shell2, resolution=0.01):
            return [shell2]
        if shell2.is_inside_shell(self, resolution=0.01):
            return [self]
        return []

    def union(self, shell2, tol=1e-8):
        '''
            Given Two closed shells, it returns
            a new united ClosedShell3D object
        '''

        validate_set_operation = \
            self.validate_set_operation(shell2, tol)
        if validate_set_operation:
            return validate_set_operation

        face_combinations = self.intersecting_faces_combinations(shell2, tol)

        intersecting_combinations = \
            self.dict_intersecting_combinations(face_combinations, tol)

        intersecting_faces1, intersecting_faces2 = \
            self.get_intersecting_faces(intersecting_combinations)
        intersecting_faces = intersecting_faces1 + intersecting_faces2
        faces = self.get_non_intersecting_faces(
            shell2, intersecting_faces) + shell2.get_non_intersecting_faces(
            self, intersecting_faces)
        if len(faces) == \
                len(self.faces + shell2.faces) and not intersecting_faces:
            return [self, shell2]
        new_valid_faces = self.union_faces(shell2, intersecting_faces,
                                           intersecting_combinations)
        faces += new_valid_faces
        new_shell = ClosedShell3D(faces)
        new_shell.merge_union_faces()
        return [new_shell]

    @staticmethod
    def get_faces_to_be_merged(union_faces):
        coincident_planes_faces = []
        valid_coicident_faces = []
        for i, face1 in enumerate(union_faces):
            for j, face2 in enumerate(union_faces):
                if j != i and \
                        face1.surface3d.is_coincident(face2.surface3d):
                    if face1 not in coincident_planes_faces:
                        coincident_planes_faces.append(face1)
                    coincident_planes_faces.append(face2)
            if coincident_planes_faces:
                for f1, f2 in \
                        product(coincident_planes_faces, repeat=2):
                    if f1 != f2 and f1.is_adjacent(f2):
                        if f1 not in valid_coicident_faces:
                            valid_coicident_faces.append(f1)
                        if f2 not in valid_coicident_faces:
                            valid_coicident_faces.append(f2)
                break
        return valid_coicident_faces

    def merge_union_faces(self):
        union_faces = self.faces
        finished = False
        list_new_faces = []
        count = 0
        while not finished:
            valid_coicident_faces = ClosedShell3D.get_faces_to_be_merged(
                union_faces)
            list_valid_coincident_faces = valid_coicident_faces[:]

            if valid_coicident_faces:
                list_new_faces += PlaneFace3D.merge_faces(valid_coicident_faces)
            if list_valid_coincident_faces:
                for face in list_valid_coincident_faces:
                    union_faces.remove(face)
            count += 1
            if count >= len(self.faces) and not list_valid_coincident_faces:
                finished = True

        list_new_faces += union_faces
        self.faces = list_new_faces

    def subtract(self, shell2, tol=1e-8):
        '''
            Given Two closed shells, it returns a new subtracted OpenShell3D object
        '''
        validate_set_operation = self.validate_set_operation(shell2, tol)
        if validate_set_operation:
            return validate_set_operation

        face_combinations = self.intersecting_faces_combinations(shell2, tol)

        intersecting_combinations = self.dict_intersecting_combinations(face_combinations, tol)

        if len(intersecting_combinations) == 0:
            return [self, shell2]

        intersecting_faces, _ = self.get_intersecting_faces(intersecting_combinations)

        faces = self.get_non_intersecting_faces(shell2, intersecting_faces)
        new_valid_faces = self.union_faces(shell2, intersecting_faces,
                                           intersecting_combinations)
        faces += new_valid_faces

        return [OpenShell3D(faces)]

    def subtract_to_closed_shell(self, shell2, tol=1e-8):
        '''
            Given Two closed shells, it returns a new subtracted ClosedShell3D object
        '''
        validate_set_operation = self.validate_set_operation(shell2, tol)
        if validate_set_operation:
            return validate_set_operation

        face_combinations = self.intersecting_faces_combinations(shell2, tol)

        intersecting_combinations = self.dict_intersecting_combinations(face_combinations, tol)

        if len(intersecting_combinations) == 0:
            return [self, shell2]

        intersecting_faces1, intersecting_faces2 = self.get_intersecting_faces(
            intersecting_combinations)
        intersecting_faces = intersecting_faces1 + intersecting_faces2

        faces = self.get_non_intersecting_faces(shell2, intersecting_faces)
        faces += shell2.get_non_intersecting_faces(self, intersecting_faces,
                                                   intersection_method=True)
        new_valid_faces = self.subtraction_faces(shell2, intersecting_faces,
                                                 intersecting_combinations)
        faces += new_valid_faces

        return [ClosedShell3D(faces)]

    def intersection(self, shell2, tol=1e-8):
        """
        Given two ClosedShell3D, it returns the new objet resulting
         from the intersection of the two
        """
        validate_set_operation = self.validate_set_operation(
            shell2, tol)
        if validate_set_operation:
            return validate_set_operation

        face_combinations = self.intersecting_faces_combinations(shell2, tol)

        intersecting_combinations = self.dict_intersecting_combinations(
            face_combinations, tol)

        if len(intersecting_combinations) == 0:
            return [self, shell2]

        intersecting_faces1, intersecting_faces2 = self.get_intersecting_faces(
            intersecting_combinations)
        intersecting_faces = intersecting_faces1 + intersecting_faces2
        faces = self.intersection_faces(shell2, intersecting_faces,
                                        intersecting_combinations)
        faces += self.get_non_intersecting_faces(shell2,
                                                 intersecting_faces,
                                                 intersection_method=True) + shell2.get_non_intersecting_faces(self,
                                                                                                               intersecting_faces,
                                                                                                               intersection_method=True)

        return [ClosedShell3D(faces)]
