"""
Surfaces & faces.

"""

import math
import warnings
from itertools import chain, product
from typing import List, Tuple, Union

import matplotlib.pyplot as plt
import networkx as nx
import numpy as npy
from scipy.optimize import least_squares, minimize
import triangle as triangle_lib

from geomdl import NURBS, BSpline, utilities
from geomdl.construct import extract_curves
from geomdl.fitting import approximate_surface, interpolate_surface
from geomdl.operations import split_surface_u, split_surface_v
from trimesh import Trimesh

from dessia_common.core import DessiaObject

import volmdlr.bspline_compiled
import volmdlr.core
from volmdlr.core import EdgeStyle
import volmdlr.core_compiled
import volmdlr.display as vmd
import volmdlr.edges as vme
import volmdlr.geometry
import volmdlr.grid
import volmdlr.utils.parametric as vm_parametric
import volmdlr.wires
from volmdlr.utils.parametric import array_range_search, repair_start_end_angle_periodicity, angle_discontinuity
from volmdlr.bspline_evaluators import evaluate_single
from volmdlr.core import point_in_list


def knots_vector_inv(knots_vector):
    """
    Compute knot-elements and multiplicities based on the global knot vector.

    """

    knots = sorted(set(knots_vector))
    multiplicities = [knots_vector.count(knot) for knot in knots]

    return knots, multiplicities


class Surface2D(volmdlr.core.Primitive2D):
    """
    A surface bounded by an outer contour.

    """

    def __init__(self, outer_contour: volmdlr.wires.Contour2D,
                 inner_contours: List[volmdlr.wires.Contour2D],
                 name: str = 'name'):
        self.outer_contour = outer_contour
        self.inner_contours = inner_contours
        self._area = None

        volmdlr.core.Primitive2D.__init__(self, name=name)

    def __hash__(self):
        return hash((self.outer_contour, tuple(self.inner_contours)))

    def _data_hash(self):
        return hash(self)

    def copy(self):
        """
        Copies the surface2d.

        """
        return self.__class__(outer_contour=self.outer_contour.copy(),
                              inner_contours=[c.copy() for c in self.inner_contours],
                              name=self.name)

    def area(self):
        """
        Computes the area of the surface.

        """
        if not self._area:
            self._area = self.outer_contour.area() - sum(contour.area() for contour in self.inner_contours)
        return self._area

    def second_moment_area(self, point: volmdlr.Point2D):
        """
        Computes the second moment area of the surface.

        """
        i_x, i_y, i_xy = self.outer_contour.second_moment_area(point)
        for contour in self.inner_contours:
            i_xc, i_yc, i_xyc = contour.second_moment_area(point)
            i_x -= i_xc
            i_y -= i_yc
            i_xy -= i_xyc
        return i_x, i_y, i_xy

    def center_of_mass(self):
        """
        Compute the center of mass of the 2D surface.

        :return: The center of mass of the surface.
        :rtype: :class:`volmdlr.Point2D`
        """
        center = self.outer_contour.area() * self.outer_contour.center_of_mass()
        for contour in self.inner_contours:
            center -= contour.area() * contour.center_of_mass()
        return center / self.area()

    def point_belongs(self, point2d: volmdlr.Point2D):
        """
        Check whether a point belongs to the 2D surface.

        :param point2d: The point to check.
        :type point2d: :class:`volmdlr.Point2D`
        :return: True if the point belongs to the surface, False otherwise.
        :rtype: bool
        """
        if not self.outer_contour.point_belongs(point2d):
            if self.outer_contour.point_over_contour(point2d):
                return True
            return False

        for inner_contour in self.inner_contours:
            if inner_contour.point_belongs(point2d) and not inner_contour.point_over_contour(point2d):
                return False
        return True

    def random_point_inside(self):
        """
        Generate a random point inside the 2D surface.

        Taking into account any inner contours (holes) it may have.

        :return: A random point inside the surface.
        :rtype: :class:`volmdlr.Point2D`
        """
        valid_point = False
        point_inside_outer_contour = None
        while not valid_point:
            point_inside_outer_contour = self.outer_contour.random_point_inside()
            inside_inner_contour = False
            for inner_contour in self.inner_contours:
                if inner_contour.point_belongs(point_inside_outer_contour):
                    inside_inner_contour = True
            if not inside_inner_contour and \
                    point_inside_outer_contour is not None:
                valid_point = True

        return point_inside_outer_contour

    @staticmethod
    def triangulation_without_holes(vertices, segments, points_grid, tri_opt):
        """
        Triangulates a surface without holes.

        :param vertices: vertices of the surface.
        :param segments: segments defined as tuples of vertices.
        :param points_grid: to do.
        :param tri_opt: triangulation option: "p"
        :return:
        """
        vertices_grid = [(p.x, p.y) for p in points_grid]
        vertices.extend(vertices_grid)
        tri = {'vertices': npy.array(vertices).reshape((-1, 2)),
               'segments': npy.array(segments).reshape((-1, 2)),
               }
        t = triangle_lib.triangulate(tri, tri_opt)
        triangles = t['triangles'].tolist()
        np = t['vertices'].shape[0]
        points = [vmd.Node2D(*t['vertices'][i, :]) for i in range(np)]
        return vmd.DisplayMesh2D(points, triangles=triangles)

    def triangulation(self, number_points_x: int = 15, number_points_y: int = 15):
        """
        Triangulates the Surface2D using the Triangle library.

        :param number_points_x: Number of discretization points in x direction.
        :type number_points_x: int
        :param number_points_y: Number of discretization points in y direction.
        :type number_points_y: int
        :return: The triangulated surface as a display mesh.
        :rtype: :class:`volmdlr.display.DisplayMesh2D`
        """
        area = self.bounding_rectangle().area()
        tri_opt = "p"
        if math.isclose(area, 0., abs_tol=1e-6):
            return vmd.DisplayMesh2D([], triangles=[])

        triangulates_with_grid = number_points_x > 0 or number_points_y > 0

        outer_polygon = self.outer_contour.to_polygon(angle_resolution=15, discretize_line=triangulates_with_grid)

        if not self.inner_contours and not triangulates_with_grid:
            return outer_polygon.triangulation()

        points_grid, x, y, grid_point_index = outer_polygon.grid_triangulation_points(number_points_x=number_points_x,
                                                                                      number_points_y=number_points_y)
        points = [vmd.Node2D(*point) for point in outer_polygon.points]
        vertices = [(point.x, point.y) for point in points]
        n = len(points)
        segments = [(i, i + 1) for i in range(n - 1)]
        segments.append((n - 1, 0))

        if not self.inner_contours:  # No holes
            return self.triangulation_without_holes(vertices, segments, points_grid, tri_opt)

        point_index = {p: i for i, p in enumerate(points)}
        holes = []
        for inner_contour in self.inner_contours:
            inner_polygon = inner_contour.to_polygon(angle_resolution=10, discretize_line=triangulates_with_grid)
            inner_polygon_nodes = [vmd.Node2D.from_point(p) for p in inner_polygon.points]
            for point in inner_polygon_nodes:
                if point not in point_index:
                    points.append(point)
                    vertices.append((point.x, point.y))
                    point_index[point] = n
                    n += 1

            for point1, point2 in zip(inner_polygon_nodes[:-1],
                                      inner_polygon_nodes[1:]):
                segments.append((point_index[point1], point_index[point2]))
            segments.append((point_index[inner_polygon_nodes[-1]], point_index[inner_polygon_nodes[0]]))
            rpi = inner_polygon.barycenter()
            if not inner_polygon.point_belongs(rpi, include_edge_points=False):
                rpi = inner_polygon.random_point_inside(include_edge_points=False)
            holes.append([rpi.x, rpi.y])

            if triangulates_with_grid:
                # removes with a region search the grid points that are in the inner contour
                xmin, xmax, ymin, ymax = inner_polygon.bounding_rectangle.bounds()
                x_grid_range = array_range_search(x, xmin, xmax)
                y_grid_range = array_range_search(y, ymin, ymax)
                for i in x_grid_range:
                    for j in y_grid_range:
                        point = grid_point_index.get((i, j))
                        if not point:
                            continue
                        if inner_polygon.point_belongs(point):
                            points_grid.remove(point)
                            grid_point_index.pop((i, j))

        if triangulates_with_grid:
            vertices_grid = [(p.x, p.y) for p in points_grid]
            vertices.extend(vertices_grid)

        tri = {'vertices': npy.array(vertices).reshape((-1, 2)),
               'segments': npy.array(segments).reshape((-1, 2)),
               'holes': npy.array(holes).reshape((-1, 2))
               }
        t = triangle_lib.triangulate(tri, tri_opt)
        triangles = t['triangles'].tolist()
        np = t['vertices'].shape[0]
        points = [vmd.Node2D(*t['vertices'][i, :]) for i in range(np)]
        return vmd.DisplayMesh2D(points, triangles=triangles)

    def split_by_lines(self, lines):
        """
        Returns a list of cut surfaces given by the lines provided as argument.
        """
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
        Split in n slices.
        """
        bounding_rectangle = self.outer_contour.bounding_rectangle
        lines = []
        for i in range(n - 1):
            xi = bounding_rectangle[0] + (i + 1) * (bounding_rectangle[1] - bounding_rectangle[0]) / n
            lines.append(vme.Line2D(volmdlr.Point2D(xi, 0),
                                    volmdlr.Point2D(xi, 1)))
        return self.split_by_lines(lines)

    def cut_by_line(self, line: vme.Line2D):
        """
        Returns a list of cut Surface2D by the given line.

        :param line: The line to cut the Surface2D with.
        :type line: :class:`volmdlr.edges.Line2D`
        :return: A list of 2D surfaces resulting from the cut.
        :rtype: List[:class:`volmdlr.faces.Surface2D`]
        """
        surfaces = []
        splitted_outer_contours = self.outer_contour.cut_by_line(line)
        splitted_inner_contours_table = []
        for inner_contour in self.inner_contours:
            splitted_inner_contours = inner_contour.cut_by_line(line)
            splitted_inner_contours_table.append(splitted_inner_contours)

        # First part of the external contour
        for outer_split in splitted_outer_contours:
            inner_contours = []
            for splitted_inner_contours in splitted_inner_contours_table:
                for inner_split in splitted_inner_contours:
                    inner_split.order_contour()
                    point = inner_split.random_point_inside()
                    if outer_split.point_belongs(point):
                        inner_contours.append(inner_split)

            if inner_contours:
                surface2d = self.from_contours(outer_split, inner_contours)
                surfaces.append(surface2d)
            else:
                surfaces.append(Surface2D(outer_split, []))
        return surfaces

    def line_crossings(self, line: vme.Line2D):
        """
        Find intersection points between a line and the 2D surface.

        :param line: The line to intersect with the shape.
        :type line: :class:`volmdlr.edges.Line2D`
        :return: A list of intersection points sorted by increasing abscissa
            along the line. Each intersection point is a tuple
            (point, primitive) where point is the intersection point and
            primitive is the intersected primitive.
        :rtype: List[Tuple[:class:`volmdlr.Point2D`,
            :class:`volmdlr.core.Primitive2D`]]

        """
        intersection_points = []
        for primitive in self.outer_contour.primitives:
            for p in primitive.line_crossings(line):
                if (p, primitive) not in intersection_points:
                    intersection_points.append((p, primitive))
        for inner_contour in self.inner_contours:
            for primitive in inner_contour.primitives:
                for p in primitive.line_crossings(line):
                    if (p, primitive) not in intersection_points:
                        intersection_points.append((p, primitive))
        return sorted(intersection_points, key=lambda ip: line.abscissa(ip[0]))

    def split_at_centers(self):
        """
        Split in n slices.

        # TODO: is this used ?
        """

        cutted_contours = []
        center_of_mass1 = self.inner_contours[0].center_of_mass()
        center_of_mass2 = self.inner_contours[1].center_of_mass()
        cut_line = vme.Line2D(center_of_mass1, center_of_mass2)

        iteration_contours2 = []

        surface_cut = self.cut_by_line(cut_line)

        iteration_contours2.extend(surface_cut)

        iteration_contours = iteration_contours2[:]
        cutted_contours.extend(iteration_contours)

        return cutted_contours

    def cut_by_line2(self, line):
        """
        Cuts a Surface2D with line (2).

        :param line: DESCRIPTION
        :type line: TYPE
        :raises NotImplementedError: DESCRIPTION
        :return: DESCRIPTION
        :rtype: TYPE

        """

        all_contours = []
        inner_1 = self.inner_contours[0]
        inner_2 = self.inner_contours[1]

        inner_intersections_1 = inner_1.line_intersections(line)
        inner_intersections_2 = inner_2.line_intersections(line)

        arc1, arc2 = inner_1.split(inner_intersections_1[1],
                                   inner_intersections_1[0])
        arc3, arc4 = inner_2.split(inner_intersections_2[1],
                                   inner_intersections_2[0])
        new_inner_1 = volmdlr.wires.Contour2D([arc1, arc2])
        new_inner_2 = volmdlr.wires.Contour2D([arc3, arc4])

        intersections = [(inner_intersections_1[0], arc1), (inner_intersections_1[1], arc2)]
        intersections += self.outer_contour.line_intersections(line)
        intersections.append((inner_intersections_2[0], arc3))
        intersections.append((inner_intersections_2[1], arc4))
        intersections += self.outer_contour.line_intersections(line)

        if not intersections:
            all_contours.extend([self])
        if len(intersections) < 4:
            return [self]
        if len(intersections) >= 4:
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

                primitives1 = [volmdlr.edges.LineSegment2D(intersections[6][0], intersections[1][0]),
                               new_inner_1.primitives[ip1],
                               volmdlr.edges.LineSegment2D(intersections[0][0], intersections[5][0]),
                               new_inner_2.primitives[ip5],
                               volmdlr.edges.LineSegment2D(intersections[4][0], intersections[7][0]),
                               sp44
                               ]
                primitives1.extend(self.outer_contour.primitives[ip3 + 1:ip4])
                primitives1.append(sp34)

                primitives2 = [volmdlr.edges.LineSegment2D(intersections[7][0], intersections[4][0]),
                               new_inner_2.primitives[ip6],
                               volmdlr.edges.LineSegment2D(intersections[5][0], intersections[0][0]),
                               new_inner_1.primitives[ip2],
                               volmdlr.edges.LineSegment2D(intersections[1][0], intersections[6][0]),
                               sp33
                               ]

                primitives2.extend(self.outer_contour.primitives[:ip3].reverse())
                primitives2.append(sp43)

                all_contours.extend([volmdlr.wires.Contour2D(primitives1),
                                     volmdlr.wires.Contour2D(primitives2)])

            else:
                raise NotImplementedError(
                    'Non convex contour not supported yet')

        return all_contours

    # def cut_by_line3(self, line):
    #     """
    #     Cuts a Surface2D with line (2).
    #
    #     :param line: DESCRIPTION
    #     :type line: TYPE
    #     :raises NotImplementedError: DESCRIPTION
    #     :return: DESCRIPTION
    #     :rtype: TYPE
    #
    #     """
    #
    #     # ax=self.outer_contour.plot()
    #     all_contours = []
    #     inner = self.inner_contours[0]
    #     inner_2 = self.inner_contours[1]
    #     inner_3 = self.inner_contours[2]
    #
    #     c = inner.center_of_mass()
    #     c_2 = inner_2.center_of_mass()
    #     c_3 = inner_3.center_of_mass()
    #     direction_vector = line.normal_vector()
    #     direction_line = vme.Line2D(c, volmdlr.Point2D(
    #         (direction_vector.y * c.x - direction_vector.x * c.y) / (
    #             direction_vector.y), 0))
    #     direction_line_2 = vme.Line2D(c_2, volmdlr.Point2D(
    #         (direction_vector.y * c_2.x - direction_vector.x * c_2.y) / (
    #             direction_vector.y), 0))
    #
    #     direction_line_3 = vme.Line2D(c_3, volmdlr.Point2D(
    #         (direction_vector.y * c_3.x - direction_vector.x * c_3.y) / (
    #             direction_vector.y), 0))
    #     inner_intersections = inner.line_intersections(direction_line)
    #     inner_intersections_2 = inner_2.line_intersections(direction_line_2)
    #     inner_intersections_3 = inner_3.line_intersections(direction_line_3)
    #     arc1, arc2 = inner.split(inner_intersections[1],
    #                              inner_intersections[0])
    #     arc3, arc4 = inner_2.split(inner_intersections_2[1],
    #                                inner_intersections_2[0])
    #     arc5, arc6 = inner_3.split(inner_intersections_3[1],
    #                                inner_intersections_3[0])
    #     new_inner = volmdlr.wires.Contour2D([arc1, arc2])
    #     new_inner_2 = volmdlr.wires.Contour2D([arc3, arc4])
    #     new_inner_3 = volmdlr.wires.Contour2D([arc5, arc6])
    #     intersections = [(inner_intersections[0], arc1), (inner_intersections[1], arc2)]
    #
    #     if len(self.outer_contour.line_intersections(direction_line)) > 2:
    #
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line)[0])
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line)[2])
    #     else:
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line)[0])
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line)[1])
    #     intersections.append((inner_intersections_2[0], arc3))
    #     intersections.append((inner_intersections_2[1], arc4))
    #     if len(self.outer_contour.line_intersections(direction_line_2)) > 2:
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line_2)[0])
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line_2)[2])
    #     else:
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line_2)[0])
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line_2)[1])
    #     intersections.append((inner_intersections_3[0], arc5))
    #     intersections.append((inner_intersections_3[1], arc6))
    #     if len(self.outer_contour.line_intersections(direction_line_3)) > 2:
    #
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line_3)[0])
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line_3)[2])
    #     else:
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line_3)[0])
    #         intersections.append(
    #             self.outer_contour.line_intersections(direction_line_3)[1])
    #
    #     if isinstance(intersections[0][0], volmdlr.Point2D) and \
    #             isinstance(intersections[1][0], volmdlr.Point2D):
    #         ip1, ip2 = sorted([new_inner.primitives.index(intersections[0][1]),
    #                            new_inner.primitives.index(
    #                                intersections[1][1])])
    #         ip5, ip6 = sorted(
    #             [new_inner_2.primitives.index(intersections[4][1]),
    #              new_inner_2.primitives.index(intersections[5][1])])
    #         ip7, ip8 = sorted(
    #             [new_inner_3.primitives.index(intersections[8][1]),
    #              new_inner_3.primitives.index(intersections[9][1])])
    #         ip3, ip4 = sorted(
    #             [self.outer_contour.primitives.index(intersections[2][1]),
    #              self.outer_contour.primitives.index(intersections[3][1])])
    #
    #         sp11, sp12 = intersections[2][1].split(intersections[2][0])
    #         sp21, sp22 = intersections[3][1].split(intersections[3][0])
    #         sp33, sp34 = intersections[6][1].split(intersections[6][0])
    #         sp44, sp43 = intersections[7][1].split(intersections[7][0])
    #         sp55, sp56 = intersections[10][1].split(intersections[10][0])
    #         sp66, sp65 = intersections[11][1].split(intersections[11][0])
    #
    #         primitives1 = []
    #         primitives1.append(volmdlr.edges.LineSegment2D(intersections[7][0],
    #                                                        intersections[5][0]))
    #         primitives1.append(new_inner_2.primitives[ip5])
    #         primitives1.append(volmdlr.edges.LineSegment2D(intersections[6][0],
    #                                                        intersections[4][0]))
    #         primitives1.append(sp33)
    #         primitives1.append(sp43)
    #
    #         primitives2 = []
    #         primitives2.append(volmdlr.edges.LineSegment2D(intersections[6][0],
    #                                                        intersections[4][0]))
    #         primitives2.append(new_inner_2.primitives[ip6])
    #         primitives2.append(volmdlr.edges.LineSegment2D(intersections[5][0],
    #                                                        intersections[7][0]))
    #         primitives2.append(volmdlr.edges.LineSegment2D(intersections[7][0],
    #                                                        intersections[11][0]))
    #         primitives2.append(
    #             vme.LineSegment2D(intersections[11][0],
    #                                         intersections[9][0]))
    #         primitives2.append(new_inner_3.primitives[ip7])
    #         primitives2.append(volmdlr.edges.LineSegment2D(intersections[8][0],
    #                                                        intersections[10][0]))
    #         primitives2.append(sp34)
    #
    #         primitives3 = []
    #         primitives3.append(
    #             vme.LineSegment2D(intersections[10][0],
    #                                         intersections[8][0]))
    #         primitives3.append(new_inner_3.primitives[ip8])
    #         primitives3.append(volmdlr.edges.LineSegment2D(intersections[9][0],
    #                                                        intersections[11][0]))
    #         primitives3.append(sp22)
    #         primitives3.append(volmdlr.edges.LineSegment2D(intersections[3][0],
    #                                                        intersections[1][0]))
    #         primitives3.append(new_inner.primitives[ip1])
    #         primitives3.append(vme.LineSegment2D(intersections[0][0],
    #                                                        intersections[2][
    #                                                            0]))
    #         primitives3.append(vme.LineSegment2D(intersections[2][0],
    #                                                        intersections[10][
    #                                                            0]))
    #
    #         primitives4 = [volmdlr.edges.LineSegment2D(intersections[3][0],
    #                                                    intersections[1][0])]
    #         a = volmdlr.edges.Arc2D(new_inner.primitives[ip2].end,
    #                                 new_inner.primitives[ip2].interior,
    #                                 new_inner.primitives[ip2].start)
    #         primitives4.append(a)
    #         primitives4.append(volmdlr.edges.LineSegment2D(intersections[0][0],
    #                                                        intersections[2][0]))
    #         primitives4.append(sp12)
    #         primitives4.append(sp21)
    #
    #         # Contour2D(primitives1),Contour2D(primitives2),
    #         #                      Contour2D(primitives3),
    #         all_contours.extend([volmdlr.wires.Contour2D(primitives4)])
    #
    #     else:
    #         raise NotImplementedError(
    #             f'{len(intersections)} intersections not supported yet')
    #
    #     return all_contours

    def bounding_rectangle(self):
        """
        Returns bounding rectangle.

        :return: Returns a python object with useful methods
        :rtype: :class:`volmdlr.core.BoundingRectangle
        """

        return self.outer_contour.bounding_rectangle

    @classmethod
    def from_contours(cls, outer_contour, inner_contours):
        """
        Create a Surface2D object from an outer contour and a list of inner contours.

        :param outer_contour: The outer contour that bounds the surface.
        :type outer_contour: volmdlr.wires.Contour2D
        :param inner_contours: The list of inner contours that define the holes of the surface.
        :type inner_contours : List[volmdlr.wires.Contour2D]
        :return: Surface2D defined by the given contours.
        """
        surface2d_inner_contours = []
        surface2d_outer_contour = outer_contour
        for inner_contour in inner_contours:
            if surface2d_outer_contour.shared_primitives_extremities(
                    inner_contour):
                # inner_contour will be merged with outer_contour
                merged_contours = surface2d_outer_contour.merge_with(
                    inner_contour)
                if len(merged_contours) >= 2:
                    raise NotImplementedError
                surface2d_outer_contour = merged_contours[0]
            else:
                # inner_contour will be added to the inner contours of the
                # Surface2D
                surface2d_inner_contours.append(inner_contour)
        return cls(surface2d_outer_contour, surface2d_inner_contours)

    def plot(self, ax=None, color='k', alpha=1, equal_aspect=False):

        if ax is None:
            _, ax = plt.subplots()
        self.outer_contour.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha,
                                                            equal_aspect=equal_aspect))
        for inner_contour in self.inner_contours:
            inner_contour.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha,
                                                           equal_aspect=equal_aspect))

        if equal_aspect:
            ax.set_aspect('equal')

        ax.margins(0.1)
        return ax

    def axial_symmetry(self, line):
        """
        Finds out the symmetric 2D surface according to a line.

        """

        outer_contour = self.outer_contour.axial_symmetry(line)
        inner_contours = []
        if self.inner_contours:
            inner_contours = [contour.axial_symmetry(line) for contour in self.inner_contours]

        return self.__class__(outer_contour=outer_contour,
                              inner_contours=inner_contours)

    def rotation(self, center, angle):
        """
        Surface2D rotation.

        :param center: rotation center.
        :param angle: angle rotation.
        :return: a new rotated Surface2D.
        """

        outer_contour = self.outer_contour.rotation(center, angle)
        if self.inner_contours:
            inner_contours = [contour.rotation(center, angle) for contour in self.inner_contours]
        else:
            inner_contours = []

        return self.__class__(outer_contour, inner_contours)

    def rotation_inplace(self, center, angle):
        """
        Rotate the surface inplace.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_surface2d = self.rotation(center, angle)
        self.outer_contour = new_surface2d.outer_contour
        self.inner_contours = new_surface2d.inner_contours

    def translation(self, offset: volmdlr.Vector2D):
        """
        Surface2D translation.

        :param offset: translation vector.
        :return: A new translated Surface2D.
        """
        outer_contour = self.outer_contour.translation(offset)
        inner_contours = [contour.translation(offset) for contour in self.inner_contours]
        return self.__class__(outer_contour, inner_contours)

    def translation_inplace(self, offset: volmdlr.Vector2D):
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_contour = self.translation(offset)
        self.outer_contour = new_contour.outer_contour
        self.inner_contours = new_contour.inner_contours

    def frame_mapping(self, frame: volmdlr.Frame2D, side: str):
        outer_contour = self.outer_contour.frame_mapping(frame, side)
        inner_contours = [contour.frame_mapping(frame, side) for contour in self.inner_contours]
        return self.__class__(outer_contour, inner_contours)

    def frame_mapping_inplace(self, frame: volmdlr.Frame2D, side: str):
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_contour = self.frame_mapping(frame, side)
        self.outer_contour = new_contour.outer_contour
        self.inner_contours = new_contour.inner_contours

    def geo_lines(self):  # , mesh_size_list=None):
        """
        Gets the lines that define a Surface2D in a .geo file.
        """

        i, i_p = None, None
        lines, line_surface, lines_tags = [], [], []
        point_account, line_account, line_loop_account = 0, 0, 1
        for outer_contour, contour in enumerate(list(chain(*[[self.outer_contour], self.inner_contours]))):

            if isinstance(contour, volmdlr.wires.Circle2D):
                points = [volmdlr.Point2D(contour.center.x - contour.radius, contour.center.y),
                          contour.center,
                          volmdlr.Point2D(contour.center.x + contour.radius, contour.center.y)]
                index = []
                for i, point in enumerate(points):
                    lines.append(point.get_geo_lines(tag=point_account + i + 1,
                                                     point_mesh_size=None))
                    index.append(point_account + i + 1)

                lines.append('Circle(' + str(line_account + 1) +
                             ') = {' + str(index[0]) + ', ' + str(index[1]) + ', ' + str(index[2]) + '};')
                lines.append('Circle(' + str(line_account + 2) +
                             ') = {' + str(index[2]) + ', ' + str(index[1]) + ', ' + str(index[0]) + '};')

                lines_tags.append(line_account + 1)
                lines_tags.append(line_account + 2)

                lines.append('Line Loop(' + str(outer_contour + 1) + ') = {' + str(lines_tags)[1:-1] + '};')
                line_surface.append(line_loop_account)

                point_account = point_account + 2 + 1
                line_account, line_loop_account = line_account + 1 + 1, line_loop_account + 1
                lines_tags = []

            elif isinstance(contour, (volmdlr.wires.Contour2D, volmdlr.wires.ClosedPolygon2D)):
                if not isinstance(contour, volmdlr.wires.ClosedPolygon2D):
                    contour = contour.to_polygon(1)
                for i, point in enumerate(contour.points):
                    lines.append(point.get_geo_lines(tag=point_account + i + 1,
                                                     point_mesh_size=None))

                for i_p, primitive in enumerate(contour.primitives):
                    if i_p != len(contour.primitives) - 1:
                        lines.append(primitive.get_geo_lines(tag=line_account + i_p + 1,
                                                             start_point_tag=point_account + i_p + 1,
                                                             end_point_tag=point_account + i_p + 2))
                    else:
                        lines.append(primitive.get_geo_lines(tag=line_account + i_p + 1,
                                                             start_point_tag=point_account + i_p + 1,
                                                             end_point_tag=point_account + 1))
                    lines_tags.append(line_account + i_p + 1)

                lines.append('Line Loop(' + str(outer_contour + 1) + ') = {' + str(lines_tags)[1:-1] + '};')
                line_surface.append(line_loop_account)
                point_account = point_account + i + 1
                line_account, line_loop_account = line_account + i_p + 1, line_loop_account + 1
                lines_tags = []

        lines.append('Plane Surface(' + str(1) + ') = {' + str(line_surface)[1:-1] + '};')

        return lines

    def mesh_lines(self,
                   factor: float,
                   curvature_mesh_size: int = None,
                   min_points: int = None,
                   initial_mesh_size: float = 5):
        """
        Gets the lines that define mesh parameters for a Surface2D, to be added to a .geo file.

        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A list of lines that describe mesh parameters
        :rtype: List[str]
        """

        lines = []
        if factor == 0:
            factor = 1e-3

        size = (math.sqrt(self.area()) / initial_mesh_size) * factor

        if min_points:
            primitives, primitives_length = [], []
            for _, contour in enumerate(list(chain(*[[self.outer_contour], self.inner_contours]))):
                if isinstance(contour, volmdlr.wires.Circle2D):
                    primitives.append(contour)
                    primitives.append(contour)
                    primitives_length.append(contour.length() / 2)
                    primitives_length.append(contour.length() / 2)
                else:
                    for primitive in contour.primitives:
                        if ((primitive not in primitives)
                                and (primitive.reverse() not in primitives)):
                            primitives.append(primitive)
                            primitives_length.append(primitive.length())

            for i, length in enumerate(primitives_length):
                if length < min_points * size:
                    lines.append('Transfinite Curve {' + str(i) + '} = ' + str(min_points) + ' Using Progression 1;')

        lines.append('Field[1] = MathEval;')
        lines.append('Field[1].F = "' + str(size) + '";')
        lines.append('Background Field = 1;')
        if curvature_mesh_size:
            lines.append('Mesh.MeshSizeFromCurvature = ' + str(curvature_mesh_size) + ';')

        # lines.append('Coherence;')

        return lines

    def to_geo(self, file_name: str,
               factor: float, **kwargs):
        # curvature_mesh_size: int = None,
        # min_points: int = None,
        # initial_mesh_size: float = 5):
        """
        Gets the .geo file for the Surface2D.
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        lines = self.geo_lines()
        lines.extend(self.mesh_lines(factor, kwargs['curvature_mesh_size'],
                                     kwargs['min_points'], kwargs['initial_mesh_size']))

        with open(file_name + '.geo', 'w', encoding="utf-8") as f:
            for line in lines:
                f.write(line)
                f.write('\n')
        f.close()

    def to_msh(self, file_name: str, mesh_dimension: int,
               factor: float, **kwargs):
        # curvature_mesh_size: int = 0,
        # min_points: int = None,
        # initial_mesh_size: float = 5):
        """
        Gets .msh file for the Surface2D generated by gmsh.

        :param file_name: The msh. file name
        :type file_name: str
        :param mesh_dimension: The mesh dimension (1: 1D-Edge, 2: 2D-Triangle, 3D-Tetrahedra)
        :type mesh_dimension: int
        :param factor: A float, between 0 and 1, that describes the mesh quality
        (1 for coarse mesh - 0 for fine mesh)
        :type factor: float
        :param curvature_mesh_size: Activate the calculation of mesh element sizes based on curvature
        (with curvature_mesh_size elements per 2*Pi radians), defaults to 0
        :type curvature_mesh_size: int, optional
        :param min_points: Check if there are enough points on small edges (if it is not, we force to have min_points
        on that edge), defaults to None
        :type min_points: int, optional
        :param initial_mesh_size: If factor=1, it will be initial_mesh_size elements per dimension, defaults to 5
        :type initial_mesh_size: float, optional

        :return: A txt file
        :rtype: .txt
        """

        for element in [('curvature_mesh_size', 0), ('min_points', None), ('initial_mesh_size', 5)]:
            if element[0] not in kwargs:
                kwargs[element[0]] = element[1]

        self.to_geo(file_name=file_name, mesh_dimension=mesh_dimension,
                    factor=factor, curvature_mesh_size=kwargs['curvature_mesh_size'],
                    min_points=kwargs['min_points'], initial_mesh_size=kwargs['initial_mesh_size'])

        volmdlr.core.VolumeModel.generate_msh_file(file_name, mesh_dimension)

        # gmsh.initialize()
        # gmsh.open(file_name + ".geo")

        # gmsh.model.geo.synchronize()
        # gmsh.model.mesh.generate(mesh_dimension)

        # gmsh.write(file_name + ".msh")

        # gmsh.finalize()


class Surface3D(DessiaObject):
    """
    Abstract class.

    """
    x_periodicity = None
    y_periodicity = None
    face_class = None

    def point2d_to_3d(self, point2d):
        raise NotImplementedError(f'point2d_to_3d is abstract and should be implemented in {self.__class__.__name__}')

    def point3d_to_2d(self, point3d):
        """
        Abstract method. Convert a 3D point to a 2D parametric point.

        :param point3d: The 3D point to convert, represented by 3 coordinates (x, y, z).
        :type point3d: `volmdlr.Point3D`
        :return: NotImplementedError: If the method is not implemented in the subclass.
        """
        raise NotImplementedError(f'point3d_to_2d is abstract and should be implemented in {self.__class__.__name__}')

    def face_from_contours3d(self, contours3d: List[volmdlr.wires.Contour3D], name: str = ''):
        """
        Returns the face generated by a list of contours. Finds out which are outer or inner contours.

        :param name: the name to inject in the new face
        """

        lc3d = len(contours3d)

        if lc3d == 1:
            outer_contour2d = self.contour3d_to_2d(contours3d[0])
            inner_contours2d = []
        elif lc3d > 1:
            area = -1
            inner_contours2d = []

            contours2d = [self.contour3d_to_2d(contour3d) for contour3d in contours3d]

            check_contours = [not contour2d.is_ordered() for contour2d in contours2d]
            if any(check_contours):
                outer_contour2d, inner_contours2d = self.repair_contours2d(contours2d[0], contours2d[1:])
            else:
                for contour2d in contours2d:
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
        return class_(self, surface2d=surface2d, name=name)

    def repair_primitives_periodicity(self, primitives2d):
        """
        Repairs the continuity of the 2D contour while using contour3d_to_2d on periodic surfaces.

        :param primitives2d: The primitives in parametric surface domain.
        :type primitives2d: list
        :return: A list of primitives.
        :rtype: list
        """
        x_periodicity = self.x_periodicity
        y_periodicity = self.y_periodicity
        # Search for a primitive that can be used as reference for repairing periodicity
        if x_periodicity or y_periodicity:
            pos = vm_parametric.find_index_defined_brep_primitive_on_periodical_surface(primitives2d,
                                                                                        [x_periodicity, y_periodicity])
            if pos != 0:
                primitives2d = primitives2d[pos:] + primitives2d[:pos]

        i = 1
        if x_periodicity is None:
            x_periodicity = -1
        if y_periodicity is None:
            y_periodicity = -1
        while i < len(primitives2d):
            previous_primitive = primitives2d[i - 1]
            delta = previous_primitive.end - primitives2d[i].start
            is_connected = math.isclose(delta.norm(), 0, abs_tol=1e-5)
            if not is_connected and \
                    primitives2d[i].end.is_close(primitives2d[i - 1].end, tol=1e-3) and \
                    math.isclose(primitives2d[i].length(), x_periodicity, abs_tol=1e-5):
                primitives2d[i] = primitives2d[i].reverse()
            elif not is_connected and \
                    primitives2d[i].end.is_close(primitives2d[i - 1].end, tol=1e-3) and \
                    math.isclose(primitives2d[i].length(), y_periodicity, abs_tol=1e-5):
                primitives2d[i] = primitives2d[i].reverse()
            elif not is_connected:
                primitives2d[i] = primitives2d[i].translation(delta)
            i += 1
        last_end = primitives2d[-1].end
        first_start = primitives2d[0].start
        if not last_end.is_close(first_start, tol=1e-5):
            deltax = first_start.x - last_end.x
            deltay = first_start.y - last_end.y
            if deltax and x_periodicity and abs(deltax) % x_periodicity == 0:
                primitives2d[-1] = vme.LineSegment2D(primitives2d[-1].start,
                                                     volmdlr.Point2D(primitives2d[-1].end.x + deltax,
                                                                     primitives2d[-1].end.y))
            elif deltay and y_periodicity and abs(deltay) % y_periodicity == 0:
                primitives2d[-1] = vme.LineSegment2D(primitives2d[-1].start,
                                                     volmdlr.Point2D(primitives2d[-1].end.x,
                                                                     primitives2d[-1].end.y + deltay))

        return primitives2d

    def repair_contours2d(self, outer_contour, inner_contours):
        """
        Abstract method. Repair 2D contours of a face on the parametric domain.

        :param outer_contour: Outer contour 2D.
        :type inner_contours: volmdlr.wires.Contour2D
        :param inner_contours: List of 2D contours.
        :type inner_contours: list
        :return: NotImplementedError: If the method is not implemented in the subclass.
        """
        raise NotImplementedError(f'repair_contours2d is abstract and should be implemented in '
                                  f'{self.__class__.__name__}')

    def primitives3d_to_2d(self, primitives3d):
        """
        Helper function to perform convertion of 3D primitives into B-Rep primitives.

        :param primitives3d: List of 3D primitives from a 3D contour.
        :type primitives3d: List[vme.Edge]
        :return: A list of 2D primitives on parametric domain.
        :rtype: List[vme.Edge]
        """
        primitives2d = []
        for primitive3d in primitives3d:
            method_name = f'{primitive3d.__class__.__name__.lower()}_to_2d'
            if hasattr(self, method_name):
                primitives = getattr(self, method_name)(primitive3d)

                if primitives is None:
                    continue
                primitives2d.extend(primitives)
            else:
                raise NotImplementedError(f'Class {self.__class__.__name__} does not implement {method_name}')
        return primitives2d

    def contour3d_to_2d(self, contour3d):
        """
        Transforms a Contour3D into a Contour2D in the parametric domain of the surface.

        :param contour3d: The contour to be transformed.
        :type contour3d: :class:`volmdlr.wires.Contour3D`
        :return: A 2D contour object.
        :rtype: :class:`volmdlr.wires.Contour2D`
        """
        primitives2d = self.primitives3d_to_2d(contour3d.primitives)

        wire2d = volmdlr.wires.Wire2D(primitives2d)
        delta_x = abs(wire2d.primitives[0].start.x - wire2d.primitives[-1].end.x)
        if math.isclose(delta_x, volmdlr.TWO_PI, abs_tol=1e-3) and wire2d.is_ordered():
            return volmdlr.wires.Contour2D(primitives2d)
        # Fix contour
        if self.x_periodicity or self.y_periodicity:
            primitives2d = self.repair_primitives_periodicity(primitives2d)
        return volmdlr.wires.Contour2D(primitives2d)

    def contour2d_to_3d(self, contour2d):
        """
        Transforms a Contour2D in the parametric domain of the surface into a Contour3D in Cartesian coordinate.

        :param contour2d: The contour to be transformed.
        :type contour2d: :class:`volmdlr.wires.Contour2D`
        :return: A 3D contour object.
        :rtype: :class:`volmdlr.wires.Contour3D`
        """
        primitives3d = []
        for primitive2d in contour2d.primitives:
            method_name = f'{primitive2d.__class__.__name__.lower()}_to_3d'
            if hasattr(self, method_name):
                try:
                    primitives = getattr(self, method_name)(primitive2d)
                    if primitives is None:
                        continue
                    primitives3d.extend(primitives)
                except NotImplementedError:
                    print(f'Class {self.__class__.__name__} does not implement {method_name}'
                          f'with {primitive2d.__class__.__name__}')
            else:
                raise NotImplementedError(
                    f'Class {self.__class__.__name__} does not implement {method_name}')

        return volmdlr.wires.Contour3D(primitives3d)

    def linesegment3d_to_2d(self, linesegment3d):
        """
        A line segment on a surface will be in any case a line in 2D?.

        """
        return [vme.LineSegment2D(self.point3d_to_2d(linesegment3d.start),
                                  self.point3d_to_2d(linesegment3d.end))]

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        """
        Is this right?.
        """
        n = len(bspline_curve3d.control_points)
        points = [self.point3d_to_2d(p)
                  for p in bspline_curve3d.discretization_points(number_points=n)]
        return [vme.BSplineCurve2D.from_points_interpolation(points, bspline_curve3d.degree, bspline_curve3d.periodic)]

    def bsplinecurve2d_to_3d(self, bspline_curve2d):
        """
        Is this right?.

        """
        n = len(bspline_curve2d.control_points)
        points = [self.point2d_to_3d(p)
                  for p in bspline_curve2d.discretization_points(number_points=n)]
        return [vme.BSplineCurve3D.from_points_interpolation(points, bspline_curve2d.degree, bspline_curve2d.periodic)]

    def normal_from_point2d(self, point2d):

        raise NotImplementedError('NotImplemented')

    def normal_from_point3d(self, point3d):
        """
        Evaluates the normal vector of the bspline surface at this 3D point.
        """

        return (self.normal_from_point2d(self.point3d_to_2d(point3d)))[1]

    def geodesic_distance_from_points2d(self, point1_2d: volmdlr.Point2D,
                                        point2_2d: volmdlr.Point2D, number_points: int = 50):
        """
        Approximation of geodesic distance via line segments length sum in 3D.
        """
        # points = [point1_2d]
        current_point3d = self.point2d_to_3d(point1_2d)
        distance = 0.
        for i in range(number_points):
            next_point3d = self.point2d_to_3d(point1_2d + (i + 1) / number_points * (point2_2d - point1_2d))
            distance += next_point3d.point_distance(current_point3d)
            current_point3d = next_point3d
        return distance

    def geodesic_distance(self, point1_3d: volmdlr.Point3D, point2_3d: volmdlr.Point3D):
        """
        Approximation of geodesic distance between 2 3D points supposed to be on the surface.
        """
        point1_2d = self.point3d_to_2d(point1_3d)
        point2_2d = self.point3d_to_2d(point2_3d)
        return self.geodesic_distance_from_points2d(point1_2d, point2_2d)


class Plane3D(Surface3D):
    """
    Defines a plane 3d.

    :param frame: u and v of frame describe the plane, w is the normal
    """
    face_class = 'PlaneFace3D'

    def __init__(self, frame: volmdlr.Frame3D, name: str = ''):

        self.frame = frame
        self.name = name
        Surface3D.__init__(self, name=name)

    def __hash__(self):
        return hash(('plane 3d', self.frame))

    def __eq__(self, other_plane):
        if other_plane.__class__.__name__ != self.__class__.__name__:
            return False
        return self.frame == other_plane.frame

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a Plane3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated
        :type object_dict: dict
        :return: The corresponding Plane3D object.
        :rtype: :class:`volmdlr.faces.Plane3D`
        """
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
        content += f"#{plane_id} = PLANE('{self.name}',#{frame_id});\n"
        return content, [plane_id]

    @classmethod
    def from_3_points(cls, *args):
        """
        Point 1 is used as origin of the plane.
        """
        point1, point2, point3 = args
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
    def from_plane_vectors(cls, plane_origin: volmdlr.Point3D, plane_x: volmdlr.Vector3D, plane_y: volmdlr.Vector3D):
        """
        Initializes a 3D plane object with a given plane origin and plane x and y vectors.

        :param plane_origin: A volmdlr.Point3D representing the origin of the plane.
        :param plane_x: A volmdlr.Vector3D representing the x-axis of the plane.
        :param plane_y: A volmdlr.Vector3D representing the y-axis of the plane.
        :return: A Plane3D object initialized from the provided plane origin and plane x and y vectors.
        """
        normal = plane_x.cross(plane_y)
        return cls(volmdlr.Frame3D(plane_origin, plane_x, plane_y, normal))

    @classmethod
    def from_points(cls, points):
        """
        Returns a 3D plane that goes through the 3 first points on the list.

        Why for more than 3 points we only do some check and never raise error?
        """
        if len(points) < 3:
            raise ValueError
        if len(points) == 3:
            return cls.from_3_points(volmdlr.Point3D(points[0].vector),
                                     volmdlr.Vector3D(points[1].vector),
                                     volmdlr.Vector3D(points[2].vector))
        points = [p.copy() for p in points]
        indexes_to_del = []
        for i, point in enumerate(points[1:]):
            if point.is_close(points[0]):
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
        return cls.from_3_points(origin, vector1 + origin, vector2_min + origin)

    def point_on_surface(self, point):
        """
        Return if the point belongs to the plane at a tolerance of 1e-6.

        """
        if math.isclose(self.frame.w.dot(point - self.frame.origin), 0,
                        abs_tol=1e-6):
            return True
        return False

    def point_distance(self, point3d):
        """
        Calculates the distance of a point to plane.

        :param point3d: point to verify distance.
        :return: a float, point distance to plane.
        """
        coefficient_a, coefficient_b, coefficient_c, coefficient_d = self.equation_coefficients()
        return abs(self.frame.w.dot(point3d) + coefficient_d) / math.sqrt(coefficient_a ** 2 +
                                                                          coefficient_b ** 2 + coefficient_c ** 2)

    def line_intersections(self, line):
        """
        Find the intersection with a line.

        :param line: Line to evaluate the intersection
        :type line: :class:`vme.Line`
        :return: ADD DESCRIPTION
        :rtype: List[volmdlr.Point3D]
        """
        u_vector = line.point2 - line.point1
        w_vector = line.point1 - self.frame.origin
        if math.isclose(self.frame.w.dot(u_vector), 0, abs_tol=1e-08):
            return []
        intersection_abscissea = - self.frame.w.dot(w_vector) / self.frame.w.dot(u_vector)
        return [line.point1 + intersection_abscissea * u_vector]

    def linesegment_intersections(self, linesegment: vme.LineSegment3D) \
            -> List[volmdlr.Point3D]:
        u_vector = linesegment.end - linesegment.start
        w_vector = linesegment.start - self.frame.origin
        normaldotu = self.frame.w.dot(u_vector)
        if math.isclose(normaldotu, 0, abs_tol=1e-08):
            return []
        intersection_abscissea = - self.frame.w.dot(w_vector) / normaldotu
        if intersection_abscissea < 0 or intersection_abscissea > 1:
            return []
        return [linesegment.start + intersection_abscissea * u_vector]

    def fullarc_intersections(self, fullarc: vme.FullArc3D):
        """
        Calculates the intersections between a Plane 3D and a FullArc 3D.

        :param fullarc: full arc to verify intersections.
        :return: list of intersections: List[volmdlr.Point3D].
        """
        fullarc_plane = Plane3D(fullarc.frame)
        plane_intersections = self.plane_intersection(fullarc_plane)
        if not plane_intersections:
            return []
        fullarc2d = fullarc.to_2d(fullarc.center, fullarc_plane.frame.u, fullarc_plane.frame.v)
        line2d = plane_intersections[0].to_2d(fullarc.center, fullarc_plane.frame.u, fullarc_plane.frame.v)
        fullarc2d_inters_line2d = fullarc2d.line_intersections(line2d)
        intersections = []
        for inter in fullarc2d_inters_line2d:
            intersections.append(inter.to_3d(fullarc.center, fullarc_plane.frame.u, fullarc_plane.frame.v))
        return intersections

    def equation_coefficients(self):
        """
        Returns the a,b,c,d coefficient from equation ax+by+cz+d = 0.

        """
        a, b, c = self.frame.w
        d = -self.frame.origin.dot(self.frame.w)
        return round(a, 12), round(b, 12), round(c, 12), round(d, 12)

    def plane_intersection(self, other_plane):
        """
        Computes intersection points between two Planes 3D.

        """
        if self.is_parallel(other_plane):
            return []
        line_direction = self.frame.w.cross(other_plane.frame.w)

        if line_direction.norm() < 1e-6:
            return None

        a1, b1, c1, d1 = self.equation_coefficients()
        a2, b2, c2, d2 = other_plane.equation_coefficients()
        if not math.isclose(a1 * b2 - a2 * b1, 0.0, abs_tol=1e-10):
            x0 = (b1 * d2 - b2 * d1) / (a1 * b2 - a2 * b1)
            y0 = (a2 * d1 - a1 * d2) / (a1 * b2 - a2 * b1)
            point1 = volmdlr.Point3D(x0, y0, 0)
        elif a2 * c1 != a1 * c2:
            x0 = (c2 * d1 - c1 * d2) / (a2 * c1 - a1 * c2)
            z0 = (a1 * d2 - a2 * d1) / (a2 * c1 - a1 * c2)
            point1 = volmdlr.Point3D(x0, 0, z0)
        elif c1 * b2 != b1 * c2:
            y0 = (- c2 * d1 + c1 * d2) / (b1 * c2 - c1 * b2)
            z0 = (- b1 * d2 + b2 * d1) / (b1 * c2 - c1 * b2)
            point1 = volmdlr.Point3D(0, y0, z0)
        else:
            raise NotImplementedError
        return [vme.Line3D(point1, point1 + line_direction)]

    def is_coincident(self, plane2):
        """
        Verifies if two planes are parallel and coincident.

        """
        if not isinstance(self, plane2.__class__):
            return False
        if self.is_parallel(plane2):
            if plane2.point_on_surface(self.frame.origin):
                return True
        return False

    def is_parallel(self, plane2):
        """
        Verifies if two planes are parallel.

        """
        if self.frame.w.is_colinear_to(plane2.frame.w):
            return True
        return False

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Plane3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Plane3D
        """
        new_frame = self.frame.rotation(center=center, axis=axis, angle=angle)
        return Plane3D(new_frame)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Plane3D rotation. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.frame.rotation_inplace(center=center, axis=axis, angle=angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Plane3D translation.

        :param offset: translation vector
        :return: A new translated Plane3D
        """
        new_frame = self.frame.translation(offset)
        return Plane3D(new_frame)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Plane3D translation. Object is updated inplace.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.frame.translation_inplace(offset)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Frame3D.

        :param frame: Frame of reference
        :type frame: `volmdlr.Frame3D`
        :param side: 'old' or 'new'
        """
        new_frame = self.frame.frame_mapping(frame, side)
        return Plane3D(new_frame, self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        :param frame: Frame of reference
        :type frame: `volmdlr.Frame3D`
        :param side: 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_frame = self.frame.frame_mapping(frame, side)
        self.frame.origin = new_frame.origin
        self.frame.u = new_frame.u
        self.frame.v = new_frame.v
        self.frame.w = new_frame.w

    def copy(self, deep=True, memo=None):
        new_frame = self.frame.copy(deep, memo)
        return Plane3D(new_frame, self.name)

    def plot(self, ax=None):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')

        self.frame.origin.plot(ax)
        self.frame.u.plot(ax, color='r')
        self.frame.v.plot(ax, color='g')
        return ax

    def point2d_to_3d(self, point2d):
        """
        Converts a 2D parametric point into a 3D point on the surface.
        """
        return point2d.to_3d(self.frame.origin, self.frame.u, self.frame.v)

    def point3d_to_2d(self, point3d):
        """
        Converts a 3D point into a 2D parametric point.
        """
        return point3d.to_2d(self.frame.origin, self.frame.u, self.frame.v)

    def contour2d_to_3d(self, contour2d):
        """
        Converts a contour 2D on parametric surface into a 3D contour.
        """
        return contour2d.to_3d(self.frame.origin, self.frame.u, self.frame.v)

    def contour3d_to_2d(self, contour3d):
        """
        Converts a contour 3D into a 2D parametric contour.
        """
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
        """
        Converts a 2D B-Spline in parametric domain into a 3D B-Spline in spatial domain.

        :param bspline_curve2d: The B-Spline curve to perform the transformation.
        :type bspline_curve2d: vme.BSplineCurve2D
        :return: A 3D B-Spline.
        :rtype: vme.BSplineCurve3D
        """
        control_points = [self.point2d_to_3d(point)
                          for point in bspline_curve2d.control_points]
        return [vme.BSplineCurve3D(
            bspline_curve2d.degree,
            control_points=control_points,
            knot_multiplicities=bspline_curve2d.knot_multiplicities,
            knots=bspline_curve2d.knots,
            weights=bspline_curve2d.weights,
            periodic=bspline_curve2d.periodic)]

    def rectangular_cut(self, x1: float, x2: float,
                        y1: float, y2: float, name: str = ''):
        """
        Cut a rectangular piece of the Plane3D object and return a PlaneFace3D object.

        """
        point1 = volmdlr.Point2D(x1, y1)
        point2 = volmdlr.Point2D(x2, y1)
        point3 = volmdlr.Point2D(x2, y2)
        point4 = volmdlr.Point2D(x1, y2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        surface = Surface2D(outer_contour, [])
        return PlaneFace3D(self, surface, name)


PLANE3D_OXY = Plane3D(volmdlr.OXYZ)
PLANE3D_OYZ = Plane3D(volmdlr.OYZX)
PLANE3D_OZX = Plane3D(volmdlr.OZXY)


class PeriodicalSurface(Surface3D):
    """
    Abstract class for surfaces with two-pi periodicity that creates some problems.
    """

    def point2d_to_3d(self, point2d):
        raise NotImplementedError(f'point2d_to_3d is abstract and should be implemented in {self.__class__.__name__}')

    def point3d_to_2d(self, point3d):
        """
        Abstract method. Convert a 3D point to a 2D parametric point.

        :param point3d: The 3D point to convert, represented by 3 coordinates (x, y, z).
        :type point3d: `volmdlr.Point3D`
        :return: NotImplementedError: If the method is not implemented in the subclass.
        """
        raise NotImplementedError(f'point3d_to_2d is abstract and should be implemented in {self.__class__.__name__}')

    def repair_contours2d(self, outer_contour, inner_contours):
        """
        Repair contours on parametric domain.

        :param outer_contour: Outer contour 2D.
        :type inner_contours: volmdlr.wires.Contour2D
        :param inner_contours: List of 2D contours.
        :type inner_contours: list
        :param name: the name to inject in the new face.
        """
        new_inner_contours = []
        point1 = outer_contour.primitives[0].start
        point2 = outer_contour.primitives[-1].end

        theta1, z1 = point1
        theta2, z2 = point2
        old_outer_contour_positioned = outer_contour
        new_outer_contour = old_outer_contour_positioned
        for inner_contour in inner_contours:
            theta3, z3 = inner_contour.primitives[0].start
            theta4, z4 = inner_contour.primitives[-1].end
            # check if inner_contour has a length of 2pi in theta.
            if math.isclose(abs(theta4 - theta3), 2 * math.pi, abs_tol=1e-3):

                outer_contour_theta = [theta1, theta2]
                inner_contour_theta = [theta3, theta4]

                # Contours are aligned
                if (math.isclose(theta1, theta3, abs_tol=1e-3) and math.isclose(theta2, theta4, abs_tol=1e-3)) \
                        or (math.isclose(theta1, theta4, abs_tol=1e-3) and math.isclose(theta2, theta3, abs_tol=1e-3)):
                    old_innner_contour_positioned = inner_contour

                else:
                    overlapping_theta, outer_contour_side, inner_contour_side = self._get_overlapping_theta(
                        outer_contour_theta,
                        inner_contour_theta)
                    line = vme.Line2D(volmdlr.Point2D(overlapping_theta, z1),
                                      volmdlr.Point2D(overlapping_theta, z3))
                    cutted_contours = inner_contour.split_by_line(line)
                    number_contours = len(cutted_contours)
                    if number_contours == 2:
                        contour1, contour2 = cutted_contours
                        increasing_theta = theta3 < theta4
                        # inner_contour_side = 0 --> left  inner_contour_side = 1 --> right
                        if (not inner_contour_side and increasing_theta) or (
                                inner_contour_side and not increasing_theta):
                            theta_offset = outer_contour_theta[outer_contour_side] - contour2.primitives[0].start.x
                            translation_vector = volmdlr.Vector2D(theta_offset, 0)
                            contour2_positionned = contour2.translation(offset=translation_vector)
                            theta_offset = contour2_positionned.primitives[-1].end.x - contour1.primitives[0].start.x
                            translation_vector = volmdlr.Vector2D(theta_offset, 0)
                            contour1_positionned = contour1.translation(offset=translation_vector)
                        else:
                            theta_offset = outer_contour_theta[outer_contour_side] - contour1.primitives[-1].end.x
                            translation_vector = volmdlr.Vector2D(theta_offset, 0)
                            contour1_positionned = contour1.translation(offset=translation_vector)
                            theta_offset = contour1_positionned.primitives[0].start.x - contour2.primitives[-1].end.x
                            translation_vector = volmdlr.Vector2D(theta_offset, 0)
                            contour2_positionned = contour2.translation(offset=translation_vector)

                        old_innner_contour_positioned = volmdlr.wires.Contour2D(contour1_positionned.primitives +
                                                                                contour2_positionned.primitives)
                        old_innner_contour_positioned.order_contour()
                    elif number_contours == 1:
                        contour = cutted_contours[0]
                        theta_offset = outer_contour_theta[outer_contour_side] - \
                                       inner_contour_theta[inner_contour_side]
                        translation_vector = volmdlr.Vector2D(theta_offset, 0)
                        old_innner_contour_positioned = contour.translation(offset=translation_vector)

                    else:
                        print(True)
                        raise NotImplementedError
                point1 = old_outer_contour_positioned.primitives[0].start
                point2 = old_outer_contour_positioned.primitives[-1].end
                point3 = old_innner_contour_positioned.primitives[0].start
                point4 = old_innner_contour_positioned.primitives[-1].end

                outer_contour_direction = point1.x < point2.x
                inner_contour_direction = point3.x < point4.x
                if outer_contour_direction == inner_contour_direction:
                    old_innner_contour_positioned = old_innner_contour_positioned.invert()
                    point3 = old_innner_contour_positioned.primitives[0].start
                    point4 = old_innner_contour_positioned.primitives[-1].end

                closing_linesegment1 = volmdlr.edges.LineSegment2D(point2, point3)
                closing_linesegment2 = volmdlr.edges.LineSegment2D(point4, point1)
                new_outer_contour_primitives = old_outer_contour_positioned.primitives + [closing_linesegment1] + \
                                               old_innner_contour_positioned.primitives + \
                                               [closing_linesegment2]
                new_outer_contour = volmdlr.wires.Contour2D(primitives=new_outer_contour_primitives)
                new_outer_contour.order_contour()
            else:
                new_inner_contours.append(inner_contour)
        return new_outer_contour, new_inner_contours

    def _get_overlapping_theta(self, outer_contour_startend_theta, inner_contour_startend_theta):
        """
        Find overlapping theta domain between two contours on periodical Surfaces.
        """
        oc_xmin_index, outer_contour_xmin = min(enumerate(outer_contour_startend_theta), key=lambda x: x[1])
        oc_xmax_index, outer_contour_xman = max(enumerate(outer_contour_startend_theta), key=lambda x: x[1])
        inner_contour_xmin = min(inner_contour_startend_theta)
        inner_contour_xmax = max(inner_contour_startend_theta)

        # check if tetha3 or theta4 is in [theta1, theta2] interval
        overlap = outer_contour_xmin <= inner_contour_xmax and outer_contour_xman >= inner_contour_xmin

        if overlap:
            if inner_contour_xmin < outer_contour_xmin:
                overlapping_theta = outer_contour_startend_theta[oc_xmin_index]
                outer_contour_side = oc_xmin_index
                side = 0
                return overlapping_theta, outer_contour_side, side
            overlapping_theta = outer_contour_startend_theta[oc_xmax_index]
            outer_contour_side = oc_xmax_index
            side = 1
            return overlapping_theta, outer_contour_side, side

        # if not direct intersection -> find intersection at periodicity
        if inner_contour_xmin < outer_contour_xmin:
            overlapping_theta = outer_contour_startend_theta[oc_xmin_index] - 2 * math.pi
            outer_contour_side = oc_xmin_index
            side = 0
            return overlapping_theta, outer_contour_side, side
        overlapping_theta = outer_contour_startend_theta[oc_xmax_index] + 2 * math.pi
        outer_contour_side = oc_xmax_index
        side = 1
        return overlapping_theta, outer_contour_side, side

    def _reference_points(self, edge):
        """
        Helper function to return points of reference on the edge to fix some parametric periodical discontinuities.
        """
        length = edge.length()
        point_after_start = self.point3d_to_2d(edge.point_at_abscissa(0.001 * length))
        point_before_end = self.point3d_to_2d(edge.point_at_abscissa(0.98 * length))
        theta3, _ = point_after_start
        theta4, _ = point_before_end
        if abs(theta3) == math.pi or abs(theta3) == 0.5 * math.pi:
            point_after_start = self.point3d_to_2d(edge.point_at_abscissa(0.002 * length))
        if abs(theta4) == math.pi or abs(theta4) == 0.5 * math.pi:
            point_before_end = self.point3d_to_2d(edge.point_at_abscissa(0.97 * length))
        return point_after_start, point_before_end

    def _verify_start_end_angles(self, edge, theta1, theta2):
        """
        Verify if there is some incoherence with start and end angles. If so, return fixed angles.
        """
        length = edge.length()
        theta3, _ = self.point3d_to_2d(edge.point_at_abscissa(0.001 * length))
        # make sure that the reference angle is not undefined
        if abs(theta3) == math.pi:
            theta3, _ = self.point3d_to_2d(edge.point_at_abscissa(0.002 * length))

        # Verify if theta1 or theta2 point should be -pi because atan2() -> ]-pi, pi]
        # And also atan2 discontinuity in 0.5 * math.pi
        if abs(theta1) == math.pi or abs(theta1) == 0.5 * math.pi:
            theta1 = repair_start_end_angle_periodicity(theta1, theta3)
        if abs(theta2) == math.pi or abs(theta2) == 0.5 * math.pi:
            theta4, _ = self.point3d_to_2d(edge.point_at_abscissa(0.98 * length))
            # make sure that the reference angle is not undefined
            if abs(theta4) == math.pi:
                theta4, _ = self.point3d_to_2d(edge.point_at_abscissa(0.97 * length))
            theta2 = repair_start_end_angle_periodicity(theta2, theta4)

        return theta1, theta2

    def _fix_angle_discontinuity_on_discretization_points(self, points, indexes_angle_discontinuity, direction):
        # angle1, angle2, angle3 = angle_list
        # if angle3 < angle1 < angle2:
        #     points = [p - volmdlr.Point2D(volmdlr.TWO_PI, 0) if p.x > 0 else p for p in points]
        #     if angle2 == 0.0:
        #         points[-1] = volmdlr.Point2D(-volmdlr.TWO_PI, points[-1].y)
        # elif angle3 > angle1 > angle2:
        #     points = [p + volmdlr.Point2D(volmdlr.TWO_PI, 0) if p.x < 0 else p for p in points]
        #     if angle2 == 0.0:
        #         points[-1] = volmdlr.Point2D(volmdlr.TWO_PI, points[-1].y)
        i = 0 if direction == "x" else 1
        if len(indexes_angle_discontinuity) == 1:
            index_angle_discontinuity = indexes_angle_discontinuity[0]
            sign = round(points[index_angle_discontinuity - 1][i] / abs(points[index_angle_discontinuity - 1][i]), 2)
            if i == 0:
                points = [p + volmdlr.Point2D(sign * volmdlr.TWO_PI, 0) if i >= index_angle_discontinuity else p
                          for i, p in enumerate(points)]
            else:
                points = [p + volmdlr.Point2D(0, sign * volmdlr.TWO_PI) if i >= index_angle_discontinuity else p
                          for i, p in enumerate(points)]
        else:
            raise NotImplementedError
        return points

    def linesegment3d_to_2d(self, linesegment3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        start = self.point3d_to_2d(linesegment3d.start)
        end = self.point3d_to_2d(linesegment3d.end)
        if start.x != end.x:
            end = volmdlr.Point2D(start.x, end.y)
        if not start.is_close(end):
            return [vme.LineSegment2D(start, end)]
        return None

    def arc3d_to_2d(self, arc3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        start = self.point3d_to_2d(arc3d.start)
        end = self.point3d_to_2d(arc3d.end)
        angle3d = arc3d.angle
        point_after_start, point_before_end = self._reference_points(arc3d)

        start, end = vm_parametric.arc3d_to_cylindrical_coordinates_verification(start, end, angle3d,
                                                                                 point_after_start.x,
                                                                                 point_before_end.x)
        return [vme.LineSegment2D(start, end)]

    def fullarc3d_to_2d(self, fullarc3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        start = self.point3d_to_2d(fullarc3d.start)
        end = self.point3d_to_2d(fullarc3d.end)

        point_after_start, point_before_end = self._reference_points(fullarc3d)

        start, end = vm_parametric.arc3d_to_cylindrical_coordinates_verification(start, end, volmdlr.TWO_PI,
                                                                                 point_after_start.x,
                                                                                 point_before_end.x)
        theta1, z1 = start
        theta2, z2 = end
        theta3, z3 = point_after_start

        if self.frame.w.is_colinear_to(fullarc3d.normal):
            if start.is_close(end):
                point1 = volmdlr.Point2D(theta1, z1)
                if theta1 > theta3:
                    point2 = volmdlr.Point2D(theta1 + volmdlr.TWO_PI, z2)
                elif theta1 < theta3:
                    point2 = volmdlr.Point2D(theta1 - volmdlr.TWO_PI, z2)
                return [vme.LineSegment2D(point1, point2)]
            return [vme.LineSegment2D(start, end)]
        # Treating one case from Revolution Surface
        if z1 > z3:
            point1 = volmdlr.Point2D(theta1, 1)
            point2 = volmdlr.Point2D(theta1, 0)
        else:
            point1 = volmdlr.Point2D(theta1, 0)
            point2 = volmdlr.Point2D(theta1, 1)
        return [vme.LineSegment2D(point1, point2)]

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        n = len(bspline_curve3d.control_points)
        points3d = bspline_curve3d.discretization_points(number_points=n)
        points = [self.point3d_to_2d(point) for point in points3d]
        theta1, z1 = points[0]
        theta2, z2 = points[-1]
        theta1, theta2 = self._verify_start_end_angles(bspline_curve3d, theta1, theta2)
        points[0] = volmdlr.Point2D(theta1, z1)
        points[-1] = volmdlr.Point2D(theta2, z2)

        theta_list = [point.x for point in points]
        theta_discontinuity, indexes_theta_discontinuity = angle_discontinuity(theta_list)
        if theta_discontinuity:
            points = self._fix_angle_discontinuity_on_discretization_points(points,
                                                                            indexes_theta_discontinuity, "x")

        return [vme.BSplineCurve2D.from_points_interpolation(points, degree=bspline_curve3d.degree,
                                                             periodic=bspline_curve3d.periodic)]

    def arcellipse3d_to_2d(self, arcellipse3d):
        """
        Transformation of a 3D arc of ellipse to a 2D primitive in a cylindrical surface.

        """
        points = [self.point3d_to_2d(p)
                  for p in arcellipse3d.discretization_points(number_points=50)]

        theta1, z1 = points[0]
        theta2, z2 = points[-1]

        theta1, theta2 = self._verify_start_end_angles(arcellipse3d, theta1, theta2)
        points[0] = volmdlr.Point2D(theta1, z1)
        points[-1] = volmdlr.Point2D(theta2, z2)

        theta_list = [point.x for point in points]
        theta_discontinuity, indexes_theta_discontinuity = angle_discontinuity(theta_list)
        if theta_discontinuity:
            points = self._fix_angle_discontinuity_on_discretization_points(points,
                                                                            indexes_theta_discontinuity, "x")

        bsplinecurve2d = vme.BSplineCurve2D.from_points_interpolation(points, degree=2)
        return [bsplinecurve2d]

    def fullarcellipse3d_to_2d(self, arcellipse3d):
        """
        Transformation of a 3D arcellipse to 2D, in a cylindrical surface.

        """
        points = [self.point3d_to_2d(p)
                  for p in arcellipse3d.discretization_points(number_points=100)]
        points.pop()
        theta1, z1 = points[0]
        theta2, z2 = points[-1]
        theta1, theta2 = self._verify_start_end_angles(arcellipse3d, theta1, theta2)
        points[0] = volmdlr.Point2D(theta1, z1)
        points[-1] = volmdlr.Point2D(theta2, z2)

        theta_list = [point.x for point in points]
        theta_discontinuity, indexes_theta_discontinuity = angle_discontinuity(theta_list)
        if theta_discontinuity:
            points = self._fix_angle_discontinuity_on_discretization_points(points,
                                                                            indexes_theta_discontinuity, "x")

        bsplinecurve2d = vme.BSplineCurve2D.from_points_interpolation(points, degree=2, periodic=True, name="ellipse")
        return [bsplinecurve2d]

    def bsplinecurve2d_to_3d(self, bspline_curve2d):
        """
        Is this right?.

        #TODO: today ArcEllipse3D are also represented as bspline curves 2D.
        """
        n = len(bspline_curve2d.control_points)
        points = [self.point2d_to_3d(p)
                  for p in bspline_curve2d.discretization_points(number_points=n)]
        return [vme.BSplineCurve3D.from_points_interpolation(points, bspline_curve2d.degree, bspline_curve2d.periodic)]

    def linesegment2d_to_3d(self, linesegment2d):
        """
        Converts a BREP line segment 2D onto a 3D primitive on the surface.
        """
        theta1, z1 = linesegment2d.start
        theta2, z2 = linesegment2d.end
        if math.isclose(theta1, theta2, abs_tol=1e-4):
            return [vme.LineSegment3D(self.point2d_to_3d(linesegment2d.start),
                                      self.point2d_to_3d(linesegment2d.end))]
        if math.isclose(z1, z2, abs_tol=1e-4):
            if math.isclose(abs(theta1 - theta2), volmdlr.TWO_PI, abs_tol=1e-4):
                return [vme.FullArc3D(center=self.frame.origin + z1 * self.frame.w,
                                      start_end=self.point2d_to_3d(linesegment2d.start),
                                      normal=self.frame.w)]

            return [vme.Arc3D(
                self.point2d_to_3d(linesegment2d.start),
                self.point2d_to_3d(volmdlr.Point2D(0.5 * (theta1 + theta2), z1)),
                self.point2d_to_3d(linesegment2d.end)
            )]
        return [vme.LineSegment3D(self.point2d_to_3d(linesegment2d.start), self.point2d_to_3d(linesegment2d.end))]


class CylindricalSurface3D(PeriodicalSurface):
    """
    The local plane is defined by (theta, z).

    :param frame: frame.w is axis, frame.u is theta=0 frame.v theta=pi/2
    :param frame:
    :param radius: Cylinder's radius
    :type radius: float
    """
    face_class = 'CylindricalFace3D'
    x_periodicity = volmdlr.TWO_PI
    y_periodicity = None

    def __init__(self, frame, radius: float, name: str = ''):
        self.frame = frame
        self.radius = radius
        PeriodicalSurface.__init__(self, name=name)

    def plot(self, ax=None, color: str = "grey", alpha: float = 0.5, z: float = 0.5):
        """
        Plot the cylindrical surface in the local frame normal direction.

        :param ax: Matplotlib Axes3D object to plot on. If None, create a new figure.
        :type ax: Axes3D or None
        :param color: color of the wireframe plot. Default is 'grey'.
        :type color: str
        :param alpha: transparency of the edge frame plot. Default is 0.5.
        :type alpha: float
        :param z: additional keyword arguments to pass the value of z to cut the surface.
        :type z: float
        :return: Matplotlib Axes3D object containing the plotted wireframe.
        :rtype: Axes3D
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        x = self.radius
        point1 = self.frame.local_to_global_coordinates(volmdlr.Point3D(x, 0, 0))
        point2 = self.frame.local_to_global_coordinates(volmdlr.Point3D(x, 0, z))
        generatrix = vme.LineSegment3D(point1, point2)
        for i in range(37):
            theta = i / 36. * volmdlr.TWO_PI
            wire = generatrix.rotation(self.frame.origin, self.frame.w, theta)
            wire.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha))
        return ax

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        """
        Coverts a parametric coordinate on the surface into a 3D spatial point (x, y, z).

        :param point2d: Point at the ToroidalSuface3D
        :type point2d: `volmdlr.`Point2D`
        """

        point = volmdlr.Point3D(self.radius * math.cos(point2d.x),
                                self.radius * math.sin(point2d.x),
                                point2d.y)
        return self.frame.local_to_global_coordinates(point)

    def point3d_to_2d(self, point3d):
        """
        Returns the cylindrical coordinates volmdlr.Point2D(theta, z) of a Cartesian coordinates point (x, y, z).

        :param point3d: Point at the CylindricalSuface3D
        :type point3d: `volmdlr.`Point3D`
        """
        x, y, z = self.frame.global_to_local_coordinates(point3d)
        # Do not delete this, mathematical problem when x and y close to zero but not 0
        if abs(x) < 1e-12:
            x = 0
        if abs(y) < 1e-12:
            y = 0

        theta = math.atan2(y, x)
        if abs(theta) < 1e-9:
            theta = 0.0

        return volmdlr.Point2D(theta, z)

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a CylindricalSurface3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated
        :type object_dict: dict
        :return: The corresponding CylindricalSurface3D object.
        :rtype: :class:`volmdlr.faces.CylindricalSurface3D`
        """

        length_conversion_factor = kwargs.get("length_conversion_factor", 1)
        frame3d = object_dict[arguments[1]]
        u_vector, w_vector = frame3d.v, -frame3d.u
        u_vector.normalize()
        w_vector.normalize()
        v_vector = w_vector.cross(u_vector)
        frame_direct = volmdlr.Frame3D(frame3d.origin, u_vector, v_vector, w_vector)
        radius = float(arguments[2]) * length_conversion_factor
        return cls(frame_direct, radius, arguments[0][1:-1])

    def to_step(self, current_id):
        """
        Translate volmdlr primitive to step syntax.
        """
        frame = volmdlr.Frame3D(self.frame.origin, self.frame.w, self.frame.u,
                                self.frame.v)
        content, frame_id = frame.to_step(current_id)
        current_id = frame_id + 1
        content += f"#{current_id} = CYLINDRICAL_SURFACE('{self.name}',#{frame_id},{round(1000 * self.radius, 3)});\n"
        return content, [current_id]

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new CylindricalSurface3D.

        :param side: 'old' or 'new'
        """
        new_frame = self.frame.frame_mapping(frame, side)
        return CylindricalSurface3D(new_frame, self.radius,
                                    name=self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        :param side: 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_frame = self.frame.frame_mapping(frame, side)
        self.frame = new_frame

    def rectangular_cut(self, theta1: float, theta2: float,
                        z1: float, z2: float, name: str = ''):
        """
        Cut a rectangular piece of the CylindricalSurface3D object and return a CylindricalFace3D object.

        """

        if theta1 == theta2:
            theta2 += volmdlr.TWO_PI

        point1 = volmdlr.Point2D(theta1, z1)
        point2 = volmdlr.Point2D(theta2, z1)
        point3 = volmdlr.Point2D(theta2, z2)
        point4 = volmdlr.Point2D(theta1, z2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        surface2d = Surface2D(outer_contour, [])
        return volmdlr.faces.CylindricalFace3D(self, surface2d, name)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        CylindricalFace3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated Plane3D.
        """
        new_frame = self.frame.rotation(center=center, axis=axis,
                                        angle=angle)
        return CylindricalSurface3D(new_frame, self.radius)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        CylindricalFace3D rotation. Object is updated inplace.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: rotation angle.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.frame.rotation_inplace(center, axis, angle)

    def translation(self, offset: volmdlr.Vector3D):
        """
        CylindricalFace3D translation.

        :param offset: translation vector.
        :return: A new translated CylindricalFace3D.
        """
        return CylindricalSurface3D(self.frame.translation(offset), self.radius)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        CylindricalFace3D translation. Object is updated inplace.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.frame.translation_inplace(offset)

    def grid3d(self, grid2d: volmdlr.grid.Grid2D):
        """
        Generate 3d grid points of a Cylindrical surface, based on a Grid2D.

        """

        points_2d = grid2d.points
        points_3d = [self.point2d_to_3d(point2d) for point2d in points_2d]

        return points_3d

    def line_intersections(self, line: vme.Line3D):
        line_2d = line.to_2d(self.frame.origin, self.frame.u, self.frame.v)
        if line_2d is None:
            return []
        origin2d = self.frame.origin.to_2d(self.frame.origin, self.frame.u, self.frame.v)
        distance_line2d_to_origin = line_2d.point_distance(origin2d)
        if distance_line2d_to_origin > self.radius:
            return []
        a_prime = line_2d.point1
        b_prime = line_2d.point2
        a_prime_minus_b_prime = a_prime - b_prime
        t_param = a_prime.dot(a_prime_minus_b_prime) / a_prime_minus_b_prime.dot(a_prime_minus_b_prime)
        k_param = math.sqrt(
            (self.radius ** 2 - distance_line2d_to_origin ** 2) / a_prime_minus_b_prime.dot(a_prime_minus_b_prime))
        intersection1 = line.point1 + (t_param + k_param) * (line.direction_vector())
        intersection2 = line.point1 + (t_param - k_param) * (line.direction_vector())
        if intersection1 == intersection2:
            return [intersection1]

        return [intersection1, intersection2]

    def linesegment_intersections(self, linesegment: vme.LineSegment3D):
        line = linesegment.to_line()
        line_intersections = self.line_intersections(line)
        linesegment_intersections = [inters for inters in line_intersections if linesegment.point_belongs(inters)]
        return linesegment_intersections

    def parallel_plane_intersection(self, plane3d):
        """
        Cylinder plane intersections when plane's normal is perpendicular with the cylinder axis.

        :param plane3d: intersecting plane
        :return: list of intersecting curves
        """
        distance_plane_cylinder_axis = plane3d.point_distance(self.frame.origin)
        if distance_plane_cylinder_axis > self.radius:
            return []
        if math.isclose(self.frame.w.dot(plane3d.frame.u), 0, abs_tol=1e-6):
            line = vme.Line3D(plane3d.frame.origin, plane3d.frame.origin + plane3d.frame.u)
        else:
            line = vme.Line3D(plane3d.frame.origin, plane3d.frame.origin + plane3d.frame.v)
        line_intersections = self.line_intersections(line)
        lines = []
        for intersection in line_intersections:
            lines.append(vme.Line3D(intersection, intersection + self.frame.w))
        return lines

    def perpendicular_plane_intersection(self, plane3d):
        """
        Cylinder plane intersections when plane's normal is parallel with the cylinder axis.

        :param plane3d: intersecting plane
        :return: list of intersecting curves
        """
        line = vme.Line3D(self.frame.origin, self.frame.origin + self.frame.w)
        center3d_plane = plane3d.line_intersections(line)[0]
        circle3d = volmdlr.wires.Circle3D(volmdlr.Frame3D(center3d_plane, plane3d.frame.u,
                                                          plane3d.frame.v, plane3d.frame.w), self.radius)
        return [circle3d]

    def concurrent_plane_intersection(self, plane3d):
        """
        Cylinder plane intersections when plane's normal is concurrent with the cylinder axis, but not orthogonal.

        Ellipse vector equation : < r*cos(t), r*sin(t), -(1 / c)*(d + a*r*cos(t) +
        b*r*sint(t)); d = - (ax_0 + by_0 + cz_0).

        :param plane3d: intersecting plane.
        :return: list of intersecting curves.
        """
        line = vme.Line3D(self.frame.origin, self.frame.origin + self.frame.w)
        center3d_plane = plane3d.line_intersections(line)[0]
        plane_coefficient_a, plane_coefficient_b, plane_coefficient_c, plane_coefficient_d = \
            plane3d.equation_coefficients()
        ellipse_0 = volmdlr.Point3D(
            self.radius * math.cos(0),
            self.radius * math.sin(0),
            - (1 / plane_coefficient_c) * (plane_coefficient_d + plane_coefficient_a * self.radius * math.cos(0) +
                                           plane_coefficient_b * self.radius * math.sin(0)))
        ellipse_pi_by_2 = volmdlr.Point3D(
            self.radius * math.cos(math.pi / 2),
            self.radius * math.sin(math.pi / 2),
            - (1 / plane_coefficient_c) * (
                    plane_coefficient_d + plane_coefficient_a * self.radius * math.cos(math.pi / 2)
                    + plane_coefficient_b * self.radius * math.sin(math.pi / 2)))
        axis_1 = center3d_plane.point_distance(ellipse_0)
        axis_2 = center3d_plane.point_distance(ellipse_pi_by_2)
        if axis_1 > axis_2:
            major_axis = axis_1
            minor_axis = axis_2
            major_dir = ellipse_0 - center3d_plane
        else:
            major_axis = axis_2
            minor_axis = axis_1
            major_dir = ellipse_pi_by_2 - center3d_plane
        return [volmdlr.wires.Ellipse3D(major_axis, minor_axis, center3d_plane, plane3d.frame.w, major_dir)]

    def plane_intersection(self, plane3d):
        """
        Cylinder intersections with a plane.

        :param plane3d: intersecting plane.
        :return: list of intersecting curves.
        """
        if math.isclose(abs(plane3d.frame.w.dot(self.frame.w)), 0, abs_tol=1e-6):
            return self.parallel_plane_intersection(plane3d)
        if math.isclose(abs(plane3d.frame.w.dot(self.frame.w)), 1, abs_tol=1e-6):
            return self.perpendicular_plane_intersection(plane3d)
        return self.concurrent_plane_intersection(plane3d)

    def is_coincident(self, surface3d):
        """
        Verifies if two CylindricalSurfaces are coincident.

        :param surface3d: surface to verify.
        :return: True if they are coincident, False otherwise.
        """
        if not isinstance(self, surface3d.__class__):
            return False
        if math.isclose(abs(self.frame.w.dot(surface3d.frame.w)), 1.0, abs_tol=1e-6) and \
                self.radius == surface3d.radius:
            return True
        return False

    def point_on_surface(self, point3d):
        """
        Verifies if a given point is on the CylindricalSurface3D.

        :param point3d: point to verify.
        :return: True if point on surface, False otherwise.
        """
        new_point = self.frame.global_to_local_coordinates(point3d)
        if math.isclose(new_point.x ** 2 + new_point.y ** 2, self.radius ** 2, abs_tol=1e-6):
            return True
        return False


class ToroidalSurface3D(PeriodicalSurface):
    """
    The local plane is defined by (theta, phi).

    Theta is the angle around the big (R) circle and phi around the small (r).

    :param frame: Tore's frame: origin is the center, u is pointing at theta=0
    :param tore_radius: Tore's radius
    :param r: Circle to revolute radius

    :See Also:

    Definitions of R and r according to https://en.wikipedia.org/wiki/Torus
    """
    face_class = 'ToroidalFace3D'
    x_periodicity = volmdlr.TWO_PI
    y_periodicity = volmdlr.TWO_PI

    def __init__(self, frame: volmdlr.Frame3D, tore_radius: float, small_radius: float, name: str = ''):
        self.frame = frame
        self.tore_radius = tore_radius
        self.small_radius = small_radius
        PeriodicalSurface.__init__(self, name=name)

        self._bbox = None

    @property
    def bounding_box(self):
        """
        Returns the surface bounding box.
        """
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    def _bounding_box(self):
        distance = self.tore_radius + self.small_radius
        point1 = self.frame.origin + \
                 self.frame.u * distance + self.frame.v * distance + self.frame.w * self.small_radius
        point2 = self.frame.origin + \
                 self.frame.u * distance + self.frame.v * distance - self.frame.w * self.small_radius
        point3 = self.frame.origin + \
                 self.frame.u * distance - self.frame.v * distance + self.frame.w * self.small_radius
        point4 = self.frame.origin + \
                 self.frame.u * distance - self.frame.v * distance - self.frame.w * self.small_radius
        point5 = self.frame.origin - \
                 self.frame.u * distance + self.frame.v * distance + self.frame.w * self.small_radius
        point6 = self.frame.origin - \
                 self.frame.u * distance + self.frame.v * distance - self.frame.w * self.small_radius
        point7 = self.frame.origin - \
                 self.frame.u * distance - self.frame.v * distance + self.frame.w * self.small_radius
        point8 = self.frame.origin - \
                 self.frame.u * distance - self.frame.v * distance - self.frame.w * self.small_radius

        return volmdlr.core.BoundingBox.from_points(
            [point1, point2, point3, point4, point5, point6, point7, point8])

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        """
        Coverts a parametric coordinate on the surface into a 3D spatial point (x, y, z).

        :param point2d: Point at the ToroidalSuface3D
        :type point2d: `volmdlr.`Point2D`
        """
        theta, phi = point2d
        x = (self.tore_radius + self.small_radius * math.cos(phi)) * math.cos(theta)
        y = (self.tore_radius + self.small_radius * math.cos(phi)) * math.sin(theta)
        z = self.small_radius * math.sin(phi)
        return self.frame.local_to_global_coordinates(volmdlr.Point3D(x, y, z))

    def point3d_to_2d(self, point3d):
        """
        Transform a 3D spatial point (x, y, z) into a 2D spherical parametric point (theta, phi).
        """
        x, y, z = self.frame.global_to_local_coordinates(point3d)
        z = min(self.small_radius, max(-self.small_radius, z))

        # Do not delete this, mathematical problem when x and y close to zero (should be zero) but not 0
        # Genarally this is related to uncertaintity of step files.

        if abs(x) < 1e-12:
            x = 0
        if abs(y) < 1e-12:
            y = 0

        zr = z / self.small_radius
        phi = math.asin(zr)
        if abs(phi) < 1e-9:
            phi = 0

        u = self.tore_radius + math.sqrt((self.small_radius ** 2) - (z ** 2))
        u1, u2 = round(x / u, 5), round(y / u, 5)
        theta = math.atan2(u2, u1)

        vector_to_tube_center = volmdlr.Vector3D(self.tore_radius * math.cos(theta),
                                                 self.tore_radius * math.sin(theta), 0)
        vector_from_tube_center_to_point = volmdlr.Vector3D(x, y, z) - vector_to_tube_center
        phi2 = volmdlr.geometry.vectors3d_angle(vector_to_tube_center, vector_from_tube_center_to_point)

        if phi >= 0 and phi2 > 0.5 * math.pi:
            phi = math.pi - phi
        elif phi < 0 and phi2 > 0.5 * math.pi:
            phi = -math.pi - phi
        if abs(theta) < 1e-9:
            theta = 0.0
        if abs(phi) < 1e-9:
            phi = 0.0
        return volmdlr.Point2D(theta, phi)

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a ToroidalSurface3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding ToroidalSurface3D object.
        :rtype: :class:`volmdlr.faces.ToroidalSurface3D`
        """

        length_conversion_factor = kwargs.get("length_conversion_factor", 1)

        frame3d = object_dict[arguments[1]]
        u_vector, w_vector = frame3d.v, -frame3d.u
        u_vector.normalize()
        w_vector.normalize()
        v_vector = w_vector.cross(u_vector)
        frame_direct = volmdlr.Frame3D(frame3d.origin, u_vector, v_vector, w_vector)
        rcenter = float(arguments[2]) * length_conversion_factor
        rcircle = float(arguments[3]) * length_conversion_factor
        return cls(frame_direct, rcenter, rcircle, arguments[0][1:-1])

    def to_step(self, current_id):
        frame = volmdlr.Frame3D(self.frame.origin, self.frame.w, self.frame.u,
                                self.frame.v)
        content, frame_id = frame.to_step(current_id)
        current_id = frame_id + 1
        content += f"#{current_id} = TOROIDAL_SURFACE('{self.name}',#{frame_id}," \
                   f"{round(1000 * self.tore_radius, 3)},{round(1000 * self.small_radius, 3)});\n"
        return content, [current_id]

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new ToroidalSurface3D.

        :param frame: The new frame to map to.
        :type frame: `volmdlr.Frame3D
        :param side: Indicates whether the frame should be mapped to the 'old' or 'new' frame.
            Acceptable values are 'old' or 'new'.
        :type side: str
        """
        new_frame = self.frame.frame_mapping(frame, side)
        return ToroidalSurface3D(new_frame, self.tore_radius, self.small_radius, name=self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        :param frame: The new frame to map to.
        :type frame: `volmdlr.Frame3D
        :param side: Indicates whether the frame should be mapped to the 'old' or 'new' frame.
            Acceptable values are 'old' or 'new'.
        :type side: str
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_frame = self.frame.frame_mapping(frame, side)
        self.frame = new_frame

    def rectangular_cut(self, theta1: float, theta2: float, phi1: float, phi2: float, name: str = ""):
        """
        Cut a rectangular piece of the ToroidalSurface3D object and return a ToroidalFace3D object.

        :param theta1: Start angle of the cut in theta direction.
        :param theta2: End angle of the cut in theta direction.
        :param phi1: Start angle of the cut in phi direction.
        :param phi2: End angle of the cut in phi direction.
        :param name: (optional) Name of the returned ToroidalFace3D object. Defaults to "".
        :return: A ToroidalFace3D object created by cutting the ToroidalSurface3D object.
        :rtype: ToroidalFace3D
        """
        if phi1 == phi2:
            phi2 += volmdlr.TWO_PI
        elif phi2 < phi1:
            phi2 += volmdlr.TWO_PI
        if theta1 == theta2:
            theta2 += volmdlr.TWO_PI
        elif theta2 < theta1:
            theta2 += volmdlr.TWO_PI

        point1 = volmdlr.Point2D(theta1, phi1)
        point2 = volmdlr.Point2D(theta2, phi1)
        point3 = volmdlr.Point2D(theta2, phi2)
        point4 = volmdlr.Point2D(theta1, phi2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        return ToroidalFace3D(self,
                              Surface2D(outer_contour, []),
                              name)

    def linesegment2d_to_3d(self, linesegment2d):
        """
        Converts the parametric boundary representation into a 3D primitive.
        """
        theta1, phi1 = linesegment2d.start
        theta2, phi2 = linesegment2d.end
        if math.isclose(theta1, theta2, abs_tol=1e-4):
            if math.isclose(abs(phi1 - phi2), volmdlr.TWO_PI, abs_tol=1e-4):
                u_vector = self.frame.u.rotation(self.frame.origin, self.frame.w, angle=theta1)
                v_vector = self.frame.u.rotation(self.frame.origin, self.frame.w, angle=theta1)
                center = self.frame.origin + self.tore_radius * u_vector
                return [vme.FullArc3D(center=center,
                                      start_end=center + self.small_radius * u_vector,
                                      normal=v_vector)]
            return [vme.Arc3D(
                self.point2d_to_3d(linesegment2d.start),
                self.point2d_to_3d(volmdlr.Point2D(theta1, 0.5 * (phi1 + phi2))),
                self.point2d_to_3d(linesegment2d.end),
            )]
        if math.isclose(phi1, phi2, abs_tol=1e-4):
            if math.isclose(abs(theta1 - theta2), volmdlr.TWO_PI, abs_tol=1e-4):
                center = self.frame.origin + self.small_radius * math.sin(phi1) * self.frame.w
                start_end = center + self.frame.u * (self.small_radius + self.tore_radius)
                return [vme.FullArc3D(center=center,
                                      start_end=start_end,
                                      normal=self.frame.w)]
            return [vme.Arc3D(
                self.point2d_to_3d(linesegment2d.start),
                self.point2d_to_3d(volmdlr.Point2D(0.5 * (theta1 + theta2), phi1)),
                self.point2d_to_3d(linesegment2d.end),
            )]
        raise NotImplementedError('Ellipse?')

    def fullarc3d_to_2d(self, fullarc3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        start = self.point3d_to_2d(fullarc3d.start)
        end = self.point3d_to_2d(fullarc3d.end)

        length = fullarc3d.length()
        angle3d = fullarc3d.angle
        point_after_start = self.point3d_to_2d(fullarc3d.point_at_abscissa(0.001 * length))
        point_before_end = self.point3d_to_2d(fullarc3d.point_at_abscissa(0.98 * length))

        start, end = vm_parametric.arc3d_to_spherical_coordinates_verification(start, end, angle3d,
                                                                               [point_after_start, point_before_end],
                                                                               [self.x_periodicity,
                                                                                self.y_periodicity])
        theta1, phi1 = start
        # theta2, phi2 = end
        theta3, phi3 = point_after_start
        # theta4, phi4 = point_before_end
        if self.frame.w.is_colinear_to(fullarc3d.normal, abs_tol=1e-4):
            point1 = start
            if theta1 > theta3:
                point2 = volmdlr.Point2D(theta1 - volmdlr.TWO_PI, phi1)
            elif theta1 < theta3:
                point2 = volmdlr.Point2D(theta1 + volmdlr.TWO_PI, phi1)
            return [vme.LineSegment2D(point1, point2)]
        point1 = start
        if phi1 > phi3:
            point2 = volmdlr.Point2D(theta1, phi1 - volmdlr.TWO_PI)
        elif phi1 < phi3:
            point2 = volmdlr.Point2D(theta1, phi1 + volmdlr.TWO_PI)
        return [vme.LineSegment2D(point1, point2)]

    def arc3d_to_2d(self, arc3d):
        start = self.point3d_to_2d(arc3d.start)
        end = self.point3d_to_2d(arc3d.end)

        length = arc3d.length()
        angle3d = arc3d.angle
        point_after_start = self.point3d_to_2d(arc3d.point_at_abscissa(0.001 * length))
        point_before_end = self.point3d_to_2d(arc3d.point_at_abscissa(0.98 * length))

        start, end = vm_parametric.arc3d_to_spherical_coordinates_verification(start, end, angle3d,
                                                                               [point_after_start, point_before_end],
                                                                               [self.x_periodicity,
                                                                                self.y_periodicity])

        return [vme.LineSegment2D(start, end)]

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        length = bspline_curve3d.length()
        theta3, phi3 = self.point3d_to_2d(bspline_curve3d.point_at_abscissa(0.001 * length))
        theta4, phi4 = self.point3d_to_2d(bspline_curve3d.point_at_abscissa(0.98 * length))
        n = len(bspline_curve3d.control_points)
        points3d = bspline_curve3d.discretization_points(number_points=n)
        points = [self.point3d_to_2d(p) for p in points3d]
        theta1, phi1 = points[0]
        theta2, phi2 = points[-1]

        # Verify if theta1 or theta2 point should be -pi because atan2() -> ]-pi, pi]
        if abs(theta1) == math.pi:
            theta1 = repair_start_end_angle_periodicity(theta1, theta3)
        if abs(theta2) == math.pi:
            theta2 = repair_start_end_angle_periodicity(theta2, theta4)

        # Verify if phi1 or phi2 point should be -pi because phi -> ]-pi, pi]
        if abs(phi1) == math.pi:
            phi1 = repair_start_end_angle_periodicity(phi1, phi3)
        if abs(phi2) == math.pi:
            phi2 = repair_start_end_angle_periodicity(phi2, phi4)

        points[0] = volmdlr.Point2D(theta1, phi1)
        points[-1] = volmdlr.Point2D(theta2, phi2)

        theta_list = [point.x for point in points]
        phi_list = [point.y for point in points]
        theta_discontinuity, indexes_theta_discontinuity = angle_discontinuity(theta_list)
        phi_discontinuity, indexes_phi_discontinuity = angle_discontinuity(phi_list)

        if theta_discontinuity:
            points = self._fix_angle_discontinuity_on_discretization_points(points,
                                                                            indexes_theta_discontinuity, "x")
        if phi_discontinuity:
            points = self._fix_angle_discontinuity_on_discretization_points(points,
                                                                            indexes_phi_discontinuity, "y")

        return [vme.BSplineCurve2D.from_points_interpolation(points, bspline_curve3d.degree, bspline_curve3d.periodic)]

    def arcellipse3d_to_2d(self, arcellipse3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        points = [self.point3d_to_2d(p) for p in arcellipse3d.discretization_points(number_points=15)]
        theta1, phi1 = self.point3d_to_2d(arcellipse3d.start)
        theta2, phi2 = self.point3d_to_2d(arcellipse3d.end)
        # TODO: create a method point_at_abscissa abscissa for ArcEllipse3D and enhance this code

        theta3, phi3 = points[1]
        theta4, phi4 = points[-2]

        # Verify if theta1 or theta2 point should be -pi because atan2() -> ]-pi, pi]
        if abs(theta1) == math.pi:
            theta1 = vm_parametric.repair_start_end_angle_periodicity(theta1, theta3)
        if abs(theta2) == math.pi:
            theta2 = vm_parametric.repair_start_end_angle_periodicity(theta2, theta4)

        # Verify if phi1 or phi2 point should be -pi because phi -> ]-pi, pi]
        if abs(phi1) == math.pi:
            phi1 = vm_parametric.repair_start_end_angle_periodicity(phi1, phi3)
        if abs(phi2) == math.pi:
            phi2 = vm_parametric.repair_start_end_angle_periodicity(phi2, phi4)

        points[0] = volmdlr.Point2D(theta1, phi1)
        points[-1] = volmdlr.Point2D(theta2, phi2)

        theta_list = [point.x for point in points]
        phi_list = [point.y for point in points]
        theta_discontinuity, indexes_theta_discontinuity = angle_discontinuity(theta_list)
        phi_discontinuity, indexes_phi_discontinuity = angle_discontinuity(phi_list)

        if theta_discontinuity:
            points = self._fix_angle_discontinuity_on_discretization_points(points,
                                                                            indexes_theta_discontinuity, "x")
        if phi_discontinuity:
            points = self._fix_angle_discontinuity_on_discretization_points(points,
                                                                            indexes_phi_discontinuity, "y")

        return [vme.BSplineCurve2D.from_points_interpolation(points=points, degree=3, periodic=True)]

    def triangulation(self):
        """
        Triangulation.

        :rtype: vmd.DisplayMesh3D
        """
        face = self.rectangular_cut(0, volmdlr.TWO_PI, 0, volmdlr.TWO_PI)
        return face.triangulation()

    def translation(self, offset: volmdlr.Vector3D):
        """
        ToroidalSurface3D translation.

        :param offset: translation vector
        :return: A new translated ToroidalSurface3D
        """
        return ToroidalSurface3D(self.frame.translation(
            offset), self.tore_radius, self.small_radius)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        ToroidalSurface3D translation. Object is updated inplace.

        :param offset: translation vector.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.frame.translation_inplace(offset)

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        ToroidalSurface3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated ToroidalSurface3D.
        """
        new_frame = self.frame.rotation(center=center, axis=axis,
                                        angle=angle)
        return self.__class__(new_frame, self.tore_radius, self.small_radius)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        ToroidalSurface3D rotation. Object is updated inplace.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: rotation angle.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.frame.rotation_inplace(center, axis, angle)


class ConicalSurface3D(PeriodicalSurface):
    """
    The local plane is defined by (theta, z).

    :param frame: Cone's frame to position it: frame.w is axis of cone frame.origin is at the angle of the cone.
    :param semi_angle: cone's semi-angle.
    """
    face_class = 'ConicalFace3D'
    x_periodicity = volmdlr.TWO_PI
    y_periodicity = None

    def __init__(self, frame: volmdlr.Frame3D, semi_angle: float,
                 name: str = ''):
        self.frame = frame
        self.semi_angle = semi_angle
        PeriodicalSurface.__init__(self, name=name)

    def plot(self, ax=None, color='grey', alpha=0.5, **kwargs):
        z = kwargs.get("z", 0.5)
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        x = z / math.sin(self.semi_angle)
        # point1 = self.frame.local_to_global_coordinates(volmdlr.Point3D(-x, 0, -z))
        point1 = self.frame.origin
        point2 = self.frame.local_to_global_coordinates(volmdlr.Point3D(x, 0, z))
        generatrix = vme.LineSegment3D(point1, point2)
        for i in range(37):
            theta = i / 36. * volmdlr.TWO_PI
            wire = generatrix.rotation(self.frame.origin, self.frame.w, theta)
            wire.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha))
        return ax

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a ConicalSurface3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding ConicalSurface3D object.
        :rtype: :class:`volmdlr.faces.ConicalSurface3D`
        """

        length_conversion_factor = kwargs.get("length_conversion_factor", 1)
        angle_conversion_factor = kwargs.get("angle_conversion_factor", 1)

        frame3d = object_dict[arguments[1]]
        u, w = frame3d.v, frame3d.u
        u.normalize()
        w.normalize()
        v = w.cross(u)
        radius = float(arguments[2]) * length_conversion_factor
        semi_angle = float(arguments[3]) * angle_conversion_factor
        origin = frame3d.origin - radius / math.tan(semi_angle) * w
        frame_direct = volmdlr.Frame3D(origin, u, v, w)
        return cls(frame_direct, semi_angle, arguments[0][1:-1])

    def to_step(self, current_id):
        frame = volmdlr.Frame3D(self.frame.origin, self.frame.w, self.frame.u,
                                self.frame.v)
        content, frame_id = frame.to_step(current_id)
        current_id = frame_id + 1
        content += f"#{current_id} = CONICAL_SURFACE('{self.name}',#{frame_id},{0.},{round(self.semi_angle, 3)});\n"
        return content, [current_id]

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new ConicalSurface3D.

        :param side: 'old' or 'new'
        """
        new_frame = self.frame.frame_mapping(frame, side)
        return ConicalSurface3D(new_frame, self.semi_angle, name=self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        :param side:'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_frame = self.frame.frame_mapping(frame, side)
        self.frame = new_frame

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        """
        Coverts a parametric coordinate on the surface into a 3D spatial point (x, y, z).

        :param point2d: Point at the ConicalSuface3D
        :type point2d: `volmdlr.`Point2D`
        """
        theta, z = point2d
        r = math.tan(self.semi_angle) * z
        new_point = volmdlr.Point3D(r * math.cos(theta),
                                    r * math.sin(theta),
                                    z)
        return self.frame.local_to_global_coordinates(new_point)

    def point3d_to_2d(self, point3d: volmdlr.Point3D):
        """
        Returns the cylindrical coordinates volmdlr.Point2D(theta, z) of a Cartesian coordinates point (x, y, z).

        :param point3d: Point at the CylindricalSuface3D.
        :type point3d: :class:`volmdlr.`Point3D`
        """
        x, y, z = self.frame.global_to_local_coordinates(point3d)
        # Do not delete this, mathematical problem when x and y close to zero (should be zero) but not 0
        # Genarally this is related to uncertainty of step files.
        if abs(x) < 1e-12:
            x = 0
        if abs(y) < 1e-12:
            y = 0
        theta = math.atan2(y, x)
        if abs(theta) < 1e-9:
            theta = 0.0
        return volmdlr.Point2D(theta, z)

    def rectangular_cut(self, theta1: float, theta2: float,
                        z1: float, z2: float, name: str = ''):
        """
        Cut a rectangular piece of the ConicalSurface3D object and return a ConicalFace3D object.

        """
        # theta1 = angle_principal_measure(theta1)
        # theta2 = angle_principal_measure(theta2)
        if theta1 == theta2:
            theta2 += volmdlr.TWO_PI

        point1 = volmdlr.Point2D(theta1, z1)
        point2 = volmdlr.Point2D(theta2, z1)
        point3 = volmdlr.Point2D(theta2, z2)
        point4 = volmdlr.Point2D(theta1, z2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        return ConicalFace3D(self, Surface2D(outer_contour, []), name)

    def linesegment3d_to_2d(self, linesegment3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        start = self.point3d_to_2d(linesegment3d.start)
        end = self.point3d_to_2d(linesegment3d.end)
        if start.x != end.x and start.is_close(volmdlr.Point2D(0, 0)):
            start = volmdlr.Point2D(end.x, 0)
        elif start.x != end.x and end == volmdlr.Point2D(0, 0):
            end = volmdlr.Point2D(start.x, 0)

        if not start.is_close(end):
            return [vme.LineSegment2D(start, end)]
        return [vme.BSplineCurve2D.from_points_interpolation([start, end], 1, False)]

    def linesegment2d_to_3d(self, linesegment2d):
        if linesegment2d.name == "construction":
            return None
        theta1, z1 = linesegment2d.start
        theta2, z2 = linesegment2d.end

        if math.isclose(theta1, theta2, abs_tol=1e-4):
            return [vme.LineSegment3D(
                self.point2d_to_3d(linesegment2d.start),
                self.point2d_to_3d(linesegment2d.end),
            )]
        if math.isclose(z1, z2, abs_tol=1e-4) and math.isclose(z1, 0.,
                                                               abs_tol=1e-6):
            return []
        if math.isclose(z1, z2, abs_tol=1e-4):
            if abs(theta1 - theta2) == volmdlr.TWO_PI:
                return [vme.FullArc3D(center=self.frame.origin + z1 * self.frame.w,
                                      start_end=self.point2d_to_3d(linesegment2d.start),
                                      normal=self.frame.w)]

            return [vme.Arc3D(
                self.point2d_to_3d(linesegment2d.start),
                self.point2d_to_3d(
                    volmdlr.Point2D(0.5 * (theta1 + theta2), z1)),
                self.point2d_to_3d(linesegment2d.end))
            ]
        raise NotImplementedError('Ellipse?')

    def contour3d_to_2d(self, contour3d):
        """
        Transforms a Contour3D into a Contour2D in the parametric domain of the surface.

        :param contour3d: The contour to be transformed.
        :type contour3d: :class:`volmdlr.wires.Contour3D`
        :return: A 2D contour object.
        :rtype: :class:`volmdlr.wires.Contour2D`
        """
        primitives2d = self.primitives3d_to_2d(contour3d.primitives)

        wire2d = volmdlr.wires.Wire2D(primitives2d)
        delta_x = abs(wire2d.primitives[0].start.x - wire2d.primitives[-1].end.x)
        if math.isclose(delta_x, volmdlr.TWO_PI, abs_tol=1e-3) and wire2d.is_ordered():
            if len(primitives2d) > 1:
                # very specific conical case due to the singularity in the point z = 0 on parametric domain.
                if primitives2d[-2].start.y == 0.0:
                    primitives2d = self.repair_primitives_periodicity(primitives2d)
            return volmdlr.wires.Contour2D(primitives2d)
        # Fix contour
        primitives2d = self.repair_primitives_periodicity(primitives2d)
        return volmdlr.wires.Contour2D(primitives2d)

    def translation(self, offset: volmdlr.Vector3D):
        """
        ConicalSurface3D translation.

        :param offset: translation vector.
        :return: A new translated ConicalSurface3D.
        """
        return self.__class__(self.frame.translation(offset),
                              self.semi_angle)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        ConicalSurface3D translation. Object is updated inplace.

        :param offset: translation vector.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.frame.translation_inplace(offset)

    def rotation(self, center: volmdlr.Point3D,
                 axis: volmdlr.Vector3D, angle: float):
        """
        ConicalSurface3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated ConicalSurface3D.
        """
        new_frame = self.frame.rotation(center=center, axis=axis, angle=angle)
        return self.__class__(new_frame, self.semi_angle)

    def rotation_inplace(self, center: volmdlr.Point3D,
                         axis: volmdlr.Vector3D, angle: float):
        """
        ConicalSurface3D rotation. Object is updated inplace.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: rotation angle.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.frame.rotation_inplace(center, axis, angle)

    def repair_primitives_periodicity(self, primitives2d):
        """
        Repairs the continuity of the 2D contour while using contour3d_to_2d on periodic surfaces.

        :param primitives2d: The primitives in parametric surface domain.
        :type primitives2d: list
        :return: A list of primitives.
        :rtype: list
        """
        # Search for a primitive that can be used as reference for repairing periodicity
        pos = vm_parametric.find_index_defined_brep_primitive_on_periodical_surface(primitives2d,
                                                                                    [self.x_periodicity,
                                                                                     self.y_periodicity])
        if pos != 0:
            primitives2d = primitives2d[pos:] + primitives2d[:pos]

        i = 1
        while i < len(primitives2d):
            previous_primitive = primitives2d[i - 1]
            delta = previous_primitive.end - primitives2d[i].start
            if not math.isclose(delta.norm(), 0, abs_tol=1e-5):
                if primitives2d[i].end.is_close(primitives2d[i - 1].end, tol=1e-4) and \
                        math.isclose(primitives2d[i].length(), volmdlr.TWO_PI, abs_tol=1e-4):
                    primitives2d[i] = primitives2d[i].reverse()
                elif delta.norm() and math.isclose(abs(previous_primitive.end.y), 0, abs_tol=1e-6):
                    primitives2d.insert(i, vme.LineSegment2D(previous_primitive.end, primitives2d[i].start,
                                                             name="construction"))
                    i += 1
                else:
                    primitives2d[i] = primitives2d[i].translation(delta)
            elif math.isclose(primitives2d[i].start.y, 0.0, abs_tol=1e-6):
                if primitives2d[i + 1].end.x < primitives2d[i].end.x:
                    primitives2d[i] = primitives2d[i].translation(volmdlr.Vector2D(volmdlr.TWO_PI, 0))
                    primitives2d.insert(i - 1, vme.LineSegment2D(previous_primitive.end, primitives2d[i].start,
                                                                 name="construction"))
                else:
                    primitives2d[i] = primitives2d[i].translation(volmdlr.Vector2D(-volmdlr.TWO_PI, 0))
                    primitives2d.insert(i - 1, vme.LineSegment2D(previous_primitive.end, primitives2d[i].start,
                                                                 name="construction"))
            i += 1
        if not primitives2d[0].start.is_close(primitives2d[-1].end) \
                and primitives2d[0].start.y == 0.0 and primitives2d[-1].end.y == 0.0:
            primitives2d.append(vme.LineSegment2D(primitives2d[-1].end, primitives2d[0].start))

        return primitives2d


class SphericalSurface3D(PeriodicalSurface):
    """
    Defines a spherical surface.

    :param frame: Sphere's frame to position it
    :type frame: volmdlr.Frame3D
    :param radius: Sphere's radius
    :type radius: float
    """
    face_class = 'SphericalFace3D'
    x_periodicity = volmdlr.TWO_PI
    y_periodicity = math.pi

    def __init__(self, frame, radius, name=''):
        self.frame = frame
        self.radius = radius
        PeriodicalSurface.__init__(self, name=name)

        # Hidden Attributes
        self._bbox = None

    @property
    def bounding_box(self):

        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    def _bounding_box(self):
        points = [self.frame.origin + volmdlr.Point3D(-self.radius,
                                                      -self.radius,
                                                      -self.radius),
                  self.frame.origin + volmdlr.Point3D(self.radius,
                                                      self.radius,
                                                      self.radius),

                  ]
        return volmdlr.core.BoundingBox.from_points(points)

    def contour2d_to_3d(self, contour2d):
        primitives3d = []
        for primitive2d in contour2d.primitives:
            method_name = f'{primitive2d.__class__.__name__.lower()}_to_3d'
            if hasattr(self, method_name):
                try:
                    primitives_list = getattr(self, method_name)(primitive2d)
                    if primitives_list:
                        primitives3d.extend(primitives_list)
                    else:
                        continue
                except NotImplementedError:
                    print(f'Class {self.__class__.__name__} does not implement {method_name}'
                          f'with {primitive2d.__class__.__name__}')
            else:
                raise NotImplementedError(f'Class {self.__class__.__name__} does not implement {method_name}')

        return volmdlr.wires.Contour3D(primitives3d)

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a SphericalSurface3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding SphericalSurface3D object.
        :rtype: :class:`volmdlr.faces.SphericalSurface3D`
        """
        length_conversion_factor = kwargs.get("length_conversion_factor", 1)

        frame3d = object_dict[arguments[1]]
        u_vector, w_vector = frame3d.v, frame3d.u
        u_vector.normalize()
        w_vector.normalize()
        v_vector = w_vector.cross(u_vector)
        frame_direct = volmdlr.Frame3D(frame3d.origin, u_vector, v_vector, w_vector)
        radius = float(arguments[2]) * length_conversion_factor
        return cls(frame_direct, radius, arguments[0][1:-1])

    def point2d_to_3d(self, point2d):
        """
        Coverts a parametric coordinate on the surface into a 3D spatial point (x, y, z).

        source: https://mathcurve.com/surfaces/sphere
        # -pi<theta<pi, -pi/2<phi<pi/2

        :param point2d: Point at the CylindricalSuface3D.
        :type point2d: `volmdlr.`Point2D`
        """
        theta, phi = point2d
        x = self.radius * math.cos(phi) * math.cos(theta)
        y = self.radius * math.cos(phi) * math.sin(theta)
        z = self.radius * math.sin(phi)
        return self.frame.local_to_global_coordinates(volmdlr.Point3D(x, y, z))

    def point3d_to_2d(self, point3d):
        """
        Tansform a 3D spatial point (x, y, z) into a 2D spherical parametric point (theta, phi).
        """
        x, y, z = self.frame.global_to_local_coordinates(point3d)
        z = min(self.radius, max(-self.radius, z))

        if z == -0.0:
            z = 0.0

        # Do not delete this, mathematical problem when x and y close to zero (should be zero) but not 0
        # Generally this is related to uncertainty of step files.
        if abs(x) < 1e-12:
            x = 0
        if abs(y) < 1e-12:
            y = 0

        theta = math.atan2(y, x)
        if abs(theta) < 1e-10:
            theta = 0

        zr = z / self.radius
        phi = math.asin(zr)
        if abs(phi) < 1e-10:
            phi = 0

        return volmdlr.Point2D(theta, phi)

    def linesegment2d_to_3d(self, linesegment2d):
        if linesegment2d.name == "construction":
            return []
        start = self.point2d_to_3d(linesegment2d.start)
        interior = self.point2d_to_3d(0.5 * (linesegment2d.start + linesegment2d.end))
        end = self.point2d_to_3d(linesegment2d.end)
        if start.is_close(end) or linesegment2d.length() == 2 * math.pi:
            u_vector = start - self.frame.origin
            u_vector.normalize()
            v_vector = interior - self.frame.origin
            v_vector.normalize()
            normal = u_vector.cross(v_vector)
            return [vme.FullArc3D(self.frame.origin, start, normal)]
        return [vme.Arc3D(start, interior, end)]

    def contour3d_to_2d(self, contour3d):
        """
        Transforms a Contour3D into a Contour2D in the parametric domain of the surface.

        :param contour3d: The contour to be transformed.
        :type contour3d: :class:`volmdlr.wires.Contour3D`
        :return: A 2D contour object.
        :rtype: :class:`volmdlr.wires.Contour2D`
        """
        primitives2d = []

        # Transform the contour's primitives to parametric domain
        for primitive3d in contour3d.primitives:
            method_name = f'{primitive3d.__class__.__name__.lower()}_to_2d'
            if hasattr(self, method_name):
                primitives = getattr(self, method_name)(primitive3d)

                if primitives is None:
                    continue
                primitives2d.extend(primitives)
            else:
                raise NotImplementedError(
                    f'Class {self.__class__.__name__} does not implement {method_name}')
        # Fix contour
        if self.x_periodicity or self.y_periodicity:
            primitives2d = self.repair_primitives_periodicity(primitives2d)
        return volmdlr.wires.Contour2D(primitives2d)

    def is_lat_long_curve(self, theta_list, phi_list):
        """
        Checks if a curve defined on the sphere is a latitude/longitude curve.

        Returns True if it is, False otherwise.
        """
        # Check if curve is a longitude curve (phi is constant)
        if all(math.isclose(theta, theta_list[0], abs_tol=1e-4) for theta in theta_list[1:]):
            return True
        # Check if curve is a latitude curve (theta is constant)
        if all(math.isclose(phi, phi_list[0], abs_tol=1e-4) for phi in phi_list[1:]):
            return True
        return False

    def arc3d_to_2d(self, arc3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        start = self.point3d_to_2d(arc3d.start)
        end = self.point3d_to_2d(arc3d.end)
        theta_i, phi_i = self.point3d_to_2d(arc3d.interior)
        theta1, phi1 = start
        theta2, phi2 = end

        point_after_start, point_before_end = self._reference_points(arc3d)
        theta3, phi3 = point_after_start
        theta4, _ = point_before_end

        # Fix sphere singularity point
        if math.isclose(abs(phi1), 0.5 * math.pi, abs_tol=1e-5) and theta1 == 0.0 \
                and math.isclose(theta3, theta_i, abs_tol=1e-6) and math.isclose(theta4, theta_i, abs_tol=1e-6):
            theta1 = theta_i
            start = volmdlr.Point2D(theta1, phi1)
        if math.isclose(abs(phi2), 0.5 * math.pi, abs_tol=1e-5) and theta2 == 0.0 \
                and math.isclose(theta3, theta_i, abs_tol=1e-6) and math.isclose(theta4, theta_i, abs_tol=1e-6):
            theta2 = theta_i
            end = volmdlr.Point2D(theta2, phi2)

        start, end = vm_parametric.arc3d_to_spherical_coordinates_verification(start, end, arc3d.angle,
                                                                               [point_after_start, point_before_end],
                                                                               [self.x_periodicity,
                                                                                self.y_periodicity])
        if start == end:  # IS THIS POSSIBLE ?
            return [vme.LineSegment2D(start, start + volmdlr.TWO_PI * volmdlr.X2D)]
        if self.is_lat_long_curve([theta1, theta_i, theta2], [phi1, phi_i, phi2]):
            return [vme.LineSegment2D(start, end)]

        return self.arc3d_to_2d_with_singularity(arc3d, start, end, point_after_start)

    def arc3d_to_2d_with_singularity(self, arc3d, start, end, point_after_start):
        # trying to treat when the arc starts at theta1 passes at the singularity at |phi| = 0.5*math.pi
        # and ends at theta2 = theta1 + math.pi
        theta1, phi1 = start
        theta2, phi2 = end
        theta3, phi3 = point_after_start

        half_pi = 0.5 * math.pi
        point_positive_singularity = self.point2d_to_3d(volmdlr.Point2D(theta1, half_pi))
        point_negative_singularity = self.point2d_to_3d(volmdlr.Point2D(theta1, -half_pi))
        positive_singularity = arc3d.point_belongs(point_positive_singularity, 1e-4)
        negative_singularity = arc3d.point_belongs(point_negative_singularity, 1e-4)
        interior = self.point3d_to_2d(arc3d.interior)
        if positive_singularity and negative_singularity:
            thetai = interior.x
            is_trigo = phi1 < phi3
            if is_trigo and abs(phi1) > half_pi:
                half_pi = 0.5 * math.pi
            elif is_trigo and abs(phi1) < half_pi:
                half_pi = - 0.5 * math.pi
            elif not is_trigo and abs(phi1) > half_pi:
                half_pi = - 0.5 * math.pi
            elif is_trigo and abs(phi1) < half_pi:
                half_pi = 0.5 * math.pi
            return [vme.LineSegment2D(volmdlr.Point2D(theta1, phi1), volmdlr.Point2D(theta1, -half_pi)),
                    vme.LineSegment2D(volmdlr.Point2D(theta1, -half_pi), volmdlr.Point2D(thetai, -half_pi),
                                      name="construction"),
                    vme.LineSegment2D(volmdlr.Point2D(thetai, -half_pi), volmdlr.Point2D(thetai, half_pi)),
                    vme.LineSegment2D(volmdlr.Point2D(thetai, half_pi), volmdlr.Point2D(theta2, half_pi),
                                      name="construction"),
                    vme.LineSegment2D(volmdlr.Point2D(theta2, half_pi), volmdlr.Point2D(theta2, phi2))
                    ]

        if (positive_singularity or negative_singularity) and \
                math.isclose(abs(theta2 - theta1), math.pi, abs_tol=1e-4):
            if abs(phi1) == 0.5 * math.pi:
                return [vme.LineSegment2D(volmdlr.Point2D(theta3, phi1), volmdlr.Point2D(theta2, phi2))]
            if theta1 == math.pi and theta2 != math.pi:
                theta1 = -math.pi
            if theta2 == math.pi and theta1 != math.pi:
                theta2 = -math.pi

            return [vme.LineSegment2D(volmdlr.Point2D(theta1, phi1), volmdlr.Point2D(theta1, half_pi)),
                    vme.LineSegment2D(volmdlr.Point2D(theta1, half_pi), volmdlr.Point2D(theta2, half_pi),
                                      name="construction"),
                    vme.LineSegment2D(volmdlr.Point2D(theta2, half_pi), volmdlr.Point2D(theta2, phi2))
                    ]

        # maybe this is incomplete and not exact
        angle3d = arc3d.angle
        number_points = math.ceil(angle3d * 50) + 1  # 50 points per radian
        number_points = max(number_points, 5)
        points3d = arc3d.discretization_points(number_points=number_points)
        points = [self.point3d_to_2d(p) for p in points3d]

        points[0] = start  # to take into account all the previous verification
        points[-1] = end  # to take into account all the previous verification

        if theta3 < theta1 < theta2:
            points = [p - volmdlr.Point2D(self.x_periodicity, 0) if p.x > 0 else p for p in points]
        elif theta3 > theta1 > theta2:
            points = [p + volmdlr.Point2D(self.x_periodicity, 0) if p.x < 0 else p for p in points]

        return [vme.BSplineCurve2D.from_points_interpolation(points, 2)]

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        length = bspline_curve3d.length()
        n = len(bspline_curve3d.control_points)
        points3d = bspline_curve3d.discretization_points(number_points=n)
        points = [self.point3d_to_2d(point) for point in points3d]
        theta1, phi1 = points[0]
        theta2, phi2 = points[-1]

        theta3, _ = self.point3d_to_2d(bspline_curve3d.point_at_abscissa(0.001 * length))
        # make sure that the reference angle is not undefined
        if abs(theta3) == math.pi:
            theta3, _ = self.point3d_to_2d(bspline_curve3d.point_at_abscissa(0.002 * length))

        # Verify if theta1 or theta2 point should be -pi because atan2() -> ]-pi, pi]
        if abs(theta1) == math.pi:
            theta1 = repair_start_end_angle_periodicity(theta1, theta3)
        if abs(theta2) == math.pi:
            theta4, _ = self.point3d_to_2d(bspline_curve3d.point_at_abscissa(0.98 * length))
            # make sure that the reference angle is not undefined
            if abs(theta4) == math.pi:
                theta4, _ = self.point3d_to_2d(bspline_curve3d.point_at_abscissa(0.97 * length))
            theta2 = repair_start_end_angle_periodicity(theta2, theta4)

        points[0] = volmdlr.Point2D(theta1, phi1)
        points[-1] = volmdlr.Point2D(theta2, phi2)

        theta_list = [point.x for point in points]
        theta_discontinuity, indexes_theta_discontinuity = angle_discontinuity(theta_list)
        if theta_discontinuity:
            points = self._fix_angle_discontinuity_on_discretization_points(points, indexes_theta_discontinuity, "x")

        return [vme.BSplineCurve2D.from_points_interpolation(points, degree=bspline_curve3d.degree,
                                                             periodic=bspline_curve3d.periodic)]

    def bsplinecurve2d_to_3d(self, bspline_curve2d):
        # TODO: this is incomplete, a bspline_curve2d can be also a bspline_curve3d
        i = round(0.5 * len(bspline_curve2d.points))
        start = self.point2d_to_3d(bspline_curve2d.points[0])
        interior = self.point2d_to_3d(bspline_curve2d.points[i])
        end = self.point2d_to_3d(bspline_curve2d.points[-1])
        arc3d = vme.Arc3D(start, interior, end)
        flag = True
        points3d = [self.point2d_to_3d(p) for p in bspline_curve2d.points]
        for point in points3d:
            if not arc3d.point_belongs(point, 1e-4):
                flag = False
                break
        if flag:
            return [arc3d]

        return [vme.BSplineCurve3D.from_points_interpolation(points3d, degree=bspline_curve2d.degree,
                                                             periodic=bspline_curve2d.periodic)]

    def fullarc3d_to_2d(self, fullarc3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        # TODO: On a spherical surface we can have fullarc3d in any plane
        length = fullarc3d.length()

        theta1, phi1 = self.point3d_to_2d(fullarc3d.start)
        theta2, phi2 = self.point3d_to_2d(fullarc3d.end)
        theta3, phi3 = self.point3d_to_2d(fullarc3d.point_at_abscissa(0.001 * length))
        theta4, _ = self.point3d_to_2d(fullarc3d.point_at_abscissa(0.98 * length))

        if self.frame.w.is_colinear_to(fullarc3d.normal):
            point1 = volmdlr.Point2D(theta1, phi1)
            if theta1 > theta3:
                point2 = volmdlr.Point2D(theta1 - volmdlr.TWO_PI, phi2)
            elif theta1 < theta3:
                point2 = volmdlr.Point2D(theta1 + volmdlr.TWO_PI, phi2)
            return [vme.LineSegment2D(point1, point2)]

        if math.isclose(self.frame.w.dot(fullarc3d.normal), 0, abs_tol=1e-4):
            if theta1 > theta3:
                theta_plus_pi = theta1 - math.pi
            else:
                theta_plus_pi = theta1 + math.pi
            if phi1 > phi3:
                half_pi = 0.5 * math.pi
            else:
                half_pi = -0.5 * math.pi
            if abs(phi1) == 0.5 * math.pi:
                return [vme.LineSegment2D(volmdlr.Point2D(theta3, phi1), volmdlr.Point2D(theta3, -half_pi)),
                        vme.LineSegment2D(volmdlr.Point2D(theta4, -half_pi), volmdlr.Point2D(theta4, phi2))]

            return [vme.LineSegment2D(volmdlr.Point2D(theta1, phi1), volmdlr.Point2D(theta1, -half_pi)),
                    vme.LineSegment2D(volmdlr.Point2D(theta_plus_pi, -half_pi),
                                      volmdlr.Point2D(theta_plus_pi, half_pi)),
                    vme.LineSegment2D(volmdlr.Point2D(theta1, half_pi), volmdlr.Point2D(theta1, phi2))]

        points = [self.point3d_to_2d(p) for p in fullarc3d.discretization_points(angle_resolution=25)]

        # Verify if theta1 or theta2 point should be -pi because atan2() -> ]-pi, pi]
        theta1 = vm_parametric.repair_start_end_angle_periodicity(theta1, theta3)
        theta2 = vm_parametric.repair_start_end_angle_periodicity(theta2, theta4)

        points[0] = volmdlr.Point2D(theta1, phi1)
        points[-1] = volmdlr.Point2D(theta2, phi2)

        if theta3 < theta1 < theta2:
            points = [p - volmdlr.Point2D(volmdlr.TWO_PI, 0) if p.x > 0 else p for p in points]
        elif theta3 > theta1 > theta2:
            points = [p + volmdlr.Point2D(volmdlr.TWO_PI, 0) if p.x < 0 else p for p in points]

        return [vme.BSplineCurve2D.from_points_interpolation(points, 2)]

    def plot(self, ax=None, color='grey', alpha=0.5):
        """Plot sphere arcs."""
        for i in range(20):
            theta = i / 20. * volmdlr.TWO_PI
            t_points = []
            for j in range(20):
                phi = j / 20. * volmdlr.TWO_PI
                t_points.append(self.point2d_to_3d(volmdlr.Point2D(theta, phi)))
            ax = volmdlr.wires.ClosedPolygon3D(t_points).plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha))

        return ax

    def rectangular_cut(self, theta1, theta2, phi1, phi2, name=''):
        """
        Cut a rectangular piece of the SphericalSurface3D object and return a SphericalFace3D object.

        """
        if phi1 == phi2:
            phi2 += volmdlr.TWO_PI
        elif phi2 < phi1:
            phi2 += volmdlr.TWO_PI
        if theta1 == theta2:
            theta2 += volmdlr.TWO_PI
        elif theta2 < theta1:
            theta2 += volmdlr.TWO_PI

        point1 = volmdlr.Point2D(theta1, phi1)
        point2 = volmdlr.Point2D(theta2, phi1)
        point3 = volmdlr.Point2D(theta2, phi2)
        point4 = volmdlr.Point2D(theta1, phi2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        return SphericalFace3D(self,
                               Surface2D(outer_contour, []),
                               name=name)

    def triangulation(self):
        face = self.rectangular_cut(0, volmdlr.TWO_PI, -0.5 * math.pi, 0.5 * math.pi)
        return face.triangulation()

    def repair_primitives_periodicity(self, primitives2d):
        """
        Repairs the continuity of the 2D contour while using contour3d_to_2d on periodic surfaces.

        :param primitives2d: The primitives in parametric surface domain.
        :type primitives2d: list
        :return: A list of primitives.
        :rtype: list
        """
        # # Search for a primitive that can be used as reference for reparing periodicity
        # pos = 0
        # x_periodicity = self.x_periodicity
        # y_periodicity = self.y_periodicity
        # if x_periodicity and not y_periodicity:
        #     for i, primitive in enumerate(primitives2d):
        #         start = primitive.start
        #         end = primitive.end
        #         if abs(start.x) != math.pi and end.x != start.x:
        #             pos = i
        #             break
        #     if pos != 0:
        #         primitives2d = primitives2d[pos:] + primitives2d[:pos]

        i = 1
        while i < len(primitives2d):
            previous_primitive = primitives2d[i - 1]
            delta = previous_primitive.end - primitives2d[i].start
            if not math.isclose(delta.norm(), 0, abs_tol=1e-5):
                if primitives2d[i].end == primitives2d[i - 1].end and \
                        primitives2d[i].length() == volmdlr.TWO_PI:
                    primitives2d[i] = primitives2d[i].reverse()
                elif math.isclose(abs(previous_primitive.end.y), 0.5 * math.pi, abs_tol=1e-6):
                    primitives2d.insert(i, vme.LineSegment2D(previous_primitive.end, primitives2d[i].start,
                                                             name="construction"))
                else:
                    primitives2d[i] = primitives2d[i].translation(delta)
            i += 1
        #     return primitives2d
        # primitives2d = repair(primitives2d)
        last_end = primitives2d[-1].end
        first_start = primitives2d[0].start
        if not last_end.is_close(first_start, tol=1e-3):
            last_end_3d = self.point2d_to_3d(last_end)
            first_start_3d = self.point2d_to_3d(first_start)
            if last_end_3d.is_close(first_start_3d, 1e-6) and \
                    not math.isclose(abs(last_end.y), 0.5 * math.pi, abs_tol=1e-5):
                if first_start.x > last_end.x:
                    half_pi = -0.5 * math.pi
                else:
                    half_pi = 0.5 * math.pi
                lines = [vme.LineSegment2D(last_end, volmdlr.Point2D(last_end.x, half_pi), name="construction"),
                         vme.LineSegment2D(volmdlr.Point2D(last_end.x, half_pi),
                                           volmdlr.Point2D(first_start.x, half_pi), name="construction"),
                         vme.LineSegment2D(volmdlr.Point2D(first_start.x, half_pi),
                                           first_start, name="construction")
                         ]
                primitives2d.extend(lines)
            else:
                primitives2d.append(vme.LineSegment2D(last_end, first_start))
        return primitives2d

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D, angle: float):
        """
        Spherical Surface 3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Spherical Surface 3D
        """
        new_frame = self.frame.rotation(center=center, axis=axis, angle=angle)
        return SphericalSurface3D(new_frame, self.radius)

    def translation(self, offset: volmdlr.Vector3D):
        """
        Spherical Surface 3D translation.

        :param offset: translation vector
        :return: A new translated Spherical Surface 3D
        """
        new_frame = self.frame.translation(offset)
        return SphericalSurface3D(new_frame, self.radius)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes Spherical Surface 3D's frame and return a new Spherical Surface 3D.

        :param frame: Frame of reference
        :type frame: `volmdlr.Frame3D`
        :param side: 'old' or 'new'
        """
        new_frame = self.frame.frame_mapping(frame, side)
        return SphericalSurface3D(new_frame, self.radius)


class RuledSurface3D(Surface3D):
    """
    Defines a ruled surface between two wires.

    :param wire1: Wire
    :type wire1: :class:`vmw.Wire3D`
    :param wire2: Wire
    :type wire2: :class:`volmdlr.wires.Wire3D`
    """
    face_class = 'RuledFace3D'

    def __init__(self, wire1: volmdlr.wires.Wire3D, wire2: volmdlr.wires.Wire3D, name: str = ''):
        self.wire1 = wire1
        self.wire2 = wire2
        self.length1 = wire1.length()
        self.length2 = wire2.length()
        Surface3D.__init__(self, name=name)

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
        """
        Cut a rectangular piece of the RuledSurface3D object and return a RuledFace3D object.

        """
        point1 = volmdlr.Point2D(x1, y1)
        point2 = volmdlr.Point2D(x2, y1)
        point3 = volmdlr.Point2D(x2, y2)
        point4 = volmdlr.Point2D(x1, y2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        surface2d = Surface2D(outer_contour, [])
        return volmdlr.faces.RuledFace3D(self, surface2d, name)


class ExtrusionSurface3D(Surface3D):
    """
    Defines a surface of revolution.

    An extrusion surface is a sufarce that is a generic cylindrical surface genarated by the linear
    extrusion of a curve, generally an Ellipse or a B-Spline curve.

    :param edge: edge.
    :type edge: Union[:class:`vmw.Wire3D`, :class:`vmw.Contour3D`]
    :param axis_point: Axis placement
    :type axis_point: :class:`volmdlr.Point3D`
    :param axis: Axis of revolution
    :type axis: :class:`volmdlr.Vector3D`
    """
    face_class = 'ExtrusionFace3D'
    y_periodicity = None

    def __init__(self, edge: Union[volmdlr.edges.FullArcEllipse3D, volmdlr.edges.BSplineCurve3D],
                 direction: volmdlr.Vector3D, name: str = ''):
        self.edge = edge
        direction.normalize()
        self.direction = direction
        self.frame = volmdlr.Frame3D.from_point_and_vector(edge.start, direction, volmdlr.Z3D)
        self._x_periodicity = False

        Surface3D.__init__(self, name=name)

    @property
    def x_periodicity(self):
        if self._x_periodicity:
            return self._x_periodicity
        start = self.edge.start
        end = self.edge.end
        if start.is_close(end, 1e-4):
            return 1
        return None

    @x_periodicity.setter
    def x_periodicity(self, value):
        self._x_periodicity = value

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        """
        Transform a parametric (u, v) point into a 3D Cartesian point (x, y, z).

        # u = [0, 1] and v = z
        """
        u, v = point2d
        if abs(u) < 1e-7:
            u = 0.0
        if abs(v) < 1e-7:
            v = 0.

        point_at_curve_global = self.edge.point_at_abscissa(u * self.edge.length())
        point_at_curve_local = self.frame.global_to_local_coordinates(point_at_curve_global)
        # x, y, z = point_at_curve_local
        point_local = point_at_curve_local.translation(volmdlr.Vector3D(0, 0, v))
        return self.frame.local_to_global_coordinates(point_local)

    def point3d_to_2d(self, point3d):
        """
        Transform a 3D Cartesian point (x, y, z) into a parametric (u, v) point.
        """
        x, y, z = self.frame.global_to_local_coordinates(point3d)
        if abs(x) < 1e-7:
            x = 0.0
        if abs(y) < 1e-7:
            y = 0.0
        if abs(z) < 1e-7:
            z = 0.0
        v = z
        point_at_curve_local = volmdlr.Point3D(x, y, 0)
        point_at_curve_global = self.frame.local_to_global_coordinates(point_at_curve_local)

        u = self.edge.abscissa(point_at_curve_global) / self.edge.length()

        u = min(u, 1.0)
        return volmdlr.Point2D(u, v)

    def rectangular_cut(self, x1: float = 0.0, x2: float = 1.0,
                        y1: float = 0.0, y2: float = 1.0, name: str = ''):
        """
        Cut a rectangular piece of the ExtrusionSurface3D object and return a ExtrusionFace3D object.

        """
        p1 = volmdlr.Point2D(x1, y1)
        p2 = volmdlr.Point2D(x2, y1)
        p3 = volmdlr.Point2D(x2, y2)
        p4 = volmdlr.Point2D(x1, y2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([p1, p2, p3, p4])
        surface2d = Surface2D(outer_contour, [])
        return volmdlr.faces.ExtrusionFace3D(self, surface2d, name)

    def plot(self, ax=None, color='grey', alpha=0.5, z: float = 0.5):
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        for i in range(21):
            step = i / 20. * 0.5
            wire = self.edge.translation(step * self.frame.w)
            wire.plot(ax=ax, color=color, alpha=alpha)

        return ax

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        name = arguments[0][1:-1]
        edge = object_dict[arguments[1]]
        if edge.__class__ is volmdlr.wires.Ellipse3D:
            fullarcellipse = vme.FullArcEllipse3D(edge.point_at_abscissa(0), edge.major_axis, edge.minor_axis,
                                                  edge.center, edge.normal, edge.major_dir, edge.name)
            edge = fullarcellipse
            direction = -object_dict[arguments[2]]
            surface = cls(edge=edge, direction=direction, name=name)
            surface.x_periodicity = 1

        else:
            direction = object_dict[arguments[2]]
            surface = cls(edge=edge, direction=direction, name=name)
        return surface

    def arcellipse3d_to_2d(self, arcellipse3d):
        """
        Transformation of an arcellipse3d to 2d, in a cylindrical surface.

        """
        points = [self.point3d_to_2d(p)
                  for p in arcellipse3d.discretization_points(number_points=15)]

        bsplinecurve2d = vme.BSplineCurve2D.from_points_interpolation(points, degree=2)
        return [bsplinecurve2d]

    def fullarcellipse3d_to_2d(self, fullarcellipse3d):
        length = fullarcellipse3d.length()
        start = self.point3d_to_2d(fullarcellipse3d.start)
        end = self.point3d_to_2d(fullarcellipse3d.end)

        u3, z3 = self.point3d_to_2d(fullarcellipse3d.point_at_abscissa(0.01 * length))
        if u3 > 0.5:
            p1 = volmdlr.Point2D(1, start.y)
            p2 = volmdlr.Point2D(0, end.y)
        elif u3 < 0.5:
            p1 = volmdlr.Point2D(0, start.y)
            p2 = volmdlr.Point2D(1, end.y)
        else:
            raise NotImplementedError
        return [vme.LineSegment2D(p1, p2)]

    def linesegment2d_to_3d(self, linesegment2d):
        """
        Converts a BREP line segment 2D onto a 3D primitive on the surface.
        """
        start3d = self.point2d_to_3d(linesegment2d.start)
        end3d = self.point2d_to_3d(linesegment2d.end)
        u1, z1 = linesegment2d.start
        u2, z2 = linesegment2d.end
        if math.isclose(u1, u2, abs_tol=1e-4):
            return [vme.LineSegment3D(start3d, end3d)]
        if math.isclose(z1, z2, abs_tol=1e-4):
            if math.isclose(abs(u1 - u2), 1.0, abs_tol=1e-6):
                primitive = self.edge.translation(self.direction * z1)
                return [primitive]
                # if primitive.__name__
                # return [vme.FullArcEllipse3D()]
            # TODO: return self.edge translated and trimmed between u1 and u2
            raise NotImplementedError
        raise NotImplementedError

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        n = len(bspline_curve3d.control_points)
        points = [self.point3d_to_2d(point)
                  for point in bspline_curve3d.discretization_points(number_points=n)]
        start = points[0]
        end = points[-1]
        if not start.is_close(end):
            linesegment = vme.LineSegment2D(start, end)
            flag = True
            for pt in points:
                if not linesegment.point_belongs(pt):
                    flag = False
                    break
            if flag:
                return [linesegment]

        # Is this always True?
        n = len(bspline_curve3d.control_points)
        points = [self.point3d_to_2d(p)
                  for p in bspline_curve3d.discretization_points(number_points=n)]
        return [vme.BSplineCurve2D.from_points_interpolation(points, bspline_curve3d.degree, bspline_curve3d.periodic)]
        # raise NotImplementedError


class RevolutionSurface3D(PeriodicalSurface):
    """
    Defines a surface of revolution.

    :param wire: Wire.
    :type wire: Union[:class:`vmw.Wire3D`, :class:`vmw.Contour3D`]
    :param axis_point: Axis placement
    :type axis_point: :class:`volmdlr.Point3D`
    :param axis: Axis of revolution
    :type axis: :class:`volmdlr.Vector3D`
    """
    face_class = 'RevolutionFace3D'
    x_periodicity = volmdlr.TWO_PI
    y_periodicity = None

    def __init__(self, wire: Union[volmdlr.wires.Wire3D, volmdlr.wires.Contour3D],
                 axis_point: volmdlr.Point3D, axis: volmdlr.Vector3D, name: str = ''):
        self.wire = wire
        self.axis_point = axis_point
        self.axis = axis

        point1 = wire.point_at_abscissa(0)
        if point1 == axis_point:
            point1 = wire.point_at_abscissa(0.1 * wire.length())
        vector1 = point1 - axis_point
        w_vector = axis
        w_vector.normalize()
        u_vector = vector1 - vector1.vector_projection(w_vector)
        u_vector.normalize()
        v_vector = w_vector.cross(u_vector)
        self.frame = volmdlr.Frame3D(origin=axis_point, u=u_vector, v=v_vector, w=w_vector)

        PeriodicalSurface.__init__(self, name=name)

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        """
        Transform a parametric (u, v) point into a 3D Cartesian point (x, y, z).

        u = [0, 2pi] and v = [0, 1] into a
        """
        u, v = point2d
        point_at_curve = self.wire.point_at_abscissa(v * self.wire.length())
        point = point_at_curve.rotation(self.axis_point, self.axis, u)
        return point

    def point3d_to_2d(self, point3d):
        """
        Transform a 3D Cartesian point (x, y, z) into a parametric (u, v) point.
        """
        x, y, _ = self.frame.global_to_local_coordinates(point3d)
        if abs(x) < 1e-12:
            x = 0
        if abs(y) < 1e-12:
            y = 0
        u = math.atan2(y, x)

        point_at_curve = point3d.rotation(self.axis_point, self.axis, -u)
        v = self.wire.abscissa(point_at_curve) / self.wire.length()
        return volmdlr.Point2D(u, v)

    def rectangular_cut(self, x1: float, x2: float,
                        y1: float, y2: float, name: str = ''):
        """
        Cut a rectangular piece of the RevolutionSurface3D object and return a RevolutionFace3D object.

        """
        point1 = volmdlr.Point2D(x1, y1)
        point2 = volmdlr.Point2D(x2, y1)
        point3 = volmdlr.Point2D(x2, y2)
        point4 = volmdlr.Point2D(x1, y2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        surface2d = Surface2D(outer_contour, [])
        return volmdlr.faces.RevolutionFace3D(self, surface2d, name)

    def plot(self, ax=None, color='grey', alpha=0.5, number_curves: int = 20):
        """
        Plot rotated Revolution surface generatrix.

        :param number_curves: Number of curves to display.
        :type number_curves: int
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
        for i in range(number_curves + 1):
            theta = i / number_curves * volmdlr.TWO_PI
            wire = self.wire.rotation(self.axis_point, self.axis, theta)
            wire.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha))

        return ax

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a RevolutionSurface3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding RevolutionSurface3D object.
        :rtype: :class:`volmdlr.faces.RevolutionSurface3D`
        """
        name = arguments[0][1:-1]
        contour3d = object_dict[arguments[1]]
        axis_point, axis = object_dict[arguments[2]]
        return cls(wire=contour3d, axis_point=axis_point, axis=axis, name=name)

    def fullarc3d_to_2d(self, fullarc3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        length = fullarc3d.length()

        start = self.point3d_to_2d(fullarc3d.start)
        end = self.point3d_to_2d(fullarc3d.end)

        theta3, _ = self.point3d_to_2d(fullarc3d.point_at_abscissa(0.001 * length))
        theta4, _ = self.point3d_to_2d(fullarc3d.point_at_abscissa(0.98 * length))

        # make sure that the references points are not undefined
        if abs(theta3) == math.pi:
            theta3, z3 = self.point3d_to_2d(fullarc3d.point_at_abscissa(0.002 * length))
        if abs(theta4) == math.pi:
            theta4, _ = self.point3d_to_2d(fullarc3d.point_at_abscissa(0.97 * length))

        start, end = vm_parametric.arc3d_to_cylindrical_coordinates_verification(start, end, volmdlr.TWO_PI, theta3,
                                                                                 theta4)

        theta1, z1 = start
        _, z2 = end

        point1 = volmdlr.Point2D(theta1, z1)
        if theta1 > theta3:
            point2 = volmdlr.Point2D(theta1 + volmdlr.TWO_PI, z2)
        elif theta1 < theta3:
            point2 = volmdlr.Point2D(theta1 - volmdlr.TWO_PI, z2)
        return [vme.LineSegment2D(point1, point2)]


class BSplineSurface3D(Surface3D):
    """
    A class representing a 3D B-spline surface.

    A B-spline surface is a smooth surface defined by a set of control points and
    a set of basis functions called B-spline basis functions. The shape of the
    surface is determined by the position of the control points and can be
    modified by moving the control points.

    :param degree_u: The degree of the B-spline curve in the u direction.
    :type degree_u: int
    :param degree_v: The degree of the B-spline curve in the v direction.
    :type degree_v: int
    :param control_points: A list of 3D control points that define the shape of
        the surface.
    :type control_points: List[`volmdlr.Point3D`]
    :param nb_u: The number of control points in the u direction.
    :type nb_u: int
    :param nb_v: The number of control points in the v direction.
    :type nb_v: int
    :param u_multiplicities: A list of multiplicities for the knots in the u direction.
        The multiplicity of a knot is the number of times it appears in the knot vector.
    :type u_multiplicities: List[int]
    :param v_multiplicities: A list of multiplicities for the knots in the v direction.
        The multiplicity of a knot is the number of times it appears in the knot vector.
    :type v_multiplicities: List[int]
    :param u_knots: A list of knots in the u direction. The knots are real numbers that
        define the position of the control points along the u direction.
    :type u_knots: List[float]
    :param v_knots: A list of knots in the v direction. The knots are real numbers that
        define the position of the control points along the v direction.
    :type v_knots: List[float]
    :param weights: (optional) A list of weights for the control points. The weights
        can be used to adjust the influence of each control point on the shape of the
        surface. Default is None.
    :type weights: List[float]
    :param name: (optional) A name for the surface. Default is an empty string.
    :type name: str
    """
    face_class = "BSplineFace3D"
    _non_serializable_attributes = ["surface", "curves"]

    def __init__(self, degree_u: int, degree_v: int, control_points: List[volmdlr.Point3D], nb_u: int, nb_v: int,
                 u_multiplicities: List[int], v_multiplicities: List[int], u_knots: List[float], v_knots: List[float],
                 weights: List[float] = None, name: str = ''):
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
        for point in control_points:
            points_row.append(point)
            if i == nb_v:
                self.control_points_table.append(points_row)
                points_row = []
                i = 1
            else:
                i += 1
        if weights is None:
            surface = BSpline.Surface()
            points = [(control_points[i][0], control_points[i][1],
                       control_points[i][2]) for i in range(len(control_points))]

        else:
            surface = NURBS.Surface()
            points = [(control_points[i][0] * weights[i], control_points[i][1] * weights[i],
                       control_points[i][2] * weights[i], weights[i]) for i in range(len(control_points))]
        surface.degree_u = degree_u
        surface.degree_v = degree_v
        surface.set_ctrlpts(points, nb_u, nb_v)
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
        self.curves = extract_curves(surface, extract_u=True, extract_v=True)
        # self.points = [volmdlr.Point3D(*p) for p in surface_points]
        Surface3D.__init__(self, name=name)

        # Hidden Attributes
        self._displacements = None
        self._grids2d = None
        self._grids2d_deformed = None
        self._bbox = None

        self._x_periodicity = False  # Use False instead of None because None is a possible value of x_periodicity
        self._y_periodicity = False

    @property
    def x_periodicity(self):
        """
        Evaluates the periodicity of the surface in u direction.
        """
        if self._x_periodicity is False:
            u = self.curves['u']
            a, b = self.surface.domain[0]
            u0 = u[0]
            point_at_a = u0.evaluate_single(a)
            point_at_b = u0.evaluate_single(b)
            if npy.linalg.norm(npy.array(point_at_b) - npy.array(point_at_a)) < 1e-6:
                self._x_periodicity = self.surface.range[0]
            else:
                self._x_periodicity = None
        return self._x_periodicity

    @property
    def y_periodicity(self):
        """
        Evaluates the periodicity of the surface in v direction.
        """
        if self._y_periodicity is False:
            v = self.curves['v']
            c, d = self.surface.domain[1]
            v0 = v[0]
            point_at_c = v0.evaluate_single(c)
            point_at_d = v0.evaluate_single(d)
            if npy.linalg.norm(npy.array(point_at_d) - npy.array(point_at_c)) < 1e-6:
                self._y_periodicity = self.surface.range[1]
            else:
                self._y_periodicity = None
        return self._y_periodicity

    @property
    def bounding_box(self):
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    def _bounding_box(self):
        """
        Computes the bounding box of the surface.

        """
        min_bounds, max_bounds = self.surface.bbox
        xmin, ymin, zmin = min_bounds
        xmax, ymax, zmax = max_bounds
        return volmdlr.core.BoundingBox(xmin, xmax, ymin, ymax, zmin, zmax)

    def control_points_matrix(self, coordinates):
        """
        Define control points like a matrix, for each coordinate: x:0, y:1, z:2.

        """

        points = npy.empty((self.nb_u, self.nb_v))
        for i in range(0, self.nb_u):
            for j in range(0, self.nb_v):
                points[i][j] = self.control_points_table[i][j][coordinates]
        return points

    # Knots_vector
    def knots_vector_u(self):
        """
        Compute the global knot vector (u direction) based on knot elements and multiplicities.

        """

        knots = self.u_knots
        multiplicities = self.u_multiplicities

        knots_vec = []
        for i, knot in enumerate(knots):
            for _ in range(0, multiplicities[i]):
                knots_vec.append(knot)
        return knots_vec

    def knots_vector_v(self):
        """
        Compute the global knot vector (v direction) based on knot elements and multiplicities.

        """

        knots = self.v_knots
        multiplicities = self.v_multiplicities

        knots_vec = []
        for i, knot in enumerate(knots):
            for _ in range(0, multiplicities[i]):
                knots_vec.append(knot)
        return knots_vec

    def basis_functions_u(self, u, k, i):
        """
        Compute basis functions Bi in u direction for u=u and degree=k.

        """

        # k = self.degree_u
        t = self.knots_vector_u()

        if k == 0:
            return 1.0 if t[i] <= u < t[i + 1] else 0.0
        if t[i + k] == t[i]:
            param_c1 = 0.0
        else:
            param_c1 = (u - t[i]) / (t[i + k] - t[i]) * self.basis_functions_u(u, k - 1, i)
        if t[i + k + 1] == t[i + 1]:
            param_c2 = 0.0
        else:
            param_c2 = (t[i + k + 1] - u) / (t[i + k + 1] - t[i + 1]) * self.basis_functions_u(u, k - 1, i + 1)
        return param_c1 + param_c2

    def basis_functions_v(self, v, k, i):
        """
        Compute basis functions Bi in v direction for v=v and degree=k.

        """

        # k = self.degree_u
        knots = self.knots_vector_v()

        if k == 0:
            return 1.0 if knots[i] <= v < knots[i + 1] else 0.0
        if knots[i + k] == knots[i]:
            param_c1 = 0.0
        else:
            param_c1 = (v - knots[i]) / (knots[i + k] - knots[i]) * self.basis_functions_v(v, k - 1, i)
        if knots[i + k + 1] == knots[i + 1]:
            param_c2 = 0.0
        else:
            param_c2 = (knots[i + k + 1] - v) / (knots[i + k + 1] - knots[i + 1]) * self.basis_functions_v(v, k - 1,
                                                                                                           i + 1)
        return param_c1 + param_c2

    def blending_vector_u(self, u):
        """
        Compute a vector of basis_functions in u direction for u=u.
        """

        blending_vect = npy.empty((1, self.nb_u))
        for j in range(0, self.nb_u):
            blending_vect[0][j] = self.basis_functions_u(u, self.degree_u, j)

        return blending_vect

    def blending_vector_v(self, v):
        """
        Compute a vector of basis_functions in v direction for v=v.

        """

        blending_vect = npy.empty((1, self.nb_v))
        for j in range(0, self.nb_v):
            blending_vect[0][j] = self.basis_functions_v(v, self.degree_v, j)

        return blending_vect

    def blending_matrix_u(self, u):
        """
        Compute a matrix of basis_functions in u direction for a vector u like [0,1].

        """

        blending_mat = npy.empty((len(u), self.nb_u))
        for i, u_i in enumerate(u):
            for j in range(self.nb_u):
                blending_mat[i][j] = self.basis_functions_u(u_i, self.degree_u, j)
        return blending_mat

    def blending_matrix_v(self, v):
        """
        Compute a matrix of basis_functions in v direction for a vector v like [0,1].

        """

        blending_mat = npy.empty((len(v), self.nb_v))
        for i, v_i in enumerate(v):
            for j in range(self.nb_v):
                blending_mat[i][j] = self.basis_functions_v(v_i, self.degree_v, j)
        return blending_mat

    def point2d_to_3d(self, point2d: volmdlr.Point2D):
        u, v = point2d
        u = min(max(u, 0), 1)
        v = min(max(v, 0), 1)
        return volmdlr.Point3D(*evaluate_single((u, v), self.surface.data,
                                                self.surface.rational,
                                                self.surface.evaluator._span_func))
        # uses derivatives for performance because it's already compiled
        # return volmdlr.Point3D(*self.derivatives(u, v, 0)[0][0])
        # return volmdlr.Point3D(*self.surface.evaluate_single((x, y)))

    def point3d_to_2d(self, point3d: volmdlr.Point3D, tol=1e-5):
        """
        Evaluates the parametric coordinates (u, v) of a 3D point (x, y, z).

        :param point3d: A 3D point to be evaluated.
        :type point3d: :class:`volmdlr.Point3D`
        :param tol: Tolerance to accept the results.
        :type tol: float
        :return: The parametric coordinates (u, v) of the point.
        :rtype: :class:`volmdlr.Point2D`
        """

        def f(x):
            return point3d.point_distance(self.point2d_to_3d(volmdlr.Point2D(x[0], x[1])))

        def fun(x):
            derivatives = self.derivatives(x[0], x[1], 1)
            r = derivatives[0][0] - point3d
            f_value = r.norm() + 1e-32
            jacobian = npy.array([r.dot(derivatives[1][0]) / f_value, r.dot(derivatives[0][1]) / f_value])
            return f_value, jacobian

        min_bound_x, max_bound_x = self.surface.domain[0]
        min_bound_y, max_bound_y = self.surface.domain[1]

        delta_bound_x = max_bound_x - min_bound_x
        delta_bound_y = max_bound_y - min_bound_y
        x0s = [((min_bound_x + max_bound_x) / 2, (min_bound_y + max_bound_y) / 2),
               ((min_bound_x + max_bound_x) / 2, min_bound_y + delta_bound_y / 10),
               ((min_bound_x + max_bound_x) / 2, max_bound_y - delta_bound_y / 10),
               ((min_bound_x + max_bound_x) / 4, min_bound_y + delta_bound_y / 10),
               (max_bound_x - delta_bound_x / 4, min_bound_y + delta_bound_y / 10),
               ((min_bound_x + max_bound_x) / 4, max_bound_y - delta_bound_y / 10),
               (max_bound_x - delta_bound_x / 4, max_bound_y - delta_bound_y / 10),
               (min_bound_x + delta_bound_x / 10, min_bound_y + delta_bound_y / 10),
               (min_bound_x + delta_bound_x / 10, max_bound_y - delta_bound_y / 10),
               (max_bound_x - delta_bound_x / 10, min_bound_y + delta_bound_y / 10),
               (max_bound_x - delta_bound_x / 10, max_bound_y - delta_bound_y / 10)]

        # Sort the initial conditions
        x0s.sort(key=f)

        # Find the parametric coordinates of the point
        results = []
        for x0 in x0s:
            res = minimize(fun, x0=npy.array(x0), jac=True,
                           bounds=[(min_bound_x, max_bound_x),
                                   (min_bound_y, max_bound_y)])
            if res.fun <= tol:
                return volmdlr.Point2D(*res.x)

            results.append((res.x, res.fun))

        return volmdlr.Point2D(*min(results, key=lambda r: r[1])[0])

    def linesegment2d_to_3d(self, linesegment2d):
        # TODO: this is a non exact method!
        lth = linesegment2d.length()
        points = [self.point2d_to_3d(
            linesegment2d.point_at_abscissa(i * lth / 20.)) for i in range(21)]

        linesegment = vme.LineSegment3D(points[0], points[-1])
        flag = True
        flag_arc = False
        for point in points:
            if not linesegment.point_belongs(point, abs_tol=1e-4):
                flag = False
                break
        if not flag:
            arc = vme.Arc3D(points[0], points[10], points[-1])
            flag_arc = all(arc.point_belongs(point, abs_tol=1e-4) for point in points)

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

        if flag and not flag_arc:
            # All the points are on the same LineSegment3D
            primitives = [linesegment]
        elif flag_arc:
            primitives = [arc]
        else:
            primitives = [vme.BSplineCurve3D.from_points_interpolation(
                points, min(self.degree_u, self.degree_v), periodic=periodic)]
        return primitives

    def linesegment3d_to_2d(self, linesegment3d):
        """
        A line segment on a BSplineSurface3D will be in any case a line in 2D?.

        """
        start = self.point3d_to_2d(linesegment3d.start)
        end = self.point3d_to_2d(linesegment3d.end)
        if self.x_periodicity:
            if start.x != end.x:
                end = volmdlr.Point2D(start.x, end.y)
            if not start.is_close(end):
                return [vme.LineSegment2D(start, end)]
            return None
        if self.y_periodicity:
            if start.y != end.y:
                end = volmdlr.Point2D(end.x, start.y)
            if not start.is_close(end):
                return [vme.LineSegment2D(start, end)]
            return None
        if start.is_close(end):
            return None
        return [vme.LineSegment2D(start, end)]

    def _repair_periodic_boundary_points(self, curve3d, points_2d, direction_periodicity):
        """
        Verifies points at boundary on a periodic BSplineSurface3D.

        :param points_2d: List of `volmdlr.Point2D` after transformation from 3D Cartesian coordinates
        :type points_2d: List[volmdlr.Point2D]
        :param direction_periodicity: should be 'x' if x_periodicity or 'y' if y periodicity
        :type direction_periodicity: str
        """
        lth = curve3d.length()
        start = points_2d[0]
        end = points_2d[-1]
        points = points_2d
        pt_after_start = self.point3d_to_2d(curve3d.point_at_abscissa(0.1 * lth))
        pt_before_end = self.point3d_to_2d(curve3d.point_at_abscissa(0.9 * lth))
        # pt_after_start = points[1]
        # pt_before_end = points[-2]

        if direction_periodicity == 'x':
            i = 0
        else:
            i = 1
        min_bound, max_bound = self.surface.domain[i]
        delta = max_bound - min_bound

        if math.isclose(start[i], min_bound, abs_tol=1e-4) and pt_after_start[i] > 0.5 * delta:
            start[i] = max_bound
        elif math.isclose(start[i], max_bound, abs_tol=1e-4) and pt_after_start[i] < 0.5 * delta:
            start[i] = min_bound

        if math.isclose(end[i], min_bound, abs_tol=1e-4) and pt_before_end[i] > 0.5 * delta:
            end[i] = max_bound
        elif math.isclose(end[i], max_bound, abs_tol=1e-4) and pt_before_end[i] < 0.5 * delta:
            end[i] = min_bound

        points[0] = start
        points[-1] = end

        boundary = [(math.isclose(p[i], max_bound, abs_tol=1e-4) or math.isclose(p[i], min_bound, abs_tol=1e-4)) for
                    p in points]
        if all(boundary):
            # if the line is at the boundary of the surface domain, we take the first point as reference
            t = max_bound if math.isclose(points[0][i], max_bound, abs_tol=1e-4) else min_bound
            if direction_periodicity == 'x':
                points = [volmdlr.Point2D(t, p[1]) for p in points]
            else:
                points = [volmdlr.Point2D(p[0], t) for p in points]

        return points

    def bsplinecurve3d_to_2d(self, bspline_curve3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        # TODO: enhance this, it is a non exact method!
        # TODO: bsplinecurve can be periodic but not around the bsplinesurface
        flag = False
        if not bspline_curve3d.points[0].is_close(bspline_curve3d.points[-1]):
            bsc_linesegment = vme.LineSegment3D(bspline_curve3d.points[0],
                                                bspline_curve3d.points[-1])
            flag = True
            for point in bspline_curve3d.points:
                if not bsc_linesegment.point_belongs(point):
                    flag = False
                    break

        if self.x_periodicity and not self.y_periodicity \
                and bspline_curve3d.periodic:
            point1 = self.point3d_to_2d(bspline_curve3d.points[0])
            p1_sup = self.point3d_to_2d(bspline_curve3d.points[0])
            new_x = point1.x - p1_sup.x + self.x_periodicity
            new_x = new_x if 0 <= new_x else 0
            reverse = False
            if new_x < 0:
                new_x = 0
            elif math.isclose(new_x, self.x_periodicity, abs_tol=1e-5):
                new_x = 0
                reverse = True

            linesegments = [
                vme.LineSegment2D(
                    volmdlr.Point2D(new_x, point1.y),
                    volmdlr.Point2D(self.x_periodicity, point1.y))]
            if reverse:
                linesegments[0] = linesegments[0].reverse()

        elif self.y_periodicity and not self.x_periodicity \
                and bspline_curve3d.periodic:
            point1 = self.point3d_to_2d(bspline_curve3d.points[0])
            p1_sup = self.point3d_to_2d(bspline_curve3d.points[0])
            new_y = point1.y - p1_sup.y + self.y_periodicity
            new_y = new_y if 0 <= new_y else 0
            reverse = False
            if new_y < 0:
                new_y = 0
            elif math.isclose(new_y, self.y_periodicity, abs_tol=1e-5):
                new_y = 0
                reverse = True

            linesegments = [
                vme.LineSegment2D(
                    volmdlr.Point2D(point1.x, new_y),
                    volmdlr.Point2D(point1.x, self.y_periodicity))]
            if reverse:
                linesegments[0] = linesegments[0].reverse()

        elif self.x_periodicity and self.y_periodicity \
                and bspline_curve3d.periodic:
            raise NotImplementedError

        if flag:
            x_perio = self.x_periodicity if self.x_periodicity is not None \
                else 1.
            y_perio = self.y_periodicity if self.y_periodicity is not None \
                else 1.

            point1 = self.point3d_to_2d(bspline_curve3d.points[0])
            point2 = self.point3d_to_2d(bspline_curve3d.points[-1])

            if point1.is_close(point2):
                print('BSplineCruve3D skipped because it is too small')
                linesegments = None
            else:
                p1_sup = self.point3d_to_2d(bspline_curve3d.points[0])
                p2_sup = self.point3d_to_2d(bspline_curve3d.points[-1])
                if self.x_periodicity and point1.point_distance(p1_sup) > 1e-5:
                    point1.x -= p1_sup.x - x_perio
                    point2.x -= p2_sup.x - x_perio
                if self.y_periodicity and point1.point_distance(p1_sup) > 1e-5:
                    point1.y -= p1_sup.y - y_perio
                    point2.y -= p2_sup.y - y_perio
                linesegments = [vme.LineSegment2D(point1, point2)]
            # How to check if end of surface overlaps start or the opposite ?
        else:
            lth = bspline_curve3d.length()
            # uses straight line length to improve performance
            # lth = bspline_curve3d.start.point_distance(bspline_curve3d.end)
            if lth > 1e-5:
                n = len(bspline_curve3d.control_points)
                points = [self.point3d_to_2d(p) for p in bspline_curve3d.discretization_points(number_points=n)]

                if self.x_periodicity:
                    points = self._repair_periodic_boundary_points(bspline_curve3d, points, 'x')
                    if bspline_curve3d.periodic and points[0].is_close(points[-1]):
                        u_min, u_max = bspline_curve3d.curve.domain
                        if math.isclose(points[0].x, u_min, abs_tol=1e-6):
                            should_be_umax = (u_max - points[1].x) < (points[1].x - u_min)
                            if should_be_umax:
                                points[0] = volmdlr.Point2D(u_max, points[0].y)
                            else:
                                points[-1] = volmdlr.Point2D(u_max, points[-1].y)
                        elif math.isclose(points[0].x, u_max, abs_tol=1e-6):
                            should_be_umin = (u_max - points[1].x) > (points[1].x - u_min)
                            if should_be_umin:
                                points[0] = volmdlr.Point2D(u_min, points[0].y)
                            else:
                                points[-1] = volmdlr.Point2D(u_min, points[-1].y)
                if self.y_periodicity:
                    points = self._repair_periodic_boundary_points(bspline_curve3d, points, 'y')
                    if bspline_curve3d.periodic and points[0].is_close(points[-1]):
                        u_min, u_max = bspline_curve3d.curve.domain
                        if math.isclose(points[0].y, u_min, abs_tol=1e-6):
                            should_be_umax = (u_max - points[1].y) < (points[1].y - u_min)
                            if should_be_umax:
                                points[0] = volmdlr.Point2D(points[0].x, u_max)
                            else:
                                points[-1] = volmdlr.Point2D(points[-1].x, u_max)
                        elif math.isclose(points[0].y, u_max, abs_tol=1e-6):
                            should_be_umin = (u_max - points[1].y) > (points[1].y - u_min)
                            if should_be_umin:
                                points[0] = volmdlr.Point2D(points[0].x, u_min)
                            else:
                                points[-1] = volmdlr.Point2D(points[-1].x, u_min)

                if not points[0].is_close(points[-1]) and not bspline_curve3d.periodic:
                    linesegment = vme.LineSegment2D(points[0], points[-1])
                    flag_line = True
                    for point in points:
                        if not linesegment.point_belongs(point, abs_tol=1e-4):
                            flag_line = False
                            break
                    if flag_line:
                        return [linesegment]

                if self.x_periodicity:
                    points = self._repair_periodic_boundary_points(bspline_curve3d, points, 'x')

                if self.y_periodicity:
                    points = self._repair_periodic_boundary_points(bspline_curve3d, points, 'y')
                # points_ = [points[0]]
                # for point in points[1:]:
                #     if not point.is_close(points[-1]):
                #         points_.append(point)
                # if len(points_) < 2:
                #     return []

                return [vme.BSplineCurve2D.from_points_interpolation(
                    points=points, degree=bspline_curve3d.degree, periodic=bspline_curve3d.periodic)]

            if 1e-6 < lth <= 1e-5:
                linesegments = [vme.LineSegment2D(
                    self.point3d_to_2d(bspline_curve3d.start),
                    self.point3d_to_2d(bspline_curve3d.end))]
            else:
                print('BSplineCruve3D skipped because it is too small')
                linesegments = None

        return linesegments

    def arc3d_to_2d(self, arc3d):
        """
        Converts the primitive from 3D spatial coordinates to its equivalent 2D primitive in the parametric space.
        """
        number_points = max(self.nb_u, self.nb_v)
        degree = max(self.degree_u, self.degree_v)
        points = [self.point3d_to_2d(point3d) for point3d in arc3d.discretization_points(number_points=number_points)]
        start = points[0]
        end = points[-1]
        min_bound_x, max_bound_x = self.surface.domain[0]
        min_bound_y, max_bound_y = self.surface.domain[1]
        if self.x_periodicity:
            points = self._repair_periodic_boundary_points(arc3d, points, 'x')
            start = points[0]
            end = points[-1]
            if start.is_close(end):
                if math.isclose(start.x, min_bound_x, abs_tol=1e-4):
                    end.x = max_bound_x
                else:
                    end.x = min_bound_x
        if self.y_periodicity:
            points = self._repair_periodic_boundary_points(arc3d, points, 'y')
            start = points[0]
            end = points[-1]
            if start.is_close(end):
                if math.isclose(start.y, min_bound_y, abs_tol=1e-4):
                    end.y = max_bound_y
                else:
                    end.y = min_bound_y
        if start.is_close(end):
            return []
        linesegment = vme.LineSegment2D(start, end)
        flag = True
        for point in points:
            if not linesegment.point_belongs(point):
                flag = False
                break
        if flag:
            return [linesegment]
        return [vme.BSplineCurve2D.from_points_interpolation(points, degree)]

    def arc2d_to_3d(self, arc2d):
        number_points = math.ceil(arc2d.angle * 7) + 1  # 7 points per radian
        length = arc2d.length()
        points = [self.point2d_to_3d(arc2d.point_at_abscissa(i * length / (number_points - 1)))
                  for i in range(number_points)]
        return [vme.BSplineCurve3D.from_points_interpolation(
            points, max(self.degree_u, self.degree_v))]

    def rectangular_cut(self, u1: float, u2: float,
                        v1: float, v2: float, name: str = ''):
        """
        Cut a rectangular piece of the BSplineSurface3D object and return a BSplineFace3D object.

        """
        point1 = volmdlr.Point2D(u1, v1)
        point2 = volmdlr.Point2D(u2, v1)
        point3 = volmdlr.Point2D(u2, v2)
        point4 = volmdlr.Point2D(u1, v2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        surface = Surface2D(outer_contour, [])
        return BSplineFace3D(self, surface, name)  # PlaneFace3D

    def rotation(self, center: volmdlr.Vector3D,
                 axis: volmdlr.Vector3D, angle: float):
        """
        BSplineSurface3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated BSplineSurface3D
        """
        new_control_points = [p.rotation(center, axis, angle)
                              for p in self.control_points]
        new_bsplinesurface3d = BSplineSurface3D(self.degree_u, self.degree_v,
                                                new_control_points, self.nb_u,
                                                self.nb_v,
                                                self.u_multiplicities,
                                                self.v_multiplicities,
                                                self.u_knots, self.v_knots,
                                                self.weights, self.name)
        return new_bsplinesurface3d

    def rotation_inplace(self, center: volmdlr.Vector3D,
                         axis: volmdlr.Vector3D, angle: float):
        """
        BSplineSurface3D rotation. Object is updated inplace.

        :param center: rotation center.
        :type center: `volmdlr.Vector3D`
        :param axis: rotation axis.
        :type axis: `volmdlr.Vector3D`
        :param angle: rotation angle.
        :type angle: float
        :return: None, BSplineSurface3D is updated inplace
        :rtype: None
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_bsplinesurface3d = self.rotation(center, axis, angle)
        self.control_points = new_bsplinesurface3d.control_points
        self.surface = new_bsplinesurface3d.surface

    def translation(self, offset: volmdlr.Vector3D):
        """
        BSplineSurface3D translation.

        :param offset: translation vector
        :return: A new translated BSplineSurface3D
        """
        new_control_points = [p.translation(offset) for p in
                              self.control_points]
        new_bsplinesurface3d = BSplineSurface3D(self.degree_u, self.degree_v,
                                                new_control_points, self.nb_u,
                                                self.nb_v,
                                                self.u_multiplicities,
                                                self.v_multiplicities,
                                                self.u_knots, self.v_knots,
                                                self.weights, self.name)

        return new_bsplinesurface3d

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        BSplineSurface3D translation. Object is updated inplace.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_bsplinesurface3d = self.translation(offset)
        self.control_points = new_bsplinesurface3d.control_points
        self.surface = new_bsplinesurface3d.surface

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new BSplineSurface3D.

        side = 'old' or 'new'
        """
        new_control_points = [p.frame_mapping(frame, side) for p in
                              self.control_points]
        new_bsplinesurface3d = BSplineSurface3D(self.degree_u, self.degree_v,
                                                new_control_points, self.nb_u,
                                                self.nb_v,
                                                self.u_multiplicities,
                                                self.v_multiplicities,
                                                self.u_knots, self.v_knots,
                                                self.weights, self.name)
        return new_bsplinesurface3d

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        new_bsplinesurface3d = self.frame_mapping(frame, side)
        self.control_points = new_bsplinesurface3d.control_points
        self.surface = new_bsplinesurface3d.surface

    def plot(self, ax=None, color='grey', alpha=0.5):
        u_curves = [vme.BSplineCurve3D.from_geomdl_curve(u) for u in self.curves['u']]
        v_curves = [vme.BSplineCurve3D.from_geomdl_curve(v) for v in self.curves['v']]
        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')
        for u in u_curves:
            u.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha))
        for v in v_curves:
            v.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha))
        for point in self.control_points:
            point.plot(ax, color=color, alpha=alpha)
        return ax

    def simplify_surface(self):
        """
        Verifies if BSplineSurface3D could be a Plane3D.

        :return: simplified a planar surface if possible, otherwise, returns self.
        """
        points = [self.control_points[0]]
        vector_list = []
        for point in self.control_points[1:]:
            vector = point - points[0]
            is_colinear = any(vector.is_colinear_to(other_vector) for other_vector in vector_list)
            if not point_in_list(point, points) and not is_colinear:
                points.append(point)
                vector_list.append(vector)
                if len(points) == 3:
                    plane3d = Plane3D.from_3_points(*points)
                    if all(plane3d.point_on_surface(point) for point in self.control_points):
                        return plane3d
                    break
        return self

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a BSplineSurface3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding BSplineSurface3D object.
        :rtype: :class:`volmdlr.faces.BSplineSurface3D`
        """
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
        if not bsplinesurface.x_periodicity and not bsplinesurface.y_periodicity:
            bsplinesurface = bsplinesurface.simplify_surface()
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
                point_ids += f'#{point_id},'
                current_id = point_id + 1
            point_ids = point_ids[:-1]
            point_ids += '),'
            point_matrix_ids += point_ids
        point_matrix_ids = point_matrix_ids[:-1]
        point_matrix_ids += ')'

        u_close = '.T.' if self.x_periodicity else '.F.'
        v_close = '.T.' if self.y_periodicity else '.F.'

        content += f"#{current_id} = B_SPLINE_SURFACE_WITH_KNOTS('{self.name}',{self.degree_u},{self.degree_v}," \
                   f"{point_matrix_ids},.UNSPECIFIED.,{u_close},{v_close},.F.,{tuple(self.u_multiplicities)}," \
                   f"{tuple(self.v_multiplicities)},{tuple(self.u_knots)},{tuple(self.v_knots)},.UNSPECIFIED.);\n"
        return content, [current_id]

    def grid3d(self, grid2d: volmdlr.grid.Grid2D):
        """
        Generate 3d grid points of a Bspline surface, based on a Grid2D.

        """

        if not self._grids2d:
            self._grids2d = grid2d

        points_2d = grid2d.points
        points_3d = [self.point2d_to_3d(point2d) for point2d in points_2d]

        return points_3d

    def grid2d_deformed(self, grid2d: volmdlr.grid.Grid2D):
        """
        Dimension and deform a Grid2D points based on a Bspline surface.

        """

        points_2d = grid2d.points
        points_3d = self.grid3d(grid2d)

        points_x, points_y = grid2d.points_xy

        # Parameters
        index_x = {}  # grid point position(i,j), x coordinates position in X(unknown variable)
        index_y = {}  # grid point position(i,j), y coordinates position in X(unknown variable)
        index_points = {}  # grid point position(j,i), point position in points_2d (or points_3d)
        k_index, p_index = 0, 0
        for i in range(0, points_x):
            for j in range(0, points_y):
                index_x.update({(j, i): k_index})
                index_y.update({(j, i): k_index + 1})
                index_points.update({(j, i): p_index})
                k_index = k_index + 2
                p_index = p_index + 1

        equation_points = []  # points combination to compute distances between 2D and 3D grid points
        for i in range(0, points_y):  # row from (0,i)
            for j in range(1, points_x):
                equation_points.append(((0, i), (j, i)))
        for i in range(0, points_x):  # column from (i,0)
            for j in range(1, points_y):
                equation_points.append(((i, 0), (i, j)))
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

        # geodesic distances between 3D grid points (based on points combination [equation_points])
        geodesic_distances = []
        for point in equation_points:
            geodesic_distances.append((self.geodesic_distance(
                points_3d[index_points[point[0]]], points_3d[index_points[point[1]]])) ** 2)

        # System of nonlinear equations
        def non_linear_equations(xparam):
            vector_f = npy.empty(len(equation_points) + 2)
            idx = 0
            for idx, point_ in enumerate(equation_points):
                vector_f[idx] = abs((xparam[index_x[point_[0]]] ** 2 +
                                     xparam[index_x[point_[1]]] ** 2 +
                                     xparam[index_y[point_[0]]] ** 2 +
                                     xparam[index_y[point_[1]]] ** 2 -
                                     2 *
                                     xparam[index_x[point_[0]]] *
                                     xparam[index_x[point_[1]]] -
                                     2 *
                                     xparam[index_y[point_[0]]] *
                                     xparam[index_y[point_[1]]] -
                                     geodesic_distances[idx]) /
                                    geodesic_distances[idx])

            vector_f[idx + 1] = xparam[0] * 1000
            vector_f[idx + 2] = xparam[1] * 1000

            return vector_f

        # Solution with "least_squares"
        x_init = []  # initial guess (2D grid points)
        for point in points_2d:
            x_init.append(point[0])
            x_init.append(point[1])
        z = least_squares(non_linear_equations, x_init)

        points_2d_deformed = [volmdlr.Point2D(z.x[i], z.x[i + 1])
                              for i in range(0, len(z.x), 2)]  # deformed 2d grid points

        grid2d_deformed = volmdlr.grid.Grid2D.from_points(points=points_2d_deformed,
                                                          points_dim_1=points_x,
                                                          direction=grid2d.direction)

        self._grids2d_deformed = grid2d_deformed

        return points_2d_deformed

    def grid2d_deformation(self, grid2d: volmdlr.grid.Grid2D):
        """
        Compute the deformation/displacement (dx/dy) of a Grid2D based on a Bspline surface.

        """

        if not self._grids2d_deformed:
            self.grid2d_deformed(grid2d)

        displacement = self._grids2d_deformed.displacement_compared_to(grid2d)
        self._displacements = displacement

        return displacement

    def point2d_parametric_to_dimension(self, point2d: volmdlr.Point3D, grid2d: volmdlr.grid.Grid2D):
        """
        Convert a point 2d from the parametric to the dimensioned frame.

        """

        # Check if the 0<point2d.x<1 and 0<point2d.y<1
        if point2d.x < 0:
            point2d.x = 0
        elif point2d.x > 1:
            point2d.x = 1
        if point2d.y < 0:
            point2d.y = 0
        elif point2d.y > 1:
            point2d.y = 1

        if self._grids2d == grid2d:
            points_2d = self._grids2d.points
        else:
            points_2d = grid2d.points
            self._grids2d = grid2d

        if self._displacements is not None:
            displacement = self._displacements
        else:
            displacement = self.grid2d_deformation(grid2d)

        points_x, points_y = grid2d.points_xy

        # Parameters
        index_points = {}  # grid point position(j,i), point position in points_2d (or points_3d)
        p_index = 0
        for i in range(0, points_x):
            for j in range(0, points_y):
                index_points.update({(j, i): p_index})
                p_index = p_index + 1

        # Form function "Finite Elements"
        def form_function(s_param, t_param):
            empty_n = npy.empty(4)
            empty_n[0] = (1 - s_param) * (1 - t_param) / 4
            empty_n[1] = (1 + s_param) * (1 - t_param) / 4
            empty_n[2] = (1 + s_param) * (1 + t_param) / 4
            empty_n[3] = (1 - s_param) * (1 + t_param) / 4
            return empty_n

        finite_elements_points = []  # 2D grid points index that define one element
        for j in range(0, points_y - 1):
            for i in range(0, points_x - 1):
                finite_elements_points.append(((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)))
        finite_elements = []  # finite elements defined with closed polygon
        for point in finite_elements_points:
            finite_elements.append(
                volmdlr.wires.ClosedPolygon2D((points_2d[index_points[point[0]]],
                                               points_2d[index_points[point[1]]],
                                               points_2d[index_points[point[2]]],
                                               points_2d[index_points[point[3]]])))
        k = 0
        for k, point in enumerate(finite_elements_points):
            if (volmdlr.wires.Contour2D(finite_elements[k].primitives).point_belongs(
                    point2d)  # finite_elements[k].point_belongs(point2d)
                    or volmdlr.wires.Contour2D(finite_elements[k].primitives).point_over_contour(point2d)
                    or ((points_2d[index_points[point[0]]][0] < point2d.x <
                         points_2d[index_points[point[1]]][0])
                        and point2d.y == points_2d[index_points[point[0]]][1])
                    or ((points_2d[index_points[point[1]]][1] < point2d.y <
                         points_2d[index_points[point[2]]][1])
                        and point2d.x == points_2d[index_points[point[1]]][0])
                    or ((points_2d[index_points[point[3]]][0] < point2d.x <
                         points_2d[index_points[point[2]]][0])
                        and point2d.y == points_2d[index_points[point[1]]][1])
                    or ((points_2d[index_points[point[0]]][1] < point2d.y <
                         points_2d[index_points[point[3]]][1])
                        and point2d.x == points_2d[index_points[point[0]]][0])):
                break

        x0 = points_2d[index_points[finite_elements_points[k][0]]][0]
        y0 = points_2d[index_points[finite_elements_points[k][0]]][1]
        x1 = points_2d[index_points[finite_elements_points[k][1]]][0]
        y2 = points_2d[index_points[finite_elements_points[k][2]]][1]
        x = point2d.x
        y = point2d.y
        s = 2 * ((x - x0) / (x1 - x0)) - 1
        t = 2 * ((y - y0) / (y2 - y0)) - 1

        n = form_function(s, t)
        dx = npy.array([displacement[index_points[finite_elements_points[k][0]]][0],
                        displacement[index_points[finite_elements_points[k][1]]][0],
                        displacement[index_points[finite_elements_points[k][2]]][0],
                        displacement[index_points[finite_elements_points[k][3]]][0]])
        dy = npy.array([displacement[index_points[finite_elements_points[k][0]]][1],
                        displacement[index_points[finite_elements_points[k][1]]][1],
                        displacement[index_points[finite_elements_points[k][2]]][1],
                        displacement[index_points[finite_elements_points[k][3]]][1]])

        return volmdlr.Point2D(point2d.x + npy.transpose(n).dot(dx), point2d.y + npy.transpose(n).dot(dy))

    def point3d_to_2d_with_dimension(self, point3d: volmdlr.Point3D, grid2d: volmdlr.grid.Grid2D):
        """
        Compute the point2d of a point3d, on a Bspline surface, in the dimensioned frame.
        """

        point2d = self.point3d_to_2d(point3d)

        point2d_with_dimension = self.point2d_parametric_to_dimension(point2d, grid2d)

        return point2d_with_dimension

    def point2d_with_dimension_to_parametric_frame(self, point2d, grid2d: volmdlr.grid.Grid2D):
        """
        Convert a point 2d from the dimensioned to the parametric frame.

        """

        if self._grids2d != grid2d:
            self._grids2d = grid2d
        if not self._grids2d_deformed:
            self.grid2d_deformed(grid2d)

        points_2d = grid2d.points
        points_2d_deformed = self._grids2d_deformed.points
        points_x, points_y = grid2d.points_xy

        # Parameters
        index_points = {}  # grid point position(j,i), point position in points_2d (or points_3d)
        p_index = 0
        for i in range(0, points_x):
            for j in range(0, points_y):
                index_points.update({(j, i): p_index})
                p_index = p_index + 1

        finite_elements_points = []  # 2D grid points index that define one element
        for j in range(0, points_y - 1):
            for i in range(0, points_x - 1):
                finite_elements_points.append(((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)))
        finite_elements = []  # finite elements defined with closed polygon  DEFORMED
        for point in finite_elements_points:
            finite_elements.append(
                volmdlr.wires.ClosedPolygon2D((points_2d_deformed[index_points[point[0]]],
                                               points_2d_deformed[index_points[point[1]]],
                                               points_2d_deformed[index_points[point[2]]],
                                               points_2d_deformed[index_points[point[3]]])))

        finite_elements_initial = []  # finite elements defined with closed polygon  INITIAL
        for point in finite_elements_points:
            finite_elements_initial.append(
                volmdlr.wires.ClosedPolygon2D((points_2d[index_points[point[0]]],
                                               points_2d[index_points[point[1]]],
                                               points_2d[index_points[point[2]]],
                                               points_2d[index_points[point[3]]])))
        k = 0
        for k, point in enumerate(finite_elements_points):
            if (finite_elements[k].point_belongs(point2d)
                    or ((points_2d_deformed[index_points[point[0]]][0] < point2d.x <
                         points_2d_deformed[index_points[point[1]]][0])
                        and point2d.y == points_2d_deformed[index_points[point[0]]][1])
                    or ((points_2d_deformed[index_points[finite_elements_points[k][1]]][1] < point2d.y <
                         points_2d_deformed[index_points[finite_elements_points[k][2]]][1])
                        and point2d.x == points_2d_deformed[index_points[point[1]]][0])
                    or ((points_2d_deformed[index_points[point[3]]][0] < point2d.x <
                         points_2d_deformed[index_points[point[2]]][0])
                        and point2d.y == points_2d_deformed[index_points[point[1]]][1])
                    or ((points_2d_deformed[index_points[point[0]]][1] < point2d.y <
                         points_2d_deformed[index_points[point[3]]][1])
                        and point2d.x == points_2d_deformed[index_points[point[0]]][0])
                    or finite_elements[k].primitives[0].point_belongs(point2d) or finite_elements[k].primitives[
                        1].point_belongs(point2d)
                    or finite_elements[k].primitives[2].point_belongs(point2d) or finite_elements[k].primitives[
                        3].point_belongs(point2d)):
                break

        frame_deformed = volmdlr.Frame2D(finite_elements[k].center_of_mass(),
                                         volmdlr.Vector2D(finite_elements[k].primitives[1].middle_point()[0] -
                                                          finite_elements[k].center_of_mass()[0],
                                                          finite_elements[k].primitives[1].middle_point()[1] -
                                                          finite_elements[k].center_of_mass()[1]),
                                         volmdlr.Vector2D(finite_elements[k].primitives[0].middle_point()[0] -
                                                          finite_elements[k].center_of_mass()[0],
                                                          finite_elements[k].primitives[0].middle_point()[1] -
                                                          finite_elements[k].center_of_mass()[1]))

        point2d_frame_deformed = volmdlr.Point2D(point2d.frame_mapping(frame_deformed, 'new')[0],
                                                 point2d.frame_mapping(frame_deformed, 'new')[1])

        frame_inital = volmdlr.Frame2D(finite_elements_initial[k].center_of_mass(),
                                       volmdlr.Vector2D(finite_elements_initial[k].primitives[1].middle_point()[0] -
                                                        finite_elements_initial[k].center_of_mass()[0],
                                                        finite_elements_initial[k].primitives[1].middle_point()[1] -
                                                        finite_elements_initial[k].center_of_mass()[1]),
                                       volmdlr.Vector2D(finite_elements_initial[k].primitives[0].middle_point()[0] -
                                                        finite_elements_initial[k].center_of_mass()[0],
                                                        finite_elements_initial[k].primitives[0].middle_point()[1] -
                                                        finite_elements_initial[k].center_of_mass()[1]))

        point2d = point2d_frame_deformed.frame_mapping(frame_inital, 'old')
        if point2d.x < 0:
            point2d.x = 0
        elif point2d.x > 1:
            point2d.x = 1
        if point2d.y < 0:
            point2d.y = 0
        elif point2d.y > 1:
            point2d.y = 1

        return point2d

    def point2d_with_dimension_to_3d(self, point2d, grid2d: volmdlr.grid.Grid2D):
        """
        Compute the point 3d, on a Bspline surface, of a point 2d define in the dimensioned frame.

        """

        point2d_01 = self.point2d_with_dimension_to_parametric_frame(point2d, grid2d)

        return self.point2d_to_3d(point2d_01)

    def linesegment2d_parametric_to_dimension(self, linesegment2d, grid2d: volmdlr.grid.Grid2D):
        """
        Convert a linesegment2d from the parametric to the dimensioned frame.

        """

        points = linesegment2d.discretization_points(number_points=20)
        points_dim = [
            self.point2d_parametric_to_dimension(
                point, grid2d) for point in points]

        return vme.BSplineCurve2D.from_points_interpolation(
            points_dim, max(self.degree_u, self.degree_v))

    def linesegment3d_to_2d_with_dimension(self, linesegment3d, grid2d: volmdlr.grid.Grid2D):
        """
        Compute the linesegment2d of a linesegment3d, on a Bspline surface, in the dimensioned frame.

        """

        linesegment2d = self.linesegment3d_to_2d(linesegment3d)
        bsplinecurve2d_with_dimension = self.linesegment2d_parametric_to_dimension(linesegment2d, grid2d)

        return bsplinecurve2d_with_dimension

    def linesegment2d_with_dimension_to_parametric_frame(self, linesegment2d):
        """
        Convert a linesegment2d from the dimensioned to the parametric frame.

        """

        try:
            linesegment2d = vme.LineSegment2D(
                self.point2d_with_dimension_to_parametric_frame(linesegment2d.start, self._grids2d),
                self.point2d_with_dimension_to_parametric_frame(linesegment2d.end, self._grids2d))
        except NotImplementedError:
            return None

        return linesegment2d

    def linesegment2d_with_dimension_to_3d(self, linesegment2d):
        """
        Compute the linesegment3d, on a Bspline surface, of a linesegment2d defined in the dimensioned frame.

        """

        linesegment2d_01 = self.linesegment2d_with_dimension_to_parametric_frame(linesegment2d)
        linesegment3d = self.linesegment2d_to_3d(linesegment2d_01)

        return linesegment3d

    def bsplinecurve2d_parametric_to_dimension(self, bsplinecurve2d, grid2d: volmdlr.grid.Grid2D):
        """
        Convert a bsplinecurve2d from the parametric to the dimensioned frame.

        """

        # check if bsplinecurve2d is in a list
        if isinstance(bsplinecurve2d, list):
            bsplinecurve2d = bsplinecurve2d[0]
        points = bsplinecurve2d.control_points
        points_dim = []

        for point in points:
            points_dim.append(self.point2d_parametric_to_dimension(point, grid2d))

        bsplinecurve2d_with_dimension = vme.BSplineCurve2D(bsplinecurve2d.degree, points_dim,
                                                           bsplinecurve2d.knot_multiplicities,
                                                           bsplinecurve2d.knots,
                                                           bsplinecurve2d.weights,
                                                           bsplinecurve2d.periodic)

        return bsplinecurve2d_with_dimension

    def bsplinecurve3d_to_2d_with_dimension(self, bsplinecurve3d, grid2d: volmdlr.grid.Grid2D):
        """
        Compute the bsplinecurve2d of a bsplinecurve3d, on a Bspline surface, in the dimensioned frame.

        """

        bsplinecurve2d_01 = self.bsplinecurve3d_to_2d(bsplinecurve3d)
        bsplinecurve2d_with_dimension = self.bsplinecurve2d_parametric_to_dimension(
            bsplinecurve2d_01, grid2d)

        return bsplinecurve2d_with_dimension

    def bsplinecurve2d_with_dimension_to_parametric_frame(self, bsplinecurve2d):
        """
        Convert a bsplinecurve2d from the dimensioned to the parametric frame.

        """

        points_dim = bsplinecurve2d.control_points
        points = []
        for point in points_dim:
            points.append(
                self.point2d_with_dimension_to_parametric_frame(point, self._grids2d))

        bsplinecurve2d = vme.BSplineCurve2D(bsplinecurve2d.degree, points,
                                            bsplinecurve2d.knot_multiplicities,
                                            bsplinecurve2d.knots,
                                            bsplinecurve2d.weights,
                                            bsplinecurve2d.periodic)
        return bsplinecurve2d

    def bsplinecurve2d_with_dimension_to_3d(self, bsplinecurve2d):
        """
        Compute the bsplinecurve3d, on a Bspline surface, of a bsplinecurve2d defined in the dimensioned frame.

        """

        bsplinecurve2d_01 = self.bsplinecurve2d_with_dimension_to_parametric_frame(bsplinecurve2d)
        bsplinecurve3d = self.bsplinecurve2d_to_3d(bsplinecurve2d_01)

        return bsplinecurve3d

    def arc2d_parametric_to_dimension(self, arc2d, grid2d: volmdlr.grid.Grid2D):
        """
        Convert an arc 2d from the parametric to the dimensioned frame.

        """

        number_points = math.ceil(arc2d.angle * 7) + 1
        length = arc2d.length()
        points = [self.point2d_parametric_to_dimension(arc2d.point_at_abscissa(
            i * length / (number_points - 1)), grid2d) for i in range(number_points)]

        return vme.BSplineCurve2D.from_points_interpolation(
            points, max(self.degree_u, self.degree_v))

    def arc3d_to_2d_with_dimension(self, arc3d, grid2d: volmdlr.grid.Grid2D):
        """
        Compute the arc 2d of a arc 3d, on a Bspline surface, in the dimensioned frame.

        """

        bsplinecurve2d = self.arc3d_to_2d(arc3d)[0]  # it's a bsplinecurve2d
        arc2d_with_dimension = self.bsplinecurve2d_parametric_to_dimension(bsplinecurve2d, grid2d)

        return arc2d_with_dimension  # it's a bsplinecurve2d-dimension

    def arc2d_with_dimension_to_parametric_frame(self, arc2d):
        """
        Convert an arc 2d from the dimensioned to the parametric frame.

        """

        number_points = math.ceil(arc2d.angle * 7) + 1
        length = arc2d.length()

        points = [self.point2d_with_dimension_to_parametric_frame(arc2d.point_at_abscissa(
            i * length / (number_points - 1)), self._grids2d) for i in range(number_points)]

        return vme.BSplineCurve2D.from_points_interpolation(points, max(self.degree_u, self.degree_v))

    def arc2d_with_dimension_to_3d(self, arc2d):
        """
        Compute the arc 3d, on a Bspline surface, of an arc 2d in the dimensioned frame.

        """

        arc2d_01 = self.arc2d_with_dimension_to_parametric_frame(arc2d)
        arc3d = self.arc2d_to_3d(arc2d_01)

        return arc3d  # it's a bsplinecurve3d

    def contour2d_parametric_to_dimension(self, contour2d: volmdlr.wires.Contour2D,
                                          grid2d: volmdlr.grid.Grid2D):
        """
        Convert a contour 2d from the parametric to the dimensioned frame.

        """

        primitives2d_dim = []

        for primitive2d in contour2d.primitives:
            method_name = f'{primitive2d.__class__.__name__.lower()}_parametric_to_dimension'

            if hasattr(self, method_name):
                primitives = getattr(self, method_name)(primitive2d, grid2d)
                if primitives:
                    primitives2d_dim.append(primitives)

            else:
                raise NotImplementedError(
                    f'Class {self.__class__.__name__} does not implement {method_name}')

        return volmdlr.wires.Contour2D(primitives2d_dim)

    def contour3d_to_2d_with_dimension(self, contour3d: volmdlr.wires.Contour3D,
                                       grid2d: volmdlr.grid.Grid2D):
        """
        Compute the Contour 2d of a Contour 3d, on a Bspline surface, in the dimensioned frame.

        """

        contour2d_01 = self.contour3d_to_2d(contour3d)

        return self.contour2d_parametric_to_dimension(contour2d_01, grid2d)

    def contour2d_with_dimension_to_parametric_frame(self, contour2d):
        """
        Convert a contour2d from the dimensioned to the parametric frame.

        """

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
        """
        Compute the contour3d, on a Bspline surface, of a contour2d define in the dimensioned frame.

        """

        contour01 = self.contour2d_with_dimension_to_parametric_frame(contour2d)

        return self.contour2d_to_3d(contour01)

    @classmethod
    def from_geomdl_surface(cls, surface):
        """
        Create a volmdlr BSpline_Surface3D from a geomdl's one.

        """

        control_points = []
        for point in surface.ctrlpts:
            control_points.append(volmdlr.Point3D(point[0], point[1], point[2]))

        (u_knots, u_multiplicities) = knots_vector_inv(surface.knotvector_u)
        (v_knots, v_multiplicities) = knots_vector_inv(surface.knotvector_v)

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
        """
        Bspline Surface interpolation through 3d points.

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

        """

        points = []
        for point in points_3d:
            points.append((point.x, point.y, point.z))

        surface = interpolate_surface(points, size_u, size_v, degree_u, degree_v)

        return cls.from_geomdl_surface(surface)

    @classmethod
    def points_approximate_into_bspline_surface(cls, points_3d, size_u, size_v, degree_u, degree_v, **kwargs):
        """
        Bspline Surface approximate through 3d points.

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

        """

        # Keyword arguments
        num_cpts_u = kwargs.get('ctrlpts_size_u', size_u - 1)  # number of datapts, r + 1 > number of ctrlpts, n + 1
        num_cpts_v = kwargs.get('ctrlpts_size_v', size_v - 1)  # number of datapts, s + 1 > number of ctrlpts, m + 1

        points = [tuple([*point]) for point in points_3d]

        surface = approximate_surface(points, size_u, size_v, degree_u, degree_v,
                                      ctrlpts_size_u=num_cpts_u, num_cpts_v=num_cpts_v)

        return cls.from_geomdl_surface(surface)

    @classmethod
    def from_cylindrical_faces(cls, cylindrical_faces, degree_u, degree_v,
                               points_x: int = 10, points_y: int = 10):
        """
        Define a bspline surface from a list of cylindrical faces.

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

        """

        if len(cylindrical_faces) == 1:
            return cls.from_cylindrical_face(cylindrical_faces[0], degree_u, degree_v, points_x=50, points_y=50)

        if len(cylindrical_faces) > 1:
            bspline_surfaces = []
            direction = cylindrical_faces[0].adjacent_direction(cylindrical_faces[1])

            if direction == 'x':
                bounding_rectangle_0 = cylindrical_faces[0].surface2d.outer_contour.bounding_rectangle
                ymin = bounding_rectangle_0[2]
                ymax = bounding_rectangle_0[3]
                for face in cylindrical_faces:
                    bounding_rectangle = face.surface2d.outer_contour.bounding_rectangle
                    ymin = min(ymin, bounding_rectangle[2])
                    ymax = max(ymax, bounding_rectangle[3])
                for face in cylindrical_faces:
                    bounding_rectangle = face.surface2d.outer_contour.bounding_rectangle

                    points_3d = face.surface3d.grid3d(
                        volmdlr.grid.Grid2D.from_properties(
                            x_limits=(bounding_rectangle[0], bounding_rectangle[1]),
                            y_limits=(ymin, ymax),
                            points_nbr=(points_x, points_y)))

                    bspline_surfaces.append(
                        cls.points_fitting_into_bspline_surface(
                            points_3d, points_x, points_y, degree_u, degree_v))

            elif direction == 'y':
                bounding_rectangle_0 = cylindrical_faces[0].surface2d.outer_contour.bounding_rectangle
                xmin = bounding_rectangle_0[0]
                xmax = bounding_rectangle_0[1]
                for face in cylindrical_faces:
                    bounding_rectangle = face.surface2d.outer_contour.bounding_rectangle
                    xmin = min(xmin, bounding_rectangle[0])
                    xmax = max(xmax, bounding_rectangle[1])
                for face in cylindrical_faces:
                    bounding_rectangle = face.surface2d.outer_contour.bounding_rectangle

                    points_3d = face.surface3d.grid3d(
                        volmdlr.grid.Grid2D.from_properties(
                            x_limits=(xmin, xmax),
                            y_limits=(bounding_rectangle[2], bounding_rectangle[3]),
                            points_nbr=(points_x, points_y)))

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
    def from_cylindrical_face(cls, cylindrical_face, degree_u, degree_v,
                              **kwargs):  # points_x: int = 50, points_y: int = 50
        """
        Define a bspline surface from a cylindrical face.

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

        """

        points_x = kwargs['points_x']
        points_y = kwargs['points_y']
        bounding_rectangle = cylindrical_face.surface2d.outer_contour.bounding_rectangle
        points_3d = cylindrical_face.surface3d.grid3d(
            volmdlr.grid.Grid2D.from_properties(x_limits=(bounding_rectangle[0],
                                                          bounding_rectangle[1]),
                                                y_limits=(bounding_rectangle[2],
                                                          bounding_rectangle[3]),
                                                points_nbr=(points_x, points_y)))

        return cls.points_fitting_into_bspline_surface(points_3d, points_x, points_x, degree_u, degree_v)

    def intersection_with(self, other_bspline_surface3d):
        """
        Compute intersection points between two Bspline surfaces.

        return u,v parameters for intersection points for both surfaces
        """

        def fun(param):
            return (self.point2d_to_3d(volmdlr.Point2D(param[0], param[1])) -
                    other_bspline_surface3d.point2d_to_3d(volmdlr.Point2D(param[2], param[3]))).norm()

        x = npy.linspace(0, 1, 10)
        x_init = []
        for xi in x:
            for yi in x:
                x_init.append((xi, yi, xi, yi))

        u1, v1, u2, v2 = [], [], [], []
        solutions = []
        for x0 in x_init:
            z = least_squares(fun, x0=x0, bounds=([0, 1]))
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

        return (u1, v1), (u2, v2)  # (uv1, uv2)

    def plane_intersection(self, plane3d):
        """
        Compute intersection points between a Bspline surface and a plane 3d.

        """

        def f(param):
            return ((self.surface.evaluate_single((param[0], param[1]))[0]) * plane3d.equation_coefficients()[0] +
                    (self.surface.evaluate_single((param[0], param[1]))[1]) * plane3d.equation_coefficients()[1] +
                    (self.surface.evaluate_single((param[0], param[1]))[2]) * plane3d.equation_coefficients()[2] +
                    plane3d.equation_coefficients()[3])

        x = npy.linspace(0, 1, 20)
        x_init = []
        for xi in x:
            for yi in x:
                x_init.append((xi, yi))

        intersection_points = []

        for x0 in x_init:
            z = least_squares(f, x0=x0, bounds=([0, 1]))
            if z.fun < 1e-20:
                solution = z.x
                intersection_points.append(volmdlr.Point3D(self.surface.evaluate_single((solution[0], solution[1]))[0],
                                                           self.surface.evaluate_single((solution[0], solution[1]))[1],
                                                           self.surface.evaluate_single((solution[0], solution[1]))[
                                                               2]))
        return intersection_points

    def error_with_point3d(self, point3d):
        """
        Compute the error/distance between the Bspline surface and a point 3d.

        """

        def f(x):
            return (point3d - self.point2d_to_3d(volmdlr.Point2D(x[0], x[1]))).norm()

        cost = []

        for x0 in [(0, 0), (0, 1), (1, 0), (1, 1), (0.5, 0.5)]:
            z = least_squares(f, x0=x0, bounds=([0, 1]))
            cost.append(z.fun)

        return min(cost)

    def error_with_edge3d(self, edge3d):
        """
        Compute the error/distance between the Bspline surface and an edge 3d.

        it's the mean of the start and end points errors'
        """

        return (self.error_with_point3d(edge3d.start) + self.error_with_point3d(edge3d.end)) / 2

    def nearest_edges3d(self, contour3d, threshold: float):
        """
        Compute the nearest edges of a contour 3d to a Bspline_surface3d based on a threshold.

        """

        nearest = []
        for primitive in contour3d.primitives:
            if self.error_with_edge3d(primitive) <= threshold:
                nearest.append(primitive)
        nearest_primitives = volmdlr.wires.Wire3D(nearest)

        return nearest_primitives

    def edge3d_to_2d_with_dimension(self, edge3d, grid2d: volmdlr.grid.Grid2D):
        """
        Compute the edge 2d of a edge 3d, on a Bspline surface, in the dimensioned frame.

        """
        method_name = f'{edge3d.__class__.__name__.lower()}_to_2d_with_dimension'

        if hasattr(self, method_name):
            edge2d_dim = getattr(self, method_name)(edge3d, grid2d)
            if edge2d_dim:
                return edge2d_dim
            raise NotImplementedError
        raise NotImplementedError(
            f'Class {self.__class__.__name__} does not implement {method_name}')

    def wire3d_to_2d(self, wire3d):
        """
        Compute the 2d of a wire 3d, on a Bspline surface.

        """

        contour = self.contour3d_to_2d(wire3d)

        return volmdlr.wires.Wire2D(contour.primitives)

    def wire3d_to_2d_with_dimension(self, wire3d):
        """
        Compute the 2d of a wire 3d, on a Bspline surface, in the dimensioned frame.

        """

        contour = self.contour3d_to_2d_with_dimension(wire3d, self._grids2d)

        return volmdlr.wires.Wire2D(contour.primitives)

    def split_surface_u(self, u: float):
        """
        Splits the surface at the input parametric coordinate on the u-direction.

        :param u: Parametric coordinate u chosen between 0 and 1
        :type u: float
        :return: Two split surfaces
        :rtype: List[:class:`volmdlr.faces.BSplineSurface3D`]
        """

        surfaces_geo = split_surface_u(self.surface, u)
        surfaces = [BSplineSurface3D.from_geomdl_surface(surface) for surface in surfaces_geo]
        return surfaces

    def split_surface_v(self, v: float):
        """
        Splits the surface at the input parametric coordinate on the v-direction.

        :param v: Parametric coordinate v chosen between 0 and 1
        :type v: float
        :return: Two split surfaces
        :rtype: List[:class:`volmdlr.faces.BSplineSurface3D`]
        """

        surfaces_geo = split_surface_v(self.surface, v)
        surfaces = [BSplineSurface3D.from_geomdl_surface(surface) for surface in surfaces_geo]
        return surfaces

    def split_surface_with_bspline_curve(self, bspline_curve3d: vme.BSplineCurve3D):
        """
        Cuts the surface into two pieces with a bspline curve.

        :param bspline_curve3d: A BSplineCurve3d used for cutting
        :type bspline_curve3d: :class:`vme.BSplineCurve3D`
        :return: Two split surfaces
        :rtype: List[:class:`volmdlr.faces.BSplineSurface3D`]
        """

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
            u_min, u_max, v_min, v_max = contour.bounding_rectangle.bounds()
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

            for line in lines:
                inter = contour.line_intersections(line)
                if inter:
                    pt_ = set()
                    for point_intersection in inter:
                        pt_.add(point_intersection[0])
                else:
                    raise NotImplementedError

                pt_ = sorted(pt_, key=pt0.point_distance)
                pt0 = pt_[0]
                edge = vme.LineSegment2D(pt_[0], pt_[1])

                points.extend(edge.discretization_points(number_points=10))

            points3d = []
            for point in points:
                points3d.append(self.point2d_to_3d(point))

            size_u, size_v, degree_u, degree_v = 10, 10, self.degree_u, self.degree_v
            surfaces.append(
                volmdlr.faces.BSplineSurface3D.points_fitting_into_bspline_surface(
                    points3d, size_u, size_v, degree_u, degree_v))

        return surfaces

    def point_belongs(self, point3d):
        """
        Check if a point 3d belongs to the bspline_surface or not.

        """

        def f(param):
            p3d = self.point2d_to_3d(volmdlr.Point2D(param[0], param[1]))
            return point3d.point_distance(p3d)

        x = npy.linspace(0, 1, 5)
        x_init = []
        for xi in x:
            for yi in x:
                x_init.append((xi, yi))

        for x0 in x_init:
            z = least_squares(f, x0=x0, bounds=([0, 1]))
            if z.fun < 1e-10:
                return True
        return False

    def is_intersected_with(self, other_bspline_surface3d):
        """
        Check if the two surfaces are intersected or not.

        return True, when there are more 50points on the intersection zone.

        """

        # intersection_results = self.intersection_with(other_bspline_surface3d)
        # if len(intersection_results[0][0]) >= 50:
        #     return True
        # else:
        #     return False

        def f(param):
            return (self.point2d_to_3d(volmdlr.Point2D(param[0], param[1])) -
                    other_bspline_surface3d.point2d_to_3d(volmdlr.Point2D(param[2], param[3]))).norm()

        x = npy.linspace(0, 1, 10)
        x_init = []
        for xi in x:
            for yi in x:
                x_init.append((xi, yi, xi, yi))

        i = 0
        for x0 in x_init:
            z = least_squares(f, x0=x0, bounds=([0, 1]))
            if z.fun < 1e-5:
                i += 1
                if i >= 50:
                    return True
        return False

    def merge_with(self, other_bspline_surface3d):
        """
        Merges two adjacent surfaces based on their faces.

        :param other_bspline_surface3d: Other adjacent surface
        :type other_bspline_surface3d: :class:`volmdlr.faces.BSplineSurface3D`

        :return: Merged surface
        :rtype: :class:`volmdlr.faces.BSplineSurface3D`
        """

        bspline_face3d = self.rectangular_cut(0, 1, 0, 1)
        other_bspline_face3d = other_bspline_surface3d.rectangular_cut(0, 1, 0, 1)

        bsplines = [self, other_bspline_surface3d]
        bsplines_new = bsplines

        center = [bspline_face3d.surface2d.outer_contour.center_of_mass(),
                  other_bspline_face3d.surface2d.outer_contour.center_of_mass()]
        grid2d_direction = (bspline_face3d.pair_with(other_bspline_face3d))[1]

        if (not bspline_face3d.outer_contour3d.is_sharing_primitives_with(other_bspline_face3d.outer_contour3d)
            and self.is_intersected_with(other_bspline_surface3d)):
            # find pimitives to split with
            contour1 = bspline_face3d.outer_contour3d
            contour2 = other_bspline_face3d.outer_contour3d

            distances = []
            for prim1 in contour1.primitives:
                dis = []
                for prim2 in contour2.primitives:
                    point1 = (prim1.start + prim1.end) / 2
                    point2 = (prim2.start + prim2.end) / 2
                    dis.append(point1.point_distance(point2))
                distances.append(dis)

            i = distances.index((min(distances)))
            j = distances[i].index(min(distances[i]))

            curves = [contour2.primitives[j], contour1.primitives[i]]

            # split surface
            for i, bspline in enumerate(bsplines):
                surfaces = bspline.split_surface_with_bspline_curve(curves[i])

                errors = []
                for surface in surfaces:
                    errors.append(surface.error_with_point3d(bsplines[i].point2d_to_3d(center[i])))

                bsplines_new[i] = surfaces[errors.index(min(errors))]

            grid2d_direction = (
                bsplines_new[0].rectangular_cut(
                    0, 1, 0, 1).pair_with(
                    bsplines_new[1].rectangular_cut(
                        0, 1, 0, 1)))[1]

        # grid3d
        nb = 10
        points3d = []
        is_true = (bspline_face3d.outer_contour3d.is_sharing_primitives_with(other_bspline_face3d.outer_contour3d)
                    or self.is_intersected_with(other_bspline_surface3d))

        for i, bspline in enumerate(bsplines_new):
            grid3d = bspline.grid3d(volmdlr.grid.Grid2D.from_properties(x_limits=(0, 1),
                                                                        y_limits=(0, 1),
                                                                        points_nbr=(nb, nb),
                                                                        direction=grid2d_direction[i]))

            if is_true and i == 1:
                points3d.extend(grid3d[nb:nb * nb])
            else:
                points3d.extend(grid3d)

        # fitting
        size_u, size_v, degree_u, degree_v = (nb * 2) - 1, nb, 3, 3

        merged_surface = BSplineSurface3D.points_fitting_into_bspline_surface(
            points3d, size_u, size_v, degree_u, degree_v)

        return merged_surface

    def xy_limits(self, other_bspline_surface3d):
        """
        Compute x, y limits to define grid2d.

        """

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

    def derivatives(self, u, v, order):
        """
        Evaluates n-th order surface derivatives at the given (u, v) parameter pair.

        :param u: Point's u coordinate.
        :type u: float
        :param v: Point's v coordinate.
        :type v: float
        :param order: Order of the derivatives.
        :type order: int
        :return: A list SKL, where SKL[k][l] is the derivative of the surface S(u,v) with respect
        to u k times and v l times
        :rtype: List[`volmdlr.Vector3D`]
        """
        if self.surface.rational:
            # derivatives = self._rational_derivatives(self.surface.data,(u, v), order)
            derivatives = volmdlr.bspline_compiled.rational_derivatives(self.surface.data, (u, v), order)
        else:
            # derivatives = self._derivatives(self.surface.data, (u, v), order)
            derivatives = volmdlr.bspline_compiled.derivatives(self.surface.data, (u, v), order)
        for i in range(order + 1):
            for j in range(order + 1):
                derivatives[i][j] = volmdlr.Vector3D(*derivatives[i][j])
        return derivatives


class BezierSurface3D(BSplineSurface3D):
    """
    A 3D Bezier surface.

    :param degree_u: The degree of the Bezier surface in the u-direction.
    :type degree_u: int
    :param degree_v: The degree of the Bezier surface in the v-direction.
    :type degree_v: int
    :param control_points: A list of lists of control points defining the Bezier surface.
    :type control_points: List[List[`volmdlr.Point3D`]]
    :param nb_u: The number of control points in the u-direction.
    :type nb_u: int
    :param nb_v: The number of control points in the v-direction.
    :type nb_v: int
    :param name: (Optional) name for the Bezier surface.
    :type name: str
    """

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
    """
    Abstract method to define 3D faces.
    """

    min_x_density = 1
    min_y_density = 1

    def __init__(self, surface3d, surface2d: Surface2D,
                 name: str = ''):
        self.surface3d = surface3d
        self.surface2d = surface2d
        self._outer_contour3d = None
        self._inner_contours3d = None
        # self.bounding_box = self._bounding_box()

        volmdlr.core.Primitive3D.__init__(self, name=name)

    def to_dict(self, *args, **kwargs):
        """Avoids storing points in memo that makes serialization slow."""
        return DessiaObject.to_dict(self, use_pointers=False)

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
        Tells you if a point is on the 3D face and inside its contour.
        """
        point2d = self.surface3d.point3d_to_2d(point3d)
        check_point3d = self.surface3d.point2d_to_3d(point2d)
        if check_point3d.point_distance(point3d) > 1e-6:
            return False

        return self.surface2d.point_belongs(point2d)

    @property
    def outer_contour3d(self):
        """
        Gives the 3d version of the outer contour of the face.
        """
        if not self._outer_contour3d:
            self._outer_contour3d = self.surface3d.contour2d_to_3d(self.surface2d.outer_contour)
        return self._outer_contour3d

    @property
    def inner_contours3d(self):
        """
        Gives the 3d version of the inner contours of the face.
        """
        if not self._inner_contours3d:
            self._inner_contours3d = [self.surface3d.contour2d_to_3d(c) for c in
                                      self.surface2d.inner_contours]
        return self._inner_contours3d

    @property
    def bounding_box(self):
        """
        Needs to be overridden if an error is raised.
        """
        raise NotImplementedError(
            f"bounding_box method must be"
            f"overloaded by {self.__class__.__name__}")

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        """
        Sets the bounding box to a new value.
        """
        raise NotImplementedError(
            f"bounding_box setter method must be"
            f"overloaded by {self.__class__.__name__}")

    def get_bounding_box(self):
        raise NotImplementedError(
            f"self.__class__.__name__"
            f"overloaded by {self.__class__.__name__}")

    def area(self):
        return self.surface2d.area()

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a Face3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding Face3D object.
        :rtype: :class:`volmdlr.faces.Face3D`
        """
        step_id = kwargs.get("step_id", "#UNKNOW_ID")
        step_name = kwargs.get("name", "ADVANCED_FACE")
        name = arguments[0][1:-1]
        contours = [object_dict[int(arg[1:])] for arg in arguments[1]]
        if any(contour is None for contour in contours):
            warnings.warn(f"Could not instantiate #{step_id} = {step_name}({arguments})"
                          f" because some of the contours are NoneType."
                          "See Face3D.from_step method")
            return None
        surface = object_dict[int(arguments[2])]
        if hasattr(surface, 'face_from_contours3d'):
            if (len(contours) == 1) and isinstance(contours[0], volmdlr.Point3D):
                return surface
            if (len(contours) == 2) and isinstance(contours[1], volmdlr.Point3D):
                raise NotImplementedError
            return surface.face_from_contours3d(contours, name)
            # except Exception:
            #     return None
        raise NotImplementedError(
            f'Not implemented :face_from_contours3d in {surface}')

    def to_step(self, current_id):
        xmin, xmax, ymin, ymax = self.surface2d.bounding_rectangle().bounds()
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
                face_content, face_id = face.to_step_without_splitting(
                    current_id)
                face_ids.append(face_id[0])
                content += face_content
                current_id = face_id[0] + 1
            return content, face_ids
        return self.to_step_without_splitting(current_id)

    def to_step_without_splitting(self, current_id):
        content, surface3d_ids = self.surface3d.to_step(current_id)
        current_id = max(surface3d_ids) + 1

        if len(surface3d_ids) != 1:
            raise NotImplementedError('What to do with more than 1 id ? with 0 id ?')
        outer_contour_content, outer_contour_id = self.outer_contour3d.to_step(
            current_id, surface_id=surface3d_ids[0], surface3d=self.surface3d)
        content += outer_contour_content
        content += f"#{outer_contour_id + 1} = FACE_BOUND('{self.name}',#{outer_contour_id},.T.);\n"
        contours_ids = [outer_contour_id + 1]
        current_id = outer_contour_id + 2
        for inner_contour3d in self.inner_contours3d:
            inner_contour_content, inner_contour_id = inner_contour3d.to_step(
                current_id)
            # surface_id=surface3d_id)
            content += inner_contour_content
            face_bound_id = inner_contour_id + 1
            content += f"#{face_bound_id} = FACE_BOUND('',#{inner_contour_id},.T.);\n"
            contours_ids.append(face_bound_id)
            current_id = face_bound_id + 1

        content += f"#{current_id} = ADVANCED_FACE('{self.name}',({volmdlr.core.step_ids_to_str(contours_ids)})" \
                   f",#{surface3d_ids[0]},.T.);\n"
        # TODO: create an ADVANCED_FACE for each surface3d_ids ?
        return content, [current_id]

    def triangulation_lines(self):
        return [], []

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        return [0, 0]

    def triangulation(self):
        number_points_x, number_points_y = self.grid_size()
        mesh2d = self.surface2d.triangulation(number_points_x, number_points_y)
        return vmd.DisplayMesh3D([vmd.Node3D(*self.surface3d.point2d_to_3d(point)) for point in mesh2d.points],
                                 mesh2d.triangles)

    def plot2d(self, ax=None, color='k', alpha=1):
        if ax is None:
            _, ax = plt.subplots()
        self.surface2d.plot(ax=ax, color=color, alpha=alpha)

    def rotation(self, center: volmdlr.Point3D,
                 axis: volmdlr.Vector3D, angle: float):
        """
        Face3D rotation.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: angle rotation
        :return: a new rotated Face3D
        """
        new_surface = self.surface3d.rotation(center=center, axis=axis,
                                              angle=angle)
        return self.__class__(new_surface, self.surface2d)

    def rotation_inplace(self, center: volmdlr.Point3D,
                         axis: volmdlr.Vector3D, angle: float):
        """
        Face3D rotation.

         Object is updated inplace.
        :param center: rotation center.
        :type center: `volmdlr.Point3D`
        :param axis: rotation axis.
        :type axis: `volmdlr.Vector3D`
        :param angle: rotation angle.
        :type angle: float
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.surface3d.rotation_inplace(center=center, axis=axis, angle=angle)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box

    def translation(self, offset: volmdlr.Vector3D):
        """
        Face3D translation.

        :param offset: Translation vector.
        :type offset: `volmdlr.Vector3D`
        :return: A new translated Face3D
        """
        new_surface3d = self.surface3d.translation(offset=offset)
        return self.__class__(new_surface3d, self.surface2d)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Face3D translation. Object is updated inplace.

        :param offset: translation vector
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.surface3d.translation_inplace(offset=offset)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Face3D.

        side = 'old' or 'new'
        """
        new_surface3d = self.surface3d.frame_mapping(frame, side)
        return self.__class__(new_surface3d, self.surface2d.copy(),
                              self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        side = 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.surface3d.frame_mapping_inplace(frame, side)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box

    def copy(self, deep=True, memo=None):
        return self.__class__(self.surface3d.copy(deep, memo), self.surface2d.copy(),
                              self.name)

    def face_inside(self, face2):
        """
        Verifies if a face is inside another one.

        It returns True if face2 is inside or False if the opposite.
        """
        if self.surface3d.is_coincident(face2.surface3d):
            self_contour2d = self.outer_contour3d.to_2d(
                self.surface3d.frame.origin, self.surface3d.frame.u, self.surface3d.frame.v)
            face2_contour2d = face2.outer_contour3d.to_2d(
                self.surface3d.frame.origin, self.surface3d.frame.u, self.surface3d.frame.v)
            if self_contour2d.is_inside(face2_contour2d):
                return True
        return False

    def edge_intersections(self, edge):
        intersections = []
        method_name = f'{edge.__class__.__name__.lower()[:-2]}_intersections'
        if hasattr(self, method_name):
            intersections = getattr(self, method_name)(edge)
        if not intersections:
            for point in [edge.start, edge.end]:
                if self.point_belongs(point):
                    if point not in intersections:
                        intersections.append(point)
        return intersections

    def line_intersections(self,
                           line: vme.Line3D,
                           ) -> List[volmdlr.Point3D]:
        intersections = []
        for intersection in self.surface3d.line_intersections(line):
            if self.point_belongs(intersection):
                intersections.append(intersection)
        if not intersections:
            for prim in self.outer_contour3d.primitives:
                intersection = prim.line_intersections(line)
                if intersection:
                    if intersection not in intersections:
                        intersections.append(intersection)

        return intersections

    def linesegment_intersections(self, linesegment: vme.LineSegment3D) -> List[volmdlr.Point3D]:
        linesegment_intersections = []
        for intersection in self.surface3d.linesegment_intersections(linesegment):
            if self.point_belongs(intersection):
                linesegment_intersections.append(intersection)
        if not linesegment_intersections:
            for prim in self.outer_contour3d.primitives:
                intersections = prim.linesegment_intersections(linesegment)
                for intersection in intersections:
                    if intersection not in linesegment_intersections:
                        linesegment_intersections.append(intersection)
        return linesegment_intersections

    def fullarc_intersections(self, fullarc: vme.FullArc3D) -> List[volmdlr.Point3D]:
        intersections = []
        for intersection in self.surface3d.fullarc_intersections(fullarc):
            if self.point_belongs(intersection):
                intersections.append(intersection)
        return intersections

    def plot(self, ax=None, color='k', alpha=1, edge_details=False):
        if not ax:
            ax = plt.figure().add_subplot(111, projection='3d')
        self.outer_contour3d.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha,
                                                              edge_ends=edge_details, edge_direction=edge_details))
        for contour3d in self.inner_contours3d:
            contour3d.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha,
                                                       edge_ends=edge_details, edge_direction=edge_details))
        return ax

    def random_point_inside(self):
        point_inside2d = self.surface2d.random_point_inside()
        return self.surface3d.point2d_to_3d(point_inside2d)

    def is_adjacent(self, face2: 'Face3D'):
        contour1 = self.outer_contour3d.to_2d(
            self.surface3d.frame.origin,
            self.surface3d.frame.u,
            self.surface3d.frame.v)
        contour2 = face2.outer_contour3d.to_2d(
            self.surface3d.frame.origin,
            self.surface3d.frame.u,
            self.surface3d.frame.v)
        if contour1.is_sharing_primitives_with(contour2):
            return True
        return False

    def geo_lines(self):  # , mesh_size_list=None):
        """
        Gets the lines that define a Face3D in a .geo file.

        """

        i_index, p_index = None, None
        lines, line_surface, lines_tags = [], [], []
        point_account, line_account, line_loop_account = 0, 0, 1
        for c_index, contour in enumerate(list(chain(*[[self.outer_contour3d], self.inner_contours3d]))):

            if isinstance(contour, volmdlr.wires.Circle2D):
                # point=[contour.radius, contour.center.y, 0]
                # lines.append('Point('+str(point_account+1)+') = {'+str(point)[1:-1]+', '+str(mesh_size)+'};')

                # point = [*contour.center, 0]
                # lines.append('Point('+str(point_account+2)+') = {'+str(point)[1:-1]+', '+str(mesh_size)+'};')

                # point=[-contour.radius, contour.center.y, 0]
                # lines.append('Point('+str(point_account+3)+') = {'+str(point)[1:-1]+', '+str(mesh_size)+'};')

                # lines.append('Circle('+str(line_account+1)+') = {'+str(point_account+1)+','+str(point_account+2) \
                #              +','+str(point_account+3)+'};')
                # lines.append('Circle('+str(line_account+2)+') = {'+str(point_account+3)+','+str(point_account+2) \
                #              + ','+str(point_account+1)+'};')

                # lines_tags.extend([line_account+1, line_account+2])

                # lines.append('Line Loop('+str(line_loop_account+1)+') = {'+str(lines_tags)[1:-1]+'};')

                # line_surface.append(line_loop_account+1)

                # lines_tags = []
                # point_account, line_account, line_loop_account = point_account+3, line_account+2, line_loop_account+1

                pass

            elif isinstance(contour, (volmdlr.wires.Contour3D, volmdlr.wires.ClosedPolygon3D)):
                if not isinstance(contour, volmdlr.wires.ClosedPolygon3D):
                    contour = contour.to_polygon(1)
                for i_index, point in enumerate(contour.points):
                    lines.append(point.get_geo_lines(tag=point_account + i_index + 1,
                                                     point_mesh_size=None))

                for p_index, primitive in enumerate(contour.primitives):
                    if p_index != len(contour.primitives) - 1:
                        lines.append(primitive.get_geo_lines(tag=line_account + p_index + 1,
                                                             start_point_tag=point_account + p_index + 1,
                                                             end_point_tag=point_account + p_index + 2))
                    else:
                        lines.append(primitive.get_geo_lines(tag=line_account + p_index + 1,
                                                             start_point_tag=point_account + p_index + 1,
                                                             end_point_tag=point_account + 1))
                    lines_tags.append(line_account + p_index + 1)

                lines.append('Line Loop(' + str(c_index + 1) + ') = {' + str(lines_tags)[1:-1] + '};')
                line_surface.append(line_loop_account)
                point_account = point_account + i_index + 1
                line_account, line_loop_account = line_account + p_index + 1, line_loop_account + 1
                lines_tags = []

        lines.append('Plane Surface(' + str(1) + ') = {' + str(line_surface)[1:-1] + '};')

        return lines

    def to_geo(self, file_name: str):  # , mesh_size_list=None):
        """
        Gets the .geo file for the Face3D.

        """

        lines = self.geo_lines()

        with open(file_name + '.geo', 'w', encoding="utf-8") as f:
            for line in lines:
                f.write(line)
                f.write('\n')
        f.close()

    def get_geo_lines(self, tag: int, line_loop_tag: List[int]):
        """
        Gets the lines that define a PlaneFace3D in a .geo file.

        """

        return 'Plane Surface(' + str(tag) + ') = {' + str(line_loop_tag)[1:-1] + '};'

    def edge3d_inside(self, edge3d):
        """
        Returns True if edge 3d is coplanar to the face.
        """
        method_name = f'{edge3d.__class__.__name__.lower()[:-2]}_inside'
        if hasattr(self, method_name):
            return getattr(self, method_name)(edge3d)
        points = edge3d.discretization_points(number_points=5)
        for point in points[1:-1]:
            if not self.point_belongs(point):
                return False
        return True

    def is_intersecting(self, face2, list_coincident_faces=None, tol: float = 1e-6):
        """
        Verifies if two face are intersecting.

        :param face2: face 2
        :param list_coincident_faces: list of coincident faces, if existent
        :param tol: tolerance for calculations
        :return: True if faces intersect, False otherwise
        """
        if list_coincident_faces is None:
            list_coincident_faces = []
        if (self.bounding_box.bbox_intersection(face2.bounding_box) or
            self.bounding_box.distance_to_bbox(face2.bounding_box) <= tol) and \
                (self, face2) not in list_coincident_faces:

            edge_intersections = []
            for prim1 in self.outer_contour3d.primitives + [prim for inner_contour in self.inner_contours3d
                                                            for prim in inner_contour.primitives]:
                edge_intersections = face2.edge_intersections(prim1)
                if edge_intersections:
                    return True
            if not edge_intersections:
                for prim2 in face2.outer_contour3d.primitives + [prim for inner_contour in face2.inner_contours3d
                                                                 for prim in inner_contour.primitives]:
                    edge_intersections = self.edge_intersections(prim2)
                    if edge_intersections:
                        return True

        return False

    def face_intersections_outer_contour(self, face2):
        """
        Returns the intersections of the face outer contour with other given face.
        """
        intersections_points = []
        for edge1 in self.outer_contour3d.primitives:
            intersection_points = face2.edge_intersections(edge1)
            if intersection_points:
                for point in intersection_points:
                    if point not in intersections_points:
                        intersections_points.append(point)

        return intersections_points

    def face_intersections(self, face2, tol=1e-6) -> List[volmdlr.wires.Wire3D]:
        """
        Calculates the intersections between two Face3D.

        """

        bbox1 = self.bounding_box
        bbox2 = face2.bounding_box
        if not bbox1.bbox_intersection(bbox2) and \
                bbox1.distance_to_bbox(bbox2) >= tol:
            return []
        if self.face_inside(face2) or face2.face_inside(self):
            return []
        face_intersections = self.get_face_intersections(face2)
        return face_intersections

    def get_face_intersections(self, face2):
        """
        Gets the intersections between two faces.

        :param face2: second face.
        :return: intersections.
        """
        method_name = f'{face2.__class__.__name__.lower()[:-2]}_intersections'
        intersections = getattr(self, method_name)(face2)
        return intersections

    def set_operations_new_faces(self, intersecting_combinations, contour_extract_inside):
        self_copy = self.copy(deep=True)
        list_cutting_contours = self_copy.get_face_cutting_contours(intersecting_combinations)
        if not list_cutting_contours:
            return [self_copy]
        return self_copy.divide_face(list_cutting_contours, contour_extract_inside)

    def get_face_cutting_contours(self, dict_intersecting_combinations):
        """
        Get all contours cutting the face, resultig from multiple faces intersections.

        :param dict_intersecting_combinations: dictionary containing as keys the combination of intersecting faces
        and as the values the resulting primitive from the intersection of these two faces
        return a list all contours cutting one particular face.
        """
        face_intersecting_primitives2d = self.select_face_intersecting_primitives(dict_intersecting_combinations)
        if not face_intersecting_primitives2d:
            return []
        list_cutting_contours = volmdlr.wires.Contour2D.contours_from_edges(face_intersecting_primitives2d[:])
        if self.surface2d.inner_contours:
            valid_cutting_contours = []
            connectig_to_outer_contour = []
            for cutting_contour in list_cutting_contours:
                if (self.surface2d.outer_contour.point_over_contour(cutting_contour.primitives[0].start) and
                    self.surface2d.outer_contour.point_over_contour(cutting_contour.primitives[-1].end)) or \
                        cutting_contour.primitives[0].start == cutting_contour.primitives[-1].end:
                    valid_cutting_contours.append(cutting_contour)
                if self.surface2d.outer_contour.contour_intersections(cutting_contour):
                    connectig_to_outer_contour.append(cutting_contour)
            if len(valid_cutting_contours) == len(list_cutting_contours):
                return valid_cutting_contours
            for cutting_contour in valid_cutting_contours:
                list_cutting_contours.remove(cutting_contour)
            new_cutting_contours, cutting_contours = self.get_inner_contours_cutting_primitives(
                list_cutting_contours, connectig_to_outer_contour)
            return valid_cutting_contours + new_cutting_contours + cutting_contours
        return list_cutting_contours

    def divide_face(self, list_cutting_contours: List[volmdlr.wires.Contour2D], inside):
        """
        Divides a Face 3D with a list of cutting contours.

        :param list_cutting_contours: list of contours cutting the face
        :param inside: when extracting a contour from another contour. It defines the extracted
        contour as being between the two points if True and outside these points if False
        return a list new faces resulting from face division
        """
        list_faces = []
        list_open_cutting_contours = []
        list_closed_cutting_contours = []
        for cutting_contour in list_cutting_contours:
            if not cutting_contour.primitives[0].start.is_close(cutting_contour.primitives[-1].end):
                list_open_cutting_contours.append(cutting_contour)
                continue
            list_closed_cutting_contours.append(cutting_contour)
        if list_open_cutting_contours:
            list_faces = self.divide_face_with_open_cutting_contours(list_open_cutting_contours, inside)
        list_faces = self.divide_face_with_closed_cutting_contours(list_closed_cutting_contours, list_faces)
        list_faces = [face for face in list_faces if not math.isclose(face.area(), 0.0, abs_tol=1e-6)]
        return list_faces

    def divide_face_with_open_cutting_contours(self, list_open_cutting_contours, inside):
        """
        Divides a face 3D with a list of closed cutting contour, that is, it will cut holes on the face.

        :param list_open_cutting_contours: list containing the open cutting contours.
        :param inside: inside portion.
        :type inside: bool.
        :return: list divided faces.
        """
        list_faces = []
        if not self.surface2d.outer_contour.edge_polygon.is_trigo:
            self.surface2d.outer_contour.invert_inplace()
        new_faces_contours = self.surface2d.outer_contour.divide(list_open_cutting_contours, inside)
        new_inner_contours = len(new_faces_contours) * [[]]
        if self.surface2d.inner_contours:
            new_faces_contours, new_inner_contours = self.get_open_contour_divided_faces_inner_contours(
                new_faces_contours)
        if isinstance(self, Triangle3D):
            class_to_instanciate = PlaneFace3D
        else:
            class_to_instanciate = self.__class__
        for contour, inner_contours in zip(new_faces_contours, new_inner_contours):
            new_face = class_to_instanciate(self.surface3d, Surface2D(contour, inner_contours))
            list_faces.append(new_face)
        return list_faces

    def divide_face_with_closed_cutting_contours(self, list_closed_cutting_contours, list_faces):
        """
        Divides a Face3D with a list of Open cutting contours.

        Contours going from one side to another of the Face, or from the outer contour to one inner contour.

        :param list_closed_cutting_contours: list containing the closed cutting contours
        :param list_faces: list of already divided faces
        :return: list divided faces
        """
        for new_contour in list_closed_cutting_contours:
            if len(new_contour.primitives) >= 3 and \
                    new_contour.primitives[0].start == new_contour.primitives[-1].end:
                inner_contours1 = [new_contour]
                inner_contours2 = []
                if list_faces:
                    new_list_faces = self.get_closed_contour_divided_faces_inner_contours(list_faces, new_contour)
                    list_faces = list_faces + new_list_faces
                    continue
                for inner_contour in self.surface2d.inner_contours:
                    if new_contour.is_inside(inner_contour):
                        inner_contours2.append(inner_contour)
                        continue
                    inner_contours1.append(inner_contour)
                if isinstance(self, Triangle3D):
                    class_to_instanciate = PlaneFace3D
                else:
                    class_to_instanciate = self.__class__
                surf3d = self.surface3d
                surf2d = Surface2D(self.surface2d.outer_contour, inner_contours1)
                new_plane = class_to_instanciate(surf3d, surf2d)
                list_faces.append(new_plane)
                list_faces.append(class_to_instanciate(surf3d, Surface2D(new_contour, inner_contours2)))
                continue
            surf3d = self.surface3d
            surf2d = Surface2D(self.surface2d.outer_contour, [])
            if isinstance(self, Triangle3D):
                class_to_instanciate = PlaneFace3D
            else:
                class_to_instanciate = self.__class__
            new_plane = class_to_instanciate(surf3d, surf2d)
            list_faces.append(new_plane)
        return list_faces

    def get_open_contour_divided_faces_inner_contours(self, new_faces_contours):
        """
        If there is any inner contour, verifies which ones belong to the new divided faces.

        :param new_faces_contours: new faces outer contour.
        :return: valid_new_faces_contours, valid_new_faces_contours.
        """
        valid_new_faces_contours = []
        valid_inner_contours = []
        for new_face_contour in new_faces_contours:
            for inner_contour in self.surface2d.inner_contours:
                if new_face_contour.is_superposing(inner_contour):
                    break
            else:
                if new_face_contour not in valid_new_faces_contours:
                    inner_contours = []
                    for inner_contour in self.surface2d.inner_contours:
                        if new_face_contour.is_inside(inner_contour):
                            inner_contours.append(inner_contour)
                    valid_new_faces_contours.append(new_face_contour)
                    valid_inner_contours.append(inner_contours)
        return valid_new_faces_contours, valid_inner_contours

    def get_closed_contour_divided_faces_inner_contours(self, list_faces, new_contour):
        """
        If there is any inner contour, verifies which ones belong to the new divided faces.

        :param list_faces: list of new faces.
        :param new_contour: current new face outer contour.
        :return: a list of new faces with its inner contours.
        """
        new_list_faces = []
        for new_face in list_faces:
            if new_face.surface2d.outer_contour.is_inside(new_contour):
                inner_contours1 = []
                inner_contours2 = []
                if not new_face.surface2d.inner_contours:
                    new_face.surface2d.inner_contours = [new_contour]
                    break
                new_contour_not_sharing_primitives = True
                for i, inner_contour in enumerate(new_face.surface2d.inner_contours):
                    if new_contour.is_inside(inner_contour):
                        if any(inner_contour.primitive_over_contour(prim)
                               for prim in new_contour.primitives):
                            new_face.surface2d.inner_contours[i] = new_contour
                            break
                        inner_contours2.append(inner_contour)
                    elif not any(inner_contour.primitive_over_contour(prim) for prim in
                                 new_contour.primitives):
                        inner_contours1.append(inner_contour)
                    else:
                        new_contour_not_sharing_primitives = False
                else:
                    surf3d = new_face.surface3d
                    if inner_contours1:
                        if new_contour_not_sharing_primitives:
                            inner_contours1.append(new_contour)
                            new_face.surface2d.inner_contours = inner_contours1
                            break
                        surf2d = Surface2D(new_face.surface2d.outer_contour, inner_contours1)
                        new_plane = self.__class__(surf3d, surf2d)
                        new_list_faces.append(new_plane)
                    if inner_contours2:
                        new_list_faces.append(
                            self.__class__(surf3d, Surface2D(new_contour, inner_contours2)))
        return new_list_faces

    def select_face_intersecting_primitives(self, dict_intersecting_combinations):
        """
        Select face intersecting primitives from a dictionary containing all intersection combinations.

        :param dict_intersecting_combinations: dictionary containing all intersection combinations
        :return: list of intersecting primitives for current face
        """
        face_intersecting_primitives2d = []
        for intersecting_combination in dict_intersecting_combinations.keys():
            if self in (intersecting_combination[0], intersecting_combination[1]):
                for intersection_wire in dict_intersecting_combinations[intersecting_combination]:
                    if len(intersection_wire.primitives) != 1:
                        raise NotImplementedError
                    # primitive2 = intersection_wire
                    primitive2_2d = self.surface3d.contour3d_to_2d(intersection_wire)
                    if not self.surface2d.outer_contour.primitive_over_contour(primitive2_2d.primitives[0], tol=1e-7):
                        face_intersecting_primitives2d.append(primitive2_2d.primitives[0])
        return face_intersecting_primitives2d

    def get_inner_contours_cutting_primitives(self, list_cutting_contours, connectig_to_outer_contour):
        """
        Gets cutting primitives connected to face inner_contours.

        :param list_cutting_contours: list of contours for resulting from intersection with other faces.
        :param connectig_to_outer_contour: list of contours from list_cutting_contours connected to the outer contour
        and not to any outer contour.
        :return: lists for final face cutting primitives.
        """
        (inner_contours_connected_cutting_contour, dict_inner_contour_intersections,
         dict_cutting_contour_intersections, list_cutting_contours) = self.dictionnaries_cutting_contours(
            list_cutting_contours, connectig_to_outer_contour)
        valid_cutting_contours = []
        list_primitives1 = []
        list_primitives2 = []
        used_cutting_contours = []
        used_inner_contour = []
        for cutting_contour, inner_contours in inner_contours_connected_cutting_contour.items():
            primitives1 = []
            primitives2 = []
            if len(inner_contours) == 1:
                if all(dict_cutting_contour_intersections[inters] in connectig_to_outer_contour for inters in
                       dict_inner_contour_intersections[inner_contours[0]]) and cutting_contour not in \
                        used_cutting_contours and inner_contours[0] not in used_inner_contour:
                    inner_contour_splitting_points = dict_inner_contour_intersections[inner_contours[0]]
                    inner_contour_splitting_points = list(sorted(
                        inner_contour_splitting_points, key=lambda point, ic=inner_contours[0]: ic.abscissa(point)))

                    point1, point2 = self.inner_contour_cutting_points(inner_contour_splitting_points, cutting_contour)
                    primitives1.extend(inner_contours[0].extract_with_points(point1, point2, True))
                    primitives2.extend(inner_contours[0].extract_with_points(point1, point2, False))
                    if sum(prim.length() for prim in primitives2) > sum(prim.length() for prim in primitives1):
                        primitives1, primitives2 = primitives2, primitives1
                    primitives1.extend(dict_cutting_contour_intersections[point2].primitives +
                                       cutting_contour.primitives[:])
                    used_cutting_contours.extend([cutting_contour,
                                                  dict_cutting_contour_intersections[point2]])
                    used_inner_contour.append(inner_contours[0])
                    list_primitives1.append(primitives1)
                    list_primitives2.append(primitives2)
                elif cutting_contour not in valid_cutting_contours:
                    valid_cutting_contours.append(cutting_contour)
            elif len(inner_contours) == 2:
                inner_contour_splitting_points1 = dict_inner_contour_intersections[inner_contours[0]]
                inner_contour_splitting_points2 = dict_inner_contour_intersections[inner_contours[1]]
                inner_contour_splitting_points1 = list(sorted(
                    inner_contour_splitting_points1, key=lambda point, ic=inner_contours[0]: ic.abscissa(point)))
                inner_contour_splitting_points2 = list(sorted(
                    inner_contour_splitting_points2, key=lambda point, ic=inner_contours[1]: ic.abscissa(point)))
                inside1, inside2 = self.is_inside_portion(cutting_contour, inner_contour_splitting_points1,
                                                          inner_contour_splitting_points2)
                primitives1.extend(cutting_contour.primitives[:])
                contour_used = False
                for inner_contour, inner_contour_splitting_points, inside in zip(
                        inner_contours, [inner_contour_splitting_points1, inner_contour_splitting_points2],
                        [inside1, inside2]):
                    if inner_contour in used_inner_contour:
                        contour_used = True
                        continue
                    point1, point2 = self.inner_contour_cutting_points(inner_contour_splitting_points, cutting_contour)
                    primitives1.extend(inner_contour.extract_with_points(point1, point2, inside))
                    primitives1.extend(dict_cutting_contour_intersections[point2].primitives)
                    primitives2.extend(inner_contour.extract_with_points(point1, point2, not inside))
                    used_cutting_contours.extend([cutting_contour, dict_cutting_contour_intersections[point2]])
                if contour_used:
                    list_primitives1 = self.get_connecting_contour(list_primitives1, primitives1)
                else:
                    list_primitives1.append(primitives1)
                list_primitives2.append(primitives2)
                used_inner_contour.extend(inner_contours)
            else:
                raise NotImplementedError
        valid_cutting_contours = [contour for contour in valid_cutting_contours
                                  if contour not in used_cutting_contours]
        new_cutting_contours = [volmdlr.wires.Contour2D(list_prim).order_contour()
                                for list_prim in list_primitives1]
        for list_prim in list_primitives2:
            new_cutting_contours.extend(volmdlr.wires.Contour2D.contours_from_edges(list_prim))
        return new_cutting_contours, valid_cutting_contours

    def dictionnaries_cutting_contours(self, list_cutting_contours, connectig_to_outer_contour):
        inner_contours_connected_cutting_contour = {}
        dict_inner_contour_intersections = {}
        dict_cutting_contour_intersections = {}
        for inner_contour in self.surface2d.inner_contours:
            if not inner_contour.edge_polygon.is_trigo:
                inner_contour.invert_inplace()
            dict_inner_contour_intersections[inner_contour] = []
            for cutting_contour in list_cutting_contours:
                inner_contour_intersections = inner_contour.contour_intersections(cutting_contour)
                if inner_contour_intersections:
                    dict_inner_contour_intersections[inner_contour].extend(inner_contour_intersections)
                    if cutting_contour not in inner_contours_connected_cutting_contour:
                        inner_contours_connected_cutting_contour[cutting_contour] = [inner_contour]
                    else:
                        inner_contours_connected_cutting_contour[cutting_contour].append(inner_contour)
                for intersection in inner_contour_intersections:
                    dict_cutting_contour_intersections[intersection] = cutting_contour
            splitting_points = dict_inner_contour_intersections[inner_contour]
            splitting_points = list(sorted(
                splitting_points, key=lambda point, ic=inner_contour: ic.abscissa(point)))
            remove_splitting_points, new_inner_contour, remove_cutting_contour = self.inner_contours_recalculation(
                inner_contour, splitting_points, dict_cutting_contour_intersections, connectig_to_outer_contour)
            (dict_cutting_contour_intersections, dict_inner_contour_intersections,
             inner_contours_connected_cutting_contour, list_cutting_contours) = \
                self.updated_dictionnaries_cutting_contours(remove_splitting_points, remove_cutting_contour,
                                                            splitting_points, dict_cutting_contour_intersections,
                                                            inner_contour, new_inner_contour, list_cutting_contours,
                                                            dict_inner_contour_intersections,
                                                            inner_contours_connected_cutting_contour)
            inner_contour = new_inner_contour
        return (inner_contours_connected_cutting_contour, dict_inner_contour_intersections,
                dict_cutting_contour_intersections, list_cutting_contours)

    @staticmethod
    def inner_contour_cutting_points(inner_contour_splitting_points, cutting_contour):
        """
        Searches the inner contour points where it must be cut.

        :param inner_contour_splitting_points: all points os intersection with this inner contour.
        :param cutting_contour: first cutting contour being used to cut inner contour.
        :return: point1, point2
        """
        if volmdlr.core.point_in_list(cutting_contour.primitives[0].start, inner_contour_splitting_points):
            index_point1 = volmdlr.core.get_point_index_in_list(cutting_contour.primitives[0].start,
                                                                inner_contour_splitting_points)
        else:
            index_point1 = volmdlr.core.get_point_index_in_list(cutting_contour.primitives[-1].end,
                                                                inner_contour_splitting_points)
        if index_point1 != len(inner_contour_splitting_points) - 1:
            index_point2 = index_point1 + 1
        else:
            index_point2 = 0
        point1 = inner_contour_splitting_points[index_point1]
        point2 = inner_contour_splitting_points[index_point2]
        return point1, point2

    @staticmethod
    def is_inside_portion(cutting_contour, inner_contour_splitting_points1, inner_contour_splitting_points2):
        """
        For multiple inner contour intersections with cutting contours, defines if we get the inside or outside portion.

        :param cutting_contour: cutting_contour cutting the two inner contours.
        :param inner_contour_splitting_points1: splitting points for contour1.
        :param inner_contour_splitting_points2: splitting points for contour1.
        :return:
        """
        if (not cutting_contour.primitives[0].start.is_close(inner_contour_splitting_points1[-1]) and
            not cutting_contour.primitives[-1].end.is_close(inner_contour_splitting_points2[-1])) or \
                (not cutting_contour.primitives[0].start.is_close(inner_contour_splitting_points2[-1]) and
                 not cutting_contour.primitives[-1].end.is_close(inner_contour_splitting_points2[-1])):
            is_inside1 = True
            is_inside2 = False
        elif (cutting_contour.primitives[0].start.is_close(inner_contour_splitting_points1[-1]) and
              cutting_contour.primitives[-1].end.is_close(inner_contour_splitting_points2[-1])):
            is_inside1 = True
            is_inside2 = False
        elif (not cutting_contour.primitives[0].start.is_close(inner_contour_splitting_points1[-1]) and
              cutting_contour.primitives[-1].end.is_close(inner_contour_splitting_points2[-1])) or \
                (cutting_contour.primitives[0].start.is_close(inner_contour_splitting_points1[-1]) and
                 not cutting_contour.primitives[-1].end.is_close(inner_contour_splitting_points2[-1])):
            is_inside1 = True
            is_inside2 = True
        else:
            raise NotImplementedError
        return is_inside1, is_inside2

    @staticmethod
    def get_connecting_contour(lists_primitives, inner_primitives):
        """
        Find which contour from resulting inner contour splitting is connected to saved cutting_contours.

        :param lists_primitives: saved cutting contours.
        :param inner_primitives: splited inner contour.
        :return: updated saved cutting contours.
        """
        if not lists_primitives:
            lists_primitives.extend(inner_primitives)
            return lists_primitives
        new_list_primitives = lists_primitives[:]
        for i, list_prim in enumerate(lists_primitives):
            if any(prim in list_prim for prim in inner_primitives):
                new_primitives = list_prim + [prim for prim in inner_primitives if prim not in list_prim]
                new_list_primitives[i] = new_primitives
                break
        lists_primitives = new_list_primitives[:]
        return lists_primitives

    def inner_contours_recalculation(self, inner_contour, splitting_points, splitting_points_and_cutting_contour,
                                     connectig_to_outer_contour):
        """
        Recalculates inner contours if a cutting contour is connected to an inner contour at two ends.

        Verifies if there is a cutting contours from face intersections connected to an inner contour at the two ends,
        if true this inner contour is updated with this cutting contour.

        :param inner_contour: inner contour.
        :param splitting_points: current inner contour splitting points.
        :param splitting_points_and_cutting_contour: dictionary containing all splitting points and
        the corresponding cutting contour.
        :param connectig_to_outer_contour: list of the cutting contours connected to the outer contour.
        :return: splitting points to be removed from list of splitting points and current inner contour updated.
        """
        j = self.surface2d.inner_contours.index(inner_contour)
        remove_splitting_points = []
        remove_cutting_contour = []
        for point1, point2 in zip(splitting_points[:-1], splitting_points[1:]):
            if splitting_points_and_cutting_contour[point1] not in connectig_to_outer_contour and \
                    splitting_points_and_cutting_contour[point2] not in connectig_to_outer_contour and \
                    splitting_points_and_cutting_contour[point1] == splitting_points_and_cutting_contour[point2]:
                remove_cutting_contour.append(splitting_points_and_cutting_contour[point1])
                remove_splitting_points.extend([point1, point2])
                primitives1 = inner_contour.extract_with_points(point1, point2, True) + \
                              splitting_points_and_cutting_contour[point1].primitives
                primitives2 = inner_contour.extract_with_points(point1, point2, False) + \
                              splitting_points_and_cutting_contour[point1].primitives
                contour1 = volmdlr.wires.Contour2D(primitives1).order_contour()
                contour2 = volmdlr.wires.Contour2D(primitives2).order_contour()
                if contour1.is_inside(inner_contour):
                    self.surface2d.inner_contours[j] = contour1
                    inner_contour = self.surface2d.inner_contours[j]
                    remove_splitting_points.extend([point1, point2])
                elif contour2.is_inside(inner_contour):
                    self.surface2d.inner_contours[j] = contour2
                    inner_contour = self.surface2d.inner_contours[j]
                    remove_splitting_points.extend([point1, point2])
        return remove_splitting_points, inner_contour, remove_cutting_contour

    @staticmethod
    def updated_dictionnaries_cutting_contours(remove_splitting_points, remove_cutting_contour, splitting_points,
                                               dict_cutting_contour_intersections, old_inner_contour,
                                               new_inner_contour, list_cutting_contours,
                                               dict_inner_contour_intersections,
                                               inner_contours_connected_cutting_contour):
        for remove_point in remove_splitting_points:
            if remove_point in splitting_points:
                splitting_points.remove(remove_point)
            if remove_point in dict_cutting_contour_intersections:
                del dict_cutting_contour_intersections[remove_point]
        del dict_inner_contour_intersections[old_inner_contour]
        dict_inner_contour_intersections[new_inner_contour] = splitting_points
        for contour in remove_cutting_contour:
            if contour in list_cutting_contours:
                list_cutting_contours.remove(contour)
            if contour in inner_contours_connected_cutting_contour:
                del inner_contours_connected_cutting_contour[contour]
        for cutting_contour, innr_cntrs in inner_contours_connected_cutting_contour.items():
            if old_inner_contour in innr_cntrs:
                inner_contours_connected_cutting_contour[cutting_contour].remove(old_inner_contour)
                inner_contours_connected_cutting_contour[cutting_contour].append(new_inner_contour)
        return (dict_cutting_contour_intersections, dict_inner_contour_intersections,
                inner_contours_connected_cutting_contour, list_cutting_contours)


class PlaneFace3D(Face3D):
    """
    Defines a PlaneFace3D class.

    :param surface3d: a plane 3d.
    :type surface3d: Plane3D.
    :param surface2d: a 2d surface to define the plane face.
    :type surface2d: Surface2D.
    """
    _standalone_in_db = False
    _generic_eq = True
    _non_serializable_attributes = ['bounding_box', 'polygon2D']
    _non_data_eq_attributes = ['name', 'bounding_box', 'outer_contour3d',
                               'inner_contours3d']
    _non_data_hash_attributes = []

    def __init__(self, surface3d: Plane3D, surface2d: Surface2D,
                 name: str = ''):
        self._bbox = None
        Face3D.__init__(self,
                        surface3d=surface3d,
                        surface2d=surface2d,
                        name=name)

    def copy(self, deep=True, memo=None):
        return PlaneFace3D(self.surface3d.copy(deep, memo), self.surface2d.copy(),
                           self.name)

    @property
    def bounding_box(self):
        """
        Returns the boundary box of a PlanFace3D.

        """
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def get_bounding_box(self):
        return self.outer_contour3d.bounding_box

    def distance_to_point(self, point, return_other_point=False):
        """
        Calculates the distance from a plane face and a point.

        :param point: point to verify.
        :param return_other_point: bool to decide if corresponding point on face should be returned.
        :return: distance to planeface3D.
        """

        projected_pt = point.plane_projection3d(self.surface3d.frame.origin,
                                                self.surface3d.frame.u,
                                                self.surface3d.frame.v)
        projection_distance = point.point_distance(projected_pt)

        if self.point_belongs(projected_pt):
            if return_other_point:
                return projection_distance, projected_pt
            return projection_distance

        point_2d = point.to_2d(self.surface3d.frame.origin, self.surface3d.frame.u,
                               self.surface3d.frame.v)

        polygon2d = self.surface2d.outer_contour.to_polygon(angle_resolution=10)
        border_distance, other_point = polygon2d.point_border_distance(point_2d, return_other_point=True)

        other_point = self.surface3d.point2d_to_3d(volmdlr.Point2D(*other_point))

        if return_other_point:
            return (projection_distance ** 2 + border_distance ** 2) ** 0.5, \
                other_point
        return (projection_distance ** 2 + border_distance ** 2) ** 0.5

    def minimum_distance_points_plane(self, other_plane_face, return_points=False):
        """
        Given two plane faces, calculates the points which corresponds to the minimal distance between these two faces.

        :param other_plane_face: Second plane face.
        :param return_points: Boolean to return corresponding points or not.
        :return: minimal distance.
        """

        min_distance = math.inf
        for edge1 in self.outer_contour3d.primitives:
            for edge2 in other_plane_face.outer_contour3d.primitives:
                dist = edge1.minimum_distance(edge2,
                                              return_points=return_points)
                if return_points:
                    if dist[0] < min_distance:
                        min_distance = dist[0]
                        point1, point2 = dist[1], dist[2]
                else:
                    if dist < min_distance:
                        min_distance = dist
        if return_points:
            return min_distance, point1, point2
        return min_distance

    def linesegment_inside(self, linesegment: vme.LineSegment3D):
        direction_vector = linesegment.unit_direction_vector()
        if not math.isclose(abs(direction_vector.dot(self.surface3d.frame.w)), 0.0, abs_tol=1e-6):
            return False
        for point in [linesegment.start, linesegment.middle_point(), linesegment.end]:
            if not self.point_belongs(point):
                return False
        return True

    def circle_inside(self, circle: volmdlr.wires.Circle3D):
        if not math.isclose(abs(circle.frame.w.dot(self.surface3d.frame.w)), 1.0, abs_tol=1e-6):
            return False
        points = circle.discretization_points(number_points=4)
        for point in points:
            if not self.point_belongs(point):
                return False
        return True

    def planeface_intersections(self, planeface):
        face2_plane_interections = planeface.surface3d.plane_intersection(self.surface3d)
        if not face2_plane_interections:
            return []
        points_intersections = []
        for contour in [self.outer_contour3d, planeface.outer_contour3d] + self.inner_contours3d + \
                       planeface.inner_contours3d:
            for intersection in contour.line_intersections(face2_plane_interections[0]):
                if intersection and not volmdlr.core.point_in_list(intersection, points_intersections):
                    points_intersections.append(intersection)
        points_intersections = face2_plane_interections[0].sort_points_along_line(points_intersections)
        planeface_intersections = []
        for point1, point2 in zip(points_intersections[:-1], points_intersections[1:]):
            linesegment3d = vme.LineSegment3D(point1, point2)
            over_self_outer_contour = self.outer_contour3d.primitive_over_contour(linesegment3d)
            over_planeface_outer_contour = planeface.outer_contour3d.primitive_over_contour(linesegment3d)
            if over_self_outer_contour and over_planeface_outer_contour:
                continue
            if self.edge3d_inside(linesegment3d) or over_self_outer_contour:
                if planeface.edge3d_inside(linesegment3d):
                    planeface_intersections.append(volmdlr.wires.Wire3D([linesegment3d]))
                elif over_planeface_outer_contour:
                    planeface_intersections.append(volmdlr.wires.Wire3D([linesegment3d]))
        return planeface_intersections

    def triangle_intersections(self, triangleface):
        return self.planeface_intersections(triangleface)

    def cylindricalface_intersections(self, cylindricalface: 'CylindricalFace3D'):
        cylindricalsurfaceface_intersections = cylindricalface.surface3d.plane_intersection(self.surface3d)
        if not isinstance(cylindricalsurfaceface_intersections[0], vme.Line3D):
            if all(self.edge3d_inside(intersection) and cylindricalface.edge3d_inside(intersection)
                   for intersection in cylindricalsurfaceface_intersections):
                return cylindricalsurfaceface_intersections
        intersections_points = self.face_intersections_outer_contour(cylindricalface)
        for point in cylindricalface.face_intersections_outer_contour(self):
            if point not in intersections_points:
                intersections_points.append(point)
        face_intersections = []
        for primitive in cylindricalsurfaceface_intersections:
            points_on_primitive = []
            for point in intersections_points:
                if primitive.point_belongs(point):
                    points_on_primitive.append(point)
            if not points_on_primitive:
                continue
            if isinstance(primitive, vme.Line3D):
                points_on_primitive = primitive.sort_points_along_line(points_on_primitive)
            else:
                points_on_primitive = primitive.sort_points_along_wire(points_on_primitive)
                points_on_primitive = points_on_primitive + [points_on_primitive[0]]
            for point1, point2 in zip(points_on_primitive[:-1], points_on_primitive[1:]):
                edge = primitive.trim(point1, point2)
                if self.edge3d_inside(edge) and cylindricalface.edge3d_inside(edge):
                    face_intersections.append(volmdlr.wires.Wire3D([edge]))
        return face_intersections

    def minimum_distance(self, other_face, return_points=False):
        """
        Returns the minimum distance between the current face and the specified other face.

        :param other_face: Face to evaluate the minimum distance.
        :type other_face: :class:`PlaneFace3D`
        :param return_points: A boolean value indicating whether to return the minimum distance as
            well as the two points on each face that are closest to each other. If True, the function
            returns a tuple of the form (distance, point1, point2). If False, the function only returns
            the distance.
        :type return_points: bool
        :return: If the return_points parameter is set to True, it also returns the two points on each face
            that are closest to each other. The return type is a float if return_points is False and
            a tuple (distance, point1, point2) if return_points is True.
        :rtype: float, Tuple[float, float, float]
        """
        if other_face.__class__ is CylindricalFace3D:
            point1, point2 = other_face.minimum_distance_points_cyl(self)
            if return_points:
                return point1.point_distance(point2), point1, point2
            return point1.point_distance(point2)

        if other_face.__class__ is PlaneFace3D:
            if return_points:
                dist, point1, point2 = self.minimum_distance_points_plane(other_face,
                                                                          return_points=return_points)
                return dist, point1, point2
            dist = self.minimum_distance_points_plane(other_face,
                                                      return_points=return_points)
            return dist

        if other_face.__class__ is ToroidalFace3D:
            point1, point2 = other_face.minimum_distance_points_plane(self)
            if return_points:
                return point1.point_distance(point2), point1, point2
            return point1.point_distance(point2)

        raise NotImplementedError

    def is_adjacent(self, face2: Face3D):
        contour1 = self.outer_contour3d.to_2d(
            self.surface3d.frame.origin,
            self.surface3d.frame.u,
            self.surface3d.frame.v)
        contour2 = face2.outer_contour3d.to_2d(
            self.surface3d.frame.origin,
            self.surface3d.frame.u,
            self.surface3d.frame.v)
        if contour1.is_sharing_primitives_with(contour2, False):
            return True
        return False

    @staticmethod
    def merge_faces(list_coincident_faces: List[Face3D]):
        valid_coicident_faces = list_coincident_faces[:]
        list_new_faces = []
        list_inner_contours = []
        merge_finished = False
        face0 = valid_coicident_faces[0]
        merged_contour = face0.outer_contour3d.to_2d(face0.surface3d.frame.origin,
                                                     face0.surface3d.frame.u,
                                                     face0.surface3d.frame.v)
        valid_coicident_faces.remove(face0)
        while not merge_finished:
            adjacent_faces = False
            list_inner_contours = []
            for face in valid_coicident_faces:
                adjacent_faces = False
                face_inside = False
                contour = face.outer_contour3d.to_2d(face0.surface3d.frame.origin,
                                                     face0.surface3d.frame.u,
                                                     face0.surface3d.frame.v)
                if contour.is_sharing_primitives_with(merged_contour):
                    merged_contour_results = merged_contour.union(contour)
                    merged_contour = merged_contour_results[0]
                    merged_inner_contours = merged_contour_results[1:]
                    list_inner_contours.extend(merged_inner_contours)
                    list_inner_contours.extend(face.surface2d.inner_contours)
                    valid_coicident_faces.remove(face)
                    adjacent_faces = True
                    break
                if merged_contour.is_inside(contour):
                    valid_coicident_faces.remove(face)
                    face_inside = True
                    break
            if not adjacent_faces and not face_inside and valid_coicident_faces:
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

    def cut_by_coincident_face(self, face):
        """
        Cuts face1 with another coincident face2.

        :param face: a face 3d.
        :type face: Face3D.
        :return: a list of faces 3d.
        :rtype: List[Face3D].
        """

        if not self.surface3d.is_coincident(face.surface3d):
            raise ValueError('The faces are not coincident')

        if self.face_inside(face):
            return self.divide_face([face.surface2d.outer_contour], True)
        # if face.is_inside(self):
        #     return face.divide_face([self.surface2d.outer_contour], True)

        outer_contour_1 = self.surface2d.outer_contour
        outer_contour_2 = self.surface3d.contour3d_to_2d(face.outer_contour3d)

        if (face.face_inside(self)
                and not outer_contour_1.contour_intersections(outer_contour_2)):
            return self.divide_face(face.surface2d.inner_contours, True)

        inner_contours = self.surface2d.inner_contours
        inner_contours.extend([self.surface3d.contour3d_to_2d(
            contour) for contour in face.inner_contours3d])

        contours = outer_contour_1.cut_by_wire(outer_contour_2)

        surfaces = []
        for contour in contours:
            inners = []
            for inner_c in inner_contours:
                if contour.is_inside(inner_c):
                    inners.append(inner_c)
            surfaces.append(Surface2D(contour, inners))

        return [self.__class__(self.surface3d, surface2d) for surface2d in surfaces]

    def check_inner_contours(self, face):
        c_inners_1 = self.surface2d.inner_contours
        c_inners_2 = [self.surface3d.contour3d_to_2d(inner) for inner in face.inner_contours3d]
        inside = set()
        for inner_contour1 in c_inners_1:
            for inner_contour2 in c_inners_2:
                if inner_contour1.is_superposing(inner_contour2):
                    inside.add(False)
                else:
                    inside.add(inner_contour2.is_inside(inner_contour1))
        return inside

    @staticmethod
    def update_faces_with_divided_faces(divided_faces, face2_2, used, list_faces):
        for d_face in divided_faces:

            if d_face.outer_contour3d.is_superposing(face2_2.outer_contour3d):
                if face2_2.surface2d.inner_contours:
                    divided_faces_d_face = []
                    for inner in face2_2.surface2d.inner_contours:

                        if True in [(((abs(inner_d.area() - inner.area()) < 1e-6)
                                      and inner.center_of_mass().is_close(inner_d.center_of_mass()))
                                     or inner_d.is_inside(inner))
                                    for inner_d in d_face.surface2d.inner_contours]:
                            divided_faces_d_face = ['', d_face]
                            continue

                        divided_faces_d_face = d_face.divide_face([inner], True)
                        divided_faces_d_face.sort(key=lambda x: x.area())

                        list_faces.append(divided_faces_d_face[0])
                        d_face = divided_faces_d_face[1]

                    if divided_faces_d_face:
                        list_faces.append(divided_faces_d_face[1])

                else:
                    list_faces.append(d_face)
            else:
                used.append(d_face)

        return used, list_faces

    def project_faces(self, faces):
        """
        Divide self based on the faces outer, and inner contours.

        :param faces: DESCRIPTION
        :type faces: TYPE
        :return: DESCRIPTION
        :rtype: TYPE
        """

        used_faces, list_faces = {}, []

        for _, face2 in enumerate(faces):
            contour1 = self.surface2d.outer_contour
            contour2 = self.surface3d.contour3d_to_2d(face2.outer_contour3d)

            inside = self.check_inner_contours(face2)
            if (self.surface3d.is_coincident(face2.surface3d)
                    and (contour1.is_overlapping(contour2)
                         or (contour1.is_inside(contour2) or True in inside))):

                if self in used_faces:
                    faces_1, face2_2 = used_faces[self][:], face2
                else:
                    faces_1, face2_2 = [self], face2

                used = []
                for face1_1 in faces_1:
                    plane3d = face1_1.surface3d
                    s2d = Surface2D(outer_contour=plane3d.contour3d_to_2d(face2_2.outer_contour3d),
                                    inner_contours=[
                                        plane3d.contour3d_to_2d(contour) for contour in face2_2.inner_contours3d])
                    face2_2 = PlaneFace3D(surface3d=plane3d, surface2d=s2d)

                    divided_faces = face1_1.cut_by_coincident_face(face2_2)

                    used, list_faces = self.update_faces_with_divided_faces(
                        divided_faces, face2_2, used, list_faces)
                used_faces[self] = used

        try:
            if isinstance(used_faces[self], list):
                list_faces.extend(used_faces[self])
            else:
                list_faces.append(used_faces[self])
        except KeyError:
            list_faces.append(self)

        return list_faces

    def get_geo_lines(self, tag: int, line_loop_tag: List[int]):
        """
        Gets the lines that define a PlaneFace3D in a .geo file.
        """

        return 'Plane Surface(' + str(tag) + ') = {' + str(line_loop_tag)[1:-1] + '};'


class Triangle3D(PlaneFace3D):
    """
    Defines a Triangle3D class.

    :param point1: The first point.
    :type point1: volmdlr.Point3D.
    :param point2: The second point.
    :type point2: volmdlr.Point3D.
    :param point3: The third point.
    :type point3: volmdlr.Point3D.
    """
    _standalone_in_db = False

    def __init__(self, point1: volmdlr.Point3D, point2: volmdlr.Point3D,
                 point3: volmdlr.Point3D, alpha=1, color=None, name: str = ''):
        self.point1 = point1
        self.point2 = point2
        self.point3 = point3
        self.points = [self.point1, self.point2, self.point3]
        self.color = color
        self.alpha = alpha
        self.name = name

        self._surface3d = None
        self._surface2d = None
        self._bbox = None
        self._outer_contour3d = None
        self._inner_contours3d = None
        # self.bounding_box = self._bounding_box()

        # DessiaObject.__init__(self, name=name)

    def _data_hash(self):
        """
        Using point approx hash to speed up.

        """
        return self.point1.approx_hash() + self.point2.approx_hash() + self.point3.approx_hash()

    def _data_eq(self, other_object):
        if other_object.__class__.__name__ != self.__class__.__name__:
            return False
        self_set = {self.point1, self.point2, self.point3}
        other_set = {other_object.point1, other_object.point2, other_object.point3}
        if self_set != other_set:
            return False
        return True

    @property
    def bounding_box(self):
        """
        Returns the surface bounding box.
        """
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bouding_box):
        self._bbox = new_bouding_box

    def get_bounding_box(self):
        return volmdlr.core.BoundingBox.from_points([self.point1,
                                                     self.point2,
                                                     self.point3])

    @property
    def surface3d(self):
        if self._surface3d is None:
            self._surface3d = Plane3D.from_3_points(self.point1, self.point2, self.point3)
        return self._surface3d

    @property
    def surface2d(self):
        if self._surface2d is None:
            plane3d = self.surface3d
            contour3d = volmdlr.wires.Contour3D([vme.LineSegment3D(self.point1, self.point2),
                                                 vme.LineSegment3D(self.point2, self.point3),
                                                 vme.LineSegment3D(self.point3, self.point1)])

            contour2d = contour3d.to_2d(plane3d.frame.origin,
                                        plane3d.frame.u, plane3d.frame.v)

            self._surface2d = Surface2D(outer_contour=contour2d, inner_contours=[])

        return self._surface2d

    def to_dict(self, *args, **kwargs):
        dict_ = {'object_class': 'volmdlr.faces.Triangle3D',
                 'point1': self.point1.to_dict(),
                 'point2': self.point2.to_dict(),
                 'point3': self.point3.to_dict()}
        if self.name:
            dict_['name'] = self.name
        return dict_

    @classmethod
    def dict_to_object(cls, dict_, *args, **kwargs):
        point1 = volmdlr.Point3D.dict_to_object(dict_['point1'])
        point2 = volmdlr.Point3D.dict_to_object(dict_['point2'])
        point3 = volmdlr.Point3D.dict_to_object(dict_['point3'])
        return cls(point1, point2, point3, dict_.get('name', ""))

    def area(self) -> float:
        """
        Calculates the area for the Triangle3D.

        :return: area triangle.
        :rtype: float.

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
        # Basis = vector point1 to point 2d
        return 2 * self.area() / self.point1.point_distance(self.point2)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Triangle3D.

        :param side: 'old' or 'new'.
        """
        np1 = self.point1.frame_mapping(frame, side)
        np2 = self.point2.frame_mapping(frame, side)
        np3 = self.point3.frame_mapping(frame, side)
        return self.__class__(np1, np2, np3, self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated in-place.

        :param side: 'old' or 'new'
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.point1.frame_mapping_inplace(frame, side)
        self.point2.frame_mapping_inplace(frame, side)
        self.point3.frame_mapping_inplace(frame, side)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box

    def copy(self, deep=True, memo=None):
        return Triangle3D(self.point1.copy(), self.point2.copy(), self.point3.copy(),
                          self.name)

    def triangulation(self):
        return vmd.DisplayMesh3D([vmd.Node3D.from_point(self.point1),
                                  vmd.Node3D.from_point(self.point2),
                                  vmd.Node3D.from_point(self.point3)],
                                 [(0, 1, 2)])

    def translation(self, offset: volmdlr.Vector3D):
        """
        Plane3D translation.

        :param offset: translation vector.
        :return: A new translated Plane3D.
        """
        new_point1 = self.point1.translation(offset)
        new_point2 = self.point2.translation(offset)
        new_point3 = self.point3.translation(offset)

        new_triangle = Triangle3D(new_point1, new_point2, new_point3,
                                  self.alpha, self.color, self.name)
        return new_triangle

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Plane3D translation. Object is updated in-place.

        :param offset: translation vector.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.point1.translation_inplace(offset)
        self.point2.translation_inplace(offset)
        self.point3.translation_inplace(offset)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Triangle3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated Triangle3D.
        """
        new_point1 = self.point1.rotation(center, axis, angle)
        new_point2 = self.point2.rotation(center, axis, angle)
        new_point3 = self.point3.rotation(center, axis, angle)
        new_triangle = Triangle3D(new_point1, new_point2, new_point3,
                                  self.alpha, self.color, self.name)
        return new_triangle

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        Triangle3D rotation. Object is updated inplace.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: rotation angle.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        self.point1.rotation_inplace(center, axis, angle)
        self.point2.rotation_inplace(center, axis, angle)
        self.point3.rotation_inplace(center, axis, angle)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box

    def subdescription(self, resolution=0.01):
        """
        Returns a list of Point3D with resolution as max between Point3D.
        """

        lengths = [self.points[0].point_distance(self.points[1]),
                   self.points[1].point_distance(self.points[2]),
                   self.points[2].point_distance(self.points[0])]
        max_length = max(lengths)

        if max_length <= resolution:
            return self.points

        pos_length_max = lengths.index(max_length)
        new_points = [self.points[-3 + pos_length_max + k] for k in range(3)]

        vector = new_points[0] - new_points[1]
        vector.normalize()
        points_0_1 = []

        for k in range(int(max_length / resolution) + 2):
            if k == 0:
                points_0_1.append(new_points[1])
            points_0_1.append(new_points[1] + vector * min(k * resolution, max_length))

        vector, length_2_1 = new_points[2] - new_points[1], new_points[2].point_distance(new_points[1])
        vector.normalize()
        points_in = []

        for p0_1 in points_0_1:
            point_on_2_1 = new_points[1] + vector * min(points_0_1[0].point_distance(p0_1) * length_2_1 / max_length,
                                                        length_2_1)

            length_2_0 = point_on_2_1.point_distance(p0_1)
            nb_int = int(length_2_0 / resolution) + 2
            if nb_int == 2:
                points_in.append(point_on_2_1)
            else:
                vector_2_0 = point_on_2_1 - p0_1
                vector_2_0.normalize()
                step_in = length_2_0 / (nb_int - 1)
                for i in range(nb_int):
                    if min(i * step_in, length_2_0) != 0:
                        points_in.append(p0_1 + vector_2_0 * min(i * step_in, length_2_0))

        return npy.unique(points_0_1 + points_in).tolist()

    def subdescription_to_triangles(self, resolution=0.01):
        """
        Returns a list of Triangle3D with resolution as max length of subtriangles side.

        """

        sub_triangles = [self.points]

        while True:
            triangles = []
            for subtri in sub_triangles:
                lengths = [subtri[0].point_distance(subtri[1]),
                           subtri[1].point_distance(subtri[2]),
                           subtri[2].point_distance(subtri[0])]
                max_length = max(lengths)

                if max_length > resolution:
                    pos_length_max = lengths.index(max_length)
                    pt_mid = (subtri[-3 + pos_length_max] + subtri[-3 + pos_length_max + 1]) / 2
                    triangles.extend([[subtri[-3 + pos_length_max], pt_mid, subtri[-3 + pos_length_max + 2]],
                                      [subtri[-3 + pos_length_max + 1], pt_mid, subtri[-3 + pos_length_max + 2]]])

                else:
                    triangles.append(subtri)

            if len(sub_triangles) == len(triangles):
                break

            sub_triangles = triangles

        return [Triangle3D(subtri[0], subtri[1], subtri[2]) for subtri in sub_triangles]

    def middle(self):
        return (self.point1 + self.point2 + self.point3) / 3

    def normal(self):
        """
        Get the normal vector to the face.

        Returns
        -------
        normal to the face

        """
        normal = self.surface3d.frame.w
        normal.normalize()
        return normal


class CylindricalFace3D(Face3D):
    """
    Defines a CylindricalFace3D class.

    :param surface3d: a cylindrical surface 3d.
    :type surface3d: CylindricalSurface3D.
    :param surface2d: a 2d surface to define the cylindrical face.
    :type surface2d: Surface2D.

    :Example:

        contours 2d is rectangular and will create a classic cylinder with x= 2*pi*radius, y=h
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
        self._bbox = None

    def copy(self, deep=True, memo=None):
        return CylindricalFace3D(self.surface3d.copy(deep, memo), self.surface2d.copy(),
                                 self.name)

    @property
    def bounding_box(self):
        """
        Returns the surface bounding box.
        """
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bouding_box):
        self._bbox = new_bouding_box

    def get_bounding_box(self):
        """
        Computes the bounding box using the contour 3d. true in this case of cylindrical face (not general).
        """
        return self.outer_contour3d.bounding_box

    def triangulation_lines(self, angle_resolution=5):
        theta_min, theta_max, zmin, zmax = self.surface2d.bounding_rectangle.bounds()
        delta_theta = theta_max - theta_min
        nlines = math.ceil(delta_theta * angle_resolution)
        lines = []
        for i in range(nlines):
            theta = theta_min + (i + 1) / (nlines + 1) * delta_theta
            lines.append(vme.Line2D(volmdlr.Point2D(theta, zmin),
                                    volmdlr.Point2D(theta, zmax)))
        return lines, []

    def point_belongs(self, point3d: volmdlr.Point3D):
        """
        Tells you if a point is on the 3D Cylindrical face and inside its contour.
        """
        point2d = self.surface3d.point3d_to_2d(point3d)
        point2d_plus_2pi = point2d.translation(volmdlr.Point2D(volmdlr.TWO_PI, 0))
        point2d_minus_2pi = point2d.translation(volmdlr.Point2D(-volmdlr.TWO_PI, 0))
        check_point3d = self.surface3d.point2d_to_3d(point2d)
        if check_point3d.point_distance(point3d) > 1e-6:
            return False

        return any(self.surface2d.point_belongs(pt2d) for pt2d in [point2d, point2d_plus_2pi, point2d_minus_2pi])

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        angle_resolution = 11
        z_resolution = 5
        theta_min, theta_max, zmin, zmax = self.surface2d.bounding_rectangle().bounds()
        delta_theta = theta_max - theta_min
        number_points_x = int(delta_theta * angle_resolution)

        delta_z = zmax - zmin
        number_points_y = int(delta_z * z_resolution)

        return number_points_x, number_points_y

    def minimum_distance(self, other_face, return_points=False):
        if other_face.__class__ is CylindricalFace3D:
            point1, point2 = self.minimum_distance_points_cyl(other_face)
            if return_points:
                return point1.point_distance(point2), point1, point2
            return point1.point_distance(point2)

        if other_face.__class__ is PlaneFace3D:
            point1, point2 = self.minimum_distance_points_plane(other_face)
            if return_points:
                return point1.point_distance(point2), point1, point2
            return point1.point_distance(point2)

        if other_face.__class__ is ToroidalFace3D:
            point1, point2 = other_face.minimum_distance_points_cyl(self)
            if return_points:
                return point1.point_distance(point2), point1, point2
            return point1.point_distance(point2)

        raise NotImplementedError

    def adjacent_direction(self, other_face3d):
        """
        Find out in which direction the faces are adjacent.

        :param other_face3d: The face to evaluation.
        :type other_face3d: volmdlr.faces.CylindricalFace3D
        """

        contour1 = self.outer_contour3d
        contour2 = other_face3d.outer_contour3d
        point1, point2 = contour1.shared_primitives_extremities(contour2)

        coord = point1 - point2
        coord = [abs(coord.x), abs(coord.y)]

        if coord.index(max(coord)) == 0:
            return 'x'
        return 'y'

    def get_geo_lines(self, tag: int, line_loop_tag: List[int]):
        """
        Gets the lines that define a CylindricalFace3D in a .geo file.

        """

        return 'Surface(' + str(tag) + ') = {' + str(line_loop_tag)[1:-1] + '};'

    def arc_inside(self, arc: vme.Arc3D):
        """
        Verifies if Arc3D is inside a CylindricalFace3D.

        :param arc: Arc3D to be verified.
        :return: True if it is inside, False otherwise.
        """
        if not math.isclose(abs(arc.frame.w.dot(self.surface3d.frame.w)), 1.0, abs_tol=1e-6):
            return False
        if not math.isclose(self.radius, arc.radius, abs_tol=1e-6):
            return False
        return self.arcellipse_inside(arc)

    def arcellipse_inside(self, arcellipse: vme.ArcEllipse3D):
        """
        Verifies if ArcEllipse3D is inside a CylindricalFace3D.

        :param arcellipse: ArcEllipse3D to be verified.
        :return: True if it is inside, False otherwise.
        """
        for point in arcellipse.points:
            if not self.point_belongs(point):
                return False
        return True

    def planeface_intersections(self, planeface: PlaneFace3D):
        planeface_intersections = planeface.cylindricalface_intersections(self)
        return planeface_intersections


class ToroidalFace3D(Face3D):
    """
    Defines a ToroidalFace3D class.

    :param surface3d: a toroidal surface 3d.
    :type surface3d: ToroidalSurface3D.
    :param surface2d: a 2d surface to define the toroidal face.
    :type surface2d: Surface2D.

    :Example:

    contours 2d is rectangular and will create a classic tore with x:2*pi, y:2*pi
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

        theta_min, theta_max, phi_min, phi_max = surface2d.outer_contour.bounding_rectangle.bounds()

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
        self._bbox = None

    def copy(self, deep=True, memo=None):
        return ToroidalFace3D(self.surface3d.copy(deep, memo), self.surface2d.copy(),
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

    @property
    def bounding_box(self):
        """
        Returns the face bounding box.
        """
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def get_bounding_box(self):
        return self.surface3d.bounding_box

    def triangulation_lines(self, angle_resolution=5):
        theta_min, theta_max, phi_min, phi_max = self.surface2d.bounding_rectangle().bounds()

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

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        theta_angle_resolution = 11
        phi_angle_resolution = 7
        theta_min, theta_max, phi_min, phi_max = self.surface2d.bounding_rectangle().bounds()

        delta_theta = theta_max - theta_min
        number_points_x = int(delta_theta * theta_angle_resolution)

        delta_phi = phi_max - phi_min
        number_points_y = int(delta_phi * phi_angle_resolution)

        return number_points_x, number_points_y


class ConicalFace3D(Face3D):
    """
    Defines a ConicalFace3D class.

    :param surface3d: a conical surface 3d.
    :type surface3d: ConicalSurface3D.
    :param surface2d: a 2d surface to define the conical face.
    :type surface2d: Surface2D.


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
        self._bbox = None

    @property
    def bounding_box(self):
        """
        Surface bounding box.
        """
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bouding_box):
        self._bbox = new_bouding_box

    def get_bounding_box(self):
        theta_min, theta_max, zmin, zmax = self.surface2d.outer_contour.bounding_rectangle.bounds()

        x_vector = (volmdlr.X3D.dot(self.surface3d.frame.u) * self.surface3d.frame.u
                    + volmdlr.X3D.dot(self.surface3d.frame.v) * self.surface3d.frame.v)
        try:
            x_vector.normalize()
        except ZeroDivisionError:
            pass
        y_vector = (volmdlr.Y3D.dot(self.surface3d.frame.u) * self.surface3d.frame.u
                    + volmdlr.Y3D.dot(self.surface3d.frame.v) * self.surface3d.frame.v)

        try:
            y_vector.normalize()
        except ZeroDivisionError:
            pass

        z_vector = (volmdlr.Z3D.dot(self.surface3d.frame.u) * self.surface3d.frame.u
                    + volmdlr.Z3D.dot(self.surface3d.frame.v) * self.surface3d.frame.v)
        try:
            z_vector.normalize()
        except ZeroDivisionError:
            pass

        lower_center = self.surface3d.frame.origin + zmin * self.surface3d.frame.w
        upper_center = self.surface3d.frame.origin + zmax * self.surface3d.frame.w
        lower_radius = math.tan(self.surface3d.semi_angle) * zmin
        upper_radius = math.tan(self.surface3d.semi_angle) * zmax

        points = [lower_center - lower_radius * x_vector,
                  lower_center + lower_radius * x_vector,
                  lower_center - lower_radius * y_vector,
                  lower_center + lower_radius * y_vector,
                  lower_center - lower_radius * z_vector,
                  lower_center + lower_radius * z_vector,
                  upper_center - upper_radius * x_vector,
                  upper_center + upper_radius * x_vector,
                  upper_center - upper_radius * y_vector,
                  upper_center + upper_radius * y_vector,
                  upper_center - upper_radius * z_vector,
                  upper_center + upper_radius * z_vector,
                  ]

        return volmdlr.core.BoundingBox.from_points(points)

    def triangulation_lines(self, angle_resolution=5):
        theta_min, theta_max, zmin, zmax = self.surface2d.bounding_rectangle().bounds()
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

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        angle_resolution = 5
        theta_min, theta_max, _, _ = self.surface2d.bounding_rectangle().bounds()
        delta_theta = theta_max - theta_min
        number_points_x = math.ceil(delta_theta * angle_resolution)

        number_points_y = 0

        return number_points_x, number_points_y


class SphericalFace3D(Face3D):
    """
    Defines a SpehericalFace3D class.

    :param surface3d: a spherical surface 3d.
    :type surface3d: SphericalSurface3D.
    :param surface2d: a 2d surface to define the spherical face.
    :type surface2d: Surface2D.


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
        self._bbox = None

    @property
    def bounding_box(self):
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bouding_box):
        self._bbox = new_bouding_box

    def get_bounding_box(self):
        # To be enhanced
        return self.surface3d.bounding_box

    def triangulation_lines(self, angle_resolution=7):
        """
        Specifies the number of subdivision when using triangulation by lines. (Old triangulation).
        """
        theta_min, theta_max, phi_min, phi_max = self.surface2d.bounding_rectangle().bounds()

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

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        angle_resolution = 11
        theta_min, theta_max, phi_min, phi_max = self.surface2d.bounding_rectangle().bounds()

        delta_theta = theta_max - theta_min
        number_points_x = int(delta_theta * angle_resolution)

        delta_phi = phi_max - phi_min
        number_points_y = int(delta_phi * angle_resolution)

        return number_points_x, number_points_y


class RuledFace3D(Face3D):
    """
    A 3D face with a ruled surface.

    This class represents a 3D face with a ruled surface, which is a surface
    formed by straight lines connecting two input curves. It is a subclass of
    the `Face3D` class and inherits all of its attributes and methods.

    :param surface3d: The 3D ruled surface of the face.
    :type surface3d: `RuledSurface3D`
    :param surface2d: The 2D projection of the face onto the parametric domain (u, v).
    :type surface2d: `Surface2D`
    :param name: The name of the face.
    :type name: str
    :param color: The color of the face.
    :type color: tuple
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
        self._bbox = None

    @property
    def bounding_box(self):
        """
        Returns the bounding box of the surface.
        """
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bouding_box):
        self._bbox = new_bouding_box

    def get_bounding_box(self):
        # To be enhance by restricting wires to cut
        # xmin, xmax, ymin, ymax = self.surface2d.outer_contour.bounding_rectangle()
        points = [self.surface3d.point2d_to_3d(volmdlr.Point2D(i / 30, 0.)) for
                  i in range(31)]
        points.extend(
            [self.surface3d.point2d_to_3d(volmdlr.Point2D(i / 30, 1.)) for i
             in range(31)])

        return volmdlr.core.BoundingBox.from_points(points)

    def triangulation_lines(self, angle_resolution=10):
        """
        Specifies the number of subdivision when using triangulation by lines. (Old triangulation).
        """
        xmin, xmax, ymin, ymax = self.surface2d.bounding_rectangle().bounds()
        delta_x = xmax - xmin
        nlines = int(delta_x * angle_resolution)
        lines = []
        for i in range(nlines):
            x = xmin + (i + 1) / (nlines + 1) * delta_x
            lines.append(vme.Line2D(volmdlr.Point2D(x, ymin),
                                    volmdlr.Point2D(x, ymax)))
        return lines, []

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        angle_resolution = 10
        xmin, xmax, _, _ = self.surface2d.bounding_rectangle().bounds()
        delta_x = xmax - xmin
        number_points_x = int(delta_x * angle_resolution)

        number_points_y = 0

        return number_points_x, number_points_y


class ExtrusionFace3D(Face3D):
    """
    A 3D face with a ruled surface.

    This class represents a 3D face with a ruled surface, which is a surface
    formed by straight lines connecting two input curves. It is a subclass of
    the `Face3D` class and inherits all of its attributes and methods.


    :param surface3d: The 3D ruled surface of the face.
    :type surface3d: `RuledSurface3D`
    :param surface2d: The 2D projection of the face onto the parametric domain (u, v).
    :type surface2d: `Surface2D`
    :param name: The name of the face.
    :type name: str
    """
    min_x_density = 50
    min_y_density = 1

    def __init__(self,
                 surface3d: RuledSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):
        Face3D.__init__(self, surface3d=surface3d,
                        surface2d=surface2d,
                        name=name)
        self._bbox = None

    @property
    def bounding_box(self):
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bouding_box):
        self._bbox = new_bouding_box

    def get_bounding_box(self):
        # To be enhanced by restricting wires to cut
        points = self.outer_contour3d.discretization_points(number_points=25)

        return volmdlr.core.BoundingBox.from_points(points)

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        angle_resolution = 36
        xmin, xmax, _, _ = self.surface2d.bounding_rectangle().bounds()
        delta_x = xmax - xmin
        number_points_x = int(delta_x * angle_resolution)

        number_points_y = number_points_x

        return number_points_x, number_points_y


class RevolutionFace3D(Face3D):
    """
    A 3D face with a ruled surface.

    This class represents a 3D face with a ruled surface, which is a surface
    formed by straight lines connecting two input curves. It is a subclass of
    the `Face3D` class and inherits all of its attributes and methods.


    :param surface3d: The 3D ruled surface of the face.
    :type surface3d: `RuledSurface3D`
    :param surface2d: The 2D projection of the face onto the parametric domain (u, v).
    :type surface2d: `Surface2D`
    :param name: The name of the face.
    :type name: str
    """
    min_x_density = 50
    min_y_density = 1

    def __init__(self,
                 surface3d: RuledSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):
        Face3D.__init__(self, surface3d=surface3d,
                        surface2d=surface2d,
                        name=name)
        self._bbox = None

    @property
    def bounding_box(self):
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bouding_box):
        self._bbox = new_bouding_box

    def get_bounding_box(self):
        # To be enhanced by restricting wires to cut
        curve_points = self.surface3d.edge.discretization_points(angle_resolution=20)
        points = []
        for i in range(37):
            angle = i * volmdlr.TWO_PI / 36
            points.extend([point.rotation(self.surface3d.axis_point, self.surface3d.axis, angle)
                           for point in curve_points])

        return volmdlr.core.BoundingBox.from_points(points)

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        angle_resolution = 10
        xmin, xmax, _, _ = self.surface2d.bounding_rectangle().bounds()
        delta_x = xmax - xmin
        number_points_x = int(delta_x * angle_resolution)

        number_points_y = number_points_x

        return number_points_x, number_points_y


class BSplineFace3D(Face3D):
    """
    A 3D face with a B-spline surface.

    This class represents a 3D face with a B-spline surface, which is a smooth
    surface defined by a set of control points and knots. It is a subclass of
    the `Face3D` class and inherits all of its attributes and methods.

    :param surface3d: The 3D B-spline surface of the face.
    :type surface3d: `BSplineSurface3D`
    :param surface2d: The 2D projection of the face onto the parametric domain (u, v).
    :type surface2d: `Surface2D`
    :param name: The name of the face.
    :type name: str
    """

    def __init__(self, surface3d: BSplineSurface3D,
                 surface2d: Surface2D,
                 name: str = ''):
        Face3D.__init__(self,
                        surface3d=surface3d,
                        surface2d=surface2d,
                        name=name)
        self._bbox = None

    @property
    def bounding_box(self):
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def get_bounding_box(self):
        return self.surface3d.bounding_box

    def triangulation_lines(self, resolution=25):
        u_min, u_max, v_min, v_max = self.surface2d.bounding_rectangle().bounds()

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

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        if self.surface3d.x_periodicity or self.surface3d.y_periodicity:
            resolution = 25
        else:
            resolution = 15
        u_min, u_max, v_min, v_max = self.surface2d.bounding_rectangle().bounds()
        delta_u = u_max - u_min
        number_points_x = int(delta_u * resolution)

        delta_v = v_max - v_min
        number_points_y = int(delta_v * resolution)

        return number_points_x, number_points_y

    def pair_with(self, other_bspline_face3d):
        """
        Finds out how the uv parametric frames are located.

        It does it by compaing to each other and also how grid 3d can be defined respected to these directions.

        :param other_bspline_face3d: BSplineFace3D
        :type other_bspline_face3d: :class:`volmdlr.faces.BSplineFace3D`
        :return: corresponding_direction, grid2d_direction
        :rtype: Tuple[?, ?]
        """

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

        if v1 == 0 and v2 == 0:
            corresponding_directions.append(('+v', '-v'))
            grid2d_direction = [['+x', '-y'], ['+x', '+y']]

        elif v1 == 1 and v2 == 1:
            if corresponding_directions == [('+u', '-u')]:
                grid2d_direction = [['+x', '+y'], ['-x', '-y']]
            else:
                grid2d_direction = [['+x', '+y'], ['+x', '-y']]
            corresponding_directions.append(('+v', '-v'))

        elif v1 == 1 and v2 == 0:
            corresponding_directions.append(('+v', '+v'))
            grid2d_direction = [['+x', '+y'], ['+x', '+y']]

        elif v1 == 0 and v2 == 1:
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

        if u1 == 0 and u2 == 0:
            corresponding_directions.append(('+u', '-v'))
            grid2d_direction = [['-y', '-x'], ['-y', '+x']]

        elif u1 == 1 and u2 == 1:
            if corresponding_directions == [('+v', '-v')]:
                grid2d_direction = [['+y', '+x'], ['-y', '-x']]
            else:
                grid2d_direction = [['+y', '+x'], ['+y', '-x']]
            corresponding_directions.append(('+u', '-u'))

        elif u1 == 0 and u2 == 1:
            corresponding_directions.append(('+u', '+u'))
            grid2d_direction = [['+y', '-x'], ['+y', '-x']]

        elif u1 == 1 and u2 == 0:
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

        if v1 == 1 and u2 == 0:
            corresponding_directions.append(('+v', '+u'))
            grid2d_direction = [['+x', '+y'], ['+y', '+x']]

        elif v1 == 0 and u2 == 1:
            corresponding_directions.append(('+v', '+u'))
            grid2d_direction = [['-x', '-y'], ['-y', '-x']]

        elif v1 == 1 and u2 == 1:
            corresponding_directions.append(('+v', '-u'))
            grid2d_direction = [['+x', '+y'], ['-y', '-x']]

        elif v1 == 0 and u2 == 0:
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

        if u1 == 1 and v2 == 0:
            corresponding_directions.append(('+u', '+v'))
            grid2d_direction = [['+y', '+x'], ['+x', '+y']]

        elif u1 == 0 and v2 == 1:
            corresponding_directions.append(('+u', '+v'))
            grid2d_direction = [['-y', '-x'], ['+x', '-y']]

        elif u1 == 0 and v2 == 0:
            corresponding_directions.append(('+u', '-v'))
            grid2d_direction = [['+y', '-x'], ['+x', '+y']]

        elif u1 == 1 and v2 == 1:
            if corresponding_directions == [('+v', '-u')]:
                grid2d_direction = [['+y', '+x'], ['-x', '-y']]
            else:
                grid2d_direction = [['+y', '+x'], ['+x', '-y']]
            corresponding_directions.append(('+u', '-v'))

        return corresponding_directions, grid2d_direction

    def extremities(self, other_bspline_face3d):
        """
        Find points extremities for nearest edges of two faces.
        """
        contour1 = self.outer_contour3d
        contour2 = other_bspline_face3d.outer_contour3d

        contour1_2d = self.surface2d.outer_contour
        contour2_2d = other_bspline_face3d.surface2d.outer_contour

        points1 = [prim.start for prim in contour1.primitives]
        points2 = [prim.start for prim in contour2.primitives]

        dis, ind = [], []
        for point_ in points1:
            point = point_.nearest_point(points2)
            ind.append(points2.index(point))
            dis.append(point_.point_distance(point))

        dis_sorted = sorted(dis)

        shared = []
        for k, point1 in enumerate(contour1.primitives):
            if dis_sorted[0] == dis_sorted[1]:
                indices = npy.where(npy.array(dis) == dis_sorted[0])[0]
                index1 = indices[0]
                index2 = indices[1]
            else:
                index1 = dis.index(dis_sorted[0])
                index2 = dis.index(dis_sorted[1])
            if ((point1.start.is_close(points1[index1]) and point1.end.is_close(points1[index2]))
                    or
                    (point1.end.is_close(points1[index1]) and point1.start.is_close(points1[index2]))):
                shared.append(point1)
                i = k

        for k, prim2 in enumerate(contour2.primitives):
            if ((prim2.start.is_close(points2[ind[index1]]) and prim2.end.is_close(points2[ind[index2]]))
                    or
                    (prim2.end.is_close(points2[ind[index1]]) and prim2.start.is_close(points2[ind[index2]]))):
                shared.append(prim2)
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
        """
        Find directions (u or v) between two faces, in the nearest edges between them.
        """

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
        """
        Find out in which direction the faces are adjacent.

        :type other_face3d: volmdlr.faces.BSplineFace3D
        :return: adjacent_direction
        """

        contour1 = self.outer_contour3d
        contour2 = other_face3d.outer_contour3d
        point1, point2 = contour1.shared_primitives_extremities(contour2)

        coord = point1 - point2
        coord = [abs(coord.x), abs(coord.y)]

        if coord.index(max(coord)) == 0:
            return 'x'
        return 'y'

    def merge_with(self, other_bspline_face3d):
        """
        Merge two adjacent faces.

        :type: other_bspline_face3d : volmdlr.faces.BSplineFace3D
        :rtype: merged_face : volmdlr.faces.BSplineFace3D
        """

        merged_surface = self.surface3d.merge_with(other_bspline_face3d.surface3d)
        contours = self.outer_contour3d.merge_with(other_bspline_face3d.outer_contour3d)
        contours.extend(self.inner_contours3d)
        contours.extend(other_bspline_face3d.inner_contours3d)
        merged_face = merged_surface.face_from_contours3d(contours)

        return merged_face


class OpenShell3D(volmdlr.core.CompositePrimitive3D):
    """
    A 3D open shell composed of multiple faces.

    This class represents a 3D open shell, which is a collection of connected
    faces with no volume. It is a subclass of the `CompositePrimitive3D` class
    and inherits all of its attributes and methods.


    :param faces: The faces of the shell.
    :type faces: List[`Face3D`]
    :param color: The color of the shell.
    :type color: Tuple[float, float, float]
    :param alpha: The transparency of the shell, should be a value in the interval (0, 1).
    :type alpha: float
    :param name: The name of the shell.
    :type name: str
    :param bounding_box: The bounding box of the shell.
    :type bounding_box: :class:`volmdlr.core.BoundingBox`
    """
    _standalone_in_db = True
    _non_serializable_attributes = ['primitives']
    _non_data_eq_attributes = ['name', 'color', 'alpha', 'bounding_box', 'primitives']
    _non_data_hash_attributes = []
    STEP_FUNCTION = 'OPEN_SHELL'

    def __init__(self, faces: List[Face3D],
                 color: Tuple[float, float, float] = None,
                 alpha: float = 1.,
                 name: str = '',
                 bounding_box: volmdlr.core.BoundingBox = None):

        self.faces = faces
        if not color:
            self.color = volmdlr.core.DEFAULT_COLOR
        else:
            self.color = color
        self.alpha = alpha

        if bounding_box:
            self._bbox = bounding_box
        else:
            self._bbox = None

        self._faces_graph = None

        volmdlr.core.CompositePrimitive3D.__init__(self,
                                                   primitives=faces, color=color, alpha=alpha,
                                                   name=name)

    def _data_hash(self):
        return len(self.faces)  # sum(face._data_hash() for face in self.faces)

    def _data_eq(self, other_object):
        if other_object.__class__.__name__ != self.__class__.__name__:
            return False
        for face1, face2 in zip(self.faces, other_object.faces):
            if not face1._data_eq(face2):
                return False

        return True

    @property
    def faces_graph(self):
        if not self._faces_graph:
            faces_graph = nx.Graph()
            for face in self.faces:
                for edge in face.outer_contour3d.primitives:
                    faces_graph.add_edge(edge.start, edge.end, edge=edge)
            self._faces_graph = faces_graph
        return self._faces_graph

    def to_dict(self, *args, **kwargs):
        """
        Serializes a 3 dimensional open shell into a dictionary.

        This method does not use pointers for faces as it has no sense
        to have duplicate faces.

        :return: A serialized version of the OpenShell3D
        :rtype: dict

        .. see also::
            How `serialization and de-serialization`_ works in dessia_common

        .. _serialization and deserialization:
        https://documentation.dessia.tech/dessia_common/customizing.html#overloading-the-dict-to-object-method

        """
        dict_ = DessiaObject.base_dict(self)
        dict_.update({'color': self.color,
                      'alpha': self.alpha,
                      'faces': [f.to_dict(use_pointers=False) for f in self.faces]})
        if self._bbox:
            dict_['bounding_box'] = self._bbox.to_dict()

        return dict_

    @classmethod
    def from_step(cls, arguments, object_dict, **kwargs):
        """
        Converts a step primitive to a Open Shell 3D.

        :param arguments: The arguments of the step primitive.
        :type arguments: list
        :param object_dict: The dictionary containing all the step primitives
            that have already been instantiated.
        :type object_dict: dict
        :return: The corresponding OpenShell3D object.
        :rtype: :class:`volmdlr.faces.OpenShell3D`
        """
        faces = [object_dict[int(face[1:])] for face in arguments[1] if object_dict[int(face[1:])] is not None]
        return cls(faces, name=arguments[0][1:-1])

    def to_step(self, current_id):
        """
        Creates step file entities from volmdlr objects.
        """
        step_content = ''
        face_ids = []
        for face in self.faces:
            if isinstance(face, (Face3D, Surface3D)):
                face_content, face_sub_ids = face.to_step(current_id)
            else:
                face_content, face_sub_ids = face.to_step(current_id)
                face_sub_ids = [face_sub_ids]
            step_content += face_content
            face_ids.extend(face_sub_ids)
            current_id = max(face_sub_ids) + 1

        shell_id = current_id
        step_content += f"#{current_id} = {self.STEP_FUNCTION}('{self.name}'," \
                        f"({volmdlr.core.step_ids_to_str(face_ids)}));\n"
        manifold_id = shell_id + 1
        step_content += f"#{manifold_id} = SHELL_BASED_SURFACE_MODEL('{self.name}',(#{shell_id}));\n"

        frame_content, frame_id = volmdlr.OXYZ.to_step(manifold_id + 1)
        step_content += frame_content
        brep_id = frame_id + 1
        step_content += f"#{brep_id} = MANIFOLD_SURFACE_SHAPE_REPRESENTATION('',(#{frame_id},#{manifold_id}),#7);\n"

        return step_content, brep_id

    def to_step_face_ids(self, current_id):
        """
        Creates step file entities from volmdlr objects.
        """
        step_content = ''
        face_ids = []
        for face in self.faces:
            if isinstance(face, Face3D):
                face_content, face_sub_ids = face.to_step(current_id)
            else:
                face_content, face_sub_ids = face.to_step(current_id)
                face_sub_ids = [face_sub_ids]
            step_content += face_content
            face_ids.extend(face_sub_ids)
            current_id = max(face_sub_ids) + 1

        shell_id = current_id
        step_content += f"#{current_id} = {self.STEP_FUNCTION}('{self.name}'," \
                        f"({volmdlr.core.step_ids_to_str(face_ids)}));\n"
        manifold_id = shell_id + 1
        step_content += f"#{manifold_id} = SHELL_BASED_SURFACE_MODEL('{self.name}',(#{shell_id}));\n"

        frame_content, frame_id = volmdlr.OXYZ.to_step(manifold_id + 1)
        step_content += frame_content
        brep_id = frame_id + 1
        step_content += f"#{brep_id} = MANIFOLD_SURFACE_SHAPE_REPRESENTATION('',(#{frame_id},#{manifold_id}),#7);\n"

        return step_content, brep_id, face_ids

    def rotation(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                 angle: float):
        """
        Open Shell 3D / Closed Shell 3D rotation.

        :param center: rotation center.
        :param axis: rotation axis.
        :param angle: angle rotation.
        :return: a new rotated OpenShell3D.
        """
        new_faces = [face.rotation(center, axis, angle) for face
                     in self.faces]
        return self.__class__(new_faces, color=self.color, alpha=self.alpha, name=self.name)

    def rotation_inplace(self, center: volmdlr.Point3D, axis: volmdlr.Vector3D,
                         angle: float):
        """
        Shell 3D rotation. Object is updated inplace.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for face in self.faces:
            face.rotation_inplace(center, axis, angle)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box

    def translation(self, offset: volmdlr.Vector3D):
        """
        Shell3D translation.

        :param offset: translation vector.
        :return: A new translated Open Shell 3D.
        """
        new_faces = [face.translation(offset) for face in
                     self.faces]
        return self.__class__(new_faces, color=self.color, alpha=self.alpha,
                              name=self.name)

    def translation_inplace(self, offset: volmdlr.Vector3D):
        """
        Open Shell 3D translation. Object is updated inplace.

        :param offset: Translation vector.
        :type offset: `volmdlr.Vector3D`.
        :return: Translate the Open Shell 3D in place.
        :rtype: None.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for face in self.faces:
            face.translation_inplace(offset)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new OpenShell3D.

        side = 'old' or 'new'.
        """
        new_faces = [face.frame_mapping(frame, side) for face in
                     self.faces]
        return self.__class__(new_faces, name=self.name)

    def frame_mapping_inplace(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and the object is updated inplace.

        side = 'old' or 'new'.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not inplace method instead.", DeprecationWarning)

        for face in self.faces:
            face.frame_mapping_inplace(frame, side)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box

    def copy(self, deep=True, memo=None):
        new_faces = [face.copy(deep=deep, memo=memo) for face in self.faces]
        return self.__class__(new_faces, color=self.color, alpha=self.alpha,
                              name=self.name)

    def union(self, shell2):
        new_faces = self.faces + shell2.faces
        new_name = self.name + ' union ' + shell2.name
        new_color = self.color
        return self.__class__(new_faces, name=new_name, color=new_color)

    def volume(self):
        """
        Does not consider holes.

        """
        volume = 0
        for face in self.faces:
            display3d = face.triangulation()
            for triangle_index in display3d.triangles:
                point1 = display3d.points[triangle_index[0]]
                point2 = display3d.points[triangle_index[1]]
                point3 = display3d.points[triangle_index[2]]

                v321 = point3[0] * point2[1] * point1[2]
                v231 = point2[0] * point3[1] * point1[2]
                v312 = point3[0] * point1[1] * point2[2]
                v132 = point1[0] * point3[1] * point2[2]
                v213 = point2[0] * point1[1] * point3[2]
                v123 = point1[0] * point2[1] * point3[2]
                volume_tetraedre = 1 / 6 * (-v321 + v231 + v312 - v132 - v213 + v123)

                volume += volume_tetraedre

        return abs(volume)

    @property
    def bounding_box(self):
        """
        Returns the boundary box.

        """
        if not self._bbox:
            self._bbox = self.get_bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def get_bounding_box(self):
        return volmdlr.core.BoundingBox.from_bounding_boxes([face.bounding_box for face in self.faces])

    def cut_by_plane(self, plane_3d: Plane3D):
        frame_block = self.bounding_box.to_frame()
        frame_block.u = 1.1 * frame_block.u
        frame_block.v = 1.1 * frame_block.v
        frame_block.w = 1.1 * frame_block.w
        block = volmdlr.primitives3d.Block(frame_block,
                                           color=(0.1, 0.2, 0.2),
                                           alpha=0.6)
        face_3d = block.cut_by_orthogonal_plane(plane_3d)
        intersection_primitives = []
        for face in self.faces:
            intersection_wires = face.face_intersections(face_3d)
            if intersection_wires:
                for intersection_wire in intersection_wires:
                    intersection_primitives.extend(intersection_wire.primitives)
        contours3d = volmdlr.wires.Contour3D.contours_from_edges(
            intersection_primitives[:])
        if not contours3d:
            return []
        contours2d = [contour.to_2d(plane_3d.frame.origin,
                                    plane_3d.frame.u,
                                    plane_3d.frame.v) for contour in contours3d]
        resulting_faces = []
        for contour2d in contours2d:
            if contour2d.area() > 1e-7:
                surface2d = Surface2D(contour2d, [])
                resulting_faces.append(PlaneFace3D(plane_3d, surface2d))
        return resulting_faces

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
        Returns a Measure object if the distance is not zero, otherwise returns None.

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
                    if distance < distance_min:
                        distance_min, point1_min, point2_min = distance, point1, point2

        return point1_min, point2_min

    def distance_to_shell(self, other_shell: 'OpenShell3D', resolution: float):
        min_dist = self.minimum_distance_points(other_shell, resolution)
        if min_dist is not None:
            point1, point2 = min_dist
            return point1.point_distance(point2)
        return 0

    def minimum_distance_point(self,
                               point: volmdlr.Point3D) -> volmdlr.Point3D:
        """
        Computes the distance of a point to a Shell3D, whether it is inside or outside the Shell3D.

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
        Aabb made of the intersection points and the points of self internal to shell2.
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
            for point in face.outer_contour3d.discretization_points(angle_resolution=resolution):
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
        Aabb made of the intersection points and the points of self external to shell2.
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
                    angle_resolution=resolution):
                if not shell2.point_belongs(point):
                    shell1_points_outside_shell2.append(point)

        if len(intersections_points + shell1_points_outside_shell2) == 0:
            return 0
        bbox = volmdlr.core.BoundingBox.from_points(
            intersections_points + shell1_points_outside_shell2)
        return bbox.volume()

    def face_on_shell(self, face):
        """
        Verifies if a face lies on the shell's surface.

        """
        for face_ in self.faces:
            if face_.face_inside(face):
                return True
        return False

    def point_on_shell(self, point: volmdlr.Point3D):
        for face in self.faces:
            if face.point_belongs(point) or face.outer_contour3d.point_over_contour(point, abs_tol=1e-7):
                return True
        return False

    def point_in_shell_face(self, point: volmdlr.Point3D):
        warnings.warn('point_in_shell_face is deprecated, please use point_on_shell instead',
                      DeprecationWarning)
        return self.point_on_shell(point)

    def triangulation(self):
        meshes = []
        for face in self.faces:
            face_mesh = face.triangulation()
            meshes.append(face_mesh)
        return vmd.DisplayMesh3D.merge_meshes(meshes)

    def plot(self, ax=None, color: str = 'k', alpha: float = 1.0):
        if ax is None:
            ax = plt.figure().add_subplot(111, projection='3d')

        for face in self.faces:
            face.plot(ax=ax, color=color, alpha=alpha)

        return ax

    def project_coincident_faces_of(self, shell):
        """
        Divides self's faces based on coincident shell's faces.

        """

        list_faces = []
        initial_faces = self.faces[:]

        for face1 in initial_faces:
            list_faces.extend(face1.project_faces(shell.faces))

        return self.__class__(list_faces)

    def get_geo_lines(self, update_data,
                      point_mesh_size: float = None):
        """
        Gets the lines that define an OpenShell3D geometry in a .geo file.

        :param update_data: Data used for VolumeModel defined with different shells
        :type update_data: dict
        :param point_mesh_size: The mesh size at a specific point, defaults to None
        :type point_mesh_size: float, optional

        :return: A list of lines that describe the geometry & the updated data
        :rtype: Tuple(List[str], dict)
        """

        primitives = []
        points = set()
        for face in self.faces:
            for _, contour in enumerate(list(chain(*[[face.outer_contour3d], face.inner_contours3d]))):
                points.update(contour.get_geo_points())
                if isinstance(contour, volmdlr.wires.Circle2D):
                    pass
                else:
                    for _, primitive in enumerate(contour.primitives):
                        if ((primitive not in primitives)
                                and (primitive.reverse() not in primitives)):
                            primitives.append(primitive)

                # if isinstance(contour, volmdlr.wires.Circle2D):
                #     points.add(volmdlr.Point3D(contour.radius, contour.center.y, 0))
                #     points.add(volmdlr.Point3D(contour.center.x, contour.center.y, 0))
                #     points.add(volmdlr.Point3D(-contour.radius, contour.center.y, 0))

                # else:
                #     for _, primitive in enumerate(contour.primitives):
                #         if isinstance(primitive, volmdlr.edges.LineSegment):
                #             points.add(primitive.start)
                #             points.add(primitive.end)

                #         if isinstance(primitive, volmdlr.edges.Arc):
                #             points.add(primitive.start)
                #             points.add(primitive.center)
                #             points.add(primitive.end)

                #         if isinstance(primitive, volmdlr.edges.BSplineCurve3D):
                #             # for point in primitive.control_points:
                #             # points.add(point)
                #             for point in primitive.discretization_points():
                #                 points.add(point)

                #         if ((primitive not in primitives)
                #                 and (primitive.reverse() not in primitives)):
                #             primitives.append(primitive)

        indices_check = len(primitives) * [None]

        point_account = update_data['point_account']
        line_account, line_loop_account = update_data['line_account'] + 1, update_data['line_loop_account']
        lines, line_surface, lines_tags = [], [], []

        points = list(points)
        for p_index, point in enumerate(points):
            lines.append(point.get_geo_lines(tag=p_index + point_account + 1,
                                             point_mesh_size=point_mesh_size))

        for f_index, face in enumerate(self.faces):
            line_surface = []
            for _, contour in enumerate(list(chain(*[[face.outer_contour3d], face.inner_contours3d]))):
                lines_tags = []
                if isinstance(contour, volmdlr.wires.Circle2D):
                    pass
                else:
                    for _, primitive in enumerate(contour.primitives):

                        try:
                            # line_account += 1
                            # print(line_account)
                            index = primitives.index(primitive)
                            if isinstance(primitive, volmdlr.edges.BSplineCurve3D):
                                discretization_points = primitive.discretization_points()
                                start_point_tag = points.index(discretization_points[0]) + 1
                                end_point_tag = points.index(discretization_points[1]) + 1
                                primitive_linesegments = volmdlr.edges.LineSegment3D(
                                    discretization_points[0], discretization_points[1])
                                lines.append(primitive_linesegments.get_geo_lines(tag=line_account,
                                                                                  start_point_tag=start_point_tag
                                                                                                  + point_account,
                                                                                  end_point_tag=end_point_tag
                                                                                                + point_account))

                            if isinstance(primitive, volmdlr.edges.LineSegment):
                                start_point_tag = points.index(primitive.start) + 1
                                end_point_tag = points.index(primitive.end) + 1
                                lines.append(primitive.get_geo_lines(tag=line_account,
                                                                     start_point_tag=start_point_tag + point_account,
                                                                     end_point_tag=end_point_tag + point_account))
                            elif isinstance(primitive, volmdlr.edges.Arc):
                                start_point_tag = points.index(primitive.start) + 1
                                center_point_tag = points.index(primitive.center) + 1
                                end_point_tag = points.index(primitive.end) + 1
                                lines.append(primitive.get_geo_lines(tag=line_account,
                                                                     start_point_tag=start_point_tag + point_account,
                                                                     center_point_tag=center_point_tag + point_account,
                                                                     end_point_tag=end_point_tag + point_account))

                            lines_tags.append(line_account)
                            indices_check[index] = line_account
                            line_account += 1

                        except ValueError:
                            index = primitives.index(primitive.reverse())
                            lines_tags.append(-indices_check[index])

                    lines.append(contour.get_geo_lines(line_loop_account + 1, lines_tags))

                    line_surface.append(line_loop_account + 1)
                    line_loop_account += 1
                    lines_tags = []

            lines.append(face.get_geo_lines((f_index + 1 + update_data['surface_account']),
                                            line_surface))

            line_surface = []

        lines.append('Surface Loop(' + str(1 + update_data['surface_loop_account']) + ') = {'
                     + str(list(range(update_data['surface_account'] + 1,
                                      update_data['surface_account'] +
                                      len(self.faces) + 1)))[1:-1] + '};')

        update_data['point_account'] += len(points)
        update_data['line_account'] += line_account - 1
        update_data['line_loop_account'] += line_loop_account
        update_data['surface_account'] += len(self.faces)
        update_data['surface_loop_account'] += 1

        return lines, update_data

    def get_mesh_lines_with_transfinite_curves(self, min_points, size):

        lines, primitives, primitives_length = [], [], []
        for face in self.faces:
            for _, contour in enumerate(list(chain(*[[face.outer_contour3d], face.inner_contours3d]))):
                if isinstance(contour, volmdlr.wires.Circle2D):
                    primitives.append(contour)
                    primitives.append(contour)
                    primitives_length.append(contour.length() / 2)
                    primitives_length.append(contour.length() / 2)
                else:
                    for _, primitive_c in enumerate(contour.primitives):
                        if ((primitive_c not in primitives)
                                and (primitive_c.reverse() not in primitives)):
                            primitives.append(primitive_c)
                            primitives_length.append(primitive_c.length())

        for i, length in enumerate(primitives_length):
            if length < min_points * size:
                lines.append('Transfinite Curve {' + str(i) + '} = ' +
                             str(min_points) + ' Using Progression 1;')
        return lines


class ClosedShell3D(OpenShell3D):
    """
    A 3D closed shell composed of multiple faces.

    This class represents a 3D closed shell, which is a collection of connected
    faces with a volume. It is a subclass of the `OpenShell3D` class and
    inherits all of its attributes and methods. In addition, it has a method
    to check whether a face is inside the shell.

    :param faces: The faces of the shell.
    :type faces: List[`Face3D`]
    :param color: The color of the shell.
    :type color: Tuple[float, float, float]
    :param alpha: The transparency of the shell, should be a value in the range (0, 1).
    :type alpha: float
    :param name: The name of the shell.
    :type name: str
    """

    STEP_FUNCTION = 'CLOSED_SHELL'

    def is_face_inside(self, face: Face3D):
        for point in face.outer_contour3d.discretization_points(angle_resolution=0.1):
            point_inside_shell = self.point_belongs(point)
            point_in_shells_faces = self.point_on_shell(point)
            if (not point_inside_shell) and (not point_in_shells_faces):
                return False
        return True

    def shell_intersection(self, shell2: 'OpenShell3D', resolution: float):
        """
        Return None if disjointed.

        Return (1, 0) or (0, 1) if one is inside the other
        Return (n1, n2) if intersection

        4 cases :
            (n1, n2) with face intersection             => (n1, n2)
            (0, 0) with face intersection               => (0, 0)
            (0, 0) with no face intersection            => None
            (1, 0) or (0, 1) with no face intersection  => 1
        """
        # Check if boundary boxes don't intersect
        if not self.bounding_box.bbox_intersection(shell2.bounding_box):
            # print("No intersection of shells' BBox")
            return None

        # Check if any point of the first shell is in the second shell
        points1 = []
        for face in self.faces:
            points1.extend(
                face.outer_contour3d.discretization_points(angle_resolution=resolution))
        points2 = []
        for face in shell2.faces:
            points2.extend(
                face.outer_contour3d.discretization_points(angle_resolution=resolution))

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

    def point_belongs(self, point3d: volmdlr.Point3D, **kwargs):
        """
        Ray Casting algorithm.

        Returns True if the point is inside the Shell, False otherwise
        """
        nb_rays = kwargs.get("nb_rays", 1)  # TODO: remove nb_rays argument in the future as it shouldn't be necessary

        bbox = self.bounding_box
        if not bbox.point_belongs(point3d):
            return False

        min_ray_length = 2 * max((bbox.xmax - bbox.xmin,
                                  bbox.ymax - bbox.ymin,
                                  bbox.zmax - bbox.zmin))
        two_min_ray_length = 2 * min_ray_length

        rays = []
        for _ in range(0, nb_rays):
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
            for _, point_inters in self.linesegment_intersections(ray):
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
        """
        Verifies if a given point belongs to some shell face.

        :param point: The point to check.
        :type point: volmdlr.Point3D
        :return: True if point belongs to some shell face. False otherwise.
        :rtype: bool
        """
        for face in self.faces:
            if (face.surface3d.point_on_surface(point) and face.point_belongs(point)) or \
                    face.outer_contour3d.point_over_contour(point, abs_tol=1e-7):
                return True
        return False

    def is_inside_shell(self, shell2, resolution: float):
        """
        Returns True if all the points of self are inside shell2 and no face are intersecting.

        This method is not exact.
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
        """
        Verifies and returns a Boolean if two shells are disjointed or not.

        """
        disjoint = True
        if self.bounding_box.bbox_intersection(shell2.bounding_box) or \
                self.bounding_box.distance_to_bbox(shell2.bounding_box) <= tol:
            return False
        return disjoint

    def intersecting_faces_combinations(self, shell2, list_coincident_faces, tol=1e-8):
        """
        Gets intersecting faces combinations.

        :param shell2: ClosedShell3D
            for two closed shells, it calculates and return a list of face
            combinations (list = [(face_shell1, face_shell2),...])
            for intersecting faces. if two faces can not be intersected,
            there is no combination for those
        :param tol: Corresponde to the tolerance to consider two faces as intersecting faces
        :param shell2:
        :param list_coincident_faces:
        :param tol:
        :return:
        """
        face_combinations = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                if face1.is_intersecting(face2, list_coincident_faces, tol):
                    face_combinations.append((face1, face2))
        return face_combinations

    @staticmethod
    def dict_intersecting_combinations(intersecting_faces_combinations, tol=1e-8):
        """
        Gets a Dictionary with the intersecting combinations.

        :param intersecting_faces_combinations: list of face combinations (list = [(face_shell1, face_shell2),...])
        for intersecting faces.
        :type intersecting_faces_combinations: list of face objects combinatons
        :param tol: tolerance
        returns a dictionary containing as keys the combination of intersecting faces
        and as the values the resulting primitive from the two intersecting faces.
        It is done so it is not needed to calculate the same intersecting primitive twice.
        """
        intersecting_combinations = {}
        for combination in intersecting_faces_combinations:
            face_intersections = combination[0].face_intersections(combination[1], tol)
            combination_face_intersections = []
            for face_intersection in face_intersections:
                for contour1 in [combination[0].outer_contour3d] + combination[0].inner_contours3d:
                    if contour1.is_superposing(face_intersection):
                        for contour2 in [combination[1].outer_contour3d] + combination[1].inner_contours3d:
                            if contour2.is_superposing(face_intersection):
                                break
                        else:
                            continue
                        break
                else:
                    combination_face_intersections.append(face_intersection)
            if combination_face_intersections:
                intersecting_combinations[combination] = combination_face_intersections
        return intersecting_combinations

    @staticmethod
    def get_intersecting_faces(dict_intersecting_combinations):
        """
        Gets intersecting faces.

        :param dict_intersecting_combinations: dictionary containing as keys the combination of intersecting faces
        and as the values the resulting primitive from the two intersecting faces

        returns two lists. One for the intersecting faces in shell1 and the other for the shell2
        """
        intersecting_faces_shell1 = []
        intersecting_faces_shell2 = []
        for face in list(dict_intersecting_combinations.keys()):
            if face[0] not in intersecting_faces_shell1:
                intersecting_faces_shell1.append(face[0])
            if face[1] not in intersecting_faces_shell2:
                intersecting_faces_shell2.append(face[1])
        return intersecting_faces_shell1, intersecting_faces_shell2

    def get_non_intersecting_faces(self, shell2, intersecting_faces, intersection_method=False):
        """
        Gets lists of faces that never intersect with any of the shell2's faces.

        :param shell2: ClosedShell3D.
        :param intersecting_faces:
        :param intersection_method: determines if running for intersection operation.
        returns a list of all the faces that never intersect any
        face of the other shell.
        """
        non_intersecting_faces = []

        for face in self.faces:
            if (face not in intersecting_faces) and (face not in non_intersecting_faces):
                if not intersection_method:
                    if not face.bounding_box.is_inside_bbox(shell2.bounding_box) or not shell2.is_face_inside(face):
                        for face2 in shell2.faces:
                            if face.surface3d.is_coincident(face2.surface3d) and \
                                    face.bounding_box.is_inside_bbox(face2.bounding_box):
                                break
                        else:
                            non_intersecting_faces.append(face)
                else:
                    if face.bounding_box.is_inside_bbox(shell2.bounding_box) and shell2.is_face_inside(face):
                        non_intersecting_faces.append(face)

        return non_intersecting_faces

    def get_coincident_and_adjacent_faces(self, shell2):
        coincident_and_adjacent_faces = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                if face1.surface3d.is_coincident(face2.surface3d) and \
                        face1.is_adjacent(face2):
                    coincident_and_adjacent_faces.append((face1, face2))

        return coincident_and_adjacent_faces

    def get_coincident_faces(self, shell2):
        """
        Finds all pairs of faces that are coincident faces, that is, faces lying on the same plane.

        returns a List of tuples with the face pairs.
        """
        list_coincident_faces = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                if isinstance(face1, face2.__class__) and face1.surface3d.is_coincident(face2.surface3d):
                    contour1 = face1.outer_contour3d.to_2d(
                        face1.surface3d.frame.origin,
                        face1.surface3d.frame.u,
                        face1.surface3d.frame.v)
                    contour2 = face2.outer_contour3d.to_2d(
                        face1.surface3d.frame.origin,
                        face1.surface3d.frame.u,
                        face1.surface3d.frame.v)
                    inters = contour1.contour_intersections(contour2)
                    if len(inters) >= 2:
                        list_coincident_faces.append((face1, face2))

        return list_coincident_faces

    def two_shells_intersecting_contour(self, shell2,
                                        list_coincident_faces: List[Face3D],
                                        dict_intersecting_combinations=None):
        """
        Computes intersecting_contour between two shells.

        :param shell2: ClosedShell3D
        :type shell2: :class:`volmdlr.faces.ClosedShell3D`
        :type list_coincident_faces: List[:class:`volmdlr.faces.Face3D`]
        :param dict_intersecting_combinations: dictionary containing as keys
            the combination of intersecting faces and as the values the
            resulting primitive from the two intersecting faces
        :returns: intersecting contour for two intersecting shells
        """
        if dict_intersecting_combinations is None:
            face_combinations = self.intersecting_faces_combinations(
                shell2, list_coincident_faces)
            dict_intersecting_combinations = \
                self.dict_intersecting_combinations(face_combinations)
        intersecting_wires = list(dict_intersecting_combinations.values())
        intersecting_contour = \
            volmdlr.wires.Contour3D([wire.primitives[0] for
                                     wires in intersecting_wires for wire in wires])
        return intersecting_contour

    def reference_shell(self, shell2, face):
        if face in shell2.faces:
            contour_extract_inside = True
            reference_shell = self
        else:
            contour_extract_inside = False
            reference_shell = shell2
        return contour_extract_inside, reference_shell

    def set_operations_valid_exterior_faces(self, new_faces: List[Face3D], valid_faces: List[Face3D],
                                            list_coincident_faces: List[Face3D], shell2, reference_shell):
        for new_face in new_faces:
            inside_reference_shell = reference_shell.point_belongs(new_face.random_point_inside())
            if self.set_operations_exterior_face(new_face, valid_faces, inside_reference_shell,
                                                 list_coincident_faces, shell2):
                valid_faces.append(new_face)
        return valid_faces

    def union_faces(self, shell2, intersecting_faces,
                    intersecting_combinations,
                    list_coincident_faces):
        faces = []
        for face in intersecting_faces:
            contour_extract_inside, reference_shell = self.reference_shell(shell2, face)
            new_faces = face.set_operations_new_faces(intersecting_combinations, contour_extract_inside)
            faces = self.set_operations_valid_exterior_faces(new_faces, faces, list_coincident_faces,
                                                             shell2, reference_shell)
        return faces

    def get_subtraction_valid_faces(self, new_faces, valid_faces, reference_shell, shell2, keep_interior_faces):
        faces = []
        for new_face in new_faces:
            inside_reference_shell = reference_shell.point_belongs(new_face.random_point_inside())
            if keep_interior_faces:
                if self.set_operations_interior_face(new_face, valid_faces, inside_reference_shell):
                    faces.append(new_face)
            elif self.set_operations_exterior_face(new_face, faces, inside_reference_shell, [], shell2):
                faces.append(new_face)
        return faces

    def validate_intersection_substractions_faces(self, faces):
        """
        Final validation of new faces created during intersections or subtractions of two closedshells.

        :param faces: new faces.
        :return: valid faces.
        """
        valid_faces = []
        finished = False
        while not finished:
            for face in valid_faces:
                if face.face_inside(faces[0]):
                    faces.remove(faces[0])
                    break
            else:
                valid_faces.append(faces[0])
                faces.remove(faces[0])
            if not faces:
                finished = True
        return valid_faces

    def subtraction_faces(self, shell2, intersecting_faces, intersecting_combinations):
        faces = []
        for face in intersecting_faces:
            keep_interior_faces = False
            if face in shell2.faces:
                keep_interior_faces = True
            contour_extract_inside, reference_shell = self.reference_shell(shell2, face)
            new_faces = face.set_operations_new_faces(intersecting_combinations, contour_extract_inside)
            valid_faces = self.get_subtraction_valid_faces(new_faces, faces, reference_shell,
                                                           shell2, keep_interior_faces)
            faces.extend(valid_faces)

        valid_faces = self.validate_intersection_substractions_faces(faces)

        return valid_faces

    def valid_intersection_faces(self, new_faces, valid_faces,
                                 reference_shell, shell2):
        faces = []
        for new_face in new_faces:
            inside_reference_shell = reference_shell.point_belongs(
                new_face.random_point_inside())
            if (inside_reference_shell or (self.face_on_shell(new_face) and shell2.face_on_shell(new_face))) \
                    and new_face not in valid_faces:
                faces.append(new_face)

        return faces

    def intersection_faces(self, shell2, intersecting_faces,
                           intersecting_combinations):
        faces = []
        for face in intersecting_faces:
            contour_extract_inside, reference_shell = \
                self.reference_shell(shell2, face)
            new_faces = face.set_operations_new_faces(
                intersecting_combinations, contour_extract_inside)
            valid_faces = self.valid_intersection_faces(
                new_faces, faces, reference_shell, shell2)
            faces.extend(valid_faces)

        valid_faces = self.validate_intersection_substractions_faces(faces)
        return valid_faces

    def set_operations_interior_face(self, new_face, faces, inside_reference_shell):
        if inside_reference_shell and new_face not in faces:
            return True
        if self.face_on_shell(new_face):
            return True
        return False

    def is_face_between_shells(self, shell2, face):
        if face.surface2d.inner_contours:
            normal_0 = face.surface2d.outer_contour.primitives[0].normal_vector()
            middle_point_0 = face.surface2d.outer_contour.primitives[0].middle_point()
            point1 = middle_point_0 + 0.0001 * normal_0
            point2 = middle_point_0 - 0.0001 * normal_0
            points = [point1, point2]
        else:
            points = [face.surface2d.outer_contour.center_of_mass()]

        for point in points:
            point3d = face.surface3d.point2d_to_3d(point)
            if face.point_belongs(point3d):
                normal1 = point3d - 0.00001 * face.surface3d.frame.w
                normal2 = point3d + 0.00001 * face.surface3d.frame.w
                if (self.point_belongs(normal1) and
                    shell2.point_belongs(normal2)) or \
                        (shell2.point_belongs(normal1) and
                         self.point_belongs(normal2)):
                    return True
        return False

    def set_operations_exterior_face(self, new_face, valid_faces, inside_reference_shell,
                                     list_coincident_faces, shell2):
        if new_face.area() < 1e-8:
            return False
        if new_face not in valid_faces and not inside_reference_shell:
            if list_coincident_faces:
                if self.is_face_between_shells(shell2, new_face):
                    return False
            return True
        return False

    def validate_set_operation(self, shell2, tol):
        """
        Verifies if two shells are valid for union or subtractions operations.

        Its Verfies if they are disjointed or if one is totally inside the other.

        If it returns an empty list, it means the two shells are valid to continue the
        operation.
        """
        if self.is_disjoint_from(shell2, tol):
            return [self, shell2]
        if self.is_inside_shell(shell2, resolution=0.01):
            return [shell2]
        if shell2.is_inside_shell(self, resolution=0.01):
            return [self]
        return []

    def is_clean(self):
        """
        Verifies if closed shell\'s faces are clean or if it is needed to be cleaned.

        :return: True if clean and False Otherwise
        """
        for face1, face2 in product(self.faces, repeat=2):
            if face1 != face2 and \
                    face1.surface3d.is_coincident(face2.surface3d) and \
                    face1.is_adjacent(face2):
                return False
        return True

    def union(self, shell2: 'ClosedShell3D', tol: float = 1e-8):
        """
        Given Two closed shells, it returns a new united ClosedShell3D object.

        """

        validate_set_operation = \
            self.validate_set_operation(shell2, tol)
        if validate_set_operation:
            return validate_set_operation
        list_coincident_faces = self.get_coincident_faces(shell2)
        face_combinations = self.intersecting_faces_combinations(shell2, list_coincident_faces, tol)
        intersecting_combinations = self.dict_intersecting_combinations(face_combinations, tol)
        intersecting_faces1, intersecting_faces2 = self.get_intersecting_faces(intersecting_combinations)
        intersecting_faces = intersecting_faces1 + intersecting_faces2
        faces = self.get_non_intersecting_faces(shell2, intersecting_faces) + \
                shell2.get_non_intersecting_faces(self, intersecting_faces)
        if len(faces) == len(self.faces + shell2.faces) and not intersecting_faces:
            return [self, shell2]
        new_valid_faces = self.union_faces(shell2, intersecting_faces,
                                           intersecting_combinations, list_coincident_faces)
        faces += new_valid_faces
        new_shell = ClosedShell3D(faces)
        return [new_shell]

    @staticmethod
    def get_faces_to_be_merged(union_faces):
        coincident_planes_faces = []
        for i, face1 in enumerate(union_faces):
            for j, face2 in enumerate(union_faces):
                if j != i and face1.surface3d.is_coincident(face2.surface3d):
                    if face1 not in coincident_planes_faces:
                        coincident_planes_faces.append(face1)
                    coincident_planes_faces.append(face2)
            if coincident_planes_faces:
                break
        return coincident_planes_faces

    @staticmethod
    def clean_faces(union_faces, list_new_faces):
        list_remove_faces = []
        if union_faces:
            for face1 in union_faces:
                for face2 in list_new_faces:
                    if face1.face_inside(face2):
                        list_remove_faces.append(face2)
                    elif face2.face_inside(face1):
                        list_remove_faces.append(face1)
        list_new_faces += union_faces
        for face in list_remove_faces:
            list_new_faces.remove(face)
        return list_new_faces

    def merge_faces(self):
        """
        Merges all shells' adjacent faces into one.

        """
        union_faces = self.faces
        finished = False
        list_new_faces = []
        count = 0
        while not finished:
            valid_coicident_faces = ClosedShell3D.get_faces_to_be_merged(union_faces)
            list_valid_coincident_faces = valid_coicident_faces[:]
            if valid_coicident_faces:
                list_new_faces += PlaneFace3D.merge_faces(valid_coicident_faces)
            for face in list_valid_coincident_faces:
                union_faces.remove(face)
            count += 1
            if count >= len(self.faces) and not list_valid_coincident_faces:
                finished = True

        list_new_faces = self.clean_faces(union_faces, list_new_faces)

        self.faces = list_new_faces

    def subtract(self, shell2, tol=1e-8):
        """
        Given Two closed shells, it returns a new subtracted OpenShell3D.

        """
        validate_set_operation = self.validate_set_operation(shell2, tol)
        if validate_set_operation:
            return validate_set_operation

        list_coincident_faces = self.get_coincident_faces(shell2)
        face_combinations = self.intersecting_faces_combinations(
            shell2, list_coincident_faces, tol)

        intersecting_combinations = self.dict_intersecting_combinations(
            face_combinations, tol)

        if len(intersecting_combinations) == 0:
            return [self, shell2]

        intersecting_faces, _ = self.get_intersecting_faces(
            intersecting_combinations)

        faces = self.get_non_intersecting_faces(shell2, intersecting_faces)
        new_valid_faces = self.union_faces(shell2, intersecting_faces,
                                           intersecting_combinations,
                                           list_coincident_faces
                                           )
        faces += new_valid_faces
        return [OpenShell3D(faces)]

    def subtract_to_closed_shell(self, shell2: OpenShell3D, tol: float = 1e-8):
        """
        Given Two closed shells, it returns a new subtracted ClosedShell3D.

        :param shell2:
        :param tol:
        :return:
        """

        validate_set_operation = self.validate_set_operation(shell2, tol)
        if validate_set_operation:
            return validate_set_operation

        list_coincident_faces = self.get_coincident_faces(shell2)
        face_combinations = self.intersecting_faces_combinations(
            shell2, list_coincident_faces, tol)
        intersecting_combinations = self.dict_intersecting_combinations(
            face_combinations, tol)

        if len(intersecting_combinations) == 0:
            return [self, shell2]

        intersecting_faces1, intersecting_faces2 = self.get_intersecting_faces(
            intersecting_combinations)
        intersecting_faces = intersecting_faces1 + intersecting_faces2

        faces = self.get_non_intersecting_faces(shell2, intersecting_faces)
        faces += shell2.get_non_intersecting_faces(self, intersecting_faces, intersection_method=True)
        new_valid_faces = self.subtraction_faces(shell2, intersecting_faces, intersecting_combinations)
        faces += new_valid_faces
        new_shell = ClosedShell3D(faces)
        # new_shell.eliminate_not_valid_closedshell_faces()
        return [new_shell]

    def intersection(self, shell2, tol=1e-8):
        """
        Given two ClosedShell3D, it returns the new object resulting from the intersection of the two.

        """
        validate_set_operation = self.validate_set_operation(
            shell2, tol)
        if validate_set_operation:
            return validate_set_operation
        list_coincident_faces = self.get_coincident_faces(shell2)
        face_combinations = self.intersecting_faces_combinations(shell2, list_coincident_faces, tol)
        intersecting_combinations = self.dict_intersecting_combinations(face_combinations, tol)

        if not intersecting_combinations:
            return [self, shell2]

        intersecting_faces1, intersecting_faces2 = self.get_intersecting_faces(intersecting_combinations)
        intersecting_faces = intersecting_faces1 + intersecting_faces2
        faces = self.intersection_faces(shell2, intersecting_faces, intersecting_combinations)
        faces += self.get_non_intersecting_faces(shell2, intersecting_faces, intersection_method=True) + \
                 shell2.get_non_intersecting_faces(self, intersecting_faces, intersection_method=True)
        new_shell = ClosedShell3D(faces)
        new_shell.eliminate_not_valid_closedshell_faces()
        return [new_shell]

    def eliminate_not_valid_closedshell_faces(self):
        nodes_with_2degrees = [node for node, degree in list(self.faces_graph.degree()) if degree <= 2]
        for node in nodes_with_2degrees:
            neighbors = nx.neighbors(self.faces_graph, node)
            for neighbor_node in neighbors:
                for face in self.faces:
                    if self.faces_graph.edges[(node, neighbor_node)]['edge'] in face.outer_contour3d.primitives:
                        self.faces.remove(face)
                        break
        self._faces_graph = None


class OpenTriangleShell3D(OpenShell3D):
    """
    A 3D open shell composed of multiple triangle faces.

    This class represents a 3D open shell, which is a collection of connected
    triangle faces with no volume. It is a subclass of the `OpenShell3D` class
    and inherits all of its attributes and methods.

    :param faces: The triangle faces of the shell.
    :type faces: List[`Triangle3D`]
    :param color: The color of the shell.
    :type color: Tuple[float, float, float]
    :param alpha: The transparency of the shell, should be a value in the range (0, 1).
    :type alpha: float
    :param name: The name of the shell.
    :type name: str
    """

    def __init__(self, faces: List[Triangle3D],
                 color: Tuple[float, float, float] = None,
                 alpha: float = 1., name: str = ''):
        OpenShell3D.__init__(self, faces=faces, color=color, alpha=alpha, name=name)

    def to_dict(self):
        dict_ = self.base_dict()
        dict_['faces'] = [t.to_dict() for t in self.faces]
        dict_['alpha'] = self.alpha
        dict_['color'] = self.color
        return dict_

    def point_on_shell(self, point: volmdlr.Point3D):
        for face in self.faces:
            # TODO: why the first check?
            if (face.surface3d.point_on_plane(point) and face.point_belongs(point)) or \
                    face.outer_contour3d.point_over_contour(point, abs_tol=1e-7):
                return True
        return False

    def to_mesh_data(self):
        positions = npy.zeros((3 * len(self.faces), 3))
        faces = npy.zeros((len(self.faces), 3))
        for i, triangle_face in enumerate(self.faces):
            i1 = 3 * i
            i2 = i1 + 1
            i3 = i1 + 2
            positions[i1, 0] = triangle_face.points[0].x
            positions[i1, 1] = triangle_face.points[0].y
            positions[i1, 2] = triangle_face.points[0].z
            positions[i2, 0] = triangle_face.points[1].x
            positions[i2, 1] = triangle_face.points[1].y
            positions[i2, 2] = triangle_face.points[1].z
            positions[i3, 0] = triangle_face.points[2].x
            positions[i3, 1] = triangle_face.points[2].y
            positions[i3, 2] = triangle_face.points[2].z

            faces[i, 0] = i1
            faces[i, 1] = i2
            faces[i, 2] = i3

        return positions, faces

    @classmethod
    def from_mesh_data(cls, positions, faces):
        triangles = []
        points = [volmdlr.Point3D(px, py, pz) for px, py, pz in positions]
        for i1, i2, i3 in faces:
            triangles.append(Triangle3D(points[i1], points[i2], points[i3]))
        return cls(triangles)

    def to_trimesh(self):
        return Trimesh(*self.to_mesh_data())

    @classmethod
    def from_trimesh(cls, trimesh):
        return cls.from_mesh_data(trimesh.vertices.tolist(), trimesh.faces.tolist())

    def triangulation(self):
        points = []
        triangles = []
        for i, triangle in enumerate(self.faces):
            points.append(vmd.Node3D.from_point(triangle.point1))
            points.append(vmd.Node3D.from_point(triangle.point2))
            points.append(vmd.Node3D.from_point(triangle.point3))
            triangles.append((3 * i, 3 * i + 1, 3 * i + 2))
        return vmd.DisplayMesh3D(points, triangles)


class ClosedTriangleShell3D(ClosedShell3D, OpenTriangleShell3D):
    """
        A 3D closed shell composed of multiple triangle faces.

    This class represents a 3D closed shell, which is a collection of connected
    triangle faces with a volume. It is a subclass of both the `ClosedShell3D`
    and `OpenTriangleShell3D` classes and inherits all of their attributes and
    methods.

    :param faces: The triangle faces of the shell.
    :type faces: List[`Triangle3D`]
    :param color: The color of the shell.
    :type color: Tuple[float, float, float]
    :param alpha: The transparency of the shell, should be a value in the range (0, 1).
    :type alpha: float
    :param name: The name of the shell.
    :type name: str
    """

    def __init__(self, faces: List[Triangle3D],
                 color: Tuple[float, float, float] = None,
                 alpha: float = 1., name: str = ''):
        OpenTriangleShell3D.__init__(self, faces=faces, color=color, alpha=alpha, name=name)
