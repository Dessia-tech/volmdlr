# -*- coding: utf-8 -*-
"""
Cloud of points classes.

"""

import math
from typing import List, Tuple
import random

import dessia_common.core as dc
import matplotlib.pyplot as plt
from trimesh.proximity import closest_point

import volmdlr as vm
import volmdlr.faces as vmf
from volmdlr import shells as vmshells
from volmdlr import surfaces
import volmdlr.primitives3d as p3d
import volmdlr.step as vstep
import volmdlr.stl as vmstl
# import volmdlr.core
import volmdlr.wires as vmw


class PointCloud3D(dc.DessiaObject):
    """
    Point Cloud3D, a list of points.

    :param points: a list of points.
    """

    def __init__(self, points: List[vm.Point3D], name: str = ''):
        self.points = points
        self.__bounding_box = None
        dc.DessiaObject.__init__(self, name=name)

    @classmethod
    def from_stl(cls, file_path, name: str = 'from_stl'):
        """
        Creates a point cloud 3d from a stl file.

        :param file_path: path to stl file.
        :param name: object's name.
        :return: point cloud 3d object.
        """
        list_points = vmstl.Stl.from_file(file_path).extract_points_bis()

        return cls(list_points, name=name)

    def _bounding_box(self):
        """
        Computes the bounding box of the point cloud.

        :return: Bounding box object.
        """
        return vm.core.BoundingBox.from_points(self.points)

    @property
    def bounding_box(self):
        """
        Property to get the bounding box of the point cloud.

        :return: Bounding box object.
        """
        if not self.__bounding_box:
            self.__bounding_box = self._bounding_box()
        return self.__bounding_box

    def to_2d(self, plane_origin, x, y):
        """
        Projects the point cloud on a 2D plane.

        :param plane_origin: Origin of the plane in 3D.
        :param x: X-axis vector in 3D.
        :param y: Y-axis vector in 3D.
        :return: PointCloud2D object.
        """
        list_points2d = [pt3d.to_2d(plane_origin, x, y) for pt3d in self.points]
        return PointCloud2D(list_points2d, name='3d_to_2d')

    def extract(self, u, umin, umax):  # -> List[PointCloud3D] :
        """
        Extracts a subset of points within a given range from a plane.

        :param u: Normal vector of the plane.
        :param umin: Minimum distance to the plane.
        :param umax: Maximum distance to the plane.
        :return: PointCloud3D object containing the extracted points.
        """
        extracted_points = []
        for points in self.points:
            dist_to_plane = points.dot(u)
            if umin < dist_to_plane < umax:
                extracted_points.append(points)
        return PointCloud3D(extracted_points)

    def determine_extrusion_vector(self):
        """
        Determines the extrusion vector based on the bounding box.

        :return: Tuple containing position, normal, vec1, and vec2.
        """
        bbox = self._bounding_box()
        xyz_bbox = [[bbox.xmin, bbox.xmax], [bbox.ymin, bbox.ymax], [bbox.zmin, bbox.zmax]]
        xyz_list = [l[1] - l[0] for l in xyz_bbox]
        absxyz_list, xyz_vect = [abs(length) for length in xyz_list], [vm.X3D, vm.Y3D, vm.Z3D]
        posmax = xyz_list.index(max(absxyz_list))
        normal = xyz_vect[posmax]
        vec1, vec2 = xyz_vect[posmax - 2], xyz_vect[posmax - 1]

        return posmax, normal, vec1, vec2

    def position_plane(self, posmax, resolution):
        """
        Calculates the position of planes.

        :param posmax: Position index.
        :param resolution: Resolution of the plane.
        :return: Tuple containing distance between planes and position of planes.
        """
        bbox = self._bounding_box()
        xyz_bbox = [[bbox.xmin, bbox.xmax], [bbox.ymin, bbox.ymax], [bbox.zmin, bbox.zmax]]
        dist_between_plane = (xyz_bbox[posmax][1] - xyz_bbox[posmax][0]) / (resolution - 1)
        position_plane = [xyz_bbox[posmax][0] + n * dist_between_plane for n in range(resolution)]

        return dist_between_plane, position_plane

    @staticmethod
    def check_area_polygon(initial_polygons2d, position_plane,
                           normal, vec1, vec2):
        """
        Checks and processes area of polygons.

        :param initial_polygons2d: List of 2D polygons.
        :param position_plane: Position of planes.
        :param normal: Normal vector.
        :param vec1: Vector 1.
        :param vec2: Vector 2.
        :return: List of processed 3D polygons.
        """
        areas = [0] * len(initial_polygons2d)
        for n, poly in enumerate(initial_polygons2d):
            if poly is not None:
                areas[n] = poly.area()
        avg_area = sum(areas) / len(areas)

        polygon2d, polygons_3d = [], []
        for n, poly in enumerate(initial_polygons2d):
            if poly is None or (poly.area() < avg_area / 10) and (n not in [0, len(initial_polygons2d) - 1]):
                continue
            if poly.area() < avg_area / 10:
                new_poly = vmw.ClosedPolygon2D.concave_hull(poly.points, -1, 0.000005)
                new_polygon = new_poly.to_3d(position_plane[n] * normal, vec1,
                                             vec2)
                polygons_3d.append(new_polygon)
            else:

                polygon2d.append(poly)
                new_polygon = poly.to_3d(position_plane[n] * normal, vec1, vec2)
                polygons_3d.append(new_polygon)

        return polygons_3d

    def to_subcloud2d(self, pos_normal, vec1, vec2):
        """
        Converts the point cloud to a simplified 2D sub-cloud.

        :param pos_normal: Position and normal vector.
        :param vec1: Vector 1.
        :param vec2: Vector 2.
        :return: Simplified PointCloud2D object.
        """
        subcloud2d_tosimp = self.to_2d(pos_normal, vec1, vec2)
        subcloud2d = subcloud2d_tosimp.simplify(resolution=5)
        return subcloud2d

    def to_shell(self, resolution: int = 10, normal=None, offset: float = 0):
        """Creates a Shell from a Cloud of points 3D."""
        if normal is None:
            posmax, normal, vec1, vec2 = self.determine_extrusion_vector()
        else:
            posmax = 0
            for n, vect in enumerate([vm.X3D, vm.Y3D, vm.Z3D]):
                if vect == normal:
                    posmax = n
            vec1, vec2 = [vm.X3D, vm.Y3D, vm.Z3D][posmax - 2], [vm.X3D, vm.Y3D, vm.Z3D][posmax - 1]

        dist_between_plane, position_plane = self.position_plane(posmax=posmax,
                                                                 resolution=resolution)
        sub_clouds3d = [self.extract(normal, pos_plane - 0.5 * dist_between_plane,
                                     pos_plane + 0.5 * dist_between_plane) for pos_plane in position_plane]
        sub_clouds2d = [sub_clouds3d[n].to_subcloud2d((position_plane[n] * normal).to_point(), vec1, vec2)
                        for n in range(resolution)]

        # Offsetting
        if offset != 0:
            initial_polygons2d = [cloud2d.to_polygon(convex=True) for cloud2d in sub_clouds2d]
            position_plane, initial_polygons2d = self.offset_to_shell(position_plane, initial_polygons2d, offset)
        else:
            initial_polygons2d = [cloud2d.to_polygon() for cloud2d in sub_clouds2d]
        polygons_3d = self.check_area_polygon(initial_polygons2d=initial_polygons2d,
                                              position_plane=position_plane,
                                              normal=normal,
                                              vec1=vec1, vec2=vec2)

        return self.generate_shell(polygons_3d, normal, vec1, vec2)

    @classmethod
    def generate_shell(cls, polygons_3d: List[vm.wires.ClosedPolygon3D],
                       normal: vm.Vector3D, vec1: vm.Vector3D, vec2: vm.Vector3D, name: str = ''):
        """
        Generates a shell from a list of polygon 3d, using a sewing algorithm.

        :param polygons_3d: list of polygon 3d to be sewed.
        :param normal: normal to the sewing plane.
        :param vec1: u vector in the sewing plane.
        :param vec2: v vector in the sewing plane.
        :param name: object's name.
        :return: return a shell.
        """
        position_plane = [p.points[0].dot(normal) for p in polygons_3d]
        resolution = len(polygons_3d)
        faces = []
        for n in range(resolution):
            poly1 = polygons_3d[n]
            poly1_simplified = cls._helper_simplify_polygon(poly1, position_plane[n], normal, vec1, vec2)

            if n in (resolution - 1, 0):
                vec2_ = vec2
                if n == 0:
                    vec2_ = -vec2
                faces.append(
                    vmf.PlaneFace3D(
                        surface3d=surfaces.Plane3D.from_plane_vectors((position_plane[n] * normal).to_point(),
                                                                      vec1, vec2_),
                        surface2d=cls._poly_to_surf2d(poly1_simplified, position_plane[n], normal, vec1, vec2_)))

            if n != resolution - 1:
                poly2 = polygons_3d[n + 1]
                poly2_simplified = cls._helper_simplify_polygon(poly2, position_plane[n + 1], normal, vec1, vec2)

                list_triangles_points = cls._helper_sew_polygons(poly1_simplified, poly2_simplified, vec1, vec2)
                list_faces = [vmf.Triangle3D(*triangle_points, alpha=0.9, color=(1, 0.1, 0.1))
                              for triangle_points in list_triangles_points]
                faces.extend(list_faces)
        return vmshells.ClosedShell3D(faces, name=name)

    @staticmethod
    def _helper_simplify_polygon(polygon, position_plane, normal, vec1, vec2):
        poly_2d = polygon.to_2d((position_plane * normal).to_point(), vec1, vec2)
        poly_2d_simplified = poly_2d.simplify_polygon(0.01, 1)
        if 1 - poly_2d_simplified.area() / poly_2d.area() > 0.3:
            poly_2d_simplified = poly_2d
        return poly_2d_simplified.to_3d(position_plane * normal, vec1, vec2)

    @staticmethod
    def _poly_to_surf2d(polygon, position_plane, normal, vec1, vec2):
        return surfaces.Surface2D(polygon.to_2d((position_plane * normal).to_point(), vec1, vec2), [])

    @staticmethod
    def _helper_sew_polygons(poly1, poly2, vec1, vec2):
        return poly1.sewing(poly2, vec1, vec2)

    def shell_distances(self, shells: vmshells.OpenTriangleShell3D) -> Tuple['PointCloud3D', List[float], List[int]]:
        """
        Computes distance of point to shell for each point in self.points.

        :return: The point cloud of points projection on nearest triangle, their distances and the corresponding
        triangles index
        :rtype: Tuple[PointCloud3D, List[float], List[int]]
        """
        nearest_coords, distances, triangles_idx = self.shell_distances_ndarray(shells)
        return (PointCloud3D([vm.Point3D(*coords) for coords in nearest_coords]),
                distances.tolist(),
                triangles_idx.tolist())

    def shell_distances_ndarray(self, shells: vmshells.OpenTriangleShell3D):
        """
        Computes distance of point to shell for each point in self.points in a numpy formatted data.

        :return: The point cloud of points projection on nearest triangle, their distances and the corresponding
        triangles index
        :rtype: Tuple[numpy.ndarray(float), numpy.ndarray(float), numpy.ndarray(int)]
        """
        shells_trimesh = shells.to_trimesh()
        return closest_point(shells_trimesh, self.to_coord_matrix())

    def to_coord_matrix(self) -> List[List[float]]:
        """Generate an n_points x 3 matrix of coordinates."""
        return [point.coordinates() for point in self.points]

    @classmethod
    def from_step(cls, step_file: str, name: str = ''):
        """
        Creates a cloud of points from a step file.

        :param step_file: step file.
        :param name: object's name.
        :return: Point Cloud 3D.
        """
        step = vstep.Step(step_file)
        points = step.to_points()
        return cls(points, name=name)

    def plot(self, ax=None, color='k'):
        """
        Plot the cloud 3d.

        """

        if self.points:
            ax = self.points[0].plot(ax=ax)
            if len(self.points) > 1000:
                plot_points = random.choices(self.points, k=1000)
            else:
                plot_points = self.points
            for point in plot_points:
                point.plot(ax=ax, color=color)

        return ax

    def extended_cloud(self, distance_extended: float):
        # it works if distance_extended >= 0
        spheres, extended_points = [], []
        for point in self.points:
            extended_zone = p3d.Sphere(point, distance_extended)

            spheres.append(extended_zone)

            extended_points.extend(extended_zone.triangulation().points)

        for sphere in spheres:
            clean_extended_zone = []
            for point in extended_points:
                if not sphere.point_belongs(point):
                    clean_extended_zone.append(point)
            extended_points = clean_extended_zone

        return extended_points

    @staticmethod
    def offset_to_shell(positions_plane: List[surfaces.Plane3D],
                        polygons2d: List[vmw.ClosedPolygon2D], offset: float):

        origin_f, origin_l = positions_plane[0], positions_plane[-1]

        new_position_plane = [origin_f - offset] + positions_plane[1:-1] + [origin_l + offset]
        polyconvexe = [vmw.ClosedPolygon2D.points_convex_hull(poly.points) for poly in polygons2d]
        new_poly = [poly.offset(offset) for poly in polyconvexe]

        return new_position_plane, new_poly


class PointCloud2D(dc.DessiaObject):
    """
    Point Cloud2D class.

    :param points: list of points for point cloud.
    """

    def __init__(self, points, name: str = ''):
        self.points = points
        dc.DessiaObject.__init__(self, name=name)

    def plot(self, ax=None, color='k'):
        """Plot a point cloud 2d using Matplotlib."""
        if ax is None:
            _, ax = plt.subplots()
        for point in self.points:
            point.plot(ax=ax, color=color)
        return ax

    def to_polygon(self, convex=False):
        """
        Use a Cloud point 2d to create a polygon.

        :param convex: if True, it will return a convex polygon. If false, it will search for a concave polygon.
        :return: closed polygon 2d.
        """
        if not self.points:
            return None

        # polygon = vmw.ClosedPolygon2D.convex_hull_points(self.points)
        if convex:
            polygon = vmw.ClosedPolygon2D.convex_hull_points(self.points)
        else:
            polygon = vmw.ClosedPolygon2D.concave_hull(self.points, -0.2, 0.000005)
        polygon = polygon.simplify(min_distance=0.002)
        if polygon is None or math.isclose(polygon.area(), 0, abs_tol=1e-6):
            return None
        return polygon

    def bounding_rectangle(self):
        """
        Calculates the bounding rectangle for the point cloud.

        :return: bounds for bounding rectangle.
        """
        x_list, y_list = [p.x for p in self.points], [p.y for p in self.points]
        return min(x_list), max(x_list), min(y_list), max(y_list)

    def simplify(self, resolution=5):
        """Simplify cloud point 2d."""
        if not self.points:
            return PointCloud2D(self.points, name=self.name + '_none')

        xy_extr = list(self.bounding_rectangle())  # xmin, xmax, ymin, ymax

        x_slide = [xy_extr[0] + n * (xy_extr[1] - xy_extr[0]) / (resolution - 1) for n in range(resolution)]
        y_slide = [xy_extr[2] + n * (xy_extr[3] - xy_extr[2]) / (resolution - 1) for n in range(resolution)]

        points = []
        for x1, x2 in zip(x_slide, x_slide[1:] + [x_slide[0]]):
            for y1, y2 in zip(y_slide, y_slide[1:] + [y_slide[0]]):
                box_points = []
                for point in self.points:
                    if point.x >= x1 and point.x <= x2:
                        if point.y >= y1 and point.y <= y2:
                            box_points.append(point)
                points.append(box_points)

        polys = [vmw.ClosedPolygon2D.points_convex_hull(pts) for pts in points]
        clean_points = []
        for poly in polys:
            if poly is not None:
                clean_points += poly.points
        return PointCloud2D(clean_points, name=self.name + '_clean')

    def to_coord_matrix(self) -> List[List[float]]:
        """
        Generate an n_points x 2 matrix of coordinates.
        """
        return [point.coordinates() for point in self.points]
