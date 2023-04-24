# -*- coding: utf-8 -*-
"""
Cloud of points classes.

"""

import math
from typing import List, Tuple

import dessia_common.core as dc
import matplotlib.pyplot as plt
from trimesh.proximity import closest_point

import volmdlr as vm
import volmdlr.faces as vmf
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
        dc.DessiaObject.__init__(self, name=name)

    @classmethod
    def from_stl(cls, file_path):
        """
        Creates a point cloud 3d from an stl file.

        :param file_path: path to stl file.
        :return: point cloud 3d object.
        """
        list_points = vmstl.Stl.from_file(file_path).extract_points_BIS()

        return cls(list_points, name='from_stl')

    def _bounding_box(self):
        return vm.core.BoundingBox.from_points(self.points)

    def to_2d(self, plane_origin, x, y):
        list_points2d = [pt3d.to_2d(plane_origin, x, y) for pt3d in self.points]
        return PointCloud2D(list_points2d, name='3d_to_2d')

    def extract(self, u, umin, umax):  # -> List[PointCloud3D] :
        extracted_points = []
        for points in self.points:
            dist_to_plane = points.dot(u)
            if dist_to_plane > umin and dist_to_plane < umax:
                extracted_points.append(points)
        return PointCloud3D(extracted_points)

    def determine_extrusion_vector(self):
        bbox = self._bounding_box()
        xyz_bbox = [[bbox.xmin, bbox.xmax], [bbox.ymin, bbox.ymax], [bbox.zmin, bbox.zmax]]
        xyz_list = [l[1] - l[0] for l in xyz_bbox]
        absxyz_list, xyz_vect = [abs(length) for length in xyz_list], [vm.X3D, vm.Y3D, vm.Z3D]
        posmax = xyz_list.index(max(absxyz_list))
        normal = xyz_vect[posmax]
        vec1, vec2 = xyz_vect[posmax - 2], xyz_vect[posmax - 1]

        return posmax, normal, vec1, vec2

    def position_plane(self, posmax, resolution):
        bbox = self._bounding_box()
        xyz_bbox = [[bbox.xmin, bbox.xmax], [bbox.ymin, bbox.ymax], [bbox.zmin, bbox.zmax]]
        dist_between_plane = (xyz_bbox[posmax][1] - xyz_bbox[posmax][0]) / (resolution - 1)
        position_plane = [xyz_bbox[posmax][0] + n * dist_between_plane for n in range(resolution)]

        return dist_between_plane, position_plane

    @staticmethod
    def check_area_polygon(initial_polygon2d, position_plane,
                           normal, vec1, vec2):
        areas = [0] * len(initial_polygon2d)
        for n, poly in enumerate(initial_polygon2d):
            if poly is not None:
                areas[n] = poly.area()
        avg_area = sum(areas) / len(areas)

        polygon2d, polygon3d = [], []
        for n, poly in enumerate(initial_polygon2d):
            if poly is None or (poly.area() < avg_area / 10) and (n not in [0, len(initial_polygon2d) - 1]):
                continue
            if poly.area() < avg_area / 10:
                new_poly = vmw.ClosedPolygon2D.concave_hull(poly.points, -1, 0.000005)
                new_polygon = new_poly.to_3d(position_plane[n] * normal, vec1,
                                             vec2)
                polygon3d.append(new_polygon)
            else:

                polygon2d.append(poly)
                new_polygon = poly.to_3d(position_plane[n] * normal, vec1, vec2)
                polygon3d.append(new_polygon)

        return polygon3d

    def to_subcloud2d(self, pos_normal, vec1, vec2):
        subcloud2d_tosimp = self.to_2d(pos_normal, vec1, vec2)
        subcloud2d = subcloud2d_tosimp.simplify(resolution=5)
        return subcloud2d

    def to_shell(self, resolution: int = 10, normal=None, offset: float = 0):
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
        subcloud3d = [
            self.extract(
                normal,
                pos_plane -
                dist_between_plane /
                2,
                pos_plane +
                dist_between_plane /
                2) for pos_plane in position_plane]
        subcloud2d = [subcloud3d[n].to_subcloud2d(position_plane[n] * normal, vec1, vec2) for n in range(resolution)]

        # Offsetting
        if offset != 0:
            initial_polygon2d = [cloud2d.to_polygon(convex=True) for cloud2d in subcloud2d]
            position_plane, initial_polygon2d = self.offset_to_shell(position_plane, initial_polygon2d, offset)
        else:
            initial_polygon2d = [cloud2d.to_polygon() for cloud2d in subcloud2d]
        polygon3d = self.check_area_polygon(initial_polygon2d=initial_polygon2d,
                                            position_plane=position_plane,
                                            normal=normal,
                                            vec1=vec1, vec2=vec2)

        return self.generate_shell(polygon3d, normal, vec1, vec2)

    @classmethod
    def generate_shell(cls, polygon3d: List[vm.wires.ClosedPolygon3D],
                       normal: vm.Vector3D, vec1: vm.Vector3D, vec2: vm.Vector3D):
        """
        Generates a shell from a list of polygon 3d, using a sewing algorithm.

        :param polygon3d: list of polygon 3d to be sewed.
        :param normal: normal to the sewing plane.
        :param vec1: u vector in the sewing plane.
        :param vec2: v vector in the sewing plane.
        :return: return a shell.
        """
        position_plane = [p.points[0].dot(normal) for p in polygon3d]
        resolution = len(polygon3d)

        faces = []
        for n in range(resolution):
            poly1 = polygon3d[n]
            poly1_simplified = cls._helper_simplify_polygon(poly1, position_plane[n], normal, vec1, vec2)

            if n in (resolution - 1, 0):
                faces.append(
                    vmf.PlaneFace3D(surface3d=vmf.Plane3D.from_plane_vectors(position_plane[n] * normal, vec1, vec2),
                                    surface2d=cls._poly_to_surf2d(poly1_simplified, position_plane[n],
                                                                  normal, vec1, vec2)))

            if n != resolution - 1:
                poly2 = polygon3d[n + 1]
                poly2_simplified = cls._helper_simplify_polygon(poly2, position_plane[n + 1], normal, vec1, vec2)

                list_triangles_points = cls._helper_sew_polygons(poly1_simplified, poly2_simplified, vec1, vec2)
                list_faces = [vmf.Triangle3D(*triangle_points, alpha=0.9, color=(1, 0.1, 0.1))
                              for triangle_points in list_triangles_points]
                faces.extend(list_faces)
        return vmf.ClosedShell3D(faces)

    @staticmethod
    def _helper_simplify_polygon(polygon, position_plane, normal, vec1, vec2):
        poly_2d = polygon.to_2d(position_plane * normal, vec1, vec2)
        poly_2d_simplified = poly_2d.simplify_polygon(0.01, 1)
        if 1 - poly_2d_simplified.area() / poly_2d.area() > 0.3:
            poly_2d_simplified = poly_2d
        return poly_2d_simplified.to_3d(position_plane * normal, vec1, vec2)

    @staticmethod
    def _poly_to_surf2d(polygon, position_plane, normal, vec1, vec2):
        return vmf.Surface2D(polygon.to_2d(position_plane * normal, vec1, vec2), [])

    @staticmethod
    def _helper_sew_polygons(poly1, poly2, vec1, vec2):
        return poly1.sewing(poly2, vec1, vec2)

    def shell_distances(self, shells: vmf.OpenTriangleShell3D) -> Tuple['PointCloud3D', List[float], List[int]]:
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

    def shell_distances_ndarray(self, shells: vmf.OpenTriangleShell3D):
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
    def from_step(cls, step_file: str):
        step = vstep.Step(step_file)
        points = step.to_points()
        return cls(points)

    def plot(self, ax=None, color='k'):
        """
        Plot the cloud 3d.

        """
        ax = self.points[0].plot(ax=ax)
        for point in self.points[1::100]:
            point.plot(ax=ax, color=color)

        return ax

    def extended_cloud(self, distance_extended: float):
        # it works if distance_extended >= 0
        spheres, extended_points = [], []
        for point in self.points:
            extended_zone = p3d.Sphere(point, distance_extended)
            sphere_primitive = extended_zone.shell_faces[0]

            spheres.append(vmf.ClosedShell3D([sphere_primitive]))

            extended_points.extend(sphere_primitive.triangulation().points)

        for sphere in spheres:
            clean_extended_zone = []
            for point in extended_points:
                if not sphere.point_belongs(point):
                    clean_extended_zone.append(point)
            extended_points = clean_extended_zone

        return extended_points

    @staticmethod
    def offset_to_shell(positions_plane: List[vmf.Plane3D],
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
        if ax is None:
            _, ax = plt.subplots()
        for point in self.points:
            point.plot(ax=ax, color=color)
        return ax

    def to_polygon(self, convex=False):
        if not self.points:
            return None

        # polygon = vmw.ClosedPolygon2D.convex_hull_points(self.points)
        if convex:
            polygon = vmw.ClosedPolygon2D.points_convex_hull(self.points)
        else:
            polygon = vmw.ClosedPolygon2D.concave_hull(self.points, -0.2, 0.000005)

        if polygon is None or math.isclose(polygon.area(), 0, abs_tol=1e-6):
            return None
        return polygon

    def bounding_rectangle(self):
        x_list, y_list = [p.x for p in self.points], [p.y for p in self.points]
        return min(x_list), max(x_list), min(y_list), max(y_list)

    def simplify(self, resolution=5):
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
        Generate a n_points x 2 matrix of coordinates.
        """
        return [point.coordinates() for point in self.points]
