"""
Class for voxel representation of volmdlr models
"""
import warnings
from typing import List, Set, Tuple, Iterable, Dict

import numpy as np
from dessia_common.core import PhysicalObject
from tqdm import tqdm
from volmdlr import Point3D, Vector3D, Point2D
from volmdlr.core import VolumeModel
from volmdlr.faces import Triangle3D, PlaneFace3D
from volmdlr.surfaces import Surface2D, PLANE3D_OYZ, PLANE3D_OXY, PLANE3D_OXZ
from volmdlr.wires import ClosedPolygon2D
from volmdlr.shells import ClosedShell3D, ClosedTriangleShell3D

# Typings
Point = Tuple[float, ...]
Triangle = Tuple[Point, ...]
Segment = Tuple[Point, ...]


class Voxelization(PhysicalObject):
    """Class for creation and manipulation of voxelization of volmdlr geometry."""

    def __init__(self, voxels_centers: Set[Point], voxel_size: float, octree_root: "OctreeNode" = None, name: str = ""):
        """
        Initialize the Voxelization.

        :param voxels_centers: The set of points representing voxel centers.
        :type voxels_centers: set[tuple[float, float, float]]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param octree_root: The octree root used to create the voxelization, if so.
        :type octree_root: OctreeNode, optional
        :param name: The name of the Voxelization.
        :type name: str, optional
        """
        self.voxels_centers = voxels_centers
        self.voxel_size = voxel_size
        self.octree_root = octree_root

        PhysicalObject.__init__(self, name=name)

    @classmethod
    def from_closed_triangle_shell(
        cls, closed_triangle_shell: ClosedTriangleShell3D, voxel_size: float, name: str = ""
    ) -> "Voxelization":
        """
        Create a Voxelization from a ClosedTriangleShell3D.

        :param closed_triangle_shell: The ClosedTriangleShell3D to voxelize.
        :type closed_triangle_shell: ClosedTriangleShell3D
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param name: The name of the Voxelization.
        :type name: str, optional

        :return: The created Voxelization.
        :rtype: Voxelization
        """
        triangles = cls._closed_triangle_shell_to_triangles(closed_triangle_shell)
        voxels = cls._triangles_to_voxels(triangles, voxel_size)
        return cls(voxels, voxel_size, name=name)

    @classmethod
    def from_closed_shell(cls, closed_shell: ClosedShell3D, voxel_size: float, name: str = "") -> "Voxelization":
        """
        Create a Voxelization from a ClosedShell3D.

        :param closed_shell: The ClosedShell3D to voxelize.
        :type closed_shell: ClosedShell3D
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param name: The name of the Voxelization.
        :type name: str, optional

        :return: The created Voxelization.
        :rtype: Voxelization
        """
        triangles = cls._closed_shell_to_triangles(closed_shell)
        voxels = cls._triangles_to_voxels(triangles, voxel_size)
        return cls(voxels, voxel_size, name=name)

    @classmethod
    def from_volume_model(cls, volume_model: VolumeModel, voxel_size: float, name: str = "") -> "Voxelization":
        """
        Create a Voxelization from a VolumeModel.

        :param volume_model: The VolumeModel to voxelize.
        :type volume_model: VolumeModel
        :param voxel_size: The voxel edges size.
        :type voxel_size: float
        :param name: The name of the Voxelization.
        :type name: str, optional

        :return: The created Voxelization.
        :rtype: Voxelization
        """
        triangles = cls._volume_model_to_triangles(volume_model)
        voxels = cls._triangles_to_voxels(triangles, voxel_size)
        return cls(voxels, voxel_size, name=name)

    @staticmethod
    def _closed_triangle_shell_to_triangles(closed_triangle_shell: ClosedTriangleShell3D) -> List[Triangle]:
        """
        Helper method to convert a ClosedTriangleShell3D to a list of triangles.

        :param closed_triangle_shell: The ClosedTriangleShell3D to convert to triangles.
        :type closed_triangle_shell: ClosedTriangleShell3D

        :return: The list of triangles extracted from the ClosedTriangleShell3D.
        :rtype: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        """
        return list(
            tuple(tuple([point.x, point.y, point.z]) for point in [tr.point1, tr.point2, tr.point3])
            for tr in closed_triangle_shell.primitives
        )

    @staticmethod
    def _closed_shell_to_triangles(closed_shell: ClosedShell3D) -> List[Triangle]:
        """
        Helper method to convert a ClosedShell3D to a list of triangles.
        It uses the "triangulation" method to triangulate the ClosedShell3D.

        :param closed_shell: The ClosedShell3D to convert to triangles.
        :type closed_shell: ClosedShell3D

        :return: The list of triangles extracted from the triangulated ClosedShell3D.
        :rtype: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        """
        triangulation = closed_shell.triangulation()
        return list(
            tuple(
                tuple([point.x, point.y, point.z])
                for point in [
                    triangulation.points[triangle[0]],
                    triangulation.points[triangle[1]],
                    triangulation.points[triangle[2]],
                ]
            )
            for triangle in triangulation.triangles
        )

    @staticmethod
    def _volume_model_to_triangles(volume_model: VolumeModel) -> List[Triangle]:
        """
        Helper method to convert a VolumeModel to a list of triangles.
        It uses the "triangulation" method to triangulate the primitives of the VolumeModel.

        :param volume_model: The VolumeModel to convert to triangles.
        :type volume_model: VolumeModel

        :return: The list of triangles extracted from the triangulated primitives of the VolumeModel.
        :rtype: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        """
        triangles = []
        for primitive in volume_model.primitives:
            triangulation = primitive.triangulation()
            triangles += list(
                tuple(
                    tuple([point.x, point.y, point.z])
                    for point in [
                        triangulation.points[triangle[0]],
                        triangulation.points[triangle[1]],
                        triangulation.points[triangle[2]],
                    ]
                )
                for triangle in triangulation.triangles
            )
        return triangles

    @staticmethod
    def triangle_voxel_intersection(triangle: Triangle, voxel_center: Point, voxel_extents: List[float]) -> bool:
        """
        Helper method to compute if there is an intersection between a 3D triangle and a voxel.

        :param triangle: The triangle to check if it intersects with the voxel.
        :type: triangle: tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]
        :param voxel_center: The center point of the voxel.
        :type voxel_center: tuple[float, float, float]
        :param voxel_extents: The extents of the voxel in each direction.
        :type voxel_extents: list[float, float, float]

        :return: True if there is an intersection, False otherwise.
        :rtype: bool
        """
        # Convert to numpy ndarray
        triangle = np.array(triangle)
        voxel_center = np.array(voxel_center)
        voxel_extents = np.array(voxel_extents)

        # Translate triangle as conceptually moving AABB to origin
        v = triangle - voxel_center

        # Compute edge vectors for triangle
        f = np.empty_like(triangle)
        f[0] = triangle[1] - triangle[0]
        f[1] = triangle[2] - triangle[1]
        f[2] = triangle[0] - triangle[2]

        for i in range(3):
            a = np.zeros(3)
            # axis a0i
            a[np.mod(i + 1, 3)] = -f[i][np.mod(i + 2, 3)]
            a[np.mod(i + 2, 3)] = f[i][np.mod(i + 1, 3)]
            p = np.dot(v, a)
            r = voxel_extents[np.mod(i + 1, 3)] * abs(f[i][np.mod(i + 2, 3)]) + voxel_extents[np.mod(i + 2, 3)] * abs(
                f[i][np.mod(i + 1, 3)]
            )
            if (max(-max(p), min(p))) > r:
                return False

            # axis a1i
            a[np.mod(i + 2, 3)] = -f[i][np.mod(i, 3)]
            a[np.mod(i, 3)] = f[i][np.mod(i + 2, 3)]
            p = np.dot(v, a)
            r = voxel_extents[np.mod(i, 3)] * abs(f[i][np.mod(i + 2, 3)]) + voxel_extents[np.mod(i + 2, 3)] * abs(
                f[i][np.mod(i, 3)]
            )
            if (max(-max(p), min(p))) > r:
                return False

            # axis a2i
            a[np.mod(i, 3)] = -f[i][np.mod(i + 1, 3)]
            a[np.mod(i + 1, 3)] = f[i][np.mod(i, 3)]
            p = np.dot(v, a)
            r = voxel_extents[np.mod(i, 3)] * abs(f[i][np.mod(i + 1, 3)]) + voxel_extents[np.mod(i + 1, 3)] * abs(
                f[i][np.mod(i, 3)]
            )
            if (max(-max(p), min(p))) > r:
                return False

        # Test the three axes corresponding to the face normals of AABB b (category 1)
        if np.any(np.max(v, axis=0) < -voxel_extents) or np.any(np.min(v, axis=0) > voxel_extents):
            return False

        # Test separating axis corresponding to triangle face normal (category 2)
        plane_normal = np.cross(f[0], f[1])
        plane_distance = np.dot(plane_normal, v[0])

        # Compute the projection interval radius of b onto L(t) = b.c + t * p.n
        r = (
            voxel_extents[0] * abs(plane_normal[0])
            + voxel_extents[1] * abs(plane_normal[1])
            + voxel_extents[2] * abs(plane_normal[2])
        )

        # Intersection occurs when plane distance falls within [-r,+r] interval
        if abs(plane_distance) > r:
            return False

        return True

    @staticmethod
    def _aabb_intersecting_boxes(min_point: Point, max_point: Point, voxel_size: float) -> List[Point]:
        """
        Helper method to compute the center of the voxels that intersect with a given axis aligned
        bounding box (defined by 2 points).

        :param min_point: The minimum point of the bounding box.
        :type min_point: tuple[float, float, float]
        :param max_point: The maximum point of the bounding box.
        :type max_point: tuple[float, float, float]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float

        :return: A list of the centers of the intersecting voxels.
        :rtype: list[tuple[float, float, float]]
        """
        # Calculate the indices of the cubes that intersect with the bounding box
        x_indices = range(int(min_point[0] / voxel_size) - 1, int(max_point[0] / voxel_size) + 1)
        y_indices = range(int(min_point[1] / voxel_size) - 1, int(max_point[1] / voxel_size) + 1)
        z_indices = range(int(min_point[2] / voxel_size) - 1, int(max_point[2] / voxel_size) + 1)

        # Create a list of the centers of all the intersecting voxels
        centers = []
        for x in x_indices:
            for y in y_indices:
                for z in z_indices:
                    center = tuple(round((_ + 1 / 2) * voxel_size, 6) for _ in [x, y, z])
                    centers.append(center)

        return centers

    @staticmethod
    def _voxels_intersecting_voxels(voxel_centers_array: np.ndarray, voxel_size: float) -> Set[Point]:
        """
        Helper method to compute the center of the voxels that intersect with a given array of voxels.
        The returned voxels are part of an implicit 3D grid defined by the voxel size, whereas the given voxels centers
        are not.

        This method is used when translating or rotating the voxels get back to the implicit 3D grid and perform
        efficient boolean operations (thanks to voxels defined in the same grid).

        :param voxel_centers_array: The array of voxel centers.
        :type voxel_centers_array: numpy.ndarray
        :param voxel_size: The voxel edges size.
        :type voxel_size: float

        :return: A set of the centers of the intersecting voxels.
        :rtype: set[tuple[float, float, float]]
        """
        # Compute the voxel indices
        indices = np.floor(voxel_centers_array / voxel_size).astype(int)

        # Compute unique indices to avoid duplicates
        unique_indices = np.unique(indices, axis=0)

        # Convert back to voxel centers
        centers = np.around((unique_indices + 0.5) * voxel_size, 6)

        return set(tuple(center) for center in centers)

    @staticmethod
    def _triangles_to_voxels(triangles: List[Triangle], voxel_size: float) -> Set[Point]:
        """
        Helper method to compute all the voxels intersecting with a given list of triangles.

        :param triangles: The triangles to compute the intersecting voxels.
        :type triangles: list[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float

        :return: The centers of the voxels that intersect with the triangles.
        :rtype: set[tuple[float, float, float]]
        """
        bbox_centers = set()

        for triangle in tqdm(triangles):
            min_point = tuple(min(p[i] for p in triangle) for i in range(3))
            max_point = tuple(max(p[i] for p in triangle) for i in range(3))

            for bbox_center in Voxelization._aabb_intersecting_boxes(min_point, max_point, voxel_size):
                if bbox_center not in bbox_centers:
                    if Voxelization.triangle_voxel_intersection(
                        triangle,
                        bbox_center,
                        [voxel_size for _ in range(3)],
                    ):
                        bbox_centers.add(bbox_center)

        return bbox_centers

    @staticmethod
    def _voxel_triangular_faces(voxel_center: Point, voxel_size: float) -> List[Triangle]:
        """
        Helper method to compute the 12 triangular faces that compose a voxel, for visualization.

        :param voxel_center: The voxel center point.
        :type voxel_center: tuple[float, float, float]
        :param voxel_size: The voxel edges size.
        :type voxel_size: float

        :return: The 12 triangles representing the 6 faces of the given voxel.
        """
        x, y, z = voxel_center
        sx, sy, sz = voxel_size, voxel_size, voxel_size
        hx, hy, hz = sx / 2, sy / 2, sz / 2
        faces = [
            # Front face
            tuple(
                [
                    (x - hx, y - hy, z + hz),
                    (x - hx, y + hy, z + hz),
                    (x + hx, y + hy, z + hz),
                ]
            ),
            tuple(
                [
                    (x - hx, y - hy, z + hz),
                    (x + hx, y + hy, z + hz),
                    (x + hx, y - hy, z + hz),
                ]
            ),
            # Back face
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x - hx, y + hy, z - hz),
                    (x + hx, y + hy, z - hz),
                ]
            ),
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x + hx, y + hy, z - hz),
                    (x + hx, y - hy, z - hz),
                ]
            ),
            # Left face
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x - hx, y + hy, z - hz),
                    (x - hx, y + hy, z + hz),
                ]
            ),
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x - hx, y + hy, z + hz),
                    (x - hx, y - hy, z + hz),
                ]
            ),
            # Right face
            tuple(
                [
                    (x + hx, y - hy, z - hz),
                    (x + hx, y + hy, z - hz),
                    (x + hx, y + hy, z + hz),
                ]
            ),
            tuple(
                [
                    (x + hx, y - hy, z - hz),
                    (x + hx, y + hy, z + hz),
                    (x + hx, y - hy, z + hz),
                ]
            ),
            # Top face
            tuple(
                [
                    (x - hx, y + hy, z - hz),
                    (x + hx, y + hy, z - hz),
                    (x + hx, y + hy, z + hz),
                ]
            ),
            tuple(
                [
                    (x - hx, y + hy, z - hz),
                    (x + hx, y + hy, z + hz),
                    (x - hx, y + hy, z + hz),
                ]
            ),
            # Bottom face
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x + hx, y - hy, z - hz),
                    (x + hx, y - hy, z + hz),
                ]
            ),
            tuple(
                [
                    (x - hx, y - hy, z - hz),
                    (x + hx, y - hy, z + hz),
                    (x - hx, y - hy, z + hz),
                ]
            ),
        ]

        triangular_faces = [tuple(tuple(round(_, 6) for _ in point) for point in face) for face in faces]
        return triangular_faces

    @staticmethod
    def _triangles_to_segments(
        triangles: Set[Triangle],
    ) -> Tuple[Dict[float, List[Segment]], Dict[float, List[Segment]], Dict[float, List[Segment]]]:
        """
        Helper method to extract the segments of a given list of triangle representing a voxelization.
        The segments are sorted by plane: a plane is defined by its normal vector (X, Y or Z) and abscissa.

        Only the segments representing the relevant contours of a faces are returned.
        These relevant segments are the one that are only present in the list of triangles, because if a segment is
        present twice, this means that it's inside the face and is not interesting for defining the face contour.

        :param triangles: The triangles defining the voxelization.
        :type triangles: Set[Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]]

        :return: The extracted relevant segments, sorted by plane.
        :rtype: tuple[dict[float, list[Segment]], dict[float, list[Segment]], dict[float, list[Segment]]]
        """
        segments_x = {}
        segments_y = {}
        segments_z = {}

        for triangle in triangles:
            if triangle[0][0] == triangle[1][0] == triangle[2][0]:
                segments = segments_x.get(triangle[0][0], set())
                triangle_segments = [
                    (tuple(triangle[0][1:]), tuple(triangle[1][1:])),
                    (tuple(triangle[1][1:]), tuple(triangle[2][1:])),
                    (tuple(triangle[2][1:]), tuple(triangle[0][1:])),
                ]
                for segment in triangle_segments:
                    if segment not in segments and segment[::-1] not in segments:
                        segments.add(segment)
                    else:
                        if segment in segments:
                            segments.remove(segment)
                        else:
                            segments.remove(segment[::-1])

                segments_x[triangle[0][0]] = segments

            elif triangle[0][1] == triangle[1][1] == triangle[2][1]:
                segments = segments_y.get(triangle[0][1], set())
                triangle_segments = [
                    (
                        tuple([triangle[0][0], triangle[0][2]]),
                        tuple([triangle[1][0], triangle[1][2]]),
                    ),
                    (
                        tuple([triangle[1][0], triangle[1][2]]),
                        tuple([triangle[2][0], triangle[2][2]]),
                    ),
                    (
                        tuple([triangle[2][0], triangle[2][2]]),
                        tuple([triangle[0][0], triangle[0][2]]),
                    ),
                ]
                for segment in triangle_segments:
                    if segment not in segments and segment[::-1] not in segments:
                        segments.add(segment)
                    else:
                        if segment in segments:
                            segments.remove(segment)
                        else:
                            segments.remove(segment[::-1])

                segments_y[triangle[0][1]] = segments

            else:
                segments = segments_z.get(triangle[0][2], set())
                triangle_segments = [
                    (tuple(triangle[0][:2]), tuple(triangle[1][:2])),
                    (tuple(triangle[1][:2]), tuple(triangle[2][:2])),
                    (tuple(triangle[2][:2]), tuple(triangle[0][:2])),
                ]
                for segment in triangle_segments:
                    if segment not in segments and segment[::-1] not in segments:
                        segments.add(segment)
                    else:
                        if segment in segments:
                            segments.remove(segment)
                        else:
                            segments.remove(segment[::-1])

                segments_z[triangle[0][2]] = segments

        return segments_x, segments_y, segments_z

    @staticmethod
    def _order_segments(segments: Iterable[Segment]) -> List[List[Segment]]:
        """
        Helper method to order a given set of segments defined in the same plane.
        The segments may define several contours, so this method define a list of lists of ordered segments,
        representing the multiple contours with ordered segments.

        :param segments: The segments to order.
        :type segments: Iterable[tuple[tuple[float, float, float], tuple[float, float, float]]]

        :return: The ordered segments lists.
        :rtype: list[list[tuple[tuple[float, float, float], tuple[float, float, float]]]]
        """
        segments = list(segments)
        ordered_segments_list = []

        while segments:
            # Start a new contour with the first segment
            ordered_segments = [segments.pop()]
            ordered_segments_list.append(ordered_segments)

            while True:
                last_segment = ordered_segments[-1]

                # Find the next segment
                for i, s in enumerate(segments):
                    if s[0] == last_segment[1]:
                        ordered_segments.append(segments.pop(i))
                        break

                    elif s[1] == last_segment[1]:
                        # If the segment is in the opposite direction, reverse it
                        ordered_segments.append((s[1], s[0]))
                        segments.pop(i)
                        break
                else:
                    # No more segments connect to this contour, start a new one
                    break

        return ordered_segments_list

    @staticmethod
    def _triangles_to_closed_polygons(triangles: Set[Triangle]) -> List[Dict[float, List[ClosedPolygon2D]]]:
        """
        Helper method to convert the triangles defined the faces of the voxelization to ClosedPolygon2D, sorted
        by plane (a plane is defined by its normal vector (X, Y or Z) and its abscissa).

        :param triangles: The triangles defining the voxelization.
        :type triangles: Set[Tuple[Tuple[float, float, float], Tuple[float, float, float], Tuple[float, float, float]]]

        :return: All the ClosedPolygon2D defining the faces of the voxelization, sorted by plane.
        :rtype: list[dict[float, list[ClosedPolygon2D]]]
        """
        segments_x, segments_y, segments_z = Voxelization._triangles_to_segments(triangles)
        polygons = [{}, {}, {}]

        for i, segments_i in enumerate([segments_x, segments_y, segments_z]):
            for abscissa, segments in segments_i.items():
                # Order the segments lists
                ordered_segments_list = Voxelization._order_segments(segments)

                # Split the self-intersecting polygons
                split_ordered_segments_list = []
                for ordered_segments in ordered_segments_list:
                    split_ordered_segments_list.extend(Voxelization._split_ordered_segments(ordered_segments))

                polygons[i][abscissa] = [
                    Voxelization._ordered_segments_to_closed_polygon_2d(split_ordered_segments)
                    for split_ordered_segments in split_ordered_segments_list
                ]

        return polygons

    @staticmethod
    def _split_ordered_segments(ordered_segments: List[Segment]) -> List[List[Segment]]:
        """
        Helper method to split the self-crossing polygons into multiple polygons to avoid triangulation errors.

        :param ordered_segments: The ordered segments defining the polygon.
        :type ordered_segments: list[tuple[tuple[float, float, float], tuple[float, float, float]]]

        :return: The split polygon into multiple polygons, defined by segment list.
        :rtype: list[list[tuple[tuple[float, float, float], tuple[float, float, float]]]]
        """
        points = [segment[0] for segment in ordered_segments]

        if len(points) == len(set(points)):
            # If no duplicate points, the polygon is not self-intersecting
            return [ordered_segments]

        graph = Graph(ordered_segments)
        ordered_segments_list = []

        for cycle in graph.find_cycles():
            ordered_segments_list.extend(Voxelization._order_segments(cycle))
            # print(Voxelization._order_segments(cycle))

        polygons = {}
        for ordered_segments in ordered_segments_list:
            polygons[tuple(ordered_segments)] = Voxelization._ordered_segments_to_closed_polygon_2d(ordered_segments)

        # Check for polygon inclusion
        children = {}
        parents = {}
        for polygon in polygons.values():
            children[polygon] = []
            parents[polygon] = []

        for polygon_1 in polygons.values():
            for polygon_2 in polygons.values():
                if polygon_1 != polygon_2:
                    if polygon_1.is_inside(polygon_2):
                        children[polygon_1].append(polygon_2)
                        parents[polygon_2].append(polygon_1)
                    elif polygon_2.is_inside(polygon_1):
                        children[polygon_2].append(polygon_1)
                        parents[polygon_1].append(polygon_2)

        # Actual split polygons have no children in the cycle-defined polygons from the self-intersecting polygon
        candidate_polygons = {}
        for ordered_segments, polygon in polygons.items():
            if len(children[polygon]) == 0:
                candidate_polygons[ordered_segments] = polygon

        # If the polygon has all its segments part of the other polygons segments, it is an inner polygon
        ordered_segments_list = []
        for ordered_segments in candidate_polygons.keys():
            ordered_segments_list.append(list(ordered_segments))

        ordered_segments_list_valid = []
        for i in range(len(ordered_segments_list)):
            all_current_polygon_segments = []
            all_other_polygons_segments = []

            for segment in ordered_segments_list[i]:
                all_current_polygon_segments.append(segment)
                all_current_polygon_segments.append(segment[::-1])

            for segments in ordered_segments_list[:i] + ordered_segments_list[i + 1 :]:
                for segment in segments:
                    all_other_polygons_segments.append(segment)
                    all_other_polygons_segments.append(segment[::-1])

            for segment in all_current_polygon_segments:
                if segment not in all_other_polygons_segments:
                    ordered_segments_list_valid.append(ordered_segments_list[i])
                    break

        return ordered_segments_list_valid

    @staticmethod
    def _simplify_polygon_points(points: List[Point]) -> List[Point]:
        """
        Helper method to simplify the list of points defining a polygon, by removing the points that are not relevant
        for the polygon definition, i.e. the points that are not in a corner of the polygon.

        :param points: The points defining the polygon.
        :type points: list[tuple[float, float, float]]

        :return: The simplified points.
        :rtype: list[tuple[float, float, float]]
        """
        simplified_point = []

        for i in range(len(points) - 1):
            if not (
                points[i - 1][0] == points[i][0] == points[i + 1][0]
                or points[i - 1][1] == points[i][1] == points[i + 1][1]
            ):
                simplified_point.append(points[i])

        if not (points[-2][0] == points[-1][0] == points[0][0] or points[-2][1] == points[-1][1] == points[0][1]):
            simplified_point.append(points[-1])

        return simplified_point

    @staticmethod
    def _ordered_segments_to_closed_polygon_2d(ordered_segments: Iterable[Segment]) -> ClosedPolygon2D:
        """
        Helper method to convert an iterable of ordered segments to a ClosedPolygon2D.

        :param ordered_segments: The segments that compose the ClosedPolygon2D.
        :type ordered_segments: Iterable[tuple[tuple[float, float, float], tuple[float, float, float]]]

        :return: The created ClosedPolygon2D.
        :rtype: ClosedPolygon2D
        """
        return ClosedPolygon2D(
            points=[
                Point2D(*point)
                for point in Voxelization._simplify_polygon_points([segment[0] for segment in ordered_segments])
            ]
        )

    def to_triangles(self) -> Set[Triangle]:
        """
        Convert the voxelization to triangles for display purpose.

        Only the relevant faces are returned (i.e. the faces that are not at the interface of two different voxel,
        i.e. the faces that are only present once in the list of triangles representing the triangulated voxels).

        :return: The triangles representing the voxelization.
        :rtype: set[tuple[tuple[float, float, float], tuple[float, float, float], tuple[float, float, float]]]
        """
        triangles = set()

        for voxel in self.voxels_centers:
            for triangle in self._voxel_triangular_faces(voxel, self.voxel_size):
                if triangle not in triangles:
                    triangles.add(triangle)
                else:
                    triangles.remove(triangle)

        return triangles

    def to_closed_shell(self) -> ClosedShell3D:
        """
        Convert the voxelization to a ClosedShell3D for display purpose.

        This ClosedShell3D is only composed of PlaneFace3D that represent the faces of the voxelized geometry, and
        are made to be the lightest geometrical representation of the voxelization (in terms of number of triangles
        in the babylonjs model).

        :return: The created ClosedShell3D with only PlaneFace3D.
        :rtype: ClosedShell3D
        """
        polygons = self._triangles_to_closed_polygons(self.to_triangles())
        planes = [PLANE3D_OYZ, PLANE3D_OXZ, PLANE3D_OXY]
        faces = []

        for i, polygons_i in enumerate(polygons):
            for abscissa, polygons in polygons_i.items():
                offset = [0, 0, 0]
                offset[i] = abscissa
                offset = Vector3D(*offset)

                plane = planes[i].translation(offset)

                if len(polygons) == 1:
                    faces.append(PlaneFace3D(plane, Surface2D(polygons[0], [])))
                else:
                    # Check for polygon inclusion
                    children = {}
                    parents = {}
                    for polygon in polygons:
                        children[polygon] = []
                        parents[polygon] = []

                    for polygon_1 in polygons:
                        for polygon_2 in polygons:
                            if polygon_1 != polygon_2:
                                if polygon_1.is_inside(polygon_2):
                                    children[polygon_1].append(polygon_2)
                                    parents[polygon_2].append(polygon_1)
                                elif polygon_2.is_inside(polygon_1):
                                    children[polygon_2].append(polygon_1)
                                    parents[polygon_1].append(polygon_2)

                    for polygon in polygons:
                        if len(parents[polygon]) // 2 == 0:
                            inner_polygons = [
                                child
                                for child in children[polygon]
                                if set(children[polygon]).intersection(children[child]) == set()
                            ]
                            faces.append(PlaneFace3D(plane, Surface2D(polygon, inner_polygons)))

        return ClosedShell3D(faces, name=self.name)

    def to_closed_triangle_shell(self) -> ClosedTriangleShell3D:
        """
        Convert the voxelization to a ClosedTriangleShell3D for display purpose.

        To create the ClosedTriangleShell3D, this method triangulates each voxel and convert it to Triangle3D.
        The method is robust and fast but creates over-triangulated geometry.

        :return: The created ClosedShell3D with only Triangle3D.
        :rtype: ClosedTriangleShell3D
        """
        triangles3d = [
            Triangle3D(Point3D(*triangle[0]), Point3D(*triangle[1]), Point3D(*triangle[2]))
            for triangle in self.to_triangles()
        ]
        shell = ClosedTriangleShell3D(triangles3d, name=self.name)

        return shell

    def volmdlr_primitives(self, **kwargs):
        """
        Return a list of volmdlr primitives to build up volume model.
        It uses the simplified representation of the voxelization given by the "to_closed_shell" method.
        """
        if len(self.voxels_centers) == 0:
            warnings.warn("Empty voxelization.")
            return []

        return [self.to_closed_shell()]

    def intersection(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Create a voxelization that is the boolean intersection of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the boolean intersection with.
        :type other_voxelization: Voxelization

        :return: The created voxelization resulting from the boolean intersection.
        :rtype: Voxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform intersection.")

        return Voxelization(self.voxels_centers.intersection(other_voxelization.voxels_centers), self.voxel_size)

    def union(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Create a voxelization that is the boolean union of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the boolean union with.
        :type other_voxelization: Voxelization

        :return: The created voxelization resulting from the boolean union.
        :rtype: Voxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform union.")

        return Voxelization(self.voxels_centers.union(other_voxelization.voxels_centers), self.voxel_size)

    def difference(self, other_voxelization: "Voxelization") -> "Voxelization":
        """
        Create a voxelization that is the boolean difference of two voxelization.
        Both voxelization must have same voxel size.

        :param other_voxelization: The other voxelization to compute the boolean difference with.
        :type other_voxelization: Voxelization

        :return: The created voxelization resulting from the boolean difference.
        :rtype: Voxelization
        """
        if self.voxel_size != other_voxelization.voxel_size:
            raise ValueError("Both voxelizations must have same voxel_size to perform difference.")

        return Voxelization(self.voxels_centers.difference(other_voxelization.voxels_centers), self.voxel_size)

    @staticmethod
    def _rotation_matrix(axis: Vector3D, angle: float) -> np.array:
        """
        Helper method that compute a rotation matrix from an axis and a radians angle.

        :param axis: The rotation axis.
        :type axis: Vector3D
        :param angle: The rotation angle.
        :type angle: float

        :return: The computed rotation matrix.
        :rtype: numpy.array
        """
        axis = np.array([axis.x, axis.y, axis.z])

        axis = axis / np.linalg.norm(axis)
        a = np.cos(angle / 2.0)

        b, c, d = -axis * np.sin(angle / 2.0)
        aa, bb, cc, dd = a * a, b * b, c * c, d * d
        bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d

        return np.array(
            [
                [aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc],
            ]
        )

    def rotation(self, center: Point3D, axis: Vector3D, angle: float):
        rotation_matrix = self._rotation_matrix(axis, angle)
        voxel_array = np.array(list(self.voxels_centers)) - np.array([center.x, center.y, center.z])
        rotated_voxels = np.dot(voxel_array, rotation_matrix.T)
        rotated_voxels += np.array([center.x, center.y, center.z])

        intersecting_voxels = self._voxels_intersecting_voxels(rotated_voxels, self.voxel_size)

        return Voxelization(intersecting_voxels, self.voxel_size)

    def translation(self, offset: Vector3D):
        voxel_array = np.array(list(self.voxels_centers))
        translated_voxels = voxel_array + np.array([offset.x, offset.y, offset.z])

        intersecting_voxels = self._voxels_intersecting_voxels(translated_voxels, self.voxel_size)

        return Voxelization(intersecting_voxels, self.voxel_size)


class OctreeNode:
    """Class representing an octree node for octree voxelization purpose."""
    def __init__(self, center: Point, size: float, depth: int, max_depth: int):
        self.children = []
        self.center = center
        self.size = size
        self.depth = depth
        self.max_depth = max_depth

    def subdivide(self, triangles: List[Triangle]):
        children = []
        if self.depth < self.max_depth:
            half_size = self.size / 2
            for i in range(2):
                for j in range(2):
                    for k in range(2):
                        # calculate the center of the child node
                        child_center = (
                            self.center[0] + (i - 0.5) * half_size,
                            self.center[1] + (j - 0.5) * half_size,
                            self.center[2] + (k - 0.5) * half_size,
                        )

                        # create a new OctreeNode for the child
                        child_node = OctreeNode(
                            child_center,
                            half_size,
                            self.depth + 1,
                            self.max_depth,
                        )

                        # check if the child node intersects with the mesh
                        if self.depth != 0:
                            # add the child node to the list of children
                            self.children.extend(self.process_node(child_node, triangles))
                        else:
                            children.append(child_node)

            if self.depth == 0:
                for node in tqdm(self.process_node(node, triangles) for node in children):
                    self.children.extend(node)

        else:
            return  # reached max depth, do not subdivide further.

    def intersecting_triangles(self, triangles: List[Triangle]) -> List[Triangle]:
        intersecting_triangles = [
            triangle
            for triangle in triangles
            if Voxelization.triangle_voxel_intersection(triangle, self.center, [self.size for _ in range(3)])
        ]

        return intersecting_triangles

    def get_leaf_centers(self):
        if self.depth == self.max_depth:  # if max depth reached, it is a leaf node
            return [self.center]
        else:
            centers = []
            for child in self.children:
                centers += child.get_leaf_centers()
            return centers

    @staticmethod
    def process_node(node, mesh):
        triangles = node.intersecting_triangles(mesh)
        if len(triangles) > 0:
            # recursively subdivide the child node
            node.subdivide(triangles)
            return [node]
        return []


class Graph:
    """Helper class for subdividing self-crossing polygons using graph theory algorithm."""

    def __init__(self, edges):
        """Initialize the graph from edges."""
        self.graph = self.build_graph(edges)
        self.visited = {node: False for node in self.graph}

    @staticmethod
    def build_graph(edges):
        """Builds adjacency list representation of graph from edge list."""
        graph = {}
        for edge in edges:
            if edge[0] not in graph:
                graph[edge[0]] = []
            if edge[1] not in graph:
                graph[edge[1]] = []
            graph[edge[0]].append(edge[1])
            graph[edge[1]].append(edge[0])  # assuming undirected graph
        return graph

    def find_cycles(self):
        """Finds all unique cycles in the graph."""
        cycles = []
        for node in self.graph:
            self._dfs(node, node, [], cycles)
        return self._convert_to_edge_cycles(cycles)

    def _dfs(self, node, start, path, cycles):
        """Performs depth-first search to find cycles."""
        self.visited[node] = True
        path.append(node)
        for neighbor in self.graph[node]:
            if not self.visited[neighbor]:
                self._dfs(neighbor, start, path, cycles)
            elif neighbor == start and len(path) >= 3:
                cycles.append(path[:])
        path.pop()
        self.visited[node] = False

    @staticmethod
    def _convert_to_edge_cycles(cycles):
        """Converts node-based cycles to edge-based cycles, removing duplicates."""
        edge_cycles_set = set()
        for cycle in cycles:
            edge_cycle = []
            for i in range(len(cycle)):
                edge_cycle.append(tuple(sorted([cycle[i], cycle[(i + 1) % len(cycle)]])))
            edge_cycle.sort()  # sort the edges in the cycle
            edge_cycles_set.add(tuple(edge_cycle))  # add it to the set
        return [list(cycle) for cycle in edge_cycles_set]
