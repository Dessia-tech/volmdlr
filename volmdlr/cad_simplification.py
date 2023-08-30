"""
volmdlr cad simplification module.
"""
from abc import ABC

import numpy as np
import pyfqmr
from dessia_common.core import DessiaObject
from typing import Tuple

import volmdlr
from volmdlr.cloud import PointCloud3D
from volmdlr.core import VolumeModel
from volmdlr.discrete_representation import MatrixBasedVoxelization
from volmdlr.primitives3d import ExtrudedProfile
from volmdlr.shells import ClosedTriangleShell3D
from volmdlr.wires import Contour2D


class Simplify(ABC, DessiaObject):
    """CAD simplification abstract class."""

    def __init__(self, volume_model: VolumeModel, name: str = ""):
        """
        Initialize an instance of the Simplify class.

        :param volume_model: The volume model to simplify.
        :type volume_model: VolumeModel
        :param name: A name for the TripleExtrusionSimplify instance, optional.
        :type name: str
        """
        self.volume_model = volume_model

        DessiaObject.__init__(self, name=name)


class TripleExtrusionSimplify(Simplify):
    """CAD simplification based on 'triple extrusion' method."""

    def simplify(self) -> VolumeModel:
        """
        Simplify the volume model using the 'triple extrusion' method, and return it.

        :return: The simplified volume model.
        :rtype: VolumeModel
        """
        points = []
        for primitive in self.volume_model.primitives:
            tri = primitive.triangulation()
            points.extend(tri.points)

        point_cloud3d = PointCloud3D(points)
        simplified_volume_model = VolumeModel(
            [self.extrusion_union_cloud_simplifier(point_cloud3d)], name=f"{self.volume_model.name} voxel simplified"
        )

        return simplified_volume_model

    @staticmethod
    def extrusion_union_cloud_simplifier(point_cloud3d: PointCloud3D) -> ExtrudedProfile:
        """
        Simplify a point cloud using extrusion and union operations.

        :param point_cloud3d: The 3D point cloud to simplify.
        :type point_cloud3d: PointCloud3D

        :return: The simplified shell.
        :rtype: ExtrudedProfile
        """
        simplified_shell = None
        list_shells = []

        bbox = point_cloud3d.bounding_box
        lx = bbox.xmax - bbox.xmin
        ly = bbox.ymax - bbox.ymin
        lz = bbox.zmax - bbox.zmin

        for w_vector, length, center in zip([volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D], [lx, ly, lz], bbox.center):
            u_vector = w_vector.random_unit_normal_vector()
            v_vector = w_vector.cross(u_vector)

            cloud2d = point_cloud3d.to_2d(volmdlr.O3D, u_vector, v_vector)
            polygon2d = cloud2d.to_polygon(convex=True)
            contour2d = Contour2D(polygon2d.line_segments)

            frame = volmdlr.Frame3D((center - 0.5 * length) * w_vector, u_vector, v_vector, w_vector)
            dir_shell = ExtrudedProfile(frame, contour2d, [], length)
            dir_shell.merge_faces()

            list_shells.append(dir_shell)
            if simplified_shell is None:
                simplified_shell = dir_shell
            else:
                list_shells = simplified_shell.intersection(dir_shell)
                simplified_shell = list_shells[0]

        return simplified_shell


class VoxelizationSimplify(Simplify):
    """CAD simplification using 'voxelization' method."""

    def simplify(self, voxel_size: float, fill: bool = True) -> VolumeModel:
        """
        Simplify the volume model using the 'voxelization' method, and return it.

        :param voxel_size: The size of the voxels used for simplification.
        :type voxel_size: float
        :param fill: Control the filling of the voxelization, default is True.
        :type fill: bool

        :return: The simplified volume model.
        :rtype: VolumeModel
        """
        voxelization = MatrixBasedVoxelization.from_volume_model(
            self.volume_model, voxel_size, name=f"{self.volume_model.name} voxel simplified"
        )
        if fill:
            voxelization = voxelization.fill_enclosed_voxels()

        return VolumeModel(voxelization.volmdlr_primitives())


class TriangleDecimationSimplify(Simplify):
    """CAD simplification based on 'triangle decimation' method."""

    def simplify(
        self,
        target_ratio: float,
        update_rate: int = 5,
        aggressiveness: float = 7.0,
        max_iterations: int = 100,
        verbose: bool = False,
        lossless: bool = False,
        threshold_lossless: float = 1e-3,
        alpha: float = 1e-9,
        k: int = 3,
        preserve_border: bool = True,
    ):
        """
        Simplify the VolumeModel using the 'triangle decimation' method, and return it.

        To do so, it triangulates each shell of the given volume model and then decimates it at the given target ratio.

        Note: threshold = alpha * pow(iteration + k, aggressiveness)

        :param target_ratio: Target number of triangles. Not used if `lossless` is True.
        :type target_ratio: float
        :param update_rate: Number of iterations between each update. If `lossless` flag is set to True, rate is 1.
        :type update_rate: int
        :param aggressiveness: Parameter controlling the growth rate of the threshold at each iteration when `lossless` is False.
        :type aggressiveness: float
        :param max_iterations: Maximal number of iterations.
        :type max_iterations: int
        :param verbose: Control verbosity.
        :type verbose: bool
        :param lossless: Use the lossless simplification method.
        :type lossless: bool
        :param threshold_lossless: Maximal error after which a vertex is not deleted. Only for `lossless` method.
        :type threshold_lossless: float
        :param alpha: Parameter for controlling the threshold growth.
        :type alpha: float
        :param k: Parameter for controlling the threshold growth.
        :type k: int
        :param preserve_border: Flag for preserving vertices on open border.
        :type preserve_border: bool

        :return: The decimated VolumeModel.
        :rtype: VolumeModel
        """
        decimated_shells = []
        simplifier = pyfqmr.Simplify()

        for shell in self.volume_model.get_shells():
            triangulation = shell.triangulation()

            # vertices = np.array([(point.x, point.y, point.z) for point in triangulation.points])
            # triangles = np.array(triangulation.triangles)

            vertices, triangles = self.get_vertices_and_faces(ClosedTriangleShell3D(triangulation.faces))

            simplifier.setMesh(vertices, triangles)
            simplifier.simplify_mesh(
                target_count=int(target_ratio * len(triangles)),
                update_rate=update_rate,
                aggressiveness=aggressiveness,
                max_iterations=max_iterations,
                verbose=verbose,
                lossless=lossless,
                threshold_lossless=threshold_lossless,
                alpha=alpha,
                K=k,
                preserve_border=preserve_border,
            )

            vertices, faces, _ = simplifier.getMesh()

            # valid_faces = []
            # for face1, fa in faces:
            #     if face[0] == face[1] or face[1] == face[2] or face[0] == face[2]:
            #         print("found")
            #     else:
            #         valid_faces.append(faces)

            decimated_shells.append(ClosedTriangleShell3D.from_mesh_data(vertices, faces))

        return VolumeModel(decimated_shells)

    @staticmethod
    def get_vertices_and_faces(
        closed_triangle_shell, round_vertices: bool = False, decimal_places: int = 6
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get the vertices and faces of the shell as numpy arrays.
        :param round_vertices: Whether to round vertices to prevent numerical errors.
        :type round_vertices: bool
        :param decimal_places: Number of decimal places to round vertices.
        :type decimal_places: int
        :return: A tuple containing two numpy arrays - vertices and faces.
        :rtype: Tuple[np.ndarray, np.ndarray]
        """
        vertices = np.array([(point.x, point.y, point.z) for face in closed_triangle_shell.faces for point in face.points])

        if round_vertices:
            vertices = np.round(vertices, decimal_places)

        vertices = np.unique(vertices, axis=0)  # Remove duplicate vertices

        from tqdm import tqdm

        vertex_indices = []
        for face in tqdm(closed_triangle_shell.faces):
            if round_vertices:
                indices = [
                    np.where(np.all(vertices == np.round((point.x, point.y, point.z), decimal_places), axis=1))[0][0]
                    for point in face.points
                ]
            else:
                indices = [
                    np.where(np.all(vertices == (point.x, point.y, point.z), axis=1))[0][0] for point in face.points
                ]

            vertex_indices.append(indices)

        faces = np.array(vertex_indices)

        return vertices, faces
