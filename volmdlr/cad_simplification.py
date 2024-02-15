"""
volmdlr cad simplification module.
"""
import math
from typing import Union

import numpy as np
from CGAL.CGAL_Alpha_wrap_3 import alpha_wrap_3
from CGAL.CGAL_Kernel import Point_3
from CGAL.CGAL_Polyhedron_3 import Polyhedron_3
from dessia_common.core import DessiaObject

import volmdlr
from volmdlr import Point3D
from volmdlr.cloud import PointCloud3D
from volmdlr.core import VolumeModel
from volmdlr.discrete_representation import MatrixBasedVoxelization
from volmdlr.display import Mesh3D
from volmdlr.faces import Triangle3D
from volmdlr.primitives3d import ExtrudedProfile
from volmdlr.shells import ClosedTriangleShell3D, OpenTriangleShell3D
from volmdlr.wires import Contour2D


class Simplify(DessiaObject):
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
        points = [Point3D(*point) for point in self.volume_model.to_mesh().vertices]

        point_cloud3d = PointCloud3D([volmdlr.Point3D(*point) for point in points])
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

            frame = volmdlr.Frame3D(((center - 0.5 * length) * w_vector).to_point(), u_vector, v_vector, w_vector)
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
        preserve_shells: bool = True,
    ):
        """
        Simplify the VolumeModel using the 'triangle decimation' method, and return it.

        To do so, it triangulates each shell of the given volume model and then decimates it at the given target ratio.

        Note: threshold = alpha * pow(iteration + k, aggressiveness)

        :param target_ratio: Target number of triangles. Not used if `lossless` is True.
        :type target_ratio: float
        :param update_rate: Number of iterations between each update. If `lossless` flag is set to True, rate is 1.
        :type update_rate: int
        :param aggressiveness: Parameter controlling the growth rate of the threshold at each iteration when `lossless`
            is False.
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
        :param preserve_shells: Argument to chose if you want to keep the shell structures (same number of shells in the
          returned volume model), or if you want a volume model with only one shell including the entire input geometry.
        :type preserve_shells: bool

        :return: The decimated VolumeModel.
        :rtype: VolumeModel
        """
        # pylint: disable=too-many-arguments

        decimated_meshes = []

        if preserve_shells:
            meshes = self.volume_model.to_mesh_list()
        else:
            meshes = [self.volume_model.to_mesh()]

        for mesh in meshes:
            decimated_meshes.append(
                mesh.round_vertices()
                .merge_vertices()
                .decimate(
                    target_count=int(target_ratio * mesh.n_triangles),
                    update_rate=update_rate,
                    aggressiveness=aggressiveness,
                    max_iterations=max_iterations,
                    verbose=verbose,
                    lossless=lossless,
                    threshold_lossless=threshold_lossless,
                    alpha=alpha,
                    k=k,
                    preserve_border=preserve_border,
                )
            )

            decimated_meshes[-1].name = mesh.name
            decimated_meshes[-1].color = mesh.color
            decimated_meshes[-1].alpha = mesh.alpha

        return VolumeModel(decimated_meshes)


class AlphaWrapSimplify(Simplify):
    """CAD simplification using CGAL 'alpha wrap' method."""

    def simplify(
        self,
        relative_alpha: int = 10,
        relative_offset: int = 300,
        preserve_shells: bool = False,
        equal_alpha: bool = True,
    ) -> VolumeModel:
        """
        Simplify the volume model using the CGAL 'alpha wrap' method, and return it.

        :param relative_alpha: Control the output complexity, by defining the diameter of cavities to be explored
          relatively to the size of the geometry.
        :type relative_alpha: int
        :param relative_offset: Control how close to the input geometry the simplification is.
        :type relative_offset: int
        :param preserve_shells: Argument to chose if you want to keep the shell structures (same number of shells in the
          returned volume model), or if you want a volume model with only one shell including the entire input geometry.
        :type preserve_shells: bool
        :param equal_alpha: Control if the alpha is relative to the entire volume model bbox (equal for all shells) or
          relative to each shell bbox. No effect if preserve_shells is False.
        :type equal_alpha: bool

        :return: The simplified volume model.
        :rtype: VolumeModel
        """
        wrapped_shells = []

        bbox, alpha, offset = None, None, None

        if preserve_shells:
            meshes = self.volume_model.to_mesh_list()
        else:
            meshes = [self.volume_model.to_mesh()]

            bbox = meshes[0].bounding_box

        if equal_alpha:
            if not bbox:
                bbox = self.volume_model.bounding_box
            diag_length = math.sqrt(
                (bbox.xmax - bbox.xmin) ** 2 + (bbox.ymax - bbox.ymin) ** 2 + (bbox.zmax - bbox.zmin) ** 2
            )
            alpha = diag_length / relative_alpha
            offset = diag_length / relative_offset

        for mesh in meshes:
            vertices = [Point_3(x, y, z) for (x, y, z) in mesh.vertices]
            triangles = mesh.triangles.tolist()

            if not equal_alpha:
                bbox = mesh.bounding_box
                diag_length = math.sqrt(
                    (bbox.xmax - bbox.xmin) ** 2 + (bbox.ymax - bbox.ymin) ** 2 + (bbox.zmax - bbox.zmin) ** 2
                )
                alpha = diag_length / relative_alpha
                offset = diag_length / relative_offset

            wrap = Polyhedron_3()
            alpha_wrap_3(vertices, triangles, alpha, offset, wrap)

            wrapped_shells.append(self._polyhedron_to_mesh(wrap))
            wrapped_shells[-1].name = mesh.name
            wrapped_shells[-1].color = mesh.color
            wrapped_shells[-1].alpha = mesh.alpha

        return VolumeModel(wrapped_shells)

    @staticmethod
    def _polyhedron_to_shell(polyhedron: Polyhedron_3) -> Union[ClosedTriangleShell3D, OpenTriangleShell3D]:
        """Convert a CGAL Polyhedron_3 to a volmdlr shell."""
        triangles = []

        for facet in polyhedron.facets():  # Iterate over each face
            vertices = []
            halfedge = facet.halfedge()  # Starting half-edge
            start_halfedge = halfedge

            while True:
                point = halfedge.vertex().point()
                vertices.append(Point3D(point.x(), point.y(), point.z()))
                halfedge = halfedge.next()  # Move to the next half-edge
                if halfedge == start_halfedge:
                    break  # Completed one loop around the facet

            if len(vertices) == 3:  # Check if it's a triangle
                try:
                    triangles.append(Triangle3D(vertices[0], vertices[1], vertices[2]))
                except ZeroDivisionError:
                    pass

        if polyhedron.is_closed():
            return ClosedTriangleShell3D(triangles)
        return OpenTriangleShell3D(triangles)

    @staticmethod
    def _polyhedron_to_mesh(polyhedron: Polyhedron_3) -> Mesh3D:
        """Convert a CGAL Polyhedron_3 to a Mesh3D."""
        points = []
        faces = []
        point_map = {}  # Map CGAL points to their triangles

        for facet in polyhedron.facets():
            face_triangles = []
            halfedge = facet.halfedge()
            start_halfedge = halfedge

            while True:
                point = halfedge.vertex().point()
                point_tuple = (point.x(), point.y(), point.z())

                # Add point if it's not already in the list
                if point_tuple not in point_map:
                    points.append(point_tuple)
                    point_map[point_tuple] = len(points) - 1

                face_triangles.append(point_map[point_tuple])
                halfedge = halfedge.next()
                if halfedge == start_halfedge:
                    break

            faces.append(face_triangles)

        points_array = np.array(points)
        faces_array = np.array(faces)

        return Mesh3D(points_array, faces_array, name="")
