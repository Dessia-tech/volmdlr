"""
volmdlr cad simplification module.
"""
from abc import ABC, abstractmethod

from dessia_common.core import DessiaObject

import volmdlr
from volmdlr.cloud import PointCloud3D
from volmdlr.core import VolumeModel
from volmdlr.primitives3d import ExtrudedProfile
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

    @abstractmethod
    def simplify(self, *args, **kwargs) -> VolumeModel:
        """
        Simplify the volume model and return it.

        :return: The simplified volume model.
        :rtype: VolumeModel
        """


class TripleExtrusionSimplify(Simplify):
    """CAD simplification based on 'triple extrusion' method."""

    def simplify(self, *args, **kwargs) -> VolumeModel:
        """
        Simplify the volume model using the 'triple extrusion' method.

        :return: The simplified volume model.
        :rtype: VolumeModel
        """
        points = []
        for primitive in self.volume_model.primitives:
            tri = primitive.triangulation()
            points.extend(tri.points)

        point_cloud3d = PointCloud3D(points)
        simplified_volume_model = VolumeModel([self.extrusion_union_cloud_simplifier(point_cloud3d)])

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
