"""
volmdlr cad simplification module.
"""
from typing import List

from dessia_common.core import DessiaObject

import volmdlr
from volmdlr.cloud import PointCloud3D
from volmdlr.core import VolumeModel
from volmdlr.primitives3d import Block, ExtrudedProfile
from volmdlr.shells import ClosedShell3D, union_list_of_shells
from volmdlr.wires import Contour2D


class OctreeBlockSimplify(DessiaObject):
    """CAD simplification based on octree voxelization."""

    def __init__(self, closed_shell: ClosedShell3D, name: str = ""):
        """
        Initialize an instance of OctreeBlockSimplify.

        :param closed_shell: A closed shell object representing the geometry to simplify.
        :param name: Optional. A name for the OctreeBlockSimplify instance.
        """
        self.closed_shell = closed_shell
        self.block = Block.from_bounding_box(self.closed_shell.bounding_box)

        DessiaObject.__init__(self, name=name)

    def get_outside_blocks(self, depth: int = 3) -> List[List[Block]]:
        """
        Gets octree blocks outside the given closed shell.

        :param depth: Optional. The depth of octree subdivision (default: 3).

        :return: A list of octree blocks that are outside the given closed shell.
        """
        block_volume = self.block.volume()
        min_block_volume = block_volume / (8 ** depth)
        dict_faces_to_box = {box: [] for box in self.block.octree()}
        for box, list_faces in dict_faces_to_box.items():
            for face in self.closed_shell.faces:
                if (
                    box.is_face_intersecting(face)
                    or box.is_face_inside(face)
                    or box.bounding_box.distance_to_bbox(face.bounding_box) <= 1e-6
                ):
                    list_faces.append(face)

        not_intersecting_boxes = []
        inside_boxes = []

        for octant in self.block.octree():
            dict_faces_to_box_ = {octant: dict_faces_to_box[octant]}
            octant_not_intersecting_boxes = []

            while dict_faces_to_box_:
                box, list_faces = next(iter(dict_faces_to_box_.items()))
                if not list_faces:
                    if not box.is_inside_shell(self.closed_shell):
                        octant_not_intersecting_boxes.append(box)
                    else:
                        inside_boxes.append(box)
                    del dict_faces_to_box_[box]
                    continue

                divided_boxes = box.octree()
                if divided_boxes[0].volume() >= min_block_volume:
                    for divided_box in divided_boxes:
                        new_list_faces = [
                            face
                            for face in list_faces
                            if (
                                divided_box.is_face_intersecting(face)
                                or divided_box.is_face_inside(face)
                                or divided_box.bounding_box.distance_to_bbox(face.bounding_box) <= 1e-6
                            )
                        ]
                        dict_faces_to_box_[divided_box] = new_list_faces

                del dict_faces_to_box_[box]

            not_intersecting_boxes.append(octant_not_intersecting_boxes)

        return not_intersecting_boxes

    def simplify(self, depth: int) -> ClosedShell3D:
        """
        Simplify the given closed shell.

        :param depth: The depth of the octree simplification (i.e. the number of subdivision)

        :return: The simplified closed shell.
        """
        outside_blocks = self.get_outside_blocks(depth=depth)
        closed_shells = []

        for octant in outside_blocks:
            if not octant:
                continue
            list_shells = union_list_of_shells(octant)
            closed_shells.extend(list_shells)

        closed_shells_ = union_list_of_shells(closed_shells.copy())
        closed_shells_ = sorted(closed_shells_, key=lambda shell: shell.volume(), reverse=True)
        simplified_closed_shell = self.block.subtract_to_closed_shell(closed_shells_[0])[0]

        return simplified_closed_shell


class TripleExtrusionSimplify(DessiaObject):
    """CAD simplification based on 'triple extrusion' method."""

    def __init__(self, volume_model: VolumeModel, name: str = ""):
        """
        Initialize an instance of TripleExtrusionSimplify.

        :param volume_model: The volume model to simplify.
        :param name (str): Optional. A name for the TripleExtrusionSimplify instance.
        """
        self.volume_model = volume_model

        DessiaObject.__init__(self, name=name)

    def simplify(self) -> VolumeModel:
        """
        Simplify the volume model using the 'triple extrusion' method.

        :return: The simplified volume model.
        """
        points = []
        for primitive in self.volume_model.primitives:
            tri = primitive.triangulation()
            points.extend(tri.points)

        cloud = PointCloud3D(points)
        model = VolumeModel([self.extrusion_union_cloud_simplifier(cloud)])

        return model

    @staticmethod
    def extrusion_union_cloud_simplifier(cloud3d: PointCloud3D) -> ExtrudedProfile:
        """
        Simplify a point cloud using extrusion and union operations.

        :param cloud3d: The 3D point cloud to simplify.

        :return: The simplified shell.
        """
        simplified_shell = None
        list_shells = []

        bbox = cloud3d.bounding_box
        lx = bbox.xmax - bbox.xmin
        ly = bbox.ymax - bbox.ymin
        lz = bbox.zmax - bbox.zmin

        for w_vector, length, center in zip([volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D], [lx, ly, lz], bbox.center):
            u_vector = w_vector.random_unit_normal_vector()
            v_vector = w_vector.cross(u_vector)

            cloud2d = cloud3d.to_2d(volmdlr.O3D, u_vector, v_vector)
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
