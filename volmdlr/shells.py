"""volmdlr shells module."""
import traceback
import warnings
from itertools import chain
from typing import List, Tuple

import matplotlib.pyplot as plt
import networkx as nx
import numpy as npy
from trimesh import Trimesh

from dessia_common.core import DessiaObject
import volmdlr.bspline_compiled
import volmdlr.core_compiled
import volmdlr.core
from volmdlr import display, edges, wires, surfaces
import volmdlr.faces
import volmdlr.geometry
from volmdlr.core import point_in_list, edge_in_list, get_edge_index_in_list, get_point_index_in_list


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

    def __init__(self, faces: List[volmdlr.faces.Face3D],
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
        self._bbox = None
        if bounding_box:
            self._bbox = bounding_box

        self._faces_graph = None
        self._vertices_points = None

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
    def vertices_points(self):
        """Gets the shell's vertices points. """
        if self._vertices_points is None:
            vertices_points = []
            for face in self.faces:
                for contour in [face.outer_contour3d] + face.inner_contours3d:
                    for edge in contour.primitives:
                        if not volmdlr.core.point_in_list(edge.start, vertices_points):
                            vertices_points.append(edge.start)
                        if not volmdlr.core.point_in_list(edge.end, vertices_points):
                            vertices_points.append(edge.end)
            self._vertices_points = vertices_points
        return self._vertices_points

    @property
    def faces_graph(self):
        """
        Gets the shells faces graph using networkx.

        :return: return a networkx graph for a shell faces' vertices.
        """
        if not self._faces_graph:
            faces_graph = nx.Graph()
            for face in self.faces:
                for edge in face.outer_contour3d.primitives:
                    edge_start_index = volmdlr.core.get_point_index_in_list(edge.start, self.vertices_points)
                    edge_end_index = volmdlr.core.get_point_index_in_list(edge.end, self.vertices_points)
                    faces_graph.add_edge(edge_start_index, edge_end_index, edge=edge)
            self._faces_graph = faces_graph
        return self._faces_graph

    def to_dict(self, *args, **kwargs):
        """
        Serializes a 3-dimensional open shell into a dictionary.

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
        # Quick fix:
        # ----------------------------------
        name = arguments[0][1:-1]
        if isinstance(arguments[-1], int):
            product = object_dict[arguments[-1]]
            name = product[1:-1]
        # ----------------------------------
        faces = [object_dict[int(face[1:])] for face in arguments[1] if object_dict[int(face[1:])]]
        return cls(faces, name=name)

    def to_step(self, current_id):
        """
        Creates step file entities from volmdlr objects.
        """
        step_content = ''
        face_ids = []
        for face in self.faces:
            if isinstance(face, (volmdlr.faces.Face3D, surfaces.Surface3D)):
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
            if isinstance(face, volmdlr.faces.Face3D):
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
        Shell 3D rotation. Object is updated in-place.

        :param center: rotation center
        :param axis: rotation axis
        :param angle: rotation angle
        """
        warnings.warn("'inplace' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

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
        Open Shell 3D translation. Object is updated in-place.

        :param offset: Translation vector.
        :type offset: `volmdlr.Vector3D`.
        :return: Translate the Open Shell 3D in place.
        :rtype: None.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

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
        Changes frame_mapping and the object is updated in-place.

        side = 'old' or 'new'.
        """
        warnings.warn("'inplace' methods are deprecated. Use a not in-place method instead.", DeprecationWarning)

        for face in self.faces:
            face.frame_mapping_inplace(frame, side)
        new_bounding_box = self.get_bounding_box()
        self.bounding_box = new_bounding_box

    def copy(self, deep=True, memo=None):
        """
        Copy of Shell3D.

        :return: return a copy a shell 3D.
        """
        new_faces = [face.copy(deep=deep, memo=memo) for face in self.faces]
        return self.__class__(new_faces, color=self.color, alpha=self.alpha,
                              name=self.name)

    def union(self, shell2):
        """
        Combine two shells faces.

        :return: a new OpenShell3D with the combined faces.
        """
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
        """Gets the Shell bounding box."""
        bounding_boxes = []
        for face in self.faces:
            if face.outer_contour3d.primitives:
                bounding_boxes.append(face.bounding_box)
        return volmdlr.core.BoundingBox.from_bounding_boxes(bounding_boxes)

    def cut_by_plane(self, plane_3d: surfaces.Plane3D):
        """
        Cut Shell3D by plane 3d.

        :param plane_3d: plane 3d o cut shell.
        :return: return a list of faces containing the shell's sections at the plane 3d given.
        """
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
        contours3d = wires.Contour3D.contours_from_edges(
            intersection_primitives[:])
        if not contours3d:
            return []
        contours2d = [contour.to_2d(plane_3d.frame.origin,
                                    plane_3d.frame.u,
                                    plane_3d.frame.v) for contour in contours3d]
        resulting_faces = []
        for contour2d in contours2d:
            if contour2d.area() > 1e-7:
                surface2d = surfaces.Surface2D(contour2d, [])
                resulting_faces.append(volmdlr.faces.PlaneFace3D(plane_3d, surface2d))
        return resulting_faces

    def linesegment_intersections(self, linesegment3d: edges.LineSegment3D) \
            -> List[Tuple[volmdlr.faces.Face3D, List[volmdlr.Point3D]]]:
        """
        Gets the intersections of a Shell3D with a Line Segment 3D.

        :param linesegment3d: other line segment.
        :return: List of tuples (face, intersections).
        """
        intersections = []
        for face in self.faces:
            face_intersections = face.linesegment_intersections(linesegment3d)
            if face_intersections:
                intersections.append((face, face_intersections))
        return intersections

    def line_intersections(self,
                           line3d: edges.Line3D) \
            -> List[Tuple[volmdlr.faces.Face3D, List[volmdlr.Point3D]]]:
        """
        Gets the intersections of a Shell3D with a Line Segment 3D.

        :param line3d: other line segment.
        :return: List of tuples (face, intersections).
        """
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
        """
        Gets the distance between two shells.

        :param other_shell: other shell.
        :param resolution: resolution used.
        :return: return distance between faces.
        """
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
        """
        Verify if a point is on the shell (on one of the shell's faces).

        :param point: point to be verified.
        :return: return True or False.
        """
        for face in self.faces:
            if face.point_belongs(point) or face.outer_contour3d.point_over_contour(point, abs_tol=1e-7):
                return True
        return False

    def point_in_shell_face(self, point: volmdlr.Point3D):
        warnings.warn('point_in_shell_face is deprecated, please use point_on_shell instead',
                      DeprecationWarning)
        return self.point_on_shell(point)

    def triangulation(self):
        """
        Triangulation of a Shell3D.

        """
        meshes = []
        for i, face in enumerate(self.faces):
            try:
                face_mesh = face.triangulation()

            except Exception:
                warnings.warn(f"Could not triangulate {face.__class__.__name__} with index {i} in the shell "
                              f"{self.name} faces. Probabaly because topology error in contour2d.")
                print(traceback.format_exc())
                continue
            if face_mesh:
                meshes.append(face_mesh)
            else:
                warnings.warn(f"Could not triangulate {face.__class__.__name__} with index {i} in the shell "
                              f"{self.name} faces. Probabaly because topology error in contour2d.")
        return display.DisplayMesh3D.merge_meshes(meshes)

    def babylon_meshes(self, merge_meshes=True):
        """
        Returns the babylonjs mesh.
        """
        if merge_meshes:
            return super().babylon_meshes()
        babylon_meshes = []
        for face in self.faces:
            face_babylon_meshes = face.babylon_meshes()
            if not face_babylon_meshes:
                continue
            if face_babylon_meshes[0]['positions']:
                babylon_meshes.extend(face.babylon_meshes())
        babylon_mesh = {'primitives_meshes': babylon_meshes}
        babylon_mesh.update(self.babylon_param())
        return [babylon_mesh]

    def plot(self, ax=None, color: str = 'k', alpha: float = 1.0):
        """
        Plot a Shell 3D using Matplotlib.

        """
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
        points = []

        for face in self.faces:
            for _, contour in enumerate(list(chain(*[[face.outer_contour3d], face.inner_contours3d]))):
                for point_contour in contour.get_geo_points():
                    if not point_in_list(point_contour, points):
                        points.append(point_contour)

                if isinstance(contour, volmdlr.wires.Circle2D):
                    pass
                else:
                    for _, primitive in enumerate(contour.primitives):
                        if (not edge_in_list(primitive, primitives)
                                and not edge_in_list(primitive.reverse(), primitives)):
                            primitives.append(primitive)

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
                        index = get_edge_index_in_list(primitive, primitives)

                        if primitives[index].is_close(primitive):

                            if isinstance(primitive, volmdlr.edges.BSplineCurve3D):
                                discretization_points = primitive.discretization_points()

                                start_point_tag = get_point_index_in_list(discretization_points[0], points) + 1
                                end_point_tag = get_point_index_in_list(discretization_points[1], points) + 1

                                primitive_linesegments = volmdlr.edges.LineSegment3D(
                                    discretization_points[0], discretization_points[1])
                                lines.append(primitive_linesegments.get_geo_lines(tag=line_account,
                                                                                  start_point_tag=start_point_tag
                                                                                  + point_account,
                                                                                  end_point_tag=end_point_tag
                                                                                  + point_account))

                            if isinstance(primitive, volmdlr.edges.LineSegment):

                                start_point_tag = get_point_index_in_list(primitive.start, points) + 1
                                end_point_tag = get_point_index_in_list(primitive.end, points) + 1

                                lines.append(primitive.get_geo_lines(tag=line_account,
                                                                     start_point_tag=start_point_tag + point_account,
                                                                     end_point_tag=end_point_tag + point_account))
                            elif isinstance(primitive, volmdlr.edges.Arc):

                                start_point_tag = get_point_index_in_list(primitive.start, points) + 1
                                center_point_tag = get_point_index_in_list(primitive.center, points) + 1
                                end_point_tag = get_point_index_in_list(primitive.end, points) + 1

                                lines.append(primitive.get_geo_lines(tag=line_account,
                                                                     start_point_tag=start_point_tag + point_account,
                                                                     center_point_tag=center_point_tag + point_account,
                                                                     end_point_tag=end_point_tag + point_account))

                            lines_tags.append(line_account)
                            indices_check[index] = line_account
                            line_account += 1

                        if primitives[index].is_close(primitive.reverse()):

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
        """Gets Shells' mesh lines with transfinite curves."""
        lines = []
        for face in self.faces:
            lines.extend(face.surface2d.get_mesh_lines_with_transfinite_curves(
                [[face.outer_contour3d], face.inner_contours3d], min_points, size))
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

    def is_face_inside(self, face: volmdlr.faces.Face3D):
        for point in face.outer_contour3d.discretization_points(angle_resolution=1):
            point_inside_shell = self.point_belongs(point)
            if not point_inside_shell:
                point_in_shells_faces = self.point_on_shell(point)
                if not point_in_shells_faces:
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
            rays.append(edges.LineSegment3D(
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

        for ray in rays[:nb_rays]:
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

    def is_inside_shell(self, shell2):
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
        :param tol: Corresponds to the tolerance to consider two faces as intersecting faces
        :param shell2:
        :param list_coincident_faces:
        :param tol:
        :return:
        """
        # todo: delete this method if not used three months from now (25/04/2023)
        face_combinations = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                if face1.is_intersecting(face2, list_coincident_faces, tol):
                    face_combinations.append((face1, face2))
        return face_combinations

    def intersecting_faces_combinations2(self, shell2, tol=1e-8):
        """
        Gets intersecting faces combinations.

        :param shell2: ClosedShell3D
        :param tol: Corresponds to the tolerance to consider two faces as intersecting faces

        :return: returns a dictionary containing as keys the combination of intersecting faces
        and as the values the resulting primitive from the two intersecting faces.
        It is done so it is not needed to calculate the same intersecting primitive twice.
        """
        face_combinations = {}
        for face1 in self.faces:
            for face2 in shell2.faces:
                face_intersections = face1.face_intersections(face2, tol)
                if face_intersections:
                    face_combinations[(face1, face2)] = face_intersections
        return face_combinations

    @staticmethod
    def dict_intersecting_combinations(intersecting_faces_combinations, tol=1e-8):
        """
        Gets a Dictionary with the intersecting combinations.

        :param intersecting_faces_combinations: list of face combinations (list = [(face_shell1, face_shell2),...])
        for intersecting faces.
        :type intersecting_faces_combinations: list of face objects combinations
        :param tol: tolerance
        returns a dictionary containing as keys the combination of intersecting faces
        and as the values the resulting primitive from the two intersecting faces.
        It is done so it is not needed to calculate the same intersecting primitive twice.
        """
        # todo: delete this method if not used three months from now (25/04/2023)
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
                    if face.bounding_box.is_inside_bbox(shell2.bounding_box):
                        if shell2.is_face_inside(face):
                            non_intersecting_faces.append(face)

        return non_intersecting_faces

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
                    if contour1.bounding_rectangle.b_rectangle_intersection(contour2.bounding_rectangle):
                        list_coincident_faces.append((face1, face2))

        return list_coincident_faces

    def two_shells_intersecting_contour(self, shell2, dict_intersecting_combinations=None):
        """
        Computes intersecting_contour between two shells.

        :param shell2: ClosedShell3D
        :type shell2: :class:`volmdlr.faces.ClosedShell3D`.
        :param dict_intersecting_combinations: dictionary containing as keys
            the combination of intersecting faces and as the values the
            resulting primitive from the two intersecting faces
        :returns: intersecting contour for two intersecting shells
        """
        if dict_intersecting_combinations is None:
            # face_combinations = self.intersecting_faces_combinations(
            #     shell2, list_coincident_faces)
            # dict_intersecting_combinations = self.dict_intersecting_combinations(face_combinations)
            dict_intersecting_combinations = self.intersecting_faces_combinations2(shell2)
        intersecting_wires = list(dict_intersecting_combinations.values())
        intersecting_contour = \
            wires.Contour3D([wire.primitives[0] for
                             wires_ in intersecting_wires for wire in wires_])
        return intersecting_contour

    def reference_shell(self, shell2, face):
        """Reference shell used during bool operations, to help decide if a new divided face should be saved or not."""
        if face in shell2.faces:
            reference_shell = self
        else:
            reference_shell = shell2
        return reference_shell

    def set_operations_valid_exterior_faces(self, new_faces: List[volmdlr.faces.Face3D],
                                            valid_faces: List[volmdlr.faces.Face3D],
                                            list_coincident_faces: List[volmdlr.faces.Face3D],
                                            shell2, reference_shell):
        """
        Select the valid faces from the new faces created during Boolean operations.

        :param new_faces: list of new divided faces.
        :param valid_faces: list of already validated faces.
        :param list_coincident_faces: if of coincident faces.
        :param shell2: shell2, used in the Boolean operation.
        :param reference_shell: reference shell, to help decide if a new divided face should be saved or not.
        :return:
        """
        for new_face in new_faces:
            if self.set_operations_exterior_face(new_face, valid_faces, reference_shell,
                                                 list_coincident_faces, shell2):
                valid_faces.append(new_face)
        return valid_faces

    def union_faces(self, shell2, intersecting_faces, intersecting_combinations, list_coincident_faces):
        """
        Gets new faces for union Boolean operation between two closed shell 3d.

        :param shell2: other shell
        :param intersecting_faces: list of all intersecting faces.
        :param intersecting_combinations: Dictionary containing all combination of faces intersection,\
        with corresponding intersections.
        :param list_coincident_faces: list of coincident faces.
        :return: list of new faces for union of two closed shell3.
        """
        faces = []
        for face in intersecting_faces:
            reference_shell = self.reference_shell(shell2, face)
            new_faces = face.set_operations_new_faces(intersecting_combinations)
            faces = self.set_operations_valid_exterior_faces(new_faces, faces, list_coincident_faces,
                                                             shell2, reference_shell)
        if list_coincident_faces:
            faces = self.validate_set_operations_faces(faces)
        return faces

    def get_subtraction_valid_faces(self, new_faces, valid_faces, reference_shell, shell2, keep_interior_faces):
        """
        Gets valid faces for subtraction Boolean operations.

        :param new_faces: list of new divided faces.
        :param valid_faces: list of already validated faces.
        :param reference_shell: reference shell, to help decide if a new divided face should be saved or not.
        :param shell2: other shell.
        :param keep_interior_faces: Boolean to decide to keep interior faces on reference shell or not.
        :return: return a list a valid faces.
        """
        faces = []
        for new_face in new_faces:
            if keep_interior_faces:
                if self.set_operations_interior_face(new_face, valid_faces, reference_shell):
                    faces.append(new_face)
            elif self.set_operations_exterior_face(new_face, faces, reference_shell, [], shell2):
                faces.append(new_face)
        return faces

    @staticmethod
    def validate_set_operations_faces(faces):
        """
        Final validation of new faces created during intersections or subtractions of two closed shells.

        :param faces: new faces.
        :return: valid faces.
        """
        valid_faces = []
        while True:
            if not faces:
                break
            for face in valid_faces:
                if face.face_inside(faces[0]):
                    faces.remove(faces[0])
                    break
            else:
                valid_faces.append(faces[0])
                faces.remove(faces[0])
        return valid_faces

    def subtraction_faces(self, shell2, intersecting_faces, intersecting_combinations):
        """
        Gets new faces for subtraction Boolean operation between two closed shell 3d.

        :param shell2: other shell
        :param intersecting_faces: list of all intersecting faces.
        :param intersecting_combinations: Dictionary containing all combination of faces intersection,\
        with corresponding intersections.
        :return: list of new faces for subtraction of two closed shells 3.
        """
        faces = []
        for face in intersecting_faces:
            keep_interior_faces = False
            if face in shell2.faces:
                keep_interior_faces = True
            reference_shell = self.reference_shell(shell2, face)
            new_faces = face.set_operations_new_faces(intersecting_combinations)
            valid_faces = self.get_subtraction_valid_faces(new_faces, faces, reference_shell,
                                                           shell2, keep_interior_faces)
            faces.extend(valid_faces)

        valid_faces = self.validate_set_operations_faces(faces)

        return valid_faces

    def valid_intersection_faces(self, new_faces, valid_faces, reference_shell, shell2):
        """
        Validate Boolean intersection operation new faces.

        :param new_faces: list of new divided faces.
        :param valid_faces: list of already validated faces.
        :param reference_shell: reference shell, to help decide if a new divided face should be saved or not.
        :param shell2: other shell.
        :return:
        """
        faces = []
        for new_face in new_faces:
            inside_reference_shell = reference_shell.point_belongs(
                new_face.random_point_inside())
            if (inside_reference_shell or (self.face_on_shell(new_face) and shell2.face_on_shell(new_face))) \
                    and new_face not in valid_faces:
                faces.append(new_face)

        return faces

    def intersection_faces(self, shell2, intersecting_faces, intersecting_combinations):
        """
        Gets new faces for intersection Boolean operation between two closed shell 3d.

        :param shell2: other shell
        :param intersecting_faces: list of all intersecting faces.
        :param intersecting_combinations: Dictionary containing all combination of faces intersection,\
        with corresponding intersections.
        :return: list of new faces for intersection of two closed shells 3d.
        """
        faces = []
        for face in intersecting_faces:
            reference_shell = self.reference_shell(shell2, face)
            new_faces = face.set_operations_new_faces(intersecting_combinations)
            valid_faces = self.valid_intersection_faces(
                new_faces, faces, reference_shell, shell2)
            faces.extend(valid_faces)

        valid_faces = self.validate_set_operations_faces(faces)
        return valid_faces

    def set_operations_interior_face(self, new_face, faces, reference_shell):
        """
        Verify if new face is inside reference shell for Boolean operations.

        :param new_face: new divided face.
        :param faces: list of already validated faces.
        :param reference_shell: reference shell, to help decide if a new divided face should be saved or not.
        """
        inside_reference_shell = reference_shell.point_belongs(new_face.random_point_inside())
        if inside_reference_shell and new_face not in faces:
            return True
        if self.face_on_shell(new_face):
            return True
        return False

    def is_face_between_shells(self, shell2, face):
        """
        Verify if face is between the two shells.

        :param shell2: other shell
        :param face: face to be verified.
        :return:
        """
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

    def set_operations_exterior_face(self, new_face, valid_faces, reference_shell,
                                     list_coincident_faces, shell2):
        """
        Selects exterior faces during bool operations, like union or subtraction.

        :param new_face: divided faces.
        :param valid_faces: list of already validated faces.
        :param reference_shell: reference shell, to help decide if a new divided face should be saved or not.
        :param shell2: other shell.
        :param list_coincident_faces: list of coincident faces.
        :return:
        """
        if new_face.area() < 1e-8:
            return False
        if new_face not in valid_faces:
            inside_reference_shell = reference_shell.point_belongs(new_face.random_point_inside())
            face_on_reference_shell = reference_shell.face_on_shell(new_face)
            if not inside_reference_shell or face_on_reference_shell:
                if list_coincident_faces:
                    if self.is_face_between_shells(shell2, new_face):
                        return False
                return True
        return False

    def validate_set_operation(self, shell2, tol):
        """
        Verifies if two shells are valid for union or subtractions operations.

        Its Verifies if they are disjointed or if one is totally inside the other.

        If it returns an empty list, it means the two shells are valid to continue the
        operation.
        """
        if self.is_disjoint_from(shell2, tol):
            return [self, shell2]
        if self.is_inside_shell(shell2):
            return [shell2]
        if shell2.is_inside_shell(self):
            return [self]
        return []

    def union(self, shell2: 'ClosedShell3D', tol: float = 1e-8):
        """
        Given Two closed shells, it returns a new united ClosedShell3D object.

        """

        validate_set_operation = self.validate_set_operation(shell2, tol)
        if validate_set_operation:
            return validate_set_operation
        list_coincident_faces = self.get_coincident_faces(shell2)
        # face_combinations = self.intersecting_faces_combinations(shell2, list_coincident_faces, tol)
        # intersecting_combinations = self.dict_intersecting_combinations(face_combinations, tol)
        intersecting_combinations = self.intersecting_faces_combinations2(shell2, tol)
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
        """Gets faces that are adjacent, and sharing the same surface, so they can be merged."""
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
                list_new_faces += volmdlr.faces.PlaneFace3D.merge_faces(valid_coicident_faces)
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

        # list_coincident_faces = self.get_coincident_faces(shell2)
        # face_combinations = self.intersecting_faces_combinations(
        #     shell2, list_coincident_faces, tol)
        #
        # intersecting_combinations = self.dict_intersecting_combinations(
        #     face_combinations, tol)

        list_coincident_faces = self.get_coincident_faces(shell2)
        intersecting_combinations = self.intersecting_faces_combinations2(shell2, tol)

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

        # list_coincident_faces = self.get_coincident_faces(shell2)
        # face_combinations = self.intersecting_faces_combinations(
        #     shell2, list_coincident_faces, tol)
        # intersecting_combinations = self.dict_intersecting_combinations(
        #     face_combinations, tol)
        # list_coincident_faces = self.get_coincident_faces(shell2)
        intersecting_combinations = self.intersecting_faces_combinations2(shell2, tol)

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

    def validate_intersection_operation(self, shell2, tol):
        """
        Verifies if two shells are valid for union or subtractions operations.

        Its Verifies if they are disjointed or if one is totally inside the other.
        If it returns an empty list, it means the two shells are valid to continue the
        operation.
        """
        if self.is_inside_shell(shell2):
            return [self]
        if shell2.is_inside_shell(self):
            return [shell2]
        return []

    def intersection(self, shell2, tol=1e-8):
        """
        Given two ClosedShell3D, it returns the new object resulting from the intersection of the two.

        """
        if self.is_disjoint_from(shell2, tol):
            return []
        validate_set_operation = self.validate_intersection_operation(
            shell2, tol)
        if validate_set_operation:
            return validate_set_operation
        # list_coincident_faces = self.get_coincident_faces(shell2)
        # face_combinations = self.intersecting_faces_combinations(shell2, list_coincident_faces, tol)
        # intersecting_combinations = self.dict_intersecting_combinations(face_combinations, tol)
        intersecting_combinations = self.intersecting_faces_combinations2(shell2, tol)

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

    def __init__(self, faces: List[volmdlr.faces.Triangle3D],
                 color: Tuple[float, float, float] = None,
                 alpha: float = 1., name: str = ''):
        OpenShell3D.__init__(self, faces=faces, color=color, alpha=alpha, name=name)

    def to_dict(self, *args, **kwargs):
        dict_ = self.base_dict()
        dict_['faces'] = [t.to_dict() for t in self.faces]
        dict_['alpha'] = self.alpha
        dict_['color'] = self.color
        return dict_

    def to_mesh_data(self):
        """To mesh data for Open Triangle Shell."""
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
        """Creates an Open Triangle Shell 3D from mesh data."""
        triangles = []
        points = [volmdlr.Point3D(px, py, pz) for px, py, pz in positions]
        for i1, i2, i3 in faces:
            triangles.append(volmdlr.faces.Triangle3D(points[i1], points[i2], points[i3]))
        return cls(triangles)

    def to_trimesh(self):
        """Creates a Trimesh from Open Triangle Shell 3D."""
        return Trimesh(*self.to_mesh_data())

    @classmethod
    def from_trimesh(cls, trimesh):
        """Creates an Open Triangle Shell 3D from Trimesh."""
        return cls.from_mesh_data(trimesh.vertices.tolist(), trimesh.faces.tolist())

    def triangulation(self):
        """Triangulation of an Open Triangle Shell 3D."""
        points = []
        triangles = []
        for i, triangle in enumerate(self.faces):
            points.append(display.Node3D.from_point(triangle.point1))
            points.append(display.Node3D.from_point(triangle.point2))
            points.append(display.Node3D.from_point(triangle.point3))
            triangles.append((3 * i, 3 * i + 1, 3 * i + 2))
        return display.DisplayMesh3D(points, triangles)


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

    def __init__(self, faces: List[volmdlr.faces.Triangle3D],
                 color: Tuple[float, float, float] = None,
                 alpha: float = 1., name: str = ''):
        OpenTriangleShell3D.__init__(self, faces=faces, color=color, alpha=alpha, name=name)
        ClosedShell3D.__init__(self, faces)
