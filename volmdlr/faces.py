"""
Surfaces & faces.

"""

import math
import warnings
from itertools import chain, product
from typing import List
import matplotlib.pyplot as plt
import numpy as npy

from dessia_common.core import DessiaObject

import volmdlr.core
from volmdlr.core import EdgeStyle
import volmdlr.core_compiled
import volmdlr.display as vmd
import volmdlr.edges as vme
import volmdlr.curves as volmdlr_curves
import volmdlr.geometry
import volmdlr.grid
from volmdlr import surfaces
from volmdlr.utils.parametric import array_range_search
import volmdlr.wires

warnings.simplefilter("once")


def octree_decomposition(bbox, faces):
    """Decomposes a list of faces into eight Bounding boxes subdivided boxes."""
    decomposition = {octant: [] for octant in bbox.octree()}
    for face in faces:
        center = face.bounding_box.center
        for octant in bbox.octree():
            if octant.point_belongs(center):
                decomposition[octant].append(face)
                break
    decomposed = {octant: faces for octant, faces in decomposition.items() if faces}
    return decomposed


def octree_face_decomposition(face):
    """
    Decomposes the face discretization triangle faces inside eight boxes from a bounding box octree structure.

    :param face: given face.
    :return: returns a dictionary containing bounding boxes as keys and as values, a list of faces
    inside that bounding box.
    """
    triangulation = face.triangulation()
    triangulation_faces = triangulation.faces
    return octree_decomposition(face.bounding_box, triangulation_faces)


class Face3D(volmdlr.core.Primitive3D):
    """
    Abstract method to define 3D faces.
    """

    min_x_density = 1
    min_y_density = 1

    def __init__(self, surface3d, surface2d: surfaces.Surface2D,
                 name: str = ''):
        self.surface3d = surface3d
        self.surface2d = surface2d
        self._outer_contour3d = None
        self._inner_contours3d = None
        self._face_octree_decomposition = None
        # self.bounding_box = self._bounding_box()

        volmdlr.core.Primitive3D.__init__(self, name=name)

    def to_dict(self, *args, **kwargs):
        """Avoids storing points in memo that makes serialization slow."""
        return DessiaObject.to_dict(self, use_pointers=False)

    def __hash__(self):
        """Computes the hash."""
        return hash(self.surface3d) + hash(self.surface2d)

    def __eq__(self, other_):
        """Computes the equality to another face."""
        if other_.__class__.__name__ != self.__class__.__name__:
            return False
        equal = (self.surface3d == other_.surface3d
                 and self.surface2d == other_.surface2d)
        return equal

    def point_belongs(self, point3d: volmdlr.Point3D, tol: float = 1e-6):
        """
        Tells you if a point is on the 3D face and inside its contour.
        """
        if not self.bounding_box.point_belongs(point3d):
            return False
        point2d = self.surface3d.point3d_to_2d(point3d)
        # check_point3d = self.surface3d.point2d_to_3d(point2d)
        # if check_point3d.point_distance(point3d) > tol:
        if not self.surface3d.point_on_surface(point3d, tol):
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

    @outer_contour3d.setter
    def outer_contour3d(self, contour3d):
        self._outer_contour3d = contour3d

    @property
    def inner_contours3d(self):
        """
        Gives the 3d version of the inner contours of the face.
        """
        if not self._inner_contours3d:
            self._inner_contours3d = [self.surface3d.contour2d_to_3d(c) for c in
                                      self.surface2d.inner_contours]
        return self._inner_contours3d

    @inner_contours3d.setter
    def inner_contours3d(self, contours3d):
        self._inner_contours3d = contours3d

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
        """General method to get the bounding box of a face 3D."""
        return self.outer_contour3d.bounding_box

    def area(self):
        """Computes the area of the surface2d."""
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
        face = globals()[surface.face_class]
        point_in_contours3d = any(isinstance(contour, volmdlr.Point3D) for contour in contours)
        if (len(contours) == 1) and isinstance(contours[0], volmdlr.Point3D):
            return face.from_surface_rectangular_cut(surface)
        if len(contours) == 2 and point_in_contours3d:
            vertex = next(contour for contour in contours if isinstance(contour, volmdlr.Point3D))
            base = next(contour for contour in contours if contour is not vertex)
            return face.from_base_and_vertex(surface, base, vertex, name)
        if point_in_contours3d:
            point = next(contour for contour in contours if isinstance(contour, volmdlr.Point3D))
            contours = [contour for contour in contours if contour is not point]
            return face.from_contours3d_and_rectangular_cut(surface, contours, point)
        return face.from_contours3d(surface, contours, name)

    @classmethod
    def from_contours3d(cls, surface, contours3d: List[volmdlr.wires.Contour3D], name: str = ''):
        """
        Returns the face generated by a list of contours. Finds out which are outer or inner contours.

        :param surface: Surface3D where the face is defined.
        :param contours3d: List of 3D contours representing the face's BREP.
        :param name: the name to inject in the new face
        """
        outer_contour2d = None
        outer_contour3d, inner_contours3d = None, []
        if len(contours3d) == 1:
            outer_contour2d = surface.contour3d_to_2d(contours3d[0])
            outer_contour3d = surface.contour2d_to_3d(outer_contour2d)
            inner_contours2d = []

        elif len(contours3d) > 1:
            area = -1
            inner_contours2d = []
            inner_contours3d = []

            contours2d = [surface.contour3d_to_2d(contour3d) for contour3d in contours3d]

            check_contours = [not contour2d.is_ordered(tol=1e-2) for contour2d in contours2d]
            if any(check_contours):
                # Not implemented yet, but connect_contours should also return outer_contour3d and inner_contours3d
                outer_contour2d, inner_contours2d = surface.connect_contours(contours2d[0], contours2d[1:])
                outer_contour3d = surface.contour2d_to_3d(outer_contour2d)
                inner_contours3d = [surface.contour2d_to_3d(contour) for contour in inner_contours2d]
            else:
                for contour2d, contour3d in zip(contours2d, contours3d):
                    # if not contour2d.is_ordered(1e-4):
                    #     contour2d = vm_parametric.contour2d_healing(contour2d)
                    inner_contours2d.append(contour2d)
                    inner_contours3d.append(contour3d)
                    contour_area = contour2d.area()
                    if contour_area > area:
                        area = contour_area
                        outer_contour2d = contour2d
                        outer_contour3d = contour3d
                inner_contours2d.remove(outer_contour2d)
                inner_contours3d.remove(outer_contour3d)
        else:
            raise ValueError('Must have at least one contour')
        if (not outer_contour2d) or (not all(outer_contour2d.primitives)) or (not outer_contour2d.is_ordered(1e-2)):
            return None
        # if outer_contour3d and outer_contour3d.primitives and not outer_contour3d.is_ordered(1e-5):
        #     outer_contour2d = contour2d_healing(outer_contour2d)
        face = cls(surface,
                   surface2d=surfaces.Surface2D(outer_contour=outer_contour2d, inner_contours=inner_contours2d),
                   name=name)
        # To improve performance while reading from step file
        face.outer_contour3d = outer_contour3d
        face.inner_contours3d = inner_contours3d
        return face

    def to_step(self, current_id):
        content, surface3d_ids = self.surface3d.to_step(current_id)
        current_id = max(surface3d_ids)

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
        """
        Specifies the number of subdivision when using triangulation by lines. (Old triangulation).
        """
        return [], []

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        return [0, 0]

    def triangulation(self, grid_size=None):
        if not grid_size:
            number_points_x, number_points_y = self.grid_size()
        else:
            number_points_x, number_points_y = grid_size
        mesh2d = self.surface2d.triangulation(number_points_x, number_points_y)
        if mesh2d is None:
            return None
        return vmd.DisplayMesh3D([vmd.Node3D(*self.surface3d.point2d_to_3d(point)) for point in mesh2d.points],
                                 mesh2d.triangles)

    def plot2d(self, ax=None, color='k', alpha=1):
        """Plot 2D of the face using matplotlib."""
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

    def translation(self, offset: volmdlr.Vector3D):
        """
        Face3D translation.

        :param offset: Translation vector.
        :type offset: `volmdlr.Vector3D`
        :return: A new translated Face3D
        """
        new_surface3d = self.surface3d.translation(offset=offset)
        return self.__class__(new_surface3d, self.surface2d)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Face3D.

        side = 'old' or 'new'
        """
        new_surface3d = self.surface3d.frame_mapping(frame, side)
        return self.__class__(new_surface3d, self.surface2d.copy(),
                              self.name)

    def copy(self, deep=True, memo=None):
        """Returns a copy of the Face3D."""
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
                if self_contour2d.is_inside(face2_contour2d):
                    for inner_contour in self.inner_contours3d:
                        inner_contour2d = inner_contour.to_2d(
                            self.surface3d.frame.origin, self.surface3d.frame.u, self.surface3d.frame.v)
                        if inner_contour2d.is_inside(face2_contour2d) or inner_contour2d.is_superposing(
                                face2_contour2d):
                            return False
                return True
            if self_contour2d.is_superposing(face2_contour2d):
                return True
        return False

    def edge_intersections(self, edge):
        """Gets the intersections of an edge and a 3D face."""
        intersections = []
        method_name = f'{edge.__class__.__name__.lower()[:-2]}_intersections'
        if hasattr(self, method_name):
            intersections = getattr(self, method_name)(edge)
        elif hasattr(self.surface3d, method_name):
            edge_surface_intersections = getattr(self.surface3d, method_name)(edge)
            for intersection in edge_surface_intersections:
                if self.point_belongs(intersection) and not volmdlr.core.point_in_list(intersection, intersections):
                    intersections.append(intersection)
        if not intersections:
            for point in [edge.start, edge.end]:
                if self.point_belongs(point):
                    if point not in intersections:
                        intersections.append(point)
        return intersections

    def line_intersections(self, line: volmdlr_curves.Line3D) -> List[volmdlr.Point3D]:
        """
        Get intersections between a face 3d and a Line 3D.

        :param line: other line.
        :return: a list of intersections.
        """
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
        """
        Get intersections between a face 3d and a Line Segment 3D.

        :param linesegment: other linesegment.
        :return: a list of intersections.
        """
        linesegment_intersections = []
        if not self.bounding_box.is_intersecting(linesegment.bounding_box):
            return []
        for intersection in self.surface3d.linesegment_intersections(linesegment):
            if self.point_belongs(intersection):
                linesegment_intersections.append(intersection)
        return linesegment_intersections

    def fullarc_intersections(self, fullarc: vme.FullArc3D) -> List[volmdlr.Point3D]:
        """
        Get intersections between a face 3d and a Full Arc 3D.

        :param fullarc: other fullarc.
        :return: a list of intersections.
        """
        intersections = []
        for intersection in self.surface3d.fullarc_intersections(fullarc):
            if self.point_belongs(intersection):
                intersections.append(intersection)
        return intersections

    def plot(self, ax=None, color='k', alpha=1, edge_details=False):
        """Plots the face."""
        if not ax:
            _, ax = plt.subplots(subplot_kw={"projection": "3d"})

        self.outer_contour3d.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha,
                                                              edge_ends=edge_details, edge_direction=edge_details))
        for contour3d in self.inner_contours3d:
            contour3d.plot(ax=ax, edge_style=EdgeStyle(color=color, alpha=alpha,
                                                       edge_ends=edge_details, edge_direction=edge_details))
        return ax

    def random_point_inside(self):
        """Gets a random point on the face."""
        point_inside2d = self.surface2d.random_point_inside()
        return self.surface3d.point2d_to_3d(point_inside2d)

    def is_adjacent(self, face2: 'Face3D'):
        """
        Verifies if two faces are adjacent or not.

        :param face2: other face.
        :return: True or False.
        """
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

            if isinstance(contour, volmdlr_curves.Circle2D):
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

        with open(file_name + '.geo', 'w', encoding="utf-8") as file:
            for line in lines:
                file.write(line)
                file.write('\n')
        file.close()

    def get_geo_lines(self, tag: int, line_loop_tag: List[int]):
        """
        Gets the lines that define a PlaneFace3D in a .geo file.

        """

        return 'Plane Surface(' + str(tag) + ') = {' + str(line_loop_tag)[1:-1] + '};'

    def edge3d_inside(self, edge3d, abs_tol: float = 1e-6):
        """
        Returns True if edge 3d is coplanar to the face.
        """
        method_name = f'{edge3d.__class__.__name__.lower()[:-2]}_inside'
        if hasattr(self, method_name):
            return getattr(self, method_name)(edge3d)
        points = edge3d.discretization_points(number_points=10)
        for point in points[1:-1]:
            if not self.point_belongs(point, abs_tol):
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
        if self.bounding_box.is_intersecting(face2.bounding_box, tol) and \
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
                    if not volmdlr.core.point_in_list(point, intersections_points):
                        intersections_points.append(point)

        return intersections_points

    def face_intersections(self, face2, tol=1e-6) -> List[volmdlr.wires.Wire3D]:
        """
        Calculates the intersections between two Face3D.

        """

        bbox1 = self.bounding_box
        bbox2 = face2.bounding_box
        if not bbox1.is_intersecting(bbox2, tol):
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

    def set_operations_new_faces(self, intersecting_combinations):
        self_copy = self.copy(deep=True)
        list_cutting_contours = self_copy.get_face_cutting_contours(intersecting_combinations)
        if not list_cutting_contours:
            return [self_copy]
        return self_copy.divide_face(list_cutting_contours)

    def split_inner_contour_intersecting_cutting_contours(self, list_cutting_contours):
        """
        Given a list contours cutting the face, it calculates inner contours intersections with these contours.

        Then, these inner contours were split at the found intersecting points.
        :param list_cutting_contours: list of contours cutting face.
        :return:
        """
        list_split_inner_contours = []
        for inner_contour in self.surface2d.inner_contours:
            list_intersecting_points_with_inner_contour = []
            for cutting_contour in list_cutting_contours:
                contour_intersection_points = inner_contour.intersection_points(cutting_contour)
                if not contour_intersection_points:
                    continue
                list_intersecting_points_with_inner_contour.extend(contour_intersection_points)
            inner_contour_intersections_with_outer_contour = inner_contour.intersection_points(
                self.surface2d.outer_contour)
            if list_intersecting_points_with_inner_contour and inner_contour_intersections_with_outer_contour:
                list_intersecting_points_with_inner_contour.extend(inner_contour_intersections_with_outer_contour)
            sorted_intersections_points_along_inner_contour = inner_contour.sort_points_along_wire(
                list_intersecting_points_with_inner_contour)
            if sorted_intersections_points_along_inner_contour:
                list_split_inner_contours.extend(inner_contour.split_with_sorted_points(
                    sorted_intersections_points_along_inner_contour))
        return list_split_inner_contours

    def get_face_cutting_contours(self, dict_intersecting_combinations):
        """
        Get all contours cutting the face, resulting from multiple faces intersections.

        :param dict_intersecting_combinations: dictionary containing as keys the combination of intersecting faces
        and as the values the resulting primitive from the intersection of these two faces
        return a list all contours cutting one particular face.
        """
        face_intersecting_primitives2d = self.select_face_intersecting_primitives(dict_intersecting_combinations)
        if not face_intersecting_primitives2d:
            return []
        list_cutting_contours = volmdlr.wires.Contour2D.contours_from_edges(face_intersecting_primitives2d[:])

        if not self.surface2d.inner_contours:
            return list_cutting_contours
        list_split_inner_contours = self.split_inner_contour_intersecting_cutting_contours(list_cutting_contours)

        valid_cutting_contours = []

        # remove split_inner_contour connected to a cutting_contour at two points.
        connected_at_two_ends = []
        for cutting_contour in list_cutting_contours:
            for split_contour in list_split_inner_contours:
                if split_contour.point_over_wire(cutting_contour.primitives[0].start) and \
                        split_contour.point_over_wire(cutting_contour.primitives[-1].end):
                    connected_at_two_ends.append(split_contour)
                    break
        list_split_inner_contours = [split_contour for split_contour in list_split_inner_contours
                                     if split_contour not in connected_at_two_ends]
        while list_cutting_contours:
            for i, cutting_contour in enumerate(list_cutting_contours[:]):
                if (self.surface2d.outer_contour.point_over_contour(cutting_contour.primitives[0].start) and
                    self.surface2d.outer_contour.point_over_contour(cutting_contour.primitives[-1].end)) or \
                        cutting_contour.primitives[0].start.is_close(cutting_contour.primitives[-1].end):
                    valid_cutting_contours.append(cutting_contour)
                    list_cutting_contours.pop(i)
                    break
                list_cutting_contours.pop(i)
                while True:
                    connecting_split_contour = cutting_contour.get_connected_wire(list_split_inner_contours)
                    list_split_inner_contours.remove(connecting_split_contour)
                    new_contour = volmdlr.wires.Contour2D.contours_from_edges(
                        cutting_contour.primitives + connecting_split_contour.primitives)[0]

                    if self.surface2d.outer_contour.are_extremity_points_touching(new_contour) or \
                            new_contour.is_contour_closed():
                        valid_cutting_contours.append(new_contour)
                        break

                    connecting_cutting_contour = new_contour.get_connected_wire(list_cutting_contours)
                    if not connecting_cutting_contour:
                        if any(self.surface2d.outer_contour.point_over_contour(point)
                               for point in [new_contour.primitives[0].start, new_contour.primitives[-1].end]) and \
                                any(valid_contour.point_over_contour(point) for valid_contour in
                                    valid_cutting_contours for
                                    point in [new_contour.primitives[0].start, new_contour.primitives[-1].end]):
                            valid_cutting_contours.append(new_contour)
                        break

                    new_contour = volmdlr.wires.Contour2D.contours_from_edges(
                        new_contour.primitives + connecting_cutting_contour.primitives)[0]
                    list_cutting_contours.remove(connecting_cutting_contour)

                    if self.surface2d.outer_contour.are_extremity_points_touching(new_contour):
                        valid_cutting_contours.append(new_contour)
                        break

                    cutting_contour = new_contour
                break
        return valid_cutting_contours + list_split_inner_contours

    def divide_face(self, list_cutting_contours: List[volmdlr.wires.Contour2D]):
        """
        Divides a Face 3D with a list of cutting contours.

        :param list_cutting_contours: list of contours cutting the face.
        """
        list_faces = []
        list_open_cutting_contours = []
        list_closed_cutting_contours = []
        face_inner_contours = self.surface2d.inner_contours[:]
        list_cutting_contours_ = []
        while list_cutting_contours:
            cutting_contour = list_cutting_contours[0]
            if not cutting_contour.primitives[0].start.is_close(cutting_contour.primitives[-1].end):
                list_cutting_contours_.append(cutting_contour)
                list_cutting_contours.remove(cutting_contour)
                continue
            for inner_contour in face_inner_contours:
                if cutting_contour.is_inside(inner_contour):
                    if cutting_contour.is_sharing_primitives_with(inner_contour):
                        merged_contours = cutting_contour.merge_with(inner_contour)
                        list_cutting_contours.remove(cutting_contour)
                        cutting_contour = merged_contours[0]
                    # list_cutting_contours = merged_contours + list_cutting_contours
                    # break
                    continue
                if cutting_contour.is_sharing_primitives_with(inner_contour):
                    merged_contours = cutting_contour.merge_with(inner_contour)
                    list_cutting_contours.remove(cutting_contour)
                    face_inner_contours.remove(inner_contour)
                    list_cutting_contours = merged_contours + list_cutting_contours
                    break
            else:
                list_cutting_contours_.append(cutting_contour)
                if cutting_contour in list_cutting_contours:
                    list_cutting_contours.remove(cutting_contour)

        list_cutting_contours = list_cutting_contours_
        list_cutting_contours_ = []
        for cutting_contour in list_cutting_contours:
            contour_intersections = self.surface2d.outer_contour.wire_intersections(cutting_contour)
            if len(contour_intersections) >= 2:
                contour_intersections = sorted(contour_intersections, key=cutting_contour.abscissa)
                list_cutting_contours_.extend(cutting_contour.split_with_sorted_points(contour_intersections))
                continue
            list_cutting_contours_.append(cutting_contour)
        list_cutting_contours = list_cutting_contours_
        for cutting_contour in list_cutting_contours:
            if not cutting_contour.primitives[0].start.is_close(cutting_contour.primitives[-1].end):
                list_open_cutting_contours.append(cutting_contour)
                continue
            list_closed_cutting_contours.append(cutting_contour)
        if list_open_cutting_contours:
            list_faces = self.divide_face_with_open_cutting_contours(list_open_cutting_contours)
        list_faces = self.divide_face_with_closed_cutting_contours(list_closed_cutting_contours, list_faces)
        list_faces = [face for face in list_faces if not math.isclose(face.area(), 0.0, abs_tol=1e-08)]
        return list_faces

    def divide_face_with_open_cutting_contours(self, list_open_cutting_contours):
        """
        Divides a face 3D with a list of closed cutting contour, that is, it will cut holes on the face.

        :param list_open_cutting_contours: list containing the open cutting contours.
        :return: list divided faces.
        """
        list_faces = []
        if not self.surface2d.outer_contour.edge_polygon.is_trigo:
            self.surface2d.outer_contour = self.surface2d.outer_contour.invert()
        new_faces_contours = self.surface2d.outer_contour.divide(list_open_cutting_contours)
        new_inner_contours = len(new_faces_contours) * [[]]
        if self.surface2d.inner_contours:
            new_faces_contours, new_inner_contours = self.get_open_contour_divided_faces_inner_contours(
                new_faces_contours)
        if isinstance(self, Triangle3D):
            class_to_instanciate = PlaneFace3D
        else:
            class_to_instanciate = self.__class__
        for contour, inner_contours in zip(new_faces_contours, new_inner_contours):
            new_face = class_to_instanciate(self.surface3d, surfaces.Surface2D(contour, inner_contours))
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
        for closed_cutting_contour in list_closed_cutting_contours:
            if closed_cutting_contour.primitives[0].start.is_close(closed_cutting_contour.primitives[-1].end):
                inner_contours1 = []
                inner_contours2 = []
                if list_faces:
                    new_list_faces = self.get_closed_contour_divided_faces_inner_contours(
                        list_faces, closed_cutting_contour)
                    list_faces = list_faces + new_list_faces
                    continue
                new_contour_adjacent_to_inner_contour = False
                for inner_contour in self.surface2d.inner_contours:
                    if closed_cutting_contour.is_inside(inner_contour):
                        inner_contours2.append(inner_contour)
                        continue
                    if closed_cutting_contour.is_sharing_primitives_with(inner_contour):
                        new_contour_adjacent_to_inner_contour = True
                        inner_contours1.extend(closed_cutting_contour.merge_with(inner_contour))
                    else:
                        inner_contours1.append(inner_contour)
                if not new_contour_adjacent_to_inner_contour:
                    inner_contours1.append(closed_cutting_contour)
                if isinstance(self, Triangle3D):
                    class_to_instanciate = PlaneFace3D
                else:
                    class_to_instanciate = self.__class__
                surf3d = self.surface3d
                surf2d = surfaces.Surface2D(self.surface2d.outer_contour, inner_contours1)
                new_plane = class_to_instanciate(surf3d, surf2d)
                list_faces.append(new_plane)
                list_faces.append(class_to_instanciate(surf3d, surfaces.Surface2D(
                    closed_cutting_contour, inner_contours2)))
                continue
            surf3d = self.surface3d
            surf2d = surfaces.Surface2D(self.surface2d.outer_contour, [])
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
        new_faces_contours_ = []
        for new_contour in new_faces_contours:
            for inner_contour in self.surface2d.inner_contours:
                if new_contour.is_superposing(inner_contour):
                    break
            else:
                new_faces_contours_.append(new_contour)
        new_faces_contours = new_faces_contours_
        while new_faces_contours:
            new_face_contour = new_faces_contours[0]
            if new_face_contour in valid_new_faces_contours:
                new_faces_contours.remove(new_face_contour)
                continue
            inner_contours = []
            for inner_contour in self.surface2d.inner_contours:
                if not new_face_contour.is_inside(inner_contour):
                    continue
                if new_face_contour.is_sharing_primitives_with(inner_contour):
                    merged_new_face_contours = new_face_contour.merge_with(inner_contour)
                    if merged_new_face_contours:
                        new_faces_contours.remove(new_face_contour)
                        new_faces_contours = merged_new_face_contours + new_faces_contours
                        break
                else:
                    inner_contours.append(inner_contour)
            else:
                valid_new_faces_contours.append(new_face_contour)
                valid_inner_contours.append(inner_contours)
                new_faces_contours.remove(new_face_contour)
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
                            new_list_faces.append(self.__class__(surf3d, surfaces.Surface2D(new_contour, [])))
                            break
                        surf2d = surfaces.Surface2D(new_face.surface2d.outer_contour, inner_contours1)
                        new_plane = self.__class__(surf3d, surf2d)
                        new_list_faces.append(new_plane)
                    if inner_contours2:
                        new_list_faces.append(
                            self.__class__(surf3d, surfaces.Surface2D(new_contour, inner_contours2)))
        return new_list_faces

    def select_face_intersecting_primitives(self, dict_intersecting_combinations):
        """
        Select face intersecting primitives from a dictionary containing all intersection combinations.

        :param dict_intersecting_combinations: dictionary containing all intersection combinations
        :return: list of intersecting primitives for current face
        """
        face_intersecting_primitives2d = []
        for intersecting_combination, intersections in dict_intersecting_combinations.items():
            if self in intersecting_combination:
                for intersection_wire in intersections:
                    if len(intersection_wire.primitives) != 1:
                        raise NotImplementedError
                    primitive2_2d = self.surface3d.contour3d_to_2d(intersection_wire).primitives[0]
                    if volmdlr.core.edge_in_list(primitive2_2d, face_intersecting_primitives2d) or \
                            volmdlr.core.edge_in_list(primitive2_2d.reverse(), face_intersecting_primitives2d):
                        continue
                    if not self.surface2d.outer_contour.primitive_over_contour(primitive2_2d, tol=1e-7):
                        face_intersecting_primitives2d.append(primitive2_2d)
        return face_intersecting_primitives2d

    def _is_linesegment_intersection_possible(self, linesegment: vme.LineSegment3D):
        """
        Verifies if intersection of face with line segment is possible or not.

        :param linesegment: other line segment.
        :return: returns True if possible, False otherwise.
        """
        if not self.bounding_box.is_intersecting(linesegment.bounding_box):
            return False
        if math.isclose(self.area(), 0.0, abs_tol=1e-10):
            return False
        bbox_block_faces = volmdlr.primitives3d.Block.from_bounding_box(self.bounding_box).faces
        if not any(bbox_face.line_intersections(linesegment.line) for bbox_face in bbox_block_faces):
            return False
        return True

    def _get_linesegment_intersections_approximation(self, linesegment: vme.LineSegment3D):
        """Generator line segment intersections approximation."""
        if self.__class__ == PlaneFace3D:
            faces_triangulation = [self]
        else:
            triangulation = self.triangulation()
            faces_triangulation = triangulation.triangular_faces()
        for face in faces_triangulation:
            inters = face.linesegment_intersections(linesegment)
            yield inters

    def is_linesegment_crossing(self, linesegment):
        """
        Verify if a face 3d is being crossed by a line segment 3d.
        """
        intersections = self.linesegment_intersections(linesegment)
        if intersections:
            return True
        return False

    def linesegment_intersections_approximation(self, linesegment: vme.LineSegment3D) -> List[volmdlr.Point3D]:
        """Approximation of intersections face 3D and a line segment 3D."""
        if not self._is_linesegment_intersection_possible(linesegment):
            return []
        linesegment_intersections = []
        for inters in self._get_linesegment_intersections_approximation(linesegment):
            for point in inters:
                if not volmdlr.core.point_in_list(point, linesegment_intersections):
                    linesegment_intersections.append(point)

        return linesegment_intersections

    def face_decomposition(self):
        if not self._face_octree_decomposition:
            self._face_octree_decomposition = octree_face_decomposition(self)
        return self._face_octree_decomposition

    def face_minimum_distance(self, other_face, return_points: bool = False):
        """
        Gets the minimum distance between two faces.

        :param other_face: second face to search for minimum distance.
        :param return_points: return corresponding point or not.
        :return:
        """
        method_name = f'{other_face.__class__.__name__.lower()[:-2]}_minimum_distance'
        if hasattr(self, method_name):
            return getattr(self, method_name)(other_face, return_points)
        face_decomposition1 = self.face_decomposition()
        face_decomposition2 = other_face.face_decomposition()
        list_set_points1 = [{point for face in faces1 for point in face.points}
                            for _, faces1 in face_decomposition1.items()]
        list_set_points1 = [npy.array([(point[0], point[1], point[2]) for point in sets_points1])
                            for sets_points1 in list_set_points1]
        list_set_points2 = [{point for face in faces2 for point in face.points} for _, faces2 in
                            face_decomposition2.items()]
        list_set_points2 = [npy.array([(point[0], point[1], point[2]) for point in sets_points2]) for sets_points2 in
                            list_set_points2]

        minimum_distance = math.inf
        index1, index2 = None, None
        for sets_points1, sets_points2 in product(list_set_points1, list_set_points2):
            distances = npy.linalg.norm(sets_points2[:, npy.newaxis] - sets_points1, axis=2)
            sets_min_dist = npy.min(distances)
            if sets_min_dist < minimum_distance:
                minimum_distance = sets_min_dist
                index1 = next((i for i, x in enumerate(list_set_points1) if npy.array_equal(x, sets_points1)), -1)
                index2 = next((i for i, x in enumerate(list_set_points2) if npy.array_equal(x, sets_points2)), -1)
        faces1 = list(face_decomposition1.values())[index1]
        faces2 = list(face_decomposition2.values())[index2]

        minimum_distance = math.inf
        best_distance_points = None

        for face1, face2 in product(faces1, faces2):
            distance, point1, point2 = face1.planeface_minimum_distance(face2, True)
            if distance < minimum_distance:
                minimum_distance = distance
                best_distance_points = [point1, point2]
        if return_points:
            return minimum_distance, *best_distance_points
        return minimum_distance

    def plane_intersections(self, plane3d: surfaces.Plane3D):
        surfaces_intersections = self.surface3d.plane_intersections(plane3d)
        outer_contour_intersections_with_plane = plane3d.contour_intersections(self.outer_contour3d)
        plane_intersections = []
        for plane_intersection in surfaces_intersections:
            points_on_curve = []
            for point in outer_contour_intersections_with_plane:
                if plane_intersection.point_belongs(point):
                    points_on_curve.append(point)
            points_on_primitive = plane_intersection.sort_points_along_curve(points_on_curve)
            if not isinstance(plane_intersection, volmdlr_curves.Line3D):
                points_on_primitive = points_on_primitive + [points_on_primitive[0]]
            for point1, point2 in zip(points_on_primitive[:-1], points_on_primitive[1:]):
                edge = plane_intersection.trim(point1, point2)
                if self.edge3d_inside(edge):
                    plane_intersections.append(volmdlr.wires.Wire3D([edge]))
        return plane_intersections

    def split_by_plane(self, plane3d: surfaces.Plane3D):
        """Split face with a plane."""
        intersections_with_plane = self.plane_intersections(plane3d)
        intersections_with_plane2d = [self.surface3d.contour3d_to_2d(intersection_wire)
                                      for intersection_wire in intersections_with_plane]
        while True:
            for i, intersection2d in enumerate(intersections_with_plane2d):
                if not self.surface2d.outer_contour.is_inside(intersection2d):
                    for translation in [volmdlr.Vector2D(-2 * math.pi, 0), volmdlr.Vector2D(2 * math.pi, 0)]:
                        translated_contour = intersection2d.translation(translation)
                        if not self.surface2d.outer_contour.is_inside(translated_contour):
                            continue
                        intersections_with_plane2d.pop(i)
                        intersections_with_plane2d.append(translated_contour)
                    break
            else:
                break
        return self.divide_face(intersections_with_plane2d)

    def _get_face_decomposition_set_closest_to_point(self, point):
        """
        Searches for the faces decomposition's set closest to given point.

        :param point: other point.
        :return: list of triangular faces, corresponding to area of the face closest to point.
        """
        face_decomposition1 = self.face_decomposition()
        list_set_points1 = [{point for face in faces1 for point in face.points}
                            for _, faces1 in face_decomposition1.items()]
        list_set_points1 = [npy.array([(point[0], point[1], point[2]) for point in sets_points1])
                            for sets_points1 in list_set_points1]
        list_set_points2 = [npy.array([(point[0], point[0], point[0])])]

        minimum_distance = math.inf
        index1 = None
        for sets_points1, sets_points2 in product(list_set_points1, list_set_points2):
            distances = npy.linalg.norm(sets_points2[:, npy.newaxis] - sets_points1, axis=2)
            sets_min_dist = npy.min(distances)
            if sets_min_dist < minimum_distance:
                minimum_distance = sets_min_dist
                index1 = next((i for i, x in enumerate(list_set_points1) if npy.array_equal(x, sets_points1)), -1)
        return list(face_decomposition1.values())[index1]

    def point_distance(self, point, return_other_point: bool = False):
        """
        Calculates the distance from a face 3d and a point.

        :param point: point to verify.
        :param return_other_point: bool to decide if corresponding point on face should be returned.
        :return: distance to face3D.
        """

        faces1 = self._get_face_decomposition_set_closest_to_point(point)
        minimum_distance = math.inf
        best_distance_point = None
        for face1 in faces1:
            distance, point1 = face1.point_distance(point, True)
            if distance < minimum_distance:
                minimum_distance = distance
                best_distance_point = point1
        if return_other_point:
            return minimum_distance, best_distance_point
        return minimum_distance


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

    def __init__(self, surface3d: surfaces.Plane3D, surface2d: surfaces.Surface2D,
                 name: str = ''):
        self._bbox = None
        Face3D.__init__(self,
                        surface3d=surface3d,
                        surface2d=surface2d,
                        name=name)

    def copy(self, deep=True, memo=None):
        """Returns a copy of the PlaneFace3D."""
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

    def point_distance(self, point, return_other_point=False):
        """
        Calculates the distance from a plane face and a point.

        :param point: point to verify.
        :param return_other_point: bool to decide if corresponding point on face should be returned.
        :return: distance to planeface3D.
        """

        projected_pt = self.surface3d.point_projection(point)
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

    def distance_to_point(self, point, return_other_point=False):
        warnings.warn('distance_to_point is deprecated, please use point_distance', category=DeprecationWarning)
        return self.point_distance(point, return_other_point)

    def minimum_distance_points_plane(self, other_plane_face, return_points=False):
        """
        Given two plane faces, calculates the points which corresponds to the minimal distance between these two faces.

        :param other_plane_face: Second plane face.
        :param return_points: Boolean to return corresponding points or not.
        :return: minimal distance.
        """
        for edge in other_plane_face.outer_contour3d.primitives:
            edge_intersections = self.edge_intersections(edge)
            if edge_intersections:
                if return_points:
                    return 0.0, edge_intersections[0], edge_intersections[0]
                return 0.0
        min_distance = math.inf
        for edge1 in self.outer_contour3d.primitives:
            for edge2 in other_plane_face.outer_contour3d.primitives:
                if hasattr(edge1, 'minimum_distance'):
                    dist = edge1.minimum_distance(edge2, return_points=return_points)
                elif hasattr(edge2, 'minimum_distance'):
                    dist = edge2.minimum_distance(edge1, return_points=return_points)
                else:
                    raise AttributeError(f'Neither {edge1} nor {edge2} has a minimum_distance method.')
                if return_points:
                    if dist[0] < min_distance:
                        min_distance = dist[0]
                        point1, point2 = dist[1], dist[2]
                elif dist < min_distance:
                    min_distance = dist
        if return_points:
            return min_distance, point1, point2
        return min_distance

    def linesegment_inside(self, linesegment: vme.LineSegment3D):
        """
        Verifies if a line segment 3D is completely inside the plane face.

        :param linesegment: the line segment to verify.
        :return: True if circle is inside False otherwise.
        """
        direction_vector = linesegment.unit_direction_vector()
        if not math.isclose(abs(direction_vector.dot(self.surface3d.frame.w)), 0.0, abs_tol=1e-6):
            return False
        for point in [linesegment.start, linesegment.middle_point(), linesegment.end]:
            if not self.point_belongs(point):
                return False
        return True

    def circle_inside(self, circle: volmdlr_curves.Circle3D):
        """
        Verifies if a circle 3D is completely inside the plane face.

        :param circle: the circle to verify.
        :return: True if circle is inside False otherwise.
        """
        if not math.isclose(abs(circle.frame.w.dot(self.surface3d.frame.w)), 1.0, abs_tol=1e-6):
            return False
        points = circle.discretization_points(number_points=4)
        for point in points:
            if not self.point_belongs(point):
                return False
        return True

    def planeface_intersections(self, planeface):
        """
        Calculates the intersections between two plane faces.

        :param planeface: the other Plane Face 3D to verify intersections with Plane Face 3D.
        :return: list of intersecting wires.
        """
        face2_plane_interections = planeface.surface3d.plane_intersections(self.surface3d)
        if not face2_plane_interections:
            return []
        points_intersections = []
        for contour in [self.outer_contour3d, planeface.outer_contour3d] + self.inner_contours3d + \
                planeface.inner_contours3d:
            for intersection in contour.line_intersections(face2_plane_interections[0]):
                if intersection and not volmdlr.core.point_in_list(intersection, points_intersections):
                    points_intersections.append(intersection)
        points_intersections = face2_plane_interections[0].sort_points_along_curve(points_intersections)
        planeface_intersections = []
        for point1, point2 in zip(points_intersections[:-1], points_intersections[1:]):
            linesegment3d = vme.LineSegment3D(point1, point2)
            over_self_outer_contour = self.outer_contour3d.primitive_over_contour(linesegment3d)
            over_planeface_outer_contour = planeface.outer_contour3d.primitive_over_contour(linesegment3d)
            over_self_inner_contour = any(inner_contour.primitive_over_contour(linesegment3d)
                                          for inner_contour in self.inner_contours3d)
            over_planeface_inner_contour = any(inner_contour.primitive_over_contour(linesegment3d)
                                               for inner_contour in planeface.inner_contours3d)
            if over_self_inner_contour and over_planeface_outer_contour:
                continue
            if over_planeface_inner_contour and over_self_outer_contour:
                continue
            if over_self_outer_contour and over_planeface_outer_contour:
                continue
            if self.edge3d_inside(linesegment3d) or over_self_outer_contour:
                if planeface.edge3d_inside(linesegment3d):
                    planeface_intersections.append(volmdlr.wires.Wire3D([linesegment3d]))
                elif over_planeface_outer_contour:
                    planeface_intersections.append(volmdlr.wires.Wire3D([linesegment3d]))
        return planeface_intersections

    def triangle_intersections(self, triangleface):
        """
        Gets the intersections between a Plane Face3D and a Triangle3D.

        :param triangleface: the other triangle face.
        :return:
        """
        return self.planeface_intersections(triangleface)

    def cylindricalface_intersections(self, cylindricalface: 'CylindricalFace3D'):
        """
        Calculates the intersections between a plane face 3D and Cylindrical Face3D.

        :param cylindricalface: the Cylindrical Face 3D to verify intersections with Plane Face 3D.
        :return: list of intersecting wires.
        """
        cylindricalsurfaceface_intersections = cylindricalface.surface3d.plane_intersections(self.surface3d)
        if not isinstance(cylindricalsurfaceface_intersections[0], volmdlr_curves.Line3D):
            if all(self.edge3d_inside(intersection) and cylindricalface.edge3d_inside(intersection)
                   for intersection in cylindricalsurfaceface_intersections):
                if isinstance(cylindricalsurfaceface_intersections[0], volmdlr_curves.Circle3D):
                    contour3d = volmdlr.wires.Contour3D([volmdlr.edges.FullArc3D.from_curve(
                        cylindricalsurfaceface_intersections[0])])
                else:
                    contour3d = volmdlr.wires.Contour3D([volmdlr.edges.FullArcEllipse3D.from_curve(
                        cylindricalsurfaceface_intersections[0])])
                return [contour3d]
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
            points_on_primitive = primitive.sort_points_along_curve(points_on_primitive)
            if not isinstance(primitive, volmdlr_curves.Line3D):
                points_on_primitive = points_on_primitive + [points_on_primitive[0]]
            for point1, point2 in zip(points_on_primitive[:-1], points_on_primitive[1:]):
                edge = primitive.trim(point1, point2)
                if self.edge3d_inside(edge) and cylindricalface.edge3d_inside(edge):
                    face_intersections.append(volmdlr.wires.Wire3D([edge]))
        return face_intersections

    def conicalface_intersections(self, conical_face: 'ConicalFace3D'):
        """
        Calculates the intersections between a plane face 3D and Conical Face3D.

        :param conical_face: the Conical Face 3D to verify intersections with Plane Face 3D.
        :return: list of intersecting wires.
        """
        surface_intersections = self.surface3d.surface_intersections(conical_face.surface3d)
        if isinstance(surface_intersections[0], volmdlr_curves.Circle3D):
            if self.edge3d_inside(surface_intersections[0]) and conical_face.edge3d_inside(surface_intersections[0]):
                contour3d = volmdlr.wires.Contour3D([volmdlr.edges.FullArc3D.from_curve(
                    surface_intersections[0])])
                return [contour3d]
        if isinstance(surface_intersections[0], volmdlr_curves.Ellipse3D):
            if self.edge3d_inside(surface_intersections[0]) and conical_face.edge3d_inside(surface_intersections[0]):
                contour3d = volmdlr.wires.Contour3D([volmdlr.edges.FullArcEllipse3D.from_curve(
                    surface_intersections[0])])
                return [contour3d]
        intersections_points = self.face_intersections_outer_contour(conical_face)
        for point in conical_face.face_intersections_outer_contour(self):
            if not volmdlr.core.point_in_list(point, intersections_points):
                intersections_points.append(point)
        face_intersections = []
        for primitive in surface_intersections:
            points_on_primitive = []
            for point in intersections_points:
                if primitive.point_belongs(point):
                    points_on_primitive.append(point)
            if not points_on_primitive:
                continue
            points_on_primitive = primitive.sort_points_along_curve(points_on_primitive)
            if isinstance(primitive, volmdlr_curves.ClosedCurve):
            # if isinstance(primitive, volmdlr_curves.Ellipse3D) or isinstance(primitive, volmdlr_curves.Circle3D):
                points_on_primitive = points_on_primitive + [points_on_primitive[0]]
            for point1, point2 in zip(points_on_primitive[:-1], points_on_primitive[1:]):
                edge = primitive.trim(point1, point2)
                if self.edge3d_inside(edge) and conical_face.edge3d_inside(edge):
                    face_intersections.append(volmdlr.wires.Wire3D([edge]))
        return face_intersections

    def toroidalface_intersections(self, toroidal_face):
        """
        Calculates the intersections between a plane face 3D and Conical Face3D.

        :param toroidal_face: the Toroidal Face 3D to verify intersections with Plane Face 3D.
        :return: list of intersecting wires.
        """
        surface_intersections = self.surface3d.surface_intersections(toroidal_face.surface3d)
        intersections_points = self.face_intersections_outer_contour(toroidal_face)
        for point in toroidal_face.face_intersections_outer_contour(self):
            if not volmdlr.core.point_in_list(point, intersections_points):
                intersections_points.append(point)
        face_intersections = []
        for primitive in surface_intersections:
            points_on_primitive = []
            for point in intersections_points:
                if primitive.point_belongs(point, 1e-5):
                    points_on_primitive.append(point)
            if not points_on_primitive:
                continue
            points_on_primitive = primitive.sort_points_along_curve(points_on_primitive)
            if primitive.periodic:
                # if isinstance(primitive, volmdlr_curves.Ellipse3D) or isinstance(primitive, volmdlr_curves.Circle3D):
                points_on_primitive = points_on_primitive + [points_on_primitive[0]]
            for point1, point2 in zip(points_on_primitive[:-1], points_on_primitive[1:]):
                edge = primitive.trim(point1, point2)
                if self.edge3d_inside(edge) and toroidal_face.edge3d_inside(edge, 1e-3):
                    face_intersections.append(volmdlr.wires.Wire3D([edge]))
        return face_intersections

    def planeface_minimum_distance(self, planeface: 'PlaneFace3D', return_points: bool = False):
        dist, point1, point2 = self.minimum_distance_points_plane(planeface, return_points=True)
        if not return_points:
            return dist
        return dist, point1, point2

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
        """Merges faces from a list of faces in the same plane, if any are adjacent to one another."""
        list_coincident_faces = sorted(list_coincident_faces, key=lambda face_: face_.area())
        current_face = list_coincident_faces[0]
        list_coincident_faces.remove(current_face)
        list_merged_faces = []
        while True:
            for face in list_coincident_faces:
                if current_face.outer_contour3d.is_sharing_primitives_with(face.outer_contour3d):
                    merged_contours = current_face.outer_contour3d.merge_with(face.outer_contour3d)
                    merged_contours2d = [contour.to_2d(face.surface3d.frame.origin, face.surface3d.frame.u,
                                                       face.surface3d.frame.v) for contour in merged_contours]
                    merged_contours2d = sorted(merged_contours2d, key=lambda contour: contour.area(), reverse=True)
                    if not merged_contours2d and current_face.outer_contour3d.is_superposing(face.outer_contour3d):
                        merged_contours2d = [current_face.surface2d.outer_contour]
                    new_outer_contour = merged_contours2d[0]
                    inner_contours = [contour.to_2d(face.surface3d.frame.origin, face.surface3d.frame.u,
                                                    face.surface3d.frame.v)
                                      for contour in current_face.inner_contours3d]
                    inner_contours += merged_contours2d[1:] + face.surface2d.inner_contours
                    new_face = PlaneFace3D(face.surface3d, surfaces.Surface2D(
                        new_outer_contour, inner_contours))
                    current_face = new_face
                    list_coincident_faces.remove(face)
                    break
                if current_face.face_inside(face):
                    list_coincident_faces.remove(face)
                    break
                new_inner_contours = []
                inner_contour_merged = False
                for inner_contour3d in face.inner_contours3d:
                    if current_face.outer_contour3d.is_sharing_primitives_with(inner_contour3d):
                        merged_inner_contours = current_face.outer_contour3d.merge_with(inner_contour3d)
                        if len(merged_inner_contours) >= 2:
                            raise NotImplementedError
                        new_inner_contours.extend(merged_inner_contours)
                        inner_contour_merged = True
                if inner_contour_merged:
                    list_coincident_faces.remove(face)
                    inner_contours2d = [inner_contour.to_2d(face.surface3d.frame.origin,
                                                            face.surface3d.frame.u,
                                                            face.surface3d.frame.v) for inner_contour in
                                        new_inner_contours]
                    current_face = PlaneFace3D(face.surface3d, surfaces.Surface2D(
                        face.surface2d.outer_contour, inner_contours2d))
                    break
            else:
                list_merged_faces.append(current_face)
                if not list_coincident_faces:
                    break
                current_face = list_coincident_faces[0]
                list_coincident_faces.remove(current_face)
        return list_merged_faces

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
            return self.divide_face([face.surface2d.outer_contour])

        outer_contour_1 = self.surface2d.outer_contour
        outer_contour_2 = self.surface3d.contour3d_to_2d(face.outer_contour3d)

        if (face.face_inside(self)
                and not outer_contour_1.intersection_points(outer_contour_2)):
            return self.divide_face(face.surface2d.inner_contours)

        inner_contours = self.surface2d.inner_contours
        inner_contours.extend([self.surface3d.contour3d_to_2d(
            contour) for contour in face.inner_contours3d])

        contours = outer_contour_1.cut_by_wire(outer_contour_2)

        list_surfaces = []
        for contour in contours:
            inners = []
            for inner_c in inner_contours:
                if contour.is_inside(inner_c):
                    inners.append(inner_c)
            list_surfaces.append(surfaces.Surface2D(contour, inners))

        return [self.__class__(self.surface3d, surface2d) for surface2d in list_surfaces]

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

                        divided_faces_d_face = d_face.divide_face([inner])
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

        for face2 in faces:
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
                    s2d = surfaces.Surface2D(outer_contour=plane3d.contour3d_to_2d(face2_2.outer_contour3d),
                                             inner_contours=[
                                                 plane3d.contour3d_to_2d(contour) for contour in
                                                 face2_2.inner_contours3d])
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

    @classmethod
    def from_surface_rectangular_cut(cls, plane3d, x1: float, x2: float,
                                     y1: float, y2: float, name: str = ''):
        """
        Cut a rectangular piece of the Plane3D object and return a PlaneFace3D object.

        """
        point1 = volmdlr.Point2D(x1, y1)
        point2 = volmdlr.Point2D(x2, y1)
        point3 = volmdlr.Point2D(x2, y2)
        point4 = volmdlr.Point2D(x1, y2)
        outer_contour = volmdlr.wires.Contour2D.from_points([point1, point2, point3, point4])
        surface = surfaces.Surface2D(outer_contour, [])
        return cls(plane3d, surface, name)


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

        PlaneFace3D.__init__(self, self.surface3d, self.surface2d)

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
        """General method to get the bounding box."""
        return volmdlr.core.BoundingBox.from_points([self.point1,
                                                     self.point2,
                                                     self.point3])

    @property
    def surface3d(self):
        if self._surface3d is None:
            self._surface3d = surfaces.Plane3D.from_3_points(self.point1, self.point2, self.point3)
        return self._surface3d

    @surface3d.setter
    def surface3d(self, new_surface3d):
        self._surface3d = new_surface3d

    @property
    def surface2d(self):
        """Boundary representation of the face."""
        if self._surface2d is None:
            plane3d = self.surface3d
            contour3d = volmdlr.wires.Contour3D([vme.LineSegment3D(self.point1, self.point2),
                                                 vme.LineSegment3D(self.point2, self.point3),
                                                 vme.LineSegment3D(self.point3, self.point1)])
            contour2d = contour3d.to_2d(plane3d.frame.origin,
                                        plane3d.frame.u, plane3d.frame.v)

            self._surface2d = surfaces.Surface2D(outer_contour=contour2d, inner_contours=[])

        return self._surface2d

    @surface2d.setter
    def surface2d(self, new_surface2d):
        self._surface2d = new_surface2d

    def to_dict(self, *args, **kwargs):
        """
        Creates a Dictionary with the object's instance attributes.

        """
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
        """
        Gets Triangle height.

        """
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

    def copy(self, deep=True, memo=None):
        """Returns a copy of the Triangle3D."""
        return Triangle3D(self.point1.copy(), self.point2.copy(), self.point3.copy(),
                          self.name)

    def triangulation(self):
        """Computes the triangulation of the Triangle3D, basically returns itself."""
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

    @staticmethod
    def get_subdescription_points(new_points, resolution, max_length):
        """Gets sub-description points."""
        vector = new_points[0] - new_points[1]
        vector = vector.unit_vector()
        points_0_1 = []

        for k in range(int(max_length / resolution) + 2):
            if k == 0:
                points_0_1.append(new_points[1])
            points_0_1.append(new_points[1] + vector * min(k * resolution, max_length))

        vector, length_2_1 = new_points[2] - new_points[1], new_points[2].point_distance(new_points[1])
        vector = vector.unit_vector()
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
                vector_2_0 = vector_2_0.unit_vector()
                step_in = length_2_0 / (nb_int - 1)
                for i in range(nb_int):
                    if min(i * step_in, length_2_0) != 0:
                        points_in.append(p0_1 + vector_2_0 * min(i * step_in, length_2_0))

        return npy.unique(points_0_1 + points_in).tolist()

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
        return self.get_subdescription_points(new_points, resolution, max_length)

    def subdescription_to_triangles(self, resolution=0.01):
        """
        Returns a list of Triangle3D with resolution as max length of sub triangles side.

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
        normal = normal.unit_vector()
        return normal

    def triangle_minimum_distance(self, triangle_face, return_points=False):
        """
        Gets the minimum distance between two triangle.
        """
        return self.planeface_minimum_distance(triangle_face, return_points)


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
                 surface3d: surfaces.CylindricalSurface3D,
                 surface2d: surfaces.Surface2D,
                 name: str = ''):

        self.radius = surface3d.radius
        self.center = surface3d.frame.origin
        self.normal = surface3d.frame.w
        Face3D.__init__(self, surface3d=surface3d,
                        surface2d=surface2d,
                        name=name)
        self._bbox = None

    def copy(self, deep=True, memo=None):
        """Returns a copy of the CylindricalFace3D."""
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

    def triangulation_lines(self, angle_resolution=5):
        """
        Specifies the number of subdivision when using triangulation by lines. (Old triangulation).
        """
        theta_min, theta_max, zmin, zmax = self.surface2d.bounding_rectangle().bounds()
        delta_theta = theta_max - theta_min
        nlines = math.ceil(delta_theta * angle_resolution)
        lines = []
        for i in range(nlines):
            theta = theta_min + (i + 1) / (nlines + 1) * delta_theta
            lines.append(volmdlr_curves.Line2D(volmdlr.Point2D(theta, zmin),
                                               volmdlr.Point2D(theta, zmax)))
        return lines, []

    def point_belongs(self, point3d: volmdlr.Point3D, tol: float = 1e-6):
        """
        Tells you if a point is on the 3D Cylindrical face and inside its contour.
        """
        point2d = self.surface3d.point3d_to_2d(point3d)
        point2d_plus_2pi = point2d.translation(volmdlr.Point2D(volmdlr.TWO_PI, 0))
        point2d_minus_2pi = point2d.translation(volmdlr.Point2D(-volmdlr.TWO_PI, 0))
        # check_point3d = self.surface3d.point2d_to_3d(point2d)
        # if check_point3d.point_distance(point3d) > tol:
        if not self.surface3d.point_on_surface(point3d, tol):
            return False
        return any(self.surface2d.point_belongs(pt2d) for pt2d in [point2d, point2d_plus_2pi, point2d_minus_2pi])

    def parametrized_grid_size(self, angle_resolution, z_resolution):
        # angle_resolution = 5
        # z_resolution = 0
        theta_min, theta_max, zmin, zmax = self.surface2d.bounding_rectangle().bounds()
        delta_theta = theta_max - theta_min
        number_points_x = max(angle_resolution, int(delta_theta * angle_resolution))

        delta_z = zmax - zmin
        number_points_y = min(int(delta_z * z_resolution), 20)
        return number_points_x, number_points_y

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        angle_resolution = 5
        z_resolution = 2
        return self.parametrized_grid_size(angle_resolution, z_resolution)

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
        if not math.isclose(abs(arc.circle.frame.w.dot(self.surface3d.frame.w)), 1.0, abs_tol=1e-6):
            return False
        if not math.isclose(self.radius, arc.circle.radius, abs_tol=1e-6):
            return False
        return self.arcellipse_inside(arc)

    def arcellipse_inside(self, arcellipse: vme.ArcEllipse3D):
        """
        Verifies if ArcEllipse3D is inside a CylindricalFace3D.

        :param arcellipse: ArcEllipse3D to be verified.
        :return: True if it is inside, False otherwise.
        """
        for point in [arcellipse.start, arcellipse.middle_point(), arcellipse.end]:
            if not self.point_belongs(point):
                return False
        return True

    def planeface_intersections(self, planeface: PlaneFace3D):
        planeface_intersections = planeface.cylindricalface_intersections(self)
        return planeface_intersections

    def conicalface_intersections(self, conical_face: 'ConicalFace3D'):
        """
        Calculates the intersections between a plane face 3D and Conical Face3D.

        :param conical_face: the Conical Face 3D to verify intersections with Plane Face 3D.
        :return: list of intersecting wires.
        """
        surface_intersections = self.surface3d.surface_intersections(conical_face.surface3d)
        intersections_points = self.face_intersections_outer_contour(conical_face)
        for point in conical_face.face_intersections_outer_contour(self):
            if not volmdlr.core.point_in_list(point, intersections_points):
                intersections_points.append(point)
        face_intersections = []
        for primitive in surface_intersections:
            points_on_primitive = []
            for point in intersections_points:
                if primitive.point_belongs(point, 1e-4):
                    points_on_primitive.append(point)
            if not points_on_primitive:
                continue
            points_on_primitive = primitive.sort_points_along_curve(points_on_primitive)
            points_on_primitive = points_on_primitive + [points_on_primitive[0]]
            for point1, point2 in zip(points_on_primitive[:-1], points_on_primitive[1:]):
                edge = primitive.trim(point1, point2)
                if self.edge3d_inside(edge, 1e-4) and conical_face.edge3d_inside(edge, 1e-3):
                    face_intersections.append(volmdlr.wires.Wire3D([edge]))
        return face_intersections

    @classmethod
    def from_surface_rectangular_cut(cls, cylindrical_surface, theta1: float, theta2: float,
                                     param_z1: float, param_z2: float, name: str = ''):
        """
        Cut a rectangular piece of the CylindricalSurface3D object and return a CylindricalFace3D object.

        """

        if theta1 == theta2:
            theta2 += volmdlr.TWO_PI

        point1 = volmdlr.Point2D(theta1, param_z1)
        point2 = volmdlr.Point2D(theta2, param_z1)
        point3 = volmdlr.Point2D(theta2, param_z2)
        point4 = volmdlr.Point2D(theta1, param_z2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        surface2d = surfaces.Surface2D(outer_contour, [])
        return cls(cylindrical_surface, surface2d, name)

    def neutral_fiber(self):
        """
        Returns the faces' neutral fiber.
        """
        _, _, zmin, zmax = self.surface2d.outer_contour.bounding_rectangle.bounds()

        point1 = self.surface3d.frame.origin + self.surface3d.frame.w * zmin
        point2 = self.surface3d.frame.origin + self.surface3d.frame.w * zmax
        return volmdlr.wires.Wire3D([vme.LineSegment3D(point1, point2)])


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

    def __init__(self, surface3d: surfaces.ToroidalSurface3D,
                 surface2d: surfaces.Surface2D,
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
        """Returns a copy of the ToroidalFace3D."""
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

    def triangulation_lines(self, angle_resolution=5):
        """
        Specifies the number of subdivision when using triangulation by lines. (Old triangulation).
        """
        theta_min, theta_max, phi_min, phi_max = self.surface2d.bounding_rectangle().bounds()

        delta_theta = theta_max - theta_min
        nlines_x = int(delta_theta * angle_resolution)
        lines_x = []
        for i in range(nlines_x):
            theta = theta_min + (i + 1) / (nlines_x + 1) * delta_theta
            lines_x.append(volmdlr_curves.Line2D(volmdlr.Point2D(theta, phi_min),
                                                 volmdlr.Point2D(theta, phi_max)))
        delta_phi = phi_max - phi_min
        nlines_y = int(delta_phi * angle_resolution)
        lines_y = []
        for i in range(nlines_y):
            phi = phi_min + (i + 1) / (nlines_y + 1) * delta_phi
            lines_y.append(volmdlr_curves.Line2D(volmdlr.Point2D(theta_min, phi),
                                                 volmdlr.Point2D(theta_max, phi)))
        return lines_x, lines_y

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        theta_angle_resolution = 5
        phi_angle_resolution = 2.3
        theta_min, theta_max, phi_min, phi_max = self.surface2d.bounding_rectangle().bounds()

        delta_theta = theta_max - theta_min
        number_points_x = max(theta_angle_resolution, int(delta_theta * theta_angle_resolution))

        delta_phi = phi_max - phi_min
        number_points_y = max(math.ceil(phi_angle_resolution), int(delta_phi * phi_angle_resolution))

        return number_points_x, number_points_y

    @classmethod
    def from_surface_rectangular_cut(cls, toroidal_surface3d, theta1: float = 0.0, theta2: float = volmdlr.TWO_PI,
                                     phi1: float = 0.0, phi2: float = volmdlr.TWO_PI, name: str = ""):
        """
        Cut a rectangular piece of the ToroidalSurface3D object and return a ToroidalFace3D object.

        :param toroidal_surface3d: surface 3d,
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
        return cls(toroidal_surface3d, surfaces.Surface2D(outer_contour, []), name)

    def neutral_fiber(self):
        """
        Returns the faces' neutral fiber.
        """
        theta_min, theta_max, _, _ = self.surface2d.outer_contour.bounding_rectangle.bounds()
        circle = volmdlr_curves.Circle3D(self.surface3d.frame, self.surface3d.tore_radius)
        point1, point2 = [circle.center + circle.radius * math.cos(theta) * circle.frame.u +
                          circle.radius * math.sin(theta) * circle.frame.v for theta in
                          [theta_min, theta_max]]
        return volmdlr.wires.Wire3D([circle.trim(point1, point2)])

    def planeface_intersections(self, planeface: PlaneFace3D):
        """
        Gets intersections between a Toroidal Face 3D and a Plane Face 3D.

        :param planeface: other plane face.
        :return: intersections.
        """
        planeface_intersections = planeface.toroidalface_intersections(self)
        return planeface_intersections


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

    def __init__(self, surface3d: surfaces.ConicalSurface3D,
                 surface2d: surfaces.Surface2D,
                 name: str = ''):
        surface2d_br = surface2d.bounding_rectangle()
        if surface2d_br[0] < 0:
            surface2d = surface2d.translation(volmdlr.Vector2D(2 * math.pi, 0))
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

    def triangulation_lines(self, angle_resolution=5):
        """
        Specifies the number of subdivision when using triangulation by lines. (Old triangulation).
        """
        theta_min, theta_max, zmin, zmax = self.surface2d.bounding_rectangle().bounds()
        delta_theta = theta_max - theta_min
        nlines = int(delta_theta * angle_resolution)
        lines_x = []
        for i in range(nlines):
            theta = theta_min + (i + 1) / (nlines + 1) * delta_theta
            lines_x.append(volmdlr_curves.Line2D(volmdlr.Point2D(theta, zmin),
                                                 volmdlr.Point2D(theta, zmax)))

        if zmin < 1e-9:
            delta_z = zmax - zmin
            lines_y = [volmdlr_curves.Line2D(volmdlr.Point2D(theta_min, zmin + 0.1 * delta_z),
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

    @classmethod
    def from_surface_rectangular_cut(cls, conical_surface3d, theta1: float, theta2: float,
                                     z1: float, z2: float, name: str = ''):
        """
        Cut a rectangular piece of the ConicalSurface3D object and return a ConicalFace3D object.

        """
        # theta1 = angle_principal_measure(theta1)
        # theta2 = angle_principal_measure(theta2)
        if theta1 < 0:
            theta1, theta2 = theta1 + 2 * math.pi, theta2 + 2 * math.pi
        if theta1 == theta2:
            theta2 += volmdlr.TWO_PI

        point1 = volmdlr.Point2D(theta1, z1)
        point2 = volmdlr.Point2D(theta2, z1)
        point3 = volmdlr.Point2D(theta2, z2)
        point4 = volmdlr.Point2D(theta1, z2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        return cls(conical_surface3d, surfaces.Surface2D(outer_contour, []), name)

    @classmethod
    def from_base_and_vertex(cls, conical_surface3d, contour: volmdlr.wires.Contour3D,
                             vertex: volmdlr.Point3D, name: str = ''):
        """
        Returns the conical face defined by the contour of the base and the cone vertex.

        :param conical_surface3d: surface 3d.
        :param contour: Cone, contour base.
        :type contour: volmdlr.wires.Contour3D
        :type vertex: volmdlr.Point3D
        :param name: the name to inject in the new face
        :return: Conical face.
        :rtype: ConicalFace3D
        """
        contour2d = conical_surface3d.contour3d_to_2d(contour)
        start_contour2d = contour2d.primitives[0].start
        end_contour2d = contour2d.primitives[-1].end
        linesegment2d_1 = vme.LineSegment2D(end_contour2d, volmdlr.Point2D(end_contour2d.x, 0))
        linesegment2d_2 = vme.LineSegment2D(volmdlr.Point2D(end_contour2d.x, 0), volmdlr.Point2D(start_contour2d.x, 0))
        linesegment2d_3 = vme.LineSegment2D(volmdlr.Point2D(start_contour2d.x, 0), start_contour2d)

        primitives2d = contour2d.primitives + [linesegment2d_1, linesegment2d_2, linesegment2d_3]
        outer_contour2d = volmdlr.wires.Contour2D(primitives2d)

        surface2d = surfaces.Surface2D(outer_contour=outer_contour2d, inner_contours=[])
        return cls(conical_surface3d, surface2d=surface2d, name=name)

    def neutral_fiber(self):
        """
        Returns the faces' neutral fiber.
        """
        _, _, zmin, zmax = self.surface2d.outer_contour.bounding_rectangle.bounds()

        point1 = self.surface3d.frame.origin + self.surface3d.frame.w * zmin
        point2 = self.surface3d.frame.origin + self.surface3d.frame.w * zmax
        return volmdlr.wires.Wire3D([vme.LineSegment3D(point1, point2)])

    def point_belongs(self, point3d: volmdlr.Point3D, tol: float = 1e-6):
        """
        Tells you if a point is on the 3D conical face and inside its contour.
        """
        if not self.bounding_box.point_belongs(point3d):
            return False
        x, y, z = self.surface3d.frame.global_to_local_coordinates(point3d)
        radius = z * math.tan(self.surface3d.semi_angle)
        point2d = volmdlr.Point2D(0, z)
        if radius != 0.0:
            theta = volmdlr.geometry.sin_cos_angle(x / radius, y / radius)
            if abs(theta) < 1e-9:
                theta = 0.0
            point2d = volmdlr.Point2D(theta, z)

        point2d_plus_2pi = point2d.translation(volmdlr.Point2D(volmdlr.TWO_PI, 0))
        check_point3d = self.surface3d.point2d_to_3d(self.surface3d.point3d_to_2d(point3d))
        if check_point3d.point_distance(point3d) > tol:
            return False

        return self.surface2d.point_belongs(point2d) or self.surface2d.point_belongs(point2d_plus_2pi)

    def circle_inside(self, circle: volmdlr_curves.Circle3D):
        """
        Verifies if a circle 3D lies completely on the Conical face.

        :param circle: Circle to be verified.
        :return: True if circle inside face. False otherwise.
        """
        if not math.isclose(abs(circle.frame.w.dot(self.surface3d.frame.w)), 1.0, abs_tol=1e-6):
            return False
        points = circle.discretization_points(number_points=10)
        for point in points:
            if not self.point_belongs(point):
                return False
        return True


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

    def __init__(self, surface3d: surfaces.SphericalSurface3D,
                 surface2d: surfaces.Surface2D,
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
            lines_x.append(volmdlr_curves.Line2D(volmdlr.Point2D(theta, phi_min),
                                                 volmdlr.Point2D(theta, phi_max)))
        delta_phi = phi_max - phi_min
        nlines_y = int(delta_phi * angle_resolution)
        lines_y = []
        for i in range(nlines_y):
            phi = phi_min + (i + 1) / (nlines_y + 1) * delta_phi
            lines_y.append(volmdlr_curves.Line2D(volmdlr.Point2D(theta_min, phi),
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

    @classmethod
    def from_surface_rectangular_cut(cls, spherical_surface, theta1: float = 0.0, theta2: float = volmdlr.TWO_PI,
                                     phi1: float = - 0.5 * math.pi, phi2: float = 0.5 * math.pi, name=''):
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
        return cls(spherical_surface, surfaces.Surface2D(outer_contour, []), name=name)

    @classmethod
    def from_contours3d_and_rectangular_cut(cls, surface3d, contours: List[volmdlr.wires.Contour3D],
                                            point: volmdlr.Point3D, name: str = ''):
        """
        Face defined by contours and a point indicating the portion of the parametric domain that should be considered.

        :param surface3d: surface 3d.
        :param contours: Cone, contour base.
        :type contours: List[volmdlr.wires.Contour3D]
        :type point: volmdlr.Point3D
        :param name: the name to inject in the new face
        :return: Spherical face.
        :rtype: SphericalFace3D
        """
        inner_contours = []
        point1 = volmdlr.Point2D(-math.pi, -0.5 * math.pi)
        point2 = volmdlr.Point2D(math.pi, -0.5 * math.pi)
        point3 = volmdlr.Point2D(math.pi, 0.5 * math.pi)
        point4 = volmdlr.Point2D(-math.pi, 0.5 * math.pi)
        surface_rectangular_cut = volmdlr.wires.Contour2D.from_points([point1, point2, point3, point4])
        contours2d = [surface3d.contour3d_to_2d(contour) for contour in contours]
        point2d = surface3d.point3d_to_2d(point)
        for contour in contours2d:
            if not contour.point_belongs(point2d):
                inner_contours.append(contour)

        surface2d = surfaces.Surface2D(outer_contour=surface_rectangular_cut,
                                       inner_contours=inner_contours)
        return cls(surface3d, surface2d=surface2d, name=name)


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
                 surface3d: surfaces.RuledSurface3D,
                 surface2d: surfaces.Surface2D,
                 name: str = ''):
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
        """
        General method to get the bounding box.

        To be enhanced by restricting wires to cut
        """
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
            lines.append(volmdlr_curves.Line2D(volmdlr.Point2D(x, ymin),
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

    @classmethod
    def from_surface_rectangular_cut(cls, ruled_surface3d, x1: float = 0.0, x2: float = 1.0,
                                     y1: float = 0.0, y2: float = 1.0, name: str = ''):
        """
        Cut a rectangular piece of the RuledSurface3D object and return a RuledFace3D object.

        """
        point1 = volmdlr.Point2D(x1, y1)
        point2 = volmdlr.Point2D(x2, y1)
        point3 = volmdlr.Point2D(x2, y2)
        point4 = volmdlr.Point2D(x1, y2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        surface2d = surfaces.Surface2D(outer_contour, [])
        return cls(ruled_surface3d, surface2d, name)


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
                 surface3d: surfaces.ExtrusionSurface3D,
                 surface2d: surfaces.Surface2D,
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

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        angle_resolution = 11
        xmin, xmax, _, _ = self.surface2d.bounding_rectangle().bounds()
        delta_x = xmax - xmin
        number_points_x = int(delta_x * angle_resolution)

        number_points_y = number_points_x

        return number_points_x, number_points_y

    @classmethod
    def from_surface_rectangular_cut(cls, extrusion_surface3d, x1: float = 0.0, x2: float = 1.0,
                                     y1: float = 0.0, y2: float = 1.0, name: str = ''):
        """
        Cut a rectangular piece of the ExtrusionSurface3D object and return a ExtrusionFace3D object.

        """
        p1 = volmdlr.Point2D(x1, y1)
        p2 = volmdlr.Point2D(x2, y1)
        p3 = volmdlr.Point2D(x2, y2)
        p4 = volmdlr.Point2D(x1, y2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([p1, p2, p3, p4])
        surface2d = surfaces.Surface2D(outer_contour, [])
        return cls(extrusion_surface3d, surface2d, name)


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
                 surface3d: surfaces.RevolutionSurface3D,
                 surface2d: surfaces.Surface2D,
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

    @classmethod
    def from_surface_rectangular_cut(cls, revolution_surface3d, x1: float, x2: float,
                                     y1: float, y2: float, name: str = ''):
        """
        Cut a rectangular piece of the RevolutionSurface3D object and return a RevolutionFace3D object.

        """
        point1 = volmdlr.Point2D(x1, y1)
        point2 = volmdlr.Point2D(x2, y1)
        point3 = volmdlr.Point2D(x2, y2)
        point4 = volmdlr.Point2D(x1, y2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        surface2d = surfaces.Surface2D(outer_contour, [])
        return cls(revolution_surface3d, surface2d, name)


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

    def __init__(self, surface3d: surfaces.BSplineSurface3D,
                 surface2d: surfaces.Surface2D,
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
        """Creates a bounding box from the face mesh."""
        try:
            number_points_x, number_points_y = self.grid_size()
            if number_points_x >= number_points_y:
                number_points_x, number_points_y = 5, 3
            else:
                number_points_x, number_points_y = 3, 5
            outer_polygon = self.surface2d.outer_contour.to_polygon(angle_resolution=15, discretize_line=True)
            points_grid, x, y, grid_point_index = outer_polygon.grid_triangulation_points(number_points_x,
                                                                                          number_points_y,
                                                                                          include_edge_points=False)
            if self.surface2d.inner_contours:
                points_grid = self._get_bbox_inner_contours_points(points_grid, x, y, grid_point_index)
            points3d = [self.surface3d.point2d_to_3d(point) for point in points_grid]
        except ZeroDivisionError:
            points3d = []
        if not points3d:
            return self.outer_contour3d.bounding_box
        return volmdlr.core.BoundingBox.from_bounding_boxes([volmdlr.core.BoundingBox.from_points(points3d),
                                                            self.outer_contour3d.bounding_box])

    def _get_bbox_inner_contours_points(self, points_grid, x, y, grid_point_index):
        """Helper function to get_bounding_box."""
        for inner_contour in self.surface2d.inner_contours:
            inner_polygon = inner_contour.to_polygon(angle_resolution=5, discretize_line=True)
            # removes with a region search the grid points that are in the inner contour
            xmin, xmax, ymin, ymax = inner_polygon.bounding_rectangle.bounds()
            for i in array_range_search(x, xmin, xmax):
                for j in array_range_search(y, ymin, ymax):
                    point = grid_point_index.get((i, j))
                    if not point:
                        continue
                    if inner_polygon.point_belongs(point):
                        points_grid.remove(point)
                        grid_point_index.pop((i, j))
        return points_grid

    def triangulation_lines(self, resolution=25):
        """
        Specifies the number of subdivision when using triangulation by lines. (Old triangulation).
        """
        u_min, u_max, v_min, v_max = self.surface2d.bounding_rectangle().bounds()

        delta_u = u_max - u_min
        nlines_x = int(delta_u * resolution)
        lines_x = []
        for i in range(nlines_x):
            u = u_min + (i + 1) / (nlines_x + 1) * delta_u
            lines_x.append(volmdlr_curves.Line2D(volmdlr.Point2D(u, v_min),
                                                 volmdlr.Point2D(u, v_max)))
        delta_v = v_max - v_min
        nlines_y = int(delta_v * resolution)
        lines_y = []
        for i in range(nlines_y):
            v = v_min + (i + 1) / (nlines_y + 1) * delta_v
            lines_y.append(volmdlr_curves.Line2D(volmdlr.Point2D(v_min, v),
                                                 volmdlr.Point2D(v_max, v)))
        return lines_x, lines_y

    def grid_size(self):
        """
        Specifies an adapted size of the discretization grid used in face triangulation.
        """
        u_min, u_max, v_min, v_max = self.surface2d.bounding_rectangle().bounds()
        delta_u = u_max - u_min
        resolution_u = self.surface3d.nb_u
        resolution_v = self.surface3d.nb_v
        if self.surface2d.inner_contours:
            resolution_u = max(15, resolution_u)
            resolution_v = max(15, resolution_v)
        delta_v = v_max - v_min
        number_points_x = int(delta_u * resolution_u)
        number_points_y = int(delta_v * resolution_v)

        return number_points_x, number_points_y

    def pair_with(self, other_bspline_face3d):
        """
        Finds out how the uv parametric frames are located.

        It does it by comparing to each other and also how grid 3d can be defined respected to these directions.

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
        """Returns the side of the faces that are adjacent."""
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
        """Returns the side of the faces that are adjacent."""
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
        merged_face = self.from_contours3d(merged_surface, contours)

        return merged_face

    @classmethod
    def from_surface_rectangular_cut(cls, bspline_surface3d, u1: float, u2: float,
                                     v1: float, v2: float, name: str = ''):
        """
        Cut a rectangular piece of the BSplineSurface3D object and return a BSplineFace3D object.

        """
        point1 = volmdlr.Point2D(u1, v1)
        point2 = volmdlr.Point2D(u2, v1)
        point3 = volmdlr.Point2D(u2, v2)
        point4 = volmdlr.Point2D(u1, v2)
        outer_contour = volmdlr.wires.ClosedPolygon2D([point1, point2, point3, point4])
        surface = surfaces.Surface2D(outer_contour, [])
        return BSplineFace3D(bspline_surface3d, surface, name)

    def to_planeface3d(self, plane3d: surfaces.Plane3D = None):
        """
        Converts a Bspline face3d to a Plane face3d (using or without a reference Plane3D).

        :param plane3d: A reference Plane3D, defaults to None
        :type plane3d: Plane3D, optional
        :return: A Plane face3d
        :rtype: PlaneFace3D
        """

        if not plane3d:
            plane3d = self.surface3d.to_plane3d()
        surface2d = surfaces.Surface2D(
            outer_contour=plane3d.contour3d_to_2d(self.outer_contour3d),
            inner_contours=[plane3d.contour3d_to_2d(contour) for contour in self.inner_contours3d])

        return PlaneFace3D(surface3d=plane3d, surface2d=surface2d)

    @staticmethod
    def approximate_with_arc(edge):
        """
        Returns an arc that approximates the given edge.

        :param edge: curve to be approximated by an arc.
        :return: An arc if possible, otherwise None.
        """
        if edge.start.is_close(edge.end):
            start = edge.point_at_abscissa(0.25 * edge.length())
            interior = edge.point_at_abscissa(0.5 * edge.length())
            end = edge.point_at_abscissa(0.75 * edge.length())
            vector1 = interior - start
            vector2 = interior - end
            if vector1.is_colinear_to(vector2) or vector1.norm() == 0 or vector2.norm() == 0:
                return None
            return vme.Arc3D.from_3_points(start, interior, end)
        interior = edge.point_at_abscissa(0.5 * edge.length())
        vector1 = interior - edge.start
        vector2 = interior - edge.end
        if vector1.is_colinear_to(vector2) or vector1.norm() == 0 or vector2.norm() == 0:
            return None
        return vme.Arc3D.from_3_points(edge.start, interior, edge.end)

    def get_approximating_arc_parameters(self, curve_list):
        """
        Approximates the given curves with arcs and returns the arcs, radii, and centers.

        :param curve_list: A list of curves to approximate.
        :type curve_list: list
        :returns: A tuple containing the radius and centers of the approximating arcs.
        :rtype: tuple
        """
        radius = []
        centers = []
        for curve in curve_list:
            if curve.simplify.__class__.__name__ in ("Arc3D", "FullArc3D"):
                arc = curve.simplify
            else:
                arc = self.approximate_with_arc(curve)
            if arc:
                radius.append(arc.circle.radius)
                centers.append(arc.circle.center)
        return radius, centers

    def neutral_fiber_points(self):
        """
        Calculates the neutral fiber points of the face.

        :returns: The neutral fiber points if they exist, otherwise None.
        :rtype: Union[list, None]
        """
        surface_curves = self.surface3d.surface_curves
        u_curves = surface_curves['u']
        v_curves = surface_curves['v']
        u_curves = [primitive.simplify
                    for primitive in u_curves if not isinstance(primitive.simplify, vme.LineSegment3D)]
        v_curves = [primitive.simplify
                    for primitive in v_curves if not isinstance(primitive.simplify, vme.LineSegment3D)]
        u_radius, u_centers = self.get_approximating_arc_parameters(u_curves)
        v_radius, v_centers = self.get_approximating_arc_parameters(v_curves)

        if not u_radius and not v_radius:
            return None
        if v_radius and not u_radius:
            return v_centers
        if u_radius and not v_radius:
            return u_centers
        u_mean = npy.mean(u_radius)
        v_mean = npy.mean(v_radius)
        if u_mean > v_mean:
            return v_centers
        if u_mean < v_mean:
            return u_centers
        return None

    def neutral_fiber(self):
        """
        Returns the faces' neutral fiber.
        """
        neutral_fiber_points = self.neutral_fiber_points()
        is_line = False
        neutral_fiber = None
        if not neutral_fiber_points[0].is_close(neutral_fiber_points[-1]):
            neutral_fiber = vme.LineSegment3D(neutral_fiber_points[0], neutral_fiber_points[-1])
            is_line = all(neutral_fiber.point_belongs(point) for point in neutral_fiber_points)
        if not is_line:
            neutral_fiber = vme.BSplineCurve3D.from_points_interpolation(neutral_fiber_points,
                                                                         min(self.surface3d.degree_u,
                                                                             self.surface3d.degree_v))
        umin, umax, vmin, vmax = self.surface2d.outer_contour.bounding_rectangle.bounds()
        min_bound_u, max_bound_u, min_bound_v, max_bound_v = self.surface3d.domain
        if not math.isclose(umin, min_bound_u, rel_tol=0.01) or not math.isclose(vmin, min_bound_v, rel_tol=0.01):
            point3d_min = self.surface3d.point2d_to_3d(volmdlr.Point2D(umin, vmin))
            point1 = neutral_fiber.point_projection(point3d_min)[0]
        else:
            point1 = neutral_fiber.start
        if not math.isclose(umax, max_bound_u, rel_tol=0.01) or not math.isclose(vmax, max_bound_v, rel_tol=0.01):
            point3d_max = self.surface3d.point2d_to_3d(volmdlr.Point2D(umax, vmax))
            point2 = neutral_fiber.point_projection(point3d_max)[0]
        else:
            point2 = neutral_fiber.end
        return volmdlr.wires.Wire3D([neutral_fiber.trim(point1, point2)])

    def linesegment_intersections(self, linesegment: vme.LineSegment3D) -> List[volmdlr.Point3D]:
        """
        Get intersections between a BSpline face 3d and a Line Segment 3D.

        :param linesegment: other linesegment.
        :return: a list of intersections.
        """
        return self.linesegment_intersections_approximation(linesegment)
