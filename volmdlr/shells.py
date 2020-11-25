#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
from typing import List, Tuple
import volmdlr.core
import volmdlr.display
import volmdlr.faces

class Shell3D(volmdlr.core.CompositePrimitive3D):
    _standalone_in_db = True
    _non_serializable_attributes = ['bounding_box']
    _non_eq_attributes = ['name', 'color', 'alpha' 'bounding_box']
    _non_hash_attributes = []

    def __init__(self, faces:List[volmdlr.faces.Face3D],
                 color:Tuple[float, float, float]=None,
                 alpha:float=1., name:str=''):
        self.faces = faces
        self.name = name
        self.color = color
        self.alpha = alpha
        self.bounding_box = self._bounding_box()

    def __hash__(self):
        return sum([hash(f) for f in self.faces])

    def __eq__(self, other_):
        if self.__class__ != other_.__class__:
            return False
        equal = True
        for face, other_face in zip(self.faces, other_.faces):
            equal = (equal and face == other_face)
        return equal

    @classmethod
    def from_step(cls, arguments, object_dict):
        faces = []
        for face in arguments[1]:
            faces.append(object_dict[int(face[1:])])
        return cls(faces, name=arguments[0][1:-1])

    def rotation(self, center, axis, angle, copy=True):
        if copy:
            new_faces = [face.rotation(center, axis, angle, copy=True) for face
                         in self.faces]
            return Shell3D(new_faces, name=self.name)
        else:
            for face in self.faces:
                face.rotation(center, axis, angle, copy=False)
            self.bounding_box = self._bounding_box()

    def translation(self, offset, copy=True):
        if copy:
            new_faces = [face.translation(offset, copy=True) for face in
                         self.faces]
            return Shell3D(new_faces, name=self.name)
        else:
            for face in self.faces:
                face.translation(offset, copy=False)
            self.bounding_box = self._bounding_box()

    def frame_mapping(self, frame, side, copy=True):
        """
        side = 'old' or 'new'
        """
        if copy:
            new_faces = [face.frame_mapping(frame, side, copy=True) for face in
                         self.faces]
            return Shell3D(new_faces, name=self.name)
        else:
            for face in self.faces:
                face.frame_mapping(frame, side, copy=False)
            self.bounding_box = self._bounding_box()

    def copy(self):
        new_faces = [face.copy() for face in self.faces]
        return Shell3D(new_faces, name=self.name)

    def union(self, shell2):
        new_faces = [face for face in self.faces + shell2.faces]
        new_name = self.name + ' union ' + shell2.name
        new_color = self.color
        return Shell3D(new_faces, name=new_name, color=new_color)

    def volume(self):
        """
        Does not consider holes
        """
        volume = 0
        for i, face in enumerate(self.faces):
            points_3D, triangles_indexes = face.triangulation()
            for triangle_indexes in triangles_indexes[0]:
                point1 = points_3D[triangle_indexes[0]]
                point2 = points_3D[triangle_indexes[1]]
                point3 = points_3D[triangle_indexes[2]]

                v321 = point3[0] * point2[1] * point1[2]
                v231 = point2[0] * point3[1] * point1[2]
                v312 = point3[0] * point1[1] * point2[2]
                v132 = point1[0] * point3[1] * point2[2]
                v213 = point2[0] * point1[1] * point3[2]
                v123 = point1[0] * point2[1] * point3[2]
                volume_tetraedre = 1 / 6 * (
                            -v321 + v231 + v312 - v132 - v213 + v123)

                volume += volume_tetraedre

        return abs(volume)

    def _bounding_box(self):
        """
        Returns the boundary box
        """
        points = []
        bbox = self.faces[0]._bounding_box()

        for face in self.faces:
            bbox += face._bounding_box()

        return bbox

    def point_belongs(self, point3d:volmdlr.Point3D, nb_rays:int=1):
        """
        Ray Casting algorithm
        Returns True if the point is inside the Shell, False otherwise
        """


        bbox = self.bounding_box
        if not bbox.point_belongs(point3d):
            print('bbox')
            print(bbox.xmin, point3d.x, bbox.xmax)
            print(bbox.ymin, point3d.y, bbox.ymax)
            print(bbox.zmin, point3d.z, bbox.zmax)
            print(point3d)

            return False

        min_ray_length = 2*max((bbox.xmax - bbox.xmin,
                                bbox.ymax - bbox.ymin,
                                bbox.zmax - bbox.zmin))

        rays = []
        for k in range(0, nb_rays):
            rays.append(volmdlr.edges.LineSegment3D(
                point3d,
                point3d+volmdlr.Point3D.random(0., min_ray_length,
                                               0., min_ray_length,
                                               0., min_ray_length)))
        print('rays', rays[0].start, rays[0].end)
        # rays = sorted(rays, key=lambda ray: ray.length())

        rays_intersections = []
        tests = []

        # for ray in rays[:3]:
        for ray in rays[:nb_rays]:
            #
            count = 0
            ray_intersection = []
            is_inside = True
            intersections = self.linesegment_intersections(ray)
                # if intersection_point is not None:
                #     ray_intersection.append(intersection_point)
            count += len(intersections)
            if count % 2 == 0:
                is_inside = False
                print('count', count)
                print('is_inside', is_inside)
                # ax = point3d.plot()
                # ray.plot(ax=ax)
            tests.append(is_inside)
            rays_intersections.append(ray_intersection)

        for test1, test2 in zip(tests[:-1], tests[1:]):
            if test1 != test2:
                raise ValueError

        return tests[0]

    def is_inside_shell(self, shell2, resolution:float):
        """
        Returns True if all the points of self are inside shell2 and no face \
        are intersecting
        This method is not exact
        """
        bbox1 = self.bounding_box
        bbox2 = shell2.bounding_box
        if not bbox1.is_inside_bbox(bbox2):
            return False

        points = []
        for face in self.faces:
            points.extend(face.outer_contour3d.discretization_points(resolution))
        ax = points[0].plot()
        for point in points:
            # print(point)
            point.plot(ax=ax)
            if not shell2.point_belongs(point):
                print('False')
                return False

        # Check if any faces are intersecting
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    #                    print('Two faces are intersecting :', face1, face2)
                    return False

        return True

    def linesegment_intersections(self,
                                 linesegment3d:volmdlr.edges.LineSegment3D):
        intersections = []
        for face in self.faces:
            face_intersections = face.linesegment_intersections(linesegment3d)
            print(face_intersections)
            if face_intersections:
                intersections.extend(face_intersections)
        return intersections

    def shell_intersection(self, shell2:'Shell3D', resolution:float):
        """
        Return None if disjointed
        Return (1, 0) or (0, 1) if one is inside the other
        Return (n1, n2) if intersection

        4 cases :
            (n1, n2) with face intersection             => (n1, n2)
            (0, 0) with face intersection               => (0, 0)
            (0, 0) with no face intersection            => None
            (1, 0) or (0, 1) with no face intersection  => 1
        """
        # Check if boundary boxes don't intersect
        bbox1 = self.bounding_box
        bbox2 = shell2.bounding_box
        if not bbox1.bbox_intersection(bbox2):
            #            print("No intersection of shells' BBox")
            return None

        # Check if any point of the first shell is in the second shell
        points1 = []
        for face in self.faces:
            points1.extend(face.outer_contour3d.discretization_points(resolution))
        points2 = []
        for face in shell2.faces:
            points2.extend(face.outer_contour3d.discretization_points(resolution))

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
        #        print('shell intersection')
        #        print('shell1 intersecte shell2 à', inter1*100, '%')
        #        print('shell2 intersecte shell1 à', inter2*100, '%')

        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    #                    print('Two faces are intersecting :', face1, face2)
                    #                    ax = face1.plot()
                    #                    face2.plot(ax)
                    return inter1, inter2
        if (inter1, inter2) == (0, 0):
            return None
        return 1

    def minimum_distance_points(self, shell2, resolution):
        """
        Returns a Mesure object if the distance is not zero, otherwise returns None
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
                    elif distance < distance_min:
                        distance_min, point1_min, point2_min = distance, point1, point2

        return point1_min, point2_min

    def distance_to_shell(self, other_shell:'Shell3D', resolution:float):
        p1, p2 = self.minimum_distance_points(other_shell, resolution)
        return p1.point_distance(p2)

        # mesure = volmdlr.measures.Measure3D(point1_min, point2_min)
        #
        # if add_to_volumemodel is not None:
        #     add_to_volumemodel.primitives.append(mesure)
        #
        # return mesure

    def distance_to_point(self, point, add_to_volumemodel=None):
        """
        Computes the distance of a point to a Shell3D, whether it is inside or outside the Shell3D
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

        mesure = Measure3D(point, point1_min)

        if add_to_volumemodel is not None:
            add_to_volumemodel.primitives.append(mesure)

        return mesure

    def intersection_internal_aabb_volume(self, shell2, resolution):
        """
        aabb made of the intersection points and the points of self internal to shell2
        """
        intersections_points = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    intersections_points.extend(intersection_points)

        shell1_points_inside_shell2 = []
        for face in self.faces:
            for point in face.outer_contour3d.discretization_points(resolution):
                if shell2.point_belongs(point):
                    shell1_points_inside_shell2.append(point)

        if len(intersections_points + shell1_points_inside_shell2) == 0:
            return 0
        bbox = BoundingBox.from_points(
            intersections_points + shell1_points_inside_shell2)
        return bbox.volume()

    def intersection_external_aabb_volume(self, shell2, resolution):
        """
        aabb made of the intersection points and the points of self external to shell2
        """
        intersections_points = []
        for face1 in self.faces:
            for face2 in shell2.faces:
                intersection_points = face1.face_intersection(face2)
                if intersection_points is not None:
                    intersections_points.extend(intersection_points)

        shell1_points_outside_shell2 = []
        for face in self.faces:
            for point in face.outer_contour3d.discretization_points(resolution):
                if not shell2.point_belongs(point):
                    shell1_points_outside_shell2.append(point)

        if len(intersections_points + shell1_points_outside_shell2) == 0:
            return 0
        bbox = volmdlr.core.BoundingBox.from_points(
            intersections_points + shell1_points_outside_shell2)
        return bbox.volume()

    def triangulation(self):

        mesh = volmdlr.display.DisplayMesh3D([], [])
        for i, face in enumerate(self.faces):

            face_mesh = face.triangulation()
            mesh += face_mesh
        return mesh

    def babylon_script(self, name='primitive_mesh'):
        s = 'var {} = new BABYLON.Mesh("{}", scene);\n'.format(name, name)

        mesh = self.babylon_meshes()[0]

        s += 'var positions = {};\n'.format(mesh['positions'])
        s += 'var indices = {};\n'.format(mesh['indices'])
        s += 'var normals = [];\n'
        s += 'var vertexData = new BABYLON.VertexData();\n'
        s += 'BABYLON.VertexData.ComputeNormals(positions, indices, normals);\n'
        s += 'vertexData.positions = positions;\n'
        s += 'vertexData.indices = indices;\n'
        s += 'vertexData.normals = normals;\n'
        s += 'vertexData.applyToMesh({});\n'.format(name)
        s += '{}.enableEdgesRendering(0.9);\n'.format(name)
        s += '{}.edgesWidth = 0.1;\n'.format(name)
        s += '{}.edgesColor = new BABYLON.Color4(0, 0, 0, 0.6);\n'.format(name)
        s += 'var mat = new BABYLON.StandardMaterial("mat", scene);\n'
        #        s += 'mat.diffuseColor = BABYLON.Color3.Green();\n'
        #        s += 'mat.specularColor = new BABYLON.Color3(0.5, 0.6, 0.87);\n'
        #        s += 'mat.emissiveColor = new BABYLON.Color3(1, 1, 1);\n'
        #        s += 'mat.ambientColor = new BABYLON.Color3(0.23, 0.98, 0.53);\n'
        s += 'mat.backFaceCulling = false;\n'
        s += 'mat.alpha = {};\n'.format(self.alpha)
        s += '{}.material = mat;\n'.format(name)
        if self.color is not None:
            s += 'mat.diffuseColor = new BABYLON.Color3({}, {}, {});\n'.format(
                *self.color)
        return s

    def plot(self, ax=None, equal_aspect=True, color='k', alpha=1):
        if ax is None:
            ax = Axes3D(plt.figure())


        for face in self.faces:
            face.plot(ax=ax, color=color, alpha=alpha)

        return ax
