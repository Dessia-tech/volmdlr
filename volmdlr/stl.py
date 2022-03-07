#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
# from binaryornot.check import is_binary
import io
from typing import BinaryIO
from typing import List
import struct
# import kaitaistruct
from kaitaistruct import KaitaiStruct, KaitaiStream, BytesIO


from binaryornot.check import is_binary

# from kaitaistruct import KaitaiStream

import dessia_common as dc
import volmdlr as vm
import volmdlr.faces as vmf


class Stl(dc.DessiaObject):
    """STL files are used to represent simple 3D models, defined using
    triangular 3D faces.

    Initially it was introduced as native format for 3D Systems
    Stereolithography CAD system, but due to its extreme simplicity, it
    was adopted by a wide range of 3D modelling, CAD, rapid prototyping
    and 3D printing applications as the simplest 3D model exchange
    format.

    STL is extremely bare-bones format: there are no complex headers, no
    texture / color support, no units specifications, no distinct vertex
    arrays. Whole model is specified as a collection of triangular
    faces.

    There are two versions of the format (text and binary), this spec
    describes binary version.
    """

    _dessia_methods = ['from_text_stream', 'from_text_stream', 'to_closed_shell', 'to_open_shell']

    def __init__(self, triangles: List[vmf.Triangle3D], name: str = ''):
        self.triangles = triangles
        self.name = name
        self.normals = None

    @classmethod
    def points_from_file(cls, filename: str, distance_multiplier=0.001):
        if is_binary(filename):
            with open(filename, 'rb') as file:
                stream = KaitaiStream(file)
                _ = stream.read_bytes(80).decode('utf8')
                num_triangles = stream.read_u4le()

                all_points = []
                for i in range(num_triangles):
                    if i % 5000 == 0:
                        print('reading stl',
                              round(i / num_triangles * 100, 2), '%')
                    # First is normal, unused
                    _ = vm.Vector3D(stream.read_f4le(),
                                    stream.read_f4le(),
                                    stream.read_f4le())
                    p1 = vm.Point3D(distance_multiplier * stream.read_f4le(),
                                    distance_multiplier * stream.read_f4le(),
                                    distance_multiplier * stream.read_f4le())
                    p2 = vm.Point3D(distance_multiplier * stream.read_f4le(),
                                    distance_multiplier * stream.read_f4le(),
                                    distance_multiplier * stream.read_f4le())
                    p3 = vm.Point3D(distance_multiplier * stream.read_f4le(),
                                    distance_multiplier * stream.read_f4le(),
                                    distance_multiplier * stream.read_f4le())
                    all_points.extend([p1, p2, p3])

                    stream.read_u2le()
        return all_points

    @classmethod
    def from_binary_stream(cls, stream: io.BytesIO, distance_multiplier: float = 0.001):
        stream = KaitaiStream(stream)
        name = stream.read_bytes(80).decode('utf8')
        # print(name)
        num_triangles = stream.read_u4le()
        # print(num_triangles)

        triangles = [None] * num_triangles
        invalid_triangles = []
        for i in range(num_triangles):
            if i % 5000 == 0:
                print('reading stl',
                      round(i / num_triangles * 100, 2), '%')
            _ = vm.Vector3D(stream.read_f4le(),
                            stream.read_f4le(),
                            stream.read_f4le())
            p1 = vm.Point3D(distance_multiplier * stream.read_f4le(),
                            distance_multiplier * stream.read_f4le(),
                            distance_multiplier * stream.read_f4le())
            p2 = vm.Point3D(distance_multiplier * stream.read_f4le(),
                            distance_multiplier * stream.read_f4le(),
                            distance_multiplier * stream.read_f4le())
            p3 = vm.Point3D(distance_multiplier * stream.read_f4le(),
                            distance_multiplier * stream.read_f4le(),
                            distance_multiplier * stream.read_f4le())
            try:
                triangles[i] = vmf.Triangle3D(p1, p2, p3)
            except ZeroDivisionError:
                invalid_triangles.append(i)

            stream.read_u2le()
            # print(abr)
        if invalid_triangles:
            # print('invalid_triangles number: ', len(invalid_triangles))
            for i in invalid_triangles[::-1]:
                del triangles[i]
        return cls(triangles, name=name)

    @classmethod
    def from_text_stream(cls, stream: io.StringIO, distance_multiplier: float = 0.001):
        header = stream.readline()
        name = header[6:]
        triangles = []
        points = []
        for line in stream.readlines():
            if 'vertex' in line:
                line = line.replace('vertex', '')
                line = line.lstrip(' ')
                x, y, z = line.split(' ')
                points.append(vm.Point3D(distance_multiplier * float(x),
                                         distance_multiplier * float(y),
                                         distance_multiplier * float(z)))
            if 'endfacet' in line:
                try:
                    triangles.append(vmf.Triangle3D(points[0], points[1], points[2]))
                except ZeroDivisionError:
                    pass
                points = []
        return cls(triangles, name=name)

    @classmethod
    def from_file(cls, filename: str = None, distance_multiplier: float = 0.001):
        if is_binary(filename):
            with open(filename, 'rb') as file:
                return cls.from_binary_stream(file, distance_multiplier=distance_multiplier)

        with open(filename, 'r', errors='ignore') as file:
            return cls.from_text_stream(file, distance_multiplier=distance_multiplier)

    # Commented because seemed invalid
    # @classmethod
    # def from_file_points(cls, filename:str, distance_multiplier=0.001):
    #     all_points = []
    #     if is_binary(filename):
    #         with open(filename, 'rb') as file:
    #             stream = KaitaiStream(file)
    #             name = stream.read_bytes(80).decode('utf8')
    #             # print(name)
    #             num_triangles = stream.read_u4le()
    #             # print(num_triangles)
    #             all_points = []
    #             for i in range(num_triangles):
    #                 if i % 5000 == 0:
    #                     print('reading stl', round(i/num_triangles*100, 2),
    #                           '%')
    #                 _ = vm.Vector3D(stream.read_f4le(),
    #                                      stream.read_f4le(),
    #                                      stream.read_f4le())
    #                 # print(n)
    #                 p1 = vm.Point3D(distance_multiplier*stream.read_f4le(),
    #                                 distance_multiplier*stream.read_f4le(),
    #                                 distance_multiplier*stream.read_f4le())
    #                 p2 = vm.Point3D(distance_multiplier*stream.read_f4le(),
    #                                 distance_multiplier*stream.read_f4le(),
    #                                 distance_multiplier*stream.read_f4le())
    #                 p3 = vm.Point3D(distance_multiplier*stream.read_f4le(),
    #                                 distance_multiplier*stream.read_f4le(),
    #                                 distance_multiplier*stream.read_f4le())
    #                 # print(p1, p2, p3)
    #                 all_points.extend([p1, p2, p3])
    #                 stream.read_u2le()
    #     return all_points

    def save_to_binary_file(self, filepath, distance_multiplier=1000):
        if not filepath.endswith('.stl'):
            filepath += '.stl'
            print('Adding .stl extension: ', filepath)

        with open(filepath, 'wb') as file:
            self.to_stream(file, distance_multiplier=distance_multiplier)

    def save_to_stream(self, stream, distance_multiplier=1000):
        stream.seek(0)

        BINARY_HEADER = "80sI"
        BINARY_FACET = "12fH"

        # counter = 0
        stream.write(struct.pack(BINARY_HEADER, self.name.encode('utf8'),
                                 len(self.triangles)))
        # counter += 1
        for triangle in self.triangles:
            data = [
                0., 0., 0.,
                distance_multiplier * triangle.point1.x,
                distance_multiplier * triangle.point1.y,
                distance_multiplier * triangle.point1.z,
                distance_multiplier * triangle.point2.x,
                distance_multiplier * triangle.point2.y,
                distance_multiplier * triangle.point2.z,
                distance_multiplier * triangle.point3.x,
                distance_multiplier * triangle.point3.y,
                distance_multiplier * triangle.point3.z,
                0]
        stream.write(struct.pack(BINARY_FACET, *data))

    def to_closed_shell(self):
        return vmf.ClosedShell3D(self.triangles, name=self.name)

    def to_open_shell(self):
        return vmf.OpenShell3D(self.triangles, name=self.name)

    def extract_points(self):

        points1 = [t.point1 for t in self.triangles]
        points2 = [t.point2 for t in self.triangles]
        points3 = [t.point3 for t in self.triangles]

        valid_points = vm.Vector3D.remove_duplicate(points1 + points2 + points3)
        return valid_points

    # TODO: decide which algorithm to be used (no _BIS)
    def extract_points_BIS(self, min_distance: float = 0.001):
        points = []
        for i, t in enumerate(self.triangles):
            distance12 = t.point1.point_distance(t.point2)
            distance13 = t.point1.point_distance(t.point3)
            distance23 = t.point2.point_distance(t.point3)
            if distance12 > min_distance:
                n_div = int(distance12 / min_distance)
                for n in range(n_div):
                    new_point = t.point1 + (t.point2 - t.point1) * n / n_div
                    points.append(new_point)
            if distance13 > min_distance:
                n_div = int(distance13 / min_distance)
                for n in range(n_div):
                    new_point = t.point1 + (t.point3 - t.point1) * (n + 1) / n_div
                    points.append(new_point)
            if distance23 > min_distance:
                n_div = int(distance23 / min_distance)
                for n in range(n_div):
                    new_point = t.point2 + (t.point3 - t.point2) * n / n_div
                    points.append(new_point)

        valid_points = vm.Vector3D.remove_duplicate(points)
        return valid_points

    @classmethod
    def from_display_mesh(cls, mesh):
        triangles = []
        for i1, i2, i3 in mesh.triangles:
            triangles.append(vmf.Triangle3D(mesh.points[i1],
                                            mesh.points[i2],
                                            mesh.points[i3]))
        return cls(triangles)

    def get_normals(self):
        """
        Returns
        -------
        points_normals : dictionary
            returns a diction
        """
        points_normals = {}
        normals = []
        for triangle in self.triangles:
            normal = triangle.normal()
            for point in triangle.points:
                if point in list(points_normals.keys()):
                    points_normals[point].append(normal)
                else:
                    points_normals[point] = [normal]
        for key, value in points_normals.items():
            point_normal = vm.O3D
            for point in value:
                point_normal += point
            points_normals[key] = point_normal
            try:
                point_normal.normalize()
            except ZeroDivisionError:
                point_normal = value[0]
                points_normals[key] = point_normal
                point_normal.normalize()
            normals.append(point_normal)
        self.normals = normals

        return points_normals
