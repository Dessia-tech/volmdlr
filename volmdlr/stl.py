#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
STL reader & writer.

https://en.wikipedia.org/wiki/STL_(file_format)
"""

import struct
import warnings
from typing import List

from binaryornot.check import is_binary
from kaitaistruct import KaitaiStream

import dessia_common.core as dc  # isort: skip
from dessia_common.files import BinaryFile, StringFile  # isort: skip

import volmdlr as vm
import volmdlr.core as vmc
import volmdlr.faces as vmf
from volmdlr import shells


class Stl(dc.DessiaObject):
    """
    STL files are used to represent simple 3D models, defined using triangular 3D faces.

    Initially it was introduced as native format for 3D Systems
    Stereo-lithography CAD system, but due to its extreme simplicity, it
    was adopted by a wide range of 3D modeling, CAD, rapid prototyping
    and 3D printing applications as the simplest 3D model exchange
    format.

    STL is extremely bare-bones format: there are no complex headers, no
    texture / color support, no units specifications, no distinct vertex
    arrays. Whole model is specified as a collection of triangular
    faces.

    There are two versions of the format (text and binary), this spec
    describes binary version.
    """

    _standalone_in_db = True

    _dessia_methods = ['from_text_stream', 'from_text_stream', 'to_closed_shell', 'to_open_shell']

    def __init__(self, triangles: List[vmf.Triangle3D], name: str = ''):
        self.triangles = triangles
        dc.DessiaObject.__init__(self, name=name)

        self.normals = None

    @classmethod
    def points_from_file(cls, filename: str, distance_multiplier=0.001):
        """
        Read points from an STL file and return a list of points.

        :param filename: The path to the STL file.
        :type filename: str
        :param distance_multiplier: (optional) The distance multiplier. Defaults to 0.001.
        :type distance_multiplier: float

        :return: A list of Point3D objects.
        :rtype: List[vm.Point3D]
        """
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
    def from_binary_stream(cls, stream: BinaryFile, distance_multiplier: float = 0.001):
        """
        Create an STL object from a binary stream.

        :param stream: The binary stream containing the STL data.
        :type stream: BinaryFile
        :param distance_multiplier: (optional) The distance multiplier. Defaults to 0.001.
        :type distance_multiplier: float
        :return: An instance of the Stl class.
        :rtype: Stl
        """
        stream.seek(0)

        stream = KaitaiStream(stream)
        name_slice = stream.read_bytes(80)
        try:
            name = name_slice.decode('utf-8')
        except UnicodeDecodeError:
            name = name_slice.decode('latin-1')

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
            except NotImplementedError:
                invalid_triangles.append(i)

            stream.read_u2le()
            # print(abr)
        if invalid_triangles:
            # print('invalid_triangles number: ', len(invalid_triangles))
            for i in invalid_triangles[::-1]:
                del triangles[i]
        return cls(triangles, name=name)

    @classmethod
    def from_text_stream(cls, stream: StringFile,
                         distance_multiplier: float = 0.001):
        """
        Create an STL object from a text stream.

        :param stream: The text stream containing the STL data.
        :type stream: StringFile
        :param distance_multiplier: (optional) The distance multiplier. Defaults to 0.001.
        :type distance_multiplier: float
        :return: An instance of the Stl class.
        :rtype: Stl
        """
        stream.seek(0)

        header = stream.readline()
        name = header[6:]
        triangles = []
        points = []
        for line in stream.readlines():
            if 'vertex' in line:
                line = line.replace('vertex', '')
                line = line.lstrip(' ')
                x, y, z = [i for i in line.split(' ') if i]

                points.append(vm.Point3D(distance_multiplier * float(x),
                                         distance_multiplier * float(y),
                                         distance_multiplier * float(z)))
            if 'endfacet' in line:
                try:
                    triangles.append(vmf.Triangle3D(points[0],
                                                    points[1],
                                                    points[2]))
                except (ZeroDivisionError, NotImplementedError):  # NotImplementedError comes from equal points
                    pass
                points = []
        return cls(triangles, name=name)

    @classmethod
    def from_file(cls, filename: str = None,
                  distance_multiplier: float = 0.001):
        warnings.warn("Use load_from_file instead of from_file",
                      DeprecationWarning)
        return cls.load_from_file(filename, distance_multiplier)

    @classmethod
    def load_from_file(cls, filepath: str, distance_multiplier: float = 0.001):
        """
        Load an STL object from a file.

        :param filepath: The path to the STL file.
        :type filepath: str
        :param distance_multiplier: (optional) The distance multiplier. Defaults to 0.001.
        :type distance_multiplier: float
        :return: An instance of the Stl class.
        :rtype: Stl
        """
        if is_binary(filepath):
            with open(filepath, 'rb') as file:
                return cls.from_binary_stream(
                    file, distance_multiplier=distance_multiplier)

        with open(filepath, 'r', encoding='utf-8', errors='ignore') as file:
            return cls.from_text_stream(
                file, distance_multiplier=distance_multiplier)

    def save_to_binary_file(self, filepath, distance_multiplier=1000):
        """
        Save the STL object into a binary file.

        :param filepath: The path to the STL file.
        :type filepath: str
        :param distance_multiplier: (optional) The distance multiplier. Defaults to 1000.
        :type distance_multiplier: float
        :return: An instance of the Stl class.
        :rtype: Stl
        """
        if not filepath.endswith('.stl'):
            filepath += '.stl'
            print('Adding .stl extension: ', filepath)

        with open(filepath, 'wb') as file:
            self.save_to_stream(file, distance_multiplier=distance_multiplier)

    def save_to_stream(self, stream, distance_multiplier=1000):
        """
        Save the STL object into a binary file.

        :param stream: The binary stream containing the STL data.
        :type filepath: BinaryFile
        :param distance_multiplier: (optional) The distance multiplier. Defaults to 1000.
        :type distance_multiplier: float
        :return: An instance of the Stl class.
        :rtype: Stl
        """
        stream.seek(0)

        binary_header = "80sI"
        binary_facet = "12fH"

        # counter = 0
        stream.write(struct.pack(binary_header, self.name.encode('utf8'),
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
            stream.write(struct.pack(binary_facet, *data))

    def to_closed_shell(self):
        """
        Convert the STL object to a closed triangle shell.

        :return: A closed triangle shell representation of the STL object.
        :rtype: shells.ClosedTriangleShell3D
        """
        return shells.ClosedTriangleShell3D(self.triangles, name=self.name)

    def to_open_shell(self):
        """
        Convert the STL object to an open triangle shell.

        :return: An open triangle shell representation of the STL object.
        :rtype: shells.OpenTriangleShell3D
        """
        return shells.OpenTriangleShell3D(self.triangles, name=self.name)

    def to_volume_model(self):
        """
        Convert the STL object to a volume model.

        :return: A volume model representation of the STL object.
        :rtype: vmc.VolumeModel
        """
        closed_shell = self.to_closed_shell()
        return vmc.VolumeModel([closed_shell], name=self.name)

    def extract_points(self):
        """
        Extract the unique points from the STL object.

        :return: A list of unique Point3D objects.
        :rtype: List[vm.Point3D]
        """
        points1 = [triangle.point1 for triangle in self.triangles]
        points2 = [triangle.point2 for triangle in self.triangles]
        points3 = [triangle.point3 for triangle in self.triangles]

        valid_points = vm.Vector3D.remove_duplicate(points1 + points2 + points3)
        return valid_points

    # TODO: decide which algorithm to be used (no _BIS)
    def extract_points_bis(self, min_distance: float = 0.001):
        points = []
        for triangle in self.triangles:
            distance12 = triangle.point1.point_distance(triangle.point2)
            distance13 = triangle.point1.point_distance(triangle.point3)
            distance23 = triangle.point2.point_distance(triangle.point3)
            if distance12 > min_distance:
                n_div = int(distance12 / min_distance)
                for n in range(n_div):
                    new_point = triangle.point1 + (triangle.point2 - triangle.point1) * n / n_div
                    points.append(new_point)
            if distance13 > min_distance:
                n_div = int(distance13 / min_distance)
                for n in range(n_div):
                    new_point = triangle.point1 + (triangle.point3 - triangle.point1) * (n + 1) / n_div
                    points.append(new_point)
            if distance23 > min_distance:
                n_div = int(distance23 / min_distance)
                for n in range(n_div):
                    new_point = triangle.point2 + (triangle.point3 - triangle.point2) * n / n_div
                    points.append(new_point)

        valid_points = vm.Vector3D.remove_duplicate(points)
        return valid_points

    @classmethod
    def from_display_mesh(cls, mesh: vm.display.DisplayMesh3D):
        """
        Create an STL object from a display mesh.

        :param mesh: The display mesh to convert to an STL object.
        :type mesh: vm.display.DisplayMesh3D
        :return: An instance of the Stl class.
        :rtype: Stl
        """
        triangles = []
        for i1, i2, i3 in mesh.triangles:
            triangles.append(vmf.Triangle3D(mesh.points[i1],
                                            mesh.points[i2],
                                            mesh.points[i3]))
        return cls(triangles)

    def get_normals(self):
        """
        Gets normals.

        points_normals : dictionary
            returns a diction
        """
        points_normals = {}
        normals = []
        for triangle in self.triangles:
            normal = triangle.normal()
            for point in triangle.points:
                try:
                    points_normals[point].append(normal)
                except KeyError:
                    points_normals[point] = [normal]

        for key, value in points_normals.items():
            point_normal = vm.O3D
            for point in value:
                point_normal += point
            points_normals[key] = point_normal
            try:
                point_normal = point_normal.unit_vector()
            except ZeroDivisionError:
                point_normal = value[0]
                points_normals[key] = point_normal
                point_normal = point_normal.unit_vector()
            normals.append(point_normal)
        self.normals = normals
        return points_normals

    def clean_flat_triangles(self, threshold: float = 1e-12) -> 'Stl':
        """
        Clean the STL object by removing flat triangles with an area below a threshold.

        :return: A new instance of the Stl class with the flat triangles removed.
        :rtype: Stl
        """
        invalid_triangles = []
        for index_t, triangles in enumerate(self.triangles):
            if triangles.area() < threshold:
                invalid_triangles.append(index_t)

        triangles = self.triangles[:]
        for invalid_triangle_index in invalid_triangles[::-1]:
            triangles.pop(invalid_triangle_index)
        return Stl(triangles)
