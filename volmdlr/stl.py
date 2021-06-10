#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""
import kaitaistruct
from kaitaistruct import KaitaiStruct, KaitaiStream, BytesIO

import struct
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
    def __init__(self, triangles, name=''):
        self.triangles = triangles
        self.name = name

    @classmethod
    def from_file(cls, filename:str, distance_multiplier=0.001):
        # TODO test binary
        with open(filename, 'rb') as file:
            stream = KaitaiStream(file)
            name = stream.read_bytes(80).decode('utf8')
            # print(name)
            num_triangles = stream.read_u4le()
            # print(num_triangles)
            
            triangles = [None] * (num_triangles)
            invalid_triangles = []
            for i in range(num_triangles):
                if i % 5000 == 0:
                    print('reading stl', round(i/num_triangles*100, 2), '%')
                normal = vm.Vector3D(stream.read_f4le(), stream.read_f4le(), stream.read_f4le())
                # print(n)
                p1 = vm.Point3D(distance_multiplier*stream.read_f4le(),
                                distance_multiplier*stream.read_f4le(),
                                distance_multiplier*stream.read_f4le())
                p2 = vm.Point3D(distance_multiplier*stream.read_f4le(),
                                distance_multiplier*stream.read_f4le(), 
                                distance_multiplier*stream.read_f4le())
                p3 = vm.Point3D(distance_multiplier*stream.read_f4le(),
                                distance_multiplier*stream.read_f4le(),
                                distance_multiplier*stream.read_f4le())
                # print(p1, p2, p3)
                try : 
                    triangles[i] = vmf.Triangle3D(p1, p2, p3)
                except ZeroDivisionError :
                    invalid_triangles.append(i)
                    
    
                stream.read_u2le()
                # print(abr)
        if invalid_triangles :
            print('invalid_triangles number: ', len(invalid_triangles))
            for i in invalid_triangles[::-1] :
                del triangles[i]
        return cls(triangles, name=name)

    def save_to_binary_file(self, filepath):
        BINARY_HEADER ="80sI"
        BINARY_FACET = "12fH"
        if not filepath.endswith('.stl'):
            filepath += '.stl'
            print('Adding .stl extension: ', filepath)

        with open(filepath, 'wb') as file:
            file.seek(0)
            # counter = 0
            file.write(struct.pack(BINARY_HEADER, self.name.encode('utf8'), len(self.triangles)))
            # counter += 1
            for triangle in self.triangles:
                data = [
                    0., 0., 0.,
                    1000*triangle.point1.x, 1000*triangle.point1.y, 1000*triangle.point1.z,
                    1000*triangle.point2.x, 1000*triangle.point2.y, 1000*triangle.point2.z,
                    1000*triangle.point3.x, 1000*triangle.point3.y, 1000*triangle.point3.z,
                    0
                    ]
                file.write(struct.pack(BINARY_FACET, *data))
            file.close()
            
    def to_closed_shell(self):
        return vmf.ClosedShell3D(self.triangles, name=self.name)
    
    def to_open_shell(self):
        return vmf.OpenShell3D(self.triangles, name=self.name)

    
    def extract_points(self):
        points1 = [t.point1 for t in self.triangles]
        points2 = [t.point2 for t in self.triangles]
        points3 = [t.point3 for t in self.triangles]
        
        return list(set(points1 + points2 + points3))
    
    @classmethod
    def from_display_mesh(cls, mesh):
        triangles = []
        for i1, i2, i3 in mesh.triangles:
            triangles.append(vmf.Triangle3D(mesh.points[i1],
                                            mesh.points[i2],
                                            mesh.points[i3]))
        return cls(triangles)