#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gmsh and related objects
"""

from dessia_common import DessiaObject
import volmdlr
import volmdlr.mesh

class Gmsh(DessiaObject):
    _standalone_in_db = False
    _non_serializable_attributes = []
    _non_eq_attributes = ['name']
    _non_hash_attributes = ['name']
    _generic_eq = True

    def __init__(self,
                 mesh_format: dict,
                 nodes: dict,
                 elements: dict,
                 entities: dict,
                 physical_name = None,
                 partitioned_entities = None,
                 periodic = None,
                 ghost_elements = None,
                 parametrizations = None,
                 node_data = None,
                 element_data = None,
                 element_node_data = None,
                 name: str = ''):

        self.mesh_format = mesh_format
        self.physical_name = physical_name,
        self.entities = entities,
        self.nodes = nodes,
        self.elements = elements,
        self.partitioned_entities = partitioned_entities,
        self.periodic = periodic,
        self.ghost_elements = ghost_elements,
        self.parametrizations = parametrizations,
        self.node_data = node_data,
        self.element_data = element_data,
        self.element_node_data = element_node_data,
        self.name = name

        DessiaObject.__init__(self, name=name)

    @classmethod
    def from_file(cls, file_path: str):
        """
        defines a gmsh object from .msh file
        """

        file_data = Gmsh.read_file(file_path)
        mesh_format = Gmsh.from_file_mesh_format(file_data['MeshFormat'])
        # physical_name = Gmsh.from_file_physical_name(file_data['PhysicalName'])
        entities = Gmsh.from_file_entities(file_data['Entities'])
        nodes = Gmsh.from_file_nodes(file_data['Nodes'])
        elements = Gmsh.from_file_elements(file_data['Elements'])

        return cls(mesh_format=mesh_format,
                    entities=entities,
                    nodes=nodes,
                    elements=elements)

    @staticmethod
    def from_file_elements(lines):
        """
        gets elements data from .msh file
        """

        elements = {}
        elements_types, elements_type = [], []
        step = 1
        while True:
            line = lines[step].split()
            num_elements_in_block = int(line[3])
            if num_elements_in_block:
                element_type = line[2]
                try:
                    elements['elements_type_'+ element_type]
                except KeyError:
                    elements_types.append(elements_type)
                    elements['elements_type_'+ element_type] = []
                    elements_type = []

                step = step +1
                elements_list = []
                for i in range(step, step+num_elements_in_block):
                    line = lines[i].split()
                    element = [int(index)-1 for index in line[1::]]
                    # elements['elements_type_'+ element_type].append(element)
                    elements_list.append(element)
                    elements_type.append(element)
                elements['elements_type_'+ element_type].append(elements_list)

            step = i+1 #num_nodes_in_block

            if step == len(lines):
                break

        return elements

    @staticmethod
    def from_file_entities(lines):
        """
        gets entities data from .msh file
        """

        entities = [int(lines[0].split()[0]), int(lines[0].split()[1]),
                    int(lines[0].split()[2]), int(lines[0].split()[3])]

        points_data = []
        for i in range(1, entities[0]+1):
            line = lines[i].split()
            points = {}
            points['pointTag'] = int(line[0])
            points['X'] = float(line[1])
            points['Y'] = float(line[2])
            points['Z'] = float(line[3])
            points['numPhysicalTags'] = int(line[4])
            if points['numPhysicalTags']:
                points['physicalTag'] = int(line[5])
            points_data.append(points)

        curves_data = []
        for i in range(entities[0]+1, entities[0]+entities[1]+1):
            line = lines[i].split()
            curves = {}
            curves['curveTag'] = int(line[0])
            curves['minX'] = float(line[1])
            curves['minY'] = float(line[2])
            curves['minZ'] = float(line[3])
            curves['maxX'] = float(line[4])
            curves['maxY'] = float(line[5])
            curves['maxZ'] = float(line[6])
            curves['numPhysicalTags'] = int(line[7])
            if curves['numPhysicalTags']:
                curves['physicalTag'] = int(line[8])
                curves['numBoundingPoints'] = int(line[9])
                curves['pointTag'] = [float(l) for l in line[10::]]
            else:
                curves['numBoundingPoints'] = int(line[8])
                curves['pointTag'] = [float(l) for l in line[9::]]
            curves_data.append(curves)

        surfaces_data = []
        for i in range(entities[0]+entities[1]+1, entities[0]+entities[1]+entities[2]+1):
            line = lines[i].split()
            surfaces = {}
            surfaces['surfaceTag'] = int(line[0])
            surfaces['minX'] = float(line[1])
            surfaces['minY'] = float(line[2])
            surfaces['minZ'] = float(line[3])
            surfaces['maxX'] = float(line[4])
            surfaces['maxY'] = float(line[5])
            surfaces['maxZ'] = float(line[6])
            surfaces['numPhysicalTags'] = int(line[7])
            if surfaces['numPhysicalTags']:
                surfaces['physicalTag'] = int(line[8])
                surfaces['numBoundingCurves'] = int(line[9])
                surfaces['curveTag'] = [float(l) for l in line[10::]]
            else:
                surfaces['numBoundingCurves'] = float(line[8])
                surfaces['curveTag'] = [float(l) for l in line[9::]]
            surfaces_data.append(surfaces)

        volumes_data = []
        for i in range(entities[0]+entities[1]+entities[2]+1, len(lines)):
            line = lines[i].split()
            volumes = {}
            volumes['volumeTag'] = int(line[0])
            volumes['minX'] = float(line[1])
            volumes['minY'] = float(line[2])
            volumes['minZ'] = float(line[3])
            volumes['maxX'] = float(line[4])
            volumes['maxY'] = float(line[5])
            volumes['maxZ'] = float(line[6])
            volumes['numPhysicalTags'] = int(line[7])
            if volumes['numPhysicalTags']:
                volumes['physicalTag'] = int(line[8])
                volumes['numBoundngSurfaces'] = int(line[9])
                volumes['surfaceTag'] = [float(l) for l in line[10::]]
            else:
                volumes['numBoundngSurfaces'] = float(line[8])
                volumes['surfaceTag'] = [float(l) for l in line[9::]]
            volumes_data.append(volumes)

        return {'points': points_data,
                'curves': curves_data,
                'surfaces': surfaces_data,
                'volumes': volumes_data}

    @staticmethod
    def from_file_mesh_format(lines):
        """
        gets mesh format data from .msh file
        """

        mesh_format = {}
        mesh_format['version_number'] = float(lines[0].split()[0])
        if int(lines[0].split()[1]) == 0:
            mesh_format['file_type'] = (int(lines[0].split()[1]), 'ASCII')
        else:
            mesh_format['file_type'] = (int(lines[0].split()[1]), 'Binary')
        mesh_format['data_size'] = int(lines[0].split()[2])

        return mesh_format

    @staticmethod
    def from_file_nodes(lines):
        """
        gets mesh nodes from .msh file
        """

        nodes = {}
        nodes_points = []
        step = 1
        while True:
            line = lines[step].split()
            num_nodes_in_block = int(line[3])
            if num_nodes_in_block:
                entity_dim = line[0]
                # nodes['nodes_dim_'+ entity_dim] = []
                try:
                    nodes['nodes_dim_'+ entity_dim]
                except KeyError:
                    nodes['nodes_dim_'+ entity_dim] = []

                step = step + num_nodes_in_block+1
                points = []
                for i in range(step, step+num_nodes_in_block):
                    line = lines[i].split()
                    points.append(volmdlr.Point3D(float(line[0]),
                                                  float(line[1]),
                                                  float(line[2])))

                    # nodes['nodes_dim_'+ entity_dim].append(
                    #     volmdlr.Point3D(float(line[0]),
                    #                     float(line[1]),
                    #                     float(line[2])))
                    # nodes_points.append(volmdlr.Point3D(float(line[0]),
                    #                                     float(line[1]),
                    #                                     float(line[2])))

                nodes['nodes_dim_'+ entity_dim].append(points)
                nodes_points.extend(points)

                step = step+num_nodes_in_block

            else:
                step += 1

            if len(nodes_points) == int(lines[0].split()[1]):
                break

        nodes['all_nodes'] = nodes_points

        return nodes

    @staticmethod
    def read_file(file_path: str):
        """
        gets lines from a .msh file
        """

        f = open(file_path, "r")
        lines = []
        data = {}
        while(True):
            line = f.readline().strip()
            if not line:
                break
            if line[0] == '$':
                data_type = line[1::]
                line = f.readline().strip()
                while line[0:4] != '$End':
                    lines.append(line)
                    line = f.readline().strip()
                data[data_type] = lines
                lines = []

        return data

    def define_mesh(self):
        """
        defines a mesh from a .msh file
        """

        nodes = self.nodes
        points = nodes['all_nodes']
        elements = self.elements

        triangles_elements = elements['elements_type_2']
        triangles_mesh, element_groups = [], []
        for triangles in triangles_elements:
            for triangle in triangles:
                if points[0].__class__.__name__[-2::] == '3D':
                    triangular_element = volmdlr.mesh.TriangularElement3D([points[triangle[0]],
                                                                           points[triangle[1]],
                                                                           points[triangle[2]]])
                else:
                    triangular_element = volmdlr.mesh.TriangularElement2D([points[triangle[0]],
                                                                           points[triangle[1]]])
                triangles_mesh.append(triangular_element)

            element_groups.append(volmdlr.mesh.ElementsGroup(triangles_mesh, name=''))

        return volmdlr.mesh.Mesh(element_groups)
