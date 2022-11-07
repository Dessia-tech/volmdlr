#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Gmsh and related objects
"""

from dessia_common import DessiaObject
import volmdlr
import volmdlr.mesh


class GmshParser(DessiaObject):
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
                 physical_names=None,
                 partitioned_entities=None,
                 periodic=None,
                 # ghost_elements = None,
                 parametrizations=None,
                 node_data=None,
                 element_data=None,
                 element_node_data=None,
                 # interpolation_scheme = None,
                 name: str = ''):

        self.mesh_format = mesh_format
        self.physical_names = physical_names
        self.entities = entities
        self.nodes = nodes
        self.elements = elements
        self.partitioned_entities = partitioned_entities
        self.periodic = periodic
        # self.ghost_elements = ghost_elements,
        self.parametrizations = parametrizations
        self.node_data = node_data
        self.element_data = element_data
        self.element_node_data = element_node_data
        # self.interpolation_scheme = interpolation_scheme,
        self.name = name

        DessiaObject.__init__(self, name=name)

    @classmethod
    def from_file(cls, file_path: str):
        """
        defines a gmsh object from .msh file
        """

        file_data = GmshParser.read_file(file_path)
        mesh_format = GmshParser.from_file_mesh_format(file_data['MeshFormat'])
        physical_names = GmshParser.from_file_physical_names(file_data['PhysicalNames'])
        entities = GmshParser.from_file_entities(file_data['Entities'])
        nodes = GmshParser.from_file_nodes(file_data['Nodes'])
        elements = GmshParser.from_file_elements(file_data['Elements'])
        partitioned_entities = GmshParser.from_file_partitioned_entities(file_data['PartitionedEntities'])
        periodic = GmshParser.from_file_periodic(file_data['Periodic'])
        # ghost_elements = GmshParser.from_file_ghost_elements(file_data['GhostElements'])
        parametrizations = GmshParser.from_file_parametrizations(file_data['Parametrizations'])
        node_data = GmshParser.from_file_node_data(file_data['NodeData'])
        element_data = GmshParser.from_file_element_data(file_data['ElementData'])
        element_node_data = GmshParser.from_file_element_node_data(file_data['ElementNodeData'])
        # interpolation_scheme = GmshParser.from_file_interpolation_scheme(file_data['InterpolationScheme'])

        return cls(mesh_format=mesh_format,
                   entities=entities,
                   nodes=nodes,
                   elements=elements,
                   physical_names=physical_names,
                   partitioned_entities=partitioned_entities,
                   periodic=periodic,
                   # ghost_elements=ghost_elements,
                   parametrizations=parametrizations,
                   node_data=node_data,
                   element_data=element_data,
                   element_node_data=element_node_data,
                   # interpolation_scheme=interpolation_scheme
                   name='')

    @staticmethod
    def from_file_elements(lines):
        """
        gets elements data from .msh file
        """

        if not lines:
            return {}

        elements = {}
        elements_types, elements_type = [], []
        step = 1
        while True:
            line = lines[step].split()
            num_elements_in_block = int(line[3])
            if num_elements_in_block:
                element_type = line[2]
                try:
                    elements['elements_type_' + element_type]
                except KeyError:
                    elements_types.append(elements_type)
                    elements['elements_type_' + element_type] = []
                    elements_type = []

                step = step + 1
                elements_list = []
                for i in range(step, step + num_elements_in_block):
                    line = lines[i].split()
                    element = [int(index) - 1 for index in line[1::]]
                    # elements['elements_type_'+ element_type].append(element)
                    elements_list.append(element)
                    elements_type.append(element)
                elements['elements_type_' + element_type].append(elements_list)

            step = i + 1  # num_nodes_in_block

            if step == len(lines):
                break

        return elements

    # $ElementData
    #   numStringTags(ASCII int)
    #   stringTag(string) ...
    #   numRealTags(ASCII int)
    #   realTag(ASCII double) ...
    #   numIntegerTags(ASCII int)
    #   integerTag(ASCII int) ...
    #   elementTag(int) value(double) ...
    #   ...
    # $EndElementData

    @staticmethod
    def from_file_element_data(lines):
        """
        gets mesh element_data from .msh file
        """

        if not lines:
            return {}

        element_data = {}

        return element_data

    # $ElementNodeData
    #   numStringTags(ASCII int)
    #   stringTag(string) ...
    #   numRealTags(ASCII int)
    #   realTag(ASCII double) ...
    #   numIntegerTags(ASCII int)
    #   integerTag(ASCII int) ...
    #   elementTag(int) numNodesPerElement(int) value(double) ...
    #   ...
    # $EndElementNodeData

    @staticmethod
    def from_file_element_node_data(lines):
        """
        gets mesh element_node_data from .msh file
        """

        if not lines:
            return {}

        element_node_data = {}

        return element_node_data

    @staticmethod
    def from_file_entities(lines):
        """
        gets entities data from .msh file
        """

        if not lines:
            return {}

        entities = [int(lines[0].split()[0]), int(lines[0].split()[1]),
                    int(lines[0].split()[2]), int(lines[0].split()[3])]

        points_data = []
        for i in range(1, entities[0] + 1):
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
        for i in range(entities[0] + 1, entities[0] + entities[1] + 1):
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
        for i in range(entities[0] + entities[1] + 1, entities[0] + entities[1] + entities[2] + 1):
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
        for i in range(entities[0] + entities[1] + entities[2] + 1, len(lines)):
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

    # $GhostElements
    #   numGhostElements(size_t)
    #   elementTag(size_t) partitionTag(int)
    #     numGhostPartitions(size_t) ghostPartitionTag(int) ...
    #   ...
    # $EndGhostElements

    @staticmethod
    def from_file_ghost_elements(lines):
        """
        gets mesh ghost_elements from .msh file
        """

        if not lines:
            return {}

        ghost_elements = {}

        return ghost_elements

    # $InterpolationScheme
    #   name(string)
    #   numElementTopologies(ASCII int)
    #   elementTopology
    #   numInterpolationMatrices(ASCII int)
    #   numRows(ASCII int) numColumns(ASCII int) value(ASCII double) ...
    # $EndInterpolationScheme

    @staticmethod
    def from_file_interpolation_scheme(lines):
        """
        gets mesh interpolation_scheme from .msh file
        """

        if not lines:
            return {}

        interpolation_scheme = {}

        return interpolation_scheme

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

        if not lines:
            return {}

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
                    nodes['nodes_dim_' + entity_dim]
                except KeyError:
                    nodes['nodes_dim_' + entity_dim] = []

                step = step + num_nodes_in_block + 1
                points = []
                for i in range(step, step + num_nodes_in_block):
                    line = lines[i].split()
                    points.append(volmdlr.mesh.Node3D(float(line[0]),
                                                      float(line[1]),
                                                      float(line[2])))

                    # nodes['nodes_dim_'+ entity_dim].append(
                    #     volmdlr.Point3D(float(line[0]),
                    #                     float(line[1]),
                    #                     float(line[2])))
                    # nodes_points.append(volmdlr.Point3D(float(line[0]),
                    #                                     float(line[1]),
                    #                                     float(line[2])))

                nodes['nodes_dim_' + entity_dim].append(points)
                nodes_points.extend(points)

                step = step + num_nodes_in_block

            else:
                step += 1

            if len(nodes_points) == int(lines[0].split()[1]):
                break

        if GmshParser.check_2d(nodes_points):
            for key, value in nodes.items():
                values = []
                for points in value:
                    values.append(GmshParser.to_2d(points))
            nodes[key] = values
            nodes['all_nodes'] = GmshParser.to_2d(nodes_points)
        else:
            nodes['all_nodes'] = nodes_points

        return nodes

    # $NodeData
    #   numStringTags(ASCII int)
    #   stringTag(string) ...
    #   numRealTags(ASCII int)
    #   realTag(ASCII double) ...
    #   numIntegerTags(ASCII int)
    #   integerTag(ASCII int) ...
    #   nodeTag(int) value(double) ...
    #   ...
    # $EndNodeData

    @staticmethod
    def from_file_node_data(lines):
        """
        gets mesh node_data from .msh file
        """

        if not lines:
            return {}

        node_data = {}

        return node_data

    # $Parametrizations
    #   numCurveParam(size_t) numSurfaceParam(size_t)
    #   curveTag(int) numNodes(size_t)
    #     nodeX(double) nodeY(double) nodeZ(double) nodeU(double)
    #     ...
    #   ...
    #   surfaceTag(int) numNodes(size_t) numTriangles(size_t)
    #     nodeX(double) nodeY(double) nodeZ(double)
    #       nodeU(double) nodeV(double)
    #       curvMaxX(double) curvMaxY(double) curvMaxZ(double)
    #       curvMinX(double) curvMinY(double) curvMinZ(double)
    #     ...
    #     nodeIndex1(int) nodeIndex2(int) nodeIndex3(int)
    #     ...
    #   ...
    # $EndParametrizations

    @staticmethod
    def from_file_parametrizations(lines):
        """
        gets mesh parametrizations from .msh file
        """

        if not lines:
            return {}

        parametrizations = {}

        return parametrizations

    # $PartitionedEntities
    #   numPartitions(size_t)
    #   numGhostEntities(size_t)
    #     ghostEntityTag(int) partition(int)
    #     ...
    #   numPoints(size_t) numCurves(size_t)
    #     numSurfaces(size_t) numVolumes(size_t)
    #   pointTag(int) parentDim(int) parentTag(int)
    #     numPartitions(size_t) partitionTag(int) ...
    #     X(double) Y(double) Z(double)
    #     numPhysicalTags(size_t) physicalTag(int) ...
    #   ...
    #   curveTag(int) parentDim(int) parentTag(int)
    #     numPartitions(size_t) partitionTag(int) ...
    #     minX(double) minY(double) minZ(double)
    #     maxX(double) maxY(double) maxZ(double)
    #     numPhysicalTags(size_t) physicalTag(int) ...
    #     numBoundingPoints(size_t) pointTag(int) ...
    #   ...
    #   surfaceTag(int) parentDim(int) parentTag(int)
    #     numPartitions(size_t) partitionTag(int) ...
    #     minX(double) minY(double) minZ(double)
    #     maxX(double) maxY(double) maxZ(double)
    #     numPhysicalTags(size_t) physicalTag(int) ...
    #     numBoundingCurves(size_t) curveTag(int) ...
    #   ...
    #   volumeTag(int) parentDim(int) parentTag(int)
    #     numPartitions(size_t) partitionTag(int) ...
    #     minX(double) minY(double) minZ(double)
    #     maxX(double) maxY(double) maxZ(double)
    #     numPhysicalTags(size_t) physicalTag(int) ...
    #     numBoundingSurfaces(size_t) surfaceTag(int) ...
    #   ...
    # $EndPartitionedEntities

    @staticmethod
    def from_file_partitioned_entities(lines):
        """
        gets mesh partitioned_entities from .msh file
        """

        if not lines:
            return {}

        partitioned_entities = {}

        return partitioned_entities

    # $Periodic
    #   numPeriodicLinks(size_t)
    #   entityDim(int) entityTag(int) entityTagMaster(int)
    #   numAffine(size_t) value(double) ...
    #   numCorrespondingNodes(size_t)
    #     nodeTag(size_t) nodeTagMaster(size_t)
    #     ...
    #   ...
    # $EndPeriodic

    @staticmethod
    def from_file_periodic(lines):
        """
        gets mesh periodic from .msh file
        """

        if not lines:
            return {}

        periodic = {}

        return periodic

    @staticmethod
    def from_file_physical_names(lines):
        """
        gets mesh physical_names from .msh file
        """

        if not lines:
            return {}

        physical_names = {}
        for i in range(1, int(lines[0].split()[0]) + 1):
            physical_dim = (lines[i].split()[0])
            try:
                physical_names['physical_dim_' + physical_dim]
            except KeyError:
                physical_names['physical_dim_' + physical_dim] = []

            physical_names['physical_dim_' + physical_dim].append({'physicalTag': int(lines[1][0]),
                                                                   'Name': lines[i][5:-1]})

        return physical_names

    @staticmethod
    def read_file(file_path: str):
        """
        gets lines from a .msh file
        """

        data = {'MeshFormat': [],
                'PhysicalNames': [],
                'Entities': [],
                'PartitionedEntities': [],
                'Nodes': [],
                'Elements': [],
                'Periodic': [],
                'GhostElements': [],
                'Parametrizations': [],
                'NodeData': [],
                'ElementData': [],
                'ElementNodeData': [],
                'InterpolationScheme': []}

        f = open(file_path, "r", encoding="utf-8")
        # lines = []
        while True:
            line = f.readline().strip()
            if not line:
                break
            if line[0] == '$':
                data_type = line[1::]
                line = f.readline().strip()
                while line[0:4] != '$End':
                    # lines.append(line)
                    data[data_type].append(line)
                    line = f.readline().strip()
                # data[data_type] = lines
                # lines = []

        return data

    def define_tetrahedron_element_mesh(self):
        """
        defines a volmdlr mesh with TetrahedronElement from a .msh file
        """

        # nodes = self.nodes[0]
        points = self.nodes['all_nodes']
        # elements = self.elements[0]

        tetrahedron_elements = self.elements['elements_type_4']
        element_groups = []
        for tetrahedrons in tetrahedron_elements:
            tetrahedrons_mesh = []
            for tetrahedron in tetrahedrons:
                tetrahedrons_mesh.append(volmdlr.mesh.TetrahedralElement([points[tetrahedron[0]],
                                                                          points[tetrahedron[1]],
                                                                          points[tetrahedron[2]],
                                                                          points[tetrahedron[3]]]))

            element_groups.append(volmdlr.mesh.ElementsGroup(tetrahedrons_mesh, name=''))

        mesh = volmdlr.mesh.Mesh(element_groups)
        # mesh.nodes = points #gmsh points are duplicated > not needed
        # mesh.node_to_index = {mesh.nodes[i]: i for i in range(len(mesh.nodes))}

        return mesh

    def define_triangular_element_mesh(self):
        """
        defines a volmdlr mesh with TriangularElement from a .msh file
        """

        # nodes = self.nodes[0]
        points = self.nodes['all_nodes']
        # elements = self.elements[0]

        triangles_elements = self.elements['elements_type_2']
        element_groups = []
        for triangles in triangles_elements:
            triangles_mesh = []
            for triangle in triangles:
                if points[0].__class__.__name__[-2::] == '3D':
                    triangles_mesh.append(volmdlr.mesh.TriangularElement3D([points[triangle[0]],
                                                                            points[triangle[1]],
                                                                            points[triangle[2]]]))
                else:
                    triangles_mesh.append(volmdlr.mesh.TriangularElement2D([points[triangle[0]],
                                                                            points[triangle[1]],
                                                                            points[triangle[2]]]))
            element_groups.append(volmdlr.mesh.ElementsGroup(triangles_mesh, name=''))

        mesh = volmdlr.mesh.Mesh(element_groups)
        # mesh.nodes = points #gmsh points are duplicated > not needed
        # mesh.node_to_index = {mesh.nodes[i]: i for i in range(len(mesh.nodes))}

        return mesh

    @staticmethod
    def check_2d(list_nodes):
        for node in list_nodes:
            if node[2] != 0:
                return False
        return True

    @staticmethod
    def to_2d(list_nodes):
        return [volmdlr.mesh.Node2D(node[0], node[1]) for
                node in list_nodes]
