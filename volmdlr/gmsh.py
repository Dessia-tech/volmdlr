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
                 mesh_format,
                 physical_name,
                 entities,
                 partitioned_entities,
                 nodes,
                 elements,
                 periodic,
                 ghost_elements,
                 parametrizations,
                 node_data,
                 element_data,
                 element_node_data,
                 name: str):

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

    def from_file(cls, file_path: str):
        """
        defines a gmsh object from .msh file
        """

        file_data = Gmsh.read_file(file_path)
        mesh_format = Gmsh.from_file_mesh_format(file_data['MeshFormat'])
        physical_name = Gmsh.from_file_physical_name(file_data['PhysicalName'])
        entities = Gmsh.from_file_entities(file_data['Entities'])
        nodes = Gmsh.from_file_nodes(file_data['Nodes'])
        elements = Gmsh.from_file_elements(file_data['Elements'])

        return cls(mesh_format,
                   physical_name,
                   entities,
                   nodes,
                   elements)

