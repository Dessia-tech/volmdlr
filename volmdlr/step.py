#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ISO STEP reader/writer.
"""

import time
from typing import List
from collections import deque
from dataclasses import dataclass, field

import matplotlib.pyplot as plt
import networkx as nx
import plot_data.graph

import dessia_common.core as dc  # isort: skip
from dessia_common.files import BinaryFile  # isort: skip

import volmdlr
import volmdlr.core
import volmdlr.edges
import volmdlr.faces
import volmdlr.curves
import volmdlr.primitives3d
import volmdlr.wires
from volmdlr.utils import step_reader
from volmdlr.utils.step_reader import STEP_TO_VOLMDLR, STEP_REPRESENTATION_ENTITIES


class StepFunction(dc.DessiaObject):
    """
    Abstract class defining a step function.

    """

    def __init__(self, function_id, function_name, function_arg):
        self.id = function_id
        self.name = function_name
        self.arg = function_arg

        # TODO : modify this continuation and simplify
        if self.name == "":
            if self.arg[1][0] == 'B_SPLINE_SURFACE':
                self.simplify('B_SPLINE_SURFACE')
            if self.arg[1][0] == 'B_SPLINE_CURVE':
                self.simplify('B_SPLINE_CURVE')
        dc.DessiaObject.__init__(self, name=function_name)

    def simplify(self, new_name):
        """ADD DOCSTRING."""
        # ITERATE ON SUBFUNCTIONS
        args = [subfun[1] for (i, subfun) in enumerate(self.arg) if
                (len(subfun[1]) != 0 or i == 0)]
        arguments = []
        for arg in args:
            if not arg:
                arguments.append("''")
            else:
                arguments.extend(arg)
        arguments.pop()  # DELETE REPRESENTATION_ITEM('')

        self.name = new_name
        self.arg = arguments


class Step(dc.DessiaObject):
    """
    Defines the Step class.

    """

    _standalone_in_db = True

    def __init__(self, lines: List[str], name: str = ''):
        self.functions, self.connections = self.read_lines(lines)
        self._graph = None
        self.global_uncertainty = 1e-6
        self.length_conversion_factor = 1
        self.angle_conversion_factor = 1
        # self.read_diagnostic = StepReaderReport
        self._roots_nodes = None

        dc.DessiaObject.__init__(self, name=name)

    @property
    def all_connections(self):
        """Returns all pairs of connections."""
        list_connections = []
        for key, values in self.connections.items():
            for value in values:
                list_connections.append((key, value))
        return list_connections

    @property
    def root_nodes(self):
        """Returns a dictionary containing the nodes of the step file function that are used as start points."""
        if not self._roots_nodes:
            self._roots_nodes = self.get_root_nodes()
        return self._roots_nodes

    def graph(self):
        """Returns the step file networkx graph of dependencies."""
        if not self._graph:
            self._graph = self.create_graph()
        return self._graph

    @classmethod
    def from_stream(cls, stream: BinaryFile, name: str = ''):
        """Instantiate a Step object from a stream."""
        stream.seek(0)
        lines = []
        for line in stream:
            line = line.decode("ISO-8859-1")
            line = line.replace("\r", "")
            lines.append(line)
        return cls(lines, name=name)

    @classmethod
    def from_file(cls, filepath: str = None, name: str = ''):
        """Instantiate a Step object from a step file."""
        with open(filepath, "r", encoding="ISO-8859-1") as file:
            lines = []
            for line in file:
                lines.append(line)
        return cls(lines, name=name)

    def read_lines(self, lines):
        """Translate the step file into step functions objects."""
        dict_connections = {}
        previous_line = ""
        functions = {}
        for line in lines:
            # line = line.replace(" ", "")
            line = line.replace("\n", "")

            # SKIP EMPTY LINE
            if not line:
                continue

            # ASSEMBLE LINES IF THEY ARE SEPARATED
            if line[-1] != ';':
                previous_line = previous_line + line
                continue

            line = previous_line + line

            # SKIP HEADER
            if line[0] != "#":
                previous_line = str()
                continue

            function = line.split("=", maxsplit=1)
            function_id = int(function[0][1:].strip())
            function_name_arg = function[1].split("(", 1)
            function_name = function_name_arg[0].replace(" ", "")
            start_index_name = function_name_arg[1].find("'")
            if start_index_name != -1:
                end_index_name = function_name_arg[1].find("'", start_index_name + 1)
                if end_index_name != -1:
                    function_arg_string = function_name_arg[1][end_index_name + 1:]
                else:
                    function_arg_string = function_name_arg[1]
            else:
                function_arg_string = function_name_arg[1]
            function_arg = function_arg_string.split("#")
            connections = []
            for connec in function_arg[1:]:
                connec = connec.split(",")
                connec = connec[0].split(")")
                if connec[0][-1] != "'":
                    function_connection = int(connec[0])
                    connections.append(function_connection)

            previous_line = str()

            # FUNCTION ARGUMENTS
            functions, connections = self._helper_intantiate_step_functions(functions, connections,
                                                                            [function_id, function_name,
                                                                             function_name_arg])

            dict_connections[function_id] = connections

        return functions, dict_connections

    def _helper_intantiate_step_functions(self, functions, connections, function_parameters):
        """Helper function to read_lines."""
        function_id, function_name, function_name_arg = function_parameters
        function_arg = function_name_arg[1]
        arguments = step_reader.step_split_arguments(function_arg)
        new_name = ''
        new_arguments = []
        if function_name == "":
            name_arg = self.step_subfunctions(arguments)
            for name, arg in name_arg:
                new_name += name + ', '
                new_arguments.extend(arg)
            new_name = new_name[:-2]
            function_name = new_name
            arguments = new_arguments
            for arg in arguments:
                if arg[0] == '#':
                    connections.append(int(arg[1:]))

        for i, argument in enumerate(arguments):
            if argument[:2] == '(#' and argument[-1] == ')':
                arg_list = step_reader.set_to_list(argument)
                arguments[i] = arg_list

        function = StepFunction(function_id, function_name, arguments)
        functions[function_id] = function
        return functions, connections

    def not_implemented(self):
        not_implemented = []
        for _, fun in self.functions.items():
            if fun.name not in STEP_TO_VOLMDLR:
                not_implemented.append(fun.name)
        return list(set(not_implemented))

    def create_graph(self):
        """
        Step functions graph.

        :return: A graph representation the step file structure.
        :rtype: networkx.DiGraph
        """
        graph = nx.DiGraph()
        labels = {}

        for function in self.functions.values():
            if function.name == 'SHAPE_REPRESENTATION_RELATIONSHIP':
                # Create short cut from id1 to id2
                id1 = int(function.arg[2][1:])
                id2 = int(function.arg[3][1:])
                elem1 = (function.id, id1)
                elem2 = (function.id, id2)
                self.all_connections.remove(elem1)
                self.all_connections.remove(elem2)
                self.all_connections.append((elem1[1], elem2[1]))

                self.functions[id1].arg.append(f'#{id2}')

            elif function.name in STEP_TO_VOLMDLR:
                graph.add_node(function.id,
                               color='rgb(0, 0, 0)',
                               shape='.',
                               name=str(function.id))
                labels[function.id] = str(function.id) + ' ' + function.name

        # Delete connection if node not found
        node_list = list(graph.nodes())
        delete_connection = []
        for connection in self.all_connections:
            if connection[0] not in node_list \
                    or connection[1] not in node_list:
                delete_connection.append(connection)
        for delete in delete_connection:
            self.all_connections.remove(delete)

        # Create graph connections
        graph.add_edges_from(self.all_connections)

        # Remove single nodes
        delete_nodes = []
        for node in graph.nodes:
            if graph.degree(node) == 0:
                delete_nodes.append(node)
        for node in delete_nodes:
            graph.remove_node(node)
            # G.remove_node(node)
        return graph

    def draw_graph(self, graph=None, reduced=False):
        """
        Draw a graph for Step data.

        :param graph: DESCRIPTION, defaults to None
        :type graph: TYPE, optional
        :param reduced: DESCRIPTION, defaults to False
        :type reduced: TYPE, optional
        :return: DESCRIPTION
        :rtype: TYPE

        """

        delete = ['CARTESIAN_POINT', 'DIRECTION']
        if graph is None:
            new_graph = self.create_graph()
        else:
            new_graph = graph.copy()

        labels = {}
        for id_nb, function in self.functions.items():
            if id_nb in new_graph.nodes and not reduced:
                labels[id_nb] = str(id_nb) + ' ' + function.name
            elif id_nb in new_graph.nodes and reduced:
                if function.name not in delete:
                    labels[id_nb] = str(id_nb) + ' ' + function.name
                else:
                    new_graph.remove_node(id_nb)
        pos = nx.kamada_kawai_layout(new_graph)
        plt.figure()
        nx.draw_networkx_nodes(new_graph, pos)
        nx.draw_networkx_edges(new_graph, pos)
        nx.draw_networkx_labels(new_graph, pos, labels)

    @staticmethod
    def step_subfunctions(subfunctions):
        """Handles context elements from step file."""
        subfunctions = subfunctions[0].replace(" ", "")
        parenthesis_count = 0
        subfunction_names = []
        subfunction_args = []
        subfunction_name = ""
        subfunction_arg = ""
        for char in subfunctions:

            if char == "(":
                parenthesis_count += 1
                if parenthesis_count == 1:
                    subfunction_names.append(subfunction_name)
                    subfunction_name = ""
                else:
                    subfunction_arg += char

            elif char == ")":
                parenthesis_count -= 1
                if parenthesis_count == 0:
                    subfunction_args.append(subfunction_arg)
                    subfunction_arg = ""
                else:
                    subfunction_arg += char

            elif parenthesis_count == 0:
                subfunction_name += char

            else:
                subfunction_arg += char
        return [
            (subfunction_names[i], step_reader.step_split_arguments(subfunction_args[i]))
            for i in range(len(subfunction_names))]

    def parse_arguments(self, arguments):
        """Converts the arguments IDs from string to integer."""
        for i, arg in enumerate(arguments):
            if isinstance(arg, str) and arg[0] == '#':
                arguments[i] = int(arg[1:])
            elif isinstance(arg, str) and arg[0:2] == '(#':
                argument = []
                arg_id = ""
                for char in arg[1:-1]:
                    if char == ',':
                        argument.append(arg_id)
                        arg_id = ""
                        continue

                    arg_id += char
                argument.append(arg_id)
                arguments[i] = argument

    def instantiate(self, name, arguments, object_dict, step_id):
        """
        Gives the volmdlr object related to the step function.
        """
        self.parse_arguments(arguments)
        fun_name = name.replace(', ', '_')
        fun_name = fun_name.lower()
        try:
            if hasattr(step_reader, fun_name):
                volmdlr_object = getattr(step_reader, fun_name)(arguments, object_dict)

            elif name in STEP_TO_VOLMDLR and hasattr(STEP_TO_VOLMDLR[name], "from_step"):
                volmdlr_object = STEP_TO_VOLMDLR[name].from_step(
                    arguments, object_dict, name=name, step_id=step_id, global_uncertainty=self.global_uncertainty,
                    length_conversion_factor=self.length_conversion_factor,
                    angle_conversion_factor=self.angle_conversion_factor)

            else:
                raise NotImplementedError(f'Dont know how to interpret #{step_id} = {name}({arguments})')
        except (ValueError, NotImplementedError) as error:
            raise ValueError(f"Error while instantiating #{step_id} = {name}({arguments})") from error
        return volmdlr_object

    def create_node_list(self, initial_nodes):
        """
        Step functions graph as a list of nodes.

        :param initial_nodes: Initial list of shell nodes and assemblies entities.
        :type initial_nodes: List[int]
        :return: A list of nodes in the right order of dependency.
        :rtype: List[int]
        """
        list_head = []
        list_nodes = []
        visited_set = set()
        stack = deque(initial_nodes)
        while stack:
            node = stack.popleft()
            name = self.functions[node].name
            if node not in visited_set and name in STEP_TO_VOLMDLR:
                visited_set.add(node)
                if self.connections[node]:
                    list_nodes.append(node)
                    for connection in self.connections[node]:
                        if connection not in visited_set:
                            stack.append(connection)
                else:
                    # Entities without connections should be instantiated first
                    list_head.append(node)
        return list_head + list_nodes[::-1]

    def get_shell_node_from_representation_entity(self, id_representation_entity: int):
        """
        Find the shell node ID related to a given representation entity.

        :param id_representation_entity: Representation entity ID.
        :type id_representation_entity: int
        :return: Shell ID.
        :rtype: int
        """
        name_representation_entity = self.functions[id_representation_entity].name
        arg = self.functions[id_representation_entity].arg[1]
        if name_representation_entity == "MANIFOLD_SURFACE_SHAPE_REPRESENTATION":
            if self.functions[int(arg[0][1:])].name == "AXIS2_PLACEMENT_3D":
                id_solid_entity = int(arg[1][1:])
            else:
                id_solid_entity = int(arg[0][1:])
            if self.functions[id_solid_entity].name in {"CLOSED_SHELL", "OPEN_SHELL"}:
                return id_solid_entity
            id_shell = self.functions[id_solid_entity].arg[1]
            if isinstance(id_shell, list):
                return int(id_shell[0][1:])
            return int(id_shell[1:])
        if self.functions[int(arg[0][1:])].name == "AXIS2_PLACEMENT_3D":
            id_solid_entity = int(arg[1][1:])
        else:
            id_solid_entity = int(arg[0][1:])
        id_shell = self.functions[id_solid_entity].arg[1]
        if isinstance(id_shell, list):
            return int(id_shell[0][1:])
        return int(id_shell[1:])

    def get_shell_node_from_shape_representation(self, id_shape_representation: int):
        """
        Find the shell node ID related to a given shape representation.

        :param id_shape_representation: Representation entity ID.
        :type id_shape_representation: int
        :return: Shell ID.
        :rtype: int
        """
        if len(self.functions[id_shape_representation].arg) < 4:
            # From the step file, the SHAPE_REPRESENTATION entity has 3 arguments. But we add a 4th argument to
            # those SHAPE_REPRESENTATION entity that are related to a representation entity. So, if the length of arg
            # is less of 4 there is no representation entity related to it, and we return None.
            return None
        id_representation_entity = int(self.functions[id_shape_representation].arg[3][1:])
        id_solid_entity = int(self.functions[id_representation_entity].arg[1][0][1:])
        id_shell = self.functions[id_solid_entity].arg[1]
        if isinstance(id_shell, list):
            return int(id_shell[0][1:])
        return int(id_shell[1:])

    def get_frame_mapped_shell_node(self, node: int):
        """
        Find the shell node in the assembly.

        :param node: Assembly step entity node.
        :type node: int
        :return: Shell ID.
        :rtype: int
        """
        id_representation_entity = None
        arguments = self.functions[node].arg
        name_arg1 = self.functions[int(arguments[2][1:])].name
        name_arg2 = self.functions[int(arguments[3][1:])].name
        if name_arg1 in STEP_REPRESENTATION_ENTITIES:
            id_representation_entity = int(arguments[2][1:])
        elif name_arg2 in STEP_REPRESENTATION_ENTITIES:
            id_representation_entity = int(arguments[3][1:])
        if id_representation_entity:
            return self.get_shell_node_from_representation_entity(id_representation_entity)
        id_shape_representation = int(arguments[3][1:])
        if len(self.functions[id_shape_representation].arg) < 4:
            id_shape_representation = int(arguments[2][1:])
        if self.functions[id_shape_representation].name == "SHAPE_REPRESENTATION":
            return self.get_shell_node_from_shape_representation(id_shape_representation)
        id_representation_entity = int(self.functions[id_shape_representation].arg[1][1][1:])
        id_shell = self.functions[id_representation_entity].arg[1]
        if isinstance(id_shell, list):
            return int(id_shell[0][1:])
        return int(id_shell[1:])

    def shape_definition_representation_to_shell_node(self, shape_definition_representation_id):
        """Returns the ID of the shell entity related to the given shape_definition_representation ID."""
        id_representation_entity = self.functions[shape_definition_representation_id].arg[1]
        function_name = self.functions[int(id_representation_entity[1:])].name
        if function_name in STEP_REPRESENTATION_ENTITIES:
            return self.get_shell_node_from_representation_entity(int(id_representation_entity[1:]))
        if function_name == "SHAPE_REPRESENTATION":
            return self.get_shell_node_from_shape_representation(int(id_representation_entity[1:]))
        raise NotImplementedError("This shape representation has not yet been recognized by volmdlr step translator.")

    def product_definition_to_product(self, id_product_definition):
        """Returns the ID of the product entity related to the given product_definition ID."""
        if self.functions[id_product_definition].name == "NEXT_ASSEMBLY_USAGE_OCCURRENCE":
            id_product_definition = int(self.functions[id_product_definition].arg[3][1:])
        id_product_definition_formation = self.functions[id_product_definition].arg[2]
        id_product = self.functions[int(id_product_definition_formation[1:])].arg[2]
        return int(id_product[1:])

    def shape_definition_representation_to_product_node(self, shape_definition_representation_id):
        """Returns the ID of the product entity related to the given shape_definition_representation ID."""
        id_product_definition_shape = self.functions[shape_definition_representation_id].arg[0]
        id_product_definition = int(self.functions[int(id_product_definition_shape[1:])].arg[2][1:])
        return self.product_definition_to_product(id_product_definition)

    def get_root_nodes(self):
        """Returns a dictionary containing the nodes of the step file function that are used as start points."""
        next_assembly_usage_occurrence = []
        product_definitions = []
        shape_representation_relationship = []
        shape_representations = []
        shape_definition_representation = []
        shell_nodes = []
        geometric_representation_context = {}
        not_shell_nodes = []
        context_dependent_shape_representation = []
        for function in self.functions.values():
            if function.name == "NEXT_ASSEMBLY_USAGE_OCCURRENCE":
                next_assembly_usage_occurrence.append(function.id)
            elif function.name == "PRODUCT_DEFINITION":
                product_definitions.append(function.id)
            elif function.name == "SHAPE_REPRESENTATION_RELATIONSHIP":
                shape_representation_relationship.append(function.id)
                id_shape_representation = int(function.arg[3][1:])
                shape_representations.append(id_shape_representation)
            elif function.name == "SHAPE_DEFINITION_REPRESENTATION":
                id_shape_representation = int(function.arg[1][1:])
                shape_representations.append(id_shape_representation)
                id_geometric_context = int(self.functions[id_shape_representation].arg[-1][1:])
                geometric_representation_context[id_shape_representation] = id_geometric_context
                shape_definition_representation.append(function.id)
            elif function.name in {"CLOSED_SHELL", "OPEN_SHELL"}:
                shell_nodes.append(function.id)
            elif function.name == 'BREP_WITH_VOIDS':
                shell_nodes.append(function.id)
                not_shell_nodes.append(int(function.arg[1][1:]))
            elif function.name == "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION":
                context_dependent_shape_representation.append(function.id)
        for node in not_shell_nodes:
            shell_nodes.remove(node)
        return {"NEXT_ASSEMBLY_USAGE_OCCURRENCE": next_assembly_usage_occurrence,
                "CONTEXT_DEPENDENT_SHAPE_REPRESENTATION": context_dependent_shape_representation,
                "PRODUCT_DEFINITION": product_definitions,
                "SHAPE_REPRESENTATION_RELATIONSHIP": shape_representation_relationship,
                "SHAPE_REPRESENTATION": shape_representations,
                "SHAPE_DEFINITION_REPRESENTATION": shape_definition_representation,
                "GEOMETRIC_REPRESENTATION_CONTEXT": geometric_representation_context,
                "SHELLS": shell_nodes}

    def get_assembly_structure(self):
        """
        Get assembly dependency structure.
        """
        assemblies_structure = {}
        assemblies = set()
        shapes = set()
        for node in self.root_nodes["NEXT_ASSEMBLY_USAGE_OCCURRENCE"]:
            function = self.functions[node]
            assembly_product_definition = int(function.arg[3][1:])
            assembly_node = int(self.functions[assembly_product_definition].arg[4][1:])
            assemblies_structure.setdefault(assembly_node, []).append(node)
            assemblies.add(assembly_node)
            id_product_definition = int(function.arg[4][1:])
            if len(self.functions[id_product_definition].arg) > 5:
                for arg in self.functions[id_product_definition].arg[5:]:
                    shapes.add(int(arg[1:]))
        valid_entities = assemblies.union(shapes)
        return assemblies_structure, valid_entities

    def get_assembly_data(self, assembly_usage_occurence, valid_entities, assembly_frame, object_dict):
        """
        Helper function to get assembly data.
        """
        assembly_shapes = []
        assembly_positions = []
        for node in assembly_usage_occurence:
            function = self.functions[node]
            id_product_definition = int(function.arg[4][1:])
            ids_shape_definition_representation = [int(arg[1:]) for
                                                   arg in self.functions[id_product_definition].arg[4:]
                                                   if int(arg[1:]) in valid_entities]
            assembly_shapes.extend(ids_shape_definition_representation)
            id_context_dependent_shape_representation = int(function.arg[-1][1:])
            id_transformation = int(self.functions[id_context_dependent_shape_representation].arg[0][1:])
            id_item_defined_transformation = int(self.functions[id_transformation].arg[4][1:])
            id_frame1 = int(self.functions[id_item_defined_transformation].arg[2][1:])
            id_frame2 = int(self.functions[id_item_defined_transformation].arg[3][1:])
            frame1 = object_dict[id_frame1]
            # frame2 = object_dict[id_frame2]
            if frame1 == assembly_frame:
                component_frame = id_frame2
            else:
                component_frame = id_frame1
            positions = [component_frame] * len(ids_shape_definition_representation)
            assembly_positions.extend(positions)
        return assembly_shapes, assembly_positions

    def context_dependent_shape_representation_to_next_assembly_usage_occurrence(self, node):
        """
        Returns id of the next_assembly_usage_occurrence related to the given context_dependent_shape_representation.
        """
        arg = self.functions[node].arg
        id_product_definition_shape = int(arg[1][1:])
        return int(self.functions[id_product_definition_shape].arg[2][1:])

    def create_connections(self):
        """
        Create connections between step entities.
        """
        for node in self.root_nodes['SHAPE_REPRESENTATION_RELATIONSHIP']:
            # Associate each step representation entity to its SHAPE_REPRESENTATION
            function = self.functions[node]
            id_shape_representation = int(function.arg[2][1:])
            id_shape = int(function.arg[3][1:])
            self.connections[id_shape_representation].append(id_shape)
            self.functions[id_shape_representation].arg.append(f'#{id_shape}')
        for node in self.root_nodes['SHAPE_DEFINITION_REPRESENTATION']:
            # Associate each step representation entity to its SHAPE_REPRESENTATION
            function = self.functions[node]
            id_product_definition_shape = int(function.arg[0][1:])
            id_product_definition = int(self.functions[id_product_definition_shape].arg[2][1:])
            id_shape_representation = int(function.arg[1][1:])
            self.connections[id_product_definition].append(node)
            self.functions[id_product_definition].arg.append(f'#{node}')
            if self.functions[id_shape_representation].name == "SHAPE_REPRESENTATION" and \
                    len(self.functions[id_shape_representation].arg) >= 4:
                # todo: take all the "arg" starting from index 3 to end ??? needs investigation
                id_shapes = [int(arg[1:]) for arg in self.functions[id_shape_representation].arg[3:]]
                self.connections[id_product_definition].extend(id_shapes)
                for id_shape in id_shapes:
                    self.functions[id_product_definition].arg.append(f'#{id_shape}')
            elif self.functions[id_shape_representation].name in STEP_REPRESENTATION_ENTITIES:
                self.connections[id_product_definition].append(id_shape_representation)
                self.functions[id_product_definition].arg.append(f'#{id_shape_representation}')

            shell_node = self.shape_definition_representation_to_shell_node(node)
            product_node = self.shape_definition_representation_to_product_node(node)
            if shell_node:
                self.connections[shell_node].append(product_node)
                self.functions[shell_node].arg.append(f'#{product_node}')

        for node in self.root_nodes['CONTEXT_DEPENDENT_SHAPE_REPRESENTATION']:
            next_assembly_usage_occurrence = \
                self.context_dependent_shape_representation_to_next_assembly_usage_occurrence(node)
            self.connections[next_assembly_usage_occurrence].append(node)
            self.functions[next_assembly_usage_occurrence].arg.append(f'#{node}')

    def instatiate_assembly(self, object_dict):
        assemblies_structure, valid_entities = self.get_assembly_structure()

        instantiate_ids = list(assemblies_structure.keys())
        error = True
        last_error = None
        none_primitives = set()
        assembly_shape_ids = []
        while error:
            try:
                # here we invert instantiate_ids because if the code enter inside the except
                # block, we want to loop from the last KeyError to the first. This avoids an infinite loop
                for instantiate_id in reversed(instantiate_ids):
                    if instantiate_id in object_dict or instantiate_id in none_primitives:
                        instantiate_ids.pop()
                        continue
                    product_id = self.shape_definition_representation_to_product_node(instantiate_id)
                    name = self.functions[product_id].arg[0]
                    id_shape_representation = int(self.functions[instantiate_id].arg[1][1:])
                    ids_frames = self.functions[id_shape_representation].arg[1]
                    self.parse_arguments(ids_frames)
                    assembly_frame = object_dict[ids_frames[0]]

                    assembly_shape_ids, assembly_position_ids = self.get_assembly_data(
                        assemblies_structure[instantiate_id], valid_entities, assembly_frame, object_dict)
                    assembly_positions = []
                    list_primitives = []
                    for id_shape, id_frame in zip(assembly_shape_ids, assembly_position_ids):
                        if id_shape not in none_primitives:
                            assembly_positions.append(object_dict[id_frame])
                            list_primitives.append(object_dict[id_shape])

                    if not list_primitives:
                        none_primitives.add(instantiate_id)

                    volmdlr_object = volmdlr.core.Assembly(list_primitives, assembly_positions, assembly_frame,
                                                           name=name)
                    object_dict[instantiate_id] = volmdlr_object
                    last_error = None
                error = False
            except KeyError as key:
                # Sometimes the search don't instantiate the nodes of a
                # depth in the right order, leading to error
                if last_error == key.args[0]:
                    raise NotImplementedError('Error instantiating assembly') from key
                if key.args[0] in assembly_shape_ids:
                    instantiate_ids.append(key.args[0])

                last_error = key.args[0]
        return volmdlr_object

    def to_volume_model(self, show_times: bool = False):
        """
        Translate a step file into a volmdlr object.

        :param show_times: if True, displays how many times a given class has been
            instantiated and the total time of all the instantiations of this
            given class.
        :type show_times: bool
        :return: A volmdlr solid object.
        :rtype: :class:`volmdlr.core.VolumeModel`
        """
        object_dict = {}
        times = {}
        self.create_connections()
        root_nodes = self.root_nodes
        # ------------------------------------------------------
        # TODO: This isn't a 100% right. Each SHAPE_REPRESENTATION has its own geometric context
        geometric_representation_dict = root_nodes["GEOMETRIC_REPRESENTATION_CONTEXT"]
        geometric_representation_nodes = list(geometric_representation_dict.values())
        object_dict, times = self._helper_instantiate(geometric_representation_nodes[0],
                                                      object_dict, times, show_times)
        arguments = self.functions[geometric_representation_nodes[0]].arg[:]
        self.global_uncertainty = object_dict[int(arguments[1][0][1:])]
        self.length_conversion_factor = object_dict[int(arguments[2][0][1:])]
        self.angle_conversion_factor = object_dict[int(arguments[2][1][1:])]
        # ------------------------------------------------------
        shape_representations = root_nodes["SHAPE_REPRESENTATION"]
        nodes = self.create_node_list(shape_representations)
        errors = set()
        for node in nodes:
            if node is None:
                continue
            object_dict, times = self._helper_instantiate(node, object_dict, times, show_times)

            if not object_dict[node]:
                errors.add(node)

        if show_times:
            print()
            for key, value in times.items():
                print(f'| {key} : {value}')
            print()

        if self.root_nodes["NEXT_ASSEMBLY_USAGE_OCCURRENCE"]:
            return volmdlr.core.VolumeModel([self.instatiate_assembly(object_dict)])
        primitives = []
        shapes = [object_dict[shape] for shape in shape_representations]
        for shape in shapes:
            if isinstance(shape, list):
                primitives.extend(shape)
            else:
                primitives.append(shape)
        volume_model = volmdlr.core.VolumeModel(primitives)
        # volume_model = volmdlr.core.VolumeModel([object_dict[shell_node] for shell_node in shell_nodes])
        return volume_model

    def _helper_instantiate(self, node, object_dict, times, show_times):
        """
        Helper method to translate step entities into volmdlr objects.
        """
        instantiate_ids = [node]
        error = True
        while error:
            try:
                # here we invert instantiate_ids because if the code enter inside the except
                # block, we want to loop from the last KeyError to the first. This avoids an infinite loop
                for instantiate_id in instantiate_ids[::-1]:
                    t_tracker = time.time()
                    volmdlr_object = self.instantiate(
                        self.functions[instantiate_id].name,
                        self.functions[instantiate_id].arg[:], object_dict, instantiate_id)
                    t_tracker = time.time() - t_tracker
                    object_dict[instantiate_id] = volmdlr_object
                    if show_times:
                        if volmdlr_object.__class__ not in times:
                            times[volmdlr_object.__class__] = [1, t_tracker]
                        else:
                            times[volmdlr_object.__class__][0] += 1
                            times[volmdlr_object.__class__][1] += t_tracker
                error = False
            except KeyError as key:
                # Sometimes the search don't instantiate the nodes of a
                # depth in the right order, leading to error
                instantiate_ids.append(key.args[0])

        return object_dict, times

    def to_points(self):
        """Returns a list containing all the points present in a step file."""
        object_dict = {}
        points3d = []
        for stepfunction in self.functions.values():
            if stepfunction.name == 'CARTESIAN_POINT':
                # INSTANTIATION
                name = self.functions[stepfunction.id].name
                arguments = self.functions[stepfunction.id].arg[:]
                self.parse_arguments(arguments)
                if arguments[1].count(',') == 2:
                    volmdlr_object = STEP_TO_VOLMDLR[name].from_step(
                        arguments, object_dict)
                    points3d.append(volmdlr_object)

        # remove first point because it refers to origin
        return points3d[1:]

    def plot_data(self):
        graph = self.graph().copy()

        graph.remove_nodes_from([stepfunction.id for stepfunction
                                 in self.functions.values()
                                 if stepfunction.name in ['CARTESIAN_POINT', 'DIRECTION']])
        return [plot_data.graph.NetworkxGraph(graph=graph)]


@dataclass
class StepReaderReport:
    """
    Data class to save a report after translating a step file to volmdlr object.
    """
    step_name: str = " "
    total_number_of_faces: int = 0
    faces_read: int = 0
    sucess_rate: float = 0.0
    errors: list = field(default_factory=list)
