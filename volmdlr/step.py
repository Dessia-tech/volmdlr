#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ISO STEP reader/writer.
"""

import time
from typing import List
from dataclasses import dataclass, field
import numpy as npy

import matplotlib.pyplot as plt
import networkx as nx
import plot_data.graph

import dessia_common.core as dc  # isort: skip
from dessia_common.files import BinaryFile  # isort: skip

import volmdlr
import volmdlr.core
import volmdlr.edges
import volmdlr.faces
import volmdlr.primitives3d
import volmdlr.wires


def set_to_list(step_set):
    """
    Convert a string representation of a set to a list of strings.

    :param step_set: String representation of a set, e.g. "{A,B,C}"
    :type step_set: str
    :return: List of strings, e.g. ["A", "B", "C"]
    :rtype: List[str]
    """
    char_list = step_set.split(',')
    char_list[0] = char_list[0][1:]
    char_list[-1] = char_list[-1][:-1]
    return list(char_list)


def step_split_arguments(function_arg):
    """
    Split the arguments of a function that doesn't start with '(' but end with ')'.

    ex: IN: '#123,#124,#125)'
       OUT: ['#123', '#124', '#125']
    """
    if len(function_arg) > 0 and function_arg[-1] != ')':
        function_arg += ')'
    arguments = []
    argument = ""
    if len(function_arg) > 0 and function_arg[0] == "(":
        function_arg += ")"
    parenthesis = 1
    for char in function_arg:
        if char == "(":
            parenthesis += 1

        if char != "," or parenthesis > 1:
            argument += char
        else:
            arguments.append(argument)
            argument = ""

        if char == ")":
            parenthesis -= 1
            if parenthesis == 0:
                arguments.append(argument[:-1])
                argument = ""
                break
    return arguments


def uncertainty_measure_with_unit(arguments, object_dict):
    """
    Gets the global length uncertainty.

    :param arguments: step primitive arguments
    :param object_dict: dictionary containing already instantiated objects.
    :return: Global length uncertainty.
    """
    length_measure = float(arguments[0].split('(')[1][:-1])
    return length_measure * object_dict[arguments[1]]


def conversion_based_unit_length_unit_named_unit(arguments, object_dict):
    """
    Gets the conversion based unit length.

    :param arguments: step primitive arguments
    :param object_dict: dictionary containing already instantiated objects.
    :return: conversion based unit length.
    """
    return object_dict[arguments[1]]


def length_measure_with_unit(arguments, object_dict):
    """
    Calculates the step file's si unit conversion factor.

    :param arguments: step primitive arguments
    :param object_dict: dictionary containing already instantiated objects.
    :return: si unit conversion factor.
    """
    if "(" in arguments[0]:
        length_measure = float(arguments[0].split('(')[1][:-1])
    else:
        length_measure = float(arguments[0])
    length_si_unit = object_dict[arguments[1]]
    return length_measure * length_si_unit


def conversion_based_unit_named_unit_plane_angle_unit(arguments, object_dict):
    """
    Gets the conversion based plane unit angle.

    :param arguments: step primitive arguments
    :param object_dict: dictionary containing already instantiated objects.
    :return: conversion based unit length.
    """
    return object_dict[arguments[1]]


def named_unit_plane_angle_unit_si_unit(arguments, object_dict):
    """
    Returns the dimension of plane angle measure.

    :param arguments: step primitive arguments
    :param object_dict: dictionary containing already instantiated objects.
    :return: SI unit dimension.
    """
    return SI_PREFIX[arguments[1]]


def named_unit_si_unit_solid_angle_unit(arguments, object_dict):
    """
    Returns the dimension of solid angle measure.

    :param arguments: step primitive arguments
    :param object_dict: dictionary containing already instantiated objects.
    :return: SI unit dimension.
    """
    return SI_PREFIX[arguments[1]]


def plane_angle_measure_with_unit(arguments, object_dict):
    """
    Returns the angle plane measure with the right unit.

    :param arguments: step primitive arguments
    :param object_dict: dictionary containing already instantiated objects.
    :return: angle measure in SI unit.
    """
    angle_measure = float(arguments[0].split('(')[1][:-1])
    angle_si_unit = object_dict[arguments[1]]
    return angle_measure * angle_si_unit


def length_unit_named_unit_si_unit(arguments, object_dict):
    """
    Gets the length si unit.

    :param arguments: step primitive arguments
    :param object_dict: dictionary containing already instantiated objects.
    :return: length si unit
    """
    si_unit_length = SI_PREFIX[arguments[1]]
    return si_unit_length


def geometric_representation_context_global_uncertainty_assigned_context_global_unit_assigned_context_representation_context(
        arguments, object_dict):
    """
    Gets the global length uncertainty.

    :param arguments: step primitive arguments
    :param object_dict: dictionary containing already instantiated objects.
    :return: Global length uncertainty.
    """
    length_global_uncertainty = object_dict[int(arguments[1][0][1:])]
    length_conversion_factor = object_dict[int(arguments[2][0][1:])]
    angle_conversion_factor = object_dict[int(arguments[2][1][1:])]
    return length_global_uncertainty, length_conversion_factor, angle_conversion_factor


def vertex_point(arguments, object_dict):
    """
    Returns the data in case of a VERTEX.
    """
    return object_dict[arguments[1]]


def axis1_placement(arguments, object_dict):
    """
    Returns the data in case of a AXIS1_PLACEMENT.
    """
    return object_dict[arguments[1]], object_dict[arguments[2]]


def oriented_edge(arguments, object_dict):
    """
    Returns the data in case of an ORIENTED_EDGE.
    """
    if not object_dict[arguments[3]]:
        # This can happen when the is too small
        return None
    edge_orientation = arguments[4]
    if edge_orientation == '.T.':
        return object_dict[arguments[3]]
    return object_dict[arguments[3]].reverse()


def face_outer_bound(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    return object_dict[arguments[1]]


def face_bound(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    return object_dict[arguments[1]]


def surface_curve(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    return object_dict[arguments[1]]


def seam_curve(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE
    """
    return object_dict[arguments[1]]


def trimmed_curve(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE
    """

    curve = object_dict[arguments[1]]
    point1 = object_dict[int(arguments[2][0][1:])]
    point2 = object_dict[int(arguments[3][0][1:])]
    return curve.trim(point1=point1, point2=point2)


def vertex_loop(arguments, object_dict):
    """
    Returns the data in case of a VERTEX_LOOP.
    """
    return object_dict[arguments[1]]


def pcurve(arguments, object_dict):
    """
    Returns the data in case of a PCURVE.
    """
    return object_dict[arguments[1]]


def geometric_curve_set(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    sub_objects = []
    for argument in arguments[1]:
        sub_obj = object_dict[int(argument[1:])]
        sub_objects.append(sub_obj)
    return sub_objects


def shell_based_surface_model(arguments, object_dict):
    """
    Returns the data in case of a Shell3D.
    """
    return object_dict[int(arguments[1][0][1:])]


def item_defined_transformation(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    # Frame3D
    volmdlr_object1 = object_dict[arguments[2]]
    volmdlr_object2 = object_dict[arguments[3]]
    # TODO : how to frame map properly from these two Frame3D ?
    # return volmdlr_object2 - volmdlr_object1
    return [volmdlr_object1, volmdlr_object2]


def manifold_surface_shape_representation(arguments, object_dict):
    """
    Returns the data in case of a manifold_surface_shape_representation, interpreted as shell3D.
    """
    shells = []
    for arg in arguments[1]:
        if isinstance(object_dict[int(arg[1:])],
                      volmdlr.faces.OpenShell3D):
            shell = object_dict[int(arg[1:])]
            shells.append(shell)
    return shells


def faceted_brep(arguments, object_dict):
    """
    Returns the data in case of a faceted_brep entity, interpreted as shell3D.
    """
    return object_dict[arguments[1]]


def faceted_brep_shape_representation(arguments, object_dict):
    """
    Returns the data in case of a faceted_brep_shape_representation, interpreted as shell3D.
    """
    shells = []
    for arg in arguments[1]:
        if isinstance(object_dict[int(arg[1:])],
                      volmdlr.faces.OpenShell3D):
            shell = object_dict[int(arg[1:])]
            shells.append(shell)
    return shells


def manifold_solid_brep(arguments, object_dict):
    """
    Returns the data in case of a manifold_solid_brep with voids.
    """
    return object_dict[arguments[1]]


def brep_with_voids(arguments, object_dict):
    """
    Returns the data in case of a BREP with voids.
    """
    return object_dict[arguments[1]]


def shape_representation(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    # does it have the extra argument comming from
    # SHAPE_REPRESENTATION_RELATIONSHIP ? In this cas return
    # them
    if len(arguments) == 4:
        shells = object_dict[int(arguments[3])]
        return shells
    shells = []
    frames = []
    for arg in arguments[1]:
        if int(arg[1:]) in object_dict and \
                isinstance(object_dict[int(arg[1:])], list) and \
                len(object_dict[int(arg[1:])]) == 1:
            shells.append(*object_dict[int(arg[1:])])
        elif int(arg[1:]) in object_dict and \
                isinstance(object_dict[int(arg[1:])],
                           volmdlr.faces.OpenShell3D):
            shells.append(object_dict[int(arg[1:])])
        elif int(arg[1:]) in object_dict and isinstance(object_dict[int(arg[1:])], volmdlr.Frame3D):
            # TODO: Is there something to read here ?
            frame = object_dict[int(arg[1:])]
            if not all(component is None for component in [frame.u, frame.u, frame.w]):
                frames.append(frame)
        elif int(arg[1:]) in object_dict and \
                isinstance(object_dict[int(arg[1:])],
                           volmdlr.edges.Arc3D):
            shells.append(object_dict[int(arg[1:])])
        elif int(arg[1:]) in object_dict and \
                isinstance(object_dict[int(arg[1:])],
                           volmdlr.edges.BSplineCurve3D):
            shells.append(object_dict[int(arg[1:])])
        else:
            pass
    if not shells and frames:
        return frames
    return shells


def advanced_brep_shape_representation(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    shells = []
    for arg in arguments[1]:
        if isinstance(object_dict[int(arg[1:])],
                      volmdlr.faces.OpenShell3D):
            shells.append(object_dict[int(arg[1:])])
    return shells


def frame_map_closed_shell(closed_shells, item_defined_transformation_frames, shape_representation_frames):
    """
    Frame maps a closed shell in an assembly to its good position.

    :param closed_shells: DESCRIPTION
    :type closed_shells: volmdlr.faces.OpenShell3D
    :param item_defined_transformation_frames: DESCRIPTION
    :type item_defined_transformation_frames: TYPE
    :param shape_representation_frames: DESCRIPTION
    :type shape_representation_frames: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    if item_defined_transformation_frames[0] == item_defined_transformation_frames[1]:
        return closed_shells
    if shape_representation_frames[0].origin.is_close(volmdlr.O3D):
        global_frame = shape_representation_frames[0]
    else:
        global_frame = [frame for frame in item_defined_transformation_frames if frame.origin.is_close(volmdlr.O3D)][0]
    transformed_frame = [frame for frame in item_defined_transformation_frames if frame != global_frame][0]
    new_closedshells = []

    for shell3d in closed_shells:
        basis_a = global_frame.basis()
        basis_b = transformed_frame.basis()
        A = npy.array([[basis_a.vectors[0].x, basis_a.vectors[0].y, basis_a.vectors[0].z],
                       [basis_a.vectors[1].x, basis_a.vectors[1].y, basis_a.vectors[1].z],
                       [basis_a.vectors[2].x, basis_a.vectors[2].y, basis_a.vectors[2].z]])
        B = npy.array([[basis_b.vectors[0].x, basis_b.vectors[0].y, basis_b.vectors[0].z],
                       [basis_b.vectors[1].x, basis_b.vectors[1].y, basis_b.vectors[1].z],
                       [basis_b.vectors[2].x, basis_b.vectors[2].y, basis_b.vectors[2].z]])
        transfer_matrix = npy.linalg.solve(A, B)
        u_vector = volmdlr.Vector3D(*transfer_matrix[0])
        v_vector = volmdlr.Vector3D(*transfer_matrix[1])
        w_vector = volmdlr.Vector3D(*transfer_matrix[2])
        new_frame = volmdlr.Frame3D(transformed_frame.origin, u_vector,
                                    v_vector,
                                    w_vector)
        new_faces = [face.frame_mapping(new_frame, 'old') for face in shell3d.faces]
        new_closed_shell3d = volmdlr.faces.ClosedShell3D(new_faces)
        new_closedshells.append(new_closed_shell3d)
    return new_closedshells


def representation_relationship_representation_relationship_with_transformation_shape_representation_relationship(
        arguments, object_dict):
    """
    Representation relationship with transformation shape. To clarify.
    """
    if arguments[2] in object_dict:
        if isinstance(object_dict[arguments[2]], list):  # arguments = {, , [], [], item_....}
            if object_dict[arguments[2]] and not isinstance(object_dict[arguments[2]][0], volmdlr.Frame3D) \
                    and isinstance(object_dict[arguments[3]][0], volmdlr.Frame3D):
                return frame_map_closed_shell(object_dict[arguments[2]],
                                              object_dict[arguments[4]], object_dict[arguments[3]])

            if object_dict[arguments[2]] and isinstance(object_dict[arguments[2]][0], volmdlr.Frame3D) and \
                    not isinstance(object_dict[arguments[3]][0], volmdlr.Frame3D):
                return frame_map_closed_shell(object_dict[arguments[3]],
                                              object_dict[arguments[4]], object_dict[arguments[2]])
            return []
        return []
    return []


def bounded_curve_b_spline_curve_b_spline_curve_with_knots_curve_geometric_representation_item_rational_b_spline_curve_representation_item(
        arguments, object_dict):
    """
    Bounded b spline with knots curve geometric representation item. To clarify.
    """
    modified_arguments = [''] + arguments
    if modified_arguments[-1] == "''":
        modified_arguments.pop()
    return STEP_TO_VOLMDLR['BOUNDED_CURVE, '
                           'B_SPLINE_CURVE, '
                           'B_SPLINE_CURVE_WITH_KNOTS, '
                           'CURVE, GEOMETRIC_REPRESENTATION_ITEM, '
                           'RATIONAL_B_SPLINE_CURVE, '
                           'REPRESENTATION_ITEM'].from_step(
        modified_arguments, object_dict)


def bounded_surface_b_spline_surface_b_spline_surface_with_knots_geometric_representation_item_rational_b_spline_surface_representation_item_surface(
        arguments, object_dict):
    """
    Bounded b spline surface with knots curve geometric representation item. To clarify.
    """
    modified_arguments = [''] + arguments
    if modified_arguments[-1] == "''":
        modified_arguments.pop()
    return STEP_TO_VOLMDLR['BOUNDED_SURFACE, B_SPLINE_SURFACE, '
                           'B_SPLINE_SURFACE_WITH_KNOTS, '
                           'GEOMETRIC_REPRESENTATION_ITEM, '
                           'RATIONAL_B_SPLINE_SURFACE, '
                           'REPRESENTATION_ITEM, SURFACE'].from_step(
        modified_arguments, object_dict)


def bounded_surface_b_spline_surface_b_spline_surface_with_knots_surface_geometric_representation_item_rational_b_spline_surface_representation_item(
        arguments, object_dict):
    """
    Bounded b spline surface with knots curve geometric representation item. To clarify.
    """
    modified_arguments = [''] + arguments
    if modified_arguments[-1] == "''":
        modified_arguments.pop()
    return STEP_TO_VOLMDLR['BOUNDED_SURFACE, B_SPLINE_SURFACE, '
                           'B_SPLINE_SURFACE_WITH_KNOTS, '
                           'GEOMETRIC_REPRESENTATION_ITEM, '
                           'RATIONAL_B_SPLINE_SURFACE, '
                           'REPRESENTATION_ITEM, SURFACE'].from_step(
        modified_arguments, object_dict)


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
        self.functions, self.all_connections, self.connections = self.read_lines(lines)
        self._graph = None
        self.global_uncertainty = 1e-6
        self.length_conversion_factor = 1
        self.angle_conversion_factor = 1
        self.read_diagnostic = StepReaderReport

        dc.DessiaObject.__init__(self, name=name)

    def graph(self):
        if not self._graph:
            self._graph = self.create_graph()
        return self._graph

    @classmethod
    def from_stream(cls, stream: BinaryFile):
        stream.seek(0)
        lines = []
        for line in stream:
            line = line.decode("ISO-8859-1")
            line = line.replace("\r", "")
            lines.append(line)
        return cls(lines)

    @classmethod
    def from_file(cls, filepath: str = None):
        with open(filepath, "r", encoding="ISO-8859-1") as file:
            lines = []
            for line in file:
                lines.append(line)
        return cls(lines)

    def read_lines(self, lines):
        all_connections = []
        dict_connections = {}
        previous_line = ""
        functions = {}

        for line in lines:
            line = line.replace(" ", "")
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

            function = line.split("=")
            function_id = int(function[0][1:])
            function_name_arg = function[1].split("(", 1)
            function_name = function_name_arg[0]
            function_arg = function_name_arg[1].split("#")
            function_connections = []
            connections = []
            # print(function_id, function_name)
            for connec in function_arg[1:]:
                connec = connec.split(",")
                connec = connec[0].split(")")
                if connec[0][-1] != "'":
                    function_connection = int(connec[0])
                    connections.append(function_connection)
                    function_connections.append(
                        (function_id, function_connection))
            # print(function_connections)
            dict_connections[function_id] = connections
            all_connections.extend(function_connections)

            previous_line = str()

            # FUNCTION ARGUMENTS
            function_arg = function_name_arg[1]
            arguments = step_split_arguments(function_arg)
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
                        function_connections.append(
                            (function_id, int(arg[1:])))
            # print('=', function_connections)

            for i, argument in enumerate(arguments):
                if argument[:2] == '(#' and argument[-1] == ')':
                    arg_list = set_to_list(argument)
                    arguments[i] = arg_list

            function = StepFunction(function_id, function_name, arguments)
            functions[function_id] = function

        return functions, all_connections, dict_connections

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
        :rtype: nx.DiGraph
        """
        F = nx.DiGraph()
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
                F.add_node(function.id,
                           color='rgb(0, 0, 0)',
                           shape='.',
                           name=str(function.id))
                labels[function.id] = str(function.id) + ' ' + function.name

        # Delete connection if node not found
        node_list = list(F.nodes())
        delete_connection = []
        for connection in self.all_connections:
            if connection[0] not in node_list \
                    or connection[1] not in node_list:
                delete_connection.append(connection)
        for delete in delete_connection:
            self.all_connections.remove(delete)

        # Create graph connections
        F.add_edges_from(self.all_connections)

        # Remove single nodes
        delete_nodes = []
        for node in F.nodes:
            if F.degree(node) == 0:
                delete_nodes.append(node)
        for node in delete_nodes:
            F.remove_node(node)
            # G.remove_node(node)
        return F

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

    def step_subfunctions(self, subfunctions):
        subfunctions = subfunctions[0]
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
            (subfunction_names[i], step_split_arguments(subfunction_args[i]))
            for i in range(len(subfunction_names))]

    def parse_arguments(self, arguments):
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
            if hasattr(volmdlr.step, fun_name):
                volmdlr_object = getattr(volmdlr.step, fun_name)(arguments, object_dict)

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

    def create_node_list(self, stack):
        """
        Step functions graph as a list of nodes.

        :param stack: Initial list of shell nodes and assemblies entities.
        :type stack: List[int]
        :return: A list of nodes in the right order of dependency.
        :rtype: List[int]
        """
        list_head = []
        list_nodes = []
        visited_set = set()
        while stack:
            node = stack.pop(0)
            name = self.functions[node].name
            if node not in visited_set and name in STEP_TO_VOLMDLR:
                visited_set.add(node)
                if self.connections[node]:
                    list_nodes.append(node)
                    for c in self.connections[node]:
                        if c not in visited_set:
                            stack.append(c)
                else:
                    # Entities without connections should be instatiate first
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
        if len(self.functions[id_shape_representation].arg) != 4:
            id_shape_representation = int(arguments[2][1:])
        if self.functions[id_shape_representation].name == "SHAPE_REPRESENTATION":
            if len(self.functions[id_shape_representation].arg) != 4:
                return None
            id_representation_entity = int(self.functions[id_shape_representation].arg[3][1:])
            id_solid_entity = int(self.functions[id_representation_entity].arg[1][0][1:])
            id_shell = self.functions[id_solid_entity].arg[1]
            if isinstance(id_shell, list):
                return int(id_shell[0][1:])
            return int(id_shell[1:])
        id_representation_entity = int(self.functions[id_shape_representation].arg[1][1][1:])
        id_shell = self.functions[id_representation_entity].arg[1]
        if isinstance(id_shell, list):
            return int(id_shell[0][1:])
        return int(id_shell[1:])

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
        frame_mapping_nodes = []
        shell_nodes = []
        for function in self.functions.values():
            if function.name == 'SHAPE_REPRESENTATION_RELATIONSHIP':
                # Create short cut from id1 to id2
                id1 = int(function.arg[2][1:])
                id2 = int(function.arg[3][1:])
                self.connections[id1].append(id2)
                self.functions[id1].arg.append(f'#{id2}')
        not_shell_nodes = []
        frame_mapped_shell_node = []

        for node in list(self.functions.keys()):
            if self.functions[node].name == 'REPRESENTATION_RELATIONSHIP, ' \
                                            'REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION, ' \
                                            'SHAPE_REPRESENTATION_RELATIONSHIP':
                frame_mapping_nodes.append(node)
                shell_node = self.get_frame_mapped_shell_node(node)
                if shell_node:
                    frame_mapped_shell_node.append(shell_node)
            elif self.functions[node].name in {"CLOSED_SHELL", "OPEN_SHELL"}:
                shell_nodes.append(node)
            elif self.functions[node].name == 'GEOMETRIC_REPRESENTATION_CONTEXT, ' \
                                              'GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT, ' \
                                              'GLOBAL_UNIT_ASSIGNED_CONTEXT, REPRESENTATION_CONTEXT':
                object_dict, times = self._helper_instantiate(node, object_dict, times, show_times)
                arguments = self.functions[node].arg[:]
                self.global_uncertainty = object_dict[int(arguments[1][0][1:])]
                self.length_conversion_factor = object_dict[int(arguments[2][0][1:])]
                self.angle_conversion_factor = object_dict[int(arguments[2][1][1:])]

            elif self.functions[node].name == 'BREP_WITH_VOIDS':
                shell_nodes.append(node)
                not_shell_nodes.append(int(self.functions[node].arg[1][1:]))

        shell_nodes_copy = shell_nodes.copy()
        remove_nodes = list(set(frame_mapped_shell_node + not_shell_nodes))
        for node in remove_nodes:
            shell_nodes.remove(node)

        nodes = self.create_node_list(shell_nodes + frame_mapping_nodes)
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
        shells = []
        step_number_faces = 0
        faces_read = 0
        if frame_mapping_nodes:
            for node in frame_mapping_nodes:
                shells.extend(object_dict[node])
        if not shells:
            for node in shell_nodes_copy:
                volmdlr_object = object_dict[node]
                if isinstance(volmdlr_object, list):
                    shells.extend(object_dict[node])
                else:
                    shells.append(volmdlr_object)
                    faces_read += len(volmdlr_object.faces)
                    step_number_faces += len(self.functions[node].arg[1])
            if step_number_faces and faces_read:
                self.read_diagnostic = StepReaderReport(self.name, step_number_faces, faces_read,
                                                        faces_read / step_number_faces, list(errors))
        volume_model = volmdlr.core.VolumeModel(shells)
        return volume_model

    def _helper_instantiate(self, node, object_dict, times, show_times):
        """
        Helper method to translate step entities into volmdlr objects.
        """
        instanciate_ids = [node]
        error = True
        while error:
            try:
                for instanciate_id in instanciate_ids[::-1]:
                    t = time.time()

                    volmdlr_object = self.instantiate(
                        self.functions[instanciate_id].name,
                        self.functions[instanciate_id].arg[:], object_dict, instanciate_id)
                    t = time.time() - t
                    object_dict[instanciate_id] = volmdlr_object
                    if show_times:
                        if volmdlr_object.__class__ not in times:
                            times[volmdlr_object.__class__] = [1, t]
                        else:
                            times[volmdlr_object.__class__][0] += 1
                            times[volmdlr_object.__class__][1] += t
                error = False
            except KeyError as key:
                # Sometimes the bfs search don't instanciate the nodes of a
                # depth in the right order, leading to error
                instanciate_ids.append(key.args[0])

        return object_dict, times

    def to_points(self):
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


STEP_TO_VOLMDLR = {
    # GEOMETRICAL ENTITIES
    'CARTESIAN_POINT': volmdlr.Point3D,
    'DIRECTION': volmdlr.Vector3D,
    'VECTOR': volmdlr.Vector3D,

    'AXIS1_PLACEMENT': None,
    'AXIS2_PLACEMENT_2D': None,  # ??????????????????
    'AXIS2_PLACEMENT_3D': volmdlr.Frame3D,

    'LINE': volmdlr.edges.Line3D,  # LineSegment3D,
    'CIRCLE': volmdlr.wires.Circle3D,
    'ELLIPSE': volmdlr.wires.Ellipse3D,
    'PARABOLA': None,
    'HYPERBOLA': None,
    # 'PCURVE': None,
    'CURVE_REPLICA': None,
    'OFFSET_CURVE_3D': None,
    'TRIMMED_CURVE': None,  # BSplineCurve3D cannot be trimmed on FreeCAD
    'B_SPLINE_CURVE': volmdlr.edges.BSplineCurve3D,
    'B_SPLINE_CURVE_WITH_KNOTS': volmdlr.edges.BSplineCurve3D,
    'BEZIER_CURVE': volmdlr.edges.BSplineCurve3D,
    'RATIONAL_B_SPLINE_CURVE': volmdlr.edges.BSplineCurve3D,
    'UNIFORM_CURVE': volmdlr.edges.BSplineCurve3D,
    'QUASI_UNIFORM_CURVE': volmdlr.edges.BSplineCurve3D,
    'SURFACE_CURVE': None,  # TOPOLOGICAL EDGE
    'SEAM_CURVE': None,
    # LineSegment3D, # TOPOLOGICAL EDGE ############################
    'COMPOSITE_CURVE_SEGMENT': None,  # TOPOLOGICAL EDGE
    'COMPOSITE_CURVE': volmdlr.wires.Wire3D,  # TOPOLOGICAL WIRE
    'COMPOSITE_CURVE_ON_SURFACE': volmdlr.wires.Wire3D,  # TOPOLOGICAL WIRE
    'BOUNDARY_CURVE': volmdlr.wires.Wire3D,  # TOPOLOGICAL WIRE

    'PLANE': volmdlr.faces.Plane3D,
    'CYLINDRICAL_SURFACE': volmdlr.faces.CylindricalSurface3D,
    'CONICAL_SURFACE': volmdlr.faces.ConicalSurface3D,
    'SPHERICAL_SURFACE': volmdlr.faces.SphericalSurface3D,
    'TOROIDAL_SURFACE': volmdlr.faces.ToroidalSurface3D,
    'DEGENERATE_TOROIDAL_SURFACE': None,
    'B_SPLINE_SURFACE_WITH_KNOTS': volmdlr.faces.BSplineSurface3D,
    'B_SPLINE_SURFACE': volmdlr.faces.BSplineSurface3D,
    'BEZIER_SURFACE': volmdlr.faces.BSplineSurface3D,
    'OFFSET_SURFACE': None,
    'SURFACE_REPLICA': None,
    'RATIONAL_B_SPLINE_SURFACE': volmdlr.faces.BSplineSurface3D,
    'RECTANGULAR_TRIMMED_SURFACE': None,
    'SURFACE_OF_LINEAR_EXTRUSION': volmdlr.faces.ExtrusionSurface3D,
    # CAN BE A BSplineSurface3D
    'SURFACE_OF_REVOLUTION': volmdlr.faces.RevolutionSurface3D,
    'UNIFORM_SURFACE': volmdlr.faces.BSplineSurface3D,
    'QUASI_UNIFORM_SURFACE': volmdlr.faces.BSplineSurface3D,
    'RECTANGULAR_COMPOSITE_SURFACE': volmdlr.faces.PlaneFace3D,  # TOPOLOGICAL FACES
    'CURVE_BOUNDED_SURFACE': volmdlr.faces.PlaneFace3D,  # TOPOLOGICAL FACE

    # Bsplines
    'BOUNDED_SURFACE, B_SPLINE_SURFACE, B_SPLINE_SURFACE_WITH_KNOTS, GEOMETRIC_REPRESENTATION_ITEM,'
    ' RATIONAL_B_SPLINE_SURFACE, REPRESENTATION_ITEM, SURFACE': volmdlr.faces.BSplineSurface3D,
    "BOUNDED_SURFACE, B_SPLINE_SURFACE, B_SPLINE_SURFACE_WITH_KNOTS, SURFACE, GEOMETRIC_REPRESENTATION_ITEM,"
    " RATIONAL_B_SPLINE_SURFACE, REPRESENTATION_ITEM": volmdlr.faces.BSplineSurface3D,
    # TOPOLOGICAL ENTITIES
    'VERTEX_POINT': None,

    'EDGE_CURVE': volmdlr.edges.Edge,  # LineSegment3D, # TOPOLOGICAL EDGE
    'ORIENTED_EDGE': None,  # TOPOLOGICAL EDGE
    # The one above can influence the direction with their last argument
    # TODO : maybe take them into consideration

    'FACE_BOUND': None,  # TOPOLOGICAL WIRE
    'FACE_OUTER_BOUND': None,  # TOPOLOGICAL WIRE
    # Both above can influence the direction with their last argument
    # TODO : maybe take them into consideration
    'EDGE_LOOP': volmdlr.wires.Contour3D,  # TOPOLOGICAL WIRE
    'POLY_LOOP': volmdlr.wires.Contour3D,  # TOPOLOGICAL WIRE
    'VERTEX_LOOP': None,  # TOPOLOGICAL WIRE

    'ADVANCED_FACE': volmdlr.faces.Face3D,
    'FACE_SURFACE': volmdlr.faces.Face3D,

    'CLOSED_SHELL': volmdlr.faces.ClosedShell3D,
    'OPEN_SHELL': volmdlr.faces.OpenShell3D,
    #        'ORIENTED_CLOSED_SHELL': None,
    'CONNECTED_FACE_SET': volmdlr.faces.OpenShell3D,
    'GEOMETRIC_CURVE_SET': None,

    # step subfunctions
    'UNCERTAINTY_MEASURE_WITH_UNIT': None,
    'CONVERSION_BASED_UNIT, LENGTH_UNIT, NAMED_UNIT': None,
    'LENGTH_MEASURE_WITH_UNIT': None,
    'LENGTH_UNIT, NAMED_UNIT, SI_UNIT': None,
    'PLANE_ANGLE_MEASURE_WITH_UNIT': None,
    'NAMED_UNIT, PLANE_ANGLE_UNIT, SI_UNIT': None,
    'CONVERSION_BASED_UNIT, NAMED_UNIT, PLANE_ANGLE_UNIT': None,
    'GEOMETRIC_REPRESENTATION_CONTEXT, GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT, GLOBAL_UNIT_ASSIGNED_CONTEXT, REPRESENTATION_CONTEXT': None,
    'REPRESENTATION_RELATIONSHIP, REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION, SHAPE_REPRESENTATION_RELATIONSHIP': volmdlr.faces.OpenShell3D.translation,
    'SHELL_BASED_SURFACE_MODEL': None,
    'MANIFOLD_SURFACE_SHAPE_REPRESENTATION': None,
    'MANIFOLD_SOLID_BREP': None,
    'BREP_WITH_VOIDS': None,
    'SHAPE_REPRESENTATION': None,
    'ADVANCED_BREP_SHAPE_REPRESENTATION': None,
    'ITEM_DEFINED_TRANSFORMATION': None,
    'SHAPE_REPRESENTATION_RELATIONSHIP': None,

    'BOUNDED_CURVE, B_SPLINE_CURVE, B_SPLINE_CURVE_WITH_KNOTS, CURVE, GEOMETRIC_REPRESENTATION_ITEM, RATIONAL_B_SPLINE_CURVE, REPRESENTATION_ITEM': volmdlr.edges.BSplineCurve3D
}

VOLMDLR_TO_STEP = {}
for k, v in STEP_TO_VOLMDLR.items():
    if v:
        if v in VOLMDLR_TO_STEP:
            VOLMDLR_TO_STEP[v].append(k)
        else:
            VOLMDLR_TO_STEP[v] = [k]

SI_PREFIX = {'.EXA.': 1e18, '.PETA.': 1e15, '.TERA.': 1e12, '.GIGA.': 1e9, '.MEGA.': 1e6, '.KILO.': 1e3,
             '.HECTO.': 1e2, '.DECA.': 1e1, '$': 1, '.DECI.': 1e-1, '.CENTI.': 1e-2, '.MILLI.': 1e-3, '.MICRO.': 1e-6,
             '.NANO.': 1e-9, '.PICO.': 1e-12, '.FEMTO.': 1e-15, '.ATTO.': 1e-18}

STEP_REPRESENTATION_ENTITIES = {"ADVANCED_BREP_SHAPE_REPRESENTATION", "FACETED_BREP_SHAPE_REPRESENTATION",
                                "MANIFOLD_SURFACE_SHAPE_REPRESENTATION",
                                "GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION",
                                "GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION",
                                "EDGE_BASED_WIREFRAME_SHAPE_REPRESENTATION"
                                }
