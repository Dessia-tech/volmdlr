#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""

import time
from typing import BinaryIO, List

import matplotlib.pyplot as plt
import networkx as nx
import plot_data.graph

import volmdlr
import volmdlr.core
import volmdlr.edges
import volmdlr.faces
import volmdlr.primitives3d
import volmdlr.wires

# import webbrowser
# from jinja2 import Environment, PackageLoader, select_autoescape
# import os


def step_split_arguments(function_arg):
    """
    Split the arguments of a function that doesn't start with '(' but end with
    ')'
    ex: IN: '#123,#124,#125)'
       OUT: ['#123', '#124', '#125']
    """
    if len(function_arg) > 0 and function_arg[-1] != ')':
        function_arg += ')'
    arguments = []
    argument = ""
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


def vertex_point(arguments, object_dict):
    return object_dict[arguments[1]]


def oriented_edge(arguments, object_dict):
    return object_dict[arguments[3]]


def face_outer_bound(arguments, object_dict):
    return object_dict[arguments[1]]


def face_bound(arguments, object_dict):
    return object_dict[arguments[1]]


def surface_curve(arguments, object_dict):
    return object_dict[arguments[1]]


def seam_curve(arguments, object_dict):
    return object_dict[arguments[1]]


def trimmed_curve(arguments, object_dict):
    curve = object_dict[arguments[1]]
    point1 = object_dict[int(arguments[2][0][1:])]
    point2 = object_dict[int(arguments[3][0][1:])]
    return curve.trim(point1=point1, point2=point2)


def vertex_loop(arguments, object_dict):
    return object_dict[arguments[1]]


def pcurve(arguments, object_dict):
    # Pas besoin de mettre PCURVE ici s'il n'est pas dans STEP_TO_VOLMDLR
    return object_dict[arguments[1]]


def geometric_curve_set(arguments, object_dict):
    sub_objects = []
    for argument in arguments[1]:
        sub_obj = object_dict[int(argument[1:])]
        sub_objects.append(sub_obj)
    return sub_objects


def shell_base_surface_model(arguments, object_dict):
    # Shell3D
    return object_dict[int(arguments[1][0][1:])]


def item_defined_transformation(arguments, object_dict):
    # Frame3D
    # volmdlr_object1 = object_dict[arguments[2]]
    volmdlr_object2 = object_dict[arguments[3]]
    # TODO : how to frame map properly from these two Frame3D ?
    # return volmdlr_object2 - volmdlr_object1
    return volmdlr_object2


def manifold_surface_shape_representation(arguments, object_dict):
    # Shell3D
    shells = []
    for arg in arguments[1]:
        if isinstance(object_dict[int(arg[1:])],
                      volmdlr.faces.OpenShell3D):
            shell = object_dict[int(arg[1:])]
            shells.append(shell)
    return shells


def manifold_solid_brep(arguments, object_dict):
    return object_dict[arguments[1]]


def brep_with_voids(arguments, object_dict):
    return object_dict[arguments[1]]


def shape_representation(arguments, object_dict):
    # does it have the extra argument comming from
    # SHAPE_REPRESENTATION_RELATIONSHIP ? In this cas return
    # them
    if len(arguments) == 4:
        shells = object_dict[int(arguments[3])]
        return shells
    else:
        shells = []
        # frames = []
        for arg in arguments[1]:
            if int(arg[1:]) in object_dict and \
                    isinstance(object_dict[int(arg[1:])], list) and \
                    len(object_dict[int(arg[1:])]) == 1:
                shells.append(*object_dict[int(arg[1:])])
            elif int(arg[1:]) in object_dict and \
                    isinstance(object_dict[int(arg[1:])],
                               volmdlr.faces.OpenShell3D):
                shells.append(object_dict[int(arg[1:])])
            elif int(arg[1:]) in object_dict and \
                    isinstance(object_dict[int(arg[1:])],
                               volmdlr.Frame3D):
                # TODO: Is there something to read here ?
                pass
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
        return shells


def advanced_brep_shape_representation(arguments, object_dict):
    shells = []
    for arg in arguments[1]:
        if isinstance(object_dict[int(arg[1:])],
                      volmdlr.faces.OpenShell3D):
            shells.append(object_dict[int(arg[1:])])
    return shells


def representation_relationship_representation_relationship_with_transformation_shape_representation_relationship(
        arguments, object_dict):
    if arguments[2] in object_dict:
        if isinstance(object_dict[arguments[2]], list):
            for shell3d in object_dict[arguments[2]]:
                frame3d = object_dict[arguments[4]]
                shell3d.frame_mapping(frame3d, 'old', copy=False)
                # return shell3d
            return None
        else:
            shell3d = object_dict[arguments[2]]
            frame3d = object_dict[arguments[4]]
            shell3d.frame_mapping(frame3d, 'old', copy=False)
            # return shell3d
            return None
    else:
        return None


def bounded_curve_b_spline_curve_b_spline_curve_with_knots_curve_geometric_representation_item_rational_b_spline_curve_representation_item(
        arguments, object_dict):
    modified_arguments = [''] + arguments
    if modified_arguments[-1] == "''":
        modified_arguments.pop()
    return STEP_TO_VOLMDLR[name].from_step(
        modified_arguments, object_dict)


def bounded_surface_b_spline_surface_b_spline_surface_with_knots_geometric_representation_item_rational_b_spline_surface_representation_item_surface(
        arguments, object_dict):
    modified_arguments = [''] + arguments
    if modified_arguments[-1] == "''":
        modified_arguments.pop()
    return STEP_TO_VOLMDLR[name].from_step(
        modified_arguments, object_dict)


class StepFunction:
    def __init__(self, function_id, function_name, function_arg):
        self.id = function_id
        self.name = function_name
        self.arg = function_arg

        # TODO : Modifier ce qui suit et simplify
        if self.name == "":
            if self.arg[1][0] == 'B_SPLINE_SURFACE':
                self.simplify('B_SPLINE_SURFACE')
            if self.arg[1][0] == 'B_SPLINE_CURVE':
                self.simplify('B_SPLINE_CURVE')

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


class Step:
    def __init__(self, lines: List[str], name: str = ''):
        self.lines = lines
        self.functions, self.all_connections = self.read_lines()
        self._utd_graph = False
        self._graph = None
        self.name = name

    @property
    def graph(self):
        if not self._utd_graph:
            self._graph = self.create_graph()
            self._utd_graph = True
        return self._graph

    @classmethod
    def from_stream(cls, stream: BinaryIO = None):
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

    def read_lines(self):
        all_connections = []

        previous_line = ""
        functions = {}

        for line in self.lines:
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
            # print(function_id, function_name)
            for connec in function_arg[1:]:
                connec = connec.split(",")
                connec = connec[0].split(")")
                if connec[0][-1] != "'":
                    function_connection = int(connec[0])
                    function_connections.append(
                        (function_id, function_connection))
            # print(function_connections)

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
                    arg_list = volmdlr.core.set_to_list(argument)
                    arguments[i] = arg_list

            function = StepFunction(function_id, function_name, arguments)
            functions[function_id] = function

        return functions, all_connections

    def not_implemented(self):
        not_implemented = []
        for _, fun in self.functions.items():
            if fun.name not in STEP_TO_VOLMDLR:
                not_implemented.append(fun.name)
        return list(set(not_implemented))

    def create_graph(self):

        G = nx.Graph()
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

                self.functions[id1].arg.append('#{}'.format(id2))

            elif function.name in STEP_TO_VOLMDLR:
                G.add_node(function.id,
                           color='rgb(0, 0, 0)',
                           shape='.',
                           name=str(function.id))
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
        G.add_edges_from(self.all_connections)
        F.add_edges_from(self.all_connections)

        # Remove single nodes
        delete_nodes = []
        for node in F.nodes:
            if F.degree(node) == 0:
                delete_nodes.append(node)
        for node in delete_nodes:
            F.remove_node(node)
            G.remove_node(node)

        # if draw:
        #     # ----------------PLOT----------------
        #     pos = nx.kamada_kawai_layout(G)
        #     plt.figure()
        #     nx.draw_networkx_nodes(F, pos)
        #     nx.draw_networkx_edges(F, pos)
        #     nx.draw_networkx_labels(F, pos, labels)
        #     # ------------------------------------
        #
        # if html:
        #
        #     env = Environment(
        #         loader=PackageLoader('powertransmission', 'templates'),
        #         autoescape=select_autoescape(['html', 'xml']))
        #     template = env.get_template('graph_visJS.html')
        #
        #     nodes = []
        #     edges = []
        #     for label in list(labels.values()):
        #         nodes.append({'name': label, 'shape': 'circular'})
        #
        #     for edge in G.edges:
        #         edge_dict = {'inode1': int(edge[0]) - 1,
        #                      'inode2': int(edge[1]) - 1}
        #         edges.append(edge_dict)
        #
        #     options = {}
        #     s = template.render(
        #         name=self.stepfile,
        #         nodes=nodes,
        #         edges=edges,
        #         options=options)
        #
        #     with open('graph_visJS.html', 'wb') as file:
        #         file.write(s.encode('utf-8'))
        #
        #     webbrowser.open('file://' + os.path.realpath('graph_visJS.html'))

        return F

    def draw_graph(self, graph=None, reduced=False):
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

    def instanciate(self, name, arguments, object_dict):
        """
        """
        self.parse_arguments(arguments)

        fun_name = name.replace(', ', '_')
        fun_name = fun_name.lower()
        if hasattr(volmdlr.step, fun_name):
            volmdlr_object = getattr(volmdlr.step, fun_name)(arguments,
                                                             object_dict)

        elif name in STEP_TO_VOLMDLR and hasattr(
                STEP_TO_VOLMDLR[name], "from_step"):
            volmdlr_object = STEP_TO_VOLMDLR[name].from_step(
                arguments, object_dict)

        else:
            raise NotImplementedError(
                'Dont know how to interpret {} with args {}'.format(name,
                                                                    arguments))
        return volmdlr_object

    def to_volume_model(self, show_times=False):
        """
        no_bug_mode=True loops on instanciate method's KeyErrors until all
        the KeyErrors can be instanciated.
        show_times=True displays the numer of times a given class has been
        instanciated and the totatl time of all the instanciations of this
        given class.
        """

        object_dict = {}

        self.graph.add_node("#0")
        frame_mapping_nodes = []
        shell_nodes = []
        # sr_nodes = []
        not_shell_nodes = []
        for node in self.graph.nodes:
            if node != '#0' and self.functions[node].name == 'REPRESENTATION_RELATIONSHIP, REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION, SHAPE_REPRESENTATION_RELATIONSHIP':
                frame_mapping_nodes.append(node)
            if node != '#0' and (self.functions[node].name in ["CLOSED_SHELL", "OPEN_SHELL"]):
                shell_nodes.append(node)
            # if node != '#0' and self.functions[node].name == 'SHAPE_REPRESENTATION':
            #     # Really a shell node ?
            #     sr_nodes.append(node)
            if node != '#0' and self.functions[node].name == 'BREP_WITH_VOIDS':
                shell_nodes.append(node)
                not_shell_nodes.append(int(self.functions[node].arg[1][1:]))

        frame_mapped_shell_node = []
        for s_node in shell_nodes:
            for fm_node in frame_mapping_nodes:
                if nx.has_path(self.graph, source=fm_node, target=s_node):
                    frame_mapped_shell_node.append(s_node)
                    break
        shell_nodes_copy = shell_nodes.copy()
        [shell_nodes.remove(node) for node in frame_mapped_shell_node]

        [shell_nodes.remove(node) for node in not_shell_nodes]

        for node in shell_nodes + frame_mapping_nodes:
            self.graph.add_edge('#0', node)

        # self.draw_graph(self.graph, reduced=True, save=True)

        nodes = []
        i = 1
        new_nodes = True
        while new_nodes:
            new_nodes = list(nx.descendants_at_distance(
                self.graph, '#0', i))[::-1]
            nodes.extend(new_nodes)
            i += 1

        # nodes = dessia_common.graph.explore_tree_from_leaves(self.graph)

        times = {}
        for node in nodes[::-1]:
            # instanciate_ids = [edge[1]]
            instanciate_ids = [node]
            error = True
            while error:
                try:
                    for instanciate_id in instanciate_ids[::-1]:
                        t = time.time()
                        volmdlr_object = self.instanciate(
                            self.functions[instanciate_id].name,
                            self.functions[instanciate_id].arg[:],
                            object_dict)
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

        if show_times:
            print()
            for key, value in times.items():
                print(f'| {key} : {value}')
            print()

        shells = []
        for node in shell_nodes_copy:
            if isinstance(object_dict[node], list):
                shells.extend(object_dict[node])
            else:
                shells.append(object_dict[node])

        return volmdlr.core.VolumeModel(shells)

    def to_points(self):
        object_dict = {}
        points3d = []
        for stepfunction in self.functions.values():
            if stepfunction.name == 'CARTESIAN_POINT':
                # INSTANCIATION
                name = self.functions[stepfunction.id].name
                arguments = self.functions[stepfunction.id].arg[:]
                self.parse_arguments(arguments)
                # for i, arg in enumerate(arguments):
                #     if type(arg) == str and arg[0] == '#':
                #         arguments[i] = int(arg[1:])
                # print(arguments)
                if arguments[1].count(',') == 2:
                    volmdlr_object = STEP_TO_VOLMDLR[name].from_step(
                        arguments, object_dict)
                    points3d.append(volmdlr_object)

        # remove first point because it refers to origin
        return points3d[1:]

    def plot_data(self):
        graph = self.graph.copy()

        graph.remove_nodes_from([stepfunction.id for stepfunction
                                 in self.functions.values()
                                 if stepfunction.name in ['CARTESIAN_POINT', 'DIRECTION']])
        return [plot_data.graph.NetworkxGraph(graph=graph)]


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
    'SURFACE_OF_LINEAR_EXTRUSION': volmdlr.primitives3d.BSplineExtrusion,
    # CAN BE A BSplineSurface3D
    'SURFACE_OF_REVOLUTION': None,
    'UNIFORM_SURFACE': volmdlr.faces.BSplineSurface3D,
    'QUASI_UNIFORM_SURFACE': volmdlr.faces.BSplineSurface3D,
    'RECTANGULAR_COMPOSITE_SURFACE': volmdlr.faces.PlaneFace3D,  # TOPOLOGICAL FACES
    'CURVE_BOUNDED_SURFACE': volmdlr.faces.PlaneFace3D,  # TOPOLOGICAL FACE

    # added on 12/08/2021 by Mack in order to read BsplinePipe
    'BOUNDED_SURFACE, B_SPLINE_SURFACE, B_SPLINE_SURFACE_WITH_KNOTS, GEOMETRIC_REPRESENTATION_ITEM, RATIONAL_B_SPLINE_SURFACE, REPRESENTATION_ITEM, SURFACE': volmdlr.faces.BSplineSurface3D,

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
