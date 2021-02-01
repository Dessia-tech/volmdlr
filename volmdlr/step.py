#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

"""


import matplotlib.pyplot as plt
import networkx as nx
import volmdlr
import volmdlr.core
import volmdlr.primitives3d
import volmdlr.edges
import volmdlr.wires
import volmdlr.faces
import plot_data.graph

import webbrowser


def step_split_arguments(function_arg):
    """
    Split the arguments of a function that doesn't start with '(' but end with ')'
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
    def __init__(self, stepfile):
        self.stepfile = stepfile

        self.functions, self.all_connections = self.read_functions()

        self.upd_graph = False

    def read_functions(self):
        f = open(self.stepfile, "r", encoding="ISO-8859-1")

        all_connections = []

        previous_line = ""
        functions = {}

        for line in f:
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
            for connec in function_arg[1:]:
                connec = connec.split(",")
                connec = connec[0].split(")")
                if connec[0][-1] != "'":
                    function_connection = int(connec[0])
                    function_connections.append(
                        (function_id, function_connection))

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

            for i, argument in enumerate(arguments):
                if argument[:2] == '(#' and argument[-1] == ')':
                    arg_list = volmdlr.core.set_to_list(argument)
                    arguments[i] = arg_list

            if function_id == 918:
                print(function_name, arguments)
            function = StepFunction(function_id, function_name, arguments)
            functions[function_id] = function

        f.close()

        return functions, all_connections

    def create_graph(self, draw=False, html=False):

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
        #         edge_dict = {}
        #         edge_dict['inode1'] = int(edge[0]) - 1
        #         edge_dict['inode2'] = int(edge[1]) - 1
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

    def draw_graph(self):
        graph = self.create_graph()
        labels = {}
        for id_nb, function in self.functions.items():
            if id_nb in graph.nodes:
                labels[id_nb] = str(id_nb) + ' ' + function.name
        pos = nx.kamada_kawai_layout(graph)
        plt.figure()
        nx.draw_networkx_nodes(graph, pos)
        nx.draw_networkx_edges(graph, pos)
        nx.draw_networkx_labels(graph, pos, labels)

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

    def instanciate(self, name, arguments, object_dict):
        """
        """
        for i, arg in enumerate(arguments):
            if type(arg) == str and arg[0] == '#':
                arguments[i] = int(arg[1:])
            elif type(arg) == str and arg[0:2] == '(#':
                argument = []
                arg_id = ""
                for char in arg[1:-1]:
                    if char == ',':
                        argument.append(arg_id)
                        arg_id = ""
                        continue

                    arg_id += char
                argument.append(arg_id)
                arguments[i] = list(argument)

        if name == 'VERTEX_POINT':
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'ORIENTED_EDGE':
            # arguments[4] is the orientation, not taken into account
            volmdlr_object = object_dict[arguments[3]]

        elif name == 'FACE_OUTER_BOUND':
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'FACE_BOUND':
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'SURFACE_CURVE':
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'SEAM_CURVE':
            volmdlr_object = object_dict[arguments[1]]

        # elif name == 'EDGE_CURVE':
        #     object_dict[instanciate_id] = object_dict[arguments[3]]

        elif name == 'VERTEX_LOOP':
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'PCURVE':
            # TODO : Pas besoin de mettre PCURVE ici s'il n'est pas dans STEP_TO_VOLMDLR
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'SHELL_BASED_SURFACE_MODEL':
            volmdlr_object = object_dict[int(arguments[1][0][1:])]
            # Shell3D

        elif name == 'ITEM_DEFINED_TRANSFORMATION':
            volmdlr_object1 = object_dict[arguments[2]]
            volmdlr_object2 = object_dict[arguments[3]]
            # TODO : how to frame map properly from these two Frame3D ?
            # volmdlr_object = volmdlr_object2 - volmdlr_object1
            volmdlr_object = volmdlr_object2
            # Frame3D

        elif name == 'MANIFOLD_SURFACE_SHAPE_REPRESENTATION':
            shells = []
            for arg in arguments[1]:
                if isinstance(object_dict[int(arg[1:])],
                              volmdlr.faces.OpenShell3D):
                    shell = object_dict[int(arg[1:])]
                    shells.append(shell)
            volmdlr_object = shells
            # Shell3D

        elif name == 'MANIFOLD_SOLID_BREP':
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'SHAPE_REPRESENTATION':
            # does it have the extra argument comming from
            # SHAPE_REPRESENTATION_RELATIONSHIP ? In this cas return
            # them
            if len(arguments) == 4:
                shells = object_dict[int(arguments[3])]
                volmdlr_object = shells
            else:
                shells = []
                # frames = []
                for arg in arguments[1]:
                    if int(arg[1:]) in object_dict and \
                            isinstance(object_dict[int(arg[1:])],
                                       volmdlr.faces.OpenShell3D):
                        shells.append(object_dict[int(arg[1:])])
                    elif int(arg[1:]) in object_dict and \
                            isinstance(object_dict[int(arg[1:])],
                                       volmdlr.Frame3D):
                        # TODO: Is there something to read here ?
                        pass
                        # frames.append(object_dict[int(arg[1:])])
                    else:
                        pass
                volmdlr_object = shells

        elif name == 'ADVANCED_BREP_SHAPE_REPRESENTATION':
            shells = []
            for arg in arguments[1]:
                if isinstance(object_dict[int(arg[1:])],
                              volmdlr.faces.OpenShell3D):
                    shells.append(object_dict[int(arg[1:])])
            volmdlr_object = shells

        elif name == 'REPRESENTATION_RELATIONSHIP, REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION, SHAPE_REPRESENTATION_RELATIONSHIP':
            if arguments[2] in object_dict:
                if type(object_dict[arguments[2]]) is list:
                    for shell3d in object_dict[arguments[2]]:
                        frame3d = object_dict[arguments[4]]
                        shell3d.frame_mapping(frame3d, 'old', copy=False)
                        # volmdlr_object = shell3d
                    volmdlr_object = None
                else:
                    shell3d = object_dict[arguments[2]]
                    frame3d = object_dict[arguments[4]]
                    shell3d.frame_mapping(frame3d, 'old', copy=False)
                    # volmdlr_object = shell3d
                    volmdlr_object = None
            else:
                volmdlr_object = None

        elif name == 'BOUNDED_CURVE, B_SPLINE_CURVE, B_SPLINE_CURVE_WITH_KNOTS, CURVE, GEOMETRIC_REPRESENTATION_ITEM, RATIONAL_B_SPLINE_CURVE, REPRESENTATION_ITEM':
            modified_arguments = ['']+arguments
            if modified_arguments[-1] == "''":
                modified_arguments.pop()
            volmdlr_object = STEP_TO_VOLMDLR[name].from_step(
                modified_arguments, object_dict)

        elif name in STEP_TO_VOLMDLR and hasattr(
                STEP_TO_VOLMDLR[name], "from_step"):
            volmdlr_object = STEP_TO_VOLMDLR[name].from_step(
                arguments, object_dict)

        else:
            raise NotImplementedError(
                'Dont know how to interpret {} with args {}'.format(name,
                                                                    arguments))
        return volmdlr_object

    def to_volume_model(self):
        if not self.upd_graph:
            self.graph = self.create_graph()

        object_dict = {}

        self.graph.add_node("#0")
        frame_mapping_nodes = []
        shell_nodes = []
        for node in self.graph.nodes:
            if node != '#0' and self.functions[node].name == 'REPRESENTATION_RELATIONSHIP, REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION, SHAPE_REPRESENTATION_RELATIONSHIP':
                frame_mapping_nodes.append(node)
            if node != '#0' and (self.functions[node].name == "CLOSED_SHELL"
                                 or
                                 self.functions[node].name == "OPEN_SHELL"):
                shell_nodes.append(node)

        frame_mapped_shell_node = []
        for s_node in shell_nodes:
            for fm_node in frame_mapping_nodes:
                if nx.has_path(self.graph, source=fm_node, target=s_node):
                    frame_mapped_shell_node.append(s_node)
                    break
        shell_nodes_copy = shell_nodes.copy()
        [shell_nodes.remove(node) for node in frame_mapped_shell_node]

        for node in shell_nodes + frame_mapping_nodes:
            self.graph.add_edge('#0', node)

        edges = list(
            nx.algorithms.traversal.breadth_first_search.bfs_edges(self.graph,
                                                                   "#0"))[::-1]
        for edge_nb, edge in enumerate(edges):
            instanciate_id = edge[1]
            volmdlr_object = self.instanciate(
                self.functions[instanciate_id].name,
                self.functions[instanciate_id].arg[:],
                object_dict)

            object_dict[instanciate_id] = volmdlr_object

        shells = []
        for node in shell_nodes_copy:
            shells.append(object_dict[node])

        return volmdlr.core.VolumeModel(shells)

    def to_scatter_volume_model(self, name):
        object_dict = {}
        points3d = []
        for stepfunction in self.functions.values():
            if stepfunction.name == 'CARTESIAN_POINT':
                # INSTANCIATION
                name = self.functions[stepfunction.id].name
                arguments = self.functions[stepfunction.id].arg[:]
                for i, arg in enumerate(arguments):
                    if type(arg) == str and arg[0] == '#':
                        arguments[i] = int(arg[1:])
                volmdlr_object = STEP_TO_VOLMDLR[name].from_step(
                    arguments, object_dict)
                points3d.append(volmdlr_object)
        return volmdlr.core.VolumeModel(points3d)

    def plot_data(self):
        if not self.upd_graph:
            self.graph = self.create_graph()
        graph = self.graph.copy()

        graph.remove_nodes_from([stepfunction.id for stepfunction
                                 in self.functions.values()
                                 if stepfunction.name == 'CARTESIAN_POINT'
                                 or stepfunction.name == 'DIRECTION'])
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

    # step subfunctions
    'REPRESENTATION_RELATIONSHIP, REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION, SHAPE_REPRESENTATION_RELATIONSHIP': volmdlr.faces.OpenShell3D.translation,
    'SHELL_BASED_SURFACE_MODEL': None,
    'MANIFOLD_SURFACE_SHAPE_REPRESENTATION': None,
    'MANIFOLD_SOLID_BREP': None,
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
