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
            if arg == []:
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
            if function_name == "":
                arguments = self.step_subfunctions(arguments)

            for i, argument in enumerate(arguments):
                if argument[:2] == '(#' and argument[-1] == ')':
                    arg_list = volmdlr.core.set_to_list(argument)
                    arguments[i] = arg_list

            function = StepFunction(function_id, function_name, arguments)
            functions[function_id] = function

        f.close()

        return functions, all_connections

    def create_graph(self, draw=False, html=False):

        G = nx.Graph()
        F = nx.DiGraph()
        labels = {}

        for function in self.functions.values():
            if function.name in STEP_TO_VOLMDLR:
                G.add_node(function.id)
                F.add_node(function.id)
                labels[function.id] = str(function.id) + ' ' + function.name

        # Delete connection if node not found
        node_list = list(F.nodes())
        delete_connection = []
        for connection in self.all_connections:
            if connection[0] not in node_list or connection[
                1] not in node_list:
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

    def instanciate(self, instanciate_id, object_dict):
        """
        Returns None if the object was instanciate
        """

        # print('instanciate_id', instanciate_id)

        name = self.functions[instanciate_id].name
        arguments = self.functions[instanciate_id].arg[:]

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
            #            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]

        # elif name == 'LINE':
        #     pass

        elif name == 'ORIENTED_EDGE':
            #            object_dict[instanciate_id] = object_dict[arguments[3]]
            volmdlr_object = object_dict[arguments[3]]

        elif name == 'FACE_OUTER_BOUND':
            #            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'FACE_BOUND':
            #            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'SURFACE_CURVE':
            #            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'SEAM_CURVE':
            #            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]
        # elif name == 'EDGE_CURVE':
        #     object_dict[instanciate_id] = object_dict[arguments[3]]

        elif name == 'VERTEX_LOOP':
            #            object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]

        elif name == 'PCURVE':
            #     # object_dict[instanciate_id] = object_dict[arguments[1]]
            volmdlr_object = object_dict[arguments[1]]

        elif name in STEP_TO_VOLMDLR and hasattr(
                STEP_TO_VOLMDLR[name], "from_step"):
            # print(object_dict)
            volmdlr_object = STEP_TO_VOLMDLR[name].from_step(
                arguments, object_dict)

        #            object_dict[instanciate_id] = volmdlr_object
        #            if hasattr(volmdlr_object, "primitive"):
        #                primitives.append(volmdlr_object.primitive)primitives

        else:
            print('name', name)
            print('arguments', arguments)
            raise NotImplementedError

        return volmdlr_object

    def to_volume_model(self):
        self.graph = self.create_graph()

        object_dict = {}

        self.graph.add_node("#0")
        for node in self.graph.nodes:
            if node != '#0' and (self.functions[node].name == "CLOSED_SHELL" or
                                 self.functions[node].name == "OPEN_SHELL"):
                self.graph.add_edge("#0", node)

        edges = list(
            nx.algorithms.traversal.breadth_first_search.bfs_edges(self.graph,
                                                                   "#0"))[::-1]
        for edge_nb, edge in enumerate(edges):
            instanciate_id = edge[1]
            volmdlr_object = self.instanciate(instanciate_id, object_dict)

            object_dict[instanciate_id] = volmdlr_object

        shells = []
        for node in list(self.graph.nodes):
            if node != '#0' and (self.functions[node].name == 'CLOSED_SHELL' or
                                 self.functions[node].name == "OPEN_SHELL"):
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
    #        'PCURVE': None,
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
}

VOLMDLR_TO_STEP = {}
for k,v in STEP_TO_VOLMDLR.items():
    if v:
        if v in VOLMDLR_TO_STEP:
            VOLMDLR_TO_STEP[v].append(k)
        else:
            VOLMDLR_TO_STEP[v] = [k]
