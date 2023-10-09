"""
volmdlr utils for importing step files.
"""
import numpy as npy
import volmdlr
import volmdlr.shells as vmshells
from volmdlr import surfaces


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
    function_arg = function_arg.strip()
    if len(function_arg) > 0 and function_arg[-1] != ')':
        function_arg += ')'
    arguments = []
    argument = ""
    if len(function_arg) > 0 and function_arg[0] == "(":
        function_arg += ")"
    parenthesis = 1
    is_str = False
    for char in function_arg:
        if char == "(":
            parenthesis += 1

        if char == "'" and not is_str:
            is_str = True
        elif char == "'" and is_str:
            is_str = False
        # if char != "," or parenthesis > 1 or is_str:
        #     argument += char
        if not is_str and char == " ":
            continue
        if parenthesis > 1 or is_str:
            argument += char
        elif char != ",":
            argument += char
        else:
            arguments.append(argument)
            argument = ""

        if char == ")":
            parenthesis -= 1
            if parenthesis == 0:
                arguments.append(argument[:-1])
                break
    return arguments


def uncertainty_measure_with_unit(arguments, object_dict):
    """
    Gets the length uncertainty related to the shape representation.

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
    Calculates the step file's SI unit conversion factor.

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


def named_unit_plane_angle_unit_si_unit(arguments, *args, **kwargs):
    """
    Returns the dimension of plane angle measure.

    :param arguments: step primitive arguments
    :return: SI unit dimension.
    """
    return SI_PREFIX[arguments[1]]


def named_unit_si_unit_solid_angle_unit(arguments, *args, **kwargs):
    """
    Returns the dimension of solid angle measure.

    :param arguments: step primitive arguments
    :return: SI unit dimension.
    """
    return SI_PREFIX[arguments[1]]


def named_unit_length_unit_si_unit(arguments, *args, **kwargs):
    """
    Returns the dimension of solid angle measure.

    :param arguments: step primitive arguments
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


def length_unit_named_unit_si_unit(arguments, *args, **kwargs):
    """
    Gets the length si unit.

    :param arguments: step primitive arguments
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
        # This can happen when the edge is too small
        return None
    edge_orientation = arguments[4]
    if edge_orientation == '.T.':
        return object_dict[arguments[3]]
    return object_dict[arguments[3]].reverse()


def face_outer_bound(arguments, object_dict):
    """
    Returns the data in case of a FACE_OUTER_BOUND.

    :param arguments: list containing the arguments of the FACE_OUTER_BOUND entity.
    :type arguments: list
    :param object_dict: Dictionary containing the objects already instantiated that will be used as arguments to the
        face_outer_bound entity.
    :type object_dict: dict
    :return: A Contour3D representing the BREP of a face.
    :rtype: volmdlr.wires.Contour3D
    """
    if arguments[2] == '.T.':
        return object_dict[arguments[1]]
    return object_dict[arguments[1]].invert()


def face_bound(arguments, object_dict):
    """
    Returns the data in case of a FACE_BOUND.

    :param arguments: list containing the arguments of the FACE_BOUND entity.
    :type arguments: list
    :param object_dict: Dictionary containing the objects already instantiated that will be used as arguments to the
        face_outer_bound entity.
    :type object_dict: dict
    :return: A Contour3D representing the BREP of a face.
    :rtype: volmdlr.wires.Contour3D
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
    if curve.__class__.__name__ == "Line3D" and point1.is_close(point2):
        return None
    return curve.trim(point1=point1, point2=point2)


def vertex_loop(arguments, object_dict):
    """
    Returns the data in case of a VERTEX_LOOP.
    """
    return object_dict[arguments[1]]


def composite_curve_segment(arguments, object_dict):
    """
    Returns the data in case of a COMPOSITE_CURVE_SEGMENT.
    """
    # arguments[0] = trasition_code (unused)
    # The transition_code type conveys the continuity properties of a composite curve or surface.
    # The continuity referred to is geometric, not parametric continuity.
    # arguments[1] = same_sense : BOOLEAN;
    # arguments[2] = parent_curve : curve;
    edge = object_dict[arguments[2]]
    if arguments[1] == ".F.":
        edge = edge.reverse()
    return edge


def composite_curve(arguments, object_dict):
    """
    Returns the data in case of a COMPOSITE_CURVE.
    """
    name = arguments[0]
    list_primitives = [object_dict[int(arg[1:])]for arg in arguments[1]]
    first_primitive = list_primitives[0]
    last_primitive = list_primitives[-1]
    if first_primitive.start.is_close(last_primitive.end):
        return volmdlr.wires.Contour3D(list_primitives, name)
    return volmdlr.wires.Wire3D(list_primitives, name)


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


def geometric_set(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    # TODO: IS THIS RIGHT?
    primitives = [object_dict[int(node[1:])]
                  for node in arguments[1]]
    return primitives


def shell_based_surface_model(arguments, object_dict):
    """
    Returns the data in case of a Shell3D.
    """
    if len(arguments[1]) == 1:
        return object_dict[int(arguments[1][0][1:])]
    primitives = [object_dict[int(arg[1:])] for arg in arguments[1]]
    compound = volmdlr.core.Compound(primitives)
    compound.compound_type = "manifold_solid_brep"
    return compound


def oriented_closed_shell(arguments, object_dict):
    """
    Returns the data in case of a Shell3D.
    """
    # TODO: How to use the orientation (arguments[3]
    return object_dict[arguments[2]]


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
    primitives = []
    for arg in arguments[1]:
        primitive = object_dict[int(arg[1:])]
        if isinstance(primitive, vmshells.Shell3D):
            primitives.append(primitive)
        if isinstance(primitive, volmdlr.core.Compound):
            counter = 0
            for sub_prim in primitive.primitives:
                sub_prim.name = arguments[0][1:-1] + str(counter)
                counter += 1
            primitives.append(primitive)
    if len(primitives) == 1:
        return primitives[0]
    compound = volmdlr.core.Compound(primitives)
    compound.compound_type = "manifold_solid_brep"
    return compound


def faceted_brep(arguments, object_dict):
    """
    Returns the data in case of a faceted_brep entity, interpreted as shell3D.
    """
    return object_dict[arguments[1]]


def faceted_brep_shape_representation(arguments, object_dict):
    """
    Returns the data in case of a faceted_brep_shape_representation, interpreted as shell3D.
    """
    if len(arguments[1]) == 1:
        return object_dict[int(arguments[1][0][1:])]
    shells = []
    for arg in arguments[1]:
        if isinstance(object_dict[int(arg[1:])],
                      vmshells.Shell3D):
            shell = object_dict[int(arg[1:])]
            shells.append(shell)
    return volmdlr.core.Compound(shells)


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
    # does it have the extra argument coming from
    # SHAPE_REPRESENTATION_RELATIONSHIP ? In this case return
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
                           vmshells.Shell3D):
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
    primitives = []
    for arg in arguments[1]:
        primitive = object_dict[int(arg[1:])]
        if isinstance(primitive, vmshells.Shell3D):
            primitives.append(primitive)
        if isinstance(primitive, volmdlr.core.Compound):
            counter = 0
            for sub_prim in primitive.primitives:
                sub_prim.name = arguments[0][1:-1] + str(counter)
                counter += 1
            primitives.append(primitive)
    if len(primitives) == 1:
        return primitives[0]
    compound = volmdlr.core.Compound(primitives)
    compound.compound_type = "manifold_solid_brep"
    return compound


def geometrically_bounded_surface_shape_representation(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    primitives = []
    for arg in arguments[1]:
        primitives.extend(object_dict[int(arg[1:])])
    if len(primitives) > 1:
        compound = volmdlr.core.Compound(primitives, name=arguments[0])
        compound.compound_type = "geometric_curve_set"
        return compound
    return primitives[0]


def geometrically_bounded_wireframe_shape_representation(arguments, object_dict):
    """
    Returns xx.

    :param arguments: DESCRIPTION
    :type arguments: TYPE
    :param object_dict: DESCRIPTION
    :type object_dict: TYPE
    :return: DESCRIPTION
    :rtype: TYPE

    """
    primitives = []
    for arg in arguments[1]:
        primitives.extend(object_dict[int(arg[1:])])
    if len(primitives) > 1:
        compound = volmdlr.core.Compound(primitives, name=arguments[0])
        compound.compound_type = "geometric_curve_set"
        return compound
    return primitives[0]


def frame_map_closed_shell(closed_shells, item_defined_transformation_frames, shape_representation_frames):
    """
    Frame maps a closed shell in an assembly to its good position.

    :param closed_shells: DESCRIPTION
    :type closed_shells: vmshells.OpenShell3D
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
        matrix_a = npy.array([[basis_a.vectors[0].x, basis_a.vectors[0].y, basis_a.vectors[0].z],
                       [basis_a.vectors[1].x, basis_a.vectors[1].y, basis_a.vectors[1].z],
                       [basis_a.vectors[2].x, basis_a.vectors[2].y, basis_a.vectors[2].z]])
        matrix_b = npy.array([[basis_b.vectors[0].x, basis_b.vectors[0].y, basis_b.vectors[0].z],
                       [basis_b.vectors[1].x, basis_b.vectors[1].y, basis_b.vectors[1].z],
                       [basis_b.vectors[2].x, basis_b.vectors[2].y, basis_b.vectors[2].z]])
        transfer_matrix = npy.linalg.solve(matrix_a, matrix_b)
        new_frame = volmdlr.Frame3D(transformed_frame.origin, volmdlr.Vector3D(*transfer_matrix[0]),
                                    volmdlr.Vector3D(*transfer_matrix[1]), volmdlr.Vector3D(*transfer_matrix[2]))
        new_closedshells.append(shell3d.frame_mapping(new_frame, 'old'))
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


def b_spline_curve_b_spline_curve_with_knots_rational_b_spline_curve_bounded_curve_representation_item_geometric_representation_item_curve(
        arguments, object_dict):
    """
    Bounded b spline with knots curve geometric representation item. To clarify.
    """
    return bounded_curve_b_spline_curve_b_spline_curve_with_knots_curve_geometric_representation_item_rational_b_spline_curve_representation_item(arguments, object_dict)


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
    return bounded_surface_b_spline_surface_b_spline_surface_with_knots_geometric_representation_item_rational_b_spline_surface_representation_item_surface(
        arguments, object_dict)


def b_spline_surface_b_spline_surface_with_knots_rational_b_spline_surface_bounded_surface_representation_item_geometric_representation_item_surface(arguments, object_dict):
    """
    Bounded b spline surface with knots curve geometric representation item. To clarify.
    """
    return bounded_surface_b_spline_surface_b_spline_surface_with_knots_geometric_representation_item_rational_b_spline_surface_representation_item_surface(
        arguments, object_dict)

def product_definition_shape(arguments, object_dict):
    """
    Returns the data in case of a product_definition_shape.
    """
    return object_dict[arguments[2]]


def product_definition(arguments, object_dict):
    """
    Returns the data in case of a product_definition.
    """
    return object_dict[arguments[2]]


def product_definition_formation(arguments, object_dict):
    """
    Returns the data in case of a product_definition_formation.
    """
    return object_dict[arguments[2]]


def product_definition_formation_with_specified_source(arguments, object_dict):
    """
    Returns the data in case of a product_definition_formation_with_specified_source.
    """
    return object_dict[arguments[2]]


def product(arguments, *args, **kwargs):
    """
    Returns the data in case of a product.
    """
    return arguments[0]


def application_context(arguments, *args, **kwargs):
    """
    Returns the data in case of an application_context.
    """
    return arguments[0]


def product_context(arguments, *args, **kwargs):
    """
    Returns the data in case of a product_context.
    """
    return arguments

STEP_TO_VOLMDLR = {
    # GEOMETRICAL ENTITIES
    'CARTESIAN_POINT': volmdlr.Point3D,
    'DIRECTION': volmdlr.Vector3D,
    'VECTOR': volmdlr.Vector3D,

    'AXIS1_PLACEMENT': None,
    'AXIS2_PLACEMENT_2D': None,  # ??????????????????
    'AXIS2_PLACEMENT_3D': volmdlr.Frame3D,

    'LINE': volmdlr.curves.Line3D,  # LineSegment3D,
    'CIRCLE': volmdlr.curves.Circle3D,
    'ELLIPSE': volmdlr.curves.Ellipse3D,
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

    'PLANE': surfaces.Plane3D,
    'CYLINDRICAL_SURFACE': surfaces.CylindricalSurface3D,
    'CONICAL_SURFACE': surfaces.ConicalSurface3D,
    'SPHERICAL_SURFACE': surfaces.SphericalSurface3D,
    'TOROIDAL_SURFACE': surfaces.ToroidalSurface3D,
    'DEGENERATE_TOROIDAL_SURFACE': surfaces.ToroidalSurface3D,
    'B_SPLINE_SURFACE_WITH_KNOTS': surfaces.BSplineSurface3D,
    'B_SPLINE_SURFACE': surfaces.BSplineSurface3D,
    'BEZIER_SURFACE': surfaces.BSplineSurface3D,

    'OFFSET_SURFACE': None,
    'SURFACE_REPLICA': None,
    'RATIONAL_B_SPLINE_SURFACE': surfaces.BSplineSurface3D,
    'RECTANGULAR_TRIMMED_SURFACE': None,
    'SURFACE_OF_LINEAR_EXTRUSION': surfaces.ExtrusionSurface3D,
    # CAN BE A BSplineSurface3D
    'SURFACE_OF_REVOLUTION': surfaces.RevolutionSurface3D,
    'UNIFORM_SURFACE': surfaces.BSplineSurface3D,
    'QUASI_UNIFORM_SURFACE': surfaces.BSplineSurface3D,
    'RECTANGULAR_COMPOSITE_SURFACE': volmdlr.faces.PlaneFace3D,  # TOPOLOGICAL FACES
    'CURVE_BOUNDED_SURFACE': volmdlr.faces.PlaneFace3D,  # TOPOLOGICAL FACE

    # Bsplines
    'BOUNDED_SURFACE, B_SPLINE_SURFACE, B_SPLINE_SURFACE_WITH_KNOTS, GEOMETRIC_REPRESENTATION_ITEM,'
    ' RATIONAL_B_SPLINE_SURFACE, REPRESENTATION_ITEM, SURFACE': surfaces.BSplineSurface3D,
    "BOUNDED_SURFACE, B_SPLINE_SURFACE, B_SPLINE_SURFACE_WITH_KNOTS, SURFACE, GEOMETRIC_REPRESENTATION_ITEM,"
    " RATIONAL_B_SPLINE_SURFACE, REPRESENTATION_ITEM": surfaces.BSplineSurface3D,
    # TOPOLOGICAL ENTITIES
    'VERTEX_POINT': None,

    'EDGE_CURVE': volmdlr.edges.Edge,  # LineSegment3D, # TOPOLOGICAL EDGE
    'ORIENTED_EDGE': None,  # TOPOLOGICAL EDGE
    "GEOMETRIC_SET": None,
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

    'CLOSED_SHELL': vmshells.ClosedShell3D,
    'OPEN_SHELL': vmshells.OpenShell3D,
    #        'ORIENTED_CLOSED_SHELL': None,
    'CONNECTED_FACE_SET': vmshells.OpenShell3D,
    'GEOMETRIC_CURVE_SET': None,

    # step sub-functions

    'UNCERTAINTY_MEASURE_WITH_UNIT': None,
    'CONVERSION_BASED_UNIT, LENGTH_UNIT, NAMED_UNIT': None,
    'LENGTH_MEASURE_WITH_UNIT': None,
    'LENGTH_UNIT, NAMED_UNIT, SI_UNIT': None,
    'PLANE_ANGLE_MEASURE_WITH_UNIT': None,
    'NAMED_UNIT, PLANE_ANGLE_UNIT, SI_UNIT': None,
    'CONVERSION_BASED_UNIT, NAMED_UNIT, PLANE_ANGLE_UNIT': None,
    'GEOMETRIC_REPRESENTATION_CONTEXT, GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT, GLOBAL_UNIT_ASSIGNED_CONTEXT, REPRESENTATION_CONTEXT': None,
    'REPRESENTATION_RELATIONSHIP, REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION, SHAPE_REPRESENTATION_RELATIONSHIP': vmshells.OpenShell3D.translation,
    'SHELL_BASED_SURFACE_MODEL': None,
    'MANIFOLD_SURFACE_SHAPE_REPRESENTATION': None,
    'MANIFOLD_SOLID_BREP': None,
    'BREP_WITH_VOIDS': None,
    'SHAPE_REPRESENTATION': None,
    'ADVANCED_BREP_SHAPE_REPRESENTATION': None,
    "FACETED_BREP_SHAPE_REPRESENTATION": None,
    "GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION": None,
    "GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION": None,
    "EDGE_BASED_WIREFRAME_SHAPE_REPRESENTATION": None,
    'ITEM_DEFINED_TRANSFORMATION': None,
    'SHAPE_REPRESENTATION_RELATIONSHIP': None,
    "NEXT_ASSEMBLY_USAGE_OCCURRENCE": None,

    'BOUNDED_CURVE, B_SPLINE_CURVE, B_SPLINE_CURVE_WITH_KNOTS, CURVE, GEOMETRIC_REPRESENTATION_ITEM, RATIONAL_B_SPLINE_CURVE, REPRESENTATION_ITEM': volmdlr.edges.BSplineCurve3D,
    "APPLICATION_CONTEXT": None,
    "PRODUCT_DEFINITION_SHAPE": None,
    "PRODUCT_DEFINITION": None,
    "PRODUCT_DEFINITION_FORMATION": None,
    "PRODUCT": None,
}

SI_PREFIX = {'.EXA.': 1e18, '.PETA.': 1e15, '.TERA.': 1e12, '.GIGA.': 1e9, '.MEGA.': 1e6, '.KILO.': 1e3,
             '.HECTO.': 1e2, '.DECA.': 1e1, '$': 1, '.DECI.': 1e-1, '.CENTI.': 1e-2, '.MILLI.': 1e-3, '.MICRO.': 1e-6,
             '.NANO.': 1e-9, '.PICO.': 1e-12, '.FEMTO.': 1e-15, '.ATTO.': 1e-18}

STEP_REPRESENTATION_ENTITIES = {"ADVANCED_BREP_SHAPE_REPRESENTATION", "FACETED_BREP_SHAPE_REPRESENTATION",
                                "MANIFOLD_SURFACE_SHAPE_REPRESENTATION",
                                "GEOMETRICALLY_BOUNDED_WIREFRAME_SHAPE_REPRESENTATION",
                                "GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION",
                                "EDGE_BASED_WIREFRAME_SHAPE_REPRESENTATION"
                                }

VOLMDLR_TO_STEP = {}
for k, v in STEP_TO_VOLMDLR.items():
    if v:
        if v in VOLMDLR_TO_STEP:
            VOLMDLR_TO_STEP[v].append(k)
        else:
            VOLMDLR_TO_STEP[v] = [k]
