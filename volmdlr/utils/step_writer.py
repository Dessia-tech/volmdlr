"""
volmdlr utils for translating volmdlr primitives into step file entity.
"""


def product_writer(current_id, product_name):
    """
    Helper function to write the step definitions of a product.
    """
    step_content = ''
    product_definition_context_id = current_id + 1
    step_content += (f"#{product_definition_context_id} = "
                     + "PRODUCT_DEFINITION_CONTEXT('part definition',#2,'design');\n")

    product_context_id = product_definition_context_id + 1
    step_content += f"#{product_context_id} = PRODUCT_CONTEXT('',#2,'mechanical');\n"
    product_id = product_context_id + 1
    step_content += f"#{product_id} = PRODUCT('{product_name}'," \
                    f"'{product_name}','',(#{product_context_id}));\n"
    product_definition_formation_id = product_id + 1
    step_content += f"#{product_definition_formation_id} = " \
                    f"PRODUCT_DEFINITION_FORMATION('','',#{product_id});\n"
    product_definition_id = product_definition_formation_id + 1
    step_content += f"#{product_definition_id} = PRODUCT_DEFINITION('design'," \
                    f"'',#{product_definition_formation_id},#{product_definition_context_id});\n"
    product_definition_shape_id = product_definition_id + 1
    step_content += f"#{product_definition_shape_id} = PRODUCT_DEFINITION_SHAPE(''," \
                    f"'',#{product_definition_id});\n"
    shape_definition_repr_id = product_definition_shape_id + 1
    shape_representation_id = shape_definition_repr_id + 1
    step_content += f"#{shape_definition_repr_id} = SHAPE_DEFINITION_REPRESENTATION(" \
                    f"#{product_definition_shape_id},#{shape_representation_id});\n"

    return step_content, shape_definition_repr_id


def geometric_context_writer(current_id, uncertainty: float = 1e-4):
    """
    Helper function to write the step definitions of a product.
    """
    step_content = ''
    legth_unit_id = current_id + 1
    plane_angle_unit_id = legth_unit_id + 1
    solid_angle_unit = plane_angle_unit_id + 1
    uncertainty_id = solid_angle_unit + 1
    geometric_representation_context_id = uncertainty_id + 1

    step_content += f"#{legth_unit_id} = ( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );\n"
    step_content += f"#{plane_angle_unit_id} = ( NAMED_UNIT(*) PLANE_ANGLE_UNIT() SI_UNIT($,.RADIAN.) );\n"
    step_content += f"#{solid_angle_unit} = ( NAMED_UNIT(*) SI_UNIT($,.STERADIAN.) SOLID_ANGLE_UNIT() );\n"
    step_content += f"#{uncertainty_id} = UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE({uncertainty})," \
                    f"#{legth_unit_id},'distance_accuracy_value','confusion accuracy');\n"
    step_content += f"#{geometric_representation_context_id} = ( GEOMETRIC_REPRESENTATION_CONTEXT(3) " \
                    f"GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#{uncertainty_id})) " \
                    f"GLOBAL_UNIT_ASSIGNED_CONTEXT((#{legth_unit_id},#{plane_angle_unit_id}," \
                    f"#{solid_angle_unit})) " \
                    f"REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY') );\n"

    return step_content, geometric_representation_context_id


def assembly_definition_writer(current_id, assembly_data, component_data, assembly_frame_id, component_frame_id):
    """Returns component assembly data."""
    shape_representation, assembly_product_definition = assembly_data
    component_shape_representation, component_product_definition = component_data
    step_content = ''
    context_dependent_id = current_id + 1
    repr_relationship_with_transfo_id = context_dependent_id + 1
    item_defined_transfo_id = repr_relationship_with_transfo_id + 1
    product_definition_shape_id = item_defined_transfo_id + 1
    next_assembly_usage_occurrence_id = product_definition_shape_id + 1

    step_content += f"#{context_dependent_id} = CONTEXT_DEPENDENT_SHAPE_REPRESENTATION(" \
                    f"#{repr_relationship_with_transfo_id},#{product_definition_shape_id});\n"
    step_content += f"#{repr_relationship_with_transfo_id} = ( REPRESENTATION_RELATIONSHIP('',''," \
                    f"#{component_shape_representation},#{shape_representation}) " \
                    f"REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION(#{item_defined_transfo_id}) " \
                    f"SHAPE_REPRESENTATION_RELATIONSHIP() );\n"
    step_content += f"#{item_defined_transfo_id} = ITEM_DEFINED_TRANSFORMATION('','',#{assembly_frame_id}," \
                    f"#{component_frame_id});\n"
    step_content += f"#{product_definition_shape_id} = PRODUCT_DEFINITION_SHAPE('Placement'," \
                    f"'Placement of an item',#{next_assembly_usage_occurrence_id});\n"
    step_content += f"#{next_assembly_usage_occurrence_id} = NEXT_ASSEMBLY_USAGE_OCCURRENCE(" \
                    f"'','','',#{assembly_product_definition},#{component_product_definition},$);\n"
    return step_content, next_assembly_usage_occurrence_id


def step_ids_to_str(ids):
    """
    Returns a string with a '#' in front of each ID and a comma separating each-one.

    :param ids: A list of step primitives IDs
    :type ids: List[int]
    :return: A string containing all the IDs
    :rtype: str
    """
    return ','.join([f"#{i}" for i in ids])


STEP_HEADER = '''ISO-10303-21;
HEADER;
FILE_DESCRIPTION(('{name}'),'2;1');
FILE_NAME('{filename}','{timestamp}',('Author'),(''),'Volmdlr v{version}','','Unknown');
FILE_SCHEMA(('AUTOMOTIVE_DESIGN {{ 1 0 10303 214 1 1 1 1 }}'));
ENDSEC;
DATA;
#1 = APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#2);
#2 = APPLICATION_CONTEXT('core data for automotive mechanical design processes');
'''

STEP_FOOTER = '''ENDSEC;
END-ISO-10303-21;
'''
