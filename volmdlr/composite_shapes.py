"""
Define some composite shapes.
"""
from typing import List
import matplotlib.pyplot as plt

import dessia_common.core as dc
import volmdlr
import volmdlr.templates
from volmdlr.utils.step_writer import (product_writer, geometric_context_writer, assembly_definition_writer,
                                       step_ids_to_str)
from volmdlr.core import Primitive3D, BoundingBox, map_primitive_with_initial_and_final_frames
from volmdlr.model import get_babylon_data
from volmdlr import from_ocp
from volmdlr.shells import Shell3D


class Assembly(dc.PhysicalObject):
    """
    Defines an assembly.

    :param components: A list of volmdlr objects
    :type components: List[:class:`volmdlr.core.Primitive3D`]
    :param positions: A list of volmdlr.Frame3D representing the positions of each component in the assembly absolute
        frame.
    :type positions: List[:class:`volmdlr.Frame3D`]
    :param name: The Assembly's name
    :type name: str
    """
    _standalone_in_db = True
    _eq_is_data_eq = True
    _non_serializable_attributes = ['bounding_box', "primitives"]
    _non_data_eq_attributes = ['name', 'bounding_box']
    _non_data_hash_attributes = ['name', 'bounding_box']

    def __init__(self, components: List[Primitive3D], positions: List[volmdlr.Frame3D],
                 frame: volmdlr.Frame3D = volmdlr.OXYZ, name: str = ''):
        self.components = components
        self.frame = frame
        self.positions = positions
        self.primitives = [map_primitive_with_initial_and_final_frames(primitive, frame, frame_primitive)
                           for primitive, frame_primitive in zip(components, positions)]
        self._bbox = None
        dc.PhysicalObject.__init__(self, name=name)

    @property
    def bounding_box(self):
        """
        Returns the bounding box.

        """
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        self._bbox = new_bounding_box

    def _bounding_box(self) -> BoundingBox:
        """
        Computes the bounding box of the model.
        """
        bbox_list = [prim.bounding_box for prim in self.primitives if hasattr(prim, "bounding_box")]
        if not bbox_list:
            return BoundingBox.from_points(self.primitives)
        return BoundingBox.from_bounding_boxes(bbox_list)

    def babylon_data(self, merge_meshes=True):
        """
        Get babylonjs data.

        :return: Dictionary with babylon data.
        """

        return get_babylon_data(self, merge_meshes=merge_meshes)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Assembly.

        side = 'old' or 'new'
        """
        new_positions = [position.frame_mapping(frame, side)
                         for position in self.positions]
        return Assembly(self.components, new_positions, self.frame, self.name)

    def volmdlr_primitives(self):
        """ Return a list of volmdlr primitives to build up an Assembly. """
        return [self]

    def to_step(self, current_id):
        """
        Creates step file entities from volmdlr objects.
        """
        step_content = ''

        product_content, current_id, assembly_data = self.to_step_product(current_id)
        step_content += product_content
        assembly_frames = assembly_data[-1]
        for i, primitive in enumerate(self.components):
            if primitive.__class__.__name__ in ('OpenShell3D', 'ClosedShell3D') or hasattr(primitive, "shell_faces"):
                primitive_content, current_id, primitive_data = primitive.to_step_product(current_id)
                assembly_frame_id = assembly_frames[0]
                component_frame_id = assembly_frames[i + 1]
                assembly_content, current_id = assembly_definition_writer(current_id, assembly_data[:-1],
                                                                          primitive_data, assembly_frame_id,
                                                                          component_frame_id)

            else:
                primitive_content, current_id, primitive_data = primitive.to_step(current_id)
                step_content += primitive_content
                assembly_frame_id = assembly_frames[0]
                component_frame_id = assembly_frames[i + 1]
                assembly_content, current_id = assembly_definition_writer(current_id, assembly_data[:-1],
                                                                          primitive_data, assembly_frame_id,
                                                                          component_frame_id)
            step_content += primitive_content
            step_content += assembly_content

        return step_content, current_id, assembly_data[:-1]

    def plot(self, ax=None, equal_aspect=True):
        """
        Matplotlib plot of model.

        To use for debug.
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d', adjustable='box')
        for primitive in self.primitives:
            primitive.plot(ax)
        if not equal_aspect:
            # ax.set_aspect('equal')
            ax.set_aspect('auto')
        ax.margins(0.1)
        return ax

    def to_step_product(self, current_id):
        """
        Returns step product entities from volmdlr objects.
        """
        step_content = ''
        product_content, shape_definition_repr_id = product_writer(current_id, self.name)
        product_definition_id = shape_definition_repr_id - 2
        step_content += product_content
        shape_representation_id = shape_definition_repr_id + 1
        current_id = shape_representation_id
        assembly_position_content = ''
        frame_ids = []
        for frame in [self.frame] + self.positions:
            frame_content, current_id = frame.to_step(current_id + 1)
            assembly_position_content += frame_content
            frame_ids.append(current_id)

        geometric_context_content, geometric_representation_context_id = geometric_context_writer(current_id)

        step_content += f"#{shape_representation_id} = SHAPE_REPRESENTATION('',({step_ids_to_str(frame_ids)})," \
                        f"#{geometric_representation_context_id});\n"

        step_content += assembly_position_content

        step_content += geometric_context_content

        return step_content, geometric_representation_context_id, \
            [shape_representation_id, product_definition_id, frame_ids]


class Compound(dc.PhysicalObject):
    """
    A class that can be a collection of any volmdlr primitives.
    """

    def __init__(self, primitives, name: str = ""):
        self.primitives = primitives
        self._bbox = None
        self._type = None
        dc.PhysicalObject.__init__(self, name=name)

    @property
    def bounding_box(self):
        """
        Returns the bounding box.

        """
        if not self._bbox:
            self._bbox = self._bounding_box()
        return self._bbox

    @bounding_box.setter
    def bounding_box(self, new_bounding_box):
        """Bounding box setter."""
        self._bbox = new_bounding_box

    @property
    def compound_type(self):
        """
        Returns the compound type.

        """
        if not self._type:
            if all(primitive.__class__.__name__ in ('OpenShell3D', 'ClosedShell3D') or
                   hasattr(primitive, "shell_faces") for primitive in self.primitives):
                self._type = "manifold_solid_brep"
            elif all(isinstance(primitive, (volmdlr.wires.Wire3D, volmdlr.edges.Edge, volmdlr.Point3D))
                     for primitive in self.primitives):
                self._type = "geometric_curve_set"
            else:
                self._type = "shell_based_surface_model"
        return self._type

    @compound_type.setter
    def compound_type(self, value):
        """Compound typesetter."""
        self._type = value

    def _bounding_box(self) -> BoundingBox:
        """
        Computes the bounding box of the model.
        """
        bbox_list = [p.bounding_box for p in self.primitives if hasattr(p, "bounding_box")]
        if not bbox_list:
            return BoundingBox.from_points(self.primitives)
        return BoundingBox.from_bounding_boxes(bbox_list)

    def frame_mapping(self, frame: volmdlr.Frame3D, side: str):
        """
        Changes frame_mapping and return a new Compound.

        side = 'old' or 'new'
        """
        new_primitives = [primitive.frame_mapping(frame, side)
                          for primitive in self.primitives]
        return Compound(new_primitives, self.name)

    def babylon_data(self, merge_meshes=True):
        """
        Get babylonjs data.

        :return: Dictionary with babylon data.
        """

        return get_babylon_data(self, merge_meshes=merge_meshes)

    def volmdlr_primitives(self):
        """Return primitives."""
        return [self]

    def to_step(self, current_id):
        """
        Creates step file entities from volmdlr objects.
        """
        step_content = ''
        primitives_content = ''
        shape_ids = []
        product_content, current_id = product_writer(current_id, self.name)
        product_definition_id = current_id - 2
        step_content += product_content
        brep_id = current_id + 1
        frame_content, frame_id = volmdlr.OXYZ.to_step(brep_id)
        current_id = frame_id

        for primitive in self.primitives:
            primitive_content, current_id = primitive.to_step(current_id)
            primitives_content += primitive_content
            shape_ids.append(current_id)

        geometric_context_content, geometric_representation_context_id = geometric_context_writer(current_id)
        current_id = geometric_representation_context_id
        if self.compound_type == "manifold_solid_brep":
            step_content += f"#{brep_id} = MANIFOLD_SURFACE_SHAPE_REPRESENTATION(''," \
                            f"({step_ids_to_str(shape_ids)})," \
                            f"#{geometric_representation_context_id});\n"
        elif self.compound_type == "geometric_curve_set":
            current_id += 1
            step_content += f"#{brep_id} = GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION(''," \
                            f"(#{current_id})," \
                            f"#{geometric_representation_context_id});\n"

            step_content += f"#{current_id} = GEOMETRIC_SET('',({step_ids_to_str(shape_ids)}));\n"
        step_content += frame_content
        step_content += primitives_content
        step_content += geometric_context_content

        return step_content, current_id, [brep_id, product_definition_id]

    @classmethod
    def from_ocp(cls, occt_compound):
        """
        Creates a volmdlr Compound from an OCCT compound.
        """
        # TODO: Today only shells are read. Needs to implement Solid or wireframe models if needed
        shell_list = [Shell3D.from_ocp(occt_shell) for occt_shell in from_ocp.get_shells(occt_compound)]
        return cls(shell_list)
