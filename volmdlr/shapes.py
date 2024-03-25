"""volmdlr shells module."""
import base64
# pylint: disable=no-name-in-module
import math
import random
import warnings
import zlib
from io import BytesIO
from itertools import chain, product
from typing import Iterable, List, Tuple, Union, Optional, Any, Dict, overload

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from dessia_common.core import DessiaObject, PhysicalObject
from dessia_common.typings import JsonSerializable
from numpy.typing import NDArray
from trimesh import Trimesh

from OCP.BRep import BRep_Tool, BRep_Builder
from OCP.TopoDS import (TopoDS, TopoDS_Shape, TopoDS_Shell, TopoDS_Face,
                        TopoDS_Solid, TopoDS_CompSolid, TopoDS_Compound, TopoDS_Builder)
from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.GProp import GProp_PGProps
from OCP.BRepGProp import BRepGProp
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_BooleanOperation, BRepAlgoAPI_Splitter
from OCP.BOPAlgo import BOPAlgo_GlueEnum, BOPAlgo_PaveFiller
from OCP.TopTools import TopTools_ListOfShape
from OCP.ShapeFix import ShapeFix_Solid
from OCP.BRepTools import BRepTools

import volmdlr.core_compiled
from volmdlr import curves, display, edges, surfaces, wires, geometry, faces as vm_faces
from volmdlr.core import edge_in_list, get_edge_index_in_list, get_point_index_in_list, point_in_list
from volmdlr.utils.step_writer import geometric_context_writer, product_writer, step_ids_to_str
from volmdlr import from_ocp, to_ocp
from volmdlr.utils.mesh_helpers import perform_decimation

import OCP.TopAbs as topabs  # Topology type enum


def shapetype(obj: TopoDS_Shape) -> topabs.TopAbs_ShapeEnum:
    """
    Gets the shape type for a TopoDS_Shape obejct.
    """

    if obj.IsNull():
        raise ValueError("Null TopoDS_Shape object")

    return obj.ShapeType()


def downcast(obj: TopoDS_Shape) -> TopoDS_Shape:
    """
    Downcasts a TopoDS object to suitable specialized type.
    """
    downcast_LUT = {
        topabs.TopAbs_VERTEX: TopoDS.Vertex_s,
        topabs.TopAbs_EDGE: TopoDS.Edge_s,
        topabs.TopAbs_WIRE: TopoDS.Wire_s,
        topabs.TopAbs_FACE: TopoDS.Face_s,
        topabs.TopAbs_SHELL: TopoDS.Shell_s,
        topabs.TopAbs_SOLID: TopoDS.Solid_s,
        topabs.TopAbs_COMPSOLID: TopoDS.CompSolid_s,
        topabs.TopAbs_COMPOUND: TopoDS.Compound_s,
    }

    f_downcast: Any = downcast_LUT[shapetype(obj)]
    downcasted_obj = f_downcast(obj)

    return downcasted_obj


class Shape(PhysicalObject):
    """
    Represents a shape in the system. Wraps TopoDS_Shape.
    """

    wrapped: TopoDS_Shape

    def __init__(self, obj: TopoDS_Shape, name: str = ""):
        self.wrapped = downcast(obj)
        self.label = name
        PhysicalObject.__init__(self, name=name)

    @classmethod
    def cast(cls, obj: TopoDS_Shape) -> "Shape":
        """
        Returns the right type of wrapper, given a OCCT object.
        """

        tr = None

        # define the shape lookup table for casting
        constructor_LUT = {
            topabs.TopAbs_SHELL: Shell,
            topabs.TopAbs_SOLID: Solid,
            topabs.TopAbs_COMPSOLID: CompSolid,
            topabs.TopAbs_COMPOUND: Compound,
        }

        shape_type = shapetype(obj)
        # NB downcast is needed to handle TopoDS_Shape types
        tr = constructor_LUT[shape_type](downcast(obj))

        return tr


class Shell(Shape):
    """
    OCP shell wrapped.
    """

    wrapped: TopoDS_Shell

    @overload
    def __init__(self, obj: TopoDS_Shell, name: str = '') -> None:
        ...

    @overload
    def __init__(self, faces: List[TopoDS_Face], name: str = '') -> None:
        ...

    def __init__(self, faces: List[vm_faces], name: str = ''):
        self._faces = None
        if isinstance(faces[0], vm_faces.Face3D):
            obj = None
            self._faces = faces
        elif isinstance(faces[0], TopoDS_Face):
            obj = None
        else:
            raise ValueError(
                f"Provided faces value: {faces} is not a valid list of faces."
            )

        self._bbox = None
        Shape.__init__(self, obj, name=name)


class Solid(Shape):
    """
    a single solid
    """

    wrapped: TopoDS_Solid


class CompSolid(Shape):
    """
    A single compsolid.
    """

    wrapped: TopoDS_CompSolid


class Compound(Shape):
    """
    A collection of disconnected solids.
    """

    wrapped: TopoDS_Compound