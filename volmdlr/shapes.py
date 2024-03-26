"""volmdlr shapes module."""
import base64
# pylint: disable=no-name-in-module
import math
import random
import sys
import warnings
import zlib
from io import BytesIO
from itertools import chain, product
from typing import Iterable, List, Tuple, Union, Optional, Any, Dict, overload, Literal, cast as tcast

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
from dessia_common.files import BinaryFile
from dessia_common.core import DessiaObject, PhysicalObject
from dessia_common.typings import JsonSerializable
from numpy.typing import NDArray
from trimesh import Trimesh

from OCP.BRep import BRep_Tool, BRep_Builder
from OCP.TopoDS import (TopoDS, TopoDS_Shape, TopoDS_Shell, TopoDS_Face,
                        TopoDS_Solid, TopoDS_CompSolid, TopoDS_Compound, TopoDS_Builder)
from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing, BRepBuilderAPI_Copy
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.GProp import GProp_PGProps
from OCP.BRepGProp import BRepGProp
from OCP.BRepExtrema import BRepExtrema_DistShapeShape
from OCP.BRepAlgoAPI import (BRepAlgoAPI_Fuse, BRepAlgoAPI_BooleanOperation, BRepAlgoAPI_Splitter,
                             BRepAlgoAPI_Cut, BRepAlgoAPI_Common)
from OCP.BOPAlgo import BOPAlgo_GlueEnum, BOPAlgo_PaveFiller
from OCP.TopTools import TopTools_ListOfShape
from OCP.ShapeFix import ShapeFix_Solid
from OCP.BRepTools import BRepTools
from OCP.TopTools import TopTools_IndexedMapOfShape
from OCP.TopExp import TopExp

import volmdlr.core_compiled
from volmdlr import curves, display, edges, surfaces, wires, geometry, faces as vm_faces
from volmdlr.core import edge_in_list, get_edge_index_in_list, get_point_index_in_list
from volmdlr.utils.step_writer import geometric_context_writer, product_writer, step_ids_to_str
from volmdlr import from_ocp, to_ocp
from volmdlr.utils.mesh_helpers import perform_decimation

import OCP.TopAbs as topabs  # Topology type enum

shape_LUT = {
    topabs.TopAbs_VERTEX: "Vertex",
    topabs.TopAbs_EDGE: "Edge",
    topabs.TopAbs_WIRE: "Wire",
    topabs.TopAbs_FACE: "Face",
    topabs.TopAbs_SHELL: "Shell",
    topabs.TopAbs_SOLID: "Solid",
    topabs.TopAbs_COMPSOLID: "CompSolid",
    topabs.TopAbs_COMPOUND: "Compound",
}

inverse_shape_LUT = {v: k for k, v in shape_LUT.items()}

Shapes = Literal[
    "Vertex", "Edge", "Wire", "Face", "Shell", "Solid", "CompSolid", "Compound"]

# pylint: disable=no-name-in-module,invalid-name,unused-import,wrong-import-order


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
    _non_serializable_attributes = ["obj"]
    _non_data_eq_attributes = ["wrapped"]
    wrapped: TopoDS_Shape

    def __init__(self, obj: TopoDS_Shape, name: str = ""):
        self.wrapped = downcast(obj)
        self.label = name
        self._bbox = None
        PhysicalObject.__init__(self, name=name)

    def copy(self, deep=True, memo=None):
        """
        Copy of Shape.

        :return: return a copy the Shape.
        """
        return self.__class__(obj=BRepBuilderAPI_Copy(self.wrapped, True, False).Shape())
        # new_faces = [face.copy(deep=deep, memo=memo) for face in self.faces]
        # return self.__class__(new_faces, color=self.color, alpha=self.alpha,
        #                       name=self.name)

    @classmethod
    def cast(cls, obj: TopoDS_Shape, name: str = '') -> "Shape":
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

        shape_type = shapetype(obj=obj)
        # NB downcast is needed to handle TopoDS_Shape types
        tr = constructor_LUT[shape_type](obj=downcast(obj), name=name)

        return tr

    @staticmethod
    def _entities(obj, topo_type: Shapes) -> Iterable[TopoDS_Shape]:
        """Gets shape's entities (vertices, edges, faces, shells...)."""
        shape_set = TopTools_IndexedMapOfShape()
        TopExp.MapShapes_s(obj, inverse_shape_LUT[topo_type], shape_set)

        return tcast(Iterable[TopoDS_Shape], shape_set)

    def _get_vertices(self):
        """Gets shape's vertices, if there exists any."""
        return [downcast(obj=i) for i in self._entities(obj=self.wrapped, topo_type="Vertex")]

    def _get_edges(self):
        """Gets shape's edges, if there exists any."""
        return [downcast(i) for i in self._entities(obj=self.wrapped, topo_type="Edge") if
                not BRep_Tool.Degenerated_s(TopoDS.Edge_s(i))]

    def _get_faces(self):
        """Gets shape's faces, if there exists any."""
        return [downcast(i) for i in self._entities(obj=self.wrapped, topo_type="Face")]

    def get_shells(self) -> List["Shell"]:
        """
        :returns: All the shells in this Shape.
        """

        return [Shell(obj=i) for i in self._entities(obj=self.wrapped, topo_type="Shell")]

    def get_solids(self) -> List["Solid"]:
        """
        :returns: All the solids in this Shape.
        """

        return [Solid(obj=i) for i in self._entities(obj=self.wrapped, topo_type="Solid")]

    def get_compsolids(self) -> List["CompSolid"]:
        """
        :returns: All the compsolids in this Shape.
        """

        return [CompSolid(obj=i) for i in self._entities(obj=self.wrapped, topo_type="CompSolid")]

    def to_brep(self, file: Union[str, BytesIO]) -> bool:
        """
        Export this shape to a BREP file.
        """

        rv = BRepTools.Write_s(self.wrapped, file)

        return True if rv is None else rv

    @classmethod
    def from_brep(cls, file: Union[str, BytesIO], name: str = '') -> "Shape":
        """
        Import shape from a BREP file.
        """
        shape = TopoDS_Shape()
        builder = BRep_Builder()

        BRepTools.Read_s(shape, file, builder)

        if shape.IsNull():
            raise ValueError(f"Could not import {file}")
        shape = cls.cast(obj=shape, name=name)
        # for shape_type in ["shells", "solids", "compsolids"]:
        #     entity = getattr(shape, f"get_{shape_type}")()
        #     print(True)

        if cls.__name__ in ["Shell", "Solid", "CompSolid"]:
            return getattr(shape, f"get_{cls.__name__.lower()}s")()[0]
        return shape

    @classmethod
    def from_brep_stream(cls, stream: BinaryFile, name: str = "") -> "Shape":
        """
        Import shape from a BREP file stream.
        """
        return cls.from_brep(file=stream, name=name)

    def to_brep_stream(self) -> BytesIO:
        """
        Export shape from a BREP file stream.
        """
        brep_bytesio = BytesIO()
        self.to_brep(brep_bytesio)
        return brep_bytesio

    def to_dict(self,
                use_pointers: bool = True,
                memo=None,
                path: str = "#",
                id_method=True,
                id_memo=None, **kwargs) -> JsonSerializable:
        """
        Serializes a 3-dimensional Shape into a dictionary.
        """
        dict_ = self.base_dict()

        brep_content = self.to_brep_stream().getvalue()
        compressed_brep_data = zlib.compress(brep_content)
        encoded_brep_string = base64.b64encode(compressed_brep_data).decode()

        dict_["brep"] = encoded_brep_string

        return dict_

    @classmethod
    def dict_to_object(
            cls,
            dict_: JsonSerializable,
            force_generic: bool = False,
            global_dict=None,
            pointers_memo: Dict[str, Any] = None,
            path: str = "#",
    ) -> "Shape":
        """
        Creates a Shape from a dictionary.
        """
        name = dict_["name"]

        encoded_brep_string = dict_["brep"]
        decoded_brep_data = base64.b64decode(encoded_brep_string)
        decompressed_brep_data = zlib.decompress(decoded_brep_data)
        new_brep_bytesio = BytesIO(decompressed_brep_data)
        obj_class = getattr(sys.modules[__name__], dict_["object_class"][15:])
        return obj_class.from_brep(new_brep_bytesio, name)

    def bounding_box(self):
        """Gets bounding box for this shape."""
        if not self._bbox:
            tol = 1e-2
            bbox = Bnd_Box()

            mesh = BRepMesh_IncrementalMesh(self.wrapped, tol, True)
            mesh.Perform()

            BRepBndLib.Add_s(self.wrapped, bbox, True)

            xmin, ymin, zmin, xmax, ymax, zmax = bbox.Get()

            self._bbox = volmdlr.core.BoundingBox(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, zmin=zmin, zmax=zmax)
        return self._bbox

    def volume(self):
        """
        Gets the Volume of a shape.

        :return:
        """
        tol = 1e-6
        prop = GProp_PGProps()
        BRepGProp.VolumeProperties_s(self.wrapped, prop, tol)
        return abs(prop.Mass())

    @staticmethod
    def _bool_op(
            args: Iterable["Shape"],
            tools: Iterable["Shape"],
            operation: Union[BRepAlgoAPI_BooleanOperation, BRepAlgoAPI_Splitter],
            parallel: bool = True,
    ) -> "TopoDS_Shell":
        """
        Generic boolean operation
        :param parallel: Sets the SetRunParallel flag, which enables parallel execution of boolean\
        operations in OCC kernel.
        """

        arg = TopTools_ListOfShape()
        for obj in args:
            arg.Append(obj.wrapped)

        tool = TopTools_ListOfShape()
        for obj in tools:
            tool.Append(obj.wrapped)

        operation.SetArguments(arg)
        operation.SetTools(tool)

        operation.SetRunParallel(parallel)
        operation.Build()

        return operation.Shape()

    def subtraction(self, *to_subtract: "Shape", tol: Optional[float] = None) -> "Shape":
        """
        Subtract the positional arguments from this Shape.

        :param tol: Fuzzy mode tolerance
        """

        cut_op = BRepAlgoAPI_Cut()

        if tol:
            cut_op.SetFuzzyValue(tol)

        return self.__class__(self._bool_op((self,), to_subtract, cut_op))

    def union(self, *to_union: "Shape", glue: bool = False, tol: Optional[float] = None):
        """
        Fuse the positional arguments with this Shape.
        :param glue: Sets the glue option for the algorithm, which allows
            increasing performance of the intersection of the input shapes
        :param tol: Fuzzy mode tolerance
        """

        fuse_op = BRepAlgoAPI_Fuse()
        if glue:
            fuse_op.SetGlue(BOPAlgo_GlueEnum.BOPAlgo_GlueShift)
        if tol:
            fuse_op.SetFuzzyValue(tol)

        union = self._bool_op((self,), to_union, fuse_op)

        return self.__class__(union)

    def intersection(self, *to_intersect: "Shape", tol: Optional[float] = None) -> "Shape":
        """
        Intersection of the positional arguments and this Shape.

        :param tol: Fuzzy mode tolerance
        """

        intersect_op = BRepAlgoAPI_Common()

        if tol:
            intersect_op.SetFuzzyValue(tol)

        return self.__class__(self._bool_op((self,), to_intersect, intersect_op))


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

    @overload
    def __init__(self, faces: List[vm_faces.Face3D], name: str = '') -> None:
        ...

    def __init__(self, faces: List[vm_faces.Face3D] = None, name: str = '', obj=None):
        self._faces = None
        if faces:
            obj = self._from_faces(faces)
            if isinstance(faces[0], vm_faces.Face3D):
                self._faces = faces
        Shape.__init__(self, obj, name=name)

    @staticmethod
    def _from_faces(faces):
        """
        Helper method to create a TopoDS_Shell from a list of faces.

        :param faces: list of faces volmdlr or a list of faces TopoDS_Face.
        :return: TopoDS_Shell object.
        """
        if isinstance(faces[0], vm_faces.Face3D):
            faces = [face.to_ocp() for face in faces]

        shell_builder = BRepBuilderAPI_Sewing()

        for face in faces:
            shell_builder.Add(face)

        shell_builder.Perform()
        return shell_builder.SewedShape()

    @property
    def faces(self):
        """Get shell's volmdlr faces."""
        if not self._faces:
            pass
            # self._faces = [from_ocp. for face in self._get_faces(self.wrapped)]
        return self._faces

    @faces.setter
    def faces(self, faces):
        self._faces = faces


class Solid(Shape):
    """
    A single solid.
    """

    wrapped: TopoDS_Solid

    @classmethod
    def make_solid(cls, shell: Shell) -> "Solid":
        """
        Makes a solid from a single shell.
        """

        return cls(ShapeFix_Solid().SolidFromShell(shell.wrapped))


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

    @staticmethod
    def _make_compound(list_of_shapes: Iterable[TopoDS_Shape]) -> TopoDS_Compound:
        comp = TopoDS_Compound()
        comp_builder = TopoDS_Builder()
        comp_builder.MakeCompound(comp)

        for shape in list_of_shapes:
            comp_builder.Add(comp, shape)

        return comp

    def remove(self, shape: Shape):
        """
        Remove the specified shape.
        """

        comp_builder = TopoDS_Builder()
        comp_builder.Remove(self.wrapped, shape.wrapped)

    @classmethod
    def make_compound(cls, list_of_shapes: Iterable[Shape]) -> "Compound":
        """
        Create a compound out of a list of shapes.
        """

        return cls(cls._make_compound((s.wrapped for s in list_of_shapes)))
