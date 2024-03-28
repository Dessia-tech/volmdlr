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
from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeFace, BRepBuilderAPI_Copy
from OCP.BRepPrimAPI import BRepPrimAPI_MakePrism, BRepPrimAPI_MakeWedge
from OCP.TopTools import TopTools_IndexedMapOfShape
from OCP.TopExp import TopExp
from OCP.Geom import Geom_Plane
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
from OCP.BRepPrimAPI import (
    BRepPrimAPI_MakeBox,
    BRepPrimAPI_MakeCone,
    BRepPrimAPI_MakeCylinder,
    BRepPrimAPI_MakeTorus,
    BRepPrimAPI_MakeWedge,
    BRepPrimAPI_MakePrism,
    BRepPrimAPI_MakeRevol,
    BRepPrimAPI_MakeSphere,
)
from OCP.TopLoc import TopLoc_Location

from OCP.gp import gp_Ax2


import volmdlr.core_compiled
from volmdlr import curves, display, edges, surfaces, wires, geometry, faces as vm_faces
from volmdlr.core import edge_in_list, get_edge_index_in_list, get_point_index_in_list
from volmdlr.utils.step_writer import geometric_context_writer, product_writer, step_ids_to_str
from volmdlr import from_ocp, to_ocp
from volmdlr.utils.ocp_helpers import plot_edge
from volmdlr.utils.mesh_helpers import perform_decimation

import OCP.TopAbs as top_abs  # Topology type enum
from OCP.TopAbs import TopAbs_Orientation

shape_LUT = {
    top_abs.TopAbs_VERTEX: "Vertex",
    top_abs.TopAbs_EDGE: "Edge",
    top_abs.TopAbs_WIRE: "Wire",
    top_abs.TopAbs_FACE: "Face",
    top_abs.TopAbs_SHELL: "Shell",
    top_abs.TopAbs_SOLID: "Solid",
    top_abs.TopAbs_COMPSOLID: "CompSolid",
    top_abs.TopAbs_COMPOUND: "Compound",
}

inverse_shape_LUT = {v: k for k, v in shape_LUT.items()}

Shapes = Literal[
    "Vertex", "Edge", "Wire", "Face", "Shell", "Solid", "CompSolid", "Compound"]


# pylint: disable=no-name-in-module,invalid-name,unused-import,wrong-import-order


def shapetype(obj: TopoDS_Shape) -> top_abs.TopAbs_ShapeEnum:
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
        top_abs.TopAbs_VERTEX: TopoDS.Vertex_s,
        top_abs.TopAbs_EDGE: TopoDS.Edge_s,
        top_abs.TopAbs_WIRE: TopoDS.Wire_s,
        top_abs.TopAbs_FACE: TopoDS.Face_s,
        top_abs.TopAbs_SHELL: TopoDS.Shell_s,
        top_abs.TopAbs_SOLID: TopoDS.Solid_s,
        top_abs.TopAbs_COMPSOLID: TopoDS.CompSolid_s,
        top_abs.TopAbs_COMPOUND: TopoDS.Compound_s,
    }

    f_downcast: Any = downcast_LUT[shapetype(obj)]
    downcasted_obj = f_downcast(obj)

    return downcasted_obj


def _make_wedge(
        dx: float,
        dy: float,
        dz: float,
        xmin: float,
        zmin: float,
        xmax: float,
        zmax: float,
        point: volmdlr.Vector3D = volmdlr.O3D,
        direction: volmdlr.Vector3D = volmdlr.Z3D,
        x_direction=volmdlr.X3D) -> BRepPrimAPI_MakeWedge:
    """
    Make a wedge builder.

    This is a private method and should not be used directlly. Please see Solid.make_wedge
    or Shell.make_wedge for details.
    """

    return BRepPrimAPI_MakeWedge(
        gp_Ax2(to_ocp.point3d_to_ocp(point), to_ocp.vector3d_to_ocp(direction, unit_vector=True),
               to_ocp.vector3d_to_ocp(x_direction, unit_vector=True)),
        dx,
        dy,
        dz,
        xmin,
        zmin,
        xmax,
        zmax,
    )


class Shape(volmdlr.core.Primitive3D):
    """
    Represents a shape in the system. Wraps TopoDS_Shape.
    """
    _non_serializable_attributes = ["obj"]
    _non_data_eq_attributes = ['wrapped', 'name', 'color', 'alpha']
    wrapped: TopoDS_Shape

    def __init__(self, obj: TopoDS_Shape, name: str = ""):
        self.wrapped = downcast(obj)
        self.label = name
        self._bbox = None
        super().__init__(name=name)

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

        # define the shape lookup table for casting
        constructor_LUT = {
            top_abs.TopAbs_SHELL: Shell,
            top_abs.TopAbs_SOLID: Solid,
            top_abs.TopAbs_COMPSOLID: CompSolid,
            top_abs.TopAbs_COMPOUND: Compound,
        }

        shape_type = shapetype(obj=obj)
        # NB downcast is needed to handle TopoDS_Shape types
        return constructor_LUT[shape_type](obj=downcast(obj), name=name)

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

    def _get_shells(self) -> List["Shell"]:
        """
        :returns: All the shells in this Shape.
        """

        return [Shell(obj=i) for i in self._entities(obj=self.wrapped, topo_type="Shell")]

    def _get_solids(self) -> List["Solid"]:
        """
        :returns: All the solids in this Shape.
        """

        return [Solid(obj=i) for i in self._entities(obj=self.wrapped, topo_type="Solid")]

    def _get_compsolids(self) -> List["CompSolid"]:
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
            return getattr(shape, f"_get_{cls.__name__.lower()}s")()[0]
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
        if isinstance(self, Shell) and not self.is_closed:
            raise ValueError("The shell is an open shell and the volume can't be calculated."
                             " Try using a closed shell.")
        BRepGProp.VolumeProperties_s(self.wrapped, prop, tol)
        return abs(prop.Mass())

    def _bool_op(
            self,
            args: Iterable["Shape"],
            tools: Iterable["Shape"],
            operation: Union[BRepAlgoAPI_BooleanOperation, BRepAlgoAPI_Splitter],
            parallel: bool = True,
    ) -> "TopoDS_Shell":
        """
        Generic boolean operation.

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
        shape = self.__class__(operation.Shape())
        return getattr(shape, f"_get_{self.__class__.__name__.lower()}s")()

    def subtraction(self, *to_subtract: "Shape", tol: Optional[float] = None) -> "Shape":
        """
        Subtract the positional arguments from this Shape.

        :param tol: Fuzzy mode tolerance
        """

        cut_op = BRepAlgoAPI_Cut()

        if tol:
            cut_op.SetFuzzyValue(tol)

        return self._bool_op((self,), to_subtract, cut_op)

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

        return self._bool_op((self,), to_union, fuse_op)

    def intersection(self, *to_intersect: "Shape", tol: Optional[float] = None) -> "Shape":
        """
        Intersection of the positional arguments and this Shape.

        :param tol: Fuzzy mode tolerance
        """

        intersect_op = BRepAlgoAPI_Common()

        if tol:
            intersect_op.SetFuzzyValue(tol)

        return self._bool_op((self,), to_intersect, intersect_op)

    def plot(self, ax=None, edge_style=volmdlr.core.EdgeStyle()):
        """Plots a shape using matplolib."""
        shape_edges = self._get_edges()
        for edge in shape_edges:
            ax = plot_edge(edge, ax, edge_style)
        return ax

    def volmdlr_primitives(self):
        """Gets shape's volmdlr primitives."""
        return [self]

    def mesh(self, tolerance: float, angular_tolerance: float = 0.1):
        """
        Generate triangulation if none exists.
        """
        if not BRepTools.Triangulation_s(self.wrapped, tolerance):
            BRepMesh_IncrementalMesh(self.wrapped, tolerance, True, angular_tolerance)

    def tessellate(self, tolerance: float, angular_tolerance: float = 0.1) ->\
            [NDArray[float],  List[Tuple[int, int, int]]]:
        """
        Tessellates the geometry into vertices and triangles.

        :param tolerance: (float) - Tolerance value for meshing.
        :param angular_tolerance: (float, optional) - Angular tolerance value. Defaults to 0.1.

        :return: Tuple[List[List[float]], List[Tuple[int, int, int]]]: A tuple containing:
                - A list of vertices represented as lists of floats,
                where each sublist represents the (x, y, z) coordinates of a vertex.
                - A list of triangles represented as tuples of integers,
                where each tuple contains the indices of the vertices forming a triangle.
        """

        self.mesh(tolerance, angular_tolerance)

        vertices = []
        triangles = []
        offset = 0

        for face in self._get_faces():
            loc = TopLoc_Location()
            poly = BRep_Tool.Triangulation_s(face, loc)
            trsf = loc.Transformation()
            reverse = face.Orientation() == TopAbs_Orientation.TopAbs_REVERSED

            # add vertices
            vertices += [
                [v.X(), v.Y(), v.Z()] for v in (poly.Node(i).Transformed(trsf) for i in range(1, poly.NbNodes() + 1))
            ]

            # add triangles
            triangles += [
                (triangle.Value(1) + offset - 1, triangle.Value(3) + offset - 1, triangle.Value(2) + offset - 1)
                if reverse else
                (triangle.Value(1) + offset - 1, triangle.Value(2) + offset - 1, triangle.Value(3) + offset - 1)
                for triangle in poly.Triangles()
            ]
            offset += poly.NbNodes()

        return vertices, triangles

    def triangulation(self):
        """
        Gets shape triagulation.

        """
        vertices, triangles = self.tessellate(tolerance=1e-2)
        mesh = display.Mesh3D(np.array(vertices), np.array(triangles))
        return mesh

    def babylon_meshes(self, *args, **kwargs):
        """
        Returns the babylonjs mesh.

        """
        mesh = self.triangulation()
        if mesh is None:
            return []
        babylon_mesh = mesh.to_babylon()
        babylon_mesh.update({
            'alpha': self.alpha,
            'name': self.name,
            'color': list(self.color) if self.color is not None else [0.8, 0.8, 0.8]
        })
        babylon_mesh["reference_path"] = self.reference_path
        return [babylon_mesh]


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
    def is_closed(self):
        """
        Returns True if shell is a closed shell and False otherwise.
        """
        return self.wrapped.Closed()

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

    @property
    def primitives(self) -> List[vm_faces.Face3D]:
        """
        Gets shell's faces.
        """
        return [vm_faces.Face3D.from_ocp(downcast(shape)) for shape in self._get_faces()]

    @classmethod
    def make_wedge(cls,
                   dx: float,
                   dy: float,
                   dz: float,
                   xmin: float,
                   zmin: float,
                   xmax: float,
                   zmax: float,
                   local_frame_origin: volmdlr.Point3D = volmdlr.O3D,
                   local_frame_direction: volmdlr.Vector3D = volmdlr.Z3D,
                   local_frame_x_direction: volmdlr.Vector3D = volmdlr.X3D,
                   ) -> "Shell":
        """
        Creates a wedge, which can represent a pyramid or a truncated pyramid.

        The origin of the local coordinate system is the corner of the base rectangle of the wedge.
        The y-axis represents the "height" of the pyramid or truncated pyramid.

        To create a pyramid, specify xmin=xmax=dx/2 and zmin=zmax=dz/2.

        :param dx: The length of the base rectangle along the x-axis.
        :type dx: float
        :param dy: The height of the pyramid or truncated pyramid along the y-axis.
        :type dy: float
        :param dz: The width of the base rectangle along the z-axis.
        :type dz: float
        :param xmin: The x-coordinate of one corner of the top rectangle.
        :type xmin: float
        :param zmin: The z-coordinate of one corner of the top rectangle.
        :type zmin: float
        :param xmax: The x-coordinate of the opposite corner of the top rectangle.
        :type xmax: float
        :param zmax: The z-coordinate of the opposite corner of the top rectangle.
        :type zmax: float
        :param local_frame_origin: The origin of the local coordinate system for the wedge.
         Defaults to the origin (0, 0, 0).
        :type local_frame_origin: volmdlr.Point3D
        :param local_frame_direction: The main direction for the local coordinate system of the wedge.
         Defaults to the z-axis (0, 0, 1).
        :type local_frame_direction: volmdlr.Vector3D
        :param local_frame_x_direction: The x direction for the local coordinate system of the wedge.
         Defaults to the x-axis (1, 0, 0).
        :type local_frame_x_direction: volmdlr.Vector3D

        :return: The created wedge.
        :rtype: Shell

        Example:
        To create a pyramid with a square base of size 1 and where its apex is located at
        volmdlr.Point3D(0.0, 0.0, 2.0):
        >>> dx, dy, dz = 1, 2, 1
        >>> wedge = Shell.make_wedge(dx=dx, dy=dy, dz=dz, xmin=dx / 2, xmax=dx / 2, zmin=dz / 2, zmax=dz / 2,
        >>>                                 local_frame_origin=volmdlr.Point3D(-0.5, 0.5, 0.0),
        >>>                                 local_frame_direction=-volmdlr.Y3D,
        >>>                                 local_frame_x_direction=volmdlr.X3D)

        """

        return cls(obj=_make_wedge(dx=dx, dy=dy, dz=dz, xmin=xmin, zmin=zmin, xmax=xmax, zmax=zmax,
                                   point=local_frame_origin,
                                   direction=local_frame_direction,
                                   x_direction=local_frame_x_direction).Shell())

    @classmethod
    def make_extrusion(cls, wire: wires.Contour3D, extrusion_direction: volmdlr.Vector3D,
                       extrusion_length: float, name: str = '') -> "Shell":
        """
        Returns a solid generated by the extrusion of a plane face.
        """
        ocp_wire = to_ocp.contour3d_to_ocp(wire)
        extrusion_vector = to_ocp.vector3d_to_ocp(extrusion_direction * extrusion_length)

        return cls(obj=BRepPrimAPI_MakePrism(ocp_wire, extrusion_vector).Shape(), name=name)


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

    @property
    def primitives(self) -> List[Shell]:
        """
        Gets shells from solid.
        """
        shape_set = TopTools_IndexedMapOfShape()
        TopExp.MapShapes_s(self.wrapped, top_abs.TopAbs_SHELL, shape_set)
        return [Shell(obj=shape) for shape in shape_set]

    @classmethod
    def make_wedge(cls,
                   dx: float,
                   dy: float,
                   dz: float,
                   xmin: float,
                   zmin: float,
                   xmax: float,
                   zmax: float,
                   local_frame_origin: volmdlr.Point3D = volmdlr.O3D,
                   local_frame_direction: volmdlr.Vector3D = volmdlr.Z3D,
                   local_frame_x_direction: volmdlr.Vector3D = volmdlr.X3D,
                   ) -> "Solid":
        """
        Creates a wedge, which can represent a pyramid or a truncated pyramid.

        The origin of the local coordinate system is the corner of the base rectangle of the wedge.
        The y-axis represents the "height" of the pyramid or truncated pyramid.

        To create a pyramid, specify xmin=xmax=dx/2 and zmin=zmax=dz/2.

        :param dx: The length of the base rectangle along the x-axis.
        :type dx: float
        :param dy: The height of the pyramid or truncated pyramid along the y-axis.
        :type dy: float
        :param dz: The width of the base rectangle along the z-axis.
        :type dz: float
        :param xmin: The x-coordinate of one corner of the top rectangle.
        :type xmin: float
        :param zmin: The z-coordinate of one corner of the top rectangle.
        :type zmin: float
        :param xmax: The x-coordinate of the opposite corner of the top rectangle.
        :type xmax: float
        :param zmax: The z-coordinate of the opposite corner of the top rectangle.
        :type zmax: float
        :param local_frame_origin: The origin of the local coordinate system for the wedge.
         Defaults to the origin (0, 0, 0).
        :type local_frame_origin: volmdlr.Point3D
        :param local_frame_direction: The main direction for the local coordinate system of the wedge.
         Defaults to the z-axis (0, 0, 1).
        :type local_frame_direction: volmdlr.Vector3D
        :param local_frame_x_direction: The x direction for the local coordinate system of the wedge.
         Defaults to the x-axis (1, 0, 0).
        :type local_frame_x_direction: volmdlr.Vector3D

        :return: The created wedge.
        :rtype: Solid

        Example:
        To create a pyramid with a square base of size 1 and where its apex is located at
        volmdlr.Point3D(0.0, 0.0, 2.0):
        >>> dx, dy, dz = 1, 2, 1
        >>> wedge = Solid.make_wedge(dx=dx, dy=dy, dz=dz, xmin=dx / 2, xmax=dx / 2, zmin=dz / 2, zmax=dz / 2,
        >>>                                 local_frame_origin=volmdlr.Point3D(-0.5, 0.5, 0.0),
        >>>                                 local_frame_direction=-volmdlr.Y3D,
        >>>                                 local_frame_x_direction=volmdlr.X3D)

        """

        return cls(obj=_make_wedge(dx=dx, dy=dy, dz=dz, xmin=xmin, zmin=zmin, xmax=xmax, zmax=zmax,
                                   point=local_frame_origin,
                                   direction=local_frame_direction,
                                   x_direction=local_frame_x_direction).Solid())

    @classmethod
    def make_extrusion(cls, face: Union[vm_faces.PlaneFace3D, Geom_Plane],
                       extrusion_length: float, name: str = '') -> "Solid":
        """
        Returns a solid generated by the extrusion of a plane face.
        """
        ocp_face = face
        if isinstance(ocp_face, vm_faces.Face3D):
            ocp_face = face.to_ocp()
        extrusion_vector = to_ocp.vector3d_to_ocp(face.surface3d.frame.w * extrusion_length)
        solid = BRepPrimAPI_MakePrism(ocp_face, extrusion_vector)
        return cls(obj=solid.Shape(), name=name)

    @classmethod
    def make_extrusion_from_frame_and_wires(cls, frame: volmdlr.Frame3D,
                                            outer_contour2d: volmdlr.wires.Contour2D,
                                            inner_contours2d: List[volmdlr.wires.Contour2D],
                                            extrusion_length: float, name: str = '') -> "Solid":
        """
        Returns a solid generated by the extrusion of a plane face.
        """
        face = vm_faces.PlaneFace3D(surface3d=surfaces.Plane3D(frame),
                                    surface2d=surfaces.Surface2D(outer_contour2d, inner_contours2d))

        solid = Solid.make_extrusion(face=face, extrusion_length=extrusion_length)

        return cls(obj=solid.wrapped, name=name)

    @classmethod
    def make_box(
            cls,
            length: float,
            width: float,
            height: float,
            point: volmdlr.Point3D = volmdlr.Point3D(0, 0, 0),
            direction: volmdlr.Vector3D = volmdlr.Vector3D(0, 0, 1),
    ) -> "Solid":
        """
        Make a box located in point with the dimensions (length,width,height).

        By default, pnt=volmdlr.Point3D(0,0,0) and dir=volmdlr.Vector3D(0,0,1).
        """
        frame = volmdlr.Frame3D.from_point_and_normal(point, direction)
        return cls(
            BRepPrimAPI_MakeBox(
                to_ocp.frame3d_to_ocp(frame=frame, right_handed=True), length, width, height
            ).Shape()
        )

    @classmethod
    def make_cone(
            cls,
            radius1: float,
            radius2: float,
            height: float,
            point: volmdlr.Point3D = volmdlr.Point3D(0, 0, 0),
            direction: volmdlr.Vector3D = volmdlr.Vector3D(0, 0, 1),
            angle_degrees: float = 360,
    ) -> "Solid":
        """
        Make a cone with given radii and height.

        By default, pnt=volmdlr.Point3D(0,0,0), dir=volmdlr.Vector3D(0,0,1) and angle=360.
        """
        frame = volmdlr.Frame3D.from_point_and_normal(point, direction)
        return cls(
            BRepPrimAPI_MakeCone(
                to_ocp.frame3d_to_ocp(frame=frame, right_handed=True),
                radius1,
                radius2,
                height,
                math.radians(angle_degrees),
            ).Shape()
        )

    @classmethod
    def make_cylinder(
            cls,
            radius: float,
            height: float,
            point: volmdlr.Point3D = volmdlr.Point3D(0, 0, 0),
            direction: volmdlr.Vector3D = volmdlr.Vector3D(0, 0, 1),
            angle_degrees: float = 360,
    ) -> "Solid":
        """
        Make a cylinder with a given radius and height.

        By default point=volmdlr.Point3D(0,0,0),dir=voldmlr.Vector3D(0,0,1) and angle=360.

        """
        frame = volmdlr.Frame3D.from_point_and_normal(point, direction)
        return cls(
            BRepPrimAPI_MakeCylinder(to_ocp.frame3d_to_ocp(frame=frame, right_handed=True),
                                     radius, height, math.radians(angle_degrees), ).Shape()
        )

    @classmethod
    def make_torus(
            cls,
            radius1: float,
            radius2: float,
            point: volmdlr.Point3D = volmdlr.Point3D(0, 0, 0),
            direction: volmdlr.Vector3D = volmdlr.Vector3D(0, 0, 1),
            angle_degrees1: float = 0,
            angle_degrees2: float = 360,
    ) -> "Solid":
        """
        Make a torus with a given radii and angles.

        By default, point=volmdlr.Point3D(0,0,0),direction=volmdlr.Vector3D(0,0,1),angle1=0
        ,angle1=360 and angle=360.
        """
        frame = volmdlr.Frame3D.from_point_and_normal(point, direction)
        return cls(
            BRepPrimAPI_MakeTorus(
                to_ocp.frame3d_to_ocp(frame=frame, right_handed=True),
                radius1,
                radius2,
                math.radians(angle_degrees1),
                math.radians(angle_degrees2),
            ).Shape()
        )

    @classmethod
    def make_sphere(
            cls,
            radius: float,
            point: volmdlr.Point3D = volmdlr.Point3D(0, 0, 0),
            direction: volmdlr.Vector3D = volmdlr.Vector3D(0, 0, 1),
            angle_degrees1: float = 0,
            angle_degrees2: float = 90,
            angle_degrees3: float = 360,
    ) -> "Shape":
        """
        Make a sphere with a given radius.

        By default, point=volmdlr.Point3D(0,0,0),direction=volmdlr.Vector3D(0,0,1), angle1=0, angle2=90 and angle3=360
        """
        frame = volmdlr.Frame3D.from_point_and_normal(point, direction)
        return cls(
            BRepPrimAPI_MakeSphere(
                to_ocp.frame3d_to_ocp(frame=frame, right_handed=True),
                radius,
                math.radians(angle_degrees1),
                math.radians(angle_degrees2),
                math.radians(angle_degrees3),
            ).Shape()
        )


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
