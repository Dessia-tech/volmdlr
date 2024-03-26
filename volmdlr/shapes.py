"""volmdlr shapes module."""
import base64
# pylint: disable=no-name-in-module
import math
import random
import warnings
import zlib
from io import BytesIO
from itertools import chain, product
from typing import Iterable, List, Tuple, Union, Optional, Any, Dict, overload, Literal, cast as tcast

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
from OCP.BRepBuilderAPI import BRepBuilderAPI_Sewing, BRepBuilderAPI_MakeFace
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
from OCP.BRepAlgoAPI import BRepAlgoAPI_Fuse, BRepAlgoAPI_BooleanOperation, BRepAlgoAPI_Splitter
from OCP.BOPAlgo import BOPAlgo_GlueEnum, BOPAlgo_PaveFiller
from OCP.TopTools import TopTools_ListOfShape
from OCP.ShapeFix import ShapeFix_Solid
from OCP.BRepTools import BRepTools
from OCP.TopTools import TopTools_IndexedMapOfShape
from OCP.TopExp import TopExp
from OCP.gp import gp_Ax2

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


class Shape(PhysicalObject):
    """
    Represents a shape in the system. Wraps TopoDS_Shape.
    """

    wrapped: TopoDS_Shape

    def __init__(self, obj: TopoDS_Shape, name: str = ""):
        self.wrapped = downcast(obj)
        self.label = name
        self._bbox = None
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

    @staticmethod
    def _entities(obj, topo_type: Shapes) -> Iterable[TopoDS_Shape]:
        """Gets shape's entities (vertices, edges, faces, shells...)."""
        shape_set = TopTools_IndexedMapOfShape()
        TopExp.MapShapes_s(obj, inverse_shape_LUT[topo_type], shape_set)

        return tcast(Iterable[TopoDS_Shape], shape_set)

    def _get_vertices(self):
        """Gets shape's vertices, if there exists any."""
        return [downcast(i) for i in self._entities(obj=self.wrapped, topo_type="Vertex")]

    def _get_edges(self):
        """Gets shape's edges, if there exists any."""
        return [downcast(i) for i in self._entities(obj=self.wrapped, topo_type="Edge") if
                not BRep_Tool.Degenerated_s(TopoDS.Edge_s(i))]

    def _get_faces(self):
        """Gets shape's faces, if there exists any."""
        return [downcast(i) for i in self._entities(obj=self.wrapped, topo_type="Face")]

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

    def __init__(self, obj: TopoDS_Solid, name: str = '') -> None:
        super().__init__(obj=obj, name=name)

    @property
    def primitives(self) -> List[Shell]:
        """
        Gets shells from solid.
        """
        shape_set = TopTools_IndexedMapOfShape()
        TopExp.MapShapes_s(self.wrapped, topabs.TopAbs_SHELL, shape_set)
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
