import os

from OCP.TopoDS import TopoDS, TopoDS_Iterator
from OCP.TopAbs import (
    TopAbs_COMPSOLID, TopAbs_COMPOUND, TopAbs_SOLID, TopAbs_SHELL,
    TopAbs_FACE
)
from OCP.TopExp import TopExp
from OCP.TopTools import TopTools_IndexedMapOfShape

from OCP.Quantity import Quantity_Color, Quantity_ColorRGBA
from OCP.TCollection import TCollection_ExtendedString
from OCP.XCAFDoc import (
    XCAFDoc_DocumentTool,
    XCAFDoc_ColorSurf, XCAFDoc_ColorCurv, XCAFDoc_ColorGen
)
from OCP.TDF import TDF_LabelSequence, TDF_Label
from OCP.TDataStd import TDataStd_Name
from OCP.TDocStd import TDocStd_Document

from OCP.IFSelect import IFSelect_RetDone
from OCP.STEPCAFControl import STEPCAFControl_Reader

from volmdlr.shapes import Shell, Solid, CompSolid, Compound


class OCAFReader:
    """
    Generic OCAF reader for import of STEP and IGES with color and names metadata.
    """

    def load_step(self, path: str):
        """
        Alias for stp.

        :param path: path to file.
        :return:
        """
        return self.load_stp(path)

    def load_stp(self, path: str):
        """
        Load a STEP model from the specified path.

        :param path: The path to the STEP file.
        :type path: str
        :return: A list of model shapes with associated metadata.
        """
        if not os.path.exists(path):
            raise FileNotFoundError(f"File '{path}' does not exist!")
        reader = STEPCAFControl_Reader()
        reader.SetColorMode(True)
        reader.SetNameMode(True)
        reader.SetLayerMode(True)
        reader.SetSHUOMode(True)
        status = reader.ReadFile(path)
        if status != IFSelect_RetDone:
            raise RuntimeError("Cannot open STEP file")

        if not reader.NbRootsForTransfer():
            raise ValueError(f"File has no shapes: {path}")

        name = TCollection_ExtendedString(f"STEP-{id(self)}")
        doc = TDocStd_Document(name)
        if not reader.Transfer(doc):
            raise ValueError(f"Failed to transfer: {path}")

        return self._process_doc(doc)

    def _process_doc(self, doc):
        """
        Process the OCAF document, extracting its shapes and associated metadata.

        :param doc: The OCAF document.
        :return: A list of dictionaries containing shape data.
        """
        main = doc.Main()
        color_tool = XCAFDoc_DocumentTool.ColorTool_s(main)
        shape_tool = XCAFDoc_DocumentTool.ShapeTool_s(main)
        labels = TDF_LabelSequence()
        shape_tool.GetFreeShapes(labels)

        shapes = []
        for i in range(1, labels.Length() + 1):
            label = labels.Value(i)
            shape = shape_tool.GetShape_s(label)
            shapes.append(self._process_shape(shape, color_tool, shape_tool))
        return shapes

    def _process_shape(self, shape, color_tool, shape_tool=None):
        """
        Process an individual shape, creating a wrapper object and collecting metadata.

        :param shape: The shape to process.
        :param color_tool: The color tool.
        :param shape_tool: The shape tool (optional).
        :return: A dictionary containing shape data.
        """
        shape_type = shape.ShapeType()
        name = self._get_shape_name(shape=shape, shape_tool=shape_tool)
        children_names_dict, map_of_colors = self._process_metadata(shape=shape,
                                                                    color_tool=color_tool, shape_tool=shape_tool)

        if shape_type == TopAbs_COMPOUND:
            wrapping_shape = Compound(obj=shape, name=name)
        elif shape_type == TopAbs_COMPSOLID:
            wrapping_shape = CompSolid(obj=shape, name=name)
        elif shape_type == TopAbs_SOLID:
            wrapping_shape = Solid(obj=shape, name=name)
        elif shape_type == TopAbs_SHELL:
            wrapping_shape = Shell(obj=shape, name=name)
        else:
            raise NotImplementedError(f"Convertion of {shape_type} not implemented yet.")
        data_structure = {
            "shape": wrapping_shape,
            "children_names": children_names_dict,
            "shell_colors": map_of_colors,
        }
        return data_structure

    def _process_metadata(self, shape, color_tool, shape_tool):
        """
        Process metadata for a given shape, including colors and names.

        :param shape: The shape to process.
        :param color_tool: The color tool.
        :param shape_tool: The shape tool.
        :return: Maps of names and colors.
        """
        shape_type = shape.ShapeType()
        map_of_colors = {}
        map_of_names = {}
        if shape_type in (TopAbs_COMPOUND, TopAbs_COMPSOLID):
            it = TopoDS_Iterator()
            it.Initialize(shape, False, False)
            map_of_names[shape] = self._get_shape_name(shape, shape_tool)
            while it.More():
                names, colors = self._process_metadata(shape=it.Value(), color_tool=color_tool, shape_tool=shape_tool)
                map_of_names.update(names)
                map_of_colors.update(colors)
                it.Next()
        elif shape_type == TopAbs_SOLID:
            map_of_colors.update(self._process_solid(shape=shape, color_tool=color_tool))
            map_of_names[shape] = self._get_shape_name(shape, shape_tool)
        elif shape_type == TopAbs_SHELL:
            map_of_colors[shape] = self._process_shell(shape=shape, color_tool=color_tool)
            map_of_names[shape] = self._get_shape_name(shape, shape_tool)
        else:
            raise NotImplementedError
        return map_of_names, map_of_colors

    def _process_solid(self, shape, color_tool):
        """
        Process metadata for a solid shape.

        :param shape: The solid shape to process.
        :param color_tool: The color tool.
        :return: A map of colors.
        """
        map_of_colors = {}
        color = self._get_color(shape, color_tool)
        it = TopoDS_Iterator()
        it.Initialize(shape, False, False)
        while it.More():
            map_of_colors[it.Value()] = self._process_shell(it.Value(), color_tool, color)
            it.Next()
        return map_of_colors

    def _process_shell(self, shape, color_tool, default_color: tuple = (0.8, 0.8, 0.8)):
        """
        Process metadata for a shell shape.

        :param shape: The shell shape to process.
        :param color_tool: The color tool.
        :param default_color: Default color if face colors are missing.
        :return: A list of face colors.
        """
        shell_color = self._get_color(shape, color_tool) or default_color
        faces_color = []
        faces_set = TopTools_IndexedMapOfShape()
        TopExp.MapShapes_s(shape, TopAbs_FACE, faces_set)
        for face in faces_set:
            face_color = self._get_color(face, color_tool)
            if face_color:
                faces_color.append(face_color)
            else:
                faces_color.append(shell_color)
        return faces_color

    def _get_color(self, shape, color_tool):
        """
        Retrieves the color of a shape.

        :param shape: The shape to get the color for.
        :param color_tool: The color tool to use.
        :return: The color of the shape, or None if no color is found.
        """
        label = TDF_Label()
        color_tool.ShapeTool().Search(shape, label)
        return self._get_color_from_label(label=label, color_tool=color_tool)

    @staticmethod
    def _get_color_from_label(label, color_tool):
        """
        Retrieves the color from a label.

        :param label: The label to get the color from.
        :param color_tool: The color tool to use.
        :return: The color from the label, or None if no color is found.
        """
        if label.IsNull():
            return None
        color = Quantity_Color()
        # color = Quantity_ColorRGBA()
        while True:
            if color_tool.GetColor_s(label, color):
                return color.Red(), color.Green(), color.Blue()
            if color_tool.GetColor_s(label, XCAFDoc_ColorSurf, color):
                return color.Red(), color.Green(), color.Blue()
            if color_tool.GetColor_s(label, XCAFDoc_ColorCurv, color):
                return color.Red(), color.Green(), color.Blue()
            label = label.Father()
            if label.IsNull():
                break
        return None

    def _get_shape_name(self, shape, shape_tool):
        """

        """
        label = TDF_Label()
        shape_tool.Search(shape, label)
        return self._get_shape_name_from_label(label)

    @staticmethod
    def _get_shape_name_from_label(label):
        """

        :param label:
        :return:
        """
        name_attr = TDataStd_Name()
        if label.FindAttribute(TDataStd_Name.GetID_s(), name_attr):
            entity_name = name_attr.Get().ToExtString()
            return entity_name
        return ""
