"""
Utils to help reading and writing OCP objects.
"""
import os
from OCP.STEPCAFControl import STEPCAFControl_Reader
# from OCP.IGESCAFControl import IGESCAFControl_Reader
from OCP.Interface import Interface_Static

from OCP.TopAbs import (
    TopAbs_COMPSOLID, TopAbs_COMPOUND, TopAbs_SOLID, TopAbs_SHELL,
)
from OCP.Quantity import Quantity_Color
from OCP.TCollection import TCollection_ExtendedString
from OCP.XCAFDoc import (
    XCAFDoc_DocumentTool, XCAFDoc_ShapeTool, XCAFDoc_ColorTool,
    XCAFDoc_ColorSurf, XCAFDoc_ColorCurv, XCAFDoc_ColorGen
)
from OCP.TDF import TDF_LabelSequence, TDF_Label
from OCP.TDataStd import TDataStd_Name
from OCP.TDocStd import TDocStd_Document
from OCP.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity

from volmdlr import from_ocp
from volmdlr import curves, edges, surfaces


OCCT_TO_VOLMDLR = {"Geom_SphericalSurface": surfaces.SphericalSurface3D,
                   "Geom_CylindricalSurface": surfaces.CylindricalSurface3D,
                   "Geom_Plane": surfaces.Plane3D,
                   "Geom_ToroidalSurface": surfaces.ToroidalSurface3D,
                   "Geom_ConicalSurface": surfaces.ConicalSurface3D,
                   "Geom_BSplineSurface": surfaces.BSplineSurface3D,
                   'Geom_SurfaceOfLinearExtrusion': surfaces.ExtrusionSurface3D,
                   "Geom_SurfaceOfRevolution": surfaces.RevolutionSurface3D,
                   "Geom_Line": curves.Line3D,
                   "Geom_Circle": curves.Circle3D,
                   "Geom_Ellipse": curves.Ellipse3D,
                   "Geom_BSplineCurve": edges.BSplineCurve3D,
                   "Geom_Parabola": curves.Parabola3D,
                   "Geom_Hyperbola": curves.Hyperbola3D,
                   "Geom2d_Line": curves.Line2D,
                   "Geom2d_Circle": curves.Circle2D,
                   "Geom2d_Ellipse": curves.Ellipse2D,
                   "Geom2d_BSplineCurve": edges.BSplineCurve2D
                   }


class OCAFReader:
    """
    Generic OCAF reader for import of STEP and IGES with color and names metadata.
    """

    def load_step(self, path: str):
        """Alias for stp."""
        return self.load_stp(path)

    def load_stp(self, path: str):
        """Load a stp model."""
        if not os.path.exists(path):
            raise FileNotFoundError(f"File '{path}' does not exist!")
        reader = STEPCAFControl_Reader()
        Interface_Static.SetCVal_s("xstep.cascade.unit", "M")
        reader.SetColorMode(True)
        reader.SetNameMode(True)
        # reader.SetLayerMode(True)
        # reader.SetSHUOMode(True)
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

        :param doc:
        :return:
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
        shapes = []
        label = TDF_Label()
        shape_tool.Search(shape, label)
        entity_name = self.get_shape_name_from_label(label)
        from volmdlr.shells import Shell3D
        from volmdlr.composite_shapes import Compound
        shape_type = from_ocp.shapetype(shape)
        if shape_type == TopAbs_COMPOUND:
            return Compound.from_ocp(shape, color_tool=color_tool, shape_tool=shape_tool, name=entity_name)
        if shape_type in (TopAbs_COMPSOLID, TopAbs_SOLID):
            list_of_shells = from_ocp.get_shells(shape)
            for shell in list_of_shells:
                return Shell3D.from_ocp(shell, color_tool=color_tool, name=entity_name)
        elif shape_type == TopAbs_SHELL:
            Shell3D.from_ocp(shape, color_tool=color_tool, name=entity_name)
        return shapes

    @staticmethod
    def get_color(shape, color_tool):
        """
        Get shape color.
        """
        label = TDF_Label()
        color_tool.ShapeTool().Search(shape, label)
        if label.IsNull():
            return None
        color = Quantity_Color()
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

    @staticmethod
    def get_shape_name_from_label(label):
        """Get the entity name."""
        name_attr = TDataStd_Name()
        if label.FindAttribute(TDataStd_Name.GetID_s(), name_attr):
            entity_name = name_attr.Get().ToExtString()
            return entity_name
        return ""