from volmdlr.from_ocp import get_shells
from OCP.TopoDS import TopoDS_Solid, TopoDS_CompSolid
import volmdlr.core
from volmdlr.shells import Shell3D


class Shape(volmdlr.core.CompositePrimitive3D):
    """
    Abstract shape.
    """

    def get_shells(self):
        """
        Gets the shells of a Solid shape.
        """
        if self.wrapped:
            return get_shells(self.wrapped)
        return None


class Solid(Shape):
    """
    Describes a solid shape which is a list of shells that define de boundaries of the solid shape.
    """

    wrapped: TopoDS_Solid


class CompSolid(Shape):
    """
    Describes a solid shape which is a list of solid shape.
    """

    wrapped: TopoDS_CompSolid
