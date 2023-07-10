import math
import unittest

from dessia_common.core import DessiaObject
import volmdlr
from volmdlr import edges, faces, primitives3d, wires, surfaces, shells


class TestShell3D(unittest.TestCase):
    def test_from_faces(self):
        openshell = DessiaObject.load_from_file('shells/open_shell.json')
        closedshell = DessiaObject.load_from_file('shells/closed_shell.json')
        fm_closedshell = closedshell.frame_mapping(volmdlr.Frame3D(volmdlr.Point3D(0.5, 0.5, 0.5),
                                                                   volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D), 'new')
        from_faces = shells.Shell3D.from_faces(openshell.faces + fm_closedshell.faces)
        self.assertEqual(from_faces[0], openshell)
        self.assertEqual(from_faces[1], fm_closedshell)


if __name__ == '__main__':
    unittest.main()