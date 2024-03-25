import unittest
import volmdlr
from volmdlr import shapes, surfaces, faces
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox


class TestShell(unittest.TestCase):
    box = BRepPrimAPI_MakeBox(2, 2, 2).Shape()
    shell = shapes.Shell(obj=box)

    def test_init1(self):
        # test init with an ocp object TopoDS_Shell
        self.assertTrue(len(self.shell._get_faces()), 8)

    def test_init2(self):
        # test init with a list of TopoDS_Face obejcts
        shell = shapes.Shell(self.shell._get_faces())
        self.assertTrue(len(shell._get_faces()), 8)

    def test_init3(self):
        # test init with a list of volmdlr.faces.Face3d objects
        faces_list = []
        for vector in [volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D]:
            for direction in [1, -1]:
                normal = direction * vector
                center = direction * vector.to_point()
                plane = surfaces.Plane3D.from_normal(center, normal)
                face = faces.PlaneFace3D.from_surface_rectangular_cut(plane, -1, 1, -1, 1)
                faces_list.append(face)
        shell = shapes.Shell(faces_list)
        self.assertTrue(len(shell._get_faces()), 8)

    def test_volume(self):
        self.assertAlmostEqual(self.shell.volume(), 8.0)

    def test_bounding_box(self):
        bbox = self.shell.bounding_box()
        self.assertAlmostEqual(bbox.volume(), 8.0, 5)
        self.assertTrue(bbox.center, volmdlr.Point3D(1.0, 1.0, 1.0))


if __name__ == '__main__':
    unittest.main()
