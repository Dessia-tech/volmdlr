import unittest
import volmdlr
from volmdlr import shapes, wires, surfaces, faces
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

    def test_make_extrusion(self):
        length, width, height = 0.4, 0.3, 0.08
        contour = wires.Contour2D.rectangle_from_center_and_sides(volmdlr.O2D, x_length=length, y_length=width,
                                                                  is_trigo=True).to_3d(volmdlr.O3D, volmdlr.X3D,
                                                                                       -volmdlr.Y3D)
        shell = shapes.Shell.make_extrusion(contour, extrusion_direction=volmdlr.Z3D, extrusion_length=height)
        self.assertEqual(len(shell.primitives), 4)
        self.assertAlmostEqual(shell.primitives[0].area(), length * height)
        self.assertAlmostEqual(shell.primitives[1].area(), width * height)
        self.assertFalse(shell.is_closed)


if __name__ == '__main__':
    unittest.main()
