import unittest
from OCP.BRepPrimAPI import BRepPrimAPI_MakeBox
from OCP.TopoDS import TopoDS_Shell
import volmdlr
from volmdlr import curves, shapes, wires, surfaces, faces
from volmdlr.models.contours import rim_contour


class TestShell(unittest.TestCase):
    box = BRepPrimAPI_MakeBox(2, 2, 2).Shell()
    shell = shapes.Shell(obj=box)

    def test_init1(self):
        # test init with an ocp object TopoDS_Shell
        self.assertTrue(len(self.shell._get_faces()), 8)

    def test_from_faces(self):
        # test init with a list of TopoDS_Face obejcts
        shell = shapes.Shell.from_faces(faces=self.shell._get_faces())
        self.assertTrue(len(shell._get_faces()), 8)

        # test init with a list of volmdlr.faces.Face3d objects
        faces_list = []
        for vector in [volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D]:
            for direction in [1, -1]:
                normal = direction * vector
                center = direction * vector.to_point()
                plane = surfaces.Plane3D.from_normal(center, normal)
                face = faces.PlaneFace3D.from_surface_rectangular_cut(plane, -1, 1, -1, 1)
                faces_list.append(face)
        shell = shapes.Shell.from_faces(faces=faces_list)
        self.assertTrue(len(shell._get_faces()), 8)

    def test_volume(self):
        self.assertAlmostEqual(self.shell.volume(), 8.0)

    def test_bounding_box(self):
        bbox = self.shell.bounding_box()
        self.assertAlmostEqual(bbox.volume(), 8.0, 5)
        self.assertTrue(bbox.center, volmdlr.Point3D(1.0, 1.0, 1.0))

    def test_make_wedge(self):
        dx, dy, dz = 1, 2, 1
        shell = shapes.Solid.make_wedge(dx=dx, dy=dy, dz=dz, xmin=dx / 2, xmax=dx / 2, zmin=dz / 2, zmax=dz / 2,
                                        local_frame_origin=volmdlr.Point3D(-0.5, 0.5, 0.0),
                                        local_frame_direction=-volmdlr.Y3D,
                                        local_frame_x_direction=volmdlr.X3D)

        self.assertAlmostEqual(shell.volume(), (1 / 3) * dy)

        shell = shapes.Shell.make_wedge(dx=dx, dy=dy, dz=dz, xmin=dx / 4, xmax=3 * dx / 4,
                                        zmin=dz / 4, zmax=3 * dz / 4,
                                        local_frame_origin=volmdlr.Point3D(-0.5, 0.5, 0.0),
                                        local_frame_direction=-volmdlr.Y3D,
                                        local_frame_x_direction=volmdlr.X3D)

        self.assertAlmostEqual(shell.volume(), (1 / 3) * dy * (1 + 0.5 ** 2 + 0.5))

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

    def test_make_revolve(self):

        y = volmdlr.X3D.random_unit_normal_vector()
        z = volmdlr.X3D.cross(y)
        axis_point = 0.5 * volmdlr.X3D.to_point()
        frame = volmdlr.Frame3D(axis_point, volmdlr.X3D, z, y)
        wire = rim_contour.to_3d(frame.origin, frame.u, frame.v)
        revolution_shape = shapes.Shell.make_revolve(shape=wire,
                                                     axis_point=axis_point, axis=volmdlr.X3D,
                                                     angle=3.1415, name="Conical rim")
        self.assertEqual(revolution_shape.name, "Conical rim")
        self.assertEqual(len(revolution_shape.primitives), 8)
        self.assertIsInstance(revolution_shape.wrapped, TopoDS_Shell)

    def test_loft(self):
        diameter = 0.3
        circle1 = curves.Circle3D(frame=volmdlr.OXYZ, radius=diameter / 2)
        circle2 = curves.Circle3D(
            frame=volmdlr.Frame3D(volmdlr.Point3D(0.3, 0.0, 0.5), volmdlr.Y3D, volmdlr.Z3D, volmdlr.X3D),
            radius=circle1.radius / 2)
        circle3 = curves.Circle3D(
            frame=volmdlr.Frame3D(volmdlr.Point3D(0.6, 0.0, 0.3), volmdlr.Y3D, volmdlr.X3D, -volmdlr.Z3D),
            radius=circle1.radius * 0.6)
        sections = [wires.Contour3D.from_circle(circle1), wires.Contour3D.from_circle(circle2),
                    wires.Contour3D.from_circle(circle3)]
        loft = shapes.Shell.make_loft(sections=sections, ruled=True, name="loft")
        self.assertEqual(loft.name, "loft")
        self.assertEqual(len(loft.primitives), 4)
        self.assertFalse(loft.is_closed)
        self.assertIsInstance(loft.wrapped, TopoDS_Shell)


if __name__ == '__main__':
    unittest.main()
