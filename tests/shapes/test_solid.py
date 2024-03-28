import os
import unittest
from dessia_common.core import DessiaObject
import volmdlr
from volmdlr import faces, surfaces, shapes

folder = os.path.join(os.path.dirname(os.path.realpath(__file__)))
objects_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "shapes_objects")


class TestSolid(unittest.TestCase):
    faces_list = []
    ax = None
    for vector in [volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D]:
        for direction in [1, -1]:
            normal = direction * vector
            center = direction * vector.to_point()
            plane = surfaces.Plane3D.from_normal(center, normal)
            face = faces.PlaneFace3D.from_surface_rectangular_cut(plane, -1, 1, -1, 1)
            faces_list.append(face)

    solid1 = shapes.Solid.make_solid(shapes.Shell(faces=faces_list))
    faces2 = [f.translation(volmdlr.Vector3D(1, 1, 1)) for f in faces_list]
    solid2 = shapes.Solid.make_solid(shapes.Shell(faces=faces2))

    def test_check_platform(self):
        self.assertIsNone(self.solid1._check_platform())

    def test_to_dict_dict_to_object(self):
        to_dict = self.solid1.to_dict()
        dict_to_obejct = DessiaObject.dict_to_object(to_dict)
        self.assertEqual(dict_to_obejct, self.solid1)

    def test_to_brep_from_brep(self):
        self.solid1.to_brep(objects_folder+"/test_to_brep.brep")
        from_brep = shapes.Solid.from_brep(objects_folder+"/test_to_brep.brep")
        self.assertEqual(from_brep, self.solid1)

    def test_union(self):
        union = self.solid1.union(self.solid2)
        self.assertAlmostEqual(union.volume(), 15)

    def test_subtraction(self):
        subtraction = self.solid1.subtraction(self.solid2)
        self.assertAlmostEqual(subtraction.volume(), 7.0)

    def test_intersection(self):
        intersection = self.solid1.intersection(self.solid2)
        self.assertAlmostEqual(intersection.volume(), 1.0)

    def test_box(self):
        box = shapes.Solid.make_box(length=2, width=3, height=5)
        self.assertEqual(box.volume(), 2*3*5)

    def test_cone(self):
        cone = shapes.Solid.make_cone(radius1=0, radius2=5, height=5, direction=volmdlr.X3D, angle_degrees=270)
        self.assertAlmostEqual(cone.volume(), 98.1747704332957)

    def test_cylinder(self):
        cylinder = shapes.Solid.make_cylinder(radius=5, height=5, direction=volmdlr.X3D, angle_degrees=270)
        self.assertAlmostEqual(cylinder.volume(), 294.5243112989299)

    def test_sphere(self):
        sphere1 = shapes.Solid.make_sphere(radius=5, direction=volmdlr.X3D, angle_degrees1=-90,
                                           angle_degrees2=60, angle_degrees3=270)
        self.assertAlmostEqual(sphere1.volume(), 387.6486925590814)

    def test_torus(self):
        torus1 = shapes.Solid.make_torus(radius1=2, radius2=.5, direction=volmdlr.X3D, angle_degrees1=0,
                                         angle_degrees2=360)
        self.assertAlmostEqual(torus1.volume(), 9.869604401089358)


if __name__ == '__main__':
    unittest.main()
