import os
import unittest

import volmdlr
from volmdlr import shapes


class TestSolidBooleanOperations(unittest.TestCase):
    def test_sequential_boolean_operations(self):
        box = shapes.Solid.make_box(length=1, width=1, height=1, point=volmdlr.Point3D(-.5, -.5, -.5))
        sphere1 = shapes.Solid.make_sphere(radius=0.68, angle_degrees1=-90, angle_degrees2=90, angle_degrees3=360)

        cyl1 = shapes.Solid.make_cylinder(radius=0.3, height=2, point=volmdlr.Point3D(0, 0, -1))
        cyl2 = shapes.Solid.make_cylinder(radius=0.3, height=2, point=volmdlr.Point3D(-1, 0, 0), direction=volmdlr.X3D)
        cyl3 = shapes.Solid.make_cylinder(radius=0.3, height=2, point=volmdlr.Point3D(0, 1, 0), direction=-volmdlr.Y3D)

        box_intersection_sphere = box.intersection(sphere1)
        self.assertAlmostEqual(box_intersection_sphere.volume(), 0.9384398022558381)

        cylinder1_cylinder2_union = cyl1.union(cyl2)
        cyl1_cyl2_cyl3_union = cylinder1_cylinder2_union.union(cyl3)
        self.assertAlmostEqual(cyl1_cyl2_cyl3_union.volume(), 1.3909899141249489)

        substraction = box_intersection_sphere.subtraction(cyl1_cyl2_cyl3_union)
        self.assertAlmostEqual(substraction.volume(), 0.3956799040690441)


if __name__ == '__main__':
    unittest.main()
