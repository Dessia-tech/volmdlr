import unittest
import volmdlr
import volmdlr.edges as vme


class MyTestCase(unittest.TestCase):
    def test_from_3_points(self):
        point1 = volmdlr.Point3D(607.597, 102.093550000191, 0.98)
        point2 = volmdlr.Point3D(607.597, 102.07645353018361, 0.9797543234978285)
        point3 = volmdlr.Point3D(607.597, 102.08536951838659, 0.9714579885610268)
        fullarc3d = vme.FullArc3D.from_3_points(point1, point2, point3)
        self.assertAlmostEqual(fullarc3d.radius,  0.00855)
        self.assertAlmostEqual(fullarc3d.normal.cross(volmdlr.X3D).norm(), 0)
        self.assertEqual(fullarc3d.angle, volmdlr.TWO_PI)


if __name__ == '__main__':
    unittest.main()
