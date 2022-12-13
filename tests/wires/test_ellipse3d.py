import unittest
import volmdlr
from volmdlr import wires
from volmdlr.models import ellipse2d as model_ellipse2d


class TestEllipse3D(unittest.TestCase):
    ellipse = wires.Ellipse3D(4, 3, volmdlr.O3D, volmdlr.Z3D, volmdlr.X3D)

    def test_point_belongs(self):
        point_on_ellipse = volmdlr.Point3D(2.8284271247461903, 2.1213203435596424, 0)
        point_on_ellipse2 = volmdlr.Point3D(4, 0, 0)
        point_not_on_ellipse = volmdlr.Point3D(3, 3, 0)
        self.assertTrue(self.ellipse.point_belongs(point_on_ellipse))
        self.assertTrue(self.ellipse.point_belongs(point_on_ellipse2))
        self.assertFalse(self.ellipse.point_belongs(point_not_on_ellipse))

    def test_length(self):
        self.assertAlmostEqual(self.ellipse.length(), 22.1034921607)

    def test_discretization_points(self):
        discretization_points = self.ellipse.discretization_points(number_points=4)
        expected_points = [volmdlr.Point3D(4.0, 0.0, 0.0),
                           volmdlr.Point3D(0, -3, 0),
                           volmdlr.Point3D(-4, 0, 0),
                           volmdlr.Point3D(0, 3, 0)]
        for expected_point, point in zip(expected_points, discretization_points):
            self.assertEqual(expected_point, point)

    def test_to_2d(self):
        expected_ellipse2d = model_ellipse2d.ellipse_2d
        vector_2 = self.ellipse.normal.cross(self.ellipse.major_dir)
        ellipse_2d = self.ellipse.to_2d(self.ellipse.center, self.ellipse.major_dir, vector_2)
        self.assertEqual(expected_ellipse2d.center, ellipse_2d.center)
        self.assertEqual(expected_ellipse2d.major_dir, ellipse_2d.major_dir)

    def test_abscissa(self):
        point_on_ellipse = volmdlr.Point3D(2.8284271247461903, 2.1213203435596424, 0)
        point_on_ellipse2 = volmdlr.Point3D(4, 0, 0)
        point_not_on_ellipse = volmdlr.Point3D(3, 3, 0)
        self.assertAlmostEqual(self.ellipse.abscissa(point_on_ellipse), 2.022137530)
        self.assertAlmostEqual(self.ellipse.abscissa(point_on_ellipse2), 22.1034921607)
        with self.assertRaises(ValueError):
            self.ellipse.abscissa(point_not_on_ellipse)


if __name__ == '__main__':
    unittest.main()
