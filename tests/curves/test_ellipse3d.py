import math
import unittest

import volmdlr
from volmdlr import curves, edges


class TestEllipse3D(unittest.TestCase):
    ellipse = curves.Ellipse3D(4, 3, volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D,
                                                     volmdlr.Z3D.cross(volmdlr.X3D), volmdlr.Z3D))
    vector1 = volmdlr.Vector3D(1, 1, 1)
    vector1 = vector1.unit_vector()
    vector2 = vector1.deterministic_unit_normal_vector()
    vector3 = vector1.cross(vector2)
    frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)
    ellipse3d = curves.Ellipse3D(2, 1, frame)

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
            self.assertTrue(expected_point.is_close(point))

    def test_to_2d(self):
        vector_2 = self.ellipse.normal.cross(self.ellipse.major_dir)
        ellipse_2d = self.ellipse.to_2d(self.ellipse.center, self.ellipse.major_dir, vector_2)
        self.assertEqual(volmdlr.O2D, ellipse_2d.center)
        self.assertEqual(volmdlr.X2D, ellipse_2d.major_dir)

    def test_abscissa(self):
        point_on_ellipse = volmdlr.Point3D(2.8284271247461903, 2.1213203435596424, 0)
        point_on_ellipse2 = volmdlr.Point3D(4, 0, 0)
        point_not_on_ellipse = volmdlr.Point3D(3, 3, 0)
        self.assertAlmostEqual(self.ellipse.abscissa(point_on_ellipse), 2.513786016470093)
        self.assertAlmostEqual(self.ellipse.abscissa(point_on_ellipse2), 0.0)
        with self.assertRaises(ValueError):
            self.ellipse.abscissa(point_not_on_ellipse)
        point3d = volmdlr.Point3D(-0.236335934849, -1.104105421298, -1.104105421298)
        self.assertAlmostEqual(self.ellipse3d.abscissa(point3d), 3.8753709794017066)

    def test_point_at_abcissa(self):
        point_at_abscissa = self.ellipse3d.point_at_abscissa(self.ellipse3d.length() * 0.4)
        self.assertTrue(point_at_abscissa.is_close(volmdlr.Point3D(-0.236335934849, -1.104105421298, -1.104105421298)))

    def test_trim(self):
        trim = self.ellipse3d.trim(volmdlr.Point3D(-0.12975651199692162, 0.9309036597828997, 0.9309036597828997),
                                   volmdlr.Point3D(0.12975651199692095, -0.9309036597829001, -0.9309036597829001))
        self.assertAlmostEqual(trim.length(), 4.844224110273849)
        self.assertTrue(trim.point_belongs(volmdlr.Point3D(1.386280848895, 0.495371104295, 0.495371104295)))

    def test_rotation(self):
        rotated_ellipse3d = self.ellipse3d.rotation(volmdlr.O3D, self.ellipse3d.frame.v, math.pi / 2)
        rotated_ellipse3d_points = rotated_ellipse3d.discretization_points(number_points=6)
        expected_points = [volmdlr.Point3D(-3.825416359218455e-16, -1.414213562373095, 1.4142135623730954),
                           volmdlr.Point3D(-0.7071067811865475, -0.35355339059327356, 1.0606601717798214),
                           volmdlr.Point3D(-0.7071067811865472, 1.0606601717798212, -0.3535533905932735),
                           volmdlr.Point3D(2.825496434870558e-16, 1.414213562373095, -1.4142135623730954),
                           volmdlr.Point3D(0.7071067811865474, 0.3535533905932742, -1.060660171779822),
                           volmdlr.Point3D(0.7071067811865475, -1.0606601717798207, 0.3535533905932728)]
        for point, expected_point in zip(rotated_ellipse3d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_traslation(self):
        translated_ellipse3d = self.ellipse3d.translation(self.ellipse3d.frame.w)
        translated_ellipse3d_points = translated_ellipse3d.discretization_points(number_points=6)
        expected_points = [volmdlr.Point3D(1.1547005383792517, 1.8618073195657994, 0.4475937571927041),
                           volmdlr.Point3D(-0.12975651199692162, 1.6380104409694474, 0.22379687859635217),
                           volmdlr.Point3D(-1.2844570503761734, 0.48330990259019585, -0.9309036597828992),
                           volmdlr.Point3D(-1.1547005383792517, -0.44759375719270406, -1.8618073195657994),
                           volmdlr.Point3D(0.12975651199692095, -0.2237968785963525, -1.6380104409694478),
                           volmdlr.Point3D(1.2844570503761727, 0.9309036597828986, -0.48330990259019657)]
        for point, expected_point in zip(translated_ellipse3d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_frame_mapping(self):
        frame_mapped_ellipse3d = self.ellipse3d.frame_mapping(self.ellipse3d.frame, 'new')
        frame_mapped_ellipse3d_points = frame_mapped_ellipse3d.discretization_points(number_points=6)
        expected_points = [volmdlr.Point3D(1.9999999999999996, 0.0, 0.0),
                           volmdlr.Point3D(1.0, -0.8660254037844384, 0.0),
                           volmdlr.Point3D(-0.9999999999999993, -0.8660254037844386, 0.0),
                           volmdlr.Point3D(-1.9999999999999996, -1.224646799147353e-16, 0.0),
                           volmdlr.Point3D(-1.0000000000000007, 0.8660254037844382, 0.0),
                           volmdlr.Point3D(0.9999999999999983, 0.8660254037844388, 0.0)]
        for point, expected_point in zip(frame_mapped_ellipse3d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_from_step(self):
        arguments = ["''", 85, '5.', '3.']
        object_dict = {85: volmdlr.OXYZ}
        ellipse = curves.Ellipse3D.from_step(arguments, object_dict)
        self.assertEqual(ellipse.major_dir, volmdlr.X3D)
        self.assertEqual(ellipse.normal, volmdlr.Z3D)
        self.assertEqual(ellipse.major_axis, 5.0)
        self.assertEqual(ellipse.minor_axis, 3.0)


if __name__ == '__main__':
    unittest.main()
