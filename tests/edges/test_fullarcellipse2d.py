import unittest
import volmdlr
import volmdlr.edges as vme
from volmdlr import curves


class TestFullArcEllipse2D(unittest.TestCase):
    start_end = volmdlr.Point2D(0.0225, 0.0)
    major_axis = 0.0225
    minor_axis = 0.0075
    center = volmdlr.O2D
    major_dir = volmdlr.X2D
    ellispe2d = curves.Ellipse2D(major_axis, minor_axis, volmdlr.OXY)
    fullarcellipse = vme.FullArcEllipse2D(ellispe2d, start_end)

    def test_init(self):
        self.assertAlmostEqual(self.fullarcellipse.ellipse.major_axis, 0.0225, places=4)
        self.assertAlmostEqual(self.fullarcellipse.ellipse.minor_axis, 0.0075, places=4)
        self.assertEqual(self.fullarcellipse.theta, 0.0)

    def test_length(self):
        self.assertAlmostEqual(self.fullarcellipse.length(), 0.10023669584870037)

    def test_to_3d(self):
        plane_origin = volmdlr.Point3D(1, 1, 1)
        x = volmdlr.Y3D
        y = volmdlr.Z3D
        ellipse3d = self.fullarcellipse.to_3d(plane_origin, x, y)
        self.assertEqual(ellipse3d.ellipse.major_dir, x)
        self.assertEqual(ellipse3d.ellipse.minor_dir, y)
        self.assertEqual(ellipse3d.normal, volmdlr.X3D)
        self.assertAlmostEqual(ellipse3d.ellipse.major_axis, 0.0225, places=4)
        self.assertAlmostEqual(ellipse3d.ellipse.minor_axis, 0.0075, places=4)

    def test_reverse(self):
        reverse = self.fullarcellipse.reverse()
        self.assertEqual(self.fullarcellipse.ellipse.frame.v.dot(reverse.ellipse.frame.v), -1)

    def test_frame_mapping(self):
        new_frame = volmdlr.Frame2D(volmdlr.O2D, -volmdlr.Y2D, volmdlr.X2D)
        new_ellipse = self.fullarcellipse.frame_mapping(new_frame, 'new')
        self.assertEqual(new_ellipse.ellipse.major_dir, volmdlr.Vector2D(0.0, 1.0))
        self.assertEqual(new_ellipse.ellipse.minor_dir, volmdlr.Vector2D(-1.0, 0.0))

    def test_abscissa(self):
        point1 = volmdlr.Point2D(0, -0.0075)
        point2 = volmdlr.Point2D(0.0225, 0)
        self.assertAlmostEqual(self.fullarcellipse.abscissa(point1), 0.75*self.fullarcellipse.length())
        self.assertAlmostEqual(self.fullarcellipse.abscissa(point2), 0.0)

        ellipse = vme.FullArcEllipse2D(curves.Ellipse2D(0.000500289037421, 0.00050027520242, volmdlr.OXY),
                                       volmdlr.Point2D(0.0005002890374210534, 0))
        point = volmdlr.Point2D(-0.00018416867811365376, 0.00046514411968310123)
        self.assertAlmostEqual(ellipse.abscissa(point), 0.00098248885770749, 4)

    def test_translation(self):
        translated_ellipse = self.fullarcellipse.translation(volmdlr.X2D)
        self.assertEqual(translated_ellipse.ellipse.center, volmdlr.Point2D(1, 0))
        self.assertEqual(translated_ellipse.start_end, volmdlr.Point2D(1.0225, 0))


if __name__ == '__main__':
    unittest.main()
