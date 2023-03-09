import unittest
import volmdlr
import volmdlr.edges as vme


class TestFullArcEllipse3D(unittest.TestCase):
    start_end = volmdlr.Point3D(0.0225, 0.0, 0.0)
    major_axis = 0.0225
    minor_axis = 0.0075
    center = volmdlr.O3D
    major_dir = volmdlr.X3D
    normal = volmdlr.Z3D
    ellipse = vme.FullArcEllipse3D(start_end, major_axis, minor_axis, center, normal, major_dir)

    def test_init(self):
        self.assertAlmostEqual(self.ellipse.Gradius, 0.0225, places=4)
        self.assertEqual(self.ellipse.Sradius, 0.0075)
        self.assertAlmostEqual(self.ellipse.angle, volmdlr.TWO_PI, places=4)

    def test_reverse(self):
        self.assertEqual(self.ellipse, self.ellipse.reverse())


if __name__ == '__main__':
    unittest.main()
