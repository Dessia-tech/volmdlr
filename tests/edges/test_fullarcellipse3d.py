import math
import os
import unittest
import volmdlr
import volmdlr.edges as vme
from volmdlr import curves


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'fullarcellipse_objects')


class TestFullArcEllipse3D(unittest.TestCase):
    start_end = volmdlr.Point3D(0.0225, 0.0, 0.0)
    major_axis = 0.0225
    minor_axis = 0.0075
    ellipse3d = curves.Ellipse3D(major_axis, minor_axis, volmdlr.OXYZ)
    ellipse = vme.FullArcEllipse3D(ellipse3d, start_end)

    def test_init(self):
        self.assertAlmostEqual(self.ellipse.ellipse.major_axis, 0.0225, places=4)
        self.assertEqual(self.ellipse.ellipse.minor_axis, 0.0075)
        self.assertAlmostEqual(self.ellipse.theta, 0, places=4)

    def test_length(self):
        self.assertAlmostEqual(self.ellipse.length(), 0.10023669584870037)

    def test_to_2d(self):
        plane_origin = volmdlr.Point3D(1, 1, 0)
        x = volmdlr.Vector3D(0.5*math.sqrt(2), 0.5*math.sqrt(2), 0)
        y = volmdlr.Vector3D(-0.5*math.sqrt(2), 0.5*math.sqrt(2), 0)
        ellipse2d = self.ellipse.to_2d(plane_origin, x, y)
        self.assertTrue(ellipse2d.ellipse.major_dir.is_close(volmdlr.Vector2D(0.5*math.sqrt(2), -0.5*math.sqrt(2))))
        self.assertTrue(ellipse2d.ellipse.minor_dir.is_close(volmdlr.Vector2D(0.5*math.sqrt(2), 0.5*math.sqrt(2))))
        self.assertAlmostEqual(ellipse2d.ellipse.major_axis, 0.0225, places=4)
        self.assertAlmostEqual(ellipse2d.ellipse.minor_axis, 0.0075, places=4)

    def test_reverse(self):
        reverse = self.ellipse.reverse()
        self.assertEqual(self.ellipse.ellipse.frame.w.dot(reverse.ellipse.frame.w), -1)

    def test_frame_mapping(self):
        new_frame = volmdlr.Frame3D(volmdlr.O3D, volmdlr.Z3D, volmdlr.X3D, -volmdlr.Y3D)
        new_ellipse = self.ellipse.frame_mapping(new_frame, 'new')
        self.assertEqual(new_ellipse.ellipse.major_dir, volmdlr.Vector3D(0.0, 1.0, 0.0))
        self.assertEqual(new_ellipse.ellipse.minor_dir, volmdlr.Vector3D(0.0, 0.0, -1.0))
        self.assertEqual(new_ellipse.normal, volmdlr.Vector3D(1.0, 0.0, 0.0))

    def test_abscissa(self):
        point1 = volmdlr.Point3D(0, -0.0075, 0)
        point2 = volmdlr.Point3D(0.0225, 0, 0)
        self.assertAlmostEqual(self.ellipse.abscissa(point1), 0.75*self.ellipse.length())
        self.assertAlmostEqual(self.ellipse.abscissa(point2), 0.0)

    def test_translation(self):
        translated_ellipse = self.ellipse.translation(volmdlr.X3D)
        self.assertEqual(translated_ellipse.ellipse.center, volmdlr.Point3D(1, 0, 0))
        self.assertEqual(translated_ellipse.start_end, volmdlr.Point3D(1.0225, 0, 0))

    def test_point_belongs(self):
        fullarcellipse = vme.FullArcEllipse3D(curves.Ellipse3D(
            major_axis=0.0100500150616, minor_axis=0.009916145846950001, frame=volmdlr.Frame3D(
                volmdlr.Point3D(-0.46362762553200004, -0.509606021605, 0.509000000114),
                volmdlr.Vector3D(-0.8577170218519301, -0.016653392447698645, 0.5138522890339582),
                volmdlr.Vector3D(-0.5141198069940289, 0.030753268830964377, -0.8571668802004851),
                volmdlr.Vector3D(-0.00152790113491038,  -0.9993882633772487, -0.0349394410649187))),
            volmdlr.Point3D(-0.47224769452020265, -0.5097733884499261, 0.514164223358229))
        point = volmdlr.Point3D(-0.45404913959118903, -0.5072162164760068, 0.5060526668248432)
        self.assertTrue(fullarcellipse.point_belongs(point, 1e-5))

    def test_split(self):
        fullarcellipse = vme.FullArcEllipse3D(curves.Ellipse3D(
            major_axis=0.0100500150616, minor_axis=0.009916145846950001, frame=volmdlr.Frame3D(
                volmdlr.Point3D(-0.46362762553200004, -0.509606021605, 0.509000000114),
                volmdlr.Vector3D(-0.8577170218519301, -0.016653392447698645, 0.5138522890339582),
                volmdlr.Vector3D(-0.5141198069940289, 0.030753268830964377, -0.8571668802004851),
                volmdlr.Vector3D(-0.00152790113491038, -0.9993882633772487, -0.0349394410649187))),
            volmdlr.Point3D(-0.47224769452020265, -0.5097733884499261, 0.514164223358229))
        split_point = volmdlr.Point3D(-0.4540526537637985, -0.5095148094812338, 0.5059723061104128)
        result = fullarcellipse.split(split_point)
        self.assertTrue(result[0].start.is_close(fullarcellipse.start_end))
        self.assertTrue(result[0].end.is_close(split_point))
        self.assertTrue(result[1].start.is_close(split_point))
        self.assertTrue(result[1].end.is_close(fullarcellipse.start_end))

    def test_discretization_points(self):
        fullarcellipse = vme.FullArcEllipse3D.from_json(
            os.path.join(folder, "fullarcellipse3d_discretization_points.json"))
        discretization_points = fullarcellipse.discretization_points(number_points=9)
        self.assertEqual(len(discretization_points), 9)
        self.assertTrue(discretization_points[0].is_close(fullarcellipse.start))
        self.assertTrue(discretization_points[-1].is_close(fullarcellipse.end))

    def test_line_intersections(self):
        fullarcellipse = vme.FullArcEllipse3D.from_json(
            os.path.join(folder, "fullarcellipse3d_line_intersections.json"))
        line = curves.Line3D.from_json(os.path.join(folder, "fullarcellipse3d_line_intersections_line.json"))
        test = fullarcellipse.line_intersections(line, 1e-4)[0]
        self.assertTrue(test, volmdlr.Point3D(0.3407914925119553, -0.10964172421958009, 0.5033056993640009))


if __name__ == '__main__':
    unittest.main()
