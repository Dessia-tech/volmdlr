import unittest
import volmdlr
from volmdlr import curves


class TestHyperbola2D(unittest.TestCase):
    def test_line_intersections(self):
        # Usage
        a = 2  # semi-major axis
        b = 1  # semi-minor axis
        u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
        v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
        frame1 = volmdlr.Frame2D(volmdlr.O2D, u_vector, v_vector)
        frame2 = volmdlr.OXY
        hyperbola1 = curves.Hyperbola2D(frame1, a, b)
        hyperbola2 = curves.Hyperbola2D(frame2, a, b)
        line1 = curves.Line2D(volmdlr.Point2D(1, 0), volmdlr.Point2D(-6, 10))
        line2 = curves.Line2D(volmdlr.Point2D(-10, -5), volmdlr.Point2D(10, 5))
        line3 = curves.Line2D(volmdlr.Point2D(-10, -5), volmdlr.Point2D(5, 10))
        expected_results = [[[volmdlr.Point2D(3.414716821729969, 1.4411665257000275),
                              volmdlr.Point2D(1.4066104613103017, 3.2667177624451798)],
                             [volmdlr.Point2D(2.529822128134703, 1.2649110640673518)],
                             [volmdlr.Point2D(2.696152422706632, 7.696152422706632)]],
                            [[volmdlr.Point2D(10.884604683797278, -5.349640621633888),
                              volmdlr.Point2D(3.452312878926738, 1.40698829188478319)],
                             [],
                             []]]
        for i, hyperbola in enumerate([hyperbola1, hyperbola2]):
            for j, line in enumerate([line1, line2, line3]):
                line_intersections = hyperbola.line_intersections(line)
                for intersection, expected_result in zip(line_intersections, expected_results[i][j]):
                    self.assertTrue(intersection.is_close(expected_result))

    def test_point_belongs(self):
        hyperbola = curves.Hyperbola2D(volmdlr.OXY, 1, 1)
        point1 = volmdlr.Point2D(1.846035698527472, -1.5517241379310343)
        point2 = volmdlr.Point2D(4.4248245764887875, 4.310344827586208)
        point3 = volmdlr.Point2D(1.4248245764887875, -2.410344827586208)
        self.assertTrue(hyperbola.point_belongs(point1))
        self.assertTrue(hyperbola.point_belongs(point2))
        self.assertFalse(hyperbola.point_belongs(point3))

    def test_tangent(self):
        hyperbola = curves.Hyperbola2D(volmdlr.OXY, 1, 1)
        point1 = volmdlr.Point2D(1.846035698527472, -1.5517241379310343)
        point2 = volmdlr.Point2D(4.4248245764887875, 4.310344827586208)
        tangent_vector1 = hyperbola.tangent(point1)
        tangent_vector2 = hyperbola.tangent(point2)
        self.assertTrue(tangent_vector1.is_close(volmdlr.Vector2D(-0.8405710350936326, 1.0)))
        self.assertTrue(tangent_vector2.is_close(volmdlr.Vector2D(0.9741278446357251, 1.0)))


class TestHyperbola3D(unittest.TestCase):
    def test_line_intersections(self):
        a = 2  # semi-major axis
        b = 1  # semi-minor axis

        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        vector3 = vector1.cross(vector2)
        frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)
        hyperbola = curves.Hyperbola3D(frame, a, b)

        points = hyperbola.get_points(number_points=400)
        line3d_1 = curves.Line3D(points[20], points[250])
        line3d_2 = curves.Line3D(points[20], points[320])
        line3d_3 = curves.Line3D(points[50], volmdlr.Point3D(10, 15, -15))
        line3d_4 = curves.Line3D(volmdlr.Point3D(-20, -15, 15), volmdlr.Point3D(10, 15, -15))
        expected_results = [[volmdlr.Point3D(3.106958942196619, 14.126593248702392, 14.126593248702392),
                             volmdlr.Point3D(5.209563472536312, 2.1093320938257762, 2.1093320938257762)],
                            [volmdlr.Point3D(3.1069589421966235, 14.126593248702383, 14.126593248702383),
                             volmdlr.Point3D(12.001168238296865, 4.603586433650937, 4.603586433650937)],
                            [volmdlr.Point3D(2.6111150368875498, 11.789027732278152, 11.789027732278152)],
                            []]

        for i, line in enumerate([line3d_1, line3d_2, line3d_3, line3d_4]):
            intersections = hyperbola.line_intersections(line)
            for intersection, expected_result in zip(intersections, expected_results[i]):
                self.assertTrue(intersection.is_close(expected_result))

    def test_trim(self):
        hyperbola = curves.Hyperbola3D(volmdlr.Frame3D(volmdlr.Point3D(0, 0.25, 0),
                                                       volmdlr.Z3D, volmdlr.X3D, volmdlr.Y3D),
                                       0.4330127018922194, 0.25000000000016875)

        point_start = volmdlr.Point3D(0.4330127018922191, 0.25, 0.866025403784)
        point_end = volmdlr.Point3D(-0.4330127018922191, 0.25, 0.866025403784)
        bspline = hyperbola.trim(point_start, point_end)
        self.assertAlmostEqual(bspline.length(), 1.259818722659198, 5)

    def test_point_belongs(self):
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        vector3 = vector1.cross(vector2)
        frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)
        hyperbola3d = curves.Hyperbola3D(frame, 1, 1)
        point1 = volmdlr.Point3D(4.181191857716377, 0.5914224070858578, 0.5914224070858578)
        point2 = volmdlr.Point3D(-0.8766643972395096, 3.9800825065547225, 3.9800825065547225)
        point3 = volmdlr.Point3D(2.3259646383580663, 1.0961993491191613, -3.0961993491191613)
        self.assertTrue(hyperbola3d.point_belongs(point1))
        self.assertTrue(hyperbola3d.point_belongs(point2))
        self.assertFalse(hyperbola3d.point_belongs(point3))

    def test_tangent(self):
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        vector3 = vector1.cross(vector2)
        frame = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)
        hyperbola3d = curves.Hyperbola3D(frame, 1, 1)
        point1 = volmdlr.Point3D(4.181191857716377, 0.5914224070858578, 0.5914224070858578)
        point2 = volmdlr.Point3D(-0.8766643972395096, 3.9800825065547225, 3.9800825065547225)
        tangent_vector1 = hyperbola3d.tangent(point1)
        tangent_vector2 = hyperbola3d.tangent(point2)
        self.assertTrue(tangent_vector1.is_close(
            volmdlr.Vector3D(1.362919855421623, 0.13817498403014217, 0.13817498403014217)))
        self.assertTrue(tangent_vector2.is_close(
            volmdlr.Vector3D(0.2566720737495096, -0.9680727976427224, -0.9680727976427224)))


if __name__ == '__main__':
    unittest.main()
