import unittest
import volmdlr
from volmdlr import curves


class TestParabola2D(unittest.TestCase):
    def test_line_intersections(self):
        expected_results = [[[volmdlr.Point2D(0.43439728232830777, 4.150547925156083),
                              volmdlr.Point2D(-1407.1677289150666, 1283.7888444682421)],
                             [volmdlr.Point2D(0.0, 0.0),
                              volmdlr.Point2D(-1.2570787221094177, -0.6285393610547088)],
                             [volmdlr.Point2D(0.15914794847249425, 5.159147948472494),
                              volmdlr.Point2D(-5.159147948472494, -0.15914794847249425)]],
                            [[volmdlr.Point2D(2.8172904669025316, 1.984281393724971),
                              volmdlr.Point2D(-6.453654103266168, 10.412412821151062)],
                             [volmdlr.Point2D(2.0, 1.0), volmdlr.Point2D(0.0, 0.0)],
                             [volmdlr.Point2D(6.898979485566356, 11.898979485566356),
                              volmdlr.Point2D(-2.8989794855663558, 2.1010205144336442)]]]

        frame1 = volmdlr.Frame2D(volmdlr.O2D, volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475),
                                 volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475))
        frame2 = volmdlr.OXY
        parabola1 = curves.Parabola2D(frame1, 1)
        parabola2 = curves.Parabola2D(frame2, 1)

        line1 = curves.Line2D(volmdlr.Point2D(5, 0), volmdlr.Point2D(-6, 10))
        line2 = curves.Line2D(volmdlr.Point2D(-10, -5), volmdlr.Point2D(10, 5))
        line3 = curves.Line2D(volmdlr.Point2D(-10, -5), volmdlr.Point2D(5, 10))
        ax = parabola1.plot()
        parabola1.frame.plot(ax)
        for i, hyperbola in enumerate([parabola1, parabola2]):
            for j, line in enumerate([line1, line2, line3]):
                line_intersections = hyperbola.line_intersections(line)
                for intersection, expected_result in zip(line_intersections, expected_results[i][j]):
                    self.assertTrue(intersection.is_close(expected_result))


class TestParabola3D(unittest.TestCase):
    def test_line_intersections(self):
        expected_results = [[[volmdlr.Point3D(2.77966655993497, 0.8078693778774468, 0.8078693778774468),
                             volmdlr.Point3D(0.9512856633182674, -3.9353925603457114, -3.9353925603457114)],
                            [volmdlr.Point3D(-0.17313627385849073, -2.067615656403906, -2.067615656403906)],
                            []],
                            [[volmdlr.Point3D(2.5376884422110555, 1.6099656574328935, 0.0),
                             volmdlr.Point3D(-3.9949748743718594, 3.9899560617156133, 0.0)],
                            [volmdlr.Point3D(-2.487437185929648, 1.5468359384864017, 0.0)],
                            []]]
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        vector3 = vector1.cross(vector2)
        frame1 = volmdlr.Frame3D(volmdlr.O3D, vector1, vector2, vector3)
        frame2 = volmdlr.OXYZ
        parabola1 = curves.Parabola3D(frame1, 1)
        parabola2 = curves.Parabola3D(frame2, 1)
        for i, parabola in enumerate([parabola1, parabola2]):
            points = parabola.get_points()
            line3d_1 = curves.Line3D(points[20], points[150])
            line3d_3 = curves.Line3D(points[50], volmdlr.Point3D(1, 1.5, -1.5))
            line3d_4 = curves.Line3D(volmdlr.Point3D(-2., -1.5, 1.5), volmdlr.Point3D(1.0, 1.5, -1.5))
            for j, line in enumerate([line3d_1, line3d_3, line3d_4]):
                line_intersections = parabola.line_intersections(line)
                for intersection, expected_result in zip(line_intersections, expected_results[i][j]):
                    self.assertTrue(intersection.is_close(expected_result))

    def test_split(self):
        parabola = curves.Parabola3D(volmdlr.Frame3D(
            volmdlr.Point3D(-0.43301270189243873, 0.0, 0.7500000000003803),
            volmdlr.Vector3D(-0.0, 1.0, 0.0), volmdlr.Vector3D(0.5000000000002298, 0.0, 0.8660254037843059),
            volmdlr.Vector3D(0.8660254037844387, 0.0, -0.49999999999999994)), 0.21650635094600354)
        point_start = volmdlr.Point3D(1.6339745962174324, -1.8921223583379627, 4.330127018924)
        point_end = volmdlr.Point3D(1.6339745962174324, 1.8921223583379627, 4.330127018924)
        bspline = parabola.split(point_start, point_end)[0]
        self.assertAlmostEqual(bspline.length(), 9.425433371950165)


if __name__ == '__main__':
    unittest.main()
