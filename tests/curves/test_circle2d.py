import math
import unittest

from geomdl import utilities

import volmdlr.edges
from volmdlr import curves


circle = curves.Circle2D(volmdlr.O2D, 0.50)
line = curves.Line2D(volmdlr.O2D, volmdlr.Point2D(0, 1))


class TestCircle2D(unittest.TestCase):
    circle2d = curves.Circle2D(volmdlr.O2D, 1)

    def test_area(self):
        area = self.circle2d.area()
        self.assertAlmostEqual(area, 3.141592653589793)

    def test_second_moment_area(self):
        expected_result = (3.9269908169872414, 3.9269908169872414, -3.141592653589793)
        for i, result in enumerate(self.circle2d.second_moment_area(volmdlr.Point2D(1, 1))):
            self.assertAlmostEqual(result, expected_result[i])

    def test_center_of_mass(self):
        center_of_masss = self.circle2d.center_of_mass()
        self.assertEqual(center_of_masss, volmdlr.O2D)

    def test_length(self):
        length = self.circle2d.length()
        self.assertAlmostEqual(length, 6.283185307179586)

    def test_point_symmetric(self):
        circle_ = self.circle2d.point_symmetric(volmdlr.Point2D(2, 2))
        self.assertEqual(self.circle2d.radius, circle_.radius)
        self.assertTrue(circle_.center.is_close(volmdlr.Point2D(4.0, 4.0)))

    def test_axial_symmetry(self):
        line_ = curves.Line2D(volmdlr.Point2D(1.5, 0), volmdlr.Point2D(1.5, 1))
        axial_symentry = self.circle2d.axial_symmetry(line_)
        self.assertEqual(self.circle2d.radius, axial_symentry.radius)
        self.assertTrue(axial_symentry.center.is_close(volmdlr.Point2D(3.0, 0.0)))

    def test_point_at_abscissa(self):
        point_at_pi = self.circle2d.point_at_abscissa(math.pi)
        self.assertTrue(point_at_pi.is_close(volmdlr.Point2D(-1.0, 0.0)))

    def test_abscissa(self):
        abscissa = self.circle2d.abscissa(self.circle2d.point_at_abscissa(math.pi))
        self.assertAlmostEqual(abscissa, math.pi)

    def test_point_belongs(self):
        self.assertTrue(self.circle2d.point_belongs(volmdlr.Point2D(-1.0, 0.0)))
        self.assertFalse(self.circle2d.point_belongs(volmdlr.Point2D(-1.0, 1.0)))

    def test_cut_by_line(self):
        line_ = curves.Line2D(volmdlr.Point2D(-1.0, -1.0), volmdlr.Point2D(1.0, 1.0))
        cut_by_line = self.circle2d.cut_by_line(line_)
        expected_points1 = [volmdlr.Point2D(0.7071067811865475, 0.7071067811865475),
                            volmdlr.Point2D(0.05086297082788327, 0.9987056414172105),
                            volmdlr.Point2D(-0.6316101530579021, 0.775286150111153),
                            volmdlr.Point2D(-0.9883708684504969, 0.15206257395694156),
                            volmdlr.Point2D(-0.8354427283221042, -0.5495775174565629),
                            volmdlr.Point2D(0.3316517992283888, 0.3316517992283888),
                            volmdlr.Point2D(-0.18772749097907948, -0.18772749097907948),
                            volmdlr.Point2D(-0.7071067811865476, -0.7071067811865476)]
        expected_points2 = [volmdlr.Point2D(-0.7071067811865476, -0.7071067811865476),
                            volmdlr.Point2D(-0.050862970827883325, -0.9987056414172106),
                            volmdlr.Point2D(0.6316101530579022, -0.7752861501111531),
                            volmdlr.Point2D(0.9883708684504972, -0.15206257395694156),
                            volmdlr.Point2D(0.8354427283221044, 0.5495775174565629),
                            volmdlr.Point2D(0.3316517992283888, 0.3316517992283888),
                            volmdlr.Point2D(-0.18772749097907948, -0.18772749097907948),
                            volmdlr.Point2D(-0.7071067811865476, -0.7071067811865476)]
        for i, point in enumerate(cut_by_line[0].discretization_points(number_points=8)):
            self.assertTrue(point.is_close(expected_points1[i]))
        for i, point in enumerate(cut_by_line[1].discretization_points(number_points=8)):
            self.assertTrue(point.is_close(expected_points2[i]))

    def test_line_intersections(self):
        line_ = curves.Line2D(volmdlr.Point2D(-1.0, -1.0), volmdlr.Point2D(1.0, 1.0))
        line_intersections = self.circle2d.line_intersections(line)
        line_intersections[0].is_close(volmdlr.Point2D(0.7071067811865475, 0.7071067811865475))
        line_intersections[1].is_close(volmdlr.Point2D(-0.7071067811865476, -0.7071067811865476))

    def test_linesegment_intersections(self):
        linesegment = volmdlr.edges.LineSegment2D(volmdlr.Point2D(-2.0, -2.0),
                                                  volmdlr.Point2D(2.0, 2.0))
        linesegment_intersections = self.circle2d.linesegment_intersections(linesegment)
        linesegment_intersections[0].is_close(volmdlr.Point2D(0.7071067811865475, 0.7071067811865475))
        linesegment_intersections[1].is_close(volmdlr.Point2D(-0.7071067811865476, -0.7071067811865476))

    def test_arc_intersections(self):
        arc = volmdlr.edges.Arc2D(curves.Circle2D(volmdlr.Point2D(2, 0), 2),
                                  volmdlr.Point2D(2, 3), volmdlr.Point2D(2, -3))
        arc_intersections = self.circle2d.arc_intersections(arc)
        arc_intersections[0].is_close(volmdlr.Point2D(0.25, -0.9682458365518543))
        arc_intersections[1].is_close(volmdlr.Point2D(0.25, 0.9682458365518543))

    def test_to_3d(self):
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        circle3d = self.circle2d.to_3d(volmdlr.O3D, vector1, vector2)
        circle3d_points = circle3d.discretization_points(number_points=5)
        expected_points = [volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                           volmdlr.Point3D(0.9549454387105718, -0.20985615202546862, -0.20985615202546862),
                           volmdlr.Point3D(0.01283846933518723, -0.7070485038896306, -0.7070485038896306),
                           volmdlr.Point3D(-0.9470108282979028, -0.22712385507308536, -0.22712385507308536),
                           volmdlr.Point3D(-0.598123348937482, 0.5666782417985585, 0.5666782417985585)]
        for point, expected_point in zip(circle3d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_rotation(self):
        rotated_arc2d = self.circle2d.rotation(volmdlr.Point2D(1, 0), math.pi / 4)
        circle2d_points = rotated_arc2d.discretization_points(number_points=5)
        expected_points = [volmdlr.Point2D(1.2928932188134525, -0.7071067811865475),
                           volmdlr.Point2D(0.6019102131883999, 0.24394973510860607),
                           volmdlr.Point2D(-0.5161237755614949, -0.11932152889407421),
                           volmdlr.Point2D(-0.516123775561495, -1.2948920334790204),
                           volmdlr.Point2D(0.6019102131883997, -1.658163297481701)]
        for point, expected_point in zip(circle2d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_translation(self):
        translated_arc2d = self.circle2d.translation(volmdlr.Vector2D(1, 1))
        circle2d_points = translated_arc2d.discretization_points(number_points=5)
        expected_points = [volmdlr.Point2D(2.0, 1.0),
                           volmdlr.Point2D(1.3090169943749475, 1.9510565162951536),
                           volmdlr.Point2D(0.19098300562505266, 1.5877852522924734),
                           volmdlr.Point2D(0.19098300562505255, 0.412214747707527),
                           volmdlr.Point2D(1.3090169943749472, 0.04894348370484636)]
        for point, expected_point in zip(circle2d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_split_at_abscissa(self):
        split_at_abscissa = self.circle2d.split_at_abscissa(self.circle2d.length() * 0.5)
        self.assertTrue(split_at_abscissa[0].is_close(volmdlr.edges.Arc2D(
            self.circle2d, volmdlr.Point2D(1, 0), volmdlr.Point2D(-1, 0))))
        self.assertTrue(split_at_abscissa[1].is_close(volmdlr.edges.Arc2D(
            self.circle2d, volmdlr.Point2D(-1, 0), volmdlr.Point2D(1, 0))))

    def test_frame_mapping(self):
        u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
        v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
        frame = volmdlr.Frame2D(volmdlr.Point2D(0, 1), u_vector, v_vector)
        frame_mapped_arc2d = self.circle2d.frame_mapping(frame, 'new')
        circle2d_points = frame_mapped_arc2d.discretization_points(number_points=5)
        expected_points = [volmdlr.Point2D(0.2928932188134524, -0.7071067811865476),
                           volmdlr.Point2D(-0.3980897868116001, 0.24394973510860596),
                           volmdlr.Point2D(-1.516123775561495, -0.11932152889407432),
                           volmdlr.Point2D(-1.516123775561495, -1.2948920334790206),
                           volmdlr.Point2D(-0.39808978681160034, -1.658163297481701)]
        for point, expected_point in zip(circle2d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_split(self):
        point1 = volmdlr.Point2D(0.7071067811865475, 0.7071067811865475)
        point2 = volmdlr.Point2D(-0.7071067811865475, 0.7071067811865475)
        split = self.circle2d.split(point1, point2)
        self.assertAlmostEqual(split[0].length(), 1.5707963267948966)
        self.assertAlmostEqual(split[1].length(), 4.71238898038469)

    def test_split_by_line(self):
        list_arcs = circle.split_by_line(line=line)
        arc1_validate = volmdlr.edges.Arc2D.from_3_points(
            volmdlr.Point2D(0, -0.5), volmdlr.Point2D(0.5, 0), volmdlr.Point2D(0, 0.5))
        arc2_validate = volmdlr.edges.Arc2D.from_3_points(
            volmdlr.Point2D(0, 0.5), volmdlr.Point2D(-0.5, 0), volmdlr.Point2D(0, -0.5))
        self.assertEqual(list_arcs[0], arc1_validate)
        self.assertEqual(list_arcs[1], arc2_validate)
        # self.assertEqual(list_arcs[0].length(), arc1_validate.length())
        # self.assertEqual(list_arcs[0].start, arc1_validate.start)
        # self.assertEqual(list_arcs[0].end, arc1_validate.end)
        # self.assertEqual(list_arcs[1].length(), arc2_validate.length())
        # self.assertEqual(list_arcs[1].start, arc2_validate.start)
        # self.assertEqual(list_arcs[1].end, arc2_validate.end)

    def test_point_distance(self):
        circle_ = curves.Circle2D(volmdlr.Point2D(1.5, 0), 1)
        point1 = volmdlr.Point2D(0.5410304786421145, 0.2835091834608327)
        self.assertEqual(circle_.point_distance(point1), 0.0)
        point2 = volmdlr.Point2D(2, 1.5)
        self.assertEqual(circle_.point_distance(point2), 0.5811388300841898)

    def test_bspline_intersections(self):
        degree = 3
        points = [volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 1), volmdlr.Point2D(2, -1), volmdlr.Point2D(3, 0)]
        knotvector = utilities.generate_knot_vector(degree, len(points))
        knot_multiplicity = [1] * len(knotvector)
        bspline = volmdlr.edges.BSplineCurve2D(degree, points, knot_multiplicity, knotvector, None, False)
        circle_ = curves.Circle2D(volmdlr.Point2D(1.5, 0), 1)
        circle_intersections = circle_.bsplinecurve_intersections(bspline)
        expected_intersections = [volmdlr.Point2D(0.5410304786421145, 0.2835091834608327),
                                  volmdlr.Point2D(2.4589695213572873, -0.2835091834628551)]
        for intersection, expected_intersection in zip(circle_intersections, expected_intersections):
            self.assertTrue(intersection.is_close(expected_intersection))
        circle_2 = curves.Circle2D(volmdlr.Point2D(4.257776625181402, -2.6149184392658222), 1.1791034225674362)
        bspline = volmdlr.edges.BSplineCurve2D.load_from_file('curves/bspline.json')
        intersections = circle_2.bsplinecurve_intersections(bspline, 1e-6)
        self.assertEqual(len(intersections), 1)
        self.assertTrue(intersections[0], volmdlr.Point2D(3.218528920632699, -3.1719185197869))

    def test_circle_intersections(self):
        circle1 = curves.Circle2D(volmdlr.Point2D(0, 0), 1)
        circle2 = curves.Circle2D(volmdlr.Point2D(1, 1), 1)
        circle_intersections = circle1.circle_intersections(circle2)
        self.assertEqual(len(circle_intersections), 2)
        self.assertTrue(circle_intersections[0].is_close(volmdlr.Point2D(1.0, 0.0)))
        self.assertTrue(circle_intersections[1].is_close(volmdlr.Point2D(0.0, 1.0)))


if __name__ == '__main__':
    unittest.main()
