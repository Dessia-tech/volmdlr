import os
import math
import unittest

from dessia_common.core import DessiaObject
from geomdl import utilities

import volmdlr.edges
from volmdlr import curves, wires

circle = curves.Circle2D(volmdlr.OXY, 0.50)
line = curves.Line2D(volmdlr.O2D, volmdlr.Point2D(0, 1))
folder = os.path.join(os.path.dirname(os.path.realpath(__file__)))
folder_2 = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'circle2d_objects')


class TestCircle2D(unittest.TestCase):
    circle2d = curves.Circle2D(volmdlr.OXY, 1)

    def test_discretization_points(self):
        points = self.circle2d.discretization_points(number_points=5)
        expected_points = [volmdlr.Point2D(1.0, 0.0), volmdlr.Point2D(0, 1.0),
                           volmdlr.Point2D(-1.0, 0), volmdlr.Point2D(0, -1.0), volmdlr.Point2D(1.0, 0)]
        for point, expected_point in zip(points, expected_points):
            self.assertTrue(point.is_close(expected_point))

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

        circle_2 = curves.Circle2D(volmdlr.Frame2D(volmdlr.Point2D(0.7137093779940084, 0.0),
                                                   volmdlr.Vector2D(-1.0, 0.0),
                                                   volmdlr.Vector2D(1.5487611193520934e-13, 1.0)), 0.15231602579123288)
        point1 = volmdlr.Point2D(0.8060483808152039, -0.12113496716812525)
        abcissa1 = circle_2.abscissa(point1)
        abscissa_point = circle_2.point_at_abscissa(abcissa1)
        self.assertTrue(point1.is_close(abscissa_point))

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
        circle_ = curves.Circle2D(volmdlr.OXY.translation(volmdlr.Vector2D(0.8, -0.3)), 0.3)
        line_2 = curves.Line2D(volmdlr.Point2D(0.5599870479815988, -0.12053417263965237),
                               volmdlr.Point2D(0.5593939842065143, -0.1196126662489295))
        circle_line_intersections = circle_.line_intersections(line_2)
        self.assertTrue(circle_line_intersections[0], volmdlr.Point2D(0.864102883983119, -0.5930713569509084))
        self.assertTrue(circle_line_intersections[1], volmdlr.Point2D(0.5598081100790778, -0.12025613775092192))
        line_2 = curves.Line2D(volmdlr.Point2D(.1, 0), volmdlr.Point2D(.1, .1))
        circle_line_intersections_2 = circle.line_intersections(line_2)
        self.assertTrue(circle_line_intersections_2[0].is_close(volmdlr.Point2D(0.1, -0.4898979485566356)))
        self.assertTrue(circle_line_intersections_2[1].is_close(volmdlr.Point2D(0.1, 0.4898979485566356)))
        line_3 = curves.Line2D(volmdlr.Point2D(0, 0.1), volmdlr.Point2D(.1, .1))
        circle_line_intersections_3 = circle.line_intersections(line_3)
        self.assertTrue(circle_line_intersections_3[0].is_close(volmdlr.Point2D(0.4898979485566356, 0.1)))
        self.assertTrue(circle_line_intersections_3[1].is_close(volmdlr.Point2D(-0.4898979485566356, 0.1)))
        circle_, line2d = DessiaObject.from_json(os.path.join(folder_2, 'test_circle2d_line2d_intersections.json')).primitives
        intersections = circle_.line_intersections(line2d)
        self.assertEqual(len(intersections), 2)
        self.assertTrue(intersections[0].is_close(volmdlr.Point2D(0.031959341501134664, -0.12658445641691748)))
        self.assertTrue(intersections[1].is_close(volmdlr.Point2D(0.031959341501134664, 0.12658445641691748)))

        circle_ = curves.Circle2D(volmdlr.OXY, 1)
        line2d = curves.Line2D.from_point_and_vector(volmdlr.Point2D(1, 0), volmdlr.Y2D)

        intersections = circle_.line_intersections(line2d)
        self.assertEqual(len(intersections), 1)
        self.assertTrue(intersections[0].is_close(volmdlr.Point2D(1, 0.0)))

    def test_linesegment_intersections(self):
        linesegment = volmdlr.edges.LineSegment2D(volmdlr.Point2D(-2.0, -2.0),
                                                  volmdlr.Point2D(2.0, 2.0))
        linesegment_intersections = self.circle2d.linesegment_intersections(linesegment)
        linesegment_intersections[0].is_close(volmdlr.Point2D(0.7071067811865475, 0.7071067811865475))
        linesegment_intersections[1].is_close(volmdlr.Point2D(-0.7071067811865476, -0.7071067811865476))

    def test_arc_intersections(self):
        arc = volmdlr.edges.Arc2D(curves.Circle2D(volmdlr.OXY.translation(volmdlr.Vector2D(2, 0)), 2),
                                  volmdlr.Point2D(2, 3), volmdlr.Point2D(2, -3))
        arc_intersections = self.circle2d.arc_intersections(arc)
        arc_intersections[0].is_close(volmdlr.Point2D(0.25, -0.9682458365518543))
        arc_intersections[1].is_close(volmdlr.Point2D(0.25, 0.9682458365518543))

    def test_to_3d(self):
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        circle3d = self.circle2d.to_3d(volmdlr.O3D, vector1, vector2)
        circle3d_points = circle3d.discretization_points(number_points=6)
        expected_points = [volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                           volmdlr.Point3D(0.9549454387105718, -0.20985615202546862, -0.20985615202546862),
                           volmdlr.Point3D(0.01283846933518723, -0.7070485038896306, -0.7070485038896306),
                           volmdlr.Point3D(-0.9470108282979028, -0.22712385507308536, -0.22712385507308536),
                           volmdlr.Point3D(-0.598123348937482, 0.5666782417985585, 0.5666782417985585),
                           volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258)]
        for point, expected_point in zip(circle3d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_rotation(self):
        rotated_arc2d = self.circle2d.rotation(volmdlr.Point2D(1, 0), math.pi / 4)
        circle2d_points = rotated_arc2d.discretization_points(number_points=6)
        expected_points = [volmdlr.Point2D(1.2928932188134525, -0.7071067811865475),
                           volmdlr.Point2D(0.6019102131883999, 0.24394973510860607),
                           volmdlr.Point2D(-0.5161237755614949, -0.11932152889407421),
                           volmdlr.Point2D(-0.516123775561495, -1.2948920334790204),
                           volmdlr.Point2D(0.6019102131883997, -1.658163297481701),
                           volmdlr.Point2D(1.2928932188134525, -0.7071067811865475)]
        for point, expected_point in zip(circle2d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

        rotated_arc2d = self.circle2d.rotation(self.circle2d.center, math.pi / 4)
        circle2d_points = rotated_arc2d.discretization_points(number_points=5)
        expected_points = [volmdlr.Point2D(1 / math.sqrt(2), 1 / math.sqrt(2)),
                           volmdlr.Point2D(-1 / math.sqrt(2), 1 / math.sqrt(2)),
                           volmdlr.Point2D(-1 / math.sqrt(2), -1 / math.sqrt(2)),
                           volmdlr.Point2D(1 / math.sqrt(2), -1 / math.sqrt(2)),
                           volmdlr.Point2D(1 / math.sqrt(2), 1 / math.sqrt(2))]
        for point, expected_point in zip(circle2d_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_translation(self):
        translated_arc2d = self.circle2d.translation(volmdlr.Vector2D(1, 1))
        circle2d_points = translated_arc2d.discretization_points(number_points=5)
        expected_points = [volmdlr.Point2D(2.0, 1.0),
                           volmdlr.Point2D(1.0, 2.0),
                           volmdlr.Point2D(0.0, 1.0),
                           volmdlr.Point2D(1.0, 0.0),
                           volmdlr.Point2D(2.0, 1.0)]
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
        circle2d_points = frame_mapped_arc2d.discretization_points(number_points=6)
        expected_points = [volmdlr.Point2D(0.0, -1.4142135623730951),
                           volmdlr.Point2D(0.18389974300182038, -0.2531162814470008),
                           volmdlr.Point2D(-0.8635412462267782, 0.2805815594085902),
                           volmdlr.Point2D(-1.6947951217816852, -0.5506723161463166),
                           volmdlr.Point2D(-1.1610972809260947, -1.5981133053749152),
                           volmdlr.Point2D(0.0, -1.4142135623730954)]
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
            volmdlr.Point2D(0, 0.5), volmdlr.Point2D(-0.5, 0), volmdlr.Point2D(0, -0.5))
        arc2_validate = volmdlr.edges.Arc2D.from_3_points(
            volmdlr.Point2D(0, -0.5), volmdlr.Point2D(0.5, 0), volmdlr.Point2D(0, 0.5))
        self.assertEqual(list_arcs[0], arc1_validate)
        self.assertEqual(list_arcs[1], arc2_validate)

    def test_point_distance(self):
        circle_ = curves.Circle2D(volmdlr.OXY.translation(volmdlr.Vector2D(1.5, 0)), 1)
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
        circle_ = curves.Circle2D(volmdlr.OXY.translation(volmdlr.Vector2D(1.5, 0)), 1)
        circle_intersections = circle_.bsplinecurve_intersections(bspline)
        expected_intersections = [volmdlr.Point2D(0.5410304786421145, 0.2835091834608327),
                                  volmdlr.Point2D(2.4589695213572873, -0.2835091834628551)]
        for intersection, expected_intersection in zip(circle_intersections, expected_intersections):
            self.assertTrue(intersection.is_close(expected_intersection))
        circle_2 = curves.Circle2D(volmdlr.OXY.translation(volmdlr.Vector2D(4.257776625181402, -2.6149184392658222)),
                                   1.1791034225674362)
        bspline = volmdlr.edges.BSplineCurve2D.from_json(os.path.join(folder, 'bspline.json'))
        intersections = circle_2.bsplinecurve_intersections(bspline, 1e-6)
        self.assertEqual(len(intersections), 1)
        self.assertTrue(intersections[0], volmdlr.Point2D(3.218528920632699, -3.1719185197869))

    def test_circle_intersections(self):
        circle1 = curves.Circle2D(volmdlr.OXY.translation(volmdlr.Vector2D(0, 0)), 1)
        circle2 = curves.Circle2D(volmdlr.OXY.translation(volmdlr.Vector2D(1, 1)), 1)
        circle_intersections = circle1.circle_intersections(circle2)
        self.assertEqual(len(circle_intersections), 2)
        self.assertTrue(circle_intersections[0].is_close(volmdlr.Point2D(1.0, 0.0)))
        self.assertTrue(circle_intersections[1].is_close(volmdlr.Point2D(0.0, 1.0)))


if __name__ == '__main__':
    unittest.main()
