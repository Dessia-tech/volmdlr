import math
import os
import unittest
from itertools import product

import dessia_common.core
import volmdlr
from volmdlr import edges, wires, curves
from volmdlr.models.curves import circle3d


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "arc_objects")


class TestArc3D(unittest.TestCase):
    list_points = circle3d.discretization_points(number_points=9)
    arc3d = edges.Arc3D(circle3d, start=volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                        end=volmdlr.Point3D(-0.9855985596534886, -0.11957315586905026, -0.11957315586905026))

    arc3d_2 = edges.Arc3D(curves.Circle3D(volmdlr.OXYZ, 1), volmdlr.Point3D(-1, 0, 0), volmdlr.Point3D(1, 0, 0))

    def test_init(self):
        expected_lengths = [0.7853981, 1.5707963, 2.3561945, 3.1415926, 3.9269908, 4.712389, 5.4977871, 5.4977872,
                            0.7853982, 1.5707963, 2.3561945, 3.1415927, 3.9269908, 4.712389, 4.712389, 5.4977871,
                            0.7853982, 1.5707963, 2.3561945, 3.1415927, 3.9269908, 3.9269908, 4.712389, 5.4977871,
                            0.7853981, 1.5707963, 2.3561945, 3.1415927, 3.1415927, 3.9269908, 4.712389, 5.4977872,
                            0.7853982, 1.5707963, 2.3561945, 2.3561945, 3.1415927, 3.9269908, 4.712389, 5.4977871,
                            0.7853982, 1.5707963, 1.5707963, 2.3561945, 3.1415927, 3.9269908, 4.712389, 5.4977871,
                            0.7853982, 0.7853982, 1.5707963, 2.3561945, 3.1415927, 3.9269908, 4.712389, 5.4977871]
        list_lengths = []
        for point1, point2 in product(self.list_points[:-1], repeat=2):
            if not point1.is_close(point2):
                arc3d_ = edges.Arc3D(circle3d, start=point1, end=point2)
                list_lengths.append(round(arc3d_.length(), 7))
        for length, expected_length in zip(list_lengths, expected_lengths):
            self.assertAlmostEqual(length, expected_length, 6)

    def test_from_3_points(self):
        points = [volmdlr.Point3D(-0.2672612419124244, -0.5345224838248488, -0.8017837257372732),
                  volmdlr.Point3D(0.9636241116594316, -0.14824986333222026, -0.22237479499833038),
                  volmdlr.Point3D(0.2672612419124244, 0.5345224838248488, 0.8017837257372732)]
        arc3d_from_3_points = edges.Arc3D.from_3_points(*points)
        self.assertTrue(arc3d_from_3_points.point_belongs(points[1]))
        arc3d_from_3_points = edges.Arc3D.from_3_points(*points[::-1])
        self.assertTrue(arc3d_from_3_points.point_belongs(points[1]))
        self.assertAlmostEqual(arc3d_from_3_points.length(), 3.141592653589793)

    def test_discretization_points(self):
        expected_points = [volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                           volmdlr.Point3D(0.9996507193171729, -0.01868742048356531, -0.01868742048356531),
                           volmdlr.Point3D(0.5334021057067797, -0.5981146184585249, -0.5981146184585249),
                           volmdlr.Point3D(-0.40696605032273603, -0.6459019406553569, -0.6459019406553569),
                           volmdlr.Point3D(-0.9855985560899309, -0.11957317055561367, -0.11957317055561367)]
        discret_points = self.arc3d.discretization_points(number_points=5)
        for point, expected_point in zip(discret_points, expected_points):
            self.assertTrue(point.is_close(expected_point))
        discret_points = self.arc3d.discretization_points(angle_resolution=1)
        for point, expected_point in zip(discret_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_point_at_abscissa(self):
        abscissas = [0, math.pi/4, math.pi/2, 3 * math.pi / 4, math.pi,
                     5 * math.pi / 4, 3 * math.pi / 2, 7 * math.pi / 4]
        points_at_abscissas = []
        for abscissa in abscissas:
            try:
                points_at_abscissas.append(self.arc3d.point_at_abscissa(abscissa))
            except ValueError:
                with self.assertRaises(ValueError):
                    self.arc3d.point_at_abscissa(abscissa)
        for point, expected_point in zip(points_at_abscissas, self.list_points):
            self.assertTrue(point.is_close(expected_point))

    def test_abscissa(self):
        expected_abscissa = [0, math.pi/4, math.pi/2, 3 * math.pi / 4, math.pi, 5 * math.pi / 4, 3 * math.pi / 2,
                             7 * math.pi / 4, 0]
        abscissas = []

        for point in self.list_points[:-1]:
            if self.arc3d.point_belongs(point):
                abscissa = self.arc3d.abscissa(point)
                abscissas.append(abscissa)
            else:
                with self.assertRaises(ValueError):
                    self.arc3d_2.abscissa(point)

        for abscissa, expected_abscissa in zip(abscissas, expected_abscissa):
            self.assertAlmostEqual(abscissa, expected_abscissa)

        arc = edges.Arc3D.from_json(os.path.join(folder, "arc_abscissa_bug.json"))
        point = volmdlr.Point3D(0.6002208894332702, -0.637646466964, 0.006570128575852758)
        abscissa = arc.abscissa(point)
        self.assertAlmostEqual(abscissa, 0.019320794819579237)

    def test_direction_vector(self):
        direction_vector = self.arc3d.direction_vector(self.arc3d.abscissa(self.list_points[2]))
        self.assertTrue(direction_vector.is_close(
            volmdlr.Vector3D(-0.5773502691896258, -0.5773502691896258, -0.5773502691896258)))

    def test_normal_vector(self):
        normal_vector = self.arc3d.normal_vector(self.arc3d.abscissa(self.list_points[2]))
        self.assertTrue(normal_vector.is_close(
            volmdlr.Vector3D(-0.8164965809277261, 0.408248290463863, 0.408248290463863)))

    def test_rotation(self):
        rotated_arc3d = self.arc3d.rotation(volmdlr.O3D, self.arc3d.circle.frame.v, math.pi / 2)
        self.assertTrue(rotated_arc3d.start.is_close(
            volmdlr.Point3D(-1.9127081796092275e-16, -0.7071067811865475, 0.7071067811865477)))
        self.assertTrue(rotated_arc3d.end.is_close(
            volmdlr.Point3D(-0.5773502691896253, 0.7886751345948129, -0.21132486540518736)))

    def test_translation(self):
        translated_arc3d = self.arc3d.translation(self.arc3d.circle.frame.w)
        self.assertTrue(translated_arc3d.start.is_close(
            volmdlr.Point3D(0.5773502691896258, 1.2844570503761734, -0.12975651199692173)))
        self.assertTrue(translated_arc3d.end.is_close(
            volmdlr.Point3D(-0.9855985596534886, 0.5875336253174973, -0.8266799370555978)))

    def test_frame_mapping(self):
        frame_mapped_arc3d = self.arc3d.frame_mapping(self.arc3d.circle.frame, 'new')
        self.assertTrue(frame_mapped_arc3d.start.is_close(
            volmdlr.Point3D(0.9999999999999998, 0.0, 0.0)))
        self.assertTrue(frame_mapped_arc3d.end.is_close(
            volmdlr.Point3D(-0.7071067811865477, -0.7071067811865472, 0.0)))

    def test_to_2d(self):
        arc2d = self.arc3d.to_2d(self.arc3d.circle.center, self.arc3d.circle.frame.u, self.arc3d.circle.frame.v)
        self.assertTrue(arc2d.start.is_close(volmdlr.Point2D(1., 0.0)))
        self.assertTrue(arc2d.end.is_close(volmdlr.Point2D(-0.7071067811865477, -0.7071067811865471)))
        point2d_ = self.arc3d.middle_point().to_2d(
            self.arc3d.circle.center, self.arc3d.circle.frame.u, self.arc3d.circle.frame.v)
        self.assertTrue(arc2d.point_belongs(point2d_))
        arc2d_2 = self.arc3d.to_2d(self.arc3d.circle.center, self.arc3d.circle.frame.u, -self.arc3d.circle.frame.v)
        point2d_2 = self.arc3d.middle_point().to_2d(
            self.arc3d.circle.center, self.arc3d.circle.frame.u, -self.arc3d.circle.frame.v)
        self.assertTrue(arc2d_2.point_belongs(point2d_2))

    def test_minimum_distance_points_arc(self):
        arc3d_2 = edges.Arc3D(circle3d, self.list_points[-2], self.list_points[-1])
        arc3d_2 = arc3d_2.translation(arc3d_2.circle.frame.w)
        minimum_distance_points_arc = self.arc3d.minimum_distance_points_arc(arc3d_2)
        self.assertTrue(minimum_distance_points_arc[0].is_close(
            volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258)) or
                        minimum_distance_points_arc[0].is_close(
                            volmdlr.Point3D(-0.9855985596534887, -0.11957315586905082, -0.11957315586905082)))
        self.assertTrue(minimum_distance_points_arc[1].is_close(
            volmdlr.Point3D(0.5773502532353942, 1.2844570583532893, -0.12975650401980598)) or
                        minimum_distance_points_arc[1].is_close(
                            volmdlr.Point3D(-0.8164965809277258, 1.1153550716504106, -0.29885849072268444)))

    def test_minimum_distance_points_line(self):
        lineseg = edges.LineSegment3D(volmdlr.Point3D(0, 1.5, -1), volmdlr.Point3D(0, -.3, 0.5))
        test_minimum_distance_points_line = self.arc3d.minimum_distance(lineseg, True)
        self.assertTrue(test_minimum_distance_points_line[1].is_close(
            volmdlr.Point3D(0.0, 0.17973768284826552, 0.10021859762644536)))
        self.assertTrue(test_minimum_distance_points_line[2].is_close(
            volmdlr.Point3D(0.577350269233884, 0.5773502691674968, 0.5773502691674968)))

    def test_minimum_distance(self):
        arc3d_2 = edges.Arc3D(circle3d, self.list_points[-2], self.list_points[-1])
        arc3d_2 = arc3d_2.translation(arc3d_2.circle.frame.w)
        minimum_distance = self.arc3d.minimum_distance(arc3d_2)
        self.assertAlmostEqual(minimum_distance, 1.0)
        lineseg = edges.LineSegment3D(volmdlr.Point3D(0, 1.5, -1), volmdlr.Point3D(0, -.3, 0.5))
        minimum_distance = self.arc3d.minimum_distance(lineseg)
        self.assertAlmostEqual(minimum_distance, 0.8479880507244569)

    def test_extrusion(self):
        extrusion = self.arc3d.extrusion(self.arc3d.circle.normal)
        self.assertAlmostEqual(extrusion[0].outer_contour3d.length(), 9.853981633974481)
        extrusion_points = extrusion[0].outer_contour3d.discretization_points(number_points=10)
        expected_points = [volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                           volmdlr.Point3D(0.9902758761944115, -0.09837095360773679, -0.09837095360773679),
                           volmdlr.Point3D(0.3300336747301572, -0.6674869937100307, -0.6674869937100307),
                           volmdlr.Point3D(-0.6878679668410815, -0.5132434413579569, -0.5132434413579569),
                           volmdlr.Point3D(-0.9855985596534889, 0.2004326678135493, -0.43957897955165043),
                           volmdlr.Point3D(-0.7535426476176112, 0.2422554725826881, -1.1719580897904072),
                           volmdlr.Point3D(0.23911314674674883, 0.020511957781626844, -1.3937016045914683),
                           volmdlr.Point3D(0.9726406135015547, 0.5428353031885098, -0.8713782591845853),
                           volmdlr.Point3D(0.6521117099944177, 1.2431806183060983, -0.17103294406699687),
                           volmdlr.Point3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258)]
        for point, expected_point in zip(extrusion_points, expected_points):
            self.assertTrue(point.is_close(expected_point))

    def test_revolution(self):
        revolution_point = self.arc3d.circle.center + self.arc3d.circle.frame.u * 4 * self.arc3d.circle.radius
        revolution = self.arc3d.revolution(revolution_point, self.arc3d.frame.v, math.pi / 4)
        self.assertAlmostEqual(self.arc3d.angle * (math.pi/4), revolution[0].surface2d.area())

    def test_point_belongs(self):
        expected_results = [True, True, True, True, True, True, False, False]
        point_belongs_assertions = [self.arc3d.point_belongs(point) for point in self.list_points]
        for result, expected_result in zip(point_belongs_assertions, expected_results):
            self.assertEqual(result, expected_result)

    def test_linesegment_intersections(self):
        lineseg1 = edges.LineSegment3D(volmdlr.Point3D(2, 2, 0), volmdlr.Point3D(-2, -2, 0))
        lineseg2 = edges.LineSegment3D(volmdlr.Point3D(0, -1, -1), volmdlr.Point3D(0, 1, 0))
        lineseg3 = edges.LineSegment3D(volmdlr.Point3D(-1, -1, -1),
                                       volmdlr.Point3D(-math.sqrt(2) / 2, -math.sqrt(2) / 2, 0))
        inters1 = self.arc3d_2.linesegment_intersections(lineseg1)
        inters2 = self.arc3d_2.linesegment_intersections(lineseg2)
        inters3 = self.arc3d_2.linesegment_intersections(lineseg3)
        expected_point = volmdlr.Point3D(-0.7071067811865476, -0.7071067811865475, 0.0)
        self.assertEqual(len(inters1), 1)
        self.assertTrue(inters1[0].is_close(expected_point))
        self.assertFalse(inters2)
        self.assertEqual(len(inters3), 1)
        self.assertTrue(inters3[0].is_close(expected_point))
        point = volmdlr.Point3D(0, 1.5, -1)
        lineseg = edges.LineSegment3D(point, point + (self.list_points[-5] - point) * 2)
        lineseg2 = edges.LineSegment3D(point, point + (volmdlr.Point3D(1, 1.5, -1) - point) * 2)
        arc3d_lineseg_inters = self.arc3d.linesegment_intersections(lineseg)
        no_arc3d_lineseg_inters = self.arc3d.linesegment_intersections(lineseg2)
        self.assertEqual(arc3d_lineseg_inters[0],
                         volmdlr.Point3D(-0.5773502691896256, -0.5773502691896257, -0.5773502691896257))
        self.assertFalse(no_arc3d_lineseg_inters)

    def test_split(self):
        split1 = self.arc3d.split(self.arc3d.start)
        self.assertIsNone(split1[0])
        self.assertEqual(split1[1], self.arc3d)
        split2 = self.arc3d.split(self.arc3d.end)
        self.assertEqual(split2[0], self.arc3d)
        self.assertIsNone(split2[1])
        split3 = self.arc3d.split(self.arc3d.middle_point())
        self.assertTrue(split3[0].start.is_close(self.arc3d.start))
        self.assertTrue(split3[0].end.is_close(self.arc3d.middle_point()))
        self.assertTrue(split3[1].start.is_close(self.arc3d.middle_point()))
        self.assertTrue(split3[1].end.is_close(self.arc3d.end))

    def test_point_distance(self):
        arc = self.arc3d_2

        point1 = volmdlr.Point3D(-1, -1, 0)
        self.assertEqual(arc.point_distance(point1), math.sqrt(2) - 1)

        point2 = volmdlr.Point3D(-0.5 / math.sqrt(2), -0.5 / math.sqrt(2), 0)
        self.assertEqual(arc.point_distance(point2), 0.5)

        point3 = volmdlr.Point3D(0, 0, 0)
        self.assertEqual(arc.point_distance(point3), 1)

        point4 = volmdlr.Point3D(0, 1, 0)
        self.assertEqual(arc.point_distance(point4), math.sqrt(2))

        point = volmdlr.Point3D(0.33525, -0.51106, 0.639935)
        arc = dessia_common.core.DessiaObject.from_json(os.path.join(folder, "arc3d_test_point_distance.json"))
        self.assertAlmostEqual(arc.point_distance(point), 0.0021097842575404403)

    def test_arc_intersections(self):
        arc1 = edges.Arc3D.from_3_points(
            volmdlr.Point3D(-3, -3, 0),
            volmdlr.Point3D(6.324555320336761, -5.692099788303083, -0.8973665961010275),
            volmdlr.Point3D(3, 3, 2))
        arc2 = edges.Arc3D.from_3_points(
            volmdlr.Point3D(1.2067665579541171, -1.246879774203074, -0.4359328108960321),
            volmdlr.Point3D(-1.2905737351068276, -5.961765089244547, -0.9872550297481824),
            volmdlr.Point3D(7.33260591629263, -4.272128323147327, -0.4240427743824422))
        intersections = arc1.intersections(arc2)
        self.assertTrue(intersections[0].is_close(
            volmdlr.Point3D(-1.2905737351057338, -5.961765089245487, -0.9872550297484957)))
        self.assertTrue(intersections[1].is_close(
            volmdlr.Point3D(7.332605916292026, -4.272128323148586, -0.42404277438286175)))

    def test_minimum_distance_bspline(self):
        arc = edges.Arc3D.from_3_points(volmdlr.Point3D(-3, -3, 0),
                                        volmdlr.Point3D(6.324555320336761, -5.692099788303083, -0.8973665961010275),
                                        volmdlr.Point3D(3, 3, 2))
        bspline = edges.BSplineCurve3D.from_points_interpolation(
            [volmdlr.Point3D(1.2067665579541171, -1.246879774203074, -0.4359328108960321),
             volmdlr.Point3D(-1.2905737351068276, -5.961765089244547, -0.9872550297481824),
             volmdlr.Point3D(2.739350840642852, -5.869347626045908, -0.7880999427201254)], 2, centripetal=True)
        dist, min_dist_point1, min_dist_point2 = bspline.minimum_distance(arc, True)
        self.assertAlmostEqual(dist, 1.9450692583417427e-11)
        self.assertTrue(min_dist_point1.is_close(
            volmdlr.Point3D(-1.29057373513423, -5.961765089185557, -0.9872550297440593)))
        self.assertTrue(min_dist_point2.is_close(
            volmdlr.Point3D(-1.2905737351257018, -5.9617650891813065, -0.9872550297271026)))

    def test_move_frame_along(self):
        origin = self.arc3d.start
        w = self.arc3d.unit_direction_vector(0.)
        u = self.arc3d.unit_normal_vector(0.)
        if not u:
            u = w.deterministic_unit_normal_vector()
        v = w.cross(u)
        starting_frame = volmdlr.Frame3D(origin, u, v, w)
        ending_frame = self.arc3d.move_frame_along(starting_frame)
        self.assertTrue(ending_frame.origin.is_close(self.arc3d.end))
        self.assertTrue(ending_frame.w.is_close(self.arc3d.unit_direction_vector(self.arc3d.length())))


if __name__ == '__main__':
    unittest.main()
