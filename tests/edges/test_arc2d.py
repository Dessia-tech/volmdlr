import math
import unittest
import os
import volmdlr
from volmdlr.edges import Arc2D
from volmdlr import curves

folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'arc_objects')


class TestArc2D(unittest.TestCase):
    circle2d = curves.Circle2D(volmdlr.OXY, 1)
    arc2d = Arc2D(circle2d, volmdlr.Point2D(-1, 0), volmdlr.Point2D(1, 0))
    arc1 = Arc2D(circle2d, volmdlr.Point2D(0, -1), volmdlr.Point2D(0, 1))
    arc2 = Arc2D(circle2d, volmdlr.Point2D(1, 0), volmdlr.Point2D(-1, 0))
    arc3 = Arc2D(curves.Circle2D(volmdlr.OXY, 1.5), 1.5 * volmdlr.Point2D(0, -1), 1.5 * volmdlr.Point2D(0, 1))
    arc4 = Arc2D(circle2d, volmdlr.Point2D(0.7071067811865475, -0.7071067811865475),
                 volmdlr.Point2D(0.7071067811865475, 0.7071067811865475))
    arc5 = Arc2D(circle2d, volmdlr.Point2D(-0.7071067811865475, 0.7071067811865475),
                 volmdlr.Point2D(-0.7071067811865475, -0.7071067811865475))
    arc6 = arc4.complementary()
    circle2d_reversed = circle2d.reverse()
    arc7 = Arc2D(circle2d_reversed, volmdlr.Point2D(-0.7071067811865475, 0.7071067811865475),
                 volmdlr.Point2D(0.7071067811865475, 0.7071067811865475))
    arc8 = arc7.complementary()
    arc9 = Arc2D(circle2d_reversed, volmdlr.Point2D(-0.7071067811865475, -0.7071067811865475),
                 volmdlr.Point2D(0.7071067811865475, -0.7071067811865475))
    arc10 = arc9.complementary()
    list_points = [volmdlr.Point2D(1, 0),
                   volmdlr.Point2D(0.7071067811865475, -0.7071067811865475),
                   volmdlr.Point2D(0, -1),
                   volmdlr.Point2D(-0.7071067811865475, -0.7071067811865475),
                   volmdlr.Point2D(-1, 0),
                   volmdlr.Point2D(-0.7071067811865475, 0.7071067811865475),
                   volmdlr.Point2D(0, 1),
                   volmdlr.Point2D(0.7071067811865475, 0.7071067811865475)]

    def test_arc_angles(self):
        frame = volmdlr.Frame2D(volmdlr.Point2D(-8.300842483198, 33.937305367912), volmdlr.Vector2D(-1.0, 0.0),
                                volmdlr.Vector2D(-0.0, -1.0))

        radius = 18.985621313014995

        circle = curves.Circle2D(frame, radius)

        start = volmdlr.Point2D(-21.897737241977318, 47.18790043228133)
        end = volmdlr.Point2D(0.5512862323650758, 50.732948667420146)
        arc = Arc2D(circle, start, end)
        self.assertAlmostEqual(arc.angle_start, 2.369092548558922)
        self.assertAlmostEqual(arc.angle_end, 1.085744328897209)
        expected_angles = [2.369092548558922, 2.9246300027275742, -2.8030178502833594, -2.247480396114707,
                           -1.6919429419460545, -1.136405487777402, -0.5808680336087492, -0.02533057944009658,
                           0.5302068747285559, 1.0857443288972086]
        for i, point in enumerate(arc.discretization_points(number_points=10)):
            self.assertAlmostEqual(arc._arc_point_angle(point), expected_angles[i])

        arc1 = Arc2D(self.circle2d, volmdlr.Point2D(1, 0), volmdlr.Point2D(0, 1))
        arc2 = Arc2D(self.circle2d, volmdlr.Point2D(0, 1), volmdlr.Point2D(-1, 0))
        arc3 = Arc2D(self.circle2d, volmdlr.Point2D(-1, 0), volmdlr.Point2D(0, -1))
        arc4 = Arc2D(self.circle2d, volmdlr.Point2D(0, -1), volmdlr.Point2D(1, 0))
        circle_reversed = self.circle2d.reverse()
        arc5 = Arc2D(circle_reversed, volmdlr.Point2D(1, 0), volmdlr.Point2D(0, -1))
        arc6 = Arc2D(circle_reversed, volmdlr.Point2D(0, -1), volmdlr.Point2D(-1, 0))
        arc7 = Arc2D(circle_reversed, volmdlr.Point2D(-1, 0), volmdlr.Point2D(0, 1))
        arc8 = Arc2D(circle_reversed, volmdlr.Point2D(0, 1), volmdlr.Point2D(1, 0))

        expected_angles = [(0.0, 0.5 * math.pi), (0.5 * math.pi, math.pi), (-math.pi, -0.5 * math.pi),
                           (-0.5 * math.pi, 0.0),
                           (0.0, -0.5 * math.pi), (-0.5 * math.pi, -math.pi), (math.pi, 0.5 * math.pi),
                           (0.5 * math.pi, 0.0)
                           ]
        for arc, angles in zip([arc1, arc2, arc3, arc4, arc5, arc6, arc7, arc8], expected_angles):
            self.assertAlmostEqual(arc.angle_start, angles[0])
            self.assertAlmostEqual(arc.angle_end, angles[1])
            self.assertAlmostEqual(arc.angle, 0.5 * math.pi)

    def test_angle(self):
        arc1 = Arc2D(self.circle2d, volmdlr.Point2D(-1/math.sqrt(2), 1/math.sqrt(2)),
                     volmdlr.Point2D(-1/math.sqrt(2), -1/math.sqrt(2)))

        circle_reversed = self.circle2d.reverse()
        arc2 = Arc2D(circle_reversed, volmdlr.Point2D(-1/math.sqrt(2), -1/math.sqrt(2)),
                     volmdlr.Point2D(-1/math.sqrt(2), 1/math.sqrt(2)))


        expected_angles = [(0.75 * math.pi, -0.75 * math.pi), (-0.75 * math.pi, 0.75 * math.pi)
                           ]
        for arc, angles in zip([arc1, arc2], expected_angles):
            self.assertAlmostEqual(arc.angle_start, angles[0])
            self.assertAlmostEqual(arc.angle_end, angles[1])
            self.assertAlmostEqual(arc.angle, 0.5 * math.pi)

    def test_split(self):
        arc_split1 = self.arc2d.split(self.arc2d.start)
        self.assertIsNone(arc_split1[0])
        self.assertEqual(arc_split1[1], self.arc2d)
        arc_split2 = self.arc2d.split(self.arc2d.end)
        self.assertEqual(arc_split2[0], self.arc2d)
        self.assertIsNone(arc_split2[1])
        arc_split3 = self.arc2d.split(volmdlr.Point2D(0, -1))
        self.assertTrue(arc_split3[0].start.is_close(self.arc2d.start))
        self.assertTrue(arc_split3[0].end.is_close(volmdlr.Point2D(0, -1)))
        self.assertTrue(arc_split3[1].start.is_close(volmdlr.Point2D(0, -1)))
        self.assertTrue(arc_split3[1].end.is_close(self.arc2d.end))

    def test_arc_intersections(self):
        arc2 = Arc2D(curves.Circle2D(volmdlr.OXY.translation(volmdlr.Point2D(0, 1.5)), 1), volmdlr.Point2D(-1, 1.5),
                     volmdlr.Point2D(1, 1.5), True)
        arc_intersections = self.arc1.arc_intersections(arc2)
        self.assertEqual(len(arc_intersections), 1)
        self.assertTrue(arc_intersections[0].is_close(volmdlr.Point2D(0.6614378277661477, 0.75)))
        self.assertFalse(self.arc5.arc_intersections(arc2))

    def test_abscissa(self):
        expected_abscissa_results = [[1.5707963267948966, 0.7853981633974483, 0, 3.141592653589793, 2.356194490192345],
                                     [0, 3.141592653589793, 2.356194490192345, 1.5707963267948966, 0.7853981633974483],
                                     [], [0.7853981633974483, 0, 1.5707963267948966],
                                     [1.5707963267948966, 0.7853981633974483, 0],
                                     [4.71238898038469, 3.9269908169872414, 3.141592653589793, 2.356194490192345,
                                      1.5707963267948966, 0.7853981633974483, 0],
                                     [0, 0.7853981633974483, 1.5707963267948966],
                                     [0.7853981633974483, 1.5707963267948966, 2.356194490192345, 3.141592653589793,
                                      3.9269908169872414, 4.71238898038469, 0],
                                     [3.9269908169872414, 4.71238898038469, 0, 0.7853981633974483, 1.5707963267948966,
                                      2.356194490192345, 3.141592653589793],
                                     [0, 0.7853981633974483, 1.5707963267948966]]
        list_abscissas = []
        for arc in [self.arc1, self.arc2, self.arc3, self.arc4, self.arc5, self.arc6, self.arc7,
                    self.arc8, self.arc9, self.arc10]:
            abscissas_ = []
            for point in self.list_points:
                try:
                    abscissa = arc.abscissa(point)
                except ValueError:
                    continue
                abscissas_.append(abscissa)
            list_abscissas.append(abscissas_)
        for abscissas, expected_abscissas in zip(list_abscissas, expected_abscissa_results):
            for abscissa, expected_abscissa in zip(abscissas, expected_abscissas):
                self.assertAlmostEqual(abscissa, expected_abscissa)

    def test_point_belongs(self):
        expected_results = [[True, True, True, False, False, False, True, True],
                            [True, False, False, False, True, True, True, True],
                            [False, False, False, False, False, False, False, False],
                            [True, True, False, False, False, False, False, True],
                            [False, False, False, True, True, True, False, False],
                            [False, True, True, True, True, True, True, True],
                            [False, False, False, False, False, True, True, True],
                            [True, True, True, True, True, True, False, True],
                            [True, True, False, True, True, True, True, True],
                            [False, True, True, True, False, False, False, False]]
        list_point_belongs = []
        for arc in [self.arc1, self.arc2, self.arc3, self.arc4, self.arc5, self.arc6, self.arc7,
                    self.arc8, self.arc9, self.arc10]:
            point_belongs_ = []
            for point in self.list_points:
                point_belongs_.append(arc.point_belongs(point))
            list_point_belongs.append(point_belongs_)
        for result_list, expected_result_list in zip(list_point_belongs, expected_results):
            self.assertEqual(result_list, expected_result_list)

        arc = Arc2D.from_json(os.path.join(folder, "arc2d_bug_point_belongs.json"))
        point = volmdlr.Point2D(0.01330629098214331, 0.0032923224261096617)
        self.assertTrue(arc.point_belongs(point))

        arc = Arc2D.from_json(os.path.join(folder, "arc2d_point_belongs.json"))
        point = volmdlr.Point2D(0.0007151488183559929, 0.007258543823331798)
        self.assertTrue(arc.point_belongs(point))

    def test_get_shared_section(self):
        # =====================Sharing one end of the arc=====================#
        shared_section1 = self.arc1.get_shared_section(self.arc2)
        expected_points = [volmdlr.Point2D(1.0, 0.0), volmdlr.Point2D(0.7071067811865476, 0.7071067811865475),
                           volmdlr.Point2D(0.0, 1.0)]
        for expected_point, point in zip(expected_points, [shared_section1[0].points[0],
                                                           shared_section1[0].middle_point(),
                                                           shared_section1[0].points[1]]):
            self.assertTrue(expected_point.is_close(point))
        # =====================Two Arcs with different radius =====================#
        self.assertFalse(self.arc2.get_shared_section(self.arc3))
        # =====================Two Arcs not touching each other =====================#
        self.assertFalse(self.arc4.get_shared_section(self.arc5))
        # =====================Second Arc Over the first =====================#
        self.assertEqual(self.arc1.get_shared_section(self.arc4)[0], self.arc4)
        # =====================First Arc Over the Second =====================#
        self.assertEqual(self.arc4.get_shared_section(self.arc1)[0], self.arc4)
        # =====================Two arcs juste touching each other =====================#
        self.assertFalse(self.arc4.get_shared_section(self.arc6))

    def test_delete_shared_section(self):
        remaining_arc1 = self.arc1.delete_shared_section(self.arc2)
        self.assertEqual(remaining_arc1, [Arc2D(curves.Circle2D(volmdlr.OXY, 1), volmdlr.Point2D(0.0, -1.0),
                                                volmdlr.Point2D(1.0, 0.0), True)])
        self.assertEqual(self.arc2.delete_shared_section(self.arc3), [self.arc2])
        remaining_arc2 = self.arc1.delete_shared_section(self.arc4)
        self.assertTrue(remaining_arc2[0].start.is_close(volmdlr.Point2D(0.0, -1.0)))
        self.assertTrue(remaining_arc2[0].middle_point().is_close(volmdlr.Point2D(0.3826834323650898, -0.9238795325112867)))
        self.assertTrue(remaining_arc2[0].end.is_close(volmdlr.Point2D(0.7071067811865475, -0.7071067811865475)))
        self.assertTrue(remaining_arc2[1].start.is_close(volmdlr.Point2D(0.7071067811865475, 0.7071067811865475)))
        self.assertTrue(remaining_arc2[1].middle_point().is_close(volmdlr.Point2D(0.38268343236508984, 0.9238795325112867)))
        self.assertTrue(remaining_arc2[1].end.is_close(volmdlr.Point2D(0.0, 1.0)))
        self.assertFalse(self.arc4.delete_shared_section(self.arc1))

    def test_point_distance(self):
        arc = self.arc2

        point1 = volmdlr.Point2D(1, 1)
        self.assertEqual(arc.point_distance(point1), math.sqrt(2) - 1)

        point2 = volmdlr.Point2D(0.5/math.sqrt(2), 0.5/math.sqrt(2))
        self.assertEqual(arc.point_distance(point2), 0.5)

        point3 = volmdlr.Point2D(0, 0)
        self.assertEqual(arc.point_distance(point3), 1)

        point4 = volmdlr.Point2D(0, -1)
        self.assertEqual(arc.point_distance(point4), math.sqrt(2))

    def test_rotation(self):
        rotated_arc2d = self.arc4.rotation(volmdlr.O2D, math.pi / 4)
        self.assertTrue(rotated_arc2d.start.is_close(volmdlr.Point2D(0.9999999999999999, -1.1102230246251565e-16)))
        self.assertTrue(rotated_arc2d.end.is_close(volmdlr.Point2D(1.1102230246251565e-16, 0.9999999999999999)))

        # Test on non trigo Arc2D
        rotated_arc2d = self.arc9.rotation(volmdlr.O2D, math.pi / 4)
        rotated_arc2d = rotated_arc2d.rotation(volmdlr.O2D, -math.pi / 4)
        self.assertTrue(rotated_arc2d.is_close(self.arc9))

        # Test on trigo Arc2D
        rotated_arc2d = self.arc2.rotation(volmdlr.O2D, math.pi / 2)
        rotated_arc2d = rotated_arc2d.rotation(volmdlr.O2D, -math.pi / 2)
        self.assertTrue(rotated_arc2d.is_close(self.arc2))

        arc2d = Arc2D.from_json(os.path.join(folder, "arc2d_rotation_test.json"))
        rotated_arc2d = arc2d.rotation(volmdlr.Point2D(0.5, 0.5), math.pi / 1.5)
        self.assertEqual(arc2d.frame.u, rotated_arc2d.frame.u)
        self.assertEqual(arc2d.frame.v, rotated_arc2d.frame.v)
        self.assertTrue(rotated_arc2d.center.is_close(volmdlr.Point2D(-0.17631397208144173, 0.5714101615137755)))
        self.assertTrue(rotated_arc2d.start.is_close(volmdlr.Point2D(-0.4330127018922192, 1.1160254037844388)))
        self.assertTrue(rotated_arc2d.end.is_close(volmdlr.Point2D(-0.3660254037844387, 2.220446049250313e-16)))

    def test_translation(self):
        translated_arc2d = self.arc4.translation(volmdlr.Vector2D(1, 1))
        self.assertTrue(translated_arc2d.start.is_close(volmdlr.Point2D(1.7071067811865475, 0.29289321881345254)))
        self.assertTrue(translated_arc2d.end.is_close(volmdlr.Point2D(1.7071067811865475, 1.7071067811865475)))

        # Test on non trigo Arc2D
        translated_arc2d = self.arc9.translation(volmdlr.Vector2D(1, 1))
        translated_arc2d = translated_arc2d.translation(volmdlr.Vector2D(-1, -1))
        self.assertTrue(translated_arc2d.is_close(self.arc9))

        # Test on trigo Arc2D
        translated_arc2d = self.arc2.translation(volmdlr.Vector2D(1, 1))
        translated_arc2d = translated_arc2d.translation(volmdlr.Vector2D(-1, -1))
        self.assertTrue(translated_arc2d.is_close(self.arc2))

    def test_frame_mapping(self):
        u_vector = volmdlr.Vector2D(0.7071067811865475, 0.7071067811865475)
        v_vector = volmdlr.Vector2D(-0.7071067811865475, 0.7071067811865475)
        frame = volmdlr.Frame2D(volmdlr.O2D, u_vector, v_vector)
        frame_mapped_arc2d = self.arc4.frame_mapping(frame, 'new')
        self.assertTrue(frame_mapped_arc2d.start.is_close(volmdlr.Point2D(0.0, -1.0)))
        self.assertTrue(frame_mapped_arc2d.end.is_close(volmdlr.Point2D(1.0, 0.0)))

        # Test on non trigo Arc2D
        frame_mapped_arc2d = self.arc9.frame_mapping(frame, 'new')
        frame_mapped_arc2d = frame_mapped_arc2d.frame_mapping(frame, 'old')
        self.assertTrue(frame_mapped_arc2d.is_close(self.arc9))

        # Test on trigo Arc2D
        frame_mapped_arc2d = self.arc2.frame_mapping(frame, 'new')
        frame_mapped_arc2d = frame_mapped_arc2d.frame_mapping(frame, 'old')
        self.assertTrue(frame_mapped_arc2d.is_close(self.arc2))

    def test_reverse(self):
        reverse = self.arc4.reverse()
        self.assertTrue(self.arc4.end.is_close(reverse.start))
        self.assertTrue(self.arc4.start.is_close(reverse.end))

    def test_complementary(self):
        complementary = self.arc4.complementary()
        self.assertTrue(self.arc4.end.is_close(complementary.start))
        self.assertTrue(self.arc4.start.is_close(complementary.end))
        self.assertAlmostEqual(complementary.length(), 4.71238898038469)

    def test_center_of_mass(self):
        center_of_mass = self.arc4.center_of_mass()
        self.assertTrue(center_of_mass.is_close(volmdlr.Point2D(0.6002108774380707, -6.663679357622141e-17)))

    def test_straight_line_second_moment_area(self):
        pass

    def test_straight_line_center_of_mass(self):
        center_of_mass = self.arc4.straight_line_center_of_mass()
        self.assertTrue(center_of_mass.is_close(volmdlr.Point2D(0.8258716790243482, -9.169017534386686e-17)))

    def test_straight_line_point_belongs(self):
        straight_line_point_belongs = self.arc4.straight_line_point_belongs(volmdlr.Point2D(0.8258716790243482, -0.2))
        self.assertTrue(straight_line_point_belongs)
        self.assertFalse(self.arc4.straight_line_point_belongs(volmdlr.Point2D(0.8258716790243482, 0.8)))

    def test_to_3d(self):
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        arc3d = self.arc4.to_3d(volmdlr.O3D, vector1, vector2)
        self.assertTrue(arc3d.start.is_close(volmdlr.Point3D(-0.169101978726, 0.696923425059, 0.696923425059)))
        self.assertTrue(arc3d.end.is_close(volmdlr.Point3D(0.985598559653, 0.119573155869, 0.119573155869)))
        point3d_ = self.arc4.middle_point().to_3d(volmdlr.O3D, vector1, vector2)
        self.assertTrue(arc3d.point_belongs(point3d_))
        circle = curves.Circle2D(
            volmdlr.OXY.translation(volmdlr.Vector2D(0.2068381066975619, 0.1167563813274402)),
            0.01500000000000002).reverse()
        arc = volmdlr.edges.Arc2D(circle,
            volmdlr.Point2D(0.21783000907195643, 0.10654961483693107),
            volmdlr.Point2D(0.19291095633428304, 0.11118552118212867))
        point3d = arc.middle_point().to_3d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D)
        arc_to_3d = arc.to_3d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D)
        self.assertTrue(arc_to_3d.point_belongs(point3d))

    def test_infinite_primitive(self):
        infinite_primitive = self.arc4.infinite_primitive(0.5)
        self.assertTrue(infinite_primitive.start.is_close(volmdlr.Point2D(1.0606601717798214, -1.0606601717798214)))
        self.assertTrue(infinite_primitive.end.is_close(volmdlr.Point2D(1.0606601717798214, 1.0606601717798214)))

    def test_from_3_points(self):
        points_ = [volmdlr.Point2D(-1, 0), volmdlr.Point2D(0, -1), volmdlr.Point2D(1, 0)]
        arc_from_3_points = Arc2D.from_3_points(*points_)
        arc_from_3_points2 = Arc2D.from_3_points(*points_[::-1])
        for arc in [arc_from_3_points, arc_from_3_points2]:
            for point in points_:
                self.assertTrue(arc.point_belongs(point))

        self.assertEqual(arc_from_3_points.rotation(volmdlr.O2D, 10).rotation(volmdlr.O2D, -10), arc_from_3_points)
        self.assertEqual(arc_from_3_points2.rotation(volmdlr.O2D, 10).rotation(volmdlr.O2D, -10), arc_from_3_points2)
        self.assertEqual(arc_from_3_points.translation(volmdlr.X2D).translation(-volmdlr.X2D), arc_from_3_points)
        self.assertEqual(arc_from_3_points2.translation(volmdlr.X2D).translation(-volmdlr.X2D), arc_from_3_points2)


if __name__ == '__main__':
    unittest.main()
