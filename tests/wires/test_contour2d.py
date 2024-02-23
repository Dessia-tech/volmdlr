import math
import os
import unittest

from dessia_common.core import DessiaObject

import volmdlr
from volmdlr import edges, wires, curves
from volmdlr.models.contours import contour2d_1, contour2d_2, contour1_cut_by_wire, contour2_cut_by_wire,\
    contour2_unittest, unordered_contour2_unittest, invalid_unordered_contour2_unittest


folder = os.path.dirname(os.path.realpath(__file__))


class TestContour2D(unittest.TestCase):
    contour1 = wires.Contour2D([edges.FullArc2D(circle=curves.Circle2D(volmdlr.OXY, 0.029999999),
                                                start_end=volmdlr.Point2D(0.029999999, 0))])
    not_ordered_contour = DessiaObject.from_json(os.path.join(folder, "contour_not_ordered.json"))
    # ordered_contour = DessiaObject.from_json('wires/contour_ordered.json')
    contour_to_extract_from = contour = wires.Contour2D.from_points(
        [volmdlr.Point2D(-.15, .15), volmdlr.Point2D(-.15, -.15), volmdlr.Point2D(.15, -.15),
         volmdlr.Point2D(.15, .15)])
    point1_ = volmdlr.Point2D(0.12500000000000003, 0.15)
    point2_ = volmdlr.Point2D(0.12500000000000003, -0.15)
    point_to_extract_with = [(point1_, point2_), (volmdlr.Point2D(0.15, -0.05), volmdlr.Point2D(0.15, 0.05)),
                       (volmdlr.Point2D(-0.15, 0.15), point2_), (volmdlr.Point2D(-0.15, 0.15), point1_)]
    contour3 = contour2_unittest.rotation(volmdlr.Point2D(0.5, 0.5), math.pi / 1.5)
    contour3 = contour3.translation(volmdlr.Vector2D(-0.3, 0))

    def test_point_inside(self):
        point1 = volmdlr.Point2D(0.0144822, 0.00595264)
        point2 = volmdlr.Point2D(0.02, 0.02)
        self.assertTrue(self.contour1.point_inside(point1))
        self.assertTrue(self.contour1.point_inside(point2))

        point3 = volmdlr.Point2D(0, 0.013)
        self.assertTrue(contour2d_1.point_inside(point3))
        self.assertFalse(contour2d_1.point_inside(point1))

        point4 = volmdlr.Point2D(0.745, 0.0685)
        self.assertTrue(contour2d_2.point_inside(point4))

        contour, point = DessiaObject.from_json(os.path.join(folder, "test_contour_point_belongs.json")).primitives
        self.assertTrue(contour.point_inside(point, False))

    def test_is_ordered(self):
        # self.assertTrue(self.ordered_contour.is_ordered())
        self.assertFalse(self.not_ordered_contour.is_ordered())

    def test_order_contour(self):
        ordered_contour = self.not_ordered_contour.order_contour()
        self.assertTrue(self.not_ordered_contour.is_ordered())
        for previous_primitive, primitive in zip(ordered_contour.primitives, ordered_contour.primitives[1:] +
                                                                             [ordered_contour.primitives[0]]):
            self.assertEqual(previous_primitive.end, primitive.start)
        self.assertFalse(unordered_contour2_unittest.is_ordered())
        unordered_contour2_unittest.order_contour()
        self.assertTrue(unordered_contour2_unittest.is_ordered())
        with self.assertRaises(NotImplementedError):
            invalid_unordered_contour2_unittest.order_contour()

    def test_cut_by_wire(self):
        results = contour1_cut_by_wire.cut_by_wire(contour2_cut_by_wire)
        results1 = contour2_cut_by_wire.cut_by_wire(contour1_cut_by_wire)
        list_expected_contour_lengths = [0.735061566418825, 3.4786699386591753, 0.7350615900834909, 0.7350613283792926,
                                         0.7350615482415725, 1.2716033189138256, 3.478669938659176, 0.8996337337370333,
                                         1.2716033047094752, 0.8996336040021796]
        self.assertEqual(len(results) + len(results1), 10)
        for i, contour in enumerate(results + results1):
            self.assertAlmostEqual(contour.length(), list_expected_contour_lengths[i])
        contour1 = wires.ClosedPolygon2D([volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 0),
                                          volmdlr.Point2D(1, 1), volmdlr.Point2D(0, 1)])
        contour1 = wires.Contour2D(contour1.line_segments)
        contour2 = contour1.translation(volmdlr.Vector2D(0.3, 0))
        wire_crossings = contour1.cut_by_wire(contour2)
        self.assertEqual(len(wire_crossings), 2)
        self.assertAlmostEqual(wire_crossings[0].length(), 2.6)
        self.assertAlmostEqual(wire_crossings[1].length(), 3.4)
        contour2 = wires.Contour2D([contour2.primitives[0]] + contour2.primitives[2:])
        contour2 = contour2.order_wire()
        wire_crossings = contour1.cut_by_wire(contour2)
        self.assertEqual(len(wire_crossings), 2)
        self.assertAlmostEqual(wire_crossings[0].length(), 2.6)
        self.assertAlmostEqual(wire_crossings[1].length(), 3.4)

    def test_offset(self):
        contour_to_offset = DessiaObject.from_json(os.path.join(folder, "contour_to_offset.json"))
        stringer_contour_offset = contour_to_offset.offset(4)
        self.assertEqual(len(stringer_contour_offset.primitives), 10)
        self.assertAlmostEqual(stringer_contour_offset.area(), 546.1486677646164)

    def test_edge_crossings(self):
        points = [volmdlr.Point2D(-0.3, -0.2), volmdlr.Point2D(0.3, -0.2),
                  volmdlr.Point2D(0.2, 0.2), volmdlr.Point2D(0, 0.35), volmdlr.Point2D(-0.2, 0.2)]
        contour = wires.ClosedPolygon2D(points)
        lineseg = edges.LineSegment2D(volmdlr.Point2D(-.3, -.3), volmdlr.Point2D(.3, .3))
        lineseg1 = edges.LineSegment2D(volmdlr.Point2D(-.3, -.3), volmdlr.Point2D(.2, .2))
        edge_crossings = contour.edge_crossings(lineseg)
        edge_crossings1 = contour.edge_crossings(lineseg1)
        wire = wires.Wire2D([edges.LineSegment2D(p1, p2) for p1, p2 in (points[:2], points[1:3])])
        wire_edge_crossings = wire.edge_crossings(lineseg)
        self.assertEqual(len(edge_crossings), 2)
        self.assertTrue(edge_crossings[0].is_close(volmdlr.Point2D(-0.2, -0.2)))
        self.assertTrue(edge_crossings[1].is_close(volmdlr.Point2D(0.2, 0.2)))
        self.assertEqual(len(edge_crossings1), 1)
        self.assertTrue(edge_crossings1[0].is_close(volmdlr.Point2D(-0.2, -0.2)))
        self.assertEqual(len(wire_edge_crossings), 1)
        self.assertTrue(wire_edge_crossings[0].is_close(volmdlr.Point2D(-0.2, -0.2)))

    def test_crossings(self):
        contour_crossings = contour2_unittest.wire_crossings(self.contour3)
        expected_crossings = [volmdlr.Point2D(0.003455042474764269, 0.010365521626940934),
                              volmdlr.Point2D(-0.08592141662920678, -0.26441019657119236),
                              volmdlr.Point2D(1.3565385114925572, 0.4261540459702289),
                              volmdlr.Point2D(1.062649150008717, 1.6296941038180754),
                              volmdlr.Point2D(-0.8413589025398691, 1.0502025995770785),
                              volmdlr.Point2D(-0.16551404839994044, -0.5695209737217725),
                              volmdlr.Point2D(0.3497299962331255, -1.1558059537665162),
                              volmdlr.Point2D(0.5885392751147405, -1.116467059799847)]
        self.assertEqual(len(contour_crossings), len(expected_crossings))
        for crossing, expected_crossing in zip(contour_crossings, expected_crossings):
            self.assertTrue(crossing.is_close(expected_crossing))
        points = [volmdlr.Point2D(-0.3, -0.2), volmdlr.Point2D(0.3, -0.2),
                  volmdlr.Point2D(0.2, 0.2), volmdlr.Point2D(0, 0.3), volmdlr.Point2D(-0.2, 0.2)]
        contour = volmdlr.wires.ClosedPolygon2D(points)
        primitives = [volmdlr.edges.LineSegment2D(volmdlr.Point2D(-0.35, -0.1), volmdlr.Point2D(-0.1, 0)),
                      volmdlr.edges.LineSegment2D(volmdlr.Point2D(-0.1, 0), volmdlr.Point2D(0.2, 0.2)),
                      volmdlr.edges.LineSegment2D(volmdlr.Point2D(0.2, 0.2), volmdlr.Point2D(0.3, 0.3))]
        wire = volmdlr.wires.Wire2D(primitives)
        wire2 = volmdlr.wires.Wire2D(primitives[:-1])
        wire3 = volmdlr.wires.Wire2D(primitives[:-1] + [volmdlr.edges.LineSegment2D(
            volmdlr.Point2D(0.2, 0.2), volmdlr.Point2D(0.1, 0.0))])
        wire_crossings = contour.wire_crossings(wire)
        wire_crossings2 = contour.wire_crossings(wire2)
        wire_crossings3 = contour.wire_crossings(wire3)
        self.assertEqual(len(wire_crossings), 2)
        self.assertTrue(wire_crossings[0].is_close(volmdlr.Point2D(-0.26666666666666666, -0.0666666666666666)))
        self.assertTrue(wire_crossings[1].is_close(volmdlr.Point2D(0.2, 0.19999999999999996)))
        self.assertEqual(len(wire_crossings2), 1)
        self.assertTrue(wire_crossings2[0].is_close(volmdlr.Point2D(-0.26666666666666666, -0.0666666666666666)))
        self.assertEqual(len(wire_crossings3), 1)
        self.assertTrue(wire_crossings3[0].is_close(volmdlr.Point2D(-0.26666666666666666, -0.0666666666666666)))

    def test_split_with_two_points(self):
        expected_results = [(3, 0.3499999999999999, 3, 0.85),
                            (1, 0.1, 5, 1.0999999999999999),
                            (2, 0.575, 3, 0.625),
                            (4, 0.9249999999999998, 1, 0.275)]
        for i, (pt1, pt2) in enumerate(self.point_to_extract_with):
            inside_prims, outside_prims = self.contour_to_extract_from.split_with_two_points(pt1, pt2)
            expected_inside_ = expected_results[i][:2]
            expected_outside_ = expected_results[i][2:]
            self.assertEqual(len(inside_prims), expected_inside_[0])
            self.assertAlmostEqual(sum(prim.length() for prim in inside_prims), expected_inside_[1])
            self.assertEqual(len(outside_prims), expected_outside_[0])
            self.assertAlmostEqual(sum(prim.length() for prim in outside_prims), expected_outside_[1])
            contour1 = wires.ClosedPolygon2D([volmdlr.Point2D(0, 0), volmdlr.Point2D(1, 0),
                                              volmdlr.Point2D(1, 1), volmdlr.Point2D(0, 1)])
            contour2 = contour1.translation(volmdlr.Vector2D(0.3, 0))
            self.assertFalse(contour1.wire_crossings(contour2))

    def test_extract_with_points(self):
        list_expected_outside_params = [(3, 0.85), (5, 1.0999999999999999), (3, 0.625), (1, 0.275)]
        list_expected_inside_params = [(3, 0.3499999999999999), (1, 0.1), (2, 0.575), (4, 0.9249999999999998)]
        for i, (pt1, pt2) in enumerate(self.point_to_extract_with):
            inside_prims = self.contour_to_extract_from.extract_with_points(pt1, pt2, inside=True)
            outside_prims = self.contour_to_extract_from.extract_with_points(pt1, pt2, inside=False)
            self.assertEqual(len(inside_prims), list_expected_inside_params[i][0])
            self.assertAlmostEqual(sum(prim.length() for prim in inside_prims), list_expected_inside_params[i][1])
            self.assertEqual(len(outside_prims), list_expected_outside_params[i][0])
            self.assertAlmostEqual(sum(prim.length() for prim in outside_prims), list_expected_outside_params[i][1])

    def test_split_by_line(self):
        line = curves.Line2D(volmdlr.Point2D(volmdlr.TWO_PI, 0.1), volmdlr.Point2D(volmdlr.TWO_PI, -0.1))
        contour = wires.Contour2D.from_json(os.path.join(folder, "contour_to_split.json"))
        intersection = contour.line_intersections(line)[0][0]
        contour1, contour2 = contour.split_by_line(line)
        self.assertTrue(contour1.primitives[-1].end.is_close(intersection))
        self.assertTrue(contour2.primitives[0].start.is_close(intersection))

    def test_closest_point_to_point2(self):
        point1 = volmdlr.Point2D(1.5, -1.5)
        point2 = volmdlr.Point2D(-1, -1)
        closest_point1 = contour2_unittest.closest_point_to_point2(point1)
        self.assertEqual(closest_point1, volmdlr.Point2D(1.0, -1.0))
        closest_point2 = contour2_unittest.closest_point_to_point2(point2)
        self.assertEqual(closest_point2, volmdlr.Point2D(-2.0, 0.7))

    def test_furthest_point_to_point2(self):
        point1 = volmdlr.Point2D(1.5, -1.5)
        point2 = volmdlr.Point2D(-1, -1)
        furthest_point1 = contour2_unittest.get_furthest_point_to_point2(point1)
        self.assertEqual(furthest_point1, volmdlr.Point2D(-2.0, 1.0))
        furthest_point2 = contour2_unittest.get_furthest_point_to_point2(point2)
        self.assertEqual(furthest_point2, volmdlr.Point2D(1.5, 1.0))

    def test_intersection_contour_with(self):
        vol = DessiaObject.from_json(os.path.join(folder, "test_intersection_contour_with.json"))
        contour2 = vol.primitives[0]
        contour3 = vol.primitives[1]
        intersection_contours1 = contour2.intersection_contour_with(contour3, abs_tol=1e-5)
        self.assertTrue(len(intersection_contours1), 2)
        self.assertAlmostEqual(intersection_contours1[0].length(), 0.16514108581676357, 4)
        intersection_contours2 = contour2_unittest.intersection_contour_with(self.contour3, abs_tol=1e-6)
        self.assertTrue(len(intersection_contours1), 2)
        self.assertAlmostEqual(intersection_contours2[0].length(), 6.915890339970204, 6)
        self.assertAlmostEqual(intersection_contours2[1].length(), 2.440847693749909, 6)

    def test_contours_from_edges(self):
        source_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                     'test_contour2d_contours_from_edges_json_files')
        expected_contour_areas = [[0.032085522557644186, 0.0855613934860544],
                                  [0.06136684795257308, 0.0055788043592579495]]
        expected_contour_lengths = [[0.8277671086559999, 1.642237241348], [1.05404105005, 0.315748639158]]
        contour_areas = []
        contour_lengths = []
        for filename in [
            '1_test_contours_from_edges.json',
            '2_test_contours_from_edges2.json']:
            if '.json' not in filename:
                continue
            file_path = os.path.join(source_folder, filename)
            obj = DessiaObject.from_json(file_path)
            primitives = obj.primitives
            contours = wires.Contour2D.contours_from_edges(primitives)
            areas = []
            lengths = []
            for contour in contours:
                areas.append(contour.area())
                lengths.append(contour.length())
            contour_lengths.append(lengths)
            contour_areas.append(areas)
        for solution, expected_solution in zip(contour_areas, expected_contour_areas):
            self.assertEqual(len(solution), len(expected_solution))
            for contour_area, expected_contour_area in zip(solution, expected_solution):
                self.assertAlmostEqual(contour_area, expected_contour_area)
        for solution, expected_solution in zip(contour_lengths, expected_contour_lengths):
            for contour_length, expected_contour_length in zip(solution, expected_solution):
                self.assertAlmostEqual(contour_length, expected_contour_length)

    def test_divide(self):
        vol = DessiaObject.from_json(os.path.join(folder, "test_contour2d_divide_1.json"))
        contour, cutting_contours = vol.primitives[0], vol.primitives[1:]
        divided_contours = contour.divide(cutting_contours)
        divided_contours = sorted(divided_contours, key=lambda cntr: cntr.area())
        expected_contour_areas = [0.0006684483866419249, 0.002005345159820638, 0.002005345159845676,
                                  0.0033422419331328506]
        expected_contour_lengths = [0.10347088858400005, 0.21026602115199994, 0.210266021154, 0.3137369097360001]
        self.assertEqual(len(divided_contours), 4)
        for i, contour_ in enumerate(divided_contours):
            self.assertAlmostEqual(contour_.area(), expected_contour_areas[i])
            self.assertAlmostEqual(contour_.length(), expected_contour_lengths[i])

    def test_merge_not_adjacent_contour(self):
        contours = DessiaObject.from_json(os.path.join(folder, "test_merge_connected_contours.json")).primitives
        contour1, contour2 = contours
        merge_not_adjacent_contour = contour2.merge_not_adjacent_contour(contour1)
        self.assertAlmostEqual(merge_not_adjacent_contour.length(), 0.1589126915239475)
        merge_not_adjacent_contour = contour1.merge_not_adjacent_contour(contour2)
        self.assertAlmostEqual(merge_not_adjacent_contour.length(), 0.1589126915239475)
        contour2 = volmdlr.wires.Contour2D.from_points(
            [volmdlr.Point2D(-0.014284827066811355, 0.00644881122080076), volmdlr.Point2D(-0.008, 0.0060),
             volmdlr.Point2D(-0.01438263879891807, 0.005871523081583764)])
        contours[0] = contour2
        contour1, contour2 = contours
        merge_not_adjacent_contour1 = contour2.merge_not_adjacent_contour(contour1)
        self.assertAlmostEqual(merge_not_adjacent_contour1.length(), 0.15238337009535752)

    def test_area(self):
        contour = wires.Contour2D.from_json(os.path.join(folder, "strange_contour_from_step_file.json"))
        self.assertAlmostEqual(contour.area(), 0.00016865275423510724, 6)

    def test_cut_by_line(self):
        contour, line = wires.Contour2D.from_json(
            os.path.join(folder, 'test_contour2d_cut_by_line.json')).primitives

        cut_by_line = contour.cut_by_line(line)
        self.assertEqual(len(cut_by_line), 1)
        self.assertEqual(cut_by_line[0], contour)


if __name__ == '__main__':
    unittest.main()
