import unittest

import volmdlr
from volmdlr.edges import Arc2D


class TestArc2D(unittest.TestCase):
    arc2d = Arc2D(volmdlr.Point2D(-1, 0), volmdlr.Point2D(0, -1), volmdlr.Point2D(1, 0))
    arc1 = Arc2D(volmdlr.Point2D(0, -1), volmdlr.Point2D(1, 0), volmdlr.Point2D(0, 1))
    arc2 = Arc2D(volmdlr.Point2D(1, 0), volmdlr.Point2D(0, 1), volmdlr.Point2D(-1, 0))
    arc3 = Arc2D(1.5 * volmdlr.Point2D(0, -1), 1.5 * volmdlr.Point2D(1, 0), 1.5 * volmdlr.Point2D(0, 1))
    arc4 = Arc2D(volmdlr.Point2D(0.7071067811865475, -0.7071067811865475), volmdlr.Point2D(1, 0),
                 volmdlr.Point2D(0.7071067811865475, 0.7071067811865475))

    arc5 = Arc2D(volmdlr.Point2D(-0.7071067811865475, 0.7071067811865475), volmdlr.Point2D(-1, 0),
                 volmdlr.Point2D(-0.7071067811865475, -0.7071067811865475))
    arc6 = arc4.complementary()

    def test_split(self):
        arc_split1 = self.arc2d.split(self.arc2d.start)
        self.assertIsNone(arc_split1[0])
        self.assertEqual(arc_split1[1], self.arc2d)
        arc_split2 = self.arc2d.split(self.arc2d.end)
        self.assertEqual(arc_split2[0], self.arc2d)
        self.assertIsNone(arc_split2[1])
        arc_split3 = self.arc2d.split(self.arc2d.interior)
        self.assertTrue(arc_split3[0].start.is_close(self.arc2d.start))
        self.assertTrue(arc_split3[0].end.is_close(self.arc2d.interior))
        self.assertTrue(arc_split3[1].start.is_close(self.arc2d.interior))
        self.assertTrue(arc_split3[1].end.is_close(self.arc2d.end))

    def test_get_shared_section(self):
        # =====================Sharing one end of the arc=====================#
        shared_section1 = self.arc1.get_shared_section(self.arc2)
        expected_points = [volmdlr.Point2D(1.0, 0.0), volmdlr.Point2D(0.7071067811865476, 0.7071067811865475),
                           volmdlr.Point2D(0.0, 1.0)]
        for expected_point, point in zip(expected_points, shared_section1[0].points):
            self.assertEqual(expected_point, point)
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
        self.assertEqual(remaining_arc1, [Arc2D(volmdlr.Point2D(0.0, -1.0),
                                                volmdlr.Point2D(0.7071067811865475, -0.7071067811865476),
                                                volmdlr.Point2D(1.0, 0.0))])
        self.assertEqual(self.arc2.delete_shared_section(self.arc3), [self.arc2])
        remaining_arc2 = self.arc1.delete_shared_section(self.arc4)
        self.assertEqual(remaining_arc2[0], Arc2D(volmdlr.Point2D(0.0, -1.0),
                                                  volmdlr.Point2D(0.3826834323650898, -0.9238795325112867),
                                                  volmdlr.Point2D(0.7071067811865475, -0.7071067811865475)))
        self.assertEqual(remaining_arc2[1], Arc2D(volmdlr.Point2D(0.7071067811865475, 0.7071067811865475),
                                                  volmdlr.Point2D(0.38268343236508984, 0.9238795325112867),
                                                  volmdlr.Point2D(0.0, 1.0)))
        self.assertFalse(self.arc4.delete_shared_section(self.arc1))


if __name__ == '__main__':
    unittest.main()
