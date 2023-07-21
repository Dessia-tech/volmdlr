"""
Unit tests for volmdlr.edges.LineSegment2D
"""
import unittest

import volmdlr
from volmdlr import edges, curves


class TestLineSegment2D(unittest.TestCase):
    lineseg1 = edges.LineSegment2D(volmdlr.O2D, volmdlr.Point2D(0, 3))
    lineseg2 = edges.LineSegment2D(volmdlr.Point2D(0, 2), volmdlr.Point2D(0, 6))
    lineseg3 = edges.LineSegment2D(volmdlr.Point2D(3, 0), volmdlr.Point2D(3, 6))
    lineseg4 = edges.LineSegment2D(volmdlr.Point2D(3, 2), volmdlr.Point2D(3, 4))
    lineseg5 = edges.LineSegment2D(volmdlr.Point2D(1.5, 0), volmdlr.Point2D(1.5, 2))
    lineseg6 = edges.LineSegment2D(volmdlr.Point2D(1.5, 3), volmdlr.Point2D(1.5, 5))
    lineseg7 = edges.LineSegment2D(volmdlr.Point2D(4.5, 0), volmdlr.Point2D(4.5, 3))
    lineseg8 = edges.LineSegment2D(volmdlr.Point2D(4.5, 3), volmdlr.Point2D(4.5, 5))

    def test_line_intersections(self):
        linesegment2d = edges.LineSegment2D(volmdlr.Point2D(1, 2), volmdlr.Point2D(3, -4))
        line1 = curves.Line2D(volmdlr.Point2D(1, 2), volmdlr.Point2D(3, 4))
        line2 = curves.Line2D(volmdlr.Point2D(3, -4), volmdlr.Point2D(-4, -3))
        line3 = curves.Line2D(volmdlr.Point2D(0, 0), volmdlr.Point2D(2, 2))
        self.assertEqual(linesegment2d.line_intersections(line1)[0], volmdlr.Point2D(1.0, 2.0))
        self.assertTrue(linesegment2d.line_intersections(line2)[0].is_close(volmdlr.Point2D(3, -4)))
        self.assertTrue(linesegment2d.line_intersections(line3)[0].is_close(volmdlr.Point2D(1.25, 1.25)))

    def test_get_shared_section(self):
        shared_section1 = self.lineseg1.get_shared_section(self.lineseg2)
        self.assertEqual(shared_section1[0].start, volmdlr.Point2D(0.0, 3.0))
        self.assertEqual(shared_section1[0].end, volmdlr.Point2D(0.0, 2.0))
        self.assertFalse(self.lineseg1.get_shared_section(self.lineseg3))
        shared_section3 = self.lineseg3.get_shared_section(self.lineseg4)
        self.assertTrue(shared_section3[0], self.lineseg4)
        self.assertTrue(self.lineseg4.get_shared_section(self.lineseg3)[0], self.lineseg4)
        self.assertFalse(self.lineseg5.get_shared_section(self.lineseg6))
        self.assertFalse(self.lineseg7.get_shared_section(self.lineseg8))

    def test_delete_shared_portion(self):
        remaining_segment = self.lineseg1.delete_shared_section(self.lineseg2)
        self.assertEqual(len(remaining_segment), 1)
        self.assertEqual(remaining_segment[0].start, volmdlr.Point2D(0.0, 0.0))
        self.assertEqual(remaining_segment[0].end, volmdlr.Point2D(0.0, 2.0))
        remaining_segment2 = self.lineseg3.delete_shared_section(self.lineseg4)
        self.assertEqual(len(remaining_segment2), 2)
        self.assertEqual(remaining_segment2[0].start, volmdlr.Point2D(3.0, 0.0))
        self.assertEqual(remaining_segment2[0].end, volmdlr.Point2D(3.0, 2.0))
        self.assertEqual(remaining_segment2[1].start, volmdlr.Point2D(3.0, 4.0))
        self.assertEqual(remaining_segment2[1].end, volmdlr.Point2D(3.0, 6.0))
        self.assertEqual(self.lineseg5.delete_shared_section(self.lineseg6)[0], self.lineseg5)
        self.assertEqual(self.lineseg7.delete_shared_section(self.lineseg8)[0], self.lineseg7)
        lineseg_9 = edges.LineSegment2D(volmdlr.Point2D(6, 0), volmdlr.Point2D(6, 5))
        lineseg_10 = edges.LineSegment2D(volmdlr.Point2D(6, 0), volmdlr.Point2D(6, 5))
        self.assertFalse(lineseg_9.delete_shared_section(lineseg_10))


if __name__ == '__main__':
    unittest.main()
