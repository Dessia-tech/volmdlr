"""
Unit tests for volmdlr.edges.LineSegment2D
"""
import unittest

import volmdlr
from volmdlr import edges


class TestLineSegment2D(unittest.TestCase):

    def test_line_intersections(self):
        linesegment2d = edges.LineSegment2D(volmdlr.Point2D(1, 2), volmdlr.Point2D(3, -4))
        line1 = edges.Line2D(volmdlr.Point2D(1, 2), volmdlr.Point2D(3, 4))
        line2 = edges.Line2D(volmdlr.Point2D(3, -4), volmdlr.Point2D(-4, -3))
        line3 = edges.Line2D(volmdlr.Point2D(0, 0), volmdlr.Point2D(2, 2))
        self.assertEqual(linesegment2d.line_intersections(line1)[0], volmdlr.Point2D(1.0, 2.0))
        self.assertTrue(linesegment2d.line_intersections(line2)[0].is_close(volmdlr.Point2D(3, -4)))
        self.assertTrue(linesegment2d.line_intersections(line3)[0].is_close(volmdlr.Point2D(1.25, 1.25)))


if __name__ == '__main__':
    unittest.main()
