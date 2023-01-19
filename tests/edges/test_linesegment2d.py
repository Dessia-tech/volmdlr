"""
Unit tests for volmdlr.edges.LineSegment2D
"""
import unittest

import volmdlr
from volmdlr import edges


class TestLineSegment2D(unittest.TestCase):

    def test_to_wire(self):
        lns2d = edges.LineSegment2D(volmdlr.Point2D(0.2, 0.1),
                                    volmdlr.Point2D(0.5, 0.8))
        wire = lns2d.to_wire(2)

        self.assertAlmostEqual(len(wire.primitives), 2)

    def test_line_intersections(self):
        linesegment2d = edges.LineSegment2D(volmdlr.Point2D(1, 2), volmdlr.Point2D(3, -4))
        line1 = edges.Line2D(volmdlr.Point2D(1, 2), volmdlr.Point2D(3, 4))
        line2 = edges.Line2D(volmdlr.Point2D(3, -4), volmdlr.Point2D(-4, -3))
        line3 = edges.Line2D(volmdlr.Point2D(0, 0), volmdlr.Point2D(2, 2))
        self.assertEqual(linesegment2d.line_intersections(line1)[0], volmdlr.Point2D(1.0, 2.0))
        self.assertEqual(linesegment2d.line_intersections(line2)[0], volmdlr.Point2D(3, -4))
        self.assertEqual(linesegment2d.line_intersections(line3)[0], volmdlr.Point2D(1.25, 1.25))


if __name__ == '__main__':
    unittest.main()
