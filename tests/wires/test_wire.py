"""
Unit tests for volmdlr.wires.Wire
"""
import unittest

import volmdlr
from volmdlr import edges, wires


class TestWire(unittest.TestCase):

    def test_from_edge(self):
        lns2d = edges.LineSegment2D(volmdlr.Point2D(0.2, 0.1),
                                    volmdlr.Point2D(0.5, 0.8))
        wire = wires.Wire2D.from_edge(lns2d, 2)
        self.assertAlmostEqual(len(wire.primitives), 2)


    def test_wires_from_edges(self):
        primitives = [
            edges.LineSegment2D(volmdlr.Point2D(-0.2, -0.2),
                                volmdlr.Point2D(0.2, -0.2)),
            edges.LineSegment2D(volmdlr.Point2D(-0.1, 0.6),
                                volmdlr.Point2D(-0.1, 0.2)),
            edges.LineSegment2D(volmdlr.Point2D(0.2, -0.2),
                                volmdlr.Point2D(0.2, 0.2)),
            edges.LineSegment2D(volmdlr.Point2D(0.3, 0.2),
                                volmdlr.Point2D(0.3, 0.6)),
            edges.LineSegment2D(volmdlr.Point2D(-0.2, 0.2),
                                volmdlr.Point2D(-0.2, -0.2)),
            edges.LineSegment2D(volmdlr.Point2D(0.3, 0.6),
                                volmdlr.Point2D(-0.1, 0.6))
            ]

        list_wires = wires.Wire2D.wires_from_edges(primitives)
        self.assertAlmostEqual(len(list_wires), 2)
        for wire in list_wires:
            self.assertAlmostEqual(len(wire.primitives), 3)

        primitives_1 = [
            edges.LineSegment2D(volmdlr.Point2D(-0.2, 0.2),
                                volmdlr.Point2D(-0.2, -0.2)),
            edges.LineSegment2D(volmdlr.Point2D(-0.2, -0.2),
                                volmdlr.Point2D(0.2, -0.2)),
            edges.LineSegment2D(volmdlr.Point2D(0.2, -0.2),
                                volmdlr.Point2D(0.2, 0.2))]

        primitives_2 = [
            edges.LineSegment2D(volmdlr.Point2D(0.3, 0.2),
                                volmdlr.Point2D(0.3, 0.6)),
            edges.LineSegment2D(volmdlr.Point2D(0.3, 0.6),
                                volmdlr.Point2D(-0.1, 0.6)),
            edges.LineSegment2D(volmdlr.Point2D(-0.1, 0.6),
                                volmdlr.Point2D(-0.1, 0.2))]

        self.assertTrue(list_wires[0], primitives_1)
        self.assertTrue(list_wires[1], primitives_2)


if __name__ == '__main__':
    unittest.main()
