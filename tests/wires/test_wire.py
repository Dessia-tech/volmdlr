"""
Unit tests for volmdlr.wires.Wire
"""
import unittest

import volmdlr
from volmdlr import edges, wires
from volmdlr.models.wires import wire3d_all_edges


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
            #TODO: uncomment this code when wires methods are refactored to share vertices between edges
            # for prim1, prim2 in zip(wire.primitives[:-1], wire.primitives[1:]):
            #     self.assertIs(prim1.end, prim2.start)

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

    def test_translation(self):
        new_wire = wire3d_all_edges.translation(volmdlr.Z3D)
        for prim1, prim2 in zip(new_wire.primitives[:-1], new_wire.primitives[1:]):
            self.assertIs(prim1.end, prim2.start)


if __name__ == '__main__':
    unittest.main()
