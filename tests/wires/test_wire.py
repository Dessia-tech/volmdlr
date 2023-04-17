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


if __name__ == '__main__':
    unittest.main()
