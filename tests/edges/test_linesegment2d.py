"""
Unit tests for volmdlr.edges.LineSegment2D
"""
import unittest
import volmdlr
import volmdlr.edges

class TestLineSegment2D(unittest.TestCase):

    def test_to_wire(self):
        
        lns2d = volmdlr.edges.LineSegment2D(volmdlr.Point2D(0.2, 0.1),
                                            volmdlr.Point2D(0.5, 0.8))
        wire = lns2d.to_wire(2)

        self.assertAlmostEqual(len(wire.primitives), 2)


if __name__ == '__main__':
    unittest.main()
