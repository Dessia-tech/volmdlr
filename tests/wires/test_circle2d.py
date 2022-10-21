import unittest
import volmdlr.wires as vmw
import volmdlr.edges

circle = vmw.Circle2D(volmdlr.O2D, 0.50)
line = volmdlr.edges.Line2D(volmdlr.O2D, volmdlr.Point2D(2, 2))


class TestCircle2D(unittest.TestCase):

    def test_cut_by_line(self):
        self.assertTrue(circle.cut_by_line(line=line))


if __name__ == '__main__':
    unittest.main()
