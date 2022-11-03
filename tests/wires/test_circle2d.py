import unittest
import volmdlr.wires as vmw
import volmdlr.edges

circle = vmw.Circle2D(volmdlr.O2D, 0.50)
line = volmdlr.edges.Line2D(volmdlr.O2D, volmdlr.Point2D(0, 1))

arc1_validate = volmdlr.edges.Arc2D(volmdlr.Point2D(0, -0.5), volmdlr.Point2D(0.5, 0), volmdlr.Point2D(0, 0.5))
arc2_validate = volmdlr.edges.Arc2D(volmdlr.Point2D(0, -0.5), volmdlr.Point2D(-0.5, 0), volmdlr.Point2D(0, 0.5))


class TestCircle2D(unittest.TestCase):
    def test_split_by_line(self):
        list_arcs = circle.split_by_line(line=line)

        self.assertEqual(list_arcs[0].length(), arc1_validate.length())
        self.assertEqual(list_arcs[0].start, arc1_validate.start)
        self.assertEqual(list_arcs[0].end, arc1_validate.end)
        self.assertEqual(list_arcs[1].length(), arc2_validate.length())
        self.assertEqual(list_arcs[1].start, arc2_validate.start)
        self.assertEqual(list_arcs[1].end, arc2_validate.end)


if __name__ == '__main__':
    unittest.main()
