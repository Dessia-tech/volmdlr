import unittest

import volmdlr
from volmdlr.edges import Arc2D


class TestArc2D(unittest.TestCase):
    arc2d = Arc2D(volmdlr.Point2D(-1, 0), volmdlr.Point2D(0, -1), volmdlr.Point2D(1, 0))

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
