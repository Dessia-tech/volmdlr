import unittest

import volmdlr
from volmdlr import edges


class TestFullArc3D(unittest.TestCase):
    fullarc3d = edges.FullArc3D(
        volmdlr.Point3D(0.0, 0.0, -0.25), volmdlr.Point3D(0.15, 0.0, -0.25), volmdlr.Vector3D(0.0, 0.0, 1.0)
    )

    def test_linesegment_intersections(self):
        line3d = edges.LineSegment3D(volmdlr.Point3D(-0.2, -0.2, -0.25), volmdlr.Point3D(0.2, 0.2, -0.25))
        line3d_ = edges.LineSegment3D(volmdlr.Point3D(-0.15, 0, 0), volmdlr.Point3D(-0.15, 0.0, -0.35))
        linesegment_intesections1 = self.fullarc3d.linesegment_intersections(line3d)
        self.assertEqual(len(linesegment_intesections1), 2)
        self.assertEqual(linesegment_intesections1[0], volmdlr.Point3D(0.10606601717798213, 0.10606601717798214, -0.25))
        self.assertEqual(
            linesegment_intesections1[1], volmdlr.Point3D(-0.10606601717798213, -0.10606601717798213, -0.25)
        )
        linesegment_intesections2 = self.fullarc3d.linesegment_intersections(line3d_)
        self.assertEqual(len(linesegment_intesections2), 1)
        self.assertEqual(linesegment_intesections2[0], volmdlr.Point3D(-0.15, 0.0, -0.25))


if __name__ == "__main__":
    unittest.main()
