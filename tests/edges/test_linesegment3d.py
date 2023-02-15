import unittest

import volmdlr
from volmdlr import edges


class TestLineSegment3D(unittest.TestCase):
    linesegment1 = edges.LineSegment3D(
        volmdlr.Point3D(-0.2, -0.2, 0.2), volmdlr.Point3D(-0.2, -0.2, -0.2), name="linesegment1"
    )
    linesegment2 = edges.LineSegment3D(
        volmdlr.Point3D(-0.225, -0.2, 0.125), volmdlr.Point3D(-0.225, -0.2, 0.275), name="linesegment2"
    )
    linesegment3 = edges.LineSegment3D(
        volmdlr.Point3D(-0.225, -0.2, 0.275), volmdlr.Point3D(0.225, -0.2, 0.275), name="linesegment3"
    )
    linesegment4 = edges.LineSegment3D(
        volmdlr.Point3D(0.225, -0.2, 0.275), volmdlr.Point3D(0.225, -0.2, 0.125), name="linesegment4"
    )
    linesegment5 = edges.LineSegment3D(
        volmdlr.Point3D(0.225, -0.2, 0.125), volmdlr.Point3D(-0.225, -0.2, 0.125), name="linesegment5"
    )

    def test_line_intersection(self):
        self.assertFalse(self.linesegment1.line_intersections(self.linesegment3.to_line()))
        self.assertEqual(
            self.linesegment3.line_intersections(self.linesegment1.to_line()), [volmdlr.Point3D(-0.2, -0.2, 0.275)]
        )

    def test_linesegment_intersection(self):
        for lineseg in [self.linesegment2, self.linesegment3, self.linesegment4]:
            self.assertFalse(self.linesegment1.linesegment_intersection(lineseg))
        self.assertEqual(
            self.linesegment1.linesegment_intersection(self.linesegment5), volmdlr.Point3D(-0.2, -0.2, 0.125)
        )


if __name__ == "__main__":
    unittest.main()
