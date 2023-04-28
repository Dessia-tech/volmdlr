
import unittest

import volmdlr
from volmdlr import edges, faces
from volmdlr.models.bspline_surfaces import bspline_surface_1


class TestBSplineFace3D(unittest.TestCase):
    bspline_face = faces.BSplineFace3D.from_surface_rectangular_cut(bspline_surface_1, 0, 1, 0, 1)

    def test_is_linesegment_crossing(self):
        linesegment = edges.LineSegment3D(volmdlr.Point3D(4, 0, 0), volmdlr.Point3D(4, 2, 2))
        self.assertTrue(self.bspline_face.is_linesegment_crossing(linesegment=linesegment))


if __name__ == '__main__':
    unittest.main()
