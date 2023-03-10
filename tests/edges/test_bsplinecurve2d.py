import unittest
import dessia_common.core as dc


class TestBSplineCurve2D(unittest.TestCase):
    def test_bounding_rectangle(self):
        contour = dc.DessiaObject.load_from_file("edges\bounding_box_contour.json")
        b_rec = contour.bounding_rectangle
        self.assertAlmostEqual(b_rec.area(), 0.46995118100796823, places=2)


if __name__ == '__main__':
    unittest.main()
