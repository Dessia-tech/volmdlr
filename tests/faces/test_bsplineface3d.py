import unittest
from volmdlr import surfaces, wires, faces


class MyTestCase(unittest.TestCase):
    def test_from_contours3d(self):
        surface = surfaces.BSplineSurface3D.load_from_file(
            "faces/objects_bspline_test/bspline_surface_openned_contour.json")
        contour3d_0 = wires.Contour3D.load_from_file(
            "faces/objects_bspline_test/bspline_contour_0_openned_contour.json")
        contour3d_1 = wires.Contour3D.load_from_file(
            "faces/objects_bspline_test/bspline_contour_1_openned_contour.json")
        contours = [contour3d_0, contour3d_1]
        face = faces.BSplineFace3D.from_contours3d(surface, contours)
        self.assertAlmostEqual(face.surface2d.area(), 0.6319342194477546, 5)


if __name__ == '__main__':
    unittest.main()
