import unittest
import math
from volmdlr.faces import SphericalFace3D
from volmdlr import surfaces
from dessia_common.core import DessiaObject


class TestSphericalFace3D(unittest.TestCase):
    def test_from_contours3d_and_rectangular_cut(self):
        surface = surfaces.SphericalSurface3D.load_from_file(
            "faces/objects_spherical_test/face_from_contours3d_and_rectangular_cut_surface.json")
        point = DessiaObject.load_from_file(
            "faces/objects_spherical_test/face_from_contours3d_and_rectangular_cut_point.json")
        contour_0 = DessiaObject.load_from_file(
            "faces/objects_spherical_test/face_from_contours3d_and_rectangular_cut_contour_0.json")
        contour_1 = DessiaObject.load_from_file(
            "faces/objects_spherical_test/face_from_contours3d_and_rectangular_cut_contour_1.json")
        face = SphericalFace3D.from_contours3d_and_rectangular_cut(surface, [contour_0, contour_1], point)
        self.assertEqual(len(face.surface2d.inner_contours), 2)
        self.assertAlmostEqual(face.surface2d.outer_contour.area(), math.pi**2 * 2)


if __name__ == '__main__':
    unittest.main()
