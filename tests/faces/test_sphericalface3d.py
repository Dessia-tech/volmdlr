import unittest
import math
import os
from volmdlr.faces import SphericalFace3D
from volmdlr import surfaces
from dessia_common.core import DessiaObject


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_spherical_test')


class TestSphericalFace3D(unittest.TestCase):
    def test_from_contours3d_and_rectangular_cut(self):
        surface = surfaces.SphericalSurface3D.load_from_file(
            os.path.join(folder, "face_from_contours3d_and_rectangular_cut_surface.json"))
        point = DessiaObject.load_from_file(
            os.path.join(folder, "face_from_contours3d_and_rectangular_cut_point.json"))
        contour_0 = DessiaObject.load_from_file(
            os.path.join(folder, "face_from_contours3d_and_rectangular_cut_contour_0.json"))
        contour_1 = DessiaObject.load_from_file(
            os.path.join(folder, "face_from_contours3d_and_rectangular_cut_contour_1.json"))
        face = SphericalFace3D.from_contours3d_and_rectangular_cut(surface, [contour_0, contour_1], point)
        self.assertEqual(len(face.surface2d.inner_contours), 2)
        self.assertAlmostEqual(face.surface2d.outer_contour.area(), math.pi**2 * 2)

    def test_from_contours3d(self):
        surface = surfaces.SphericalSurface3D.load_from_file(
            os.path.join(folder, "sphericalface_disconnected_contours_surface.json"))
        contour_0 = DessiaObject.load_from_file(
            os.path.join(folder, "sphericalface_disconnected_contours_contour_0.json"))
        contour_1 = DessiaObject.load_from_file(
            os.path.join(folder, "sphericalface_disconnected_contours_contour_1.json"))
        face = SphericalFace3D.from_contours3d(surface, [contour_0, contour_1])
        self.assertAlmostEqual(face.surface2d.area(), 3.7274576655120804, 2)

        surface = surfaces.SphericalSurface3D.load_from_file(
            os.path.join(folder, "sphericalsurface_contour3d_to_2d_bug.json"))
        contour = DessiaObject.load_from_file(
            os.path.join(folder, "sphericalsurface_contour3d_to_2d_bug_contour.json"))
        face = SphericalFace3D.from_contours3d(surface, [contour])
        self.assertTrue(face.triangulation())

if __name__ == '__main__':
    unittest.main()
