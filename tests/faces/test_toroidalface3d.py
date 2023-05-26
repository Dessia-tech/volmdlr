import math
import unittest
from volmdlr import faces, surfaces, wires


class TestToroidalFace3D(unittest.TestCase):
    def test_from_contours3d(self):
        surface = surfaces.ToroidalSurface3D.load_from_file("faces/objects_toroidal_tests/surface_4.json")
        contour = wires.Contour3D.load_from_file("faces/objects_toroidal_tests/contour_4_0.json")
        face = faces.ToroidalFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(face.surface2d.area(), 0.07116351378250674, 4)

        surface = surfaces.ToroidalSurface3D.load_from_file(
            "faces/objects_toroidal_tests/repair_periodicity_toroidal_surface.json")
        contour = wires.Contour3D.load_from_file(
            "faces/objects_toroidal_tests/repair_periodicity_toroidal_surface_contour.json")
        face = faces.ToroidalFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(face.surface2d.area(), math.pi**2, 4)
        self.assertTrue(face.surface2d.outer_contour.is_ordered())


if __name__ == '__main__':
    unittest.main()
