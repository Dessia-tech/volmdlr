"""
Unit tests for volmdlr.faces.BSplineSurface3D
"""
import unittest
from volmdlr.models import bspline_surfaces
import volmdlr.grid
import volmdlr.faces as vmf


class TestBSplineSurface3D(unittest.TestCase):
    def test_contour2d_parametric_to_dimension(self):
        bspline_face = bspline_surfaces.bspline_surface_2.rectangular_cut(0, 1, 0, 1)
        contour2d = bspline_surfaces.bspline_surface_2.contour3d_to_2d(bspline_face.outer_contour3d)
        grid2d = volmdlr.grid.Grid2D.from_properties((0, 1), (0, 1), (10, 10))
        contour2d_dim = bspline_surfaces.bspline_surface_2.contour2d_parametric_to_dimension(contour2d, grid2d)
        self.assertEqual(len(contour2d_dim.primitives), 4)
        self.assertAlmostEqual(contour2d_dim.area(), 18.112582105805615, places=2)
        self.assertAlmostEqual(contour2d_dim.length(), 16.816547325087043, places=2)

    def test_periodicity(self):
        bspline_suface = vmf.BSplineSurface3D.load_from_file('faces/surface3d_8.json')
        self.assertAlmostEqual(bspline_suface.x_periodicity,  0.8888888888888888)
        self.assertFalse(bspline_suface.y_periodicity)

    def test_bbox(self):
        surface = bspline_surfaces.bspline_surface_3
        bbox = surface.bounding_box
        volume = bbox.volume()

        # Check if the bounding box volume is correct
        self.assertEqual(volume, 4.0)

if __name__ == '__main__':
    unittest.main(verbosity=0)
