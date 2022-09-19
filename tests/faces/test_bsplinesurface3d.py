"""
Unit tests for volmdlr.faces.BSplineSurface3D
"""
import unittest
from volmdlr.models import bspline_surfaces
import volmdlr.grid


class TestBSplineSurface3D(unittest.TestCase):
    def test_contour2d_parametric_to_dimension(self):
        bspline_face = bspline_surfaces.bspline_surface_2.rectangular_cut(0, 1, 0, 1)
        contour2d = bspline_surfaces.bspline_surface_2.contour3d_to_2d(bspline_face.outer_contour3d)
        grid2d = volmdlr.grid.Grid2D.from_properties((0, 1), (0, 1), (10, 10))
        contour2d_dim = bspline_surfaces.bspline_surface_2.contour2d_parametric_to_dimension(contour2d, grid2d)
        self.assertEqual(len(contour2d_dim.primitives), 4)
        self.assertAlmostEqual(contour2d_dim.area(), 18.12438798529036, places=5)
        self.assertAlmostEqual(contour2d_dim.length(), 16.816547325087043)


if __name__ == '__main__':
    unittest.main()
