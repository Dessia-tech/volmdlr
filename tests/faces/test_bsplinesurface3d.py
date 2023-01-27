"""
Unit tests for volmdlr.faces.BSplineSurface3D
"""
import unittest

import volmdlr.edges as vme
import volmdlr.faces as vmf
import volmdlr.grid
from volmdlr.models import bspline_surfaces


class TestBSplineSurface3D(unittest.TestCase):
    def test_contour2d_parametric_to_dimension(self):
        bspline_face = bspline_surfaces.bspline_surface_2.rectangular_cut(0, 1, 0, 1)
        contour2d = bspline_surfaces.bspline_surface_2.contour3d_to_2d(bspline_face.outer_contour3d)
        grid2d = volmdlr.grid.Grid2D.from_properties((0, 1), (0, 1), (10, 10))
        contour2d_dim = bspline_surfaces.bspline_surface_2.contour2d_parametric_to_dimension(contour2d, grid2d)
        self.assertEqual(len(contour2d_dim.primitives), 4)
        self.assertAlmostEqual(contour2d_dim.area(), 18.15574028200993, places=2)
        self.assertAlmostEqual(contour2d_dim.length(), 16.81606170335965, places=2)

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

    def test_arc3d_to_2d(self):
        bspline_surface = vmf.BSplineSurface3D.load_from_file('faces/BSplineSurface3D_with_Arc3D.json')
        arc = vme.Arc3D(volmdlr.Point3D(0.01, 0.018, 0.014),
                        volmdlr.Point3D(0.00970710678118655, 0.018, 0.014707106781186547),
                        volmdlr.Point3D(0.009, 0.018, 0.015))

        test = bspline_surface.arc3d_to_2d(arc3d=arc)[0]

        inv_prof = bspline_surface.linesegment2d_to_3d(test)[0]

        # Verifies the inversion operation
        self.assertIsInstance(inv_prof, vme.Arc3D)
        self.assertEqual(inv_prof.start, arc.start)
        self.assertEqual(inv_prof.end, arc.end)

    def test_bsplinecurve3d_to_2d(self):
        bspline_surface = bspline_surfaces.bspline_surface_4
        control_points = [volmdlr.Point3D(-0.012138106431296442, 0.11769707710908962, -0.10360094389690414),
         volmdlr.Point3D(-0.012153195391844274, 0.1177764571887428, -0.10360691055433219),
         volmdlr.Point3D(-0.01216612946601426, 0.11785649353385147, -0.10361063821784446),
         volmdlr.Point3D(-0.012176888504086755, 0.11793706145749239, -0.10361212108019317)]
        weights = [1.0, 0.9994807070752826, 0.9994807070752826, 1.0]
        original_bspline = vme.BSplineCurve3D(3, control_points, [4, 4], [0, 1], weights, False)
        bspline_on_parametric_domain = bspline_surface.bsplinecurve3d_to_2d(original_bspline)[0]
        bspline_after_transfomation = bspline_surface.linesegment2d_to_3d(bspline_on_parametric_domain)[0]
        original_length = original_bspline.length()
        length_after_transformation = bspline_after_transfomation.length()
        point = original_bspline.point_at_abscissa(0.5 * original_length)
        point_test = bspline_after_transfomation.point_at_abscissa(0.5 * length_after_transformation)
        self.assertAlmostEqual(original_length, length_after_transformation, places=6)
        # self.assertTrue(point.is_close(point_test, 1e-6))


if __name__ == '__main__':
    unittest.main(verbosity=0)
