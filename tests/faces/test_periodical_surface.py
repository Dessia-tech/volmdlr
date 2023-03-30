import unittest

import volmdlr
import volmdlr.wires as vmw
import volmdlr.edges as vme
import volmdlr.faces as vmf


class TestPeriodicalSurface(unittest.TestCase):
    def test_face_from_contours3d(self):
        surface = vmf.CylindricalSurface3D.load_from_file("faces/objects_periodical_surface/surface_openned_one_contour.json")
        contour3d_0 = vmw.Contour3D.load_from_file("faces/objects_periodical_surface/contour3d__openned_one_contour_0.json")
        contour3d_1 = vmw.Contour3D.load_from_file("faces/objects_periodical_surface/contour3d__openned_one_contour_1.json")

        contours = [contour3d_0, contour3d_1]
        face = surface.face_from_contours3d(contours)
        self.assertAlmostEqual(face.surface2d.area(), 0.4272566008882119, 1e-6)
        self.assertTrue(face.surface2d.outer_contour.is_ordered())

    def test_bsplinecurve3d_to_2d(self):
        surface = vmf.CylindricalSurface3D.load_from_file(
            "faces/objects_periodical_surface/periodicalsurface_with_theta_discontinuity.json")
        bspline = vme.BSplineCurve3D.load_from_file(
            "faces/objects_periodical_surface/bsplinecurve_with_theta_discontinuity.json")
        bspline2d = surface.bsplinecurve3d_to_2d(bspline)[0]
        theta1 = bspline2d.start.x
        theta2 = bspline2d.end.x
        self.assertEqual(theta1,  0.9979944870045463)
        self.assertEqual(theta2, volmdlr.TWO_PI)

if __name__ == '__main__':
    unittest.main()
