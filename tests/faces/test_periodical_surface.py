import unittest
import volmdlr.wires as vmw
import volmdlr.faces as vmf


class TestPeriodicalSurface(unittest.TestCase):
    def test_face_from_contours3d(self):
        surface = vmf.CylindricalSurface3D.load_from_file("faces/objects_periodical_surface/surface_openned_one_contour.json")
        contour3d_0 = vmw.Contour3D.load_from_file("faces/objects_periodical_surface/contour3d__openned_one_contour_0.json")
        contour3d_1 = vmw.Contour3D.load_from_file("faces/objects_periodical_surface/contour3d__openned_one_contour_1.json")

        contours = [contour3d_0, contour3d_1]
        face = surface.face_from_contours3d(contours)
        self.assertEqual(face.surface2d.area(), 0.4272566008882119)
        self.assertTrue(face.surface2d.outer_contour.is_ordered())


if __name__ == '__main__':
    unittest.main()
