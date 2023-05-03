import unittest
import math
from dessia_common.core import DessiaObject
import volmdlr
from volmdlr import faces, surfaces, wires
from volmdlr.models import conical_surfaces


class TestConicalFace3D(unittest.TestCase):
    def test_from_contours(self):
        buggy_conical_surface = DessiaObject.load_from_file(
            'faces/objects_conical_tests/conical_surface1.json')
        buggy_contours3d1 = DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_from_contours1_0.json')
        buggy_contours3d2 = DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_from_contours1_1.json')

        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface,
                                                           [buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.003769911184307754, 4)

        buggy_conical_surface = DessiaObject.load_from_file(
            'faces/objects_conical_tests/conical_surface3d_1.json')
        buggy_contours3d1 = DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_contour1.json')
        buggy_contours3d2 = DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_contour2.json')

        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface,
                                                           [buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.0016000193084354127, 4)

        buggy_conical_surface = DessiaObject.load_from_file(
            'faces/objects_conical_tests/conical_surface3d_2.json')
        buggy_contours3d1 = DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_contour3_.json')
        buggy_contours3d2 = DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_contour4_.json')

        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface,
                                             [buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.055154411016251716, 4)

        buggy_conical_surface = surfaces.ConicalSurface3D.load_from_file(
            "faces/objects_conical_tests/conical_surface_with_singularity.json")
        buggy_contours3d = wires.Contour3D.load_from_file(
            'faces/objects_conical_tests/conical_contour_with_singularity.json')
        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface, [buggy_contours3d])
        self.assertEqual(len(conical_face.surface2d.outer_contour.primitives), 5)
        self.assertAlmostEqual(conical_face.area(), 0.0009613769926732048*volmdlr.TWO_PI, 4)

    def test_from_base_and_vertex(self):
        circle = wires.Circle3D(
            volmdlr.Frame3D(volmdlr.Point3D(0, 0, 1), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D), 0.5 * math.sqrt(3)
        )
        trim_point = volmdlr.Point3D(0.5 * math.sqrt(3), 0, 1)
        fullarc = circle.trim(trim_point, trim_point)
        contour = wires.Contour3D([fullarc])
        face = faces.ConicalFace3D.from_base_and_vertex(conical_surfaces.conical_surface1, contour, volmdlr.O3D)
        self.assertEqual(face.surface2d.area(), volmdlr.TWO_PI)

if __name__ == '__main__':
    unittest.main()
