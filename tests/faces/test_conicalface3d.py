import unittest
import os
import math
from dessia_common.core import DessiaObject
import volmdlr
from volmdlr import faces, surfaces, wires, curves
from volmdlr.models import conical_surfaces

folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_conical_tests')


class TestConicalFace3D(unittest.TestCase):
    def test_from_contours(self):
        buggy_conical_surface = DessiaObject.load_from_file(os.path.join(folder, "conical_surface1.json"))
        buggy_contours3d1 = DessiaObject.load_from_file(os.path.join(folder, 'face_from_contours1_0.json'))
        buggy_contours3d2 = DessiaObject.load_from_file(os.path.join(folder, 'face_from_contours1_1.json'))

        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface,
                                                           [buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.003769911184307754, 4)

        buggy_conical_surface = DessiaObject.load_from_file(os.path.join(folder, 'conical_surface3d_1.json'))
        buggy_contours3d1 = DessiaObject.load_from_file(os.path.join(folder, 'face_contour1.json'))
        buggy_contours3d2 = DessiaObject.load_from_file(os.path.join(folder, 'face_contour2.json'))

        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface,
                                                           [buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.0016000193084354127, 4)

        buggy_conical_surface = DessiaObject.load_from_file(os.path.join(folder, 'conical_surface3d_2.json'))
        buggy_contours3d1 = DessiaObject.load_from_file(os.path.join(folder, 'face_contour3_.json'))
        buggy_contours3d2 = DessiaObject.load_from_file(os.path.join(folder, 'face_contour4_.json'))

        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface,
                                                           [buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.055154411016251716, 4)

        buggy_conical_surface = surfaces.ConicalSurface3D.load_from_file(
            os.path.join(folder, "conical_surface_with_singularity.json"))
        buggy_contours3d = wires.Contour3D.load_from_file(
            os.path.join(folder, 'conical_contour_with_singularity.json'))
        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface, [buggy_contours3d])
        self.assertEqual(len(conical_face.surface2d.outer_contour.primitives), 5)
        self.assertAlmostEqual(conical_face.area(), 0.0009613769926732048 * volmdlr.TWO_PI, 4)

    def test_from_base_and_vertex(self):
        circle = curves.Circle3D(
            volmdlr.Frame3D(volmdlr.Point3D(0, 0, 1), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D), 0.5 * math.sqrt(3)
        )
        trim_point = volmdlr.Point3D(0.5 * math.sqrt(3), 0, 1)
        fullarc = circle.trim(trim_point, trim_point)
        contour = wires.Contour3D([fullarc])
        face = faces.ConicalFace3D.from_base_and_vertex(conical_surfaces.conical_surface1, contour, volmdlr.O3D)
        self.assertEqual(face.surface2d.area(), volmdlr.TWO_PI)

    def test_neutral_fiber(self):
        surface = conical_surfaces.conical_surface1
        face = faces.ConicalFace3D.from_surface_rectangular_cut(surface, 0, math.pi, 0.5, 1)
        neutral_fiber = face.neutral_fiber()
        self.assertEqual(neutral_fiber.length(), 0.5)

    def test_point_belongs(self):
        expected_results_lists = [[0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1],
                                  [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0],
                                  [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                                  [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]]
        conical_surface = surfaces.ConicalSurface3D(volmdlr.OXYZ, math.pi / 6)
        theta_ranges = [(.3 * math.pi, math.pi), (-.3 * math.pi, math.pi), (-1.5 * math.pi, -.3 * math.pi),
                        (0.0, math.pi), (-math.pi, 0.0)]
        circle_radius = 0.5 * math.tan(conical_surface.semi_angle)
        circle = curves.Circle3D(volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0.5),
                                                 volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D), radius=circle_radius)
        points = circle.discretization_points(number_points=20)
        for i, (theta1, theta2) in enumerate(theta_ranges):
            conical_face = faces.ConicalFace3D.from_surface_rectangular_cut(
                conical_surface, theta1, theta2, 0., 1)
            for j, expected_result in enumerate(expected_results_lists[i]):
                if expected_result == 1:
                    self.assertTrue(conical_face.point_belongs(points[j]))
                    continue
                self.assertFalse(conical_face.point_belongs(points[j]))


if __name__ == '__main__':
    unittest.main()
