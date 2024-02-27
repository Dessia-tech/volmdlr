import unittest
import os
import math
from dessia_common.core import DessiaObject
import volmdlr
from volmdlr import edges, faces, surfaces, wires, curves
from volmdlr.models import conical_surfaces

folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_conical_tests')


class TestConicalFace3D(unittest.TestCase):
    conical_surface = conical_surfaces.conical_surface1
    outer_contour2d = wires.Contour2D(
        [edges.LineSegment2D(volmdlr.Point2D(-math.pi, 0.0), volmdlr.Point2D(math.pi, 0.0)),
         edges.LineSegment2D(volmdlr.Point2D(math.pi, 0.0), volmdlr.Point2D(math.pi, 1.0)),
         edges.LineSegment2D(volmdlr.Point2D(math.pi, 1.0), volmdlr.Point2D(-math.pi, 1.0)),
         edges.LineSegment2D(volmdlr.Point2D(-math.pi, 1.0), volmdlr.Point2D(-math.pi, 0.0))])
    inner_contour = wires.Contour2D(
        [edges.LineSegment2D(volmdlr.Point2D(-0.5 * math.pi, 0.4), volmdlr.Point2D(0.5 * math.pi, 0.4)),
         edges.LineSegment2D(volmdlr.Point2D(0.5 * math.pi, 0.4), volmdlr.Point2D(0.75 * math.pi, 0.6)),
         edges.LineSegment2D(volmdlr.Point2D(0.75 * math.pi, 0.6), volmdlr.Point2D(-0.5 * math.pi, 0.7)),
         edges.LineSegment2D(volmdlr.Point2D(-0.5 * math.pi, 0.7), volmdlr.Point2D(-0.5 * math.pi, 0.4))])
    surface2d = surfaces.Surface2D(outer_contour2d, [inner_contour])
    conical_face = faces.ConicalFace3D(conical_surface, surface2d)

    def test_primitives_mapping(self):
        primitives_mapping = self.conical_face.primitives_mapping
        self.assertEqual(len(primitives_mapping), 7)
        expected_pimitives = ["LineSegment3D", "FullArc3D", "LineSegment3D", "Arc3D", "BSplineCurve3D",
                              "BSplineCurve3D", "LineSegment3D"]
        self.assertIsNone(primitives_mapping.get(self.outer_contour2d.primitives[0]))
        for prim, expected in zip(self.conical_face.surface2d.outer_contour.primitives[1:]
                                  + self.conical_face.surface2d.inner_contours[0].primitives,
                                  expected_pimitives):
            self.assertEqual(primitives_mapping.get(prim).__class__.__name__, expected)

    def test_get_face_polygons(self):
        outer_polygon, inner_polygons = self.conical_face.get_face_polygons()
        self.assertAlmostEqual(outer_polygon.area(), 2 * math.pi)
        self.assertAlmostEqual(inner_polygons[0].area(), 0.9032078879070156)

    def test_from_contours(self):
        buggy_conical_surface = DessiaObject.from_json(os.path.join(folder, "conical_surface1.json"))
        buggy_contours3d1 = DessiaObject.from_json(os.path.join(folder, 'face_from_contours1_0.json'))
        buggy_contours3d2 = DessiaObject.from_json(os.path.join(folder, 'face_from_contours1_1.json'))

        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface,
                                                           [buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.003769911184307754, 4)

        buggy_conical_surface = DessiaObject.from_json(os.path.join(folder, 'conical_surface3d_1.json'))
        buggy_contours3d1 = DessiaObject.from_json(os.path.join(folder, 'face_contour1.json'))
        buggy_contours3d2 = DessiaObject.from_json(os.path.join(folder, 'face_contour2.json'))

        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface,
                                                           [buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.0016000193084354127, 4)

        buggy_conical_surface = DessiaObject.from_json(os.path.join(folder, 'conical_surface3d_2.json'))
        buggy_contours3d1 = DessiaObject.from_json(os.path.join(folder, 'face_contour3_.json'))
        buggy_contours3d2 = DessiaObject.from_json(os.path.join(folder, 'face_contour4_.json'))

        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface,
                                                           [buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.055154411016251716, 4)

        buggy_conical_surface = surfaces.ConicalSurface3D.from_json(
            os.path.join(folder, "conical_surface_with_singularity.json"))
        buggy_contours3d = wires.Contour3D.from_json(
            os.path.join(folder, 'conical_contour_with_singularity.json'))
        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface, [buggy_contours3d])
        self.assertEqual(len(conical_face.surface2d.outer_contour.primitives), 5)
        self.assertAlmostEqual(conical_face.area(), 0.0009613769926732048 * volmdlr.TWO_PI, 4)

        buggy_conical_surface = DessiaObject.from_json(
            os.path.join(folder, 'conicalsurface_openned_contours.json'))
        buggy_contours3d1 = DessiaObject.from_json(
            os.path.join(folder, 'conicalsurface_openned_contours_contour_0.json'))
        buggy_contours3d2 = DessiaObject.from_json(
            os.path.join(folder, 'conicalsurface_openned_contours_contour_1.json'))

        conical_face = faces.ConicalFace3D.from_contours3d(buggy_conical_surface,
                                                           [buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(conical_face.surface2d.inner_contours)
        self.assertAlmostEqual(conical_face.area(), 0.021682359796019755 , 3)

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
        points = circle.discretization_points(number_points=21)
        for i, (theta1, theta2) in enumerate(theta_ranges):
            conical_face = faces.ConicalFace3D.from_surface_rectangular_cut(
                conical_surface, theta1, theta2, 0., 1)
            for j, expected_result in enumerate(expected_results_lists[i]):
                if expected_result == 1:
                    self.assertTrue(conical_face.point_belongs(points[j]))
                    continue
                self.assertFalse(conical_face.point_belongs(points[j]))

    def test_triangulation(self):
        face = faces.ConicalFace3D.from_json(
            os.path.join(folder, "conicalface_segfault_with_tri_opt_set_to_pq.json"))
        mesh2d = face.triangulation()
        self.assertIsNotNone(mesh2d)

    def test_normal_at_point(self):
        face, point = DessiaObject.from_json(os.path.join(folder,'test_conicalface_normal_at_point.json')).primitives
        normal = face.normal_at_point(point)
        self.assertTrue(normal.is_close(volmdlr.Vector3D(0.4736709994871299, 0.7250074373721028, -0.4999999999999999)))


if __name__ == '__main__':
    unittest.main()
