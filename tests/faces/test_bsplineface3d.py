import unittest
import os
from time import perf_counter
import volmdlr
from volmdlr import edges, surfaces, wires, faces, core
from volmdlr.models.bspline_surfaces import bspline_surface_1


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_bspline_test')


class TestBSplineFace3D(unittest.TestCase):
    bspline_face = faces.BSplineFace3D.from_surface_rectangular_cut(bspline_surface_1, 0, 1, 0, 1)

    def test_bounding_box(self):
        face = faces.BSplineFace3D.load_from_file(os.path.join(folder, "bsplineface_bbox_test.json"))
        bbox = face.bounding_box
        self.assertAlmostEqual(bbox.volume(), 0.00018, 5)

    def test_is_linesegment_crossing(self):
        linesegment = edges.LineSegment3D(volmdlr.Point3D(4, 0, 0), volmdlr.Point3D(4, 2, 2))
        self.assertTrue(self.bspline_face.is_linesegment_crossing(linesegment=linesegment))

    def test_linesegment_intersections_approximation(self):
        bsplineface = faces.BSplineFace3D.load_from_file(os.path.join(folder, 'bspline_face1.json'))
        lineseg = edges.LineSegment3D(volmdlr.Point3D(0, 0, 00.0015), volmdlr.Point3D(0, 0.005, 0.0025))
        intersections = bsplineface.linesegment_intersections(lineseg)
        self.assertEqual(len(intersections), 1)
        self.assertTrue(intersections[0], volmdlr.Point3D(0.0, 0.002350000000000002, 0.0019700000000000004))

    def test_from_contours3d(self):
        surface = surfaces.BSplineSurface3D.load_from_file(
            os.path.join(folder, "bspline_surface_openned_contour.json"))
        contour3d_0 = wires.Contour3D.load_from_file(os.path.join(folder, "bspline_contour_0_openned_contour.json"))
        contour3d_1 = wires.Contour3D.load_from_file(os.path.join(folder, "bspline_contour_1_openned_contour.json"))
        contours = [contour3d_0, contour3d_1]
        face = faces.BSplineFace3D.from_contours3d(surface, contours)
        self.assertTrue(face.surface2d.outer_contour.is_ordered(1e-5))
        self.assertAlmostEqual(face.surface2d.area(), 0.6257540398385436, 3)

        surface = surfaces.BSplineSurface3D.load_from_file(
            os.path.join(folder, "bsplineface_closedsurface_surface.json"))
        contour3d = wires.Contour3D.load_from_file(os.path.join(folder, "bsplineface_closedsurface_contour.json"))
        face = faces.BSplineFace3D.from_contours3d(surface, [contour3d])
        self.assertTrue(face.surface2d.outer_contour.is_ordered())
        self.assertAlmostEqual(face.surface2d.area(), 0.9996, 2)

        surface = surfaces.BSplineSurface3D.load_from_file(
            os.path.join(folder, "bsplineface_closedsurface_surface_2.json"))
        contour3d = wires.Contour3D.load_from_file(os.path.join(folder, "bsplineface_closedsurface_contour_2.json"))
        face = faces.BSplineFace3D.from_contours3d(surface, [contour3d])
        self.assertTrue(face.surface2d.outer_contour.is_ordered(1e-4))
        self.assertAlmostEqual(face.surface2d.area(), 1.0, 2)

        surface = surfaces.BSplineSurface3D.load_from_file(
            os.path.join(folder, "bsplinesurface_bsplineface_with_openned_contour.json"))
        contours3d = core.VolumeModel.load_from_file(
            os.path.join(folder, "bsplinesurface_bsplineface_with_openned_contour_contours.json")).primitives
        face = faces.BSplineFace3D.from_contours3d(surface, contours3d)
        self.assertAlmostEqual(face.surface2d.area(), 0.4261703133157918, 2)

        surface = surfaces.BSplineSurface3D.load_from_file(
            os.path.join(folder, "bsplineface_periodical_spiral_surface.json"))
        contour3d = wires.Contour3D.load_from_file(os.path.join(folder, "bsplineface_periodical_spiral_contour.json"))
        face = faces.BSplineFace3D.from_contours3d(surface, [contour3d])
        self.assertTrue(face.surface2d.outer_contour.is_ordered(1e-3))
        self.assertAlmostEqual(face.surface2d.area(), 0.49941, 2)

        surface = surfaces.BSplineSurface3D.load_from_file(
            os.path.join(folder, "slender_face_surface.json"))
        contour3d = wires.Contour3D.load_from_file(os.path.join(folder, "slender_face_contour.json"))
        face = faces.BSplineFace3D.from_contours3d(surface, [contour3d])
        self.assertTrue(face.surface2d.outer_contour.is_ordered(1e-6))
        self.assertAlmostEqual(face.surface2d.area(), 1.0, 2)

    def test_neutral_fiber(self):
        face = faces.BSplineFace3D.load_from_file(os.path.join(folder, "test_neutral_fiber.json"))
        neutral_fiber = face.neutral_fiber()
        self.assertAlmostEqual(neutral_fiber.length(), 0.030801389245691566, 2)

        face = faces.BSplineFace3D.load_from_file(
            os.path.join(folder, "test_neutral_fiber_2.json"))
        neutral_fiber = face.neutral_fiber()
        self.assertAlmostEqual(neutral_fiber.length(), 0.5327006535550406, 2)

        face = faces.BSplineFace3D.load_from_file(os.path.join(folder, "test_neutral_fiber_3.json"))
        neutral_fiber = face.neutral_fiber()
        self.assertAlmostEqual(neutral_fiber.length(), 0.1464370131293568, 2)

    def test_triangulation(self):
        surface = surfaces.BSplineSurface3D.load_from_file(os.path.join(folder, "spiral_bsplineface_surface.json"))
        contour3d = wires.Contour3D.load_from_file(os.path.join(folder, "spiral_bsplineface_contour.json"))
        start = perf_counter()
        face = faces.BSplineFace3D.from_contours3d(surface, [contour3d])
        end = perf_counter()
        total_time = end - start
        mesh = face.triangulation()
        self.assertAlmostEqual(face.surface2d.area(), 1, 2)
        self.assertGreaterEqual(len(mesh.vertices), 650)
        self.assertLessEqual(len(mesh.vertices), 2750)
        self.assertLessEqual(total_time, 2)


if __name__ == '__main__':
    unittest.main()
