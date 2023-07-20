import unittest

import volmdlr
from volmdlr import edges, surfaces, wires, faces
from volmdlr.models.bspline_surfaces import bspline_surface_1


class TestBSplineFace3D(unittest.TestCase):
    bspline_face = faces.BSplineFace3D.from_surface_rectangular_cut(bspline_surface_1, 0, 1, 0, 1)

    def test_bounding_box(self):
        face = faces.BSplineFace3D.load_from_file("faces/objects_bspline_test/bsplineface_bbox_test.json")
        bbox = face.bounding_box
        self.assertAlmostEqual(bbox.volume(), 0.00018, 5)

    def test_is_linesegment_crossing(self):
        linesegment = edges.LineSegment3D(volmdlr.Point3D(4, 0, 0), volmdlr.Point3D(4, 2, 2))
        self.assertTrue(self.bspline_face.is_linesegment_crossing(linesegment=linesegment))

    def test_linesegment_intersections_approximation(self):
        bsplineface = faces.BSplineFace3D.load_from_file('faces/objects_bspline_test/bspline_face1.json')
        lineseg = edges.LineSegment3D(volmdlr.Point3D(0, 0, 00.0015), volmdlr.Point3D(0, 0.005, 0.0025))
        intersections = bsplineface.linesegment_intersections(lineseg)
        self.assertEqual(len(intersections), 1)
        self.assertTrue(intersections[0], volmdlr.Point3D(0.0, 0.002350000000000002, 0.0019700000000000004))


    def test_from_contours3d(self):
        surface = surfaces.BSplineSurface3D.load_from_file(
            "faces/objects_bspline_test/bspline_surface_openned_contour.json")
        contour3d_0 = wires.Contour3D.load_from_file(
            "faces/objects_bspline_test/bspline_contour_0_openned_contour.json")
        contour3d_1 = wires.Contour3D.load_from_file(
            "faces/objects_bspline_test/bspline_contour_1_openned_contour.json")
        contours = [contour3d_0, contour3d_1]
        face = faces.BSplineFace3D.from_contours3d(surface, contours)
        self.assertAlmostEqual(face.surface2d.area(), 0.6324561566366691, 3)

    def test_neutral_fiber(self):
        face = faces.BSplineFace3D.load_from_file("faces/objects_bspline_test/test_neutral_fiber.json")
        neutral_fiber = face.neutral_fiber()
        self.assertAlmostEqual(neutral_fiber.length(), 0.030801389245691566, 2)

        face = faces.BSplineFace3D.load_from_file("faces/objects_bspline_test/test_neutral_fiber_2.json")
        neutral_fiber = face.neutral_fiber()
        self.assertAlmostEqual(neutral_fiber.length(), 0.5327006535550406, 2)

    def test_triangulation(self):
        surface = surfaces.BSplineSurface3D.load_from_file(
            "faces/objects_bspline_test/spiral_bsplineface_surface.json")
        contour3d = wires.Contour3D.load_from_file(
            "faces/objects_bspline_test/spiral_bsplineface_contour.json")
        face = faces.BSplineFace3D.from_contours3d(surface, [contour3d])
        mesh = face.triangulation()
        self.assertAlmostEqual(face.surface2d.area(), 1, 2)
        self.assertGreaterEqual(len(mesh.points), 650)
        self.assertLessEqual(len(mesh.points), 1300)


if __name__ == '__main__':
    unittest.main()
