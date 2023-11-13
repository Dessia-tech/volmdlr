import os
import math
import os
import numpy as npy
import unittest
import volmdlr
from volmdlr import faces, surfaces, wires


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_toroidal_tests')


class TestToroidalFace3D(unittest.TestCase):
    def test_from_contours3d(self):
        surface = surfaces.ToroidalSurface3D.load_from_file(os.path.join(folder, "surface_4.json"))
        contour = wires.Contour3D.load_from_file(os.path.join(folder, "contour_4_0.json"))
        face = faces.ToroidalFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(face.surface2d.area(), 0.07116351378250674, 4)

        surface = surfaces.ToroidalSurface3D.load_from_file(
            os.path.join(folder, "repair_periodicity_toroidal_surface.json"))
        contour = wires.Contour3D.load_from_file(
            os.path.join(folder, "repair_periodicity_toroidal_surface_contour.json"))
        face = faces.ToroidalFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(face.surface2d.area(), math.pi**2, 4)
        self.assertTrue(face.surface2d.outer_contour.is_ordered())

    def test_planeface_intersections(self):
        expected_results = [[14.700000000000001], [9.388571408528668], [9.282044462349344], [9.107655321906883],
                            [8.870824383803773], [8.58245537985896], [5.0000000000008145, 5.000000000000816],
                            [3.717538057019154, 3.717538040295004], [3.325530330504112, 3.325530342893853],
                            [3.0819608468437045, 3.0819608470021023]]
        ts = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        tf = faces.ToroidalFace3D.from_surface_rectangular_cut(ts, -1.4, 3.5, 0., 2.5)

        list_expected_lenghts1 = []
        plane1 = surfaces.Plane3D(volmdlr.OXYZ)
        plane1 = plane1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi / 4)
        for i, n in enumerate(npy.linspace(0, math.pi / 4, 10)):
            plane = plane1.rotation(plane1.frame.origin, volmdlr.X3D, n)
            plane_face = faces.PlaneFace3D.from_surface_rectangular_cut(plane, 4, -4, 4, -4)
            plane_intersections = tf.face_intersections(plane_face)
            list_expected_lenghts1.append([i.length() for i in plane_intersections])
            for result, expected_result in zip(plane_intersections, expected_results[i]):
                self.assertAlmostEqual(result.length(), expected_result)

    def test_cylindricalface_intersections(self):
        expected_results = [[2.546120994711518], [2.454558505161535], [2.7679469885415657], [2.8109172462675667],
                            [1.3806998364554988, 3.028332404171969], [2.1248783089966574], [1.7368469566130957],
                            [2.558338114997606], [2.812361380094013, 1.3899450007345244], [2.4475153123576954]]
        toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        tf = faces.ToroidalFace3D.from_surface_rectangular_cut(toroidal_surface, 0, 3, 1, 3)
        frame = volmdlr.OXYZ.translation(volmdlr.Vector3D(1, 1, 0))
        for i, theta in enumerate(npy.linspace(0, math.pi * .7, 10)):
            frame = frame.rotation(frame.origin, volmdlr.Y3D, theta)
            cylindrical_surface = surfaces.CylindricalSurface3D(frame, 1.5)
            cylface = faces.CylindricalFace3D.from_surface_rectangular_cut(cylindrical_surface, 0, 4, -4, 4)
            inters = tf.face_intersections(cylface)
            for inter, expected_result in zip(inters, expected_results[i]):
                self.assertAlmostEqual(inter.length(), expected_result)


if __name__ == '__main__':
    unittest.main()
