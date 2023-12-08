import os
import math
import os
import numpy as npy
import unittest
from dessia_common.core import DessiaObject
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
        expected_results = [[14.700000000000001], [9.388571088159908], [9.282044462953378],
                            [9.10765519911981], [8.870824370954356], [8.582455381664765],
                            [4.999999999998194, 4.999999999998194], [3.717538102728644, 3.7240013105263103],
                            [3.325530350863936, 3.325501013465095], [3.0819608531081744, 3.094049612384317]]

        ts = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        tf = faces.ToroidalFace3D.from_surface_rectangular_cut(ts, -1.4, 3.5, 0., 2.5)

        list_expected_lenghts1 = []
        plane1 = surfaces.Plane3D(volmdlr.OXYZ)
        plane1 = plane1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi / 4)
        for i, n in enumerate(npy.linspace(0, math.pi / 4, 10)):
            plane = plane1.rotation(plane1.frame.origin, volmdlr.X3D, n)
            plane_face = faces.PlaneFace3D.from_surface_rectangular_cut(plane, 4, -4, 4, -4)
            planeface_intersections = tf.face_intersections(plane_face)
            list_expected_lenghts1.append([i.length() for i in planeface_intersections])
            self.assertEqual(len(planeface_intersections), len(expected_results[i]))
            for result, expected_result in zip(planeface_intersections, expected_results[i]):
                self.assertAlmostEqual(result.length(), expected_result, 6)

        planeface, toroidalface = DessiaObject.load_from_file(
            os.path.join(folder, "test_planeface_toroidialface_intersections301123.json")).primitives

        inters = planeface.face_intersections(toroidalface)
        self.assertEqual(len(inters), 1)
        self.assertAlmostEqual(inters[0].length(), 0.08139556829160953)

    def test_cylindricalface_intersections(self):
        expected_results = [[2.546120994711518], [2.454558505161535], [2.7679469885415657], [2.8109172462675667],
                            [1.3806998364554988, 3.028332404171969], [2.1248783089966574], [1.736847875568775],
                            [2.558338114997606], [2.812361380094013, 1.3899450007345244], [2.4475153123576954]]
        toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        tf = faces.ToroidalFace3D.from_surface_rectangular_cut(toroidal_surface, 0, 3, 1, 3)
        frame = volmdlr.OXYZ.translation(volmdlr.Vector3D(1, 1, 0))
        for i, theta in enumerate(npy.linspace(0, math.pi * .7, 10)):
            frame = frame.rotation(frame.origin, volmdlr.Y3D, theta)
            cylindrical_surface = surfaces.CylindricalSurface3D(frame, 1.5)
            cylface = faces.CylindricalFace3D.from_surface_rectangular_cut(cylindrical_surface, 0, 4, -4, 4)
            inters = tf.face_intersections(cylface)
            self.assertEqual(len(inters), len(expected_results[i]))
            for inter, expected_result in zip(inters, expected_results[i]):
                self.assertAlmostEqual(inter.length(), expected_result, 6)


if __name__ == '__main__':
    unittest.main()
