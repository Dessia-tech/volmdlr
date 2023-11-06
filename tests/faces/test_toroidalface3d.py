import math
import numpy as npy
import unittest
import volmdlr
from volmdlr import faces, surfaces, wires


class TestToroidalFace3D(unittest.TestCase):
    def test_from_contours3d(self):
        surface = surfaces.ToroidalSurface3D.load_from_file("faces/objects_toroidal_tests/surface_4.json")
        contour = wires.Contour3D.load_from_file("faces/objects_toroidal_tests/contour_4_0.json")
        face = faces.ToroidalFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(face.surface2d.area(), 0.07116351378250674, 4)

        surface = surfaces.ToroidalSurface3D.load_from_file(
            "faces/objects_toroidal_tests/repair_periodicity_toroidal_surface.json")
        contour = wires.Contour3D.load_from_file(
            "faces/objects_toroidal_tests/repair_periodicity_toroidal_surface_contour.json")
        face = faces.ToroidalFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(face.surface2d.area(), math.pi**2, 4)
        self.assertTrue(face.surface2d.outer_contour.is_ordered())

    def test_planeface_intersections(self):
        expected_results = [[14.700000000000001], [9.388571408528646], [9.282044462349342], [9.107655321473693],
                            [8.870824383893368], [8.582455379382449], [4.9999999999983755, 4.9999999999983755],
                            [3.7175380418398443, 3.717538025078572], [3.325530314317449, 3.3255303255139186],
                            [3.081960732019063, 3.0819607320979014]]
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
        expected_results = [[2.5461207980560485], [2.454557959936191], [2.767950454333099], [2.8109131068618605],
                            [1.3806998377604578, 3.0341016797202873], [2.1248777973309823], [1.7368914447612895],
                            [2.5583377804228014], [1.3899444850660725, 2.812800677753292], [2.4475236530517783]]
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
