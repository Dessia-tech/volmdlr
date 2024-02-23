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
    surface1 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 0.32, 0.08)
    face1 = faces.ToroidalFace3D.from_surface_rectangular_cut(surface1, -0.1, 1.3, 2, 0.3)

    def test_from_contours3d(self):
        surface = surfaces.ToroidalSurface3D.from_json(os.path.join(folder, "surface_4.json"))
        contour = wires.Contour3D.from_json(os.path.join(folder, "contour_4_0.json"))
        face = faces.ToroidalFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(face.surface2d.area(), 0.07116351378250674, 4)

        surface = surfaces.ToroidalSurface3D.from_json(
            os.path.join(folder, "repair_periodicity_toroidal_surface.json"))
        contour = wires.Contour3D.from_json(
            os.path.join(folder, "repair_periodicity_toroidal_surface_contour.json"))
        face = faces.ToroidalFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(face.surface2d.area(), math.pi**2, 4)
        self.assertTrue(face.surface2d.outer_contour.is_ordered())

        surface = surfaces.ToroidalSurface3D.from_json(
            os.path.join(folder, "face_with_inner_contour_surface.json"))
        contour0 = wires.Contour3D.from_json(os.path.join(folder, "face_with_inner_contour_contour0.json"))
        contour1 = wires.Contour3D.from_json(os.path.join(folder, "face_with_inner_contour_contour1.json"))
        face = faces.ToroidalFace3D.from_contours3d(surface, [contour0, contour1])
        self.assertAlmostEqual(face.surface2d.area(), 33.03042743115413, 2)

        surface = surfaces.ToroidalSurface3D.from_json(
            os.path.join(folder, "repair_inner_contour_periodicity_surface.json"))
        contour0 = wires.Contour3D.from_json(
            os.path.join(folder, "repair_inner_contour_periodicity_contour_0.json"))
        contour1 = wires.Contour3D.from_json(
            os.path.join(folder, "repair_inner_contour_periodicity_contour_1.json"))
        face = faces.ToroidalFace3D.from_contours3d(surface, [contour0, contour1])
        self.assertAlmostEqual(face.surface2d.area(), 36.56961010698211, 2)

    def test_planeface_intersections(self):
        expected_results = [1, 1, 1, 1, 1, 1, 2, 2, 2, 2]

        ts = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        tf = faces.ToroidalFace3D.from_surface_rectangular_cut(ts, -1.4, 3.5, 0., 2.5)

        # list_expected_lenghts1 = []
        plane1 = surfaces.Plane3D(volmdlr.OXYZ)
        plane1 = plane1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi / 4)
        for i, n in enumerate(npy.linspace(0, math.pi / 4, 10)):
            plane = plane1.rotation(plane1.frame.origin, volmdlr.X3D, n)
            plane_face = faces.PlaneFace3D.from_surface_rectangular_cut(plane, 4, -4, 4, -4)
            planeface_intersections = tf.face_intersections(plane_face)
            # list_expected_lenghts1.append([i.length() for i in planeface_intersections])
            self.assertEqual(len(planeface_intersections), expected_results[i])
            self.assertTrue(all(tf.point_belongs(p, 1e-4) and plane_face.point_belongs(p, 1e-4)
                                for i in planeface_intersections for p in i.primitives[0].points))
            # for result, expected_result in zip(planeface_intersections, expected_results[i]):
            #     self.assertAlmostEqual(result.length(), expected_result, 5)

        planeface, toroidalface = DessiaObject.from_json(
            os.path.join(folder, "test_planeface_toroidialface_intersections301123.json")).primitives

        inters = planeface.face_intersections(toroidalface)
        self.assertEqual(len(inters), 1)
        self.assertAlmostEqual(inters[0].length(), 0.08139556829160953, 5)

        planeface, toroidalface = DessiaObject.from_json(
            os.path.join(folder, 'test_planeface3d_toroidalface3d_121223.json')).primitives
        intersections = planeface.face_intersections(toroidalface)
        self.assertEqual(len(intersections), 1)
        self.assertAlmostEqual(intersections[0].length(), 0.0033804467442557404, 5)

        planeface, toroidalface = DessiaObject.from_json(
            os.path.join(folder, "test_planeface3d_toroidalface3d_131223.json")).primitives

        inters = planeface.face_intersections(toroidalface)
        self.assertEqual(len(inters), 1)
        self.assertAlmostEqual(inters[0].length(), 0.030296492908080553, 5)

    def test_cylindricalface_intersections(self):
        expected_results = [[2.5461209954222026], [2.454561591082158], [2.7679468571575105], [2.8109179729321183],
                            [1.3806998569480715, 3.0283316710422508], [2.1248782869459646], [1.7368478889595058],
                            [2.55833794579346], [2.8123613465408064, 1.3899450251554277], [2.447515630586587]]
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

    def test_conicalface_intersections(self):
        tf = faces.ToroidalFace3D.from_surface_rectangular_cut(
            surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1), -1.4,
            3.5, 0., 2.5)
        conical_surface = surfaces.ConicalSurface3D(volmdlr.OXYZ.translation(volmdlr.Vector3D(-1, 0, 0)), math.pi / 6)
        conical_face = faces.ConicalFace3D.from_surface_rectangular_cut(
            conical_surface, 0, volmdlr.TWO_PI, 0, 2)
        face_intersections = tf.face_intersections(conical_face)
        self.assertEqual(len(face_intersections), 1)
        for point in face_intersections[0].discretization_points(number_points=20):
            self.assertTrue(tf.point_belongs(point))
            self.assertTrue(conical_face.point_belongs(point))

    def test_triangulation_quality(self):
        """
        The triangle middle of triangulation should be at least at radius/20 of the surface
        """
        triangulation = self.face1.triangulation()
        for i1, i2, i3 in triangulation.triangles:
            point1 = volmdlr.Point3D(*triangulation.vertices[i1])
            point2 = volmdlr.Point3D(*triangulation.vertices[i2])
            point3 = volmdlr.Point3D(*triangulation.vertices[i3])

            triangle = faces.Triangle3D(point1, point2, point3)
            # Test orthogonality
            self.assertAlmostEqual(triangle.surface3d.frame.w.dot(point1 - point2), 0.)
            # Test distance from middle to surface

            self.assertLess(self.surface1.point_distance(triangle.middle()),
                            self.surface1.minor_radius * 0.05)

    def test_number_triangles(self):
        triangulation = self.face1.triangulation()
        triangulation.plot()
        n_triangles = len(triangulation.triangles)
        n_triangles_max = 250
        self.assertLess(n_triangles, n_triangles_max,
                        f'Too much triangles in toroidal face triangulation: {n_triangles}/{n_triangles_max}')

    def test_normal_at_point(self):
        toroidalsurface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        toroidaface = faces.ToroidalFace3D.from_surface_rectangular_cut(toroidalsurface, -1.4, 3.5, 0., 2.5)

        point = volmdlr.Point3D(1.1520373740641632, 2.008364391878863, 0.9489846193555862)

        normal = toroidaface.normal_at_point(point)
        self.assertTrue(normal.is_close(volmdlr.Vector3D(0.1568952782807088, 0.27351794069082946, 0.9489846193555862)))


if __name__ == '__main__':
    unittest.main()
