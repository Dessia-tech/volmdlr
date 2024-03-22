import unittest
import math
import os

import volmdlr
from volmdlr.faces import SphericalFace3D
from volmdlr import surfaces, wires
from dessia_common.core import DessiaObject
from OCP.GProp import GProp_GProps
from OCP.BRepGProp import BRepGProp_Face, BRepGProp  # used for mass calculation

folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_spherical_test')


class TestSphericalFace3D(unittest.TestCase):
    def test_from_contours3d_and_rectangular_cut(self):
        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "face_from_contours3d_and_rectangular_cut_surface.json"))
        point = DessiaObject.from_json(
            os.path.join(folder, "face_from_contours3d_and_rectangular_cut_point.json"))
        contour_0 = DessiaObject.from_json(
            os.path.join(folder, "face_from_contours3d_and_rectangular_cut_contour_0.json"))
        contour_1 = DessiaObject.from_json(
            os.path.join(folder, "face_from_contours3d_and_rectangular_cut_contour_1.json"))
        face = SphericalFace3D.from_contours3d_and_rectangular_cut(surface, [contour_0, contour_1], point)
        self.assertEqual(len(face.surface2d.inner_contours), 2)
        self.assertAlmostEqual(face.surface2d.outer_contour.area(), math.pi ** 2 * 2)

    def test_from_contours3d(self):
        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "sphericalface_disconnected_contours_surface.json"))
        contour_0 = DessiaObject.from_json(
            os.path.join(folder, "sphericalface_disconnected_contours_contour_0.json"))
        contour_1 = DessiaObject.from_json(
            os.path.join(folder, "sphericalface_disconnected_contours_contour_1.json"))
        face = SphericalFace3D.from_contours3d(surface, [contour_0, contour_1])
        self.assertAlmostEqual(face.surface2d.area(), 3.7274576655120804, 2)

        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "sphericalsurface_contour3d_to_2d_bug.json"))
        contour = DessiaObject.from_json(
            os.path.join(folder, "sphericalsurface_contour3d_to_2d_bug_contour.json"))
        face = SphericalFace3D.from_contours3d(surface, [contour])
        self.assertTrue(face.triangulation())

        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "sphericalface_from_contours3d_repair_primitives_periodicity_surface.json"))
        contour = DessiaObject.from_json(
            os.path.join(folder, "sphericalface_from_contours3d_repair_primitives_periodicity_contour.json"))
        face = SphericalFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(face.surface2d.area(), math.pi * math.pi)
        self.assertTrue(face.triangulation())

        surface = surfaces.SphericalSurface3D.from_json(
            os.path.join(folder, "sphericalface_from_contours3d_repair_primitives_periodicity_surface.json"))
        contour = DessiaObject.from_json(
            os.path.join(folder, "sphericalface_from_contours3d_repair_primitives_periodicity_contour.json"))
        face = SphericalFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(face.surface2d.area(), math.pi * math.pi)
        self.assertTrue(face.triangulation())

    def test_grid_points(self):
        surface3d = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)
        outer_contour2d = wires.Contour2D.rectangle(-math.pi, math.pi, -0.5 * math.pi, 0.5 * math.pi)
        inner_contour2d = wires.Contour2D.rectangle_from_center_and_sides(
            volmdlr.Point2D(0.0, 0.0), 0.5 * math.pi, 0.5 * math.pi, False)
        surface2d = surfaces.Surface2D(outer_contour2d, [inner_contour2d])
        face = SphericalFace3D(surface3d, surface2d)
        grid_points = face.grid_points([10, 10])
        self.assertEqual(len(grid_points), 518)

    def test_boudingbox(self):
        spherical_surface = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)
        spherical_face = volmdlr.faces.SphericalFace3D.from_surface_rectangular_cut(
            spherical_surface, 0, math.pi)
        self.assertAlmostEqual(spherical_face.bounding_box.volume(), 4.0, delta=0.1)
        spherical_face = volmdlr.faces.SphericalFace3D.from_surface_rectangular_cut(
            spherical_surface, 0, 2 * math.pi)
        self.assertAlmostEqual(spherical_face.bounding_box.volume(), 8.0, delta=0.1)
        spherical_surface = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)
        spherical_face = volmdlr.faces.SphericalFace3D.from_surface_rectangular_cut(
            spherical_surface, 0, 1.5 * math.pi)
        self.assertAlmostEqual(spherical_face.bounding_box.volume(), 8.0, delta=0.1)
        phi1 = math.asin(0.5 / spherical_surface.radius)
        spherical_face = volmdlr.faces.SphericalFace3D.from_surface_rectangular_cut(
            spherical_surface, 0, 2 * math.pi, phi1, 0.5 * math.pi)
        self.assertAlmostEqual(spherical_face.bounding_box.volume(),
                               0.5 * (2 * spherical_surface.radius * math.cos(phi1)) ** 2, delta=0.1)

    def test_planeface_intersections(self):
        planeface, sphericalface = DessiaObject.from_json(
            os.path.join(folder, "test_planeface_sphericalface_intersections.json")).primitives
        intersections = planeface.face_intersections(sphericalface)
        self.assertEqual(len(intersections), 1)
        self.assertAlmostEqual(intersections[0].length(), 2.8957237263187805)
        planeface, sphericalface = DessiaObject.from_json(
            os.path.join(folder, "test_planeface_sphericalface_intersections_2.json")).primitives
        intersections = planeface.face_intersections(sphericalface)
        self.assertEqual(len(intersections), 1)
        self.assertAlmostEqual(intersections[0].length(), 2.8957237263187805)

    def test_to_ocp(self):
        sphere = surfaces.SphericalSurface3D(volmdlr.OXYZ, 1)
        spherical_face = SphericalFace3D.from_surface_rectangular_cut(spherical_surface=sphere)
        ocp_face = spherical_face.to_ocp()
        properties = GProp_GProps()
        BRepGProp.SurfaceProperties_s(ocp_face, properties)
        face_area = properties.Mass()
        self.assertAlmostEqual(face_area, 4 * math.pi) # SI


if __name__ == '__main__':
    unittest.main()
