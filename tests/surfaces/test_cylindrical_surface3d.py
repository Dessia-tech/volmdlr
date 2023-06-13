"""
Unit tests for CylindriSurface3D
"""
import unittest
import math
import numpy as npy

import dessia_common.core
import volmdlr
from volmdlr import Z3D, Point2D, Point3D, edges, wires, surfaces
from volmdlr.models import cylindrical_surfaces


class TestCylindricalSurface3D(unittest.TestCase):
    cylindrical_surface = cylindrical_surfaces.cylindrical_surface
    cylindrical_surface2 = cylindrical_surfaces.cylindrical_surface2
    cylindrical_surface3 = cylindrical_surfaces.cylindrical_surface3
    cylindrical_surface4 = cylindrical_surfaces.cylindrical_surface4

    def test_line_intersections(self):
        line3d = edges.Line3D(volmdlr.O3D, volmdlr.Point3D(0.3, 0.3, .3))
        line_inters = self.cylindrical_surface.line_intersections(line3d)
        self.assertEqual(len(line_inters), 2)
        self.assertTrue(line_inters[0].is_close(volmdlr.Point3D(0.22627416, 0.22627416, 0.22627416)))
        self.assertTrue(line_inters[1].is_close(volmdlr.Point3D(-0.22627416, -0.22627416, -0.22627416)))

    def test_plane_intersections(self):
        plane_surface = surfaces.Plane3D(volmdlr.OZXY)
        parallel_plane_secant_cylinder = plane_surface.frame_mapping(volmdlr.Frame3D(
            volmdlr.O3D, volmdlr.Point3D(.5, 0, 0), volmdlr.Point3D(0, .5, 0), volmdlr.Point3D(0, 0, .2)), 'new')
        cylinder_tanget_plane = plane_surface.frame_mapping(volmdlr.Frame3D(
            volmdlr.Point3D(0, 0.32 / 2, 0), volmdlr.Point3D(.5, 0, 0),
            volmdlr.Point3D(0, .5, 0), volmdlr.Point3D(0, 0, .2)), 'new')
        not_intersecting_cylinder_parallel_plane = plane_surface.frame_mapping(volmdlr.Frame3D(
            volmdlr.Point3D(0, 0.32, 0), volmdlr.Point3D(.5, 0, 0),
            volmdlr.Point3D(0, .5, 0), volmdlr.Point3D(0, 0, .2)), 'new')
        cylinder_concurrent_plane = plane_surface.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 4)
        cylinder_perpendicular_plane = plane_surface.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 2)

        cylinder_surface_secant_parallel_plane_intersec = self.cylindrical_surface.plane_intersection(
            parallel_plane_secant_cylinder)
        self.assertEqual(len(cylinder_surface_secant_parallel_plane_intersec), 2)
        self.assertTrue(isinstance(cylinder_surface_secant_parallel_plane_intersec[0], edges.Line3D))
        self.assertTrue(isinstance(cylinder_surface_secant_parallel_plane_intersec[1], edges.Line3D))

        cylinder_surface_tangent_plane = self.cylindrical_surface.plane_intersection(
            cylinder_tanget_plane)
        self.assertEqual(len(cylinder_surface_tangent_plane), 1)
        self.assertTrue(isinstance(cylinder_surface_tangent_plane[0], edges.Line3D))

        cylinder_surface_tangent_plane_not_intersecting = self.cylindrical_surface.plane_intersection(
            not_intersecting_cylinder_parallel_plane)
        self.assertEqual(len(cylinder_surface_tangent_plane_not_intersecting), 0)

        cylinder_surface_concurrent_plane_intersec = self.cylindrical_surface.plane_intersection(
            cylinder_concurrent_plane)
        self.assertTrue(isinstance(cylinder_surface_concurrent_plane_intersec[0], wires.Ellipse3D))

        cylinder_surface_perpendicular_plane_intersec = self.cylindrical_surface.plane_intersection(
            cylinder_perpendicular_plane)
        self.assertTrue(isinstance(cylinder_surface_perpendicular_plane_intersec[0], wires.Circle3D))

    def test_is_coincident(self):
        cyl_surface1 = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 1)
        cyl_surface2 = surfaces.CylindricalSurface3D(volmdlr.OXYZ.translation(volmdlr.Vector3D(0, 0, 1)), 1)
        plane_face = surfaces.Plane3D(volmdlr.OXYZ)
        self.assertTrue(cyl_surface1.is_coincident(cyl_surface2))
        self.assertFalse(cyl_surface1.is_coincident(plane_face))
        self.assertFalse(self.cylindrical_surface.is_coincident(cyl_surface1))

    def test_point_on_surface(self):
        point = volmdlr.Point3D(0.32, 0, 1)
        point2 = volmdlr.Point3D(1, 1, 1)
        self.assertTrue(self.cylindrical_surface.point_on_surface(point))
        self.assertFalse((self.cylindrical_surface.point_on_surface(point2)))

    def test_arcellipse3d_to_2d(self):
        pass

    def test_arc3d_to_2d(self):
        arc1 = edges.Arc3D(volmdlr.Point3D(1, 0, 0), volmdlr.Point3D(1 / math.sqrt(2), 1 / math.sqrt(2), 0),
                           volmdlr.Point3D(0, 1, 0))
        arc2 = edges.Arc3D(volmdlr.Point3D(1, 0, 0), volmdlr.Point3D(1 / math.sqrt(2), -1 / math.sqrt(2), 0),
                           volmdlr.Point3D(0, -1, 0))
        arc3 = edges.Arc3D(volmdlr.Point3D(-1 / math.sqrt(2), 1 / math.sqrt(2), 0), volmdlr.Point3D(-1, 0, 0),
                           volmdlr.Point3D(-1 / math.sqrt(2), -1 / math.sqrt(2), 0))
        arc4 = edges.Arc3D(volmdlr.Point3D(0, -1, 0), volmdlr.Point3D(-1 / math.sqrt(2), 1 / math.sqrt(2), 0),
                           volmdlr.Point3D(1, 0, 0))
        test1 = self.cylindrical_surface2.arc3d_to_2d(arc3d=arc1)[0]
        test2 = self.cylindrical_surface2.arc3d_to_2d(arc3d=arc2)[0]
        test3 = self.cylindrical_surface2.arc3d_to_2d(arc3d=arc3)[0]
        test4 = self.cylindrical_surface2.arc3d_to_2d(arc3d=arc4)[0]

        inv_prof = self.cylindrical_surface2.linesegment2d_to_3d(test4)[0]

        # Assert that the returned object is an edges.LineSegment2D
        self.assertIsInstance(test1, edges.LineSegment2D)
        self.assertIsInstance(test2, edges.LineSegment2D)
        self.assertIsInstance(test3, edges.LineSegment2D)
        self.assertIsInstance(test4, edges.LineSegment2D)

        # Assert that the returned object is right on the parametric domain (take into account periodicity)
        self.assertEqual(test1.start, volmdlr.Point2D(0, 0))
        self.assertEqual(test1.end, volmdlr.Point2D(0.5 * math.pi, 0))
        self.assertEqual(test2.start, volmdlr.Point2D(0, 0))
        self.assertEqual(test2.end, volmdlr.Point2D(-0.5 * math.pi, 0))
        self.assertEqual(test3.start, volmdlr.Point2D(0.75 * math.pi, 0))
        self.assertEqual(test3.end, volmdlr.Point2D(1.25 * math.pi, 0))
        self.assertEqual(test4.start, volmdlr.Point2D(-0.5 * math.pi, 0))
        self.assertEqual(test4.end, volmdlr.Point2D(-2 * math.pi, 0))

        # Verifies the inversion operation
        self.assertIsInstance(inv_prof, edges.Arc3D)
        # self.assertEqual(inv_prof, arc4)
        self.assertTrue(inv_prof.start.is_close(arc4.start))
        self.assertTrue(inv_prof.interior.is_close(arc4.interior))
        self.assertTrue(inv_prof.end.is_close(arc4.end))

    def test_contour3d_to_2d(self):
        primitives_cylinder = [edges.LineSegment3D(Point3D(0.03, 0, 0.003), Point3D(0.03, 0, 0.013)),
                               edges.FullArc3D(Point3D(0, 0, 0.013), Point3D(0.03, 0, 0.013), Z3D),
                               edges.LineSegment3D(Point3D(0.03, 0, 0.013), Point3D(0.03, 0, 0.003)),
                               edges.FullArc3D(Point3D(0, 0, 0.003), Point3D(0.03, 0, 0.003), Z3D)
                               ]
        contour_cylinder = wires.Contour3D(primitives_cylinder)

        contour2d_cylinder = self.cylindrical_surface4.contour3d_to_2d(contour_cylinder)
        ax = contour2d_cylinder.plot()
        ax.set_aspect("auto")
        area = contour2d_cylinder.area()
        linesegment2d = contour2d_cylinder.primitives[3]
        fullarc2d = contour2d_cylinder.primitives[2]

        self.assertEqual(area, 0.02 * math.pi)
        self.assertEqual(fullarc2d.start, Point2D(volmdlr.TWO_PI, 0.003))
        self.assertEqual(fullarc2d.end, Point2D(0, 0.003))
        self.assertEqual(linesegment2d.start, Point2D(0, 0.003))
        self.assertEqual(linesegment2d.end, Point2D(0, 0.013))

        surface = dessia_common.core.DessiaObject.load_from_file(
            'surfaces/objects_cylindrical_tests/cylindrical_surface_bspline_openned_contour.json')
        contour = dessia_common.core.DessiaObject.load_from_file(
            'surfaces/objects_cylindrical_tests/cylindrical_contour_bspline_openned_contour.json')

        contour2d = surface.contour3d_to_2d(contour)
        self.assertEqual(len(contour2d.primitives), 2)
        self.assertFalse(contour2d.is_ordered())

        surface = dessia_common.core.DessiaObject.load_from_file(
            'surfaces/objects_cylindrical_tests/test_contour3d_to_2d_surface.json')
        contour = dessia_common.core.DessiaObject.load_from_file(
            'surfaces/objects_cylindrical_tests/test_contour3d_to_2d_contour.json')

        contour2d = surface.contour3d_to_2d(contour)
        self.assertAlmostEqual(contour2d.area(), 0.29361767646954695, 2)
        self.assertTrue(contour2d.is_ordered())

    def test_bsplinecurve3d_to_2d(self):
        surface = dessia_common.core.DessiaObject.load_from_file(
            'surfaces/objects_cylindrical_tests/cylindrical_surf_bug.json')
        bsplinecurve3d = dessia_common.core.DessiaObject.load_from_file(
            'surfaces/objects_cylindrical_tests/bsplinecurve3d_bug.json')
        primitive2d = surface.bsplinecurve3d_to_2d(bsplinecurve3d)[0]
        self.assertTrue(primitive2d.start.is_close(volmdlr.Point2D(-0.001540582016168617, -0.0006229082591074433)))
        self.assertTrue(primitive2d.end.is_close(volmdlr.Point2D(0.004940216577284154, -0.000847814405768888)))

        # Test to _fix_angle_discontinuity_on_discretization_points
        z = npy.linspace(0, 2 * math.pi, 50)
        theta = math.pi + 0.5 * math.pi * npy.cos(z)
        points_2d = [volmdlr.Point2D(x, y / (2 * math.pi)) for x, y in zip(theta, z)]
        cylinder = surfaces.CylindricalSurface3D(volmdlr.OXYZ, 1)
        points_3d = [cylinder.point2d_to_3d(point) for point in points_2d]
        bspline = edges.BSplineCurve3D.from_points_interpolation(points_3d, 3)
        result = cylinder.bsplinecurve3d_to_2d(bspline)[0]
        self.assertTrue(all(point.x < 0 for point in result.points))

    def test_point_projection(self):
        test_points = [Point3D(-2.0, -2.0, 0.0), Point3D(0.0, -2.0, 0.0), Point3D(2.0, -2.0, 0.0),
                       Point3D(2.0, 0.0, 0.0), Point3D(2.0, 2.0, 0.0), Point3D(0.0, 2.0, 0.0),
                       Point3D(-2.0, 2.0, 0.0), Point3D(-2.0, 0.0, 0.0)]
        expected_points = [volmdlr.Point3D(-0.5 * math.sqrt(2), -0.5 * math.sqrt(2), 0.0),
                           volmdlr.Point3D(0.0, -1.0, 0.0),
                           volmdlr.Point3D(0.5 * math.sqrt(2), -0.5 * math.sqrt(2), 0.0),
                           volmdlr.Point3D(1.0, 0.0, 0.0),
                           volmdlr.Point3D(0.5 * math.sqrt(2), 0.5 * math.sqrt(2), 0.0),
                           volmdlr.Point3D(0.0, 1.0, 0.0),
                           volmdlr.Point3D(-0.5 * math.sqrt(2), 0.5 * math.sqrt(2), 0.0),
                           volmdlr.Point3D(-1.0, 0.0, 0.0)]

        for i, point in enumerate(test_points):
            self.assertTrue(self.cylindrical_surface2.point_projection(point).is_close(expected_points[i]))

    def test_plot(self):
        ax = self.cylindrical_surface.plot()
        self.assertTrue(ax)


if __name__ == '__main__':
    unittest.main()
