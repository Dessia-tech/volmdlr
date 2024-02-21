import math
import os
import unittest
import numpy as np
from dessia_common.core import DessiaObject
import volmdlr
from volmdlr import edges, surfaces, wires, curves

folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_toroidal_tests')


class TestToroidalSurface3D(unittest.TestCase):
    toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 1, 0.1)
    frame = volmdlr.Frame3D(volmdlr.Point3D(-0.005829, 0.000765110438227, -0.0002349369830163),
                            volmdlr.Vector3D(-0.6607898454031987, 0.562158151695499, -0.4973278523210991),
                            volmdlr.Vector3D(-0.7505709694705869, -0.4949144228333324, 0.43783893597935386),
                            volmdlr.Vector3D(-0.0, 0.6625993710787045, 0.748974013865705))
    toroidal_surface2 = surfaces.ToroidalSurface3D(frame, 0.000725, 0.000125)
    linesegs = [
        edges.LineSegment3D(volmdlr.Point3D(4, 0, 0), volmdlr.Point3D(-4, 0.25, 0.25)),
        edges.LineSegment3D(volmdlr.Point3D(4, 0, 0), volmdlr.Point3D(-4, 0.25, 4)),
    ]

    def test_parametric_points_to_3d(self):
        parametric_points = np.array([[0.0, 0.0], [0.0, 0.5 * math.pi], [0.0, math.pi], [0.0, 1.5 * math.pi],
                                      [0.5 * math.pi, 0.0], [0.5 * math.pi, 0.5 * math.pi],
                                      [0.5 * math.pi, math.pi], [0.5 * math.pi, 1.5 * math.pi],
                                      [math.pi, 0.0], [math.pi, 0.5 * math.pi],
                                      [math.pi, math.pi], [math.pi, 1.5 * math.pi],
                                      [1.5 * math.pi, 0.0], [1.5 * math.pi, 0.5 * math.pi],
                                      [1.5 * math.pi, math.pi], [1.5 * math.pi, 1.5 * math.pi]])
        points3d = self.toroidal_surface.parametric_points_to_3d(parametric_points)
        expected_points = np.array([[1.1, 0.0, 0.0], [1.0, 0.0, 0.1], [0.9, 0.0, 0.0], [1.0, 0.0, -0.1],
                                    [0.0, 1.1, 0.0], [0.0, 1.0, 0.1], [0.0, 0.9, 0.0], [0.0, 1.0, -0.1],
                                    [-1.1, 0.0, 0.0], [-1.0, 0.0, 0.1], [-0.9, 0.0, 0.0], [-1.0, 0.0, -0.1],
                                    [0.0, -1.1, 0.0], [0.0, -1.0, 0.1], [0.0, -0.9, 0.0], [0.0, -1.0, -0.1]])
        for point, expected_point in zip(points3d, expected_points):
            self.assertAlmostEqual(np.linalg.norm(point - expected_point), 0.0)

    def test_arc3d_to_2d(self):
        arc1 = edges.Arc3D.from_3_points(volmdlr.Point3D(1 - 0.1 / math.sqrt(2), 0, 0.1 / math.sqrt(2)),
                                         volmdlr.Point3D(0.9, 0, 0),
                                         volmdlr.Point3D(1 - 0.1 / math.sqrt(2), 0, -0.1 / math.sqrt(2)))

        test1 = self.toroidal_surface.arc3d_to_2d(arc3d=arc1)[0]

        # Assert that the returned object is an edges.LineSegment2D
        self.assertIsInstance(test1, edges.LineSegment2D)

        # Assert that the returned object is right on the parametric domain (take into account periodicity)
        self.assertTrue(test1.start.is_close(volmdlr.Point2D(0, 0.75 * math.pi)))
        self.assertTrue(test1.end.is_close(volmdlr.Point2D(0, 1.25 * math.pi)))

        arc2 = edges.Arc3D.from_3_points(volmdlr.Point3D(-0.169132244445, 0.06508125180570001, 0.627719515715),
                                         volmdlr.Point3D(-0.169169279223, 0.064939567779, 0.628073066814),
                                         volmdlr.Point3D(-0.169258691383, 0.064597504793, 0.628219515715))

        surface2 = surfaces.ToroidalSurface3D.from_json(os.path.join(folder, "surface.json"))
        test2 = surface2.arc3d_to_2d(arc3d=arc2)[0]
        self.assertTrue(test2.start.is_close(volmdlr.Point2D(-0.2868131934235978, -math.pi)))
        self.assertTrue(test2.end.is_close(volmdlr.Point2D(-0.28681319342359773, -0.5 * math.pi)))

        surface = surfaces.ToroidalSurface3D.from_json(
            os.path.join(folder, "degenerated_toroidalsurface.json"))
        arc3d = edges.Arc3D.from_json(
            os.path.join(folder, "degenerated_toroidalsurface_arc3d_undefined_end.json"))
        brep_primitive = surface.arc3d_to_2d(arc3d)[0]
        inverse_prof = surface.linesegment2d_to_3d(brep_primitive)[0]
        self.assertAlmostEqual(brep_primitive.length(), 0.1993422098906592, 3)
        self.assertEqual(brep_primitive.end.y, -math.pi)
        self.assertAlmostEqual(arc3d.length(), inverse_prof.length(), 5)
        self.assertTrue(arc3d.start.is_close(inverse_prof.start))
        self.assertTrue(arc3d.end.is_close(inverse_prof.end))

        surface = surfaces.ToroidalSurface3D.from_json(
            os.path.join(folder, "degenerated_toroidalsurface_2.json"))
        arc3d = edges.Arc3D.from_json(
            os.path.join(folder, "degenerated_toroidalsurface_2_arc3d_undefined_end.json"))
        brep_primitive = surface.arc3d_to_2d(arc3d)[0]
        inverse_prof = surface.linesegment2d_to_3d(brep_primitive)[0]
        self.assertAlmostEqual(brep_primitive.length(), 0.5 * math.pi, 3)
        self.assertEqual(brep_primitive.end.y, math.pi)
        self.assertAlmostEqual(arc3d.length(), inverse_prof.length())
        self.assertTrue(arc3d.start.is_close(inverse_prof.start, 5e-5))
        self.assertTrue(arc3d.end.is_close(inverse_prof.end))

    def test_bsplinecurve3d_to_2d(self):
        control_points = [volmdlr.Point3D(-0.006429000000000001, 0.000765110438227, -0.0002349369830163),
                          volmdlr.Point3D(-0.006429000000000001, 0.0007527699876436001, -0.0002071780906932),
                          volmdlr.Point3D(-0.006429000000000001, 0.0007289073904888, -0.0001535567864537),
                          volmdlr.Point3D(-0.006429000000000001, 0.0006930461151679999, -7.904060141388999e-05),
                          volmdlr.Point3D(-0.006429000000000001, 0.0006567972565296, -1.323031929501e-05),
                          volmdlr.Point3D(-0.006429000000000001, 0.0006198714960685, 4.331237693835e-05),
                          volmdlr.Point3D(-0.006429000000000001, 0.0005818146300831, 9.111011896111e-05),
                          volmdlr.Point3D(-0.006429000000000001, 0.0005560693653727, 0.00011689162321630001),
                          volmdlr.Point3D(-0.006429000000000001, 0.0005431250195402, 0.00012834317593889998)]
        knot_multiplicities = [4, 1, 1, 1, 1, 1, 4]
        knots = [0.0, 0.1666666666667, 0.3333333333333, 0.5, 0.6666666666667, 0.8333333333333, 1.0]
        bspline_curve3d = edges.BSplineCurve3D(3, control_points, knot_multiplicities, knots)

        test = self.toroidal_surface2.bsplinecurve3d_to_2d(bspline_curve3d)[0]
        inv_prof = self.toroidal_surface2.bsplinecurve2d_to_3d(test)[0]

        self.assertTrue(test.start.is_close(volmdlr.Point2D(0.8489257266674017, math.pi)))
        self.assertTrue(test.end.is_close(volmdlr.Point2D(1.4449281658550368, 1.5707952339189064)))

        self.assertTrue(inv_prof.end.is_close(bspline_curve3d.end))

        surface = surfaces.ToroidalSurface3D.from_json(
            os.path.join(folder, "toroidalsurface_bsplinecurve3d_to_2d.json"))
        bspline_curve3d = edges.BSplineCurve3D.from_json(
            os.path.join(folder, "toroidalsurface_bsplinecurve3d_to_2d_curve.json"))
        brep_primitive = surface.bsplinecurve3d_to_2d(bspline_curve3d)[0]
        inverse_prof = surface.bsplinecurve2d_to_3d(brep_primitive)[0]
        self.assertAlmostEqual(brep_primitive.length(), 0.013265398542202636, 3)
        self.assertAlmostEqual(bspline_curve3d.length(), inverse_prof.length(), 5)
        self.assertTrue(bspline_curve3d.start.is_close(inverse_prof.start))
        self.assertTrue(bspline_curve3d.end.is_close(inverse_prof.end))

    def test_point_projection(self):
        test_points = [volmdlr.Point3D(-2.0, -2.0, 0.0), volmdlr.Point3D(0.0, -2.0, 0.0),
                       volmdlr.Point3D(2.0, -2.0, 0.0),
                       volmdlr.Point3D(2.0, 0.0, 0.0), volmdlr.Point3D(2.0, 2.0, 0.0), volmdlr.Point3D(0.0, 2.0, 0.0),
                       volmdlr.Point3D(-2.0, 2.0, 0.0), volmdlr.Point3D(-2.0, 0.0, 0.0),
                       ]

        expected_points = [volmdlr.Point3D(-0.55 * math.sqrt(2), -0.55 * math.sqrt(2), 0.0),
                           volmdlr.Point3D(0.0, -1.1, 0.0),
                           volmdlr.Point3D(0.55 * math.sqrt(2), -0.55 * math.sqrt(2), 0.0),
                           volmdlr.Point3D(1.1, 0.0, 0.0),
                           volmdlr.Point3D(0.55 * math.sqrt(2), 0.55 * math.sqrt(2), 0.0),
                           volmdlr.Point3D(0.0, 1.1, 0.0),
                           volmdlr.Point3D(-0.55 * math.sqrt(2), 0.55 * math.sqrt(2), 0.0),
                           volmdlr.Point3D(-1.1, 0.0, 0.0),
                           ]

        for i, point in enumerate(test_points):
            self.assertTrue(self.toroidal_surface.point_projection(point).is_close(expected_points[i]))

    def test_contour3d_to_2d(self):
        surface = surfaces.ToroidalSurface3D.from_json(os.path.join(folder, "toroidal_surface_bug_2.json"))
        contour = wires.Contour3D.from_json(os.path.join(folder, "toroidal_surface_bug_2_contour_0.json"))
        contour2d = surface.contour3d_to_2d(contour)

        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 1.3773892114076673, 2)

        surface = surfaces.ToroidalSurface3D.from_json(os.path.join(folder, "buggy_toroidalface_surface.json"))
        contour = wires.Contour3D.from_json(os.path.join(folder, "buggy_toroidalface_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)

        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 1.0990644259885822, 2)

        surface = surfaces.ToroidalSurface3D.from_json(
            os.path.join(folder, "toroidalsurface_small_bsplinecurve.json"))
        contour = wires.Contour3D.from_json(
            os.path.join(folder, "toroidalsurface_small_bsplinecurve_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)

        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.12142017346476651, 4)

        surface = surfaces.ToroidalSurface3D.from_json(os.path.join(folder, "toroidalsurface_small_arc3d.json"))
        contour = wires.Contour3D.from_json(os.path.join(folder, "toroidalsurface_small_arc3d_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)

        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 0.09543353484687866, 2)

    def test_line_intersections(self):
        expected_results = [[volmdlr.Point3D(2.9993479584651066, 0.031270376297965426, 0.031270376297965426),
                             volmdlr.Point3D(1.0000193965498871, 0.09374939385781603, 0.09374939385781603),
                             volmdlr.Point3D(-1.0001508657814862, 0.15625471455567144, 0.15625471455567144),
                             volmdlr.Point3D(-2.968027405412837, 0.21775085641915115, 0.21775085641915115)],
                            [volmdlr.Point3D(2.799597842042955, 0.03751256743615766, 0.6002010789785226),
                             volmdlr.Point3D(2.000000955072264, 0.06249997015399175, 0.999999522463868)]]

        toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        for i, lineseg in enumerate(self.linesegs):
            inters = toroidal_surface.line_intersections(lineseg.line)
            for expected_result, inter in zip(expected_results[i], inters):
                self.assertTrue(expected_result.is_close(inter))
        expected_results = [[volmdlr.Point3D(-1.0385780861224632, -2.73372290825405, -0.3815197706145943),
                             volmdlr.Point3D(-1.0385780861224632, -0.27992053908174386, -0.3815197706145943),
                             volmdlr.Point3D(-1.0385780861224632, 0.2799205390817603, -0.3815197706145943),
                             volmdlr.Point3D(-1.0385780861224632, 2.7337229082540286, -0.3815197706145943)],
                            [volmdlr.Point3D(0.02037825356907985, -1.7074248427439929, -1.7074248427439933),
                             volmdlr.Point3D(1.1557561697724674, -2.2751138008456864, -2.2751138008456877)],
                            [volmdlr.Point3D(2.8868233917320207, 2.0, -0.8590189402508468),
                             volmdlr.Point3D(3.3891999753110165, 2.0, 0.353825421244224)],
                            [volmdlr.Point3D(1.3426661840222276, -1.6569652433720718, -0.9911601091085722),
                             volmdlr.Point3D(2.6785064972435375, -1.252458418216524, -0.2905337359307829)],
                            []]
        surface1, lineseg1 = DessiaObject.from_json(os.path.join(folder, "test_torus_line_intersections.json")).primitives
        surface2, line2 = DessiaObject.from_json(os.path.join(folder, "test_torus_line_itnersections_08_11_2023.json")).primitives
        surface3, lineseg3 = DessiaObject.from_json(os.path.join(folder, "test_torus_lineseg141223.json")).primitives
        surface4, lineseg4 = DessiaObject.from_json(os.path.join(folder, "test_toroidal_surface_lineseg_intersections201223.json")).primitives
        surface5, lineseg5 = DessiaObject.from_json(os.path.join(folder, "test_toroidal_surface_line_intersections.json")).primitives
        for i, (surface, line) in enumerate([[surface1, lineseg1.line], [surface2, line2], [surface3, lineseg3.line],
                                             [surface4, lineseg4.line], [surface5, lineseg5.line]]):
            line_intersections = surface.line_intersections(line)
            self.assertEqual(len(line_intersections), len(expected_results[i]))
            for result, expected_result in zip(line_intersections, expected_results[i]):
                self.assertTrue(result.is_close(expected_result))

    def test_plane_intersections(self):
        expected_results2 = [18.007768707061828, 7.124972521656522]

        toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        # Test 1
        plane1 = surfaces.Plane3D(volmdlr.OXYZ)
        plane1 = plane1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi / 4)
        for i, n in enumerate(np.linspace(0, math.pi / 4, 10)):
            plane = plane1.rotation(plane1.frame.origin, volmdlr.X3D, n)
            plane_intersections = toroidal_surface.plane_intersections(plane)
            for intersection in plane_intersections:
                for p in intersection.discretization_points(number_points=50):
                    self.assertLess(toroidal_surface.point_distance(p), 1e-5)
                    self.assertLess(plane.point_distance(p), 1e-5)
            # for intersection, expected_result in zip(plane_intersections, expected_results1[i]):
            #     self.assertAlmostEqual(intersection.length(), expected_result, 5)

        # Test 2
        plane2 = surfaces.Plane3D(volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0.5), volmdlr.X3D,
                                                  volmdlr.Y3D, volmdlr.Z3D))
        plane_intersections = toroidal_surface.plane_intersections(plane2)
        self.assertAlmostEqual(plane_intersections[0].length(), expected_results2[0], 6)
        self.assertAlmostEqual(plane_intersections[1].length(), expected_results2[1], 6)

        # Test 3
        plane3 = surfaces.Plane3D(volmdlr.OYZX)
        for i, n in enumerate(np.linspace(0, 2, 20)):
            plane = plane3.translation(n * volmdlr.X3D)
            plane_intersections = toroidal_surface.plane_intersections(plane)
            for intersection in plane_intersections:
                for p in intersection.discretization_points(number_points=50):
                    self.assertLess(toroidal_surface.point_distance(p), 1e-5)
                    self.assertLess(plane.point_distance(p), 1e-5)
            # for intersection, expected_result in zip(plane_intersections, expected_results3[i]):
            #     self.assertAlmostEqual(intersection.length(), expected_result, 6)
        # Test 4
        plane4 = surfaces.Plane3D(volmdlr.OYZX)
        plane4 = plane4.translation(volmdlr.X3D)
        plane_intersections = toroidal_surface.plane_intersections(plane4)
        for intersection in plane_intersections:
            for p in intersection.discretization_points(number_points=50):
                self.assertLess(toroidal_surface.point_distance(p), 1e-5)
                self.assertLess(plane4.point_distance(p), 1e-5)
        # for intersection, expected_result in zip(plane_intersections, [7.415366424519409, 7.415366424519409]):
        #     self.assertAlmostEqual(intersection.length(), expected_result, 6)

        # Test 5
        plane5 = plane4.translation(volmdlr.X3D * 3.1)
        plane_intersections = toroidal_surface.plane_intersections(plane5)
        self.assertFalse(plane_intersections)

        # Test 6
        plane6 = surfaces.Plane3D(
            volmdlr.Frame3D(origin=volmdlr.Point3D(2.265348976860137, 1.0, 1.2653489768601376),
                            u=volmdlr.Vector3D(0.7071067811865476, 0.0, -0.7071067811865475),
                            v=volmdlr.Vector3D(0.0, 1.0, 0.0),
                            w=volmdlr.Vector3D(0.7071067811865475, 0.0, 0.7071067811865476)))
        plane_intersections = toroidal_surface.plane_intersections(plane6)
        self.assertFalse(plane_intersections)
        toroidalsurface, plane = DessiaObject.from_json(
            os.path.join(folder, 'test_toroidalsurface_plane3d_intersections_211223.json')).primitives
        intersections = toroidalsurface.surface_intersections(plane)
        self.assertEqual(len(intersections), 2)
        self.assertTrue(intersections[0].center.is_close(volmdlr.Point3D(3.0, 0.0, 0.0)))
        self.assertEqual(intersections[0].radius, 1)
        self.assertTrue(intersections[1].center.is_close(volmdlr.Point3D(-3.0, 0.0, 0.0)))
        self.assertEqual(intersections[1].radius, 1)

    def test_cylindrical_surface_intersections(self):
        toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)

        # Test1
        frame = volmdlr.OXYZ.translation(volmdlr.Vector3D(1, 1, 0))
        frame = frame.rotation(volmdlr.Point3D(1, 1, 0), volmdlr.Y3D, math.pi / 4)
        cylindrical_surface = surfaces.CylindricalSurface3D(frame, 1)
        inters = toroidal_surface.cylindricalsurface_intersections(cylindrical_surface)
        self.assertEqual(len(inters), 1)
        self.assertAlmostEqual(inters[0].length(),   14.655771126896285, 6)
        # Test2
        expected_results = [[9.424777944721708, 9.424777944721708], [6.283185307179586], []]
        frame = volmdlr.OXYZ
        cylindrical_surfaces = [surfaces.CylindricalSurface3D(frame, 1.5),
                                surfaces.CylindricalSurface3D(frame, 1),
                                surfaces.CylindricalSurface3D(frame, 0.9)]
        for i, surface in enumerate(cylindrical_surfaces):
            inters = toroidal_surface.cylindricalsurface_intersections(surface)
            for sol, expected_result in zip(inters, expected_results[i]):
                self.assertAlmostEqual(sol.length(), expected_result)

        #Test3
        expected_results = [[17.15507502094234], [17.44854519606042], [8.189776671441997, 11.901135669170262],
                            [9.342188106943269, 6.783371061263169, 6.6266277842571295],
                            [8.454952065863425, 11.776550916194452], [18.761719845054934],
                            [6.937795281803973, 15.192491122547677], [19.04178257950678], [19.712211179693842],
                            [9.106322135020985, 6.606873336946121, 6.606872989299915]]

        frame = volmdlr.OXYZ.translation(volmdlr.Vector3D(1, 1, 0))
        for i, theta in enumerate(np.linspace(0, math.pi * .7, 10)):
            frame = frame.rotation(frame.origin, volmdlr.Y3D, theta)
            cylindrical_surface = surfaces.CylindricalSurface3D(frame, 1.5)
            inters = toroidal_surface.cylindricalsurface_intersections(cylindrical_surface)
            for sol, expected_result in zip(inters, expected_results[i]):
                self.assertAlmostEqual(sol.length(), expected_result, 5)

    def test_circle_intersections(self):
        toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        circle = curves.Circle3D(volmdlr.Frame3D(origin=volmdlr.Point3D(1.0, 1.0, -0.8947368421052632),
                                                 u=volmdlr.Vector3D(1.0, 0.0, 0.0),
                                                 v=volmdlr.Vector3D(0.0, 1.0, 0.0),
                                                 w=volmdlr.Vector3D(0.0, 0.0, 1.0)), 1)
        circle_intersections = toroidal_surface.circle_intersections(circle)
        expected_point1 = volmdlr.Point3D(1.544982741074, 0.161552737537, -0.894736842105)
        expected_point2 = volmdlr.Point3D(0.161552737537, 1.544982741074, -0.894736842105)
        self.assertTrue(circle_intersections[0].is_close(expected_point1))
        self.assertTrue(circle_intersections[1].is_close(expected_point2))
        torus, circle = DessiaObject.from_json(os.path.join(folder,
            'test_toroidalsurface_circle_intersections211223_2.json')).primitives
        circle_intersections = torus.circle_intersections(circle)
        self.assertEqual(len(circle_intersections), 2)
        self.assertTrue(circle_intersections[0].is_close(volmdlr.Point3D(2.0000006438528177, -0.5135128860482583, -0.9978935668376178)))
        self.assertTrue(circle_intersections[1].is_close(volmdlr.Point3D(2.0000002080103414, -0.5135127741429286, 0.9978935960903826)))

    def test_ellipse_intersections(self):
        toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.Frame3D(origin=volmdlr.Point3D(1.0, 1.0, 0.0),
                                                                      u=volmdlr.Vector3D(-5.551115123125783e-17, 0.0,
                                                                                         0.9999999999999998),
                                                                      v=volmdlr.Vector3D(0.0, 0.9999999999999998, 0.0),
                                                                      w=volmdlr.Vector3D(-0.9999999999999998, 0.0,
                                                                                         -5.551115123125783e-17)), 3,
                                                      1)

        frame = volmdlr.Frame3D(origin=volmdlr.Point3D(0.0, 0.0, 0.0),
                                u=volmdlr.Vector3D(0.5773502691896258, 0.5773502691896258, 0.5773502691896258),
                                v=volmdlr.Vector3D(0.8164965809277258, -0.40824829046386313, -0.40824829046386313),
                                w=volmdlr.Vector3D(0.0, 0.7071067811865476, -0.7071067811865476))

        ellipse = curves.Ellipse3D(2, 1, frame)
        ellipse_intersections = toroidal_surface.ellipse_intersections(ellipse)
        self.assertFalse(ellipse_intersections)

        frame1 = frame.translation(volmdlr.Vector3D(3, 0.0, 0.0))
        ellipse = curves.Ellipse3D(7, 2.5, frame1)
        ellipse_intersections = toroidal_surface.ellipse_intersections(ellipse)
        self.assertFalse(ellipse_intersections)

        frame = frame.translation(volmdlr.Vector3D(3, 0.0, 0.0))
        ellipse = curves.Ellipse3D(2, 1, frame)
        ellipse_intersections = toroidal_surface.ellipse_intersections(ellipse)
        self.assertEqual(len(ellipse_intersections), 2)
        self.assertTrue(ellipse_intersections[0].is_close(
            volmdlr.Point3D(1.6865642161149017, -1.0274512473410842, -1.0274512473410844)))
        self.assertTrue(ellipse_intersections[1].is_close(
            volmdlr.Point3D(1.817953260018375, -1.1400067506585763, -1.1400067506585763)))

    def test_conicalsurface_intersections(self):
        conical_surface = surfaces.ConicalSurface3D(volmdlr.OXYZ, math.pi / 7)
        conical_surface = conical_surface.translation(volmdlr.Vector3D(2, 2, -3))
        toroidal_surface1 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 3, 1)
        list_curves = toroidal_surface1.conicalsurface_intersections(conical_surface)
        self.assertEqual(len(list_curves), 2)
        self.assertAlmostEqual(list_curves[0].length(), 7.290767246711664)
        self.assertAlmostEqual(list_curves[1].length(),  7.290781630732165)

        conical_surface = surfaces.ConicalSurface3D(volmdlr.OXYZ, math.pi / 8)
        conical_surface = conical_surface.translation(volmdlr.Vector3D(2, 2, -3))
        toroidal_surface1 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 3, 1)
        list_curves = toroidal_surface1.conicalsurface_intersections(conical_surface)
        self.assertEqual(len(list_curves), 1)
        self.assertAlmostEqual(list_curves[0].length(), 15.26648920774545, 6)

    def test_sphericalsurface_intersections(self):
        spherical_surface = surfaces.SphericalSurface3D(
            volmdlr.OXYZ.translation(volmdlr.Vector3D(0.5, 0.5, 0)), 2)
        frame = volmdlr.OXYZ
        toroidal_surface1 = surfaces.ToroidalSurface3D(frame, 2, 1)

        intersections = toroidal_surface1.sphericalsurface_intersections(spherical_surface)
        self.assertEqual(len(intersections), 2)
        self.assertAlmostEqual(intersections[0].length(), 11.364812376610685)
        self.assertAlmostEqual(intersections[1].length(), 11.364812376610685)
        frame = frame.rotation(frame.origin, volmdlr.Y3D, math.pi / 5)
        toroidal_surface2 = surfaces.ToroidalSurface3D(frame, 2, 1)
        intersections = toroidal_surface2.sphericalsurface_intersections(spherical_surface)
        self.assertEqual(len(intersections), 2)
        self.assertAlmostEqual(intersections[0].length(), 10.264046962680238)
        self.assertAlmostEqual(intersections[1].length(), 12.024102432013244)
        frame = volmdlr.OXYZ.rotation(frame.origin, volmdlr.Y3D, math.pi / 5)
        frame = frame.translation(volmdlr.X3D * 1.6)
        toroidal_surface3 = surfaces.ToroidalSurface3D(frame, 2, 1)
        intersections = toroidal_surface3.sphericalsurface_intersections(spherical_surface)
        self.assertEqual(len(intersections), 1)
        self.assertAlmostEqual(intersections[0].length(), 20.514857203053506)

    def test_toroidal_surfaces(self):
        toroidal_surface1 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        """ ======= INCLINED AND TRANSLATED ========== """
        toroidal_surface2 = surfaces.ToroidalSurface3D(volmdlr.OYZX, 2, 1)
        toroidal_surface2 = toroidal_surface2.rotation(volmdlr.O3D, volmdlr.Y3D, math.pi / 6)
        toroidal_surface2 = toroidal_surface2.translation(volmdlr.X3D * 2.5)

        inters = toroidal_surface1.surface_intersections(toroidal_surface2)
        self.assertEqual(len(inters), 1)
        for i in inters:
            for p in i.discretization_points(number_points=50):
                self.assertFalse(toroidal_surface1.point_distance(p) > 1e-4)
                self.assertFalse(toroidal_surface2.point_distance(p) > 1e-4)

        """" ========================# PARALLEL NOT INTERSECTING ========================"""
        toroidal_surface2 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        toroidal_surface2 = toroidal_surface2.translation(volmdlr.Z3D * 2.5)
        inters = toroidal_surface1.surface_intersections(toroidal_surface2)
        self.assertFalse(inters)

        """  ======================== # PARALLEL INTERSECTING  ======================== """
        surfaces2 = []
        toroidal_surface2 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 3, 1)
        surfaces2.append(toroidal_surface2)
        toroidal_surface2 = toroidal_surface2.translation(volmdlr.Z3D * 1.5)
        surfaces2.append(toroidal_surface2)
        toroidal_surface2 = toroidal_surface2.translation(volmdlr.Z3D * 0.5)
        surfaces2.append(toroidal_surface2)

        toroidal_surface2 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 5, 3)
        surfaces2.append(toroidal_surface2)
        toroidal_surface2 = toroidal_surface2.translation(volmdlr.Z3D * 1.5)
        surfaces2.append(toroidal_surface2)
        toroidal_surface2 = toroidal_surface2.translation(volmdlr.Z3D * 0.5)
        surfaces2.append(toroidal_surface2)

        toroidal_surface2 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 3, 2)
        surfaces2.append(toroidal_surface2)
        toroidal_surface2 = toroidal_surface2.translation(volmdlr.Z3D * 1.5)
        surfaces2.append(toroidal_surface2)
        toroidal_surface2 = toroidal_surface2.translation(volmdlr.Z3D * 0.5)
        surfaces2.append(toroidal_surface2)

        toroidal_surface2 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 3, 2.5)
        expected_number_sol = [2, 2, 0, 2, 2, 2, 1, 2, 2, 0]

        surfaces2.append(toroidal_surface2)
        for i, surface2 in enumerate(surfaces2):
            inters = toroidal_surface1.surface_intersections(surface2)
            self.assertEqual(len(inters), expected_number_sol[i])
            for inter in inters:
                for p in inter.discretization_points(number_points=50):
                    self.assertFalse(toroidal_surface1.point_distance(p) > 1e-6)
                    self.assertFalse(surface2.point_distance(p) > 1e-6)

        """  ======================== # PARALLEL and not coincident INTERSECTING  ======================== """
        toroidal_surface2 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        toroidal_surface2 = toroidal_surface2.translation(volmdlr.X3D * 4)

        toroidal_surface1_1 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 1, .5)
        toroidal_surface2_1 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 0.8, .5)
        toroidal_surface2_1 = toroidal_surface2_1.translation(volmdlr.X3D * (-1.8))
        toroidal_surface3_1 = toroidal_surface2_1.translation(volmdlr.X3D * (-0.1))
        expected_number_sol = [2, 2, 1]
        for i, (sf1, sf2) in enumerate([(toroidal_surface1, toroidal_surface2),
                         (toroidal_surface1_1, toroidal_surface2_1),
                         (toroidal_surface1_1, toroidal_surface3_1)]):
            inters = sf1.surface_intersections(sf2)
            self.assertEqual(len(inters), expected_number_sol[i])
            for inter in inters:
                for p in inter.discretization_points(number_points=50):
                    self.assertFalse(sf1.point_distance(p) > 1e-5)
                    self.assertFalse(sf2.point_distance(p) > 1e-5)
        """ ==================== Yvone-Villarceau circles of T1 and T2 ==========================="""

        toroidal_surface1 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 1, .5)
        toroidal_surface2 = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 1, .3)
        toroidal_surface2_1 = toroidal_surface2.translation(volmdlr.X3D * 0.8)

        toroidal_surface2_2 = toroidal_surface1.translation(volmdlr.X3D)
        toroidal_surface2_3 = toroidal_surface1.translation(volmdlr.X3D * 1.1)
        toroidal_surface2_4 = toroidal_surface1.translation(volmdlr.X3D * 1.8)
        toroidal_surface2_5 = toroidal_surface1.translation(volmdlr.X3D * 0.8)
        expected_number_sol = [4, 4, 3, 3, 4]
        expected_sols_lengths = [[3.4902240711559385, 3.4902240711559585, 2.8025110017478196, 2.802511001747819],
                                 [6.283185307179586, 6.283185307179586, 3.707738420898486, 3.707738420898486],
                                 [6.907653689757426, 5.0272056418053115, 5.027217971367434],
                                 [5.82219814019078, 3.3338714185008955, 3.3338735379661655],
                                 [3.351031375990407, 3.351031375990407, 6.088038294280911, 6.088038257995996]]
        for i, toroidal_surface2 in enumerate([toroidal_surface2_1,  toroidal_surface2_2, toroidal_surface2_3,
                                               toroidal_surface2_4, toroidal_surface2_5]):
            inters = toroidal_surface1.surface_intersections(toroidal_surface2)
            self.assertEqual(len(inters), expected_number_sol[i])
            for inter, expected_inter_length in zip(inters, expected_sols_lengths[i]):
                self.assertAlmostEqual(inter.length(), expected_inter_length, 6)


if __name__ == '__main__':
    unittest.main()
