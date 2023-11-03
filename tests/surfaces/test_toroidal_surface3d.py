import math
import os
import unittest
import numpy as np
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

    def test_arc3d_to_2d(self):
        arc1 = edges.Arc3D.from_3_points(volmdlr.Point3D(1-0.1/math.sqrt(2), 0, 0.1/math.sqrt(2)),
                           volmdlr.Point3D(0.9, 0, 0), volmdlr.Point3D(1-0.1/math.sqrt(2), 0, -0.1/math.sqrt(2)))

        test1 = self.toroidal_surface.arc3d_to_2d(arc3d=arc1)[0]

        # Assert that the returned object is an edges.LineSegment2D
        self.assertIsInstance(test1, edges.LineSegment2D)

        # Assert that the returned object is right on the parametric domain (take into account periodicity)
        self.assertTrue(test1.start.is_close(volmdlr.Point2D(0, 0.75 * math.pi)))
        self.assertTrue(test1.end.is_close(volmdlr.Point2D(0, 1.25 * math.pi)))

        arc2 = edges.Arc3D.from_3_points(volmdlr.Point3D(-0.169132244445, 0.06508125180570001, 0.627719515715),
        volmdlr.Point3D(-0.169169279223, 0.064939567779, 0.628073066814),
        volmdlr.Point3D(-0.169258691383, 0.064597504793, 0.628219515715))
        surface2 = surfaces.ToroidalSurface3D.load_from_file(os.path.join(folder, "surface.json"))
        test2 = surface2.arc3d_to_2d(arc3d=arc2)[0]
        self.assertTrue(test2.start.is_close(volmdlr.Point2D(-0.28681306221029024,  -math.pi)))
        self.assertTrue(test2.end.is_close(volmdlr.Point2D(-0.28681611209686064, -0.5 * math.pi)))

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

        self.assertTrue(test.start.is_close(volmdlr.Point2D(0.8489211153847066, math.pi)))
        self.assertTrue(test.end.is_close(volmdlr.Point2D(1.4449243890313308, 1.5707974196708867)))

        self.assertTrue(inv_prof.end.is_close(bspline_curve3d.end))

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
        surface = surfaces.ToroidalSurface3D.load_from_file(os.path.join(folder, "toroidal_surface_bug_2.json"))
        contour = wires.Contour3D.load_from_file(os.path.join(folder, "toroidal_surface_bug_2_contour_0.json"))
        contour2d = surface.contour3d_to_2d(contour)

        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 1.3773892114076673, 2)

        surface = surfaces.ToroidalSurface3D.load_from_file(os.path.join(folder, "buggy_toroidalface_surface.json"))
        contour = wires.Contour3D.load_from_file(os.path.join(folder, "buggy_toroidalface_contour.json"))
        contour2d = surface.contour3d_to_2d(contour)

        self.assertTrue(contour2d.is_ordered())
        self.assertAlmostEqual(contour2d.area(), 1.0990644259885822, 2)

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

    def test_plane_intersections(self):
        expected_results1 = [[18.84955592153876, 6.283185307179586], [18.774778972399815, 6.306324629324003],
                             [18.561749480682863, 6.382385910045063], [18.213003549947153, 6.522535326954354],
                             [17.73936325618824, 6.755617023729998], [17.162568850138605, 7.158697717832439],
                             [12.566370614359176, 12.566370614359176], [9.548698487359694, 9.548754306669936],
                             [8.513195638064289, 8.627615279189477], [7.859510527335956, 7.859510528088488]]
        expected_results2 = [18.007768707061835, 7.124972521656521]
        expected_results3 = [[6.283185307179586, 6.283185307179586], [6.287527793065231, 6.2875277930653075],
                             [6.304007940714406, 6.304008313559747], [6.402392698169013, 6.385653138168592],
                             [6.374756663452417, 6.3742078689160975], [6.432102115248554, 6.435832133353287],
                             [6.51052494933726, 6.510526578371843], [6.617587948267936, 6.6174395158125385],
                             [6.770799162158078, 6.770798442215525], [7.02765259262962, 7.027622387212875],
                             [14.078057920670332], [13.573542097891794], [13.223883969710258], [12.919853160635844],
                             [12.627496331426345], [12.329955764747575], [12.016630378969065], [11.679649471224833],
                             [11.312293366198755], [10.908108182904483]]
        toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)
        # Test 1
        plane1 = surfaces.Plane3D(volmdlr.OXYZ)
        plane1 = plane1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi / 4)
        for i, n in enumerate(np.linspace(0, math.pi / 4, 10)):
            plane = plane1.rotation(plane1.frame.origin, volmdlr.X3D, n)
            plane_intersections = toroidal_surface.plane_intersections(plane)
            for intersection, expected_result in zip(plane_intersections, expected_results1[i]):
                self.assertAlmostEqual(intersection.length(), expected_result)

        # Test 2
        plane2 = surfaces.Plane3D(volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0.5), volmdlr.X3D,
                                                  volmdlr.Y3D, volmdlr.Z3D))
        plane_intersections = toroidal_surface.plane_intersections(plane2)
        self.assertAlmostEqual(plane_intersections[0].length(), expected_results2[0])
        self.assertAlmostEqual(plane_intersections[1].length(), expected_results2[1])

        # Test 3
        plane3 = surfaces.Plane3D(volmdlr.OYZX)
        for i, n in enumerate(np.linspace(0, 2, 20)):
            plane = plane3.translation(n * volmdlr.X3D)
            plane_intersections = toroidal_surface.plane_intersections(plane)
            for intersection, expected_result in zip(plane_intersections, expected_results3[i]):
                self.assertAlmostEqual(intersection.length(), expected_result)
        # Test 4
        plane4 = surfaces.Plane3D(volmdlr.OYZX)
        plane4 = plane4.translation(volmdlr.X3D)
        plane_intersections = toroidal_surface.plane_intersections(plane4)
        for intersection, expected_result in zip(plane_intersections, [7.4155670073846425, 7.4149870554616255]):
            self.assertAlmostEqual(intersection.length(), expected_result)

        # Test 5
        plane5 = plane4.translation(volmdlr.X3D*3.1)
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

    def test_cylindrical_surface_intersections(self):
        toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 2, 1)

        # Test1
        frame = volmdlr.OXYZ.translation(volmdlr.Vector3D(1, 1, 0))
        frame = frame.rotation(volmdlr.Point3D(1, 1, 0), volmdlr.Y3D, math.pi / 4)
        cylindrical_surface = surfaces.CylindricalSurface3D(frame, 1)
        inters = toroidal_surface.cylindricalsurface_intersections(cylindrical_surface)
        self.assertEqual(len(inters), 1)
        self.assertAlmostEqual(inters[0].length(), 14.655915673071078)
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
        expected_results = [[17.155067582864444], [17.448537924734147], [8.189509579659324, 11.901146744658861],
                            [9.342209542119194, 6.783365484889515, 6.626606372190104],
                            [8.45252506578058, 11.78315465974606], [18.761705558048746],
                            [6.937783966922164, 15.192472907400315], [19.041705498182075], [19.712239534244855],
                            [9.106314389754695, 6.606874413821903, 6.606874200591648]]
        frame = volmdlr.OXYZ.translation(volmdlr.Vector3D(1, 1, 0))
        for i, theta in enumerate(np.linspace(0, math.pi * .7, 10)):
            frame = frame.rotation(frame.origin, volmdlr.Y3D, theta)
            cylindrical_surface = surfaces.CylindricalSurface3D(frame, 1.5)
            inters = toroidal_surface.cylindricalsurface_intersections(cylindrical_surface)
            for sol, expected_result in zip(inters, expected_results[i]):
                self.assertAlmostEqual(sol.length(), expected_result)

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
            volmdlr.Point3D(1.686564519293097, -1.027451632802102, -1.027451632802102)))
        self.assertTrue(ellipse_intersections[1].is_close(
            volmdlr.Point3D(1.8179453474500284, -1.1400021883520102, -1.1400021883520102)))


if __name__ == '__main__':
    unittest.main()
