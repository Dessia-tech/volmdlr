import math
import unittest

import volmdlr
from volmdlr import edges, faces, wires, surfaces


class TestToroidalSurface3D(unittest.TestCase):
    toroidal_surface = surfaces.ToroidalSurface3D(volmdlr.OXYZ, 1, 0.1)
    frame = volmdlr.Frame3D(volmdlr.Point3D(-0.005829, 0.000765110438227, -0.0002349369830163),
                            volmdlr.Vector3D(-0.6607898454031987, 0.562158151695499, -0.4973278523210991),
                            volmdlr.Vector3D(-0.7505709694705869, -0.4949144228333324, 0.43783893597935386),
                            volmdlr.Vector3D(-0.0, 0.6625993710787045, 0.748974013865705))
    toroidal_surface2 = surfaces.ToroidalSurface3D(frame, 0.000725, 0.000125)

    def test_arc3d_to_2d(self):
        arc1 = edges.Arc3D(volmdlr.Point3D(1-0.1/math.sqrt(2), 0, 0.1/math.sqrt(2)),
                           volmdlr.Point3D(0.9, 0, 0), volmdlr.Point3D(1-0.1/math.sqrt(2), 0, -0.1/math.sqrt(2)))

        test1 = self.toroidal_surface.arc3d_to_2d(arc3d=arc1)[0]

        # Assert that the returned object is an edges.LineSegment2D
        self.assertIsInstance(test1, edges.LineSegment2D)

        # Assert that the returned object is right on the parametric domain (take into account periodicity)
        self.assertTrue(test1.start.is_close(volmdlr.Point2D(0, 0.75 * math.pi)))
        self.assertTrue(test1.end.is_close(volmdlr.Point2D(0, 1.25 * math.pi)))

        arc2 = edges.Arc3D(volmdlr.Point3D(-0.169132244445, 0.06508125180570001, 0.627719515715),
        volmdlr.Point3D(-0.169169279223, 0.064939567779, 0.628073066814),
        volmdlr.Point3D(-0.169258691383, 0.064597504793, 0.628219515715))
        surface2 = surfaces.ToroidalSurface3D.load_from_file("surfaces/objects_toroidal_tests/surface.json")
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

    def test_face_from_contours3d(self):
        surface = surfaces.ToroidalSurface3D.load_from_file("faces/objects_toroidal_tests/surface_4.json")
        contour = wires.Contour3D.load_from_file("faces/objects_toroidal_tests/contour_4_0.json")
        face = surface.face_from_contours3d([contour])
        self.assertAlmostEqual(face.surface2d.area(), 0.07116351378250674, 4)

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


if __name__ == '__main__':
    unittest.main()
