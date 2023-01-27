import math
import unittest

import dessia_common.core
import volmdlr
import volmdlr.edges as vme
import volmdlr.faces as vmf
import volmdlr.wires as vmw
from volmdlr import O3D, X3D, Y3D, Z3D, Point2D, Point3D
from volmdlr.models import conical_surfaces

class TestConicalSurface3D(unittest.TestCase):
    conical_surface = conical_surfaces.conical_surface1
    conical_surface2 = conical_surfaces.conical_surface2

    def test_arc3d_to_2d(self):
        arc1 = vme.Arc3D(volmdlr.Point3D(-1 / math.sqrt(2), 1 / math.sqrt(2), 1 / math.sqrt(3)),
                         volmdlr.Point3D(-1, 0, 1 / math.sqrt(3)),
                         volmdlr.Point3D(-1 / math.sqrt(2), -1 / math.sqrt(2), 1 / math.sqrt(3)))
        arc2 = vme.Arc3D(volmdlr.Point3D(0, -1, 1 / math.sqrt(3)),
                         volmdlr.Point3D(-1 / math.sqrt(2), 1 / math.sqrt(2), 1 / math.sqrt(3)),
                         volmdlr.Point3D(1, 0, 1 / math.sqrt(3)))
        test1 = self.conical_surface.arc3d_to_2d(arc3d=arc1)[0]
        test2 = self.conical_surface.arc3d_to_2d(arc3d=arc2)[0]

        # Assert that the returned object is an edges.LineSegment2D
        self.assertIsInstance(test1, vme.LineSegment2D)
        self.assertIsInstance(test2, vme.LineSegment2D)

        # Assert that the returned object is right on the parametric domain (take into account periodicity)
        self.assertEqual(test1.start, volmdlr.Point2D(0.75 * math.pi, 0.5773502691896258))
        self.assertEqual(test1.end, volmdlr.Point2D(1.25 * math.pi, 0.5773502691896258))
        self.assertEqual(test2.start, volmdlr.Point2D(-0.5 * math.pi, 0.5773502691896258))
        self.assertEqual(test2.end, volmdlr.Point2D(-2 * math.pi, 0.5773502691896258))

    def test_contour3d_to_2d(self):
        primitives_cone = [vme.LineSegment3D(Point3D(0, 0, 0.1), Point3D(0.035, 0, 0.0)),
                           vme.FullArc3D(O3D, Point3D(0.035, 0, 0), Z3D),
                           vme.LineSegment3D(Point3D(0.035, 0, 0.0), Point3D(0, 0, 0.1))]

        primitives_demi_cone = [primitives_cone[0],
                                vme.Arc3D(Point3D(0.035, 0, 0), Point3D(0, 0.035, 0), Point3D(-0.035, 0, 0)),
                                primitives_cone[2]
                                ]

        contour_cone = vmw.Contour3D(primitives_cone)
        contour2d_cone = self.conical_surface2.contour3d_to_2d(contour_cone)

        contour_demi_cone = vmw.Contour3D(primitives_demi_cone)
        contour2d_demi_cone = self.conical_surface2.contour3d_to_2d(contour_demi_cone)

        area_cone = contour2d_cone.area()
        area_demi_cone = contour2d_demi_cone.area()
        fullarc2d = contour2d_cone.primitives[0]
        linesegment2d_cone = contour2d_cone.primitives[1]

        # Assert that the returned object is an edges.LineSegment2D
        self.assertIsInstance(fullarc2d, vme.LineSegment2D)
        self.assertIsInstance(linesegment2d_cone, vme.LineSegment2D)

        self.assertEqual(area_cone, 0.2 * math.pi)
        self.assertEqual(area_demi_cone, 0.1 * math.pi)
        self.assertEqual(fullarc2d.start, Point2D(0, 0.1))
        self.assertEqual(fullarc2d.end, Point2D(2 * math.pi, 0.1))
        self.assertEqual(fullarc2d.length(), 2 * math.pi)
        self.assertEqual(linesegment2d_cone.start, Point2D(2 * math.pi, 0.1))
        self.assertEqual(linesegment2d_cone.end, Point2D(2 * math.pi, 0.0))

    def bsplinecurve3d_to_2d(self):
        conical_surface3 = conical_surfaces.conical_surface3
        control_points = [volmdlr.Point3D(-0.00235270234694772, 0.004075, 0.000502294734194974),
                          volmdlr.Point3D(-0.00158643061573795, 0.004075, 0.000281091139051792),
                          volmdlr.Point3D(-1.38558964783719e-06, 0.004075, -4.49804036433251e-06),
                          volmdlr.Point3D(0.00158690018976903, 0.004075, 0.000281226693398416),
                          volmdlr.Point3D(0.00235270234694773, 0.004075, 0.000502294734194975)]

        bspline_curve = vme.BSplineCurve3D(3, control_points, [4, 1, 4], knots=[0.0, 0.5, 1.0])
        bspline_curve2d = conical_surface3.bsplinecurve3d_to_2d(bspline_curve)
        bspline_curve3d = conical_surface3.bsplinecurve2d_to_3d(bspline_curve2d)
        original_length = bspline_curve.length()
        length_after_transformation = bspline_curve3d.length()
        point = bspline_curve.point_at_abscissa(0.5*original_length)
        point_test = bspline_curve3d.point_at_abscissa(0.5 * length_after_transformation)
        self.assertAlmostEqual(original_length, length_after_transformation, places=6)
        self.assertTrue(point.is_close(point_test, 1e-6))

    def test_face_from_contours(self):
        buggy_conical_surface = dessia_common.core.DessiaObject.load_from_file(
            'faces/objects_conical_tests/conical_surface1.json')
        buggy_contours3d1 = dessia_common.core.DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_from_contours1_0.json')
        buggy_contours3d2 = dessia_common.core.DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_from_contours1_1.json')

        conical_face = buggy_conical_surface.face_from_contours3d([buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.003769911184307754, 4)

        buggy_conical_surface = dessia_common.core.DessiaObject.load_from_file(
            'faces/objects_conical_tests/conical_surface3d_1.json')
        buggy_contours3d1 = dessia_common.core.DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_contour1.json')
        buggy_contours3d2 = dessia_common.core.DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_contour2.json')

        conical_face = buggy_conical_surface.face_from_contours3d([buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.004709873912166911, 4)

        buggy_conical_surface = dessia_common.core.DessiaObject.load_from_file(
            'faces/objects_conical_tests/conical_surface3d_2.json')
        buggy_contours3d1 = dessia_common.core.DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_contour3_.json')
        buggy_contours3d2 = dessia_common.core.DessiaObject.load_from_file(
            'faces/objects_conical_tests/face_contour4_.json')

        conical_face = buggy_conical_surface.face_from_contours3d([buggy_contours3d1, buggy_contours3d2])
        self.assertFalse(len(conical_face.surface2d.inner_contours))
        self.assertAlmostEqual(conical_face.area(), 0.05663155229526495, 4)


if __name__ == '__main__':
    unittest.main()
