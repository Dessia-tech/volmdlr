import unittest

import volmdlr
from volmdlr.curves import (Line2D, Line3D, Circle2D, Circle3D, Ellipse2D, Ellipse3D, Hyperbola3D,
                            Parabola3D)
from volmdlr.from_ocp import (line2d_from_ocp, line3d_from_ocp, circle2d_from_ocp, circle3d_from_ocp,
                              ellipse2d_from_ocp, ellipse3d_from_ocp, hyperbola3d_from_ocp,
                              parabola3d_from_ocp)
from OCP.gp import gp_Pnt, gp_Pnt2d, gp_Dir, gp_Dir2d, gp_Ax1, gp_Ax2, gp_Ax2d, gp_Ax22d
from OCP.Geom import Geom_Line, Geom_Circle, Geom_Ellipse, Geom_Hyperbola, Geom_Parabola
from OCP.Geom2d import Geom2d_Line, Geom2d_Circle, Geom2d_Ellipse, Geom2d_BSplineCurve


class TestOCPConversion(unittest.TestCase):
    ocp_line2d = Geom2d_Line(gp_Ax2d(gp_Pnt2d(0, 0), gp_Dir2d(1, 0)))
    ocp_line3d = Geom_Line(gp_Ax1(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)))

    ocp_circle2d = Geom2d_Circle(gp_Ax22d(gp_Pnt2d(0, 0), gp_Dir2d(1, 0)), 1)
    ocp_circle3d = Geom_Circle(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), 1)

    ocp_ellipse2d = Geom2d_Ellipse(gp_Ax22d(gp_Pnt2d(0, 0), gp_Dir2d(1, 0)), 1, 0.5)
    ocp_ellipse3d = Geom_Ellipse(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), 1, 0.5)

    ocp_hyperbola3d = Geom_Hyperbola(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), 1, 0.5)

    ocp_parabola3d = Geom_Parabola(gp_Ax2(gp_Pnt(0, 0, 0), gp_Dir(1, 0, 0)), 1)

    def test_line2d_from_ocp(self):
        result = line2d_from_ocp(Line2D, self.ocp_line2d)
        self.assertIsInstance(result, Line2D)
        self.assertEqual(result.point1, volmdlr.Point2D(0, 0))
        self.assertEqual(result.point2, volmdlr.Point2D(1, 0))

    def test_line3d_from_ocp(self):
        result = line3d_from_ocp(Line3D, self.ocp_line3d)
        self.assertIsInstance(result, Line3D)
        self.assertEqual(result.point1, volmdlr.Point3D(0, 0, 0))
        self.assertEqual(result.point2, volmdlr.Point3D(1, 0, 0))

    def test_circle2d_from_ocp(self):
        result = circle2d_from_ocp(Circle2D, self.ocp_circle2d)
        self.assertIsInstance(result, Circle2D)
        self.assertEqual(result.center, volmdlr.Point2D(0, 0))
        self.assertEqual(result.radius, 1)

    def test_circle3d_from_ocp(self):
        result = circle3d_from_ocp(Circle3D, self.ocp_circle3d)
        self.assertIsInstance(result, Circle3D)
        self.assertEqual(result.center, volmdlr.Point3D(0, 0, 0))
        self.assertEqual(result.normal, volmdlr.Vector3D(1, 0, 0))
        self.assertEqual(result.radius, 1)

    def test_ellipse2d_from_ocp(self):
        result = ellipse2d_from_ocp(Ellipse2D, self.ocp_ellipse2d)
        self.assertIsInstance(result, Ellipse2D)
        self.assertTupleEqual(tuple(result.center),
                              (self.ocp_ellipse2d.Location().X(), self.ocp_ellipse2d.Location().Y()))
        self.assertTupleEqual(tuple(result.focus1), (self.ocp_ellipse2d.Focus1().X(), self.ocp_ellipse2d.Focus1().Y()))
        self.assertTupleEqual(tuple(result.focus2), (self.ocp_ellipse2d.Focus2().X(), self.ocp_ellipse2d.Focus2().Y()))
        self.assertEqual(result.major_axis, self.ocp_ellipse2d.MajorRadius())
        self.assertEqual(result.minor_axis, self.ocp_ellipse2d.MinorRadius())

    def test_ellipse3d_from_ocp(self):
        result = ellipse3d_from_ocp(Ellipse3D, self.ocp_ellipse3d)
        self.assertIsInstance(result, Ellipse3D)
        self.assertTupleEqual(tuple(result.center), (self.ocp_ellipse3d.Position().Location().X(),
                                                     self.ocp_ellipse3d.Position().Location().Y(),
                                                     self.ocp_ellipse3d.Position().Location().Z()))
        self.assertTupleEqual(tuple(result.normal),
                              (self.ocp_ellipse3d.Position().Direction().X(),
                               self.ocp_ellipse3d.Position().Direction().Y(),
                               self.ocp_ellipse3d.Position().Direction().Z()))
        self.assertTupleEqual(tuple(result.focus1), (self.ocp_ellipse3d.Focus1().X(),
                                                     self.ocp_ellipse3d.Focus1().Y(),
                                                     self.ocp_ellipse3d.Focus1().Z()))
        self.assertTupleEqual(tuple(result.focus2), (self.ocp_ellipse3d.Focus2().X(),
                                                     self.ocp_ellipse3d.Focus2().Y(),
                                                     self.ocp_ellipse3d.Focus2().Z()))
        self.assertEqual(result.major_axis, self.ocp_ellipse3d.MajorRadius())
        self.assertEqual(result.minor_axis, self.ocp_ellipse3d.MinorRadius())

    def test_hyperbola3d_from_ocp(self):
        result = hyperbola3d_from_ocp(Hyperbola3D, self.ocp_hyperbola3d)
        self.assertIsInstance(result, Hyperbola3D)
        self.assertTupleEqual(tuple(result.focus), (self.ocp_hyperbola3d.Focus1().X(),
                                                    self.ocp_hyperbola3d.Focus1().Y(),
                                                    self.ocp_hyperbola3d.Focus1().Z()))
        self.assertEqual(result.semi_major_axis, self.ocp_hyperbola3d.MajorRadius())
        self.assertEqual(result.semi_minor_axis, self.ocp_hyperbola3d.MinorRadius())

    def test_parabola3d_from_ocp(self):
        result = parabola3d_from_ocp(Parabola3D, self.ocp_parabola3d)
        self.assertIsInstance(result, Parabola3D)
        self.assertTupleEqual(tuple(result.focus), (self.ocp_parabola3d.Focus().X(),
                                                    self.ocp_parabola3d.Focus().Y(),
                                                    self.ocp_parabola3d.Focus().Z()))


if __name__ == '__main__':
    unittest.main()
