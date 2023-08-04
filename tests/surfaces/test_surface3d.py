import math
import unittest

import volmdlr
import volmdlr.edges as vme
import volmdlr.wires as vmw
from volmdlr import OXYZ, X3D, Y3D, Z3D, Point2D, Point3D
from volmdlr import surfaces, curves


class TestSurface3D(unittest.TestCase):
    cylindrical_surface = surfaces.CylindricalSurface3D(OXYZ, radius=0.03)

    def test_contour3d_to_2d(self):
        center1 = volmdlr.Point3D(0, 0, 0.013)
        start_end1 = Point3D(0.03, 0, 0.013)
        frame1 = volmdlr.Frame3D(center1, volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        circle1 = curves.Circle3D(frame1, center1.point_distance(start_end1))
        center2 = Point3D(0, 0, 0.003)
        start_end2 = Point3D(0.03, 0, 0.003)
        frame2 = volmdlr.Frame3D(center2, volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        circle2 = curves.Circle3D(frame2, start_end2.point_distance(center2))
        primitives_cylinder = [vme.LineSegment3D(Point3D(0.03, 0, 0.003), Point3D(0.03, 0, 0.013)),
                        vme.FullArc3D(circle1, Point3D(0.03, 0, 0.013)),
                        vme.LineSegment3D(Point3D(0.03, 0, 0.013), Point3D(0.03, 0, 0.003)),
                        vme.FullArc3D(circle2, Point3D(0.03, 0, 0.003))
                        ]
        contour_cylinder = vmw.Contour3D(primitives_cylinder)

        contour2d_cylinder = self.cylindrical_surface.contour3d_to_2d(contour_cylinder)

        area = contour2d_cylinder.area()
        linesegment2d = contour2d_cylinder.primitives[3]
        fullarc2d = contour2d_cylinder.primitives[2]

        self.assertEqual(area, 0.02*math.pi)
        self.assertEqual(fullarc2d.start, Point2D(volmdlr.TWO_PI, 0.003))
        self.assertEqual(fullarc2d.end, Point2D(0, 0.003))
        self.assertEqual(linesegment2d.start, Point2D(0, 0.003))
        self.assertEqual(linesegment2d.end, Point2D(0, 0.013))


if __name__ == '__main__':
    unittest.main()
