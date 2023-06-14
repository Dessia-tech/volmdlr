import math
import unittest

import volmdlr
import volmdlr.edges as vme
import volmdlr.wires as vmw
from volmdlr import OXYZ, X3D, Y3D, Z3D, Point2D, Point3D
from volmdlr import surfaces


class TestSurface3D(unittest.TestCase):
    cylindrical_surface = surfaces.CylindricalSurface3D(OXYZ, radius=0.03)

    def test_contour3d_to_2d(self):
        primitives_cylinder = [vme.LineSegment3D(Point3D(0.03, 0, 0.003), Point3D(0.03, 0, 0.013)),
                        vme.FullArc3D(Point3D(0, 0, 0.013), Point3D(0.03, 0, 0.013), Z3D),
                        vme.LineSegment3D(Point3D(0.03, 0, 0.013), Point3D(0.03, 0, 0.003)),
                        vme.FullArc3D(Point3D(0, 0, 0.003), Point3D(0.03, 0, 0.003), Z3D)
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
