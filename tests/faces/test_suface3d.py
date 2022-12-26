import unittest
import math

import volmdlr.faces as vmf
import volmdlr.edges as vme
import volmdlr.wires as vmw
from volmdlr import Point3D, OXYZ, Z3D, Point2D


class TestSurface3D(unittest.TestCase):
    cylindrical_surface = vmf.CylindricalSurface3D(OXYZ, radius=0.03)

    def test_contour3d_to_2d(self):
        primitives = [vme.LineSegment3D(Point3D(0.03, 0, 0.003), Point3D(0.03, 0, 0.013)),
                        vme.FullArc3D(Point3D(0, 0, 0.013), Point3D(0.03, 0, 0.013), Z3D),
                        vme.LineSegment3D(Point3D(0.03, 0, 0.013), Point3D(0.03, 0, 0.003)),
                        vme.FullArc3D(Point3D(0, 0, 0.003), Point3D(0.03, 0, 0.003), Z3D)
                        ]
        contour_cylinder = vmw.Contour3D(primitives)

        contour2d_cylinder = self.cylindrical_surface.contour3d_to_2d(contour_cylinder)

        area = contour2d_cylinder.area()
        fullarc2d = contour2d_cylinder.primitives[3]
        linesegment2d = contour2d_cylinder.primitives[2]

        self.assertEqual(area, 0.02*math.pi)
        self.assertEqual(fullarc2d.start, Point2D(2*math.pi, 0.003))
        self.assertEqual(fullarc2d.end, Point2D(0, 0.003))
        self.assertEqual(linesegment2d.start, Point2D(2*math.pi, 0.013))
        self.assertEqual(linesegment2d.end, Point2D(2*math.pi, 0.003))


if __name__ == '__main__':
    unittest.main()
