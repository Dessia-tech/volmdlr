import unittest
import volmdlr
from volmdlr import wires, geometry, edges, curves


class TestGeometry(unittest.TestCase):
    def test_clockwise_interior_from_circle3d(self):
        point1 = volmdlr.Point3D(-0.210848709736, 0.624097020692, 0.6439887354620001)
        point2 = volmdlr.Point3D(-0.21084862123900003, 0.624, 0.6439887354620001)
        circle = curves.Circle3D(volmdlr.Frame3D(volmdlr.Point3D(-0.20750033288800002, 0.624, 0.548106302611),
                                                 volmdlr.Vector3D(0.0348994967025025, 0.0, -0.9993908270190957),
                                                 volmdlr.Vector3D(0.0, -0.9999999999999999, -0.0),
                                                 volmdlr.Vector3D(-0.9993908270190955, 0.0, -0.034899496702503335)),
                                 0.0959408774411)
        interior = geometry.clockwise_interior_from_circle3d(point1, point2, circle)
        arc = edges.Arc3D.from_3_points(point1, interior, point2)
        self.assertAlmostEqual(circle.radius, arc.circle.radius)  # add assertion here


if __name__ == '__main__':
    unittest.main()
