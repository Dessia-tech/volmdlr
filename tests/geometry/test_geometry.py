import unittest
import volmdlr
from volmdlr import wires, geometry, edges


class TestGeometry(unittest.TestCase):
    def test_clockwise_interior_from_circle3d(self):
        point1 = volmdlr.Point3D.load_from_file("geometry/clockwise_interior_from_circle3d_start.json")
        point2 = volmdlr.Point3D.load_from_file("geometry/clockwise_interior_from_circle3d_end.json")
        circle = wires.Circle3D.load_from_file("geometry/clockwise_interior_from_circle3d_circle.json")
        interior = geometry.clockwise_interior_from_circle3d(point1, point2, circle)
        arc = edges.Arc3D(point1, interior, point2)
        self.assertAlmostEqual(circle.radius, arc.radius)  # add assertion here


if __name__ == '__main__':
    unittest.main()
