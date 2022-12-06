import unittest

import volmdlr
from volmdlr import edges


class TestArcEllipse3D(unittest.TestCase):
    arc_ellipse3d = edges.ArcEllipse3D(
        volmdlr.Point3D(-0.14790199457730885, 0.024999999999999994, 0.028874700719328038),
        volmdlr.Point3D(-0.09682458365486075, 0.11456439237413951, 0.07995211164177612),
        volmdlr.Point3D(0.0, 0.15, 0.17677669529663687),
        volmdlr.Point3D(0.0, 0.0, 0.17677669529663687),
        volmdlr.Vector3D(0.7071067811865475, 0.0, 0.7071067811865476))

    def test_discretization_points(self):
        discretization_points = self.arc_ellipse3d.discretization_points(number_points=3)
        self.assertEqual(len(discretization_points), 3)
        expected_points = [volmdlr.Point3D(-0.14790199457730882, 0.024999999999875906, 0.028874700719328045),
                           volmdlr.Point3D(-0.09682458365510088, 0.11456439237385536, 0.07995211164153597),
                           volmdlr.Point3D(0.0, 0.15, 0.17677669529663684)]
        for expected_point, point in zip(expected_points, discretization_points):
            self.assertEqual(expected_point, point)


if __name__ == '__main__':
    unittest.main()
