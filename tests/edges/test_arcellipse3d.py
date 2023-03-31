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
    arc_ellipse3d_2 = edges.ArcEllipse3D(
        volmdlr.Point3D(0.029550000005800003, -2.071669979697e-17, 0.0039499999994),
        volmdlr.Point3D(0.03070692821993751, -0.0027930717852626134, 0.0027930717852624955),
        volmdlr.Point3D(0.0335000000052, -0.0039499999994, -2.504633374077e-33),
        volmdlr.Point3D(0.0335000000052, 0.0, 3.095745883019e-18),
        volmdlr.Vector3D(-0.7071067811865476, -3.605549748844693e-15, 0.7071067811865476)
    )


    def test_discretization_points(self):
        discretization_points = self.arc_ellipse3d.discretization_points(number_points=3)
        self.assertEqual(len(discretization_points), 3)
        expected_points = [volmdlr.Point3D(-0.14790199457730882, 0.024999999999875906, 0.028874700719328045),
                           volmdlr.Point3D(-0.09682458365510088, 0.11456439237385536, 0.07995211164153597),
                           volmdlr.Point3D(0.0, 0.15, 0.17677669529663684)]
        for expected_point, point in zip(expected_points, discretization_points):
            self.assertTrue(expected_point.is_close(point))

        expected_points = [volmdlr.Point3D(0.029550000005800003, 0.0, 0.003949999999400002),
        volmdlr.Point3D(0.029610009381392667, -0.0006859103016802056, 0.003889990623807338),
        volmdlr.Point3D(0.029788214153659485, -0.001350979565931198, 0.0037117858515405226),
        volmdlr.Point3D(0.03007919966077109, -0.0019749999997000166, 0.0034208003444289163),
        volmdlr.Point3D(0.030474124455339677, -0.002539011057876173, 0.0030258755498603343),
        volmdlr.Point3D(0.030960988947323857, -0.003025875549860349, 0.0025390110578761556),
        volmdlr.Point3D(0.031525000005500015, -0.0034208003444289267, 0.0019749999996999975),
        volmdlr.Point3D(0.03214902043926884, -0.0037117858515405286, 0.0013509795659311767),
        volmdlr.Point3D(0.03281408970351983, -0.00388999062380734, 0.0006859103016801831),
        volmdlr.Point3D(0.0335000000052, -0.0039499999993999994, 3.095745883019e-18)]
        discretization_points = self.arc_ellipse3d_2.discretization_points(number_points=10)
        self.assertEqual(len(discretization_points), 10)
        for expected_point, point in zip(expected_points, discretization_points):
            self.assertTrue(expected_point.is_close(point))


if __name__ == '__main__':
    unittest.main()
