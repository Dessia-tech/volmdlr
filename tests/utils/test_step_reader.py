import unittest
import volmdlr
from volmdlr.curves import Line3D
from volmdlr.utils import step_reader


class TestStepReader(unittest.TestCase):
    def test_trimmed_curve(self):
        point1 = volmdlr.Point3D(2.2, -1.0, 0.05)
        point2 = volmdlr.Point3D(2.201, -1.0, 0.05)
        object_dict = {326: Line3D(point1, point2)}
        arguments = ["'n° 6246'", 326, '(PARAMETER_VALUE(0.0))', '(PARAMETER_VALUE(20.0))', '.T.', '.PARAMETER.']
        linesegment = step_reader.trimmed_curve(arguments, object_dict, length_conversion_factor=0.001)
        self.assertAlmostEqual(linesegment.length(), 0.02)


if __name__ == '__main__':
    unittest.main()
