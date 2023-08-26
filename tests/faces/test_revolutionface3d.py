import unittest

import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.step import Step


class TestRevolutionFace3D(unittest.TestCase):

    def test_to_step(self):
        model = VolumeModel.load_from_file("faces/objects_revolution_tests/revolutionface_export_test.json")
        model.to_step("faces/objects_revolution_tests/test_export.step")
        step_import = Step.from_file("faces/objects_revolution_tests/test_export.step")
        model2 = step_import.to_volume_model()
        revolutionface = model2.primitives[0].primitives[0]
        self.assertTrue(revolutionface.outer_contour3d.is_ordered())
        self.assertAlmostEqual(revolutionface.surface2d.area(), 0.00738824 * volmdlr.TWO_PI, 6)


if __name__ == '__main__':
    unittest.main()
