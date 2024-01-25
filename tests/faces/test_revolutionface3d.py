import unittest
import os
import volmdlr
from volmdlr.core import VolumeModel
from volmdlr.step import Step


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_revolution_tests')


class TestRevolutionFace3D(unittest.TestCase):

    def test_to_step(self):
        model = VolumeModel.from_json(os.path.join(folder, "revolutionface_export_test.json"))
        model.to_step(os.path.join(folder, "test_export.step"))
        step_import = Step.from_file(os.path.join(folder, "test_export.step"))
        model2 = step_import.to_volume_model()
        revolutionface = model2.primitives[0].primitives[0]
        self.assertTrue(revolutionface.outer_contour3d.is_ordered())
        self.assertAlmostEqual(revolutionface.surface2d.area(), 0.00738824 * volmdlr.TWO_PI, 6)


if __name__ == '__main__':
    unittest.main()
