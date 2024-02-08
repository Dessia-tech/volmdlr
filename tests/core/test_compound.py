import unittest
import os

from volmdlr import step

folder = os.path.dirname(os.path.realpath(__file__))


class TestCompound(unittest.TestCase):
    def test_to_step(self):
        step_obj = step.Step.from_file(filepath=os.path.join(folder, "compound_geometric_set.step"))
        model = step_obj.to_volume_model()
        model.to_step("compound_export.stp")
        step_obj = step.Step.from_file(filepath=os.path.join(folder, "compound_export.stp"))
        model = step_obj.to_volume_model()
        self.assertEqual(len(model.primitives), 1)
        self.assertEqual(len(model.primitives[0].primitives), 3)
        self.assertEqual(model.primitives[0].compound_type, "geometric_curve_set")


if __name__ == '__main__':
    unittest.main()
