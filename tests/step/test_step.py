import unittest
import volmdlr.step
class TestStep(unittest.TestCase):
    step = volmdlr.step.Step.from_file(filepath="step/test_conversion_factor.step")
    model = step.to_volume_model()
    def test_to_volume_model(self):
        conical_face = self.model.primitives[0].primitives[0]
        conical_surface = conical_face.surface3d
        fullarc = conical_face.outer_contour3d.primitives[0]
        self.assertAlmostEqual(conical_surface.semi_angle, 0.785398163397, places=8)
        self.assertAlmostEqual(fullarc.radius, 0.007, places=3)


if __name__ == '__main__':
    unittest.main()
