import unittest
import volmdlr.step
class TestStep(unittest.TestCase):
    step = volmdlr.step.Step.from_file(filepath="step/test_conversion_factor.step")
    model = step.to_volume_model()

    def test_to_volume_model(self):
        conical_face = self.model.primitives[0].primitives[0]
        conical_surface = conical_face.surface3d
        not_found = True
        fullarc = None
        while not_found:
            for primitive in conical_face.outer_contour3d.primitives:
                if hasattr(primitive, "circle"):
                    fullarc = primitive
                    not_found = False
                    break
        self.assertAlmostEqual(conical_surface.semi_angle, 0.785398163397, places=8)
        self.assertAlmostEqual(fullarc.circle.radius, 0.007, places=3)


if __name__ == '__main__':
    unittest.main()
