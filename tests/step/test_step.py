import unittest
import os
import volmdlr.step
from volmdlr.utils import step_reader


folder = os.path.dirname(os.path.realpath(__file__))


class TestStep(unittest.TestCase):

    def test_to_volume_model(self):
        step = volmdlr.step.Step.from_file(filepath=os.path.join(folder, "test_conversion_factor.step"))
        model = step.to_volume_model()
        conical_face = model.primitives[0].primitives[0]
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

    def test_read_lines(self):
        step = volmdlr.step.Step.from_file(filepath=os.path.join(folder, "test_names.step"))
        model = step.to_volume_model()
        self.assertEqual(model.primitives[0].name, "'cube assembly =,v1'")
        self.assertEqual(model.primitives[0].primitives[0].name, "Part 2")
        self.assertEqual(model.primitives[0].primitives[1].name, "Part 1")

        step = volmdlr.step.Step.from_file(filepath=os.path.join(folder, "test_step_read_lines.stp"))
        self.assertEqual(step.functions[3].name, "REALLY_CHALLENGING_ENTITY")
        self.assertEqual(len(step.functions[3].arg), 5)
        self.assertEqual(step.functions[3].arg[0], "'nom, spécial'")
        self.assertEqual(step.functions[3].arg[1], ['#1', '#2'])
        self.assertEqual(step.functions[3].arg[2], "((#3,#4),(#5,#6))")
        self.assertEqual(step.functions[3].arg[3], "'nom, d'un entité'")
        self.assertEqual(step.functions[3].arg[4], "(PARAMETER_VALUE(20.0))")
        self.assertEqual(step.functions[4].name, "TRIMMED_CURVE")
        self.assertEqual(len(step.functions[4].arg), 6)
        self.assertEqual(step.functions[4].arg[2], "(#9,PARAMETER_VALUE(0.))")
        self.assertEqual(step.functions[4].arg[3], "(#10,PARAMETER_VALUE(13.))")



    def test_split_arguments_special(self):
        function_args = "'nom, spécial', (#1, #2), ((#3, #4), (#5, #6)), (PARAMETER_VALUE(20.0)));"
        arguments = step_reader.step_split_arguments_special(function_args)
        self.assertEqual(len(arguments), 4)

    def test_shape_representation(self):
        step = volmdlr.step.Step.from_file(filepath=os.path.join(folder, "cone.step"))
        model = step.to_volume_model()
        self.assertEqual(len(model.primitives), 1)

    def test_create_connections(self):
        step = volmdlr.step.Step.from_file(filepath=os.path.join(folder, "test_wireframe.step"))
        _ = step.to_volume_model()
        self.assertEqual(len(step.functions[327014].arg), 9)
        self.assertEqual(step.functions[327014].arg[-1], '.UNSPECIFIED.')


if __name__ == '__main__':
    unittest.main()
