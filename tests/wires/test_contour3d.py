import unittest
from volmdlr.step import Step
from volmdlr.wires import Contour3D


class TestContour3D(unittest.TestCase):

    def test_order_contour(self):
        contour_to_order = Contour3D.load_from_file('wires/contour_order.json')
        self.assertFalse(contour_to_order.is_ordered())
        contour_to_order.order_contour()
        self.assertTrue(contour_to_order.is_ordered())

    def test_merge_with(self):
        contour1_to_merge = Contour3D.load_from_file('wires/contour3d_merge_with1.json')
        contour2_to_merge = Contour3D.load_from_file('wires/contour3d_merge_with2.json')
        expected_contour1 = Contour3D.load_from_file('wires/expected_contour_merge_with1.json')
        expected_contour2 = Contour3D.load_from_file('wires/expected_contour_merge_with2.json')
        merged_contours = contour1_to_merge.merge_with(contour2_to_merge)
        self.assertEqual(merged_contours[0], expected_contour1)
        self.assertEqual(merged_contours[1], expected_contour2)
        contour1 = Contour3D.load_from_file('wires/contour1_merge_bug.json')
        contour2 = Contour3D.load_from_file('wires/contour2_merge_bug.json')
        merged_contour1_contour2 = contour1.merge_with(contour2)
        merged_contour2_contour1 = contour2.merge_with(contour1)
        self.assertEqual(merged_contour1_contour2, merged_contour2_contour1)

    def test_is_sharing_primitives_with(self):
        contour1_sharing_primitives = Contour3D.load_from_file('wires/contour3d_sharing_primitives1.json')
        contour2_sharing_primitives = Contour3D.load_from_file('wires/contour3d_sharing_primitives2.json')

        self.assertTrue(contour1_sharing_primitives.is_sharing_primitives_with(contour2_sharing_primitives))

    def test_from_step(self):
        step = Step.from_file(filepath="wires/contour_with_repeated_edge_in_contour3d.step")
        model = step.to_volume_model()
        face = model.primitives[0].primitives[0]
        self.assertEqual(len(face.outer_contour3d.primitives), 4)


if __name__ == '__main__':
    unittest.main()
