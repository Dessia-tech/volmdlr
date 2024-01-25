import os
import unittest
import volmdlr
from volmdlr.wires import Contour3D
from volmdlr.step import Step
from volmdlr import edges, core
from volmdlr.models.contours import contour3d


folder = os.path.dirname(os.path.realpath(__file__))


class TestContour3D(unittest.TestCase):

    def test_order_contour(self):
        contour_to_order = Contour3D.from_json(os.path.join(folder, 'contour_order.json'))
        self.assertFalse(contour_to_order.is_ordered())
        contour_to_order.order_contour()
        self.assertTrue(contour_to_order.is_ordered())

    def test_merge_with(self):
        contour1_to_merge = Contour3D.from_json(os.path.join(folder, 'contour3d_merge_with1.json'))
        contour2_to_merge = Contour3D.from_json(os.path.join(folder, 'contour3d_merge_with2.json'))
        expected_contour1 = Contour3D.from_json(os.path.join(folder, 'expected_contour_merge_with1.json'))
        expected_contour2 = Contour3D.from_json(os.path.join(folder, 'expected_contour_merge_with2.json'))
        merged_contours = contour1_to_merge.merge_with(contour2_to_merge)
        self.assertEqual(merged_contours[0], expected_contour1)
        self.assertEqual(merged_contours[1], expected_contour2)
        contour1 = Contour3D.from_json(os.path.join(folder, 'contour1_merge_bug.json'))
        contour2 = Contour3D.from_json(os.path.join(folder, 'contour2_merge_bug.json'))
        merged_contour1_contour2 = contour1.merge_with(contour2)
        merged_contour2_contour1 = contour2.merge_with(contour1)
        self.assertEqual(len(merged_contour1_contour2), len(merged_contour2_contour1))
        self.assertEqual(merged_contour1_contour2[0], merged_contour2_contour1[0])

    def test_is_sharing_primitives_with(self):
        contour1_sharing_primitives = Contour3D.from_json(
            os.path.join(folder, 'contour3d_sharing_primitives1.json'))
        contour2_sharing_primitives = Contour3D.from_json(
            os.path.join(folder, 'contour3d_sharing_primitives2.json'))
        self.assertTrue(contour1_sharing_primitives.is_sharing_primitives_with(contour2_sharing_primitives))

    def test_from_step(self):
        step = Step.from_file(filepath=os.path.join(folder, "contour_with_repeated_edge_in_contour3d.step"))
        model = step.to_volume_model()
        face = model.primitives[0].primitives[0]
        self.assertEqual(len(face.outer_contour3d.primitives), 5)
        self.assertTrue(face.outer_contour3d.is_ordered())

        # todo: refactor SphericalSuface3D repair periodicity
        # step = Step.from_file(filepath="wires/sphere_with_singularity.step")
        # model = step.to_volume_model()
        # self.assertTrue(model)

        step = Step.from_file(filepath=os.path.join(folder, "contour_with_repeated_edge_in_contour3d.step"))
        model = step.to_volume_model()
        face = model.primitives[0].primitives[0]
        self.assertEqual(len(face.outer_contour3d.primitives), 5)
        self.assertTrue(face.outer_contour3d.is_ordered())

        arguments = ["", ["#2518728", "#2518729"]]
        primitives = core.VolumeModel.from_json(
                os.path.join(folder, "strange_contour_from_step_primitives.json")).primitives
        object_dict = {2518728: primitives[0], 2518729: primitives[1]}

        contour = Contour3D.from_step(arguments, object_dict)

        self.assertFalse(contour)

        arguments = ["''", ['#13203123', '#13203124', '#13203125', '#13203126', '#13203127',
                            '#13203128', '#13203129', '#13203130', '#13203131', '#13203132']]
        primitives = Contour3D.from_json(
            os.path.join(folder, "edge_loop_with_small_edges_and_gaps.json")).primitives
        object_dict = {int(arg[1:]): edge for arg, edge in zip(arguments[1], primitives)}
        contour = Contour3D.from_step(arguments, object_dict)
        self.assertFalse(contour.is_ordered())
        self.assertTrue(contour.is_ordered(5e-6))

    def test_edge_intersections(self):
        points = [volmdlr.Point3D(1.2918566581549966, 2.3839907440191492, 0.5678759590090421),
                  volmdlr.Point3D(1.2067665579541171, -1.246879774203074, -0.4359328108960321),
                  volmdlr.Point3D(-1.2905737351068276, -5.961765089244547, -0.9872550297481824),
                  volmdlr.Point3D(7.33260591629263, -4.272128323147327, -0.4240427743824422),
                  volmdlr.Point3D(7.115095014105684, 0.40888620982702983, 1.1362954032756774),
                  volmdlr.Point3D(-3.0, 1.022248896290622, 0.5746069851843745),
                  volmdlr.Point3D(2.739350840642852, -5.869347626045908, -0.7880999427201254)]
        bspline = edges.BSplineCurve3D.from_points_interpolation(points, 3)
        edge_intersections = contour3d.edge_intersections(bspline, abs_tol=1e-6)
        expected_results = [volmdlr.Point3D(1.2918566581549966, 2.3839907440191492, 0.5678759590090421),
                            volmdlr. Point3D(1.206766559907027, -1.2468797685507946, -0.4359328046874991),
                            volmdlr. Point3D(-3.0, 1.0222488954206392, 0.5746069850600913),
                            volmdlr. Point3D(-1.2905737300311637, -5.9617650927233345, -0.9872550300602736),
                            volmdlr. Point3D(7.332606025327417, -4.272128068303522, -0.42404268513457977),
                            volmdlr. Point3D(7.115095100684387, 0.408886063311686, 1.136295354437229)]
        for intersection, expected_result in zip(edge_intersections, expected_results):
            self.assertTrue(intersection.is_close(expected_result))


if __name__ == '__main__':
    unittest.main()
