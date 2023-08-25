
import unittest
import volmdlr
from volmdlr.wires import Contour3D
from volmdlr.step import Step
from volmdlr import edges
from volmdlr.models.contours import contour3d


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
        self.assertEqual(len(merged_contour1_contour2), len(merged_contour2_contour1))
        self.assertEqual(merged_contour1_contour2[0], merged_contour2_contour1[0])

    def test_is_sharing_primitives_with(self):
        contour1_sharing_primitives = Contour3D.load_from_file('wires/contour3d_sharing_primitives1.json')
        contour2_sharing_primitives = Contour3D.load_from_file('wires/contour3d_sharing_primitives2.json')

        self.assertTrue(contour1_sharing_primitives.is_sharing_primitives_with(contour2_sharing_primitives))

    def test_from_step(self):
        step = Step.from_file(filepath="wires/contour_with_repeated_edge_in_contour3d.step")
        model = step.to_volume_model()
        face = model.primitives[0].primitives[0]
        self.assertEqual(len(face.outer_contour3d.primitives), 4)

        # todo: refactor SphericalSuface3D repair periodicity
        # step = Step.from_file(filepath="wires/sphere_with_singularity.step")
        # model = step.to_volume_model()
        # self.assertTrue(model)

        step = Step.from_file(filepath="wires/contour_with_repeated_edge_in_contour3d.step")
        model = step.to_volume_model()
        face = model.primitives[0].primitives[0]
        self.assertEqual(len(face.outer_contour3d.primitives), 4)

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
                            volmdlr.Point3D(1.2067665517702832, -1.2468797560646305, -0.4359328005283978),
                            volmdlr.Point3D(-3.000003493713931, 1.022247107552729, 0.5746061300812078),
                            volmdlr.Point3D(-1.2905737300371376, -5.961765092719032, -0.9872550300598877),
                            volmdlr.Point3D(7.332606025292464, -4.272128068376155, -0.42404268515985877),
                            volmdlr.Point3D(7.115095070785982, 0.40888610210140913, 1.1362953668149238)]
        for intersection, expected_result in zip(edge_intersections, expected_results):
            self.assertTrue(intersection.is_close(expected_result))


if __name__ == '__main__':
    unittest.main()
