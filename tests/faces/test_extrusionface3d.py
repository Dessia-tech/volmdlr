import unittest
from volmdlr.core import VolumeModel
from volmdlr.step import Step
from volmdlr import faces, surfaces, wires
import os

folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'objects_extrusion_tests')


class TestExtrusionFace3D(unittest.TestCase):
    face = faces.ExtrusionFace3D.from_json(
        os.path.join(folder, "extrusionface_with_ellipse_test_boundingbox.json"))

    def test_bounding_box(self):
        bbox = self.face.bounding_box
        self.assertAlmostEqual(bbox.volume(), 4.078418559129855e-08)

    def test_from_contours3d(self):
        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusionsurface_inner_contour.json"))
        contour_0 = wires.Contour3D.from_json(
            os.path.join(folder, "extrusionsurface_inner_contour_contour_0.json"))
        contour_1 = wires.Contour3D.from_json(
            os.path.join(folder, "extrusionsurface_inner_contour_contour_1.json"))

        extrusionface = faces.ExtrusionFace3D.from_contours3d(surface, [contour_0, contour_1])
        self.assertAlmostEqual(extrusionface.surface2d.area(), 1.5451364316425414e-06, 8)
        self.assertTrue(extrusionface.surface2d.outer_contour.is_ordered())
        self.assertTrue(extrusionface.surface2d.inner_contours[0].is_ordered())

        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "periodic_extrusionsurface_small_bsplinecurve.json"))
        contour = wires.Contour3D.from_json(
            os.path.join(folder, "periodic_extrusionsurface_small_bsplinecurve_contour.json"))

        extrusionface = faces.ExtrusionFace3D.from_contours3d(surface, [contour])
        self.assertAlmostEqual(extrusionface.surface2d.area(),  6.186940971694699e-08, 9)
        
        surface = surfaces.ExtrusionSurface3D.from_json(
            os.path.join(folder, "extrusionsurface_fullarcellipse.json"))
        contour_0 = wires.Contour3D.from_json(
            os.path.join(folder, "extrusionsurface_fullarcellipse_contour.json"))
        extrusionface = faces.ExtrusionFace3D.from_contours3d(surface, [contour_0])
        self.assertAlmostEqual(extrusionface.surface2d.area(), 0.0015836412348717846, 4)

        self.assertTrue(extrusionface.surface2d.outer_contour.is_ordered())

    def test_to_step(self):
        model = VolumeModel.from_json(os.path.join(folder, "extrusionface_export_test.json"))
        model.to_step(os.path.join(folder, "test_export.step"))
        step_import = Step.from_file(os.path.join(folder, "test_export.step"))
        model2 = step_import.to_volume_model()
        extrusionface = model2.primitives[0].primitives[0]
        self.assertTrue(extrusionface.outer_contour3d.is_ordered())


if __name__ == '__main__':
    unittest.main()
