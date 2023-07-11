import unittest
from volmdlr.core import VolumeModel
from volmdlr.step import Step
from volmdlr import faces


class TestExtrusionFace3D(unittest.TestCase):
    face = faces.ExtrusionFace3D.load_from_file(
        "faces/objects_extrusion_tests/extrusionface_with_ellipse_test_boundingbox.json")

    def test_bounding_box(self):
        bbox = self.face.bounding_box
        self.assertAlmostEqual(bbox.volume(), 3.54136919512989e-08)

    def test_to_step(self):
        model = VolumeModel.load_from_file("faces/objects_extrusion_tests/extrusionface_export_test.json")
        model.to_step("faces/objects_extrusion_tests/test_export.step")
        step_import = Step.from_file("faces/objects_extrusion_tests/test_export.step")
        model2 = step_import.to_volume_model()
        extrusionface = model2.primitives[0].primitives[0]
        self.assertTrue(extrusionface.outer_contour3d.is_ordered())


if __name__ == '__main__':
    unittest.main()
