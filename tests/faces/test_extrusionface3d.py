import unittest
import volmdlr
from volmdlr import faces


class TestExtrusionFace3D(unittest.TestCase):
    face = faces.ExtrusionFace3D.load_from_file(
        "faces/objects_extrusion_tests/extrusionface_with_ellipse_test_boundingbox.json")
    def test_bounding_box(self):
        bbox = self.face.bounding_box
        self.assertAlmostEqual(bbox.volume(), 3.54136919512989e-08 )


if __name__ == '__main__':
    unittest.main()
