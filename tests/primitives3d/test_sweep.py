import unittest

import volmdlr
from volmdlr import wires, curves, edges, faces
from volmdlr.primitives3d import Sweep


class TestSweep(unittest.TestCase):
    def test_init(self):
        path = wires.Wire3D(primitives=[
            edges.Arc3D.from_3_points(volmdlr.Point3D(0, -1, 0), volmdlr.Point3D(1, 0.0, 0), volmdlr.Point3D(0, 1, 0)),
            edges.LineSegment3D(volmdlr.Point3D(0, 1, 0), volmdlr.Point3D(2, 1, 2)),
            edges.LineSegment3D(volmdlr.Point3D(2, 1, 2), volmdlr.Point3D(3, -1, 1))])
        section = wires.Contour2D([curves.Circle2D(volmdlr.O2D, 0.05)])
        sweep = Sweep(section, path)
        self.assertEqual(len(sweep.faces), 5)
        for face, expected_face_class in zip(sweep.faces, [faces.PlaneFace3D, faces.ToroidalFace3D,
                                                           faces.CylindricalFace3D, faces.CylindricalFace3D,
                                                           faces.PlaneFace3D]):
            self.assertTrue(isinstance(face, expected_face_class))


if __name__ == '__main__':
    unittest.main()
