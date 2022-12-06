"""
Tests for places faces
"""
import unittest

import dessia_common
import volmdlr
from volmdlr import edges


class TestPlaneFace3D(unittest.TestCase):
    face_with_3holes = dessia_common.DessiaObject.load_from_file('faces/face_with_3holes.json')
    face = dessia_common.DessiaObject.load_from_file('faces/face_to_cut_the_one_with_3holes.json')

    def test_area(self):
        self.assertAlmostEqual(self.face_with_3holes.area(), 0.12160000)

    def test_face_iside(self):
        face2 = self.face.frame_mapping(volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0), volmdlr.Vector3D(0.5, 0, 0),
                                        volmdlr.Vector3D(0, 0.5, 0), volmdlr.Vector3D(0, 0, 0.5)), 'old')
        self.assertEqual(self.face.face_inside(face2), True)
        self.assertEqual(face2.face_inside(self.face), False)

    def test_face_intersections_with_holes(self):
        face_intersections = self.face.face_intersections(self.face_with_3holes)
        self.assertEqual(len(face_intersections), 4)

    def test_line_intersections(self):
        line_inside_hole = edges.Line3D(volmdlr.Point3D(0.1, 0.0, -.3), volmdlr.Point3D(-0.1, 0.0, 0.3))
        line_inside_face = edges.Line3D(volmdlr.Point3D(-0.05, 0.0, -.3), volmdlr.Point3D(-0.1, 0.0, 0.3))
        self.assertEqual([], self.face_with_3holes.line_intersections(line_inside_hole))
        self.assertEqual(1, len(self.face_with_3holes.line_intersections(line_inside_face)))

    def test_divide_face(self):
        face_intersections = self.face.face_intersections(self.face_with_3holes)
        cutting_contours = self.face_with_3holes.get_face_cutting_contours(
            {(self.face, self.face_with_3holes): face_intersections})
        new_faces = self.face_with_3holes.divide_face(cutting_contours, False)
        self.assertEqual(len(new_faces), 2)


if __name__ == '__main__':
    unittest.main()
