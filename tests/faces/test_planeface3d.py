"""
Tests for places faces
"""
import math
import unittest

import dessia_common.core as dc

import volmdlr
from volmdlr import edges, wires, faces


class TestPlaneFace3D(unittest.TestCase):
    face_with_3holes = dc.DessiaObject.load_from_file('faces/face_with_3holes.json')
    face = dc.DessiaObject.load_from_file('faces/face_to_cut_the_one_with_3holes.json')
    plane_face_cylindricalface_intersec = dc.DessiaObject.load_from_file(
        'faces/plane_face_cylindrical_face_intersec.json')

    def test_area(self):
        self.assertAlmostEqual(self.face_with_3holes.area(), 0.12160000)

    def test_face_inside(self):
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

    def test_cylindricalface_intersections(self):
        R = 0.15
        cylindricalsurface = volmdlr.faces.CylindricalSurface3D(volmdlr.OXYZ, R)
        face = cylindricalsurface.rectangular_cut(0, volmdlr.TWO_PI, -.25, .25)
        """ ========== CIRCLE3D ========="""
        plane_face_3 = self.plane_face_cylindricalface_intersec.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 2)
        face_intersections = plane_face_3.face_intersections(face)
        self.assertEqual(len(face_intersections), 1)
        self.assertIsInstance(face_intersections[0], wires.Circle3D)
        self.assertEqual(face_intersections[0].center, volmdlr.O3D)
        self.assertEqual(face_intersections[0].radius, 0.15)
        """ ========== FULL ELLIPSE3D ========="""
        plane_face_3 = self.plane_face_cylindricalface_intersec.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 4)
        face_intersections = plane_face_3.face_intersections(face)
        self.assertEqual(len(face_intersections), 1)
        self.assertIsInstance(face_intersections[0].primitives[0], wires.Ellipse3D)
        self.assertEqual(face_intersections[0].primitives[0].center, volmdlr.O3D)
        self.assertAlmostEqual(face_intersections[0].primitives[0].major_axis, 0.21213203435596426)
        self.assertTrue(face_intersections[0].primitives[0].major_dir.is_close(volmdlr.Vector3D(0, 0.7071067811865475,
                                                                                                -0.7071067811865475)))
        """ ========== THREE ARC ELLIPSES ========="""
        plane_face_3 = self.plane_face_cylindricalface_intersec.rotation(volmdlr.O3D, volmdlr.X3D, math.pi / 7)
        face_intersections = plane_face_3.face_intersections(face)
        self.assertEqual(len(face_intersections), 3)
        for inter in face_intersections:
            self.assertIsInstance(inter.primitives[0], edges.ArcEllipse3D)
        self.assertEqual(face_intersections[0].primitives[0].center, volmdlr.O3D)
        self.assertEqual(face_intersections[0].primitives[0].major_dir, volmdlr.Point3D(
            2.6567716615652136e-17, 0.4338837391180807, -0.9009688679021675))
        self.assertAlmostEqual(face_intersections[0].primitives[0].Gradius, 0.3457147306439571)
        list_expected_points = [[volmdlr.Point3D(0.08947272158306664, 0.12039365470206077, -0.25),
                                 volmdlr.Point3D(0.13401661881546628, 0.0673761521702598, -0.13990802160005056),
                                 volmdlr.Point3D(0.15, 0.0, 0.0)],
                                [volmdlr.Point3D(0.15, 0.0, 0.0),
                                 volmdlr.Point3D(0.1340166188154663, -0.0673761521702598, 0.13990802160005056),
                                 volmdlr.Point3D(0.08947272158306664, -0.12039365470206077, 0.25)],
                                [volmdlr.Point3D(-0.08947272158306664, -0.12039365470206077, 0.25),
                                 volmdlr.Point3D(-0.15, 0.0, 0.0),
                                 volmdlr.Point3D(-0.08947272158306664, 0.12039365470206077, -0.25)]]
        for expected_points, wire in zip(list_expected_points, face_intersections):
            arcellipse = wire.primitives[0]
            self.assertTrue(expected_points[0].is_close(arcellipse.start))
            self.assertTrue(expected_points[1].is_close(arcellipse.interior))
            self.assertTrue(expected_points[2].is_close(arcellipse.end))
        """ ========== TWO PARALLEL LINES ========="""
        plane_face_3 = self.plane_face_cylindricalface_intersec.rotation(volmdlr.O3D, volmdlr.X3D, math.pi)
        face_intersections = plane_face_3.face_intersections(face)
        self.assertEqual(face_intersections[0].primitives[0], edges.LineSegment3D(volmdlr.Point3D(0.15, 0.0, -0.25),
                                                                                  volmdlr.Point3D(0.15, 0.0, 0.25)))
        self.assertEqual(face_intersections[1].primitives[0], edges.LineSegment3D(volmdlr.Point3D(-0.15, 0.0, -0.25),
                                                                                  volmdlr.Point3D(-0.15, 0.0, 0.25)))
        """ ========== ONE LINE ========="""
        plane_face_3 = self.plane_face_cylindricalface_intersec.translation(R * volmdlr.Y3D)
        face_intersections = plane_face_3.face_intersections(face)
        self.assertEqual(face_intersections[0].primitives[0], edges.LineSegment3D(volmdlr.Point3D(0.0, 0.15, -0.25),
                                                                                  volmdlr.Point3D(0.0, 0.15, 0.25)))

    def test_linesegment_inside(self):
        lineseg = volmdlr.edges.LineSegment3D(volmdlr.Point3D(0.2, 0, -0.2), volmdlr.Point3D(0.1, 0.0, 0.2))
        self.assertTrue(self.plane_face_cylindricalface_intersec.linesegment_inside(lineseg))
        lineseg1 = volmdlr.edges.LineSegment3D(volmdlr.Point3D(0.2, 0, -0.2), volmdlr.Point3D(0.1, 0.1, 0.2))
        self.assertFalse(self.plane_face_cylindricalface_intersec.linesegment_inside(lineseg1))

    def test_circle_inside(self):
        circle = volmdlr.wires.Circle3D(volmdlr.OZXY, 0.1)
        self.assertTrue(self.plane_face_cylindricalface_intersec.circle_inside(circle))
        circle2 = volmdlr.wires.Circle3D(volmdlr.OYZX, 0.1)
        self.assertFalse(self.plane_face_cylindricalface_intersec.circle_inside(circle2))

if __name__ == '__main__':
    unittest.main()
