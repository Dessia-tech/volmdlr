import unittest
import volmdlr
from volmdlr import primitives3d, wires, faces, edges


class TestPlaneFace3D(unittest.TestCase):
    def test_face_intersections_with_holes(self):
        contour_primitives = [edges.LineSegment2D(start, end) for start, end in
                              [(volmdlr.Point2D(-0.2, -0.2), volmdlr.Point2D(0.2, -0.2)),
                               (volmdlr.Point2D(0.2, -0.2), volmdlr.Point2D(0.2, 0.2)),
                               (volmdlr.Point2D(0.2, 0.2), volmdlr.Point2D(-0.2, 0.2)),
                               (volmdlr.Point2D(-0.2, 0.2), volmdlr.Point2D(-0.2, -0.2))]]
        contour = wires.Contour2D(contour_primitives)

        inter_contour1 = contour.frame_mapping(volmdlr.Frame2D(volmdlr.Point2D(-0.1, -0.05), volmdlr.Vector2D(0.3, 0),
                                                               volmdlr.Vector2D(0, 0.4)), 'old')
        inter_contour2 = contour.frame_mapping(volmdlr.Frame2D(volmdlr.Point2D(0.1, 0.05), volmdlr.Vector2D(0.3, 0),
                                                               volmdlr.Vector2D(0, 0.4)), 'old')
        surface3d_1 = volmdlr.faces.Plane3D(frame=volmdlr.Frame3D(origin=volmdlr.Point3D(0, 0, 0),
                                                                  u=volmdlr.Point3D(1, 0, 0),
                                                                  v=volmdlr.Point3D(0, 1, 0),
                                                                  w=volmdlr.Point3D(0, 0, 1)))
        face_to_cut = volmdlr.faces.PlaneFace3D(surface3d_1, volmdlr.faces.Surface2D(
            contour, [inter_contour1, inter_contour2]))

        surface3d_2 = volmdlr.faces.Plane3D(frame=volmdlr.Frame3D(origin=volmdlr.Point3D(0, 0, 0),
                                                                  u=-volmdlr.Point3D(0, 0, 1),
                                                                  v=volmdlr.Point3D(1, 0, 0),
                                                                  w=-volmdlr.Point3D(0, 1, 0)))
        contour_face_to_cut = contour.frame_mapping(volmdlr.Frame2D(
            volmdlr.O2D, volmdlr.Vector2D(1.5, 0), volmdlr.Vector2D(0, 1.5)), 'old')
        surface2d = volmdlr.faces.Surface2D(contour_face_to_cut, [])
        face = volmdlr.faces.PlaneFace3D(surface3d_2, surface2d)
        face = face.rotation(volmdlr.O3D, volmdlr.Z3D, 3.14 / 3)
        face_intersections = face.face_intersections(face_to_cut)
        self.assertEqual(len(face_intersections), 3)


if __name__ == '__main__':
    unittest.main()
