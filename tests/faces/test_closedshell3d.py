import math
import unittest
from volmdlr import primitives3d, wires, faces, edges
import volmdlr


class TestPlaneFace3D(unittest.TestCase):
    def test_face_intersections(self):
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


class TestClosedShell3D(unittest.TestCase):

    def test_union_adjacent_blocks(self):
        frame1 = volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        frame2 = volmdlr.Frame3D(volmdlr.Point3D(0, 1, 0), volmdlr.X3D, volmdlr.Y3D, 3 * volmdlr.Z3D)
        block1 = primitives3d.Block(frame1)
        block2 = primitives3d.Block(frame2)
        closed_shell1 = faces.ClosedShell3D(block1.faces)
        closed_shell2 = faces.ClosedShell3D(block2.faces)
        block3 = closed_shell1.union(closed_shell2)
        block3[0].merge_faces()
        self.assertEqual(len(block3[0].faces), 10)

        contour = volmdlr.wires.ClosedPolygon2D([volmdlr.Point2D(0, 0), volmdlr.Point2D(-1, 0),
                                                 volmdlr.Point2D(-1, 1), volmdlr.Point2D(1, 1),
                                                 volmdlr.Point2D(1, -1), volmdlr.Point2D(0, -1)])
        extrude1 = volmdlr.primitives3d.ExtrudedProfile(
            volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D, contour, [], -volmdlr.Z3D)
        frame1 = volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0.5), 2 * volmdlr.X3D, 2 * volmdlr.Y3D, volmdlr.Z3D)
        block1 = volmdlr.primitives3d.Block(frame1)
        union1 = extrude1.union(block1)
        union1[0].merge_faces()
        self.assertEqual(len(union1[0].faces), 9)

        contour_primitives = [edges.LineSegment2D(start, end) for start, end in [
            (volmdlr.Point2D(0.0, 0.0), volmdlr.Point2D(-0.03, 0.0)),
            (volmdlr.Point2D(-0.03, 0.0), volmdlr.Point2D(-0.03, 0.02)),
            (volmdlr.Point2D(-0.03, 0.02), volmdlr.Point2D(-0.020436, 0.029871)),
            (volmdlr.Point2D(-0.020436, 0.029871), volmdlr.Point2D(0.0, 0.029871)),
            (volmdlr.Point2D(0.0, 0.029871), volmdlr.Point2D(0.0, 0.0))]]
        extruded_prifile1 = primitives3d.ExtrudedProfile(volmdlr.O3D, volmdlr.Y3D, volmdlr.Z3D,
                                                         wires.Contour2D(contour_primitives), [],
                                                         volmdlr.Vector3D(0.01, 0, 0), (0.4, 0.1, 0.1), 0.6)
        extruded_prifile2 = extruded_prifile1.translation(volmdlr.Vector3D(0.01, 0, 0))
        union_shell1_shell2 = extruded_prifile1.union(extruded_prifile2)[0]
        union_shell1_shell2.merge_faces()
        self.assertEqual(len(union_shell1_shell2.faces), 7)

    def test_set_operations_blocks(self):

        box_red = primitives3d.Block(
            volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0), volmdlr.Vector3D(0.4, 0, 0),
                       volmdlr.Vector3D(0, 0.4, 0), volmdlr.Vector3D(0, 0, 0.4)),
            color=(0.2, 1, 0.4), alpha=0.6)
        box_green = box_red.frame_mapping(volmdlr.Frame3D(volmdlr.Point3D(-0.4, 0, -0.1), volmdlr.Vector3D(1, 0, 0),
                                                     volmdlr.Vector3D(0, 1, 0), volmdlr.Vector3D(0, 0, 1)), 'new')
        union_red_green_boxes = box_red.union(box_green)[0]
        union_red_green_boxes.merge_faces()
        self.assertEqual(len(union_red_green_boxes.faces), 10)

        box_blue = primitives3d.Block(
            volmdlr.Frame3D(volmdlr.Point3D(0.1, 0, 0), volmdlr.Vector3D(0.2, 0, 0),
                       volmdlr.Vector3D(0, 0.1, 0), volmdlr.Vector3D(0, 0, 1)),
            alpha=0.6)
        box_blue2 = box_blue.frame_mapping(volmdlr.Frame3D(volmdlr.Point3D(0.2, 0, 0), volmdlr.Vector3D(1, 0, 0),
                                                      volmdlr.Vector3D(0, 1.8, 0), volmdlr.Vector3D(0, 0, 1)), 'old')
        union_blue_blue2_boxes = box_blue.union(box_blue2)[0]
        union_blue_blue2_boxes.merge_faces()
        self.assertEqual(len(union_blue_blue2_boxes.faces), 10)

        union_box = union_red_green_boxes.union(union_blue_blue2_boxes)[0]
        union_box.merge_faces()
        self.assertEqual(len(union_box.faces), 28)
        # subtraction_box = union_red_green_boxes.subtract(union_blue_blue2_boxes)
        intersection_box = union_red_green_boxes.intersection(union_blue_blue2_boxes)[0]
        intersection_box.merge_faces()
        self.assertEqual(len(intersection_box.faces), 12)
        # subtraction_closedbox = union_red_green_boxes.subtract_to_closed_shell(union_blue_blue2_boxes)[0]
        # subtraction_closedbox.merge_faces()
        # self.assertEqual(len(subtraction_closedbox), 16)

    def test_union_equal_overlapping_blocks(self):
        frame1 = volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        frame2 = volmdlr.Frame3D(volmdlr.Point3D(0, 0.8, 0), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        block1 = volmdlr.primitives3d.Block(frame1)
        block2 = volmdlr.primitives3d.Block(frame2)
        block3 = block1.union(block2)
        block3[0].merge_faces()
        self.assertEqual(len(block3[0].faces), 6)

    def test_union_block1_inside_block1(self):
        contour = wires.ClosedPolygon2D([volmdlr.Point2D(0, 0), volmdlr.Point2D(-1, 0),
                                         volmdlr.Point2D(-1, 1), volmdlr.Point2D(1, 1),
                                         volmdlr.Point2D(1, -1), volmdlr.Point2D(0, -1)])
        extrude1 = volmdlr.primitives3d.ExtrudedProfile(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D,
                                                        contour, [], -volmdlr.Z3D)
        frame1 = volmdlr.Frame3D(volmdlr.Point3D(0, 0.5, -0.5), 2 * volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        block1 = volmdlr.primitives3d.Block(frame1)
        union1 = extrude1.union(block1)[0]
        self.assertEqual(len(union1.faces), 8)

    def test_union_two_disjoint_objects(self):
        poly1_vol1 = wires.ClosedPolygon3D([volmdlr.Point3D(-0.1, -0.05, 0), volmdlr.Point3D(-0.15, 0.1, 0),
                                            volmdlr.Point3D(0.05, 0.2, 0), volmdlr.Point3D(0.12, 0.15, 0),
                                            volmdlr.Point3D(0.1, -0.02, 0)])

        poly2_vol1 = poly1_vol1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi).translation(0.2 * volmdlr.Z3D)
        poly3_vol1 = poly2_vol1.rotation(volmdlr.O3D, volmdlr.Z3D, math.pi / 8).translation(
            0.1 * (volmdlr.Z3D + volmdlr.X3D + volmdlr.Y3D))

        shell_faces = [faces.Triangle3D(*points)
                       for points in poly1_vol1.sewing(poly2_vol1, volmdlr.X3D, volmdlr.Y3D)] + \
                      [faces.Triangle3D(*points)
                       for points in poly2_vol1.sewing(poly3_vol1, volmdlr.X3D, volmdlr.Y3D)]

        plane3d_1 = faces.Plane3D.from_plane_vectors(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D)
        surf2d_1 = faces.Surface2D(poly1_vol1.to_2d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D), [])

        plane3d_2 = faces.Plane3D.from_plane_vectors(0.3 * volmdlr.Z3D, volmdlr.X3D, volmdlr.Y3D)
        surf2d_2 = faces.Surface2D(poly3_vol1.to_2d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D), [])
        shell_faces += [faces.PlaneFace3D(plane3d_1, surf2d_1), faces.PlaneFace3D(plane3d_2, surf2d_2)]

        shell1 = faces.ClosedShell3D(shell_faces)
        shell2 = shell1.translation(volmdlr.Point3D(0, -0.28, -0.2)).rotation(volmdlr.O3D, volmdlr.X3D, math.pi)
        union_shell1_shell2 = shell1.union(shell2)
        self.assertEqual(len(union_shell1_shell2), 2)


if __name__ == '__main__':
    unittest.main()
