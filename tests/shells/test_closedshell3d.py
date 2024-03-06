import math
import unittest
import os
import numpy

import dessia_common.core
import volmdlr
from volmdlr import edges, faces, primitives3d, wires, surfaces, shells, step


folder = os.path.join(os.path.dirname(os.path.realpath(__file__)))
objects_folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "closedshell_objects")


class TestClosedShell3D(unittest.TestCase):

    def test_union(self):
        stepfile1 = step.Step.from_file(os.path.join(objects_folder, 'test_faces.step'))
        volume1 = stepfile1.to_volume_model()
        closed_shell = volume1.primitives[0]
        closed_shell1 = closed_shell.rotation(volmdlr.O3D, volmdlr.Y3D, math.pi / 7)
        union = closed_shell.union(closed_shell1)[0]
        expected_shell_faces_areas = dessia_common.core.DessiaObject.from_json(os.path.join(
            objects_folder, 'test_union_expected_faces_areas.json')).primitives
        self.assertEqual(len(union.faces), len(expected_shell_faces_areas))
        for i, area in enumerate(sorted([face.area() for face in union.faces])):
            self.assertAlmostEqual(area, expected_shell_faces_areas[i], 4)

    def test_intersection(self):
        extrusion1, extrusion2 = dessia_common.core.DessiaObject.from_json(
            os.path.join(objects_folder, 'test_closedshell_intersection.json'))
        intersection = extrusion1.intersection(extrusion2)
        expected_areas = [0.0001603762338471547, 0.0001603762338598547, 0.0002258060000078384, 0.0005865250644030277,
                          0.0005865250644053065, 0.0005865256929828994, 0.0005865256929851797, 0.0006507280562579344,
                          0.0006507281997243776, 0.0006507287536446724, 0.0006507288971112704, 0.002079815726730889,
                          0.0020798157267342588, 0.0027911081188436927, 0.002791109662733529, 0.0037635988992237926,
                          0.0037635988992285145, 0.0038394588771942153, 0.0038394588771985295, 0.0046984506173449216,
                          0.004698453152620571, 0.005812094842241624, 0.005812103536328831, 0.0058121998046592235,
                          0.005812200540825445, 0.00581220473314671, 0.005812204737898536, 0.00581220602916541,
                          0.005812206548041761, 0.005812209253899912, 0.005812210963215569, 0.005812210964111839,
                          0.005812210965759321, 0.010283114506497748, 0.010283123293653551, 0.0102831564624514,
                          0.010283183708478005, 0.012309686730356234, 0.012309686730356262, 0.021026086509749957,
                          0.02775624156108652, 0.027756241561086553, 0.027756279461254495, 0.0277562794612545,
                          0.04320287693191287, 0.04320287693191288, 0.04710646280458197, 0.047106508562932636,
                          0.07318715514507368, 0.0797964534011808, 0.07979645340118091, 0.14314330722704466]
        self.assertEqual(len(intersection.faces), len(expected_areas))
        for i, area in enumerate(sorted([face.area() for face in intersection.faces])):
            self.assertAlmostEqual(area, expected_areas[i])

    def test_consecutive_boolean_operations(self):
        block = primitives3d.Block(volmdlr.OXYZ, color=(1, 0.3, 0.3), alpha=0.6)
        sphere = primitives3d.Sphere(volmdlr.O3D, 0.68, color=(0.2, 1, 0.3), alpha=.6)

        block_inters_sphere = block.intersection(sphere)[0]
        expected_areas1 = [0.6672300326944091, 0.6672300326944091, 0.6672300326944091, 0.6672300326944091,
                           0.6672300326944091, 0.6672742796224723, 3.2330061945128783]
        self.assertEqual(len(block_inters_sphere.faces), len(expected_areas1))
        for i, face_area in enumerate(sorted([f.area() for f in block_inters_sphere.faces])):
            self.assertAlmostEqual(face_area, expected_areas1[i])

        cyl1 = primitives3d.Cylinder(volmdlr.OXYZ, 0.3, 2, color=(1, 0, 0))
        cyl2 = primitives3d.Cylinder(volmdlr.OYZX, 0.3, 2, color=(1, 1, 0))
        cyl3 = primitives3d.Cylinder(volmdlr.OZXY, 0.3, 2, color=(0, 1, 1))

        cyl1_union_cyl2 = cyl1.union(cyl2)[0]
        union_cyl3 = cyl1_union_cyl2.union(cyl3)[0]

        expected_areas2 = [0.28272459012474943, 0.28272459012474943, 0.28272459012474943, 0.28272459012474943,
                           0.28272459012474943, 0.28272459012474943, 4.586132621569709, 4.58613262156971,
                           4.586135300242545, 4.586137982316728, 4.5861380209112195, 4.586140703735797]
        self.assertEqual(len(union_cyl3.faces), len(expected_areas2))
        for i, face_area in enumerate(sorted([f.area() for f in union_cyl3.faces])):
            self.assertAlmostEqual(face_area, expected_areas2[i])

        substraction = block_inters_sphere.subtract_to_closed_shell(union_cyl3)[0]
        expected_areas3 = [0.38450544256965963, 0.38450544256965963, 0.38450544256965963, 0.38450544256965963,
                           0.38450544256965963, 0.3845496894977229, 1.4445399679799156, 1.4445399679799158,
                           1.4445426466527518, 1.444545328726934, 1.4445453673214272, 1.4445480501460044,
                           3.2330061945128783]

        self.assertEqual(len(substraction.faces), len(expected_areas3))
        for i, face_area in enumerate(sorted([f.area() for f in substraction.faces])):
            self.assertAlmostEqual(face_area, expected_areas3[i])

    def test_union_adjacent_blocks(self):
        frame1 = volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0), volmdlr.X3D, volmdlr.Y3D, volmdlr.Z3D)
        frame2 = volmdlr.Frame3D(volmdlr.Point3D(0, 1, 0), volmdlr.X3D, volmdlr.Y3D, 3 * volmdlr.Z3D)
        block1 = primitives3d.Block(frame1)
        block2 = primitives3d.Block(frame2)
        closed_shell1 = shells.ClosedShell3D(block1.faces)
        closed_shell2 = shells.ClosedShell3D(block2.faces)
        block3 = closed_shell1.union(closed_shell2)
        block3[0].merge_faces()
        self.assertEqual(len(block3[0].faces), 10)

        contour = volmdlr.wires.ClosedPolygon2D([volmdlr.Point2D(0, 0), volmdlr.Point2D(-1, 0),
                                                 volmdlr.Point2D(-1, 1), volmdlr.Point2D(1, 1),
                                                 volmdlr.Point2D(1, -1), volmdlr.Point2D(0, -1)])
        extrude1 = volmdlr.primitives3d.ExtrudedProfile(
            volmdlr.OXYZ, contour, [], -1)
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
        extruded_prifile1 = primitives3d.ExtrudedProfile(volmdlr.OYZX,
                                                         wires.Contour2D(contour_primitives), [],
                                                         0.01, (0.4, 0.1, 0.1), 0.6)
        extruded_prifile2 = extruded_prifile1.translation(volmdlr.Vector3D(0.01, 0, 0))
        union_shell1_shell2 = extruded_prifile1.union(extruded_prifile2)[0]
        union_shell1_shell2.merge_faces()
        self.assertEqual(len(union_shell1_shell2.faces), 7)
        boundary1 = primitives3d.Block(volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D, 0.3 * volmdlr.Y3D, 0.1 * volmdlr.Z3D))
        boundary2 = primitives3d.Block(
            volmdlr.Frame3D(volmdlr.O3D, 0.3 * volmdlr.X3D, 0.8 * volmdlr.Y3D, 0.2 * volmdlr.Z3D))
        boundary2 = boundary2.translation(offset=(0.5 + 0.15) * volmdlr.X3D)
        union = boundary1.union(boundary2)[0]
        self.assertEqual(len(union.faces), 11)

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
                                                           volmdlr.Vector3D(0, 1.8, 0), volmdlr.Vector3D(0, 0, 1)),
                                           'old')
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
        extrude1 = volmdlr.primitives3d.ExtrudedProfile(volmdlr.OXYZ,
                                                        contour, [], -1)
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

        plane3d_1 = surfaces.Plane3D.from_plane_vectors(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D)
        surf2d_1 = surfaces.Surface2D(poly1_vol1.to_2d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D), [])

        plane3d_2 = surfaces.Plane3D.from_plane_vectors(0.3 * volmdlr.Z3D.to_point(), volmdlr.X3D, volmdlr.Y3D)
        surf2d_2 = surfaces.Surface2D(poly3_vol1.to_2d(volmdlr.O3D, volmdlr.X3D, volmdlr.Y3D), [])
        shell_faces += [faces.PlaneFace3D(plane3d_1, surf2d_1), faces.PlaneFace3D(plane3d_2, surf2d_2)]

        shell1 = shells.ClosedShell3D(shell_faces)
        shell2 = shell1.translation(volmdlr.Point3D(0, -0.28, -0.2)).rotation(volmdlr.O3D, volmdlr.X3D, math.pi)
        union_shell1_shell2 = shell1.union(shell2)
        self.assertEqual(len(union_shell1_shell2), 2)

    def test_cut_by_plane(self):
        boundary1 = primitives3d.Block(volmdlr.Frame3D(volmdlr.O3D, volmdlr.X3D, 0.3 * volmdlr.Y3D, 0.1 * volmdlr.Z3D))
        boundary2 = primitives3d.Block(
            volmdlr.Frame3D(volmdlr.O3D, 0.4 * volmdlr.X3D, 0.8 * volmdlr.Y3D, 0.4 * volmdlr.Z3D))
        boundary2 = boundary2.translation(offset=(0.5 + 0.14) * volmdlr.X3D)
        boundary2 = boundary2.translation(offset=(0.1) * volmdlr.Z3D)
        union = boundary1.union(boundary2)[0]
        center = union.bounding_box.center
        plane = surfaces.Plane3D.from_normal(center, volmdlr.Y3D)
        cut_by_plane = union.cut_by_plane(plane)
        self.assertEqual(len(cut_by_plane), 1)
        self.assertAlmostEqual(cut_by_plane[0].area(), 0.254)

    def test_intersection(self):
        box1 = primitives3d.Block(
            volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0), volmdlr.Vector3D(0.6, 0, 0),
                            volmdlr.Vector3D(0, 0.6, 0), volmdlr.Vector3D(0, 0, 0.6)), color=(1, 0.2, 0.2),
            alpha=0.6)
        box2 = primitives3d.Block(
            volmdlr.Frame3D(volmdlr.Point3D(0, 0, 0), volmdlr.Vector3D(0.3, 0, 0),
                            volmdlr.Vector3D(0, 0.3, 0), volmdlr.Vector3D(0, 0, 0.3)), color=(.1, 0.2, 1),
            alpha=0.6)
        self.assertEqual(box1.intersection(box2)[0], box2)
        box3 = primitives3d.Block(
            volmdlr.Frame3D(volmdlr.Point3D(.3, 0, 0), volmdlr.Vector3D(0.3, 0, 0),
                            volmdlr.Vector3D(0, 0.3, 0), volmdlr.Vector3D(0, 0, 0.3)), color=(.1, 0.2, 1), alpha=0.6)
        expected_box = primitives3d.Block(
            volmdlr.Frame3D(volmdlr.Point3D(.15 * 3/2, 0, 0), volmdlr.Vector3D(0.15, 0, 0),
                            volmdlr.Vector3D(0, 0.3, 0), volmdlr.Vector3D(0, 0, 0.3)), color=(.1, 0.2, 1), alpha=0.6)
        self.assertTrue(all(expected_box.face_on_shell(face) for face in box1.intersection(box3)[0].faces))
        box4 = primitives3d.Block(
            volmdlr.Frame3D(volmdlr.Point3D(.5, 0, 0), volmdlr.Vector3D(0.3, 0, 0),
                            volmdlr.Vector3D(0, 0.3, 0), volmdlr.Vector3D(0, 0, 0.3)), color=(.1, 0.2, 1), alpha=0.6)
        self.assertFalse(box1.intersection(box4))

    def test_point_belongs(self):
        closed_shell = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, 'test_closed_shell_point_belongs2.json')).primitives[0]
        points = [volmdlr.Point3D(-.2, -0.6, 0.08), volmdlr.Point3D(-0.340920128805, -0.418071198223, 0.007036661148),
                  volmdlr.Point3D(-0.287522562519, -0.574786328164, 0.157256628036),
                  volmdlr.Point3D(-0.314221345662, -0.522547951517, 0.057109983444)]
        expected_results = [True, True, False, True]
        for i, expected_result in enumerate(expected_results):
            self.assertEqual(closed_shell.point_inside(points[i]), expected_result)

    def test_minimum_distance(self):
        closed_shell = dessia_common.core.DessiaObject.from_json(
            os.path.join(folder, 'test_shells_distance2.json'))
        u_vector = volmdlr.Vector3D(-0.5773502691896258, -0.5773502691896258, -0.5773502691896258)
        v_vector = volmdlr.Vector3D(0.8164965809277258, -0.40824829046386313, -0.40824829046386313)
        w_vector = volmdlr.Vector3D(0.0, -0.7071067811865476, 0.7071067811865476)
        frame = volmdlr.Frame3D(volmdlr.Point3D(-0.01661984584195119, -0.04221251977732219, -0.04351622102493058),
                                u_vector, v_vector, w_vector)
        fm_shell = closed_shell.frame_mapping(frame, 'new')
        min_distance = closed_shell.minimum_distance(fm_shell, False)
        self.assertAlmostEqual(min_distance, 0.022811959708641426, 6)
        frame = volmdlr.Frame3D(volmdlr.Point3D(0.011516851705803667, 0.012859651289434018, 0.015147046170848444),
                                u_vector, v_vector, w_vector)
        fm_shell = closed_shell.frame_mapping(frame, 'new')
        min_distance, point1, point2 = closed_shell.minimum_distance(fm_shell, True)
        self.assertEqual(min_distance, 0.0)

    def test_volume(self):
        closed_shell = dessia_common.core.DessiaObject.from_json(os.path.join(folder, 'test_shell_volume.json'))
        closed_shell2 = closed_shell.rotation(volmdlr.O3D - 0.95 * volmdlr.X3D, volmdlr.Z3D, numpy.pi / 2)
        closed_shell2 = closed_shell2.translation(-0.95 * volmdlr.X3D + 0.45 * volmdlr.Z3D + 0.2 * volmdlr.Y3D)
        self.assertAlmostEqual(closed_shell.volume(), closed_shell2.volume())


if __name__ == '__main__':
    unittest.main()
