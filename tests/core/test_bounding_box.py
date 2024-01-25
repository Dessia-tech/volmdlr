import unittest
import volmdlr
from volmdlr.core import BoundingBox
from volmdlr.faces import Triangle3D


class TestBoundingBox(unittest.TestCase):
    def setUp(self):
        self.bbox1 = BoundingBox(0.0, 2.0, 0.0, 2.0, 0.0, 2.0)
        self.bbox2 = BoundingBox(1.0, 3.0, 0.0, 5.0, -1.0, 1.0)
        self.bbox3 = BoundingBox(4.0, 5.0, 4.0, 5.0, 4.0, 5.0)
        self.bbox4 = BoundingBox(0.5, 1.5, 0.5, 1.5, 0.5, 1.5)

    def test_center(self):
        self.assertEqual(self.bbox1.center, volmdlr.Point3D(1.0, 1.0, 1.0))
        self.assertEqual(self.bbox2.center, volmdlr.Point3D(2.0, 2.5, 0.0))

    def test_add(self):
        self.assertEqual(self.bbox1 + self.bbox2, BoundingBox(0.0, 3.0, 0.0, 5.0, -1.0, 2.0))

    def test_to_dict(self):
        self.assertEqual(
            self.bbox1.to_dict(),
            {
                "object_class": "volmdlr.core.BoundingBox",
                "name": "",
                "xmin": 0,
                "xmax": 2,
                "ymin": 0,
                "ymax": 2,
                "zmin": 0,
                "zmax": 2,
            },
        )

    def test_points(self):
        self.assertListEqual(
            self.bbox1.points,
            [
                volmdlr.Point3D(0, 0, 0),
                volmdlr.Point3D(2, 0, 0),
                volmdlr.Point3D(2, 2, 0),
                volmdlr.Point3D(0, 2, 0),
                volmdlr.Point3D(0, 0, 2),
                volmdlr.Point3D(2, 0, 2),
                volmdlr.Point3D(2, 2, 2),
                volmdlr.Point3D(0, 2, 2),
            ],
        )

    def test_from_bounding_boxes(self):
        self.assertEqual(
            BoundingBox(0.0, 5.0, 0.0, 5.0, -1.0, 5.0),
            BoundingBox.from_bounding_boxes([self.bbox1, self.bbox2, self.bbox3]),
        )

    def test_from_points(self):
        self.assertEqual(
            BoundingBox.from_points(
                [
                    volmdlr.Point3D(4.0, 4.0, 4.0),
                    volmdlr.Point3D(4.1, 4.0, 4.0),
                    volmdlr.Point3D(4.0, 4.5, 4.0),
                    volmdlr.Point3D(4.0, 5.0, 4.9),
                    volmdlr.Point3D(5.0, 5.0, 5.0),
                ]
            ),
            self.bbox3,
        )

    def test_to_frame(self):
        frame = self.bbox1.to_frame()
        self.assertEqual(frame.u, volmdlr.Vector3D(2.0, 0.0, 0.0))
        self.assertEqual(frame.v, volmdlr.Vector3D(0.0, 2.0, 0.0))
        self.assertEqual(frame.w, volmdlr.Vector3D(0.0, 0.0, 2.0))
        self.assertEqual(frame.origin, volmdlr.Point3D(1.0, 1.0, 1.0))

    def test_volume(self):
        self.assertEqual(self.bbox1.volume(), 8.0)
        self.assertEqual(self.bbox2.volume(), 20.0)
        self.assertEqual(self.bbox3.volume(), 1.0)

    def test_bbox_intersection(self):
        self.assertTrue(self.bbox1.is_intersecting(self.bbox2))
        self.assertTrue(self.bbox2.is_intersecting(self.bbox1))
        self.assertTrue(self.bbox1.is_intersecting(self.bbox1))
        self.assertFalse(self.bbox1.is_intersecting(self.bbox3))

    def test_is_inside_bbox(self):
        self.assertFalse(self.bbox2.is_inside_bbox(self.bbox1))
        self.assertFalse(self.bbox3.is_inside_bbox(self.bbox1))
        self.assertTrue(self.bbox4.is_inside_bbox(self.bbox1))

    def test_intersection_volume(self):
        self.assertEqual(self.bbox1.intersection_volume(self.bbox2), 2.0)
        self.assertEqual(self.bbox1.intersection_volume(self.bbox3), 0.0)
        self.assertEqual(self.bbox1.intersection_volume(self.bbox4), 1.0)

    def test_distance_to_bbox(self):
        self.assertEqual(self.bbox1.distance_to_bbox(self.bbox2), 0.0)
        self.assertEqual(self.bbox1.distance_to_bbox(self.bbox3), 12**0.5)
        self.assertEqual(self.bbox1.distance_to_bbox(self.bbox4), 0.0)
        bbox5 = BoundingBox(-2.0, -1.0, 3.0, 4.0, -2.0, -1.0)
        self.assertEqual(self.bbox1.distance_to_bbox(bbox5), 3**0.5)
        self.assertEqual(bbox5.distance_to_bbox(self.bbox1), 3**0.5)

    def test_point_belongs(self):
        self.assertTrue(self.bbox1.point_inside(volmdlr.Point3D(1.0, 1.0, 1.0)))
        self.assertTrue(self.bbox1.point_inside(volmdlr.Point3D(1.0, 1.0, 0.0)))
        self.assertFalse(self.bbox1.point_inside(volmdlr.Point3D(3.0, 3.0, 3.0)))

    def test_distance_to_point(self):
        p0 = volmdlr.O3D
        self.assertEqual(self.bbox1.distance_to_point(p0), 0.0)
        p1 = volmdlr.Point3D(-1.0, -1.0, -1.0)
        self.assertEqual(self.bbox1.distance_to_point(p1), 3**0.5)
        p2 = volmdlr.Point3D(3.0, 3.0, 3.0)
        self.assertEqual(self.bbox1.distance_to_point(p2), 3**0.5)
        p3 = volmdlr.Point3D(1.0, 1.0, 1.0)
        self.assertEqual(self.bbox1.distance_to_point(p3), 1.0)

    def test_plot(self):
        ax = self.bbox1.plot()
        self.assertEqual(len(ax.lines), 12)

    def test_is_close(self):
        self.assertTrue(self.bbox1.is_close(self.bbox1))

        bbox_almost_equal_1 = BoundingBox(0, 1, 0, 1, 0, 1)
        bbox_almost_equal_2 = BoundingBox(0 + 1e-9, 1 + 1e-9, 0 + 1e-9, 1 + 1e-9, 0 + 1e-9, 1 + 1e-9)

        self.assertTrue(bbox_almost_equal_1.is_close(bbox_almost_equal_2))
        self.assertFalse(bbox_almost_equal_1.is_close(bbox_almost_equal_2, tol=1e-10))
        self.assertNotEqual(bbox_almost_equal_1, bbox_almost_equal_2)

    def test_scale(self):
        scaled_bbox = self.bbox1.scale(2)
        self.assertEqual(BoundingBox(-1.0, 3.0, -1.0, 3.0, -1.0, 3.0), scaled_bbox)

    def test_is_intersecting_triangle(self):
        # Triangle intersecting the bounding box
        triangle1 = Triangle3D(volmdlr.Point3D(5, -3, 0.5), volmdlr.Point3D(-3, 5, 0.5), volmdlr.Point3D(5, 5, 0.5))
        self.assertTrue(self.bbox1.is_intersecting_triangle(triangle1))

        # Triangle completely inside the bounding box
        triangle2 = Triangle3D(
            volmdlr.Point3D(0.5, 0.5, 0.5), volmdlr.Point3D(1.5, 0.5, 0.5), volmdlr.Point3D(1, 1.5, 1)
        )
        self.assertTrue(self.bbox1.is_intersecting_triangle(triangle2))

        # Triangle with one vertex inside the bounding box
        triangle3 = Triangle3D(volmdlr.Point3D(1, 1, 1), volmdlr.Point3D(3, 0, 0), volmdlr.Point3D(3, 3, 0))
        self.assertTrue(self.bbox1.is_intersecting_triangle(triangle3))

        # Triangle with an edge passing through the bounding box
        triangle4 = Triangle3D(volmdlr.Point3D(1, -1, 1), volmdlr.Point3D(1, 3, 1), volmdlr.Point3D(3, 0, 3))
        self.assertTrue(self.bbox1.is_intersecting_triangle(triangle4))

        # Triangle that just touches one face of the bounding box
        triangle5 = Triangle3D(volmdlr.Point3D(0, 2, 1), volmdlr.Point3D(0, 3, 2), volmdlr.Point3D(0, 3, 0))
        self.assertTrue(self.bbox1.is_intersecting_triangle(triangle5))

        # Triangle that is completely outside and does not intersect
        triangle6 = Triangle3D(volmdlr.Point3D(-1, -1, -1), volmdlr.Point3D(-2, -1, -1), volmdlr.Point3D(-1, -2, -1))
        self.assertFalse(self.bbox1.is_intersecting_triangle(triangle6))

        # Edge case: Triangle is coplanar with one of the bounding box faces but outside
        triangle7 = Triangle3D(volmdlr.Point3D(3, 3, 2), volmdlr.Point3D(4, 3, 2), volmdlr.Point3D(3, 4, 2))
        self.assertFalse(self.bbox1.is_intersecting_triangle(triangle7))

        # Edge case: Triangle intersects at exactly one vertex of the bounding box
        triangle8 = Triangle3D(volmdlr.Point3D(2, 2, 2), volmdlr.Point3D(3, 3, 3), volmdlr.Point3D(2, 3, 3))
        self.assertTrue(self.bbox1.is_intersecting_triangle(triangle8))


if __name__ == "__main__":
    unittest.main()
