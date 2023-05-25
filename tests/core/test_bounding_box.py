import unittest
import volmdlr
from volmdlr.core import BoundingBox


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
            BoundingBox.from_bounding_boxes([self.bbox1, self.bbox2, self.bbox3]),
            BoundingBox(0.0, 5.0, 0.0, 5.0, -1.0, 5.0),
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
        self.assertTrue(self.bbox1.bbox_intersection(self.bbox2))
        self.assertTrue(self.bbox2.bbox_intersection(self.bbox1))
        self.assertTrue(self.bbox1.bbox_intersection(self.bbox1))
        self.assertFalse(self.bbox1.bbox_intersection(self.bbox3))

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
        self.assertTrue(self.bbox1.point_belongs(volmdlr.Point3D(1.0, 1.0, 1.0)))
        self.assertTrue(self.bbox1.point_belongs(volmdlr.Point3D(1.0, 1.0, 0.0)))
        self.assertFalse(self.bbox1.point_belongs(volmdlr.Point3D(3.0, 3.0, 3.0)))

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


if __name__ == "__main__":
    unittest.main()
