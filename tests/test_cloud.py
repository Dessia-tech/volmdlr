import unittest

import numpy as npy

import volmdlr as vm
from volmdlr import cloud, faces, shells


class TestCloud(unittest.TestCase):
    points_cloud = cloud.PointCloud3D([vm.Point3D(0, 0, 0), vm.Point3D(10, 10, 10),vm.Point3D(-2.5, 1.2, 0)])

    point1 = vm.Point3D(1, 1, 1)
    point2 = vm.Point3D(1, 2, 1)
    point3 = vm.Point3D(1, 1, 2)
    point4 = vm.Point3D(2, 1, 1)
    point5 = vm.Point3D(2, 2, 1)
    point6 = vm.Point3D(1, 2, 2)
    point7 = vm.Point3D(2, 1, 2)
    point8 = vm.Point3D(2, 2, 2)
    open_triangle_shell = shells.OpenTriangleShell3D([faces.Triangle3D(point1, point2, point3),
                                                     faces.Triangle3D(point3, point4, point5),
                                                     faces.Triangle3D(point5, point6, point7),
                                                     faces.Triangle3D(point7, point8, point6)])

    def test_shell_distances(self):
        images, distances, nearest_shells = self.points_cloud.shell_distances(self.open_triangle_shell)

        self.assertIsInstance(images, cloud.PointCloud3D)
        self.assertIsInstance(distances, list)
        self.assertIsInstance(nearest_shells, list)

        self.assertEqual(nearest_shells, [0, 3, 0])
        self.assertEqual(round(distances[0], 2), 1.73)
        self.assertEqual(round(distances[1], 2), 13.86)
        self.assertEqual(round(distances[2], 2), 3.64)

    def test_shell_distances_ndarray(self):
        images, distances, nearest_shells = self.points_cloud.shell_distances_ndarray(self.open_triangle_shell)

        self.assertIsInstance(images, npy.ndarray)
        self.assertIsInstance(distances, npy.ndarray)
        self.assertIsInstance(nearest_shells, npy.ndarray)

        self.assertTrue(npy.all(nearest_shells == npy.array([0, 3, 0])))
        self.assertEqual(round(distances[0], 2), 1.73)
        self.assertEqual(round(distances[1], 2), 13.86)
        self.assertEqual(round(distances[2], 2), 3.64)

if __name__ == '__main__':
    unittest.main()
