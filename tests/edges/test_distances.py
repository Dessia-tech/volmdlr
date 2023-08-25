import unittest
from itertools import product

import volmdlr
import volmdlr.edges as vme
import volmdlr.models.edges as edges_models


class TestEdgesDistances(unittest.TestCase):
    def test_edges_distances_2d(self):
        expected_distances = [3.0, 3.171929282184431, 1.413379747337626, 2.9154759474226504, 3.144603017460846,
                              3.2934193036985913, 1.4282856857085704, 2.8735269769962803, 3.183195634502654,
                              3.3500730304603152, 1.4456832294800965, 2.8575751838382266, 1.0000000009680357,
                              1.2027923774744735, 0.0, 1.2597969896641665]

        vector1 = volmdlr.Vector2D(1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        distances = []
        for edge1, edge2 in product([edges_models.bspline1, edges_models.lineseg,
                                     edges_models.arc, edges_models.arc_ellipse2d], repeat=2):
            edge2 = edge2.frame_mapping(volmdlr.Frame2D(volmdlr.Point2D(3, 3), vector1, vector2), 'new')
            dist, min_dist_point_arc, min_dist_point_lineseg = edge1.minimum_distance(edge2, True)
            distances.append(dist)
        for distance, expected_distance in zip(distances, expected_distances):
            self.assertAlmostEqual(distance, expected_distance)

    def test_edges_distances_3d(self):
        expected_distances = [3.297515976380515, 0.6600032792406747, 0.01931857722664732, 0.6383061365421937,
                              0.6600032792406747, 6.117055924652082, 2.6674164255110373, 1.2247448713912275,
                              0.019318579007616275, 2.6674164255110373, 5.084297068558137, 0.8731527317517747,
                              0.6383061445002628, 1.2247448713912275, 0.8731527317517747, 4.289502509060332]
        vector1 = volmdlr.Vector3D(1, 1, 1)
        vector1 = vector1.unit_vector()
        vector2 = vector1.deterministic_unit_normal_vector()
        vector3 = vector1.cross(vector2)

        distances = []
        for edge1, edge2 in product([edges_models.bspline_curve3d(), edges_models.linesegment3d(),
                                     edges_models.arc3d(), edges_models.arc_ellipse3d()], repeat=2):
            if edge1 == edge2:
                edge2 = edge2.frame_mapping(volmdlr.Frame3D(volmdlr.Point3D(2, 5, -3), vector1, vector2, vector3),
                                            'new')
            dist, min_dist_point_arc, min_dist_point_lineseg = edge1.minimum_distance(edge2, True)
            distances.append(dist)
        for distance, expected_distance in zip(distances, expected_distances):
            self.assertAlmostEqual(distance, expected_distance)


if __name__ == '__main__':
    unittest.main()
