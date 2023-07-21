import unittest
from itertools import product

import volmdlr
from volmdlr.models.edges import bspline1, lineseg, arc, arc_ellipse2d


class TestEdge2DIntersections(unittest.TestCase):

    def test_edge_intersections(self):
        expected_results = [[volmdlr.Point2D(0.21547678763106481, 0.17126976164919136),
                             volmdlr.Point2D(1.4999999999999996, 5.551115123125783e-17),
                             volmdlr.Point2D(2.7845232122344106, -0.17126976163125474)],
                            [volmdlr.Point2D(0.13933375799254094, 0.1205208441552475),
                             volmdlr.Point2D(1.7686124077855423, -0.1299992524853488)],
                            [volmdlr.Point2D(1.3169500217886863, 0.09016194906901831)],
                            [volmdlr.Point2D(0.21547678763106481, 0.17126976164919136),
                             volmdlr.Point2D(1.4999999999999996, 5.551115123125783e-17),
                             volmdlr.Point2D(2.7845232122344106, -0.17126976163125474)],
                            [volmdlr.Point2D(1.8893801268948387, -0.051917350252645156),
                             volmdlr.Point2D(0.07997689764965771, 0.18933641364671233)],
                            [volmdlr.Point2D(1.2821273688601553, 0.029049684151979283)],
                            [volmdlr.Point2D(0.13933375799254094, 0.1205208441552475),
                             volmdlr.Point2D(1.7686124077855423, -0.1299992524853488)],
                            [volmdlr.Point2D(1.8893801268948387, -0.051917350252645156),
                             volmdlr.Point2D(0.07997689764965771, 0.18933641364671233)],
                            [volmdlr.Point2D(1.0591654855246853, -0.30366407228860615)],
                            [volmdlr.Point2D(1.3169500217886863, 0.09016194906901831)],
                            [volmdlr.Point2D(1.2821273688601553, 0.029049684151979283)],
                            [volmdlr.Point2D(1.0591654855246853, -0.30366407228860615)]]

        intersection_results = []
        for edge1, edge2 in product([bspline1, lineseg, arc, arc_ellipse2d], repeat=2):
            if edge1 == edge2:
                continue
            intersections = edge1.intersections(edge2)
            intersection_results.append(intersections)
        for intersections, expected_intersections in zip(intersection_results, expected_results):
            for intersection, expected_intersection in zip(intersections, expected_intersections):
                self.assertTrue(intersection.is_close(expected_intersection))


if __name__ == '__main__':
    unittest.main()
