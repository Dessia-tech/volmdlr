import unittest
import volmdlr
from volmdlr import edges, wires


class TestEdgeCollection3D(unittest.TestCase):
    def test_to_dict(self):
        pt_1 = volmdlr.Point3D(1, 2, 3)
        pt_2 = volmdlr.Point3D(2, 3, 4)
        edge = edges.LineSegment3D(pt_1, pt_2)
        collection = wires.EdgeCollection3D([edge])
        objt_dict = collection.to_dict()

        self.assertEqual(len(objt_dict["primitives"]), 1)
        self.assertEqual(objt_dict["primitives"][0], {'object_class': 'volmdlr.edges.LineSegment3D',
                                                      'name': '',
                                                      'start': {'object_class': 'volmdlr.Point3D',
                                                                'x': 1.0, 'y': 2.0, 'z': 3.0},
                                                      'end': {'object_class': 'volmdlr.Point3D',
                                                              'x': 2.0, 'y': 3.0, 'z': 4.0},
                                                      'reference_path': "#"})


if __name__ == '__main__':
    unittest.main()
