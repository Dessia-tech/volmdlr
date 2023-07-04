import unittest

import volmdlr


class TestFrame3D(unittest.TestCase):
    local_frame = volmdlr.Frame3D(
        volmdlr.Point3D(0.1, 0.3, 0.6),
        volmdlr.Y3D,
        volmdlr.Z3D,
        volmdlr.X3D,
    )

    def test_global_to_local_coordinates(self):
        vector_global = volmdlr.Vector3D(3, 4, 5)
        vector_local = self.local_frame.global_to_local_coordinates(vector_global)

        # Check that the converted vector has the expected coordinates
        self.assertEqual(vector_local.x, 4 - 0.3)
        self.assertEqual(vector_local.y, 5 - 0.6)
        self.assertEqual(vector_local.z, 3 - 0.1)

    def test_local_to_global_coordinates(self):
        vector_local = volmdlr.Vector3D(4 - 0.3, 5 - 0.6, 3 - 0.1)
        vector_global = self.local_frame.local_to_global_coordinates(vector_local)

        # Check that the converted vector has the expected coordinates
        self.assertEqual(vector_global.x, 3)
        self.assertEqual(vector_global.y, 4)
        self.assertEqual(vector_global.z, 5)

    def test_from_step(self):
        arguments = [' ', 14, '$', '$']
        object_dict = {14: volmdlr.O3D}
        frame = volmdlr.Frame3D.from_step(arguments, object_dict)
        self.assertEqual(frame, volmdlr.OXYZ)

        arguments = [' ', 43, 44, 45]
        object_dict = {43: volmdlr.Point3D(-5.829000000001e0, -9.909144910505e-1, 8.766383164265e-1),
                       44: volmdlr.Vector3D(0.e0, 6.625993710787e-1, 7.489740138657e-1),
                       45: volmdlr.Vector3D(-5.829000000001e0, -9.909144910505e-1, 8.766383164265e-1)
                       }
        frame = volmdlr.Frame3D.from_step(arguments, object_dict)
        self.assertEqual(frame.w, object_dict[44])
        self.assertEqual(frame.origin, object_dict[43])

        arguments = [' ', 20, 21, '$']
        object_dict = {20: volmdlr.O3D,
                       21: volmdlr.Y3D
                       }
        frame = volmdlr.Frame3D.from_step(arguments, object_dict)
        self.assertTrue(frame.w.is_close(object_dict[21]))

    def test_to_step(self):
        step_content, _ = volmdlr.OXYZ.to_step(10)
        expected_result = "#11 = CARTESIAN_POINT('',(0.000000,0.000000,0.000000));\n" \
                          "#12 = DIRECTION('',(0.000000,0.000000,1.000000));\n" \
                          "#13 = DIRECTION('',(1.000000,0.000000,0.000000));\n" \
                          "#14 = AXIS2_PLACEMENT_3D('',#11,#12,#13);\n"
        self.assertEqual(step_content, expected_result)

        step_content, _ = volmdlr.OYZX.to_step(10)
        expected_result = "#11 = CARTESIAN_POINT('',(0.000000,0.000000,0.000000));\n" \
                          "#12 = DIRECTION('',(1.000000,0.000000,0.000000));\n" \
                          "#13 = DIRECTION('',(0.000000,1.000000,0.000000));\n" \
                          "#14 = AXIS2_PLACEMENT_3D('',#11,#12,#13);\n"
        self.assertEqual(step_content, expected_result)


if __name__ == "__main__":
    unittest.main()
