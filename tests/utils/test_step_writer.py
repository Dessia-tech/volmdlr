import unittest
from volmdlr.utils import step_writer


class TestStepIdsToStr(unittest.TestCase):
    ids_0 = [0]
    ids_1 = [1, 2, 3, 4, 5]
    ids_2 = [11, 22, 33, 44, 55]

    def test_step_ids_to_str(self):
        self.assertEqual(step_writer.step_ids_to_str(self.ids_0), "#0")
        self.assertEqual(step_writer.step_ids_to_str(self.ids_1), "#1,#2,#3,#4,#5")
        self.assertEqual(step_writer.step_ids_to_str(self.ids_2), "#11,#22,#33,#44,#55")


if __name__ == '__main__':
    unittest.main()
