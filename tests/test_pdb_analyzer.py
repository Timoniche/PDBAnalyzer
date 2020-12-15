import unittest

from PDBAnalyzer import PDBAnalyzer


class TestCase(unittest.TestCase):
    file_name = '3DMax_chr1.pdb'
    bins_cnt = 1000

    def test_curve_length(self):
        # 0 2 3
        # 2 0 7
        # 3 7 0
        test_dist_matrix = [[0, 2, 3], [2, 0, 7], [3, 7, 0]]
        curve_length = PDBAnalyzer.pdb_curve_length(test_dist_matrix, 3)
        self.assertEqual(9, curve_length)


if __name__ == '__main__':
    unittest.main()