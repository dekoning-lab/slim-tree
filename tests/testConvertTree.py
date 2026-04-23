import unittest
import sys
import io
from contextlib import redirect_stderr
from utils import convertTree


class testConvertTree(unittest.TestCase):

    def test_single_branch(self):
        # A tree with one branch; after conversion the branch should be a positive integer
        result = convertTree.substitutions_to_generations("(A:0.5);", 100, 2.5e-6)
        # Extract the number from the result
        import re
        nums = re.findall(r'\d+', result)
        self.assertEqual(len(nums), 1)
        self.assertGreater(int(nums[0]), 0)

    def test_result_is_integers(self):
        # All branch lengths must be integers (no decimal points) after conversion
        result = convertTree.substitutions_to_generations("(A:0.1,B:0.2)C:0.3;", 100, 2.5e-6)
        import re
        self.assertFalse(re.search(r'\d+\.\d+', result), "Branch lengths should be integers after conversion")

    def test_proportions_preserved(self):
        # Branch of length 0.2 should be exactly twice the length of 0.1
        result = convertTree.substitutions_to_generations("(A:0.1,B:0.2);", 500, 1e-6)
        import re
        nums = [int(x) for x in re.findall(r'\d+', result)]
        self.assertEqual(len(nums), 2)
        # Allow ±1 rounding error
        self.assertAlmostEqual(nums[1] / nums[0], 2.0, delta=0.1)

    def test_larger_population_gives_longer_branches(self):
        # More generations expected with larger N (longer coalescent time)
        import re
        r1 = convertTree.substitutions_to_generations("(A:0.1,B:0.1);", 100, 2.5e-6)
        r2 = convertTree.substitutions_to_generations("(A:0.1,B:0.1);", 10000, 2.5e-6)
        nums1 = sum(int(x) for x in re.findall(r'\d+', r1))
        nums2 = sum(int(x) for x in re.findall(r'\d+', r2))
        self.assertGreater(nums2, nums1)

    def test_higher_mutation_rate_gives_shorter_branches(self):
        # Higher mutation rate means fewer generations per substitution
        import re
        r1 = convertTree.substitutions_to_generations("(A:0.1,B:0.1);", 100, 1e-7)
        r2 = convertTree.substitutions_to_generations("(A:0.1,B:0.1);", 100, 1e-4)
        nums1 = sum(int(x) for x in re.findall(r'\d+', r1))
        nums2 = sum(int(x) for x in re.findall(r'\d+', r2))
        self.assertGreater(nums1, nums2)

    def test_newick_structure_preserved(self):
        # Parentheses, commas, semicolons, and leaf names must be unchanged
        result = convertTree.substitutions_to_generations("(A:0.1,B:0.2)C:0.3;", 100, 2.5e-6)
        self.assertIn("(", result)
        self.assertIn(",", result)
        self.assertIn(";", result)
        self.assertIn("A:", result)
        self.assertIn("B:", result)
        self.assertIn("C:", result)

    def test_no_decimal_branch_lengths_exits(self):
        # A tree with no decimal branch lengths should print an error and exit
        with redirect_stderr(io.StringIO()) as serr:
            with self.assertRaises(SystemExit) as cm:
                convertTree.substitutions_to_generations("(A:1,B:2);", 100, 2.5e-6)
        self.assertEqual(cm.exception.code, 1)
        self.assertIn("decimal branch lengths", serr.getvalue())


if __name__ == '__main__':
    unittest.main()
