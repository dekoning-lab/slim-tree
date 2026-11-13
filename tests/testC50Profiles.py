#Regression check for the C50 profiles shipped with SLiM-Tree. The bundled amino-acid table and every codon-level C50 table
#(examples and the test fixture) must reproduce Table 4 of the supplement to Le, Gascuel & Lartillot (2008), Bioinformatics
#24:2317 (tests/testFiles/le2008_c50_table4.csv). Until October 2026 these files carried the table's valine column in the
#tryptophan row and vice versa, so every simulation from them used valine-poor, tryptophan-rich profiles.
import unittest
import glob, os
import numpy as np
import pandas as pd
from Bio.Seq import Seq

AMINO_ACIDS = list("ACDEFGHIKLMNPQRSTVWY")

#The bundled table was re-derived from fitted fitness profiles and differs from the published rows by less than 2e-5;
#a valine/tryptophan transposition moves profiles by up to 0.43
TOLERANCE = 1e-4

class testC50Profiles(unittest.TestCase):

    def setUp(self):
        self.repo_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
        self.table4_counts = pd.read_csv(os.path.join(self.repo_path, "tests", "testFiles", "le2008_c50_table4.csv"),
                                         comment = "#", index_col = 0)
        self.table4 = self.table4_counts.div(self.table4_counts.sum(axis = 1), axis = 0)
        self.bundled = pd.read_csv(os.path.join(self.repo_path, "fitnessDataFiles", "table_stationary_distributions.csv"), header = None)
        self.codon_files = (sorted(glob.glob(os.path.join(self.repo_path, "Examples", "*", "stationary_distributions.csv"))) +
                            [os.path.join(self.repo_path, "tests", "testFiles", "table_stationary_dists_full.csv")])


    def test_table4_fixture(self):
        #50 profiles of integer percentages in alphabetical one-letter order, as printed
        self.assertEqual(self.table4_counts.shape, (50, 20))
        self.assertEqual(list(self.table4_counts.columns), AMINO_ACIDS)
        self.assertTrue(self.table4_counts.sum(axis = 1).between(95, 103).all())


    def test_bundled_amino_acid_table(self):
        #Rows are the 20 amino acids in alphabetical one-letter order followed by stops, columns are profiles
        self.assertEqual(self.bundled.shape, (21, 50))
        np.testing.assert_allclose(self.bundled.iloc[:20].transpose().values, self.table4.values, rtol = 0, atol = TOLERANCE)
        self.assertTrue((self.bundled.iloc[20] == 0).all())


    def test_codon_tables(self):
        #Codon tables are matched through their codon labels, so the check does not depend on row order
        self.assertEqual(len(self.codon_files), 11)
        for codon_file in self.codon_files:
            with self.subTest(codon_file = os.path.relpath(codon_file, self.repo_path)):
                codons = pd.read_csv(codon_file, header = None, index_col = 0)
                amino_acids = [*str(Seq("".join(codons.index)).translate())]
                profiles = codons.groupby(amino_acids).sum().loc[AMINO_ACIDS].transpose().values
                np.testing.assert_allclose(profiles, self.table4.values, rtol = 0, atol = TOLERANCE)


    def test_valine_tryptophan(self):
        #Profile 2 is I 41, L 8, V 43, W 0 (of 99); with valine and tryptophan transposed it is a 43% tryptophan profile
        self.assertEqual(list(self.table4_counts.loc[2, ["I", "L", "V", "W"]]), [41, 8, 43, 0])
        self.assertAlmostEqual(self.bundled.iloc[AMINO_ACIDS.index("V"), 1], 43 / 99, places = 4)
        self.assertLess(self.bundled.iloc[AMINO_ACIDS.index("W"), 1], 1e-9)
        for codon_file in self.codon_files:
            with self.subTest(codon_file = os.path.relpath(codon_file, self.repo_path)):
                codons = pd.read_csv(codon_file, header = None, index_col = 0)
                self.assertAlmostEqual(codons.loc[["GTA", "GTC", "GTG", "GTT"], 2].sum(), 43 / 99, places = 4)
                self.assertLess(codons.loc["TGG", 2], 1e-9)


if __name__ == '__main__':
    unittest.main()
