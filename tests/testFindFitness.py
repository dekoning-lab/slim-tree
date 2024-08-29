import unittest
from unittest import mock
from utils import findFitness
from contextlib import redirect_stdout
import pandas as pd
import os, io

class testCladeReader(unittest.TestCase):

    def setUp(self):
        self.test_file_path = os.path.dirname(os.path.realpath(__file__))  + "/testFiles/"
    
        self.fit = findFitness.findFitness(self.test_file_path + "table_stationary_dists_full.csv", False)
        self.fit_neut = findFitness.findFitness(None, True)
       
    
    def test_validify_stationary_distributions(self):
        #Test a matrix that works
        self.fit.validify_stationary_distribution(pd.read_csv(self.test_file_path + "table_stationary_dists_full.csv", header = None, index_col = 0))
        
        #Test a matrix that does not have codon names in first row
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.fit.validify_stationary_distribution(pd.read_csv(self.test_file_path + "table_stationary_dists_wrong_names.csv", header = None, index_col = 0))
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("Please ensure the first row of your stationary distributions is the codon names. Exiting.\n"))
        sout.close()
        
        
        #Test a matrix that has invalid codons
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.fit.validify_stationary_distribution(pd.read_csv(self.test_file_path + "table_stationary_dists_bad_codons.csv", header = None, index_col = 0))
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("Please ensure that your stationary distribution only has valid codons. Exiting.\n"))
        sout.close()
        
        
        #Test a matrix wil stop codons
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.fit.validify_stationary_distribution(pd.read_csv(self.test_file_path + "table_stationary_dists_stop_codons.csv", header = None, index_col = 0))
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("Do not include stop codons in your stationary distribution. Exiting.\n"))
        sout.close()
        
        
        #Test a matrix without all codonswith self.assertRaises(SystemExit) as cm:
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.fit.validify_stationary_distribution(pd.read_csv(self.test_file_path + "table_stationary_dists_partial.csv", header = None, index_col = 0))
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("Please ensure that every codon is represented in your stationary distributions. Exiting.\n"))
        sout.close()
    
    
    
    #Tests both process_existing_fitness_file and validify_fitness_file
    def test_process_existing_fitness_file(self):
        
        #Test a fitness file that works
        good_fitness_file = self.test_file_path + "table_fitness_dists_full.csv"
        correct_fitness_mat = pd.read_csv(good_fitness_file, header = None, index_col = 0).values.tolist()
        self.fit.process_existing_fitness_file(self.test_file_path +  "table_fitness_dists_wrong_order.csv")
        self.assertEqual(self.fit.fitness_mat.values.tolist(), correct_fitness_mat)
        
        #Test a fitness file with too many profiles
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.fit.process_existing_fitness_file(self.test_file_path + "table_fitness_dists_too_few.csv")
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("The same number of fitness profiles and stationary distributions must be provided. Exiting.\n"))
        sout.close()
        
        
        #Test a fitness file with too few profiles
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.fit.process_existing_fitness_file(self.test_file_path + "table_fitness_dists_too_many.csv")
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("The same number of fitness profiles and stationary distributions must be provided. Exiting.\n"))
        sout.close()
        
        
        #Test a fitness file with the incorrect number of amino acids
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.fit.process_existing_fitness_file(self.test_file_path + "table_fitness_dists_too_few_amino.csv")
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("Fitness data files must be in terms of amino acids. There should be the 20 amino acids and stops (ie. 21 rows). Exiting.\n"))
        sout.close()
        
        
        #Test a fitness file with an incorrect amino acid
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.fit.process_existing_fitness_file(self.test_file_path + "table_fitness_dists_bad_amino.csv")
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("Fitness data files must be in terms of amino acids. There should be the 20 amino acids and stops (ie. 21 rows). Exiting.\n"))
        sout.close()
        
    
    def test_find_optimal_fitness(self):
        #Run test with a smaller fitness distribution