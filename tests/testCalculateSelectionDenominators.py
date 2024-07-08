
import unittest, os
import pandas as pd
import numpy
from unittest import mock
from utils import calculateSelectionDenominators

class testCalculateSelectionDenominators(unittest.TestCase):

    def setUp(self):
        self.test_file_path = os.path.dirname(os.path.realpath(__file__))  + "/testFiles"
        
        #Read in associated test files
        self.partial_stationary_dist = pd.read_csv(self.test_file_path + "/table_stationary_dists_partial.csv", header = None, index_col = 0)
       
        #Set up class, input is just to set up, will test functions individually
        self.denom_calculator = calculateSelectionDenominators.calculateSelectionDenominators(
                self.partial_stationary_dist, [1], 2.5e-06, None)
        
    def test_find_syn_codons(self):
        
        #Test whether the synonymous codons that are found are correct
        syn_codons = self.denom_calculator.find_syn_codons(self.partial_stationary_dist)
        correct_syn_codons = [[True,False,False,False,False],
                                [False,True,True,True,False],
                                [False,True,True,True,False],
                                [False,True,True,True,False],
                                [False,False,False,False,True]]
        self.assertEqual(syn_codons, correct_syn_codons)
        
        
    def test_find_rate_mut(self):
        #Test with a Jukes-Cantor matrix where all mutation rates are the same, want 0 when >1 mutation and no mutations
        mu_mat = numpy.array([[0.0e+00, 3.5e-06, 3.5e-06, 3.5e-06],
                    [3.5e-06, 0.0e+00, 3.5e-06, 3.5e-06], 
                    [3.5e-06, 3.5e-06, 0.0e+00, 3.5e-06], 
                    [3.5e-06, 3.5e-06, 3.5e-06, 0.0e+00]])
                    
        
        calculated_mut_rates = self.denom_calculator.find_rate_mut(mu_mat)
        correct_mut_rates = [[0.0, 0.0, 0.0, 3.5e-6, 0.0],
                            [0.0, 0.0, 3.5e-6, 3.5e-6, 0.0],
                            [0.0, 3.5e-6, 0.0, 3.5e-6, 0.0],
                            [3.5e-6, 3.5e-6, 3.5e-6, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0]]
                            
        self.assertEqual(calculated_mut_rates.tolist(), correct_mut_rates)
        
        #Test with a matrix where mutation rates differ
        mu_mat2 = numpy.array([[0.0e+00, 1e-06, 2e-06, 3e-06],
                    [4e-06, 0.0e+00, 5e-06, 6e-06], 
                    [7e-06, 8e-06, 0.0e+00, 9e-06], 
                    [1e-05, 2e-05, 3e-05, 0.0e+00]])
                    
                    
        calculated_mut_rates2 = self.denom_calculator.find_rate_mut(mu_mat2)
        correct_mut_rates2 = [[0.0, 0.0, 0.0, 2e-6, 0.0],
                            [0.0, 0.0, 5e-6, 6e-6, 0.0],
                            [0.0, 8e-6, 0.0, 9e-6, 0.0],
                            [7e-6, 2e-5, 3e-5, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0]]
                            
                          
        self.assertEqual(calculated_mut_rates2.tolist(), correct_mut_rates2)
        
        
    def test_get_dist_ds_dn(self):
        correct_syn_codons = [[True,False,False,False,False],
                                [False,True,True,True,False],
                                [False,True,True,True,False],
                                [False,True,True,True,False],
                                [False,False,False,False,True]]
        
        #Test calculation of dn/ds denominator on one profile used only once
        test_dn_ds = self.denom_calculator.get_dist_ds_dn(0,1, correct_syn_codons)
        self.assertEqual(test_dn_ds,  [1.7676577833333335e-07, 6.733935916666667e-08])
        
        # Test calculation of dn/ds denominator on one profile used 5 times
        test_dn_ds = self.denom_calculator.get_dist_ds_dn(0,5, correct_syn_codons)
        self.assertEqual(test_dn_ds, [8.838288916666669e-07, 3.3669679583333337e-07])
        
        # Test calculation of dn/ds when mutation rate varies
        self.denom_calculator.mu_mat = numpy.array(
                            [[0.0, 0.0, 0.0, 2e-6, 0.0],
                            [0.0, 0.0, 5e-6, 6e-6, 0.0],
                            [0.0, 8e-6, 0.0, 9e-6, 0.0],
                            [7e-6, 2e-5, 3e-5, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0]])
        test_dn_ds = self.denom_calculator.get_dist_ds_dn(0,1, correct_syn_codons)
        self.assertEqual(test_dn_ds, [2.303004757e-06, 2.87875642e-07])
        
    
    
    
    def test_calculate_selection_denominators(self):

        # Test a situation where we have the same profile 5 times
        correct_syn_codons = [[True,False,False,False,False],
                                [False,True,True,True,False],
                                [False,True,True,True,False],
                                [False,True,True,True,False],
                                [False,False,False,False,True]]
                                
        test_dn_ds = self.denom_calculator.calculate_selection_denominators([0,0,0,0,0], correct_syn_codons)
        self.assertEqual(list(test_dn_ds), [8.838288916666669e-07, 3.3669679583333337e-07])
        
        
        # Test with 2 different profiles
        test_dn_ds = self.denom_calculator.calculate_selection_denominators([0,1,1,1,0], correct_syn_codons)
        self.assertEqual(list(test_dn_ds),[3.535315566667927e-07, 1.3467871833337535e-07])
        
        
        # Test with all 3 profiles
        test_dn_ds = self.denom_calculator.calculate_selection_denominators([0,0,1,1,2,2], correct_syn_codons)
        self.assertEqual(list(test_dn_ds), [3.5353155666683473e-07, 1.3467871833338936e-07])
        
        # Test when mutation rate varies
        # Test calculation of dn/ds when mutation rate varies
        self.denom_calculator.mu_mat = numpy.array(
                            [[0.0, 0.0, 0.0, 2e-6, 0.0],
                            [0.0, 0.0, 5e-6, 6e-6, 0.0],
                            [0.0, 8e-6, 0.0, 9e-6, 0.0],
                            [7e-6, 2e-5, 3e-5, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0]])
        test_dn_ds = self.denom_calculator.calculate_selection_denominators([0,0,0,0,0], correct_syn_codons)
        self.assertEqual(list(test_dn_ds), [2.303004757e-06*5, 2.87875642e-07*5])
    