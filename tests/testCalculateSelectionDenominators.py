
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
        
    
    
    
    #Numbers of synonymous and nonsynonymous sites (Goldman & Yang 1994) in n_codons codons whose summed synonymous and
    #nonsynonymous mutational opportunities are given: 3 * n_codons, split in proportion to the opportunities
    def site_numbers(self, opportunity, n_codons):
        opportunity = numpy.array(opportunity, dtype = float)
        return 3 * n_codons * opportunity / opportunity.sum()


    def test_calculate_selection_denominators(self):

        # Test a situation where we have the same profile 5 times
        correct_syn_codons = [[True,False,False,False,False],
                                [False,True,True,True,False],
                                [False,True,True,True,False],
                                [False,True,True,True,False],
                                [False,False,False,False,True]]
                                
        test_dn_ds = self.denom_calculator.calculate_selection_denominators([0,0,0,0,0], correct_syn_codons)
        numpy.testing.assert_allclose(test_dn_ds, self.site_numbers([8.838288916666669e-07, 3.3669679583333337e-07], 5), rtol = 1e-12)
        self.assertAlmostEqual(sum(test_dn_ds), 15.0, places = 10) #3 sites per codon
        
        
        # Test with 2 different profiles
        test_dn_ds = self.denom_calculator.calculate_selection_denominators([0,1,1,1,0], correct_syn_codons)
        numpy.testing.assert_allclose(test_dn_ds, self.site_numbers([3.535315566667927e-07, 1.3467871833337535e-07], 5), rtol = 1e-12)
        
        
        # Test with all 3 profiles
        test_dn_ds = self.denom_calculator.calculate_selection_denominators([0,0,1,1,2,2], correct_syn_codons)
        numpy.testing.assert_allclose(test_dn_ds, self.site_numbers([3.5353155666683473e-07, 1.3467871833338936e-07], 6), rtol = 1e-12)
        
        # Test calculation of dn/ds when mutation rate varies
        self.denom_calculator.mu_mat = numpy.array(
                            [[0.0, 0.0, 0.0, 2e-6, 0.0],
                            [0.0, 0.0, 5e-6, 6e-6, 0.0],
                            [0.0, 8e-6, 0.0, 9e-6, 0.0],
                            [7e-6, 2e-5, 3e-5, 0.0, 0.0],
                            [0.0, 0.0, 0.0, 0.0, 0.0]])
        test_dn_ds = self.denom_calculator.calculate_selection_denominators([0,0,0,0,0], correct_syn_codons)
        numpy.testing.assert_allclose(test_dn_ds, self.site_numbers([2.303004757e-06, 2.87875642e-07], 5), rtol = 1e-12)


    def test_calculate_selection_denominators_unused_profile(self):
        correct_syn_codons = [[True,False,False,False,False],
                                [False,True,True,True,False],
                                [False,True,True,True,False],
                                [False,True,True,True,False],
                                [False,False,False,False,True]]

        # Sites must be matched to their own profile when a lower-numbered profile is not used at all
        ds_dn = [self.denom_calculator.get_dist_ds_dn(profile, 1, correct_syn_codons) for profile in range(3)]

        test_dn_ds = self.denom_calculator.calculate_selection_denominators([2, 2, 2], correct_syn_codons)
        numpy.testing.assert_allclose(test_dn_ds, self.site_numbers(ds_dn[2], 3), rtol = 1e-12)

        test_dn_ds = self.denom_calculator.calculate_selection_denominators([0, 2, 2], correct_syn_codons)
        numpy.testing.assert_allclose(test_dn_ds, self.site_numbers(numpy.array(ds_dn[0]) + 2 * numpy.array(ds_dn[2]), 3), rtol = 1e-12)

        # The order in which sites are listed does not matter
        test_dn_ds = self.denom_calculator.calculate_selection_denominators([2, 0, 1, 2], correct_syn_codons)
        numpy.testing.assert_allclose(test_dn_ds,
                self.site_numbers(numpy.array(ds_dn[0]) + numpy.array(ds_dn[1]) + 2 * numpy.array(ds_dn[2]), 4), rtol = 1e-12)

        # No counted sites gives no mutational opportunity
        self.assertEqual(list(self.denom_calculator.calculate_selection_denominators([], correct_syn_codons)), [0.0, 0.0])


    def test_sites_restrict_denominators(self):
        # Only the listed codon positions contribute: here the first and last codons are left out
        all_sites = calculateSelectionDenominators.calculateSelectionDenominators(
                self.partial_stationary_dist, [1, 2], 2.5e-06, None)
        inner_sites = calculateSelectionDenominators.calculateSelectionDenominators(
                self.partial_stationary_dist, [0, 1, 2, 0], 2.5e-06, None, sites = [1, 2])
        self.assertEqual(inner_sites.get_dn(), all_sites.get_dn())
        self.assertEqual(inner_sites.get_ds(), all_sites.get_ds())


    def test_selection_sites(self):
        selection_sites = calculateSelectionDenominators.selection_sites

        # One gene over the whole genome: its start and terminal stop codons are left out
        self.assertEqual(selection_sites([[0, 5]], 6), [1, 2, 3, 4])

        # Two genes with non-coding codons between and after them
        self.assertEqual(selection_sites(numpy.array([[0, 3], [6, 9]]), 12), [1, 2, 7, 8])

        # A gene with no codon between its start and stop has nothing to count
        self.assertEqual(selection_sites([[0, 1]], 2), [])

        # Neutral runs have no genes, so every codon counts
        self.assertEqual(selection_sites([[0, 5]], 6, neutral_evolution = True), [0, 1, 2, 3, 4, 5])

