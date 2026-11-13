import unittest
from unittest import mock
from utils import findFitness
from contextlib import redirect_stderr
import pandas as pd
import numpy as np
import os, io, pathlib

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
            with redirect_stderr(io.StringIO()) as serr:
                self.fit.validify_stationary_distribution(pd.read_csv(self.test_file_path + "table_stationary_dists_wrong_names.csv", header = None, index_col = 0))
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("Please ensure the first row of your stationary distributions is the codon names. Exiting.\n"))
        serr.close()
        
        
        #Test a matrix that has invalid codons
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.fit.validify_stationary_distribution(pd.read_csv(self.test_file_path + "table_stationary_dists_bad_codons.csv", header = None, index_col = 0))
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("Please ensure that your stationary distribution only has valid codons. Exiting.\n"))
        serr.close()
        
        
        #Test a matrix wil stop codons
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.fit.validify_stationary_distribution(pd.read_csv(self.test_file_path + "table_stationary_dists_stop_codons.csv", header = None, index_col = 0))
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("Do not include stop codons in your stationary distribution. Exiting.\n"))
        serr.close()
        
        
        #Test a matrix without all codonswith self.assertRaises(SystemExit) as cm:
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.fit.validify_stationary_distribution(pd.read_csv(self.test_file_path + "table_stationary_dists_partial.csv", header = None, index_col = 0))
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("Please ensure that every codon is represented in your stationary distributions. Exiting.\n"))
        serr.close()
    
    
    
    def test_apply_stop_model(self):
        # site model: stop value = lowest amino-acid value - S_stop / (4N - 2), derived from the 20 amino-acid rows only,
        # so applying it twice gives the same table
        self.fit.process_existing_fitness_file(self.test_file_path + "table_fitness_dists_full.csv")
        aa_rows = [a for a in self.fit.fitness_mat.index if a != "X"]
        expected = self.fit.fitness_mat.loc[aa_rows].min(axis = 0) - 100.0 / (4 * 100 - 2)
        self.fit.apply_stop_model("site", 100.0, 100)
        self.fit.apply_stop_model("site", 100.0, 100)
        self.assertTrue(np.allclose(self.fit.fitness_mat.loc["X"].values, expected.values))

        # gene model: a stop value above a profile's lowest amino-acid value is refused
        self.fit.process_existing_fitness_file(self.test_file_path + "table_fitness_dists_full.csv")
        self.fit.fitness_mat.loc["X"] = self.fit.fitness_mat.loc[aa_rows].max(axis = 0)
        with self.assertRaises(SystemExit):
            with redirect_stderr(io.StringIO()):
                self.fit.apply_stop_model("gene", 100.0, 100)

    #Tests both process_existing_fitness_file and validify_fitness_file
    def test_process_existing_fitness_file(self):
        
        #Test a fitness file that works
        good_fitness_file = self.test_file_path + "table_fitness_dists_full.csv"
        correct_fitness_mat = pd.read_csv(good_fitness_file, header = None, index_col = 0).values.tolist()
        self.fit.process_existing_fitness_file(self.test_file_path +  "table_fitness_dists_wrong_order.csv")
        self.assertEqual(self.fit.fitness_mat.values.tolist(), correct_fitness_mat)
        
        #Test a fitness file with too many profiles
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.fit.process_existing_fitness_file(self.test_file_path + "table_fitness_dists_too_few.csv")
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("The same number of fitness profiles and stationary distributions must be provided. Exiting.\n"))
        serr.close()
        
        
        #Test a fitness file with too few profiles
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.fit.process_existing_fitness_file(self.test_file_path + "table_fitness_dists_too_many.csv")
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("The same number of fitness profiles and stationary distributions must be provided. Exiting.\n"))
        serr.close()
        
        
        #Test a fitness file with the incorrect number of amino acids
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.fit.process_existing_fitness_file(self.test_file_path + "table_fitness_dists_too_few_amino.csv")
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("Fitness data files must be in terms of amino acids. There should be the 20 amino acids and stops (ie. 21 rows). Exiting.\n"))
        serr.close()
        
        
        #Test a fitness file with an incorrect amino acid
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.fit.process_existing_fitness_file(self.test_file_path + "table_fitness_dists_bad_amino.csv")
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("Fitness data files must be in terms of amino acids. There should be the 20 amino acids and stops (ie. 21 rows). Exiting.\n"))
        serr.close()
        
      
      
    def test_process_fitness_dists(self):
    
        #Set up with the correct fitness distribution and fitness file
        self.fit.stationary_dist_file = self.test_file_path + "/table_stationary_dists_full.csv"
        good_fitness_file = self.test_file_path + "table_fitness_dists_full.csv"
        correct_fitness_mat = pd.read_csv(good_fitness_file, header = None, index_col = 0).values.tolist()
        self.fit.process_existing_fitness_file(good_fitness_file)
        
        #Compare output of process_fitness_dists
        profile_process_output = self.fit.process_fitness_dists()
        self.assertEqual(profile_process_output[0]['A'], 
                    [0.993863393422138, 0.99229128101363, 0.9919459893795, 1.0, 0.945899647493657, 0.990795072439949, 0.997310399791655, 
                        0.994918819337198, 0.998167004338997, 0.993620676778232, 0.993086442469101, 0.996315685918207, 0.991139770732561, 0.995793891396369,
                        0.991627372437062, 0.997289366704523, 1.0, 0.997697631580694, 1.0, 0.998029709598833, 0.995040054604082, 0.997184312204647, 
                        0.996645719887739, 0.999631500297592, 0.997239493334031, 0.997907180213274, 0.996735376249632, 0.998476787369641, 0.997239683109686, 
                        1.00086616162111, 0.997838593411436, 0.994092002221287, 0.99792272174356, 0.997479051848255, 1.0, 0.995110789341356, 0.997943938004739, 1.0, 
                        0.994535488030139, 0.995498047028917, 0.998645502514295, 1.00087497092797, 1.0, 0.995620646428902, 0.996381143314001, 0.995233133222628, 
                        0.995716711753162, 0.997396673692362, 0.998258242840597, 1.0, 1.0])
        
        self.assertEqual(round(profile_process_output[1],9), 0.945899647)


    def test_define_fitness_profiles(self):
        
        #Set up with the correct fitness distribution and fitness file
        self.fit.stationary_dist_file = self.test_file_path + "/table_stationary_dists_full.csv"
        good_fitness_file = self.test_file_path + "table_fitness_dists_full.csv"
        correct_fitness_mat = pd.read_csv(good_fitness_file, header = None, index_col = 0).values.tolist()
        self.fit.process_existing_fitness_file(self.test_file_path +  "table_fitness_dists_full.csv")
        profile_process_output = self.fit.process_fitness_dists()
        
        
        #Test with random fitness profiles all coding one coding area
        assigned_fit_profiles = self.fit.define_fitness_profiles( True, np.array([[0,29]]), 30)
        self.assertEqual(len(assigned_fit_profiles), 30)
        self.assertEqual(assigned_fit_profiles[0],50) #50 is the neutral codon
        self.assertTrue(all(x < 50 for x in assigned_fit_profiles[1:29]))
        self.assertEqual(assigned_fit_profiles[29], 50)
        
        
        #Test with random fitness profiles with multiple coding areas
        assigned_fit_profiles = self.fit.define_fitness_profiles( True, np.array([[0,9],
        [19,29]]), 30)
        self.assertEqual(len(assigned_fit_profiles), 30)
        self.assertEqual(assigned_fit_profiles[0],50) #50 is the neutral codon
        self.assertEqual(assigned_fit_profiles[9], 50)
        self.assertTrue(all(x < 50 for x in assigned_fit_profiles[1:9]))
        self.assertTrue(all(x == 50 for x in assigned_fit_profiles[10:19]))
        self.assertEqual(assigned_fit_profiles[19], 50)
        self.assertTrue(all(x < 50 for x in assigned_fit_profiles[20:29]))
        self.assertEqual(assigned_fit_profiles[29], 50)

        
        #Test with a user provided sequence that is the same length as the number of profiles
        assigned_fit_profiles = self.fit.define_fitness_profiles( False, np.array([[0,49]]), 50)
        self.assertEqual(list(range(50)), assigned_fit_profiles)
        
        #Test with a shorter sequence with self.assertRaises(SystemExit) as cm:
        with self.assertRaises(SystemExit) as cm1:    
            with redirect_stderr(io.StringIO()) as serr1:
                self.fit.define_fitness_profiles( False, np.array([[0,49]]), 30)
        
        self.assertEqual(cm1.exception.code, 1)
        self.assertEqual(serr1.getvalue(), ("Please ensure that when using a fasta file, the same number of fitness profiles are provided as the length " +
                "of the genome in the fasta file. Exiting.\n"))
        serr1.close()
        
        
        # Test with a longer sequence
        with self.assertRaises(SystemExit) as cm2:
            with redirect_stderr(io.StringIO()) as serr2:
                self.fit.define_fitness_profiles( False, np.array([[0,49]]), 60)
        
        self.assertEqual(cm2.exception.code, 1)
        self.assertEqual(serr2.getvalue(), ("Please ensure that when using a fasta file, the same number of fitness profiles are provided as the length " +
                "of the genome in the fasta file. Exiting.\n"))
        serr2.close()
    

    
    def test_find_anscestral(self):
        #Set up with the correct fitness distribution and fitness file
        self.fit.stationary_dist_file = self.test_file_path + "/table_stationary_dists_full.csv"
        good_fitness_file = self.test_file_path + "table_fitness_dists_full.csv"
        correct_fitness_mat = pd.read_csv(good_fitness_file, header = None, index_col = 0).values.tolist()
        self.fit.process_existing_fitness_file(self.test_file_path +  "table_fitness_dists_full.csv")
        profile_process_output = self.fit.process_fitness_dists()
        
    
        #Test finding ancestral sequence with one coding region
        assigned_fit_profiles = self.fit.define_fitness_profiles( True, np.array([[0,29]]), 30)
        ans_seq = self.fit.find_ancestral(np.array([[0,29]]), assigned_fit_profiles)
        self.assertEqual(ans_seq[0], "14")
        self.assertTrue(ans_seq[29] in ["48", "50", "56"])
        self.assertTrue(all(x not in ["48", "50", "56"] for x in ans_seq[1:29]))
        
        #Test finding ancestral sequence with two coding regions
        assigned_fit_profiles = self.fit.define_fitness_profiles( True, np.array([[0,9],
        [19,29]]), 30)
        ans_seq = self.fit.find_ancestral( np.array([[0,9], [19,29]]), assigned_fit_profiles)
        self.assertEqual(ans_seq[0], "14")
        self.assertTrue(ans_seq[9] in ["48", "50", "56"])
        self.assertTrue(all(x not in ["48", "50", "56"] for x in ans_seq[1:9]))
        self.assertEqual(ans_seq[19], "14")
        self.assertTrue(ans_seq[29] in ["48", "50", "56"])
        self.assertTrue(all(x not in ["48", "50", "56"] for x in ans_seq[20:29]))
        
    
    
    def test_find_anscestral_neutral(self):
        #Test to ensure correct length and values
        ans_seq = self.fit.find_ancestral_neutral(30)
        self.assertEqual(len(ans_seq), 30)
        self.assertTrue(all (int(x) <= 63 for x in ans_seq ))
        
    
    def test_find_anscestral_fasta(self):
        #Test with a good sequence -> AA seq = LINLH
        ans_seq, genome_length = self.fit.find_ancestral_fasta(self.test_file_path + "good_seq.fas")
        self.assertEqual(genome_length, 5)
        #L = 28,29, 30, 31, 60, 62
        #I = 12,13,15
        #N = 1, 3
        #H = 17,19
        self.assertTrue(str(ans_seq[0]) in ["28", "29", "30", "31", "60", "62"])
        self.assertTrue(str(ans_seq[1]) in ["12","13","15"])
        self.assertTrue(str(ans_seq[2]) in ["1", "3"])
        self.assertTrue(str(ans_seq[3]) in ["28","29", "30", "31", "60", "62"])
        self.assertTrue(str(ans_seq[4]) in ["17", "19"])
        
        #Test with a sequence that is not in fasta format
        with self.assertRaises(SystemExit) as cm1:    
            with redirect_stderr(io.StringIO()) as serr1:
                self.fit.find_ancestral_fasta(self.test_file_path + "not_fasta.fas")
        
        self.assertEqual(cm1.exception.code, 1)
        self.assertEqual(serr1.getvalue(), ("Please provide ancestral sequence file in fasta format. Exiting.\n"))
        serr1.close()
        
        
        #Test with a sequence in terms of amino acids
        with self.assertRaises(SystemExit) as cm1:    
            with redirect_stderr(io.StringIO()) as serr1:
                self.fit.find_ancestral_fasta(self.test_file_path + "amino.fas")
        
        self.assertEqual(cm1.exception.code, 1)
        self.assertEqual(serr1.getvalue(), ("Please ensure that your fasta file is in terms of nucleotides, not amino acids. Exiting.\n"))
        serr1.close()
        
        
        #Test a fasta file with more than one sequence
        with self.assertRaises(SystemExit) as cm1:    
            with redirect_stderr(io.StringIO()) as serr1:
                self.fit.find_ancestral_fasta(self.test_file_path + "2_seq.fas")
        
        self.assertEqual(cm1.exception.code, 1)
        self.assertEqual(serr1.getvalue(), ("You have provided more than one sequence in your fasta file. Please provide only one sequence. Exiting.\n"))
        serr1.close()
        
        
        
    def test_find_fitness_scaling(self):
        #Set up with the correct fitness distribution and fitness file
        self.fit.stationary_dist_file = self.test_file_path + "/table_stationary_dists_full.csv"
        good_fitness_file = self.test_file_path + "table_fitness_dists_full.csv"
        correct_fitness_mat = pd.read_csv(good_fitness_file, header = None, index_col = 0).values.tolist()
        self.fit.process_existing_fitness_file(good_fitness_file)
        profile_process_output = self.fit.process_fitness_dists()
        
        
        #Test with only one profile (the first profile)
        self.assertEqual(round(0.9979253977255378,8), round(self.fit.find_fitness_scaling([0], False),8))
        
        #Test with one profile for a longer genome
        self.assertEqual(round(0.9979253977255378*5,8), round(self.fit.find_fitness_scaling([0,0,0,0,0], False),8))
        
        #Test with 1 gene of two profiles
        self.assertEqual(round(0.9979253977255378 + 0.999351298565126,8), round(self.fit.find_fitness_scaling([1,0], False),8))
        
        #Test with multiple genes
        self.assertEqual(round(1+0.9979253977255378 + 0.999351298565126+1+1+1+1+0.9979253977255378 + 0.999351298565126+1,8), round(self.fit.find_fitness_scaling([50,1,0,50,50,50,50,1,0,50], False),8))
        
        