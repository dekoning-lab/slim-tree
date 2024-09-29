import unittest
from unittest import mock
from utils import findFitness
from contextlib import redirect_stdout
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
        
    
    # def test_find_optimal_fitness(self):
        
        # #Run test with a smaller stationary distribution that totals 1
        # self.fit.stationary_dist_file = self.test_file_path + "table_stationary_dists_partial_total_1.csv"
        # self.fit.find_optimal_fitnesses(8e-4, 100, False, None, None)
        
        # #Test that stationary distribution file is created
        # path = pathlib.Path(os.getcwd() + "/table_fitness_dists.csv")
        # self.assertEqual((str(path), path.is_file()), (str(path), True))
        # os.remove(path)
        
        # #Test that desired output is obtained
        # self.assertEqual(self.fit.fitness_mat.round(6).values.tolist(),  
                    # [[0.986429, 0.837947, 0.837947, 0.994507],
                    # [0.989663, 0.843440, 0.843440, 1.000000],
                    # [1.000000, 1.000000, 1.000000, 0.994507],
                    # [0.986429, 0.837947, 0.837947, 0.994507]])

    
    
    def test_find_optimal_fitness_mu_mat(self):
    
        #Run test with a smaller stationary distribution that totals 1, with mutation matrix
        mutation_matrix = pd.read_csv(self.test_file_path + "/mut_mat_different.csv", header = None)
        ave_mu = self.fit.find_optimal_fitnesses_mu_mat(mutation_matrix, 5, False, None, None, True)
        self.assertEqual(1.07875e-2, ave_mu)
        
      
      
    def test_process_fitness_dists(self):
    
        #Set up with the correct fitness distribution and fitness file
        self.fit.stationary_dist_file = self.test_file_path + "/table_stationary_dists_full.csv"
        good_fitness_file = self.test_file_path + "table_fitness_dists_full.csv"
        correct_fitness_mat = pd.read_csv(good_fitness_file, header = None, index_col = 0).values.tolist()
        self.fit.process_existing_fitness_file(self.test_file_path +  "table_fitness_dists_full.csv")
        
        #Compare output of process_fitness_dists
        profile_process_output = self.fit.process_fitness_dists()
        self.assertEqual(profile_process_output[0]['A'], 
                    [0.98772671022561, 0.977616132881315, 0.983892402711877, 1.0, 0.848840724523026, 0.981590226576898, 0.994620506976017, 
                        0.989837430922411, 0.996333910360753, 0.987241361388325, 0.98617271805462, 0.985665017801387, 0.982279583261853, 0.991587326200583,
                        0.983254343103396, 0.994578628643976, 1.0, 0.995395066625556, 1.0, 0.996059148953552, 0.990079894431839, 0.994368350336289, 
                        0.993291302786467, 0.999262810460322, 0.994478816655827, 0.995814365642572, 0.993470535463265, 0.99695275234019, 0.994479280421373, 
                        1.00173197443177, 0.99567698287822, 0.988184007763769, 0.995845153305058, 0.994957926138398, 1.0, 0.990221441426225, 0.988921472309081, 1.0, 
                        0.989070994223239, 0.990996003301216, 0.997290358765916, 1.00174964010603, 1.0, 0.991241038550079, 0.99276226457365, 0.990465979336448, 
                        0.991433398314745, 0.994793137690509, 0.996515820049071, 1.0, 1.0])
        
        self.assertEqual(round(profile_process_output[1],9), 0.830820542)


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
            with redirect_stdout(io.StringIO()) as sout1:
                self.fit.define_fitness_profiles( False, np.array([[0,49]]), 30)
        
        self.assertEqual(cm1.exception.code, 0)
        self.assertEqual(sout1.getvalue(), ("Please ensure that when using a fasta file, the same number of fitness profiles are provided as the length " +
                "of the genome in the fasta file. Exiting.\n"))
        sout1.close()
        
        
        # Test with a longer sequence
        with self.assertRaises(SystemExit) as cm2:
            with redirect_stdout(io.StringIO()) as sout2:
                self.fit.define_fitness_profiles( False, np.array([[0,49]]), 60)
        
        self.assertEqual(cm2.exception.code, 0)
        self.assertEqual(sout2.getvalue(), ("Please ensure that when using a fasta file, the same number of fitness profiles are provided as the length " +
                "of the genome in the fasta file. Exiting.\n"))
        sout2.close()
    

    
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
            with redirect_stdout(io.StringIO()) as sout1:
                self.fit.find_ancestral_fasta(self.test_file_path + "not_fasta.fas")
        
        self.assertEqual(cm1.exception.code, 0)
        self.assertEqual(sout1.getvalue(), ("Please provide ancestral sequence file in fasta format. Exiting.\n"))
        sout1.close()
        
        
        #Test with a sequence in terms of amino acids
        with self.assertRaises(SystemExit) as cm1:    
            with redirect_stdout(io.StringIO()) as sout1:
                self.fit.find_ancestral_fasta(self.test_file_path + "amino.fas")
        
        self.assertEqual(cm1.exception.code, 0)
        self.assertEqual(sout1.getvalue(), ("Please ensure that your fasta file is in terms of nucleotides, not amino acids. Exiting.\n"))
        sout1.close()
        
        
        #Test a fasta file with more than one sequence
        with self.assertRaises(SystemExit) as cm1:    
            with redirect_stdout(io.StringIO()) as sout1:
                self.fit.find_ancestral_fasta(self.test_file_path + "2_seq.fas")
        
        self.assertEqual(cm1.exception.code, 0)
        self.assertEqual(sout1.getvalue(), ("You have provided more than one sequence in your fasta file. Please provide only one sequence. Exiting.\n"))
        sout1.close()
        
        
        
    def test_find_fitness_scaling(self):
        #Set up with the correct fitness distribution and fitness file
        self.fit.stationary_dist_file = self.test_file_path + "/table_stationary_dists_full.csv"
        good_fitness_file = self.test_file_path + "table_fitness_dists_full.csv"
        correct_fitness_mat = pd.read_csv(good_fitness_file, header = None, index_col = 0).values.tolist()
        self.fit.process_existing_fitness_file(self.test_file_path +  "table_fitness_dists_full.csv")
        profile_process_output = self.fit.process_fitness_dists()
        
        
        #Test with only one profile (the first profile)
        self.assertEqual(0.99585076, round(self.fit.find_fitness_scaling([0], False),8))
        
        #Test with one profile for a longer genome
        self.assertEqual(round(0.99585076*5,8), round(self.fit.find_fitness_scaling([0,0,0,0,0], False),8))
        
        #Test with 1 gene of two profiles
        self.assertEqual(round(0.99585076 + 0.994762035,8), round(self.fit.find_fitness_scaling([1,0], False),8))
        
        #Test with multiple genes
        self.assertEqual(round(1+0.99585076 + 0.994762035+1+1+1+1+0.99585076 + 0.994762035+1,7), round(self.fit.find_fitness_scaling([50,1,0,50,50,50,50,1,0,50], False),8))
        
        