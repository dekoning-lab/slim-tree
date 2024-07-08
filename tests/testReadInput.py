
import unittest
from unittest import mock
from utils import readInput
from contextlib import redirect_stdout
import os, io
import numpy
import builtins
import copy
import filecmp


class testReadInput(unittest.TestCase):

    def setUp(self):
        self.input_reader = readInput.readInput()
        self.test_file_path = os.path.dirname(os.path.realpath(__file__))  + "/testFiles"
        
        
        #Set arguments to defaults with test tree
        class arguments:
            def __init__(self, file_path):
                self.input_tree = file_path + "/test_tree.txt"
                self.codon_stationary_distributions = file_path + "/table_stationary_distributions"
                self.aa_fitness_distributions = None
                self.high_performance_computing = False
                self.partition = 'apophis'
                self.time = '12:00:00'
                self.nonWF = True
                self.population_size = 100
                self.burn_in_multiplier = 10
                self.recombination_rate = 2.5e-8
                self.mutation_rate = 2.5e-6
                self.mutation_matrix = None
                self.tree_data_file = None
                self.gene_count = 1
                self.genome_length = 500
                self.coding_ratio = 1.0
                self.fasta_file = None
                self.sample_size = 'all'
                self.split_ratio = 0.5
                self.count_subs = False
                self.output_gens = False
                self.backup = False
                self.polymorphisms = False
                self.calculate_selection = False
                self.neutral_evolution = False
            
       
        self.arguments = arguments(self.test_file_path)
    
    
    def test_make_mutation_matrix(self):
       
        #Check that mutation matrices with different values can be run
        diff_vals = self.input_reader.make_mutation_matrix(self.test_file_path + "/mut_mat_different.csv")
        self.assertEqual(diff_vals[0].tolist(), 
        [[0.0e+00, 3.5e-07, 3.0e-06, 2.5e-05],[4.0e-08, 0.0e+00, 3.0e-06, 2.5e-05], [4.0e-08, 3.5e-07, 0.0e+00, 2.5e-05], [4.0e-08, 3.5e-07, 3.0e-06, 0.0e+00]])
       
        self.assertEqual(diff_vals[1], 
        'matrix(c(0.0, 3.5e-07, 3e-06, 2.5e-05, 4e-08, 0.0, 3e-06, 2.5e-05, 4e-08, 3.5e-07, 0.0, 2.5e-05, 4e-08, 3.5e-07, 3e-06, 0.0), ncol = 4, byrow = T)')
       
       
        #Check that mutation matrices with same values can be processed
        same_vals = self.input_reader.make_mutation_matrix(self.test_file_path + "/mut_mat_same.csv")
        self.assertEqual(same_vals[0].tolist(), 
        [[0.0e+00, 3.5e-07, 3.5e-07, 3.5e-07],[3.5e-07, 0.0e+00,3.5e-07, 3.5e-07], [3.5e-07, 3.5e-07, 0.0e+00, 3.5e-07], [3.5e-07, 3.5e-07, 3.5e-07, 0.0e+00]])
       
        self.assertEqual(same_vals[1], 
        'matrix(c(0.0, 3.5e-07, 3.5e-07, 3.5e-07, 3.5e-07, 0.0, 3.5e-07, 3.5e-07, 3.5e-07, 3.5e-07, 0.0, 3.5e-07, 3.5e-07, 3.5e-07, 3.5e-07, 0.0), ncol = 4, byrow = T)')
       
        #Check that system correctly closes when mutational matrix is not 4 by 4
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.input_reader.make_mutation_matrix(self.test_file_path + "/mut_mat_3by4.csv")

        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), 'Mutational matrices must be 4 by 4. Representing mutations from nucleotide to nucleotide.\n')
        sout.close()
        
        
        with self.assertRaises(SystemExit) as cm1:
            with redirect_stdout(io.StringIO()) as sout1:
                self.input_reader.make_mutation_matrix(self.test_file_path + "/mut_mat_4by3.csv")
            
        self.assertEqual(cm1.exception.code, 0)
        self.assertEqual(sout1.getvalue(), 'Mutational matrices must be 4 by 4. Representing mutations from nucleotide to nucleotide.\n')
        sout1.close()
        
        # Check that the system closes correctly when mutations to self are not 0
        with self.assertRaises(SystemExit) as cm2:
            with redirect_stdout(io.StringIO()) as sout2:
                self.input_reader.make_mutation_matrix(self.test_file_path + "/mut_mat_nonzero.csv")
            
        self.assertEqual(cm2.exception.code, 0)
        self.assertEqual(sout2.getvalue(), 'All mutations from a nucleotide to itself must be 0.\n')
        sout2.close()


    def test_check_arguments(self):
                
                
        #Test to make sure return is true when arguments good without hpc
        return_val = self.input_reader.check_arguments(self.arguments)
        self.assertTrue(return_val)
        
        #Test to make sure return is true when arguments good with fasta_file
        self.arguments.fasta_file = "my_fasta.fa"
        return_val = self.input_reader.check_arguments(self.arguments)
        self.assertTrue(return_val)
        self.arguments.fasta_file = None
        
        #Test to make sure return is true when arguments good with hpc
        self.arguments.high_performance_computing = True
        return_val = self.input_reader.check_arguments(self.arguments)
        self.assertTrue(return_val)

        # Test that arguments with hpc close if time and partition are not included
        self.arguments.partition = 'apophis'
        self.arguments.time = None
        with redirect_stdout(io.StringIO()) as sout:
            return_val = self.input_reader.check_arguments(self.arguments)
        self.assertFalse(return_val)
            
        self.assertEqual(sout.getvalue(), 'When using high performance computing, partition and time data must be provided. Closing program.\n')
        sout.close()
        

        self.arguments.partition = None
        self.arguments.time ='10:00:00'
        with redirect_stdout(io.StringIO()) as sout2:
            return_val = self.input_reader.check_arguments(self.arguments)
        self.assertFalse(return_val)
            
        self.assertEqual(sout2.getvalue(), 'When using high performance computing, partition and time data must be provided. Closing program.\n')
        sout2.close()
        
        #Test that program does not run if gene count incorrect
        self.arguments.high_performance_computing = False
        self.arguments.gene_count = -0.05
        with redirect_stdout(io.StringIO()) as sout:
            return_val = self.input_reader.check_arguments(self.arguments)
        self.assertFalse(return_val)
            
        self.assertEqual(sout.getvalue(), "Number of genes must be greater than 0 and less than the length of the genome. Closing program.\n")
        sout.close()
        
        self.arguments.gene_count = 501
        with redirect_stdout(io.StringIO()) as sout:
            return_val = self.input_reader.check_arguments(self.arguments)
        self.assertFalse(return_val)
            
        self.assertEqual(sout.getvalue(), "Number of genes must be greater than 0 and less than the length of the genome. Closing program.\n")
        sout.close()
        
        
        #Test that program does not run if coding ratio incorrect
        self.arguments.gene_count = 1
        self.arguments.coding_ratio = -0.05
        with redirect_stdout(io.StringIO()) as sout:
            return_val = self.input_reader.check_arguments(self.arguments)
        self.assertFalse(return_val)
            
        self.assertEqual(sout.getvalue(), "Coding ratio must be greater than 0 and less than or equal to 1. Please re-enter as a ratio (0, 1]. Closing program.\n")
        sout.close()
        
        self.arguments.gene_count = 1.1
        with redirect_stdout(io.StringIO()) as sout:
            return_val = self.input_reader.check_arguments(self.arguments)
        self.assertFalse(return_val)
            
        self.assertEqual(sout.getvalue(), "Coding ratio must be greater than 0 and less than or equal to 1. Please re-enter as a ratio (0, 1]. Closing program.\n")
        sout.close()
        
        
        #Test the program if there is a fasta file
        self.arguments.fasta_file = "my_fasta.fa"
        self.arguments.gene_count = 1
        self.arguments.coding_ratio = 0.5
        with redirect_stdout(io.StringIO()) as sout:
            return_val = self.input_reader.check_arguments(self.arguments)
        self.assertFalse(return_val)
            
        self.assertEqual(sout.getvalue(), "When specifying an ancestral sequence with a fasta file, the sequence of only one fully coding gene should be provided. Closing program.\n")
        sout.close()

        
        self.arguments.gene_count = 2
        self.arguments.coding_ratio = 1.0
        with redirect_stdout(io.StringIO()) as sout:
            return_val = self.input_reader.check_arguments(self.arguments)
        self.assertFalse(return_val)
            
        self.assertEqual(sout.getvalue(), "When specifying an ancestral sequence with a fasta file, the sequence of only one fully coding gene should be provided. Closing program.\n")
        sout.close()
        
        
        self.arguments.gene_count = 1
        self.arguments.coding_ratio = 1.0
        with redirect_stdout(io.StringIO()) as sout:
            return_val = self.input_reader.check_arguments(self.arguments)
        self.assertTrue(return_val)
            
        sout.close()


    def test_process_filenames(self):
        
        tree_name = "tests/testFiles/test_tree"
        
        #Test process with new directories, no backup, ensure folders created
        processed_output = self.input_reader.process_filenames(tree_name + ".txt", False)
        correct_output = (self.test_file_path + "/slimScripts/test_tree", self.test_file_path + "/test_tree", None)
        self.assertTrue(processed_output == correct_output)
        self.assertTrue(os.path.exists(self.test_file_path + "/slimScripts/"))
        self.assertTrue(os.path.exists(self.test_file_path + "/aa_FASTA/"))
        self.assertTrue(os.path.exists(self.test_file_path + "/nuc_FASTA/"))
        
        #Test process with existing directories, no backup
        with mock.patch.object(builtins, 'input', lambda _: 'y'):
            with redirect_stdout(io.StringIO()) as sout:
                processed_output2 = self.input_reader.process_filenames(tree_name + ".txt", False)
 
            sout_correct = sout.getvalue() == "using the same nuc_FASTA folder\nusing the same  aa_FASTA folder\n"
            sout.close()
            correct = processed_output2 == correct_output and sout_correct
            assert correct
        
        
        #Remove directories to test with a backup
        os.rmdir(self.test_file_path + "/slimScripts/")
        os.rmdir(self.test_file_path + "/aa_FASTA/")
        os.rmdir(self.test_file_path + "/nuc_FASTA/") 

        #Test process with new directories, include backup, ensure folders created
        processed_output3 = self.input_reader.process_filenames(tree_name + ".txt", True)
        correct_output2 = (self.test_file_path + "/slimScripts/test_tree", self.test_file_path + "/test_tree", self.test_file_path + "/backupFiles")
        self.assertTrue(processed_output3 == correct_output2)
        self.assertTrue(os.path.exists(self.test_file_path + "/slimScripts/"))
        self.assertTrue(os.path.exists(self.test_file_path + "/aa_FASTA/"))
        self.assertTrue(os.path.exists(self.test_file_path + "/nuc_FASTA/"))
        self.assertTrue(os.path.exists(self.test_file_path + "/backupFiles/"))
        
        
        # Test process with existing directories, no backup
        with mock.patch.object(builtins, 'input', lambda _: 'y'):
            with redirect_stdout(io.StringIO()) as sout:
                processed_output4 = self.input_reader.process_filenames(tree_name + ".txt", True)
 
            sout_correct = sout.getvalue() == "using the same nuc_FASTA folder\nusing the same  aa_FASTA folder\nusing same backup folder\n"
            sout_val = sout.getvalue()
            sout.close()
            correct = processed_output4 == correct_output2 and sout_correct
            assert correct
        
        # Remove directories
        os.rmdir(self.test_file_path + "/slimScripts")
        os.rmdir(self.test_file_path + "/aa_FASTA/")
        os.rmdir(self.test_file_path + "/nuc_FASTA/")
        os.rmdir(self.test_file_path + "/backupFiles/")
        
    def test_make_param_dict(self):
        
        #Set up comparison dict
        dict_arguments = copy.deepcopy (vars(self.arguments))
        dict_arguments["burn_in"] = self.arguments.burn_in_multiplier * self.arguments.population_size
        del dict_arguments["burn_in_multiplier"]
        dict_arguments["jukes_cantor"] = self.arguments.mutation_matrix == None
    
        
        # Check with default settings
        self.assertTrue(self.input_reader.make_param_dict(self.arguments) == dict_arguments)
            
        
        #Check with different burn in multiplier
        self.arguments.burn_in_multiplier = 50
        dict_arguments["burn_in"] = 50 * self.arguments.population_size
        self.assertTrue(self.input_reader.make_param_dict(self.arguments) == dict_arguments)
            
       
        #Check with different population size
        self.arguments.population_size = 1000
        dict_arguments["burn_in"] = self.arguments.burn_in_multiplier *1000
        dict_arguments["population_size"] = 1000
        self.assertTrue(self.input_reader.make_param_dict(self.arguments) == dict_arguments)
         
        
        #Check with a mutation_matrix
        self.arguments.mutation_matrix = ([[0.0e+00, 3.5e-07, 3.0e-06, 2.5e-05],
                    [4.0e-08, 0.0e+00, 3.0e-06, 2.5e-05], 
                    [4.0e-08, 3.5e-07, 0.0e+00, 2.5e-05], 
                    [4.0e-08, 3.5e-07, 3.0e-06, 0.0e+00]])
       
        dict_arguments["jukes_cantor"] = False
        dict_arguments["mutation_matrix"] = self.arguments.mutation_matrix
        self.assertTrue(self.input_reader.make_param_dict(self.arguments) == dict_arguments)
        
        
        
    def test_save_input(self):
        #Make parameter dict and file outputs - already tested
        param_dict = self.input_reader.make_param_dict(self.arguments)
        param_dict["filenames"] = self.input_reader.process_filenames("tests/testFiles/test_tree.txt", False)
        
        #Test default
        self.input_reader.save_input(param_dict)
        self.assertTrue(filecmp.cmp(self.test_file_path + "/test_tree_parameters.yaml",self.test_file_path + "/correct_parameters.yaml"))
            
        #Test with different mutation rate
        old_mu = param_dict["mutation_rate"]
        param_dict["mutation_rate"] = 1.0e-4
        self.input_reader.save_input(param_dict)
        self.assertTrue(filecmp.cmp(self.test_file_path + "/test_tree_parameters.yaml",self.test_file_path + "/correct_parameters_change_mu.yaml"))
        param_dict["mutation_rate"] = old_mu
        
        
        #Test with different population size
        old_pop = param_dict["population_size"]
        param_dict["population_size"] = 50
        self.input_reader.save_input(param_dict)
        self.assertTrue(filecmp.cmp(self.test_file_path + "/test_tree_parameters.yaml",self.test_file_path + "/correct_parameters_change_pop.yaml"))
        param_dict["population_size"] = old_pop
        
        #Test with mutation matrix - same value
        param_dict["mutation_matrix"] =  [numpy.array([[0.0e+00, 3.5e-06, 3.5e-06, 3.5e-06],
                    [3.5e-06, 0.0e+00, 3.5e-06, 3.5e-06], 
                    [3.5e-06, 3.5e-06, 0.0e+00, 3.5e-06], 
                    [3.5e-06, 3.5e-06, 3.5e-06, 0.0e+00]]), 
                    'matrix(c(0.0, 3.5e-06, 3.5e-06, 3.5e-06, 3.5e-06, 0.0, 3.5e-06, 3.5e-06, 3.5e-06, 3.5e-06, 0.0, 3.5e-06, 3.5e-06, 3.5e-06, 3.5e-06, 0.0), ncol = 4, byrow = T)']
        param_dict["jukes_cantor"] = False
        self.input_reader.save_input(param_dict)
        self.assertTrue(filecmp.cmp(self.test_file_path + "/test_tree_parameters.yaml",self.test_file_path + "/correct_parameters_same_mu_matrix.yaml"))
        
        
        #Test with mutation matrix - different values
        param_dict["mutation_matrix"] =  [numpy.array([[0.0e+00, 2.5e-06, 3.5e-06, 3.5e-04],
                    [6.5e-06, 0.0e+00, 3.5e-03, 1.5e-06], 
                    [0.5e-06, 9.0e-06, 0.0e+00, 3.5e-02], 
                    [3.5e-07, 5.5e-02, 3.0e-06, 0.0e+00]]), 
                    'matrix(c(0.0, 2.5e-06, 3.5e-06, 3.5e-04, 6.5e-06, 0.0, 3.5e-03, 1.5e-06, 0.5e-06, 9.0e-06, 0.0, 3.5e-02, 3.5e-07, 3.5e-06, 3.0e-06, 0.0), ncol = 4, byrow = T)']
        param_dict["jukes_cantor"] = False
        self.input_reader.save_input(param_dict)
        self.assertTrue(filecmp.cmp(self.test_file_path + "/test_tree_parameters.yaml",self.test_file_path + "/correct_parameters_diff_mu_matrix.yaml"))
        
        
        
        # Remove directories
        os.rmdir(self.test_file_path + "/slimScripts")
        os.rmdir(self.test_file_path + "/aa_FASTA/")
        os.rmdir(self.test_file_path + "/nuc_FASTA/")
        os.remove(self.test_file_path + "/test_tree_parameters.yaml")