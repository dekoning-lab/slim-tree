import unittest
import os, io, yaml, copy
from unittest import mock
from utils import cladeReader, findFitness
from contextlib import redirect_stderr
from Bio import Phylo
import numpy as np

class testCladeReader(unittest.TestCase):

    def setUp(self):
        self.test_file_path = os.path.dirname(os.path.realpath(__file__))  + "/testFiles/"
        with open(self.test_file_path + "correct_parameters.yaml", 'r') as yaml_file:
            starting_params = yaml.safe_load(yaml_file)
        yaml_file.close()
        starting_params["input_tree"] = self.test_file_path + starting_params["input_tree"]
        self.input_tree = starting_params["input_tree"]
        
        #Set up for a fitness profile shift
        starting_params['fitness_profile_nums'] = [50, 6, 12, 35, 1, 5, 30, 42, 6, 50]
        starting_params['coding_seqs'] = [[0, 9]]
        starting_params['fitness_finder'] = findFitness.findFitness(self.test_file_path + "table_stationary_dists_full.csv", False)
        good_fitness_file = self.test_file_path + "table_fitness_dists_full.csv"
        starting_params['fitness_finder'].process_existing_fitness_file(self.test_file_path +  "table_fitness_dists_full.csv")
        starting_params['fitness_finder'].process_fitness_dists()
        starting_params["coding_ratio"] = 1.0
        starting_params["genome_length"] = len(starting_params['fitness_profile_nums'])

        self.read_clades = cladeReader.cladeReader(starting_params)
        

    def test_read_clade_data(self):
 
        #Test failure if non-yaml data file is used
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("Please make sure your changes are in yaml format. " +
                        "For more information on yaml format visit: " +
                        "https://en.wikipedia.org/wiki/YAML. Program closing.\n"))
        serr.close()
        
        #Test successful yaml file with abbreviated names
        yaml_output = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_abbr.yaml"])
        self.assertEqual({"A":{"population_size" : 50}}, yaml_output)
        
        #Test successful yaml file with full names
        yaml_output = self.read_clades.read_clade_data([self.test_file_path + "good_yaml.yaml"])
        self.assertEqual({"A":{"population_size" : 50}}, yaml_output)
        
        #Test failure of yaml when no branch specified
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_no_branch.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("Please check the formatting of your yaml file, you likely did not include the branch name. Closing program.\n"))
        serr.close()
        
        #Test failure with yaml file with changes that cannot be made
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_wrong_vals.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("When using slim-tree without HPC, only the population size may be modified for specific branches. Exiting\n"))
        serr.close()
        
        #Test successful hpc yaml file
        self.read_clades.start_params["high_performance_computing"] = True
        yaml_output = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_abbr_hpc.yaml"])
        self.assertEqual({"A":{"population_size" : 50, "mutation_rate" : 0.0025, "jukes_cantor" : True, "partition": "anothernode", "time":7200, "recombination_rate": 0.05, "sample_size": 30, "split_ratio":0.4}}, yaml_output)
        
        
        #Test hpc yaml file with changes that cannot be made
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_abbr_hpc.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("When using slim-tree with HPC only the following parameters may be modified for specific branches:\n" + 
                        "\n".join(["partition", "time", "population_size", "recombination_rate", "mutation_rate", "mutation_matrix", "sample_size", "split_ratio", "profile_shift\n"] )))
        serr.close()
        
        #Test successful profile shift yaml with good changes
        yaml_output = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_profile_shift.yaml"])
        self.assertEqual({"A":{'fitness_profile_nums': [50, 6, 15, 40, 1, 5, 30, 42, 6, 50],
                    "profile_shift":{"profile_positions" : [2,3], "new_profile_nums": [15,40]},'scaling_value': 9.967312178527093,
                    'shift': True}}, yaml_output)
        
        
        #Test profile shift yaml with changes that cannot be made
        #Incorrect profile names
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_profile_shift_wrong_profile_names.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), "Please include new profile numbers as a list of integer new profile numbers. Exiting.\n")
        serr.close()
        
        #Incorrect names of profiles to shift
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_profile_shift_wrong_shift_names.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), "Please include profiles to shift as a list of integer profile numbers. Exiting.\n")
        serr.close()
        
        #Profiles to shift outside the genome
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_profile_shift_outside_genome.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), "Please ensure that all your profiles to shift are within the genome. Exiting.\n")
        serr.close()
        
        #Profiles to shift are start or stop codon
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_profile_shift_stop_codon.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("Please ensure that all your profiles to shift are within coding regions of your given " +
                        "genome and are not in start or stop codon positions. Exiting.\n"))
        serr.close()
        
        #Different number of profiles as shifts
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_profile_shift_diff_num.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), "Please ensure that you provide the same number of new profile numbers as profiles. Exiting.\n")
        serr.close()
        
        #shifts != integers
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_profile_shift_profile_not_int.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), "Please include profiles to shift as a list of integer profile numbers. Exiting.\n")
        serr.close()
        
        #Profiles != integers
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_profile_shift_shift_not_int.yaml"])
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), "Please include new profile numbers as a list of integer new profile numbers. Exiting.\n")
        serr.close()
        
        self.read_clades.start_params["high_performance_computing"] = False
        
        
        
        
    def test_read_input_trees(self):
    
        #If something fails in this test, can also be a failure of recurse through clades or get_clade_data
    
        #Check a correct tree
        tree_output = copy.deepcopy(self.read_clades.read_input_tree())
        self.assertEqual(str([Phylo.BaseTree.Clade(branch_length=200.0, name='population_A')]), str(tree_output[0]['child_clades']))
        tree_output[0].pop("child_clades")
        
        self.assertEqual([{'pop_name': 'p1', 'population_size': 100, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 
                'split_ratio': 0.5, 'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False, 
                'backup': False, 'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1000,
                'fitness_profile_nums': [50, 6, 12, 35, 1, 5, 30, 42, 6, 50], 'scaling_value': None, 
                'mutation_rate': 2.5e-06, 'parent_pop_name': None, 'clade_name': 'unnamed_population_p1', 'dist_from_start': 0, 
                'pop_end': 1000, 'terminal_clade': False, 'last_child_clade': False}, {'pop_name': 'p2', 'child_clades': [], 
                'population_size': 100, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 'split_ratio': 0.5, 
                'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False, 
                'backup': False, 'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1200.0,
                'fitness_profile_nums': [50, 6, 12, 35, 1, 5, 30, 42, 6, 50], 'scaling_value': None, 'backup': False,
                'mutation_rate': 2.5e-06, 'parent_pop_name': 'p1', 'clade_name': 'population_A', 'dist_from_start': 1000,
                'pop_end': 1200.0, 'terminal_clade': True, 'last_child_clade': True}], tree_output)
                
                
        #Check a tree with a changed value - not hpc
        self.read_clades.start_params["tree_data_file"] = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_abbr.yaml"])
        tree_output = copy.deepcopy(self.read_clades.read_input_tree())
        self.assertEqual(str([Phylo.BaseTree.Clade(branch_length=200.0, name='population_A')]), str(tree_output[0]['child_clades']))
        tree_output[0].pop("child_clades")
        self.assertEqual([{'pop_name': 'p1', 'population_size': 100, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 
                'split_ratio': 0.5, 'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False,
                'fitness_profile_nums': [50, 6, 12, 35, 1, 5, 30, 42, 6, 50], 'scaling_value': None, 
                'backup': False, 'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1000, 
                'mutation_rate': 2.5e-06, 'parent_pop_name': None,  'clade_name': 'unnamed_population_p1','dist_from_start': 0, 
                'pop_end': 1000, 'terminal_clade': False, 'last_child_clade': False}, {'pop_name': 'p2', 'child_clades': [], 
                'population_size': 50, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 'split_ratio': 0.5, 
                'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False, 
                'fitness_profile_nums': [50, 6, 12, 35, 1, 5, 30, 42, 6, 50], 'scaling_value': None, 'backup': False, 
                'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1200.0, 
                'mutation_rate': 2.5e-06, 'parent_pop_name': 'p1',  'clade_name': 'population_A', 'dist_from_start': 1000, 
                'pop_end': 1200.0, 'terminal_clade': True, 'last_child_clade': True}], tree_output)
        
        #Check a tree that is not in nexus format
        self.read_clades.start_params["input_tree"] = self.test_file_path + "nexus_tree.txt"
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_input_tree()
        self.assertEqual(cm.exception.code, 1)
          
        #Check a bad tree
        self.read_clades.start_params["input_tree"] = self.test_file_path + "bad_tree.txt"
        with self.assertRaises(SystemExit) as cm:
            with redirect_stderr(io.StringIO()) as serr:
                self.read_clades.read_input_tree()
        self.assertEqual(cm.exception.code, 1)
        
        self.assertEqual(cm.exception.code, 1)
        self.assertEqual(serr.getvalue(), ("Please make sure your input tree is in Newick format. Program closing\n"))
        serr.close()
        self.read_clades.start_params["input_tree"] = self.test_file_path + "test_tree.txt"


        #Check a tree with all possible changed values mu, not mutation matrix - hpc
        with open(self.test_file_path + "correct_parameters_hpc.yaml", 'r') as yaml_file:
            starting_params = yaml.safe_load(yaml_file)
            starting_params["input_tree"] = self.input_tree
        yaml_file.close()
        self.read_clades.start_params = starting_params
        self.read_clades.start_params["tree_data_file"] = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_abbr_hpc.yaml"])
        
        tree_output = copy.deepcopy(self.read_clades.read_input_tree())
        self.assertEqual(str([Phylo.BaseTree.Clade(branch_length=200.0, name='population_A')]), str(tree_output[0]['child_clades']))
        tree_output[0].pop("child_clades")
        self.assertEqual([{'pop_name': 'p1', 'population_size': 100, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 
                'split_ratio': 0.5, 'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False,
                'fitness_profile_nums': [50, 6, 12, 35, 1, 5, 30, 42, 6, 50], 'scaling_value': None, 
                'backup': False, 'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1000, 
                'mutation_rate': 2.5e-06, 'parent_pop_name': None,  'clade_name': 'unnamed_population_p1','dist_from_start': 0,
                'pop_end': 1000, 'terminal_clade': False, 'last_child_clade': False}, {'pop_name': 'p2', 'child_clades': [], 
                'partition': 'anothernode', 'population_size': 50, 'recombination_rate': 0.05, 'sample_size': '30', 'split_ratio': 0.4, 
                'time': '7200', 'count_subs': False, 'output_gens': False, 
                'fitness_profile_nums': [50, 6, 12, 35, 1, 5, 30, 42, 6, 50], 'scaling_value': None, 'backup': False, 
                'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1200.0, 
                'mutation_rate': 0.0025, 'parent_pop_name': 'p1', 'clade_name': 'population_A', 'dist_from_start': 1000, 'pop_end': 1200.0, 
                'terminal_clade': True, 'last_child_clade': True}], tree_output)
        
        self.maxDiff = None
        
        
        #Check a tree with all possible changed values mutation matrix, not mu - hpc
        self.read_clades.start_params["tree_data_file"] = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_abbr_hpc_mut_mat.yaml"])
        tree_output = copy.deepcopy(self.read_clades.read_input_tree())
        self.assertEqual(str([Phylo.BaseTree.Clade(branch_length=200.0, name='population_A')]), str(tree_output[0]['child_clades']))
        tree_output[0].pop("child_clades")
        mut_mat = tree_output[1]['mutation_matrix']
        self.assertEqual([0.0e+00, 3.5e-07, 3.5e-07, 3.5e-07,3.5e-07, 0.0e+00, 3.5e-07, 3.5e-07,
                        3.5e-07, 3.5e-07, 0.0e+00, 3.5e-07, 3.5e-07, 3.5e-07, 3.5e-07, 0.0e+00], mut_mat[0].flatten().tolist())
        self.assertEqual('matrix(c(0.0e+00, 3.5e-07, 3.5e-07, 3.5e-07, 3.5e-07, 0.0e+00, ' +
                        '3.5e-07, 3.5e-07,\n 3.5e-07, 3.5e-07, 0.0e+00, 3.5e-07, ' +
                        '3.5e-07, 3.5e-07, 3.5e-07, 0.0e+00), ncol = 4, byrow = T)', mut_mat[1])
        tree_output[1].pop("mutation_matrix")
        self.assertEqual([{'pop_name': 'p1', 'population_size': 100, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 
                'split_ratio': 0.5, 'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False,
                'backup': False, 'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1000,
                'fitness_profile_nums': [50, 6, 12, 35, 1, 5, 30, 42, 6, 50], 'scaling_value': None,  
                'mutation_rate': 2.5e-06, 'parent_pop_name': None, 'clade_name': 'unnamed_population_p1', 'dist_from_start': 0, 'pop_end': 1000, 
                'terminal_clade': False, 'last_child_clade': False}, {'pop_name': 'p2', 'child_clades': [], 
                'partition': 'anothernode', 'population_size': 50, 'recombination_rate': 0.05, 'sample_size': '30', 'split_ratio': 0.4, 
                'time': '7200', 'count_subs': False, 'output_gens': False, 
                'fitness_profile_nums': [50, 6, 12, 35, 1, 5, 30, 42, 6, 50], 'scaling_value': None, 'backup': False, 
                'polymorphisms': False, 'jukes_cantor': False, 'calculate_selection': False, 'end_dist': 1200.0,
                'parent_pop_name': 'p1','clade_name': 'population_A', 'dist_from_start': 1000, 'pop_end': 1200.0,
                'terminal_clade': True, 'last_child_clade': True}], tree_output)

        
