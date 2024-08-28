import unittest
import os, io, yaml, copy
from unittest import mock
from utils import cladeReader
from contextlib import redirect_stdout
from Bio import Phylo
import numpy as np

class testCladeReader(unittest.TestCase):

    def setUp(self):
        self.test_file_path = os.path.dirname(os.path.realpath(__file__))  + "/testFiles/"
        
        with open(self.test_file_path + "correct_parameters.yaml", 'r') as yaml_file:
            starting_params = yaml.safe_load(yaml_file)
        yaml_file.close()
            
        self.read_clades = cladeReader.cladeReader(starting_params)
        

    def test_read_clade_data(self):
 
        #Test failure if non-yaml data file is used
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml.yaml"])
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("Please make sure your changes are in yaml format. " +
                        "For more information on yaml format visit: " +
                        "https://en.wikipedia.org/wiki/YAML. Program closing.\n"))
        sout.close()
        
        #Test successful yaml file with abbreviated names
        yaml_output = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_abbr.yaml"])
        self.assertEqual({"A":{"population_size" : 50}}, yaml_output)
        
        #Test successful yaml file with full names
        yaml_output = self.read_clades.read_clade_data([self.test_file_path + "good_yaml.yaml"])
        self.assertEqual({"A":{"population_size" : 50}}, yaml_output)
        
        #Test failure of yaml when no branch specified
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_no_branch.yaml"])
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("Please check the formatting of your yaml file, you likely did not include the branch name. Closing program.\n"))
        sout.close()
        
        #Test failure with yaml file with changes that cannot be made
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_wrong_vals.yaml"])
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("When using slim-tree without HPC, only the population size may be modified for specific branches. Exiting\n"))
        sout.close()
        
        #Test successful hpc yaml file
        self.read_clades.start_params["high_performance_computing"] = True
        yaml_output = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_abbr_hpc.yaml"])
        self.assertEqual({"A":{"population_size" : 50, "mutation_rate" : 0.0025, "jukes_cantor" : True, "partition": "anothernode", "time":7200, "recombination_rate": 0.05, "sample_size": 30, "split_ratio":0.4}}, yaml_output)
        
        
        #Test hpc yaml file with changes that cannot be made
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.read_clades.read_clade_data([self.test_file_path + "bad_yaml_abbr_hpc.yaml"])
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("When using slim-tree with HPC only the following parameters may be modified for specific branches:\n" + 
                        "\n".join(["partition", "time", "population_size", "recombination_rate", "mutation_rate", "mutation_matrix", "sample_size", "split_ratio\n"] )))
        sout.close()
        self.read_clades.start_params["high_performance_computing"] = False
      
        
        
    def test_read_input_trees(self):
    
        #If something fails in this test, can also be a failure of recurse through clades or get_clade_data
    
        #Check a correct tree
        tree_output = copy.deepcopy(self.read_clades.read_input_tree())
        self.assertEqual(str([Phylo.BaseTree.Clade(branch_length=200.0, name='p2: A')]), str(tree_output[0]['child_clades']))
        tree_output[0].pop("child_clades")
        self.assertEqual([{'pop_name': 'p1', 'population_size': 100, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 
                'split_ratio': 0.5, 'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False,
                'backup': False, 'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1000, 
                'mutation_rate': 2.5e-06, 'parent_pop_name': None, 'dist_from_start': 0, 'pop_end': 1000, 
                'terminal_clade': False, 'last_child_clade': False}, {'pop_name': 'p2', 'child_clades': [], 
                'population_size': 100, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 'split_ratio': 0.5, 
                'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False, 'backup': False, 
                'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1200.0, 
                'mutation_rate': 2.5e-06, 'parent_pop_name': 'p1', 'dist_from_start': 1000, 'pop_end': 1200.0, 
                'terminal_clade': True, 'last_child_clade': True}], tree_output)
        
        #Check a tree with a changed value - not hpc
        self.read_clades.start_params["tree_data_file"] = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_abbr.yaml"])
        tree_output = copy.deepcopy(self.read_clades.read_input_tree())
        self.assertEqual(str([Phylo.BaseTree.Clade(branch_length=200.0, name='p2: A')]), str(tree_output[0]['child_clades']))
        tree_output[0].pop("child_clades")
        self.assertEqual([{'pop_name': 'p1', 'population_size': 100, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 
                'split_ratio': 0.5, 'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False,
                'backup': False, 'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1000, 
                'mutation_rate': 2.5e-06, 'parent_pop_name': None, 'dist_from_start': 0, 'pop_end': 1000, 
                'terminal_clade': False, 'last_child_clade': False}, {'pop_name': 'p2', 'child_clades': [], 
                'population_size': 50, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 'split_ratio': 0.5, 
                'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False, 'backup': False, 
                'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1200.0, 
                'mutation_rate': 2.5e-06, 'parent_pop_name': 'p1', 'dist_from_start': 1000, 'pop_end': 1200.0, 
                'terminal_clade': True, 'last_child_clade': True}], tree_output)
        
        #Check a tree that is not in nexus format
        self.read_clades.start_params["input_tree"] = self.test_file_path + "nexus_tree.txt"
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.read_clades.read_input_tree()
        self.assertEqual(cm.exception.code, 0)
          
        #Check a bad tree
        self.read_clades.start_params["input_tree"] = self.test_file_path + "bad_tree.txt"
        with self.assertRaises(SystemExit) as cm:
            with redirect_stdout(io.StringIO()) as sout:
                self.read_clades.read_input_tree()
        self.assertEqual(cm.exception.code, 0)
        
        self.assertEqual(cm.exception.code, 0)
        self.assertEqual(sout.getvalue(), ("Please make sure your input tree is in Newick format. Program closing\n"))
        sout.close()
        self.read_clades.start_params["input_tree"] = self.test_file_path + "test_tree.txt"


        #Check a tree with all possible changed values mu, not mutation matrix - hpc
        with open(self.test_file_path + "correct_parameters_hpc.yaml", 'r') as yaml_file:
            starting_params = yaml.safe_load(yaml_file)
        yaml_file.close()
        self.read_clades.start_params = starting_params
        self.read_clades.start_params["tree_data_file"] = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_abbr_hpc.yaml"])
        
        tree_output = copy.deepcopy(self.read_clades.read_input_tree())
        self.assertEqual(str([Phylo.BaseTree.Clade(branch_length=200.0, name='p2: A')]), str(tree_output[0]['child_clades']))
        tree_output[0].pop("child_clades")
        self.assertEqual([{'pop_name': 'p1', 'population_size': 100, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 
                'split_ratio': 0.5, 'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False,
                'backup': False, 'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1000, 
                'mutation_rate': 2.5e-06, 'parent_pop_name': None, 'dist_from_start': 0, 'pop_end': 1000, 
                'terminal_clade': False, 'last_child_clade': False}, {'pop_name': 'p2', 'child_clades': [], 
                'partition': 'anothernode', 'population_size': 50, 'recombination_rate': 0.05, 'sample_size': '30', 'split_ratio': 0.4, 
                'time': '7200', 'count_subs': False, 'output_gens': False, 'backup': False, 
                'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1200.0, 
                'mutation_rate': 0.0025, 'parent_pop_name': 'p1', 'dist_from_start': 1000, 'pop_end': 1200.0, 
                'terminal_clade': True, 'last_child_clade': True}], tree_output)
        
        self.maxDiff = None
        #Check a tree with all possible changed values mutation matrix, not mu - hpc
        self.read_clades.start_params["tree_data_file"] = self.read_clades.read_clade_data([self.test_file_path + "good_yaml_abbr_hpc_mut_mat.yaml"])
        tree_output = copy.deepcopy(self.read_clades.read_input_tree())
        self.assertEqual(str([Phylo.BaseTree.Clade(branch_length=200.0, name='p2: A')]), str(tree_output[0]['child_clades']))
        tree_output[0].pop("child_clades")
        mut_mat = tree_output[1]['mutation_matrix']
        self.assertEqual([0.0e+00, 3.5e-07, 3.5e-07, 3.5e-07,3.5e-07, 0.0e+00, 3.5e-07, 3.5e-07,
                        3.5e-07, 3.5e-07, 0.0e+00, 3.5e-07, 3.5e-07, 3.5e-07, 3.5e-07, 0.0e+00], mut_mat[0].flatten().tolist())
        self.assertEqual('matrix(c(0.0, 3.5e-07, 3.5e-07, 3.5e-07, 3.5e-07, 0.0, ' +
                        '3.5e-07, 3.5e-07, 3.5e-07, 3.5e-07, 0.0, 3.5e-07, ' +
                        '3.5e-07, 3.5e-07, 3.5e-07, 0.0), ncol = 4, byrow = T)', mut_mat[1])
        tree_output[1].pop("mutation_matrix")
        self.assertEqual([{'pop_name': 'p1', 'population_size': 100, 'recombination_rate': 2.5e-08, 'sample_size': 'all', 
                'split_ratio': 0.5, 'partition': 'apophis', 'time': '12:00:00', 'count_subs': False, 'output_gens': False,
                'backup': False, 'polymorphisms': False, 'jukes_cantor': True, 'calculate_selection': False, 'end_dist': 1000, 
                'mutation_rate': 2.5e-06, 'parent_pop_name': None, 'dist_from_start': 0, 'pop_end': 1000, 
                'terminal_clade': False, 'last_child_clade': False}, {'pop_name': 'p2', 'child_clades': [], 
                'partition': 'anothernode', 'population_size': 50, 'recombination_rate': 0.05, 'sample_size': '30', 'split_ratio': 0.4, 
                'time': '7200', 'count_subs': False, 'output_gens': False, 'backup': False, 
                'polymorphisms': False, 'jukes_cantor': False, 'calculate_selection': False, 'end_dist': 1200.0,
                'parent_pop_name': 'p1', 'dist_from_start': 1000, 'pop_end': 1200.0,
                'terminal_clade': True, 'last_child_clade': True}], tree_output)
        
