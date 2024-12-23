#Script to parse through input and find data for each clade in a tree

from Bio import Phylo
import copy
import math
import yaml
import sys
from utils import readInput
import numpy as np


class cladeReader:
    
    def __init__(self, start_params):     
        self.start_params = start_params
        
        if (self.start_params["tree_data_file"] != None):
            self.start_params["tree_data_file"] = self.read_clade_data(self.start_params["tree_data_file"])
        
        self.clade_dict_list = self.read_input_tree()



    #Read individual data from the file - given in YAML formatting
    def read_clade_data(self, data_file):
        data_file = data_file[0]
        
        #Read yaml data file of changes in tree
        with open(data_file, 'r') as yaml_file:
            yaml_data = yaml.safe_load(yaml_file)
            
            #Check to make sure yaml file is loaded
            if (not isinstance(yaml_data, dict)):
                print("Please make sure your changes are in yaml format. " +
                        "For more information on yaml format visit: " +
                        "https://en.wikipedia.org/wiki/YAML. Program closing.")
                sys.exit(0)
            
        yaml_file.close()
        
        #Dictionary to change key names for abbreviations
        data_translation = {
            'p': 'partition',
            't': 'time',
            'n': 'population_size',
            'r': 'recombination_rate',
            'v': 'mutation_rate',
            'm': 'mutation_matrix',
            'k': 'sample_size',
            'sr': 'split_ratio'
        }
        
        possible_hpc_changes = list(data_translation.keys()) + list(data_translation.values())
        possible_changes = ['n', 'population_size']
        
        #Recurse through changed branches
        for dat in yaml_data:
            #Translate keynames to full names if abbreviation used
            try:
                new_dict = {(data_translation[k] if k in data_translation else k):v  for (k,v) in yaml_data[dat].items() }
            except AttributeError:
                print("Please check the formatting of your yaml file, you likely did not include the branch name. Closing program.")
                sys.exit(0)
            
            #Check to make sure no changes are specified that cannot be handled
            if self.start_params["high_performance_computing"] and len(np.setdiff1d(list(new_dict.keys()), possible_hpc_changes)) != 0:
                print("When using slim-tree with HPC only the following parameters may be modified for specific branches:")
                print(*list(data_translation.values()), sep = "\n")
                sys.exit(0)
            elif not self.start_params["high_performance_computing"] and len(np.setdiff1d(list(new_dict.keys()), possible_changes)) != 0:
                print("When using slim-tree without HPC, only the population size may be modified for specific branches. Exiting")
                sys.exit(0)
                
            #If mutation rate or mutation matrix changed, specify Jukes-Cantor accordingly
            if('mutation_rate' in new_dict.keys()):
                new_dict['jukes_cantor'] = True
            elif('mutation_matrix' in new_dict.keys()):
                new_dict['jukes_cantor'] = False
                
            
            
            yaml_data[dat] = new_dict

            
        return(yaml_data)



    #Read the phylogenetic tree data given by the user 
    def read_input_tree(self):

        self.pop_num = 0
        
        try:
            phylogeny = Phylo.read(self.start_params["input_tree"],"newick")
        except ValueError:
            print ("Please make sure your input tree is in Newick format. Program closing")
            sys.exit(0)

        #Figure out how long the simulation is going to run for
        max_depth = int(max(phylogeny.depths().values())) + self.start_params["burn_in"] + 1
        self.start_params["num_generations"] = max_depth

        #Set up starting parameters and make list of dictionaries of variables
        starting_parameter_dict = {
            "pop_name": None,
            "child_clades" : None,
            "population_size" : self.start_params["population_size"],
            "recombination_rate" : self.start_params["recombination_rate"],
            "sample_size": self.start_params["sample_size"],
            "split_ratio": self.start_params["split_ratio"],
            "partition": self.start_params["partition"],
            "time" : self.start_params["time"],
            "count_subs" : self.start_params["count_subs"],
            "output_gens" : self.start_params["output_gens"],
            "backup" : self.start_params["backup"],
            "polymorphisms": self.start_params["polymorphisms"],
            "jukes_cantor": self.start_params["jukes_cantor"],
            "calculate_selection": self.start_params["calculate_selection"],
            "end_dist" : 0,
        }

        if(self.start_params["jukes_cantor"]):
            starting_parameter_dict["mutation_rate"] = self.start_params["mutation_rate"]
        else:
            starting_parameter_dict["mutation_matrix"] = self.start_params["mutation_matrix"]

        try:
            clade_dict_list = self.recurse_through_clades(phylogeny.get_nonterminals()[0],
                                                 starting_parameter_dict, self.start_params["tree_data_file"], phylogeny)
        except IndexError:
            print ("Please make sure your input tree is in Newick format. Program closing")
            sys.exit(0)

        #Sort clade dict list by the distance from the start of the simulation so that it works properly
        #in SLiM

        clade_dict_list = sorted(clade_dict_list, key=lambda k: k["dist_from_start"])


        return (clade_dict_list)




    #Recurses through each of the clades in the phylogeny and makes a list of dictionaries of their parameters
    #Recursion is depth first - list is made so that it can be traversed bredth first
    def recurse_through_clades(self, current_clade, parent_clade_dict, clade_data, phylogeny):
        clade_dict_list = self.get_clade_data(current_clade,parent_clade_dict,clade_data, phylogeny)
        clade_dict = clade_dict_list[0]

        #Make list of clades from data
        if (len(clade_dict["child_clades"])==0):
            return (clade_dict_list)
        else:
            #Recurse through all child clades
            list_of_child_clades = []
            for child_clade in clade_dict["child_clades"]:
                child_clade_dict = self.recurse_through_clades(child_clade, clade_dict,clade_data, phylogeny)
                list_of_child_clades = list_of_child_clades + child_clade_dict

            return clade_dict_list + list_of_child_clades




    #Set up data such as clade name, mutation rate, population size, etc. for a clade
    def get_clade_data (self, clade, parent_clade_dict, clade_data, phylogeny):
        #Set up the default parameters based on the parent dictionary
        clade_dict = copy.deepcopy(parent_clade_dict)


        #Change parameters if specified by user for a specific clade
        if(clade_data != None):
            clade_name = clade.name
            if(clade_name in clade_data.keys()):
                current_clade_data = self.start_params["tree_data_file"][clade_name]
                for keyname in current_clade_data.keys():
                    if(keyname == 'mutation_matrix'):
                        input_reader = readInput.readInput()
                        clade_dict[keyname] = input_reader.make_mutation_matrix(str(current_clade_data[keyname]))
                    elif (keyname == 'partition' or keyname == 'time' or keyname == 'sample_size'):
                        clade_dict[keyname] = str(current_clade_data[keyname])
                    else:
                        clade_dict[keyname] = current_clade_data[keyname]


        #Figure out what population name is for self and assign clade name appropriately
        self.pop_num += 1
        pop_name = "p" + str(self.pop_num)

        if(clade.name != None):
            clade.name = pop_name + ": "+ clade.name
        else:
            clade.name = pop_name

        #Figure out when the population needs to be formed
        dist_from_start = parent_clade_dict["end_dist"]
        
        #Determine when sim of pop ends - note must run for at least 1 gen
        pop_end = self.start_params["burn_in"]  + phylogeny.distance(clade)
        
        #Make sure that given tree is in generations
        if (math.floor(dist_from_start) >= math.floor(pop_end)):
            print("Please make sure that your tree is in generations not mutations. Exiting.")
            sys.exit(0)

        #Determine whether population belongs to the last child clade - allows for removal of extraneous data
        parents_children = parent_clade_dict["child_clades"]

        if(parents_children == None):
            last_child_clade = False
        else:
            last_child_clade = clade == parents_children[-1]

        #Update some values of clade
        clade_dict['parent_pop_name'] = clade_dict['pop_name']
        clade_dict['pop_name'] = pop_name 
        clade_dict['child_clades'] = clade.clades
        clade_dict['dist_from_start'] = dist_from_start
        clade_dict['pop_end'] = pop_end
        clade_dict['terminal_clade'] = clade.clades == []
        clade_dict['last_child_clade'] = last_child_clade
        clade_dict['end_dist'] = pop_end
        
        #Remove mutation rate or mutation matrix if not in this clade
        if(clade_dict['jukes_cantor'] and 'mutation_matrix' in clade_dict.keys()):
            clade_dict.pop('mutation_matrix')
        
        if(not clade_dict['jukes_cantor'] and 'mutation_rate' in clade_dict.keys()):
            clade_dict.pop('mutation_rate')

        return [clade_dict]
        
        
    #Return clade dict list to public
    def get_clade_dict_list(self):
        return self.clade_dict_list
