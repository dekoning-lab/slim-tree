#Program to take in a tree in Newick file and output SLiM code to run a simulation from that tree

#Required packages:
#sys
#argparse
#BioPython
#matplotlib
#random
#pandas
#numpy
#os
#json
#string
#math

import sys, argparse, random, os, json, string, math, pandas
import numpy as np
from Bio import Phylo
from matplotlib.pyplot import show, savefig
from writeSLiM import writeSLiM
from writeSLiMHPC import writeSLiMHPC

from getUserDefinedSequence import getUserDefinedSequence
from writeSLiMProtein import writeSLiMProtein
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import PPBuilder
from ContactMap import ContactMap
from utils import get_structure

class SLiMTree:

    #Main script to run other commands
    def __init__(self):
        self.read_user_input()
        self.find_fitness_profile()

        clade_data = self.read_clade_data()
        if (self.data_file != None):
            self.data_file.close()

        clade_dict_list = self.read_input_tree(clade_data)

        self.write_slim_code(clade_dict_list)



    #Read input parameters from the user
    def read_user_input(self):

        #Set up starting parameters dictionary
        self.starting_parameters = {}

        #Parse for arguments given by the user
        parser = argparse.ArgumentParser(description='A program to make slim simulations from newick phylogeny files')
        parser.add_argument('-i','--input_tree', nargs = 1, required = True, type = str,
                help = 'tree file in newick format specifying the clades to simulate')
        parser.add_argument('-d','--tree_data_file', nargs = 1, type=argparse.FileType('r'),
                help = 'file specifying population size, mutation rate, etc. for each node, see documentation')
        parser.add_argument('-T', '--tool', type = str, required = True,
                help = 'name of tool you would like to use. Options include SLiM-Tree, SLiM-Tree-HPC. Default = SLiM-Tree')
        parser.add_argument('-p', '--partition', type = str, help = 'partition to run SLiM-Tree HPC on')
        parser.add_argument('-t', '--time', type = str,
                help = 'maximum time to run each simulation for - suggested time is the maximum time available for a partition')


        #Default parameters are somewhat arbitrary and should be changed for each sim
        parser.add_argument('-n','--population_size', help = 'starting population size for the simulation, default = 100', type = int, default = 100)
        parser.add_argument('-v','--mutation_rate', help = 'starting mutation rate for the simulation, default = 2.5e-6', type=float, default = 2.5e-6)
        parser.add_argument('-g','--genome_length', help = 'length of the genome - amino acids, default = 500', type=int, default = 500)
        parser.add_argument('-r','--recombination_rate', help = 'recombination rate, default = 2.5e-8', type=float, default = 2.5e-8)
        parser.add_argument('-b','--burn_in_multiplier', help = 'value to multiply population size by for burn in, default = 10', type=float, default = 10)
        parser.add_argument('-k','--sample_size', help = 'size of sample obtained from each population at output. Input all for whole sample and consensus for consensus sequence. default = 10', type=str, default = "10")

        parser.add_argument('-c','--count_subs', type = self.str2bool, default = True, const=True, nargs='?',
                help = 'boolean specifying whether to count substitutions, turning off will speed up sims. default = True')
        parser.add_argument('-o','--output_gens', type = self.str2bool, default = True, const=True, nargs='?',
                help = 'boolean specifying whether to output the generation after every 100th generation. default = True')
        parser.add_argument('-B','--backup', type = self.str2bool, default = True, const=True, nargs='?',
                help = 'boolean specifying whether to backup simulations, turning off will save space. default = True')

        parser.add_argument('-w', '--wright_fisher_model', type = self.str2bool, default=True, const=True, nargs='?',
                help = 'boolean specifying whether this is a wright-fisher model or non-wright-fisher model. default = True')

        parser.add_argument('-G', '--gene_count', type = int, default = 1, help = "Number of genes in the model. Default = 1.")
        parser.add_argument('-C', '--coding_ratio', type = float, default = 1.0, help = "Ratio of the genome which are coding regions as a ratio coding/noncoding. Default = 1.0")

        parser.add_argument('-s', '--user_provided_sequence', type = self.str2bool, default = False, const = True, nargs = '?',
                help = 'boolean specifying whether user provides ancestral sequence and coding regions, Default = False')
        parser.add_argument('-f', '--fasta_file', type = str, default = None, help = 'fasta file containing ancestral sequence - please provide only 1 sequence')
        parser.add_argument('-gb', '--genbank_file', type = str, default = None, help = 'genbank file containing information about ancestral genome - please provide data for only 1 genome')

        parser.add_argument('-fc', '--fitness_profile_calc', type = self.str2bool, default = True, const = True, nargs='?',
                help = 'boolean specifying whether fitness profiles should be used to calculate fitness. If false, protein structure fitness will be calculated, and a contact map must be provided. Default = True.')

        parser.add_argument('-pdb', '--pdb_file', type = str, help = 'Path to file containing a PDB file with a valid protein structure. Either this or a contact map is required for calculating fitness based on protein structure.')
        #Add in PDB functionality later perhaps?

        parser.add_argument('-cm', '--contact_map', type = str, help = "Path to file containing a contact map showing protein contacts. This is required for calculating fitness based on protein structure.")

        #Get arguments from user
        arguments = parser.parse_args()

        #Set up tree
        self.input_file = arguments.input_tree[0]

        #Get simulation type and ensure that required arguments are given for simulation type
        self.simulationType = arguments.tool.translate(str.maketrans('', '', string.punctuation)).lower()
        if (self.simulationType == "slimtreehpc" and (arguments.partition == None or arguments.time == None)):
            print("When using SLiM-Tree-HPC, partition and time data must be provided. Closing program.")
            sys.exit(0)

        #Check to make sure gene count and coding ratio are valid
        if (arguments.gene_count < 0 or arguments.gene_count > arguments.genome_length):
            print("Number of genes must be greater than 0 and less than the length of the genome. Closing program.")
            sys.exit(0);

        if (arguments.coding_ratio < 0 or arguments.coding_ratio > 1.0):
            print("Coding ratio must be greater than 0 and less than or equal to 1. Please re-enter as a ratio (0, 1]. Closing program.")
            sys.exit(0);


        #Ensure that if users specify a user provided sequence they include fasta file and gb file
        if (arguments.user_provided_sequence and (arguments.fasta_file == None or arguments.genbank_file == None )):
            print("When specifying an ancestral sequence, a fasta file containing the sequence and a genbank file containing coding regions must " +
                    "be provided. Closing program.")
            sys.exit(0);


        #If using protein-based fitness, make sure a contact map is provided.

        if (arguments.fitness_profile_calc == False and arguments.contact_map == None and arguments.pdb_file == None):
            print("When calculating protein-based fitness, a contact map or PDB file must be provided. Closing program.")
            sys.exit(0)



        #Set up the starting parameters
        self.starting_parameters["mutation_rate"] = arguments.mutation_rate
        self.starting_parameters["population_size"] = arguments.population_size
        self.starting_parameters["recombination_rate"] = arguments.recombination_rate
        self.starting_parameters["burn_in"] = arguments.burn_in_multiplier * arguments.population_size
        self.starting_parameters["sample_size"] = arguments.sample_size

        self.starting_parameters["fitness_profile_calc"] = arguments.fitness_profile_calc

        self.starting_parameters["contact_map"] = arguments.contact_map

        if (arguments.contact_map != None):
            self.starting_parameters["contact_map"] = np.loadtxt(open(arguments.contact_map, "rb"), delimiter=",")
            #Overwrite the length of the genome with the length of the contact map.
            self.starting_parameters["genome_length"] = len(self.starting_parameters["contact_map"]) + 1
        else:
            self.starting_parameters["contact_map"] = None

        if (arguments.pdb_file != None):
            # structure = get_structure(arguments.pdb_file)
            # model = structure[0]

            pdbstructure = PDBParser().get_structure("pdb", arguments.pdb_file) #Placeholder - replace with file later
            model = pdbstructure[0]
            ppb = PPBuilder()
            aa_seq = ""
            for pp in ppb.build_peptides(pdbstructure):
                aa_seq += pp.get_sequence()

            #aa_seq = aa_seq.tomutable()
            self.starting_parameters["aa_seq"] = aa_seq
            self.starting_parameters["genome_length"] = len(aa_seq) + 1
            if (arguments.contact_map == None):
                map = ContactMap(model, threshold = 7.0)
                self.starting_parameters["contact_map"] = map.get_data().tolist() #Generates a contact map.
        else:
            self.starting_parameters["aa_seq"] = None

        self.starting_parameters["wf_model"] = arguments.wright_fisher_model

        self.starting_parameters["gene_count"] = arguments.gene_count
        self.starting_parameters["coding_ratio"] = arguments.coding_ratio
        self.starting_parameters["coding_seqs"] = self.get_coding_seqs()


        self.starting_parameters["partition"] = arguments.partition
        self.starting_parameters["time"] = arguments.time

        self.starting_parameters["count_subs"] = arguments.count_subs
        self.starting_parameters["output_gens"] = arguments.output_gens
        self.starting_parameters["backup"] = arguments.backup


        self.starting_parameters["user_provided_sequence"] = arguments.user_provided_sequence
        self.starting_parameters["fasta_file"] = arguments.fasta_file
        self.starting_parameters["genbank_file"] = arguments.genbank_file

        self.starting_parameters["wf_model"] = arguments.wright_fisher_model

        #Set up coding sequences if no user defined sequence is specified
        if (not arguments.user_provided_sequence and arguments.pdb_file == None and arguments.contact_map == None):
            self.starting_parameters["genome_length"] = int(arguments.genome_length)
            self.starting_parameters["gene_count"] = arguments.gene_count
            self.starting_parameters["coding_ratio"] = arguments.coding_ratio
            self.starting_parameters["coding_seqs"] = self.get_coding_seqs()






        #Set up the filenames for file io
        input_file_start = os.getcwd() + '/' + self.input_file.split('.')[0]
        self.starting_parameters["tree_filename"] = input_file_start + "_phylogeny.png"
        self.starting_parameters["fasta_filename"] = input_file_start

        if(arguments.tree_data_file != None):
            self.data_file = arguments.tree_data_file[0]
        else:
            self.data_file = None




        #Set up the output of scripts to a single folder
        split_starting_file = input_file_start.split('/')
        output_files_directory = "/".join(split_starting_file[0:(len(split_starting_file)-1)]) + "/slimScripts"
        backup_files_directory = "/".join(split_starting_file[0:(len(split_starting_file)-1)]) + "/backupFiles"


        try:
            os.mkdir(output_files_directory)
            os.mkdir(backup_files_directory)
        except OSError:
            print ("The directory %s already exits, program files will be overwritten" % output_files_directory)

        self.starting_parameters["output_file"] = output_files_directory + "/" + split_starting_file[-1]



        #Save starting parameters and value of theta to a file for later reference
        theta = 4*arguments.mutation_rate*arguments.population_size

        parameter_file = open(input_file_start + "_parameters.txt", 'w')
        parameter_file.write("Simulation parameters\n\n")

        for key, value in self.starting_parameters.items():
            #Don't need to record these filenames as they are not yet complete
            if(key in ['fasta_filename', 'tree_filename', 'output_file']):
                continue

            parameter_file.write('%s:%s\n' % (key, value))

        parameter_file.write("theta: " + str(theta))
        parameter_file.close()


    #Command to take input from user and convert to bool
    #From: https://stackoverflow.com/questions/15008758/parsing-boolean-values-with-argparse
    def str2bool(self, v):
        if isinstance(v, bool):
            return v
        if v.lower() in ('yes', 'true', 't', 'y', '1'):
            return True
        elif v.lower() in ('no', 'false', 'f', 'n', '0'):
            return False
        else:
            raise argparse.ArgumentTypeError('Boolean value expected.')



    #Return the range of coding sequences for the set number of genes.
    def get_coding_seqs(self):

        genome_length = self.starting_parameters["genome_length"]
        coding_ratio = self.starting_parameters["coding_ratio"]
        gene_count = self.starting_parameters["gene_count"]

        if(gene_count == 0 or coding_ratio == 0):
            return None

        percent_coding = math.ceil(int(genome_length) * coding_ratio) #Gives approximate number of coding amino acids
        avg_coding_length = math.ceil(percent_coding / gene_count) #Gives avg length of coding region
        avg_noncoding_length = 0
        if (gene_count != 1):
            avg_noncoding_length = math.floor((genome_length - percent_coding) / (gene_count - 1)) #Average length of non-coding regions by subtracting number of coding aa from total aa

        coding_regions = []
        current_aa = 0

        for i in range (gene_count):
            coding_regions.append(current_aa)
            coding_regions.append(min(current_aa + avg_coding_length, genome_length - 1)) #Ensures that you are not surpasssing the length of the genome
            current_aa = current_aa + avg_noncoding_length + avg_coding_length  #Accounts for the non-coding region + coding region added previously

            #Make sure that the genome will not be longer than the genome
            if (current_aa + avg_coding_length > genome_length - 1):
                current_aa = genome_length - avg_coding_length - 1

        coding_regions = np.stack(np.array_split(coding_regions, gene_count))
        return coding_regions


    #Read fitness profile and stationary distribution data from psi_c50 file, make fitness profiles
    def find_fitness_profile(self):
        #Open stationary and fitness effects csv
        fitness_dist = pandas.io.parsers.read_csv (os.path.dirname(os.path.realpath(__file__)) + '/fitnessDataFiles/table_fitness_profiles.csv')
        stationary_distributions = pandas.io.parsers.read_csv (os.path.dirname(os.path.realpath(__file__)) + '/fitnessDataFiles/table_stationary_distributions.csv')

        #Add an extra distribution for non-coding neutral alleles - neutral, 20 AA + stop codon
        stationary_distributions["neutral"] = 1
        fitness_dist["neutral"] = 1
        fitness_distributions = fitness_dist.values.tolist()


        #Find fitness effects for each amino acid
        amino_acids = (["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N"] +
                               ["P", "Q", "R", "S", "T", "V", "W", "Y", "X"])

        stationary_distributions.rename(index = {0 : 'A', 1 : 'C', 2 : 'D', 3 : 'E', 4 : 'F', 5 : 'G',
                    6 : 'H', 7 : 'I', 8 : 'K', 9 : 'L', 10 : 'M', 11 : 'N', 12 : 'P', 13 : 'Q',
                    14 : 'R', 15 : 'S', 16 : 'T', 17 : 'V', 18 : 'W', 19 : 'Y', 20 : 'X'}, inplace = True)


        fitness_profiles = {}

        for amino_acid_num in range(len(amino_acids)):
                aa = amino_acids[amino_acid_num]
                fitness_profiles[aa] = fitness_distributions[amino_acid_num]


        #Set up fitness profiles
        if (self.starting_parameters["user_provided_sequence"]): #Use user provided fitness profiles if given
            get_seq = getUserDefinedSequence(self.starting_parameters["genbank_file"],
                        self.starting_parameters["fasta_file"], stationary_distributions, fitness_profiles)

            coding_feats = get_seq.get_coding_features()
            ans_seq = get_seq.get_ancestral_sequence()
            fitness_profile_nums = get_seq.find_fitness_profiles(ans_seq, coding_feats)

            self.starting_parameters["coding_seqs"] = coding_feats
            self.starting_parameters["ancestral_sequence"] = ans_seq
            self.starting_parameters["genome_length"] = len(ans_seq)

        else: #Randomly set up fitness profiles if not provided by user
            fitness_profile_nums = []
            fitness_length = fitness_dist.shape[1] - 1
            coding_poses = self.starting_parameters["coding_seqs"]
            for coding_pos in range(len(coding_poses)):
                fitness_profile_nums = (fitness_profile_nums + [fitness_length] + #starting codon is correct
                            random.choices(range(fitness_length),k=coding_poses[coding_pos,1] - coding_poses[coding_pos,0] - 1))

                if (coding_pos != len(coding_poses) - 1):
                    fitness_profile_nums = fitness_profile_nums + list(np.repeat(fitness_length, coding_poses[coding_pos+1,0] - coding_poses[coding_pos,1] ))
                else:
                    fitness_profile_nums = fitness_profile_nums + list(np.repeat(fitness_length, self.starting_parameters["genome_length"] - coding_poses[coding_pos,1]))



        #Set up distributions in starting parameters
        self.starting_parameters["fitness_profile_nums"] = fitness_profile_nums
        self.starting_parameters["min_fitness"] = min(np.array(fitness_distributions).flatten())
        self.starting_parameters["stationary_distributions"] = stationary_distributions
        self.starting_parameters["fitness_profiles"] = fitness_profiles
        self.starting_parameters["amino_acids"] = amino_acids

        #Find scaling for non-wright-fisher models
        if(self.starting_parameters["wf_model"] == False):
            if(self.starting_parameters["user_provided_sequence"]): #Different fitness scaling for user defined sequence because expected value based on AA rather then seq
                self.starting_parameters["scaling_value"] = get_seq.get_fitness_scaling()
            else:
                self.find_fitness_scaling(fitness_distributions)



    #Calculates expected fitness from fitness profiles for scaling in non-wright-fisher models
    def find_fitness_scaling(self, fitness_profiles):
        fitness_profiles = np.transpose(fitness_profiles)
        stationary_dists = self.starting_parameters["stationary_distributions"]
        expected_fitness_profiles = []

        #Find the expected value for each fitness profile
        for num_dist in range(len(stationary_dists.columns)-1):
            mean = 0
            stationary = stationary_dists.iloc[:,num_dist]
            fitness = fitness_profiles[num_dist]
            for num_profile in range(len(fitness)-1):
                mean += stationary[num_profile] * fitness[num_profile]

            expected_fitness_profiles.append(mean)

        #Add fitness profile with mean of 1 to account for the neutral areas
        expected_fitness_profiles = expected_fitness_profiles + [1]
        expected_fitnesses = []

        #Find the expected value for each site in the genome based on it's respective fitness profile
        for fitness_profile in self.starting_parameters["fitness_profile_nums"]:
            expected_fitnesses.append(expected_fitness_profiles[fitness_profile])

        #Find the expected value of all sites by multiplying expected values - squared because there are 2 xsomes
        self.starting_parameters["scaling_value"] = np.prod(expected_fitnesses)**2



    #Read individaul data from the file - to add more parameters modify data_translation_dict
    def read_clade_data(self):

        data_translation_dict = {
            'v': 'mutation_rate',
            'n': 'population_size',
            'r': 'recombination_rate',
            'k': 'sample_size',
            'p': 'partition',
            't': 'time',
            'c': 'count_subs',
            'o': 'output_gens',
            'B': 'backup',
            'mutation_rate': 'mutation_rate',
            'population_size': 'population_size',
            'recombination_rate': 'recombination_rate',
            'sample_size': 'sample_size',
            'partition': 'partition',
            'time': 'time',
            'count_subs': 'count_subs',
            'output_gens': 'output_gens',
            'backup': 'backup'

        }

        if(self.data_file == None):
            return (None)
        else:
            data = {}
            line = self.data_file.readline()

            while (line != ''):
                line = line.split('\n')[0]
                if (line == ''):
                    pass
                elif(line[0] == '@'):
                    data_for_node = {}
                    data[line[1:]] = data_for_node
                elif(line[0] == '-'):
                    data_label = line[1]
                    data_value = line.split(' ')[1]

                    data_for_node[data_translation_dict[data_label]] = data_value

                line = self.data_file.readline()
            return(data)



    #Read the phylogenetic tree data given by the user
    def read_input_tree(self, clade_data):

        self.pop_num = 0

        phylogeny = Phylo.read(self.input_file,"newick")

        #Figure out how long the simulation is going to run for
        max_depth = int(max(phylogeny.depths().values())) + self.starting_parameters["burn_in"] + 1
        self.starting_parameters["num_generations"] = max_depth

        #Set up starting parameters and make list of dictionaries of variables
        starting_parameter_dict = {
            "pop_name": None,
            "child_clades" : None,
            "mutation_rate" : self.starting_parameters["mutation_rate"],
            "population_size" : self.starting_parameters["population_size"],
            "recombination_rate" : self.starting_parameters["recombination_rate"],
            "sample_size": self.starting_parameters["sample_size"],
            "partition": self.starting_parameters["partition"],
            "time" : self.starting_parameters["time"],
            "count_subs" : self.starting_parameters["count_subs"],
            "output_gens" : self.starting_parameters["output_gens"],
            "backup" : self.starting_parameters["backup"]
        }

        try:
            clade_dict_list = self.recurse_through_clades(phylogeny.get_nonterminals()[0],
                                                 starting_parameter_dict, clade_data, phylogeny)
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

        #Set up the default parameters based on the parent
        mut_rate = parent_clade_dict["mutation_rate"]
        pop_size = parent_clade_dict["population_size"]
        rec_rate = parent_clade_dict["recombination_rate"]
        samp_size = parent_clade_dict["sample_size"]
        part = parent_clade_dict["partition"]
        time = parent_clade_dict["time"]
        subs = parent_clade_dict["count_subs"]
        gens = parent_clade_dict["output_gens"]
        backup = parent_clade_dict["backup"]

        #Change parameters if specified by user
        if(clade_data != None):
            clade_name = clade.name
            if(clade_name in clade_data.keys()):
                current_clade_data = clade_data[clade_name]

                if('mutation_rate' in current_clade_data.keys()):
                    mut_rate = float(current_clade_data['mutation_rate'])
                if('population_size' in current_clade_data.keys()):
                    pop_size = int(current_clade_data['population_size'])
                if('recombination_rate' in current_clade_data.keys()):
                    rec_rate = float(current_clade_data['recombination_rate'])
                if('sample_size' in current_clade_data.keys()):
                    samp_size = int(current_clade_data['sample_size'])
                if('partition' in current_clade_data.keys()):
                    part = current_clade_data['partition']
                if('time' in current_clade_data.keys()):
                    time = current_clade_data['time']
                if('count_subs' in current_clade_data.keys()):
                    subs = self.str2bool(current_clade_data['count_subs'])
                if('output_gens' in current_clade_data.keys()):
                    gens = self.str2bool(current_clade_data['output_gens'])
                if('backup' in current_clade_data.keys()):
                    backup = self.str2bool(current_clade_data['backup'])


        #Figure out what population name is for self and assign clade name appropriately
        self.pop_num += 1
        pop_name = "p" + str(self.pop_num)

        if(clade.name != None):
            clade.name = pop_name + ": "+ clade.name
        else:
            clade.name = pop_name

        #Figure out when the population needs to be formed
        dist_from_parent = clade.branch_length
        if(dist_from_parent == None):
            dist_from_start = 0
        else:
            dist_from_start = parent_clade_dict["end_dist"]


        #Determine whether population belongs to the last child clade - allows for removal of extraneous data
        parents_children = parent_clade_dict["child_clades"]

        if(parents_children == None):
            last_child_clade = False
        else:
            last_child_clade = clade == parents_children[-1]

        #Set up the dictionary of values for the clade
        clade_dict = {
            "pop_name": pop_name,
            "parent_pop_name" : parent_clade_dict["pop_name"],
            "child_clades" : clade.clades,
            "mutation_rate" : mut_rate,
            "population_size" : pop_size,
            "recombination_rate": rec_rate,
            "dist_from_start" : dist_from_start,
            "end_dist": self.starting_parameters["burn_in"]  + phylogeny.distance(clade),
            "terminal_clade" : clade.clades == [],
            "last_child_clade" : last_child_clade,
            "sample_size": samp_size,
            "partition": part,
            "time" : time,
            "count_subs" : subs,
            "output_gens" : gens,
            "backup" : backup
        }

        return [clade_dict]


    #Using the dictionary of the data from the clades, write slim code to represent the phylogeny
    def write_slim_code (self, clade_dict_list):

        #Open SLiM writer based on tool type and write the initialize statement
        if(self.simulationType == "slimtree"):
            if (self.starting_parameters["fitness_profile_calc"]): #If this is using protein-based fitness, use writeSLiMProtein
                SLiM_Writer = writeSLiM(self.starting_parameters)
            else:
                SLiM_Writer = writeSLiMProtein(self.starting_parameters)
        elif(self.simulationType == "slimtreehpc"):
            SLiM_Writer = writeSLiMHPC(self.starting_parameters)
        else:
            print ("Invalid tool type. Please specify a tool as SLiM-Tree or SLiM-Tree-HPC. Program closing")
            sys.exit(0)

        #Write a script for each clade which will be run sequentially
        if (self.starting_parameters["wf_model"]): #If this is a Wright-Fisher model, use a different write_subpop function
            for clade in clade_dict_list:
                SLiM_Writer.write_subpop(clade)
        else:
            for clade in clade_dict_list:
                SLiM_Writer.write_subpop_nonwf(clade);

        #Start the SLiM code to run
        if(self.simulationType == "slimtree"):
            SLiM_Writer.close_file()
            os.system("slim \"" + self.starting_parameters["output_file"] + "_p1.slim\"")
        elif(self.simulationType == "slimtreehpc"):
            os.system("sbatch \"" + self.starting_parameters["output_file"] + "_p1.sh\"")



if __name__=='__main__':
    SLiMTree()
