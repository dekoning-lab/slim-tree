#Script which reads user input for slim-tree and writes to dictionary which is returned to main slim-tree

import argparse, sys, os, yaml



class readInput:
    def __init__(self):
        pass
        
    #Run all commands to get the input from the user
    def process_input(self):
        args = self.read_user_input()
        if(not self.check_arguments(args)):
            sys.exit(0)
            
        self.param_dict = self.make_param_dict(args)
        self.param_dict["filenames"] = self.process_filenames(args.input_tree, args.backup)
        
        


    #Read input parameters from the user
    def read_user_input(self):

        #Set up starting parameters dictionary
        self.starting_parameters = {}

        #Parse for arguments given by the user
        parser = argparse.ArgumentParser(description='Wrapper program to make slim scripts from a newick formatted' +
                                'phylogeny with realistic fitness effects using either fitness profiles or structure-based fitness effects.')
        
        #Main tree that is a required argument
        parser.add_argument('input_tree', type = str, help = 'newick formatted tree file with branch lengths in generations')
        
        #Table of possible stationary distributions - which we convert to fitness effects, can also input fitness effects
        parser.add_argument('codon_stationary_distributions', type = str, help = 'file containing a series of stationary codon ' +
                                'distributions from which to calculate fitness.')
        parser.add_argument('-fd', '--aa_fitness_distributions', type = str, help = 'file containing a amino acid fitnesses')
                
        #High performance computing parameters - allows for computation using Slurm
        parser.add_argument('-hpc','--high_performance_computing', action='store_true', default=False, help = 'boolean flag to turn on ' +
                                'slim-tree high performance computing. Slurm is required')
        parser.add_argument('-p', '--partition', type = str, help = 'partition to run Slurm on - required if using high performance computing')
        parser.add_argument('-t', '--time', type = str, help = 'maximum time to run each simulation for - suggested time is the maximum ' +
                                'time available for a partition - required if using high performance computing')

        #Simulation parameters
        parser.add_argument('-w', '--nonWF', action='store_true', default=False, help = 'boolean flag to specify that a non-wright-fisher ' +
                                'model should be used in lieu of a wright-fisher model.')
        #+'Note: a non-wright-fisher model is used for all simulations using structure based calculations of protein fitness.')
        # parser.add_argument('-sf', '--structure_fitness_effects', action='store_true', default=False,
                # help = 'boolean flag specifying that fitness effects are to be calculated using protein structures rather than fitness profiles.' +
                                # 'A pdb file specifying the ancestral protein structure is required. A non-wright-fisher model will be used')
        
                
        #Arguments for specifying population parameters
        parser.add_argument('-n','--population_size', help = 'starting population size for the simulation, default = 100', type = int, default = 100)
        parser.add_argument('-b','--burn_in_multiplier', type=int, default = 10, help = 'value to multiply population size by ' +
                                    'for burn in, default = 10')
        parser.add_argument('-r','--recombination_rate', help = 'recombination rate, default = 2.5e-8', type=float, default = 2.5e-8)
        parser.add_argument('-v','--mutation_rate', help = 'starting mutation rate for the simulation, default = 2.5e-6', type=float, default = 2.5e-6)
        parser.add_argument('-m', '--mutation_matrix',  type = self.make_mutation_matrix, help = 'CSV file specifying a mutation ' +
                                    'rate matrix, matrix should be either 4 by 4 or 4 by 64 specifying rates from nucleotide to nucleotide and tri-nucleotide to nucleotide respectfully. Nucleotides and tri-nucleotides should be in alphabetical order with no headers. If mutation rate matrix is supplied, mutation rate will be ignored')
        
        
        #Specify population parameters for specific branches
        parser.add_argument('-d','--tree_data_file', nargs = 1, type=argparse.FileType('r'), help = 'file to change the population size ' +
                                    'for specific branches using YAML formatting. When using HPC, other parameters may also be changed.')
                
        
        #Genome parameters
        parser.add_argument('-g','--genome_length', help = 'length of the genome - amino acids, default = 300', type=int, default = 300)
        parser.add_argument('-G', '--gene_count', type = int, default = 1, help = "number of genes to be simulated by the model, default = 1")
        parser.add_argument('-C', '--coding_ratio', type = float, default = 1.0, help = "ratio of the genome which is coding, default = 1.0")
        
        
        #Fitness profile calculations
        parser.add_argument('-f', '--fasta_file', type = str, default = None, help = 'fasta file containing ancestral sequence (amino acids), replaces random creation of ancestral sequence. Fitness profiles for each amino acid are required')
        
        
        #Tree parameters
        parser.add_argument('-k','--sample_size', help = 'size of sample obtained from each population at a  tree tip at the end of the simulations.' +
                                    'Input \'all\' for the every member of the tree tip samples and consensus for the consensus sequence of the ' + 
                                    'population at each tip. default = all', type=str, default = 'all')
        parser.add_argument('-sr', '--split_ratio', help = 'proportion of a population that goes into the first daughter branch at a tree ' +
                                    'branching point in non-wright fisher models. must be ratio between 0 and 1.0. default = 0.5', 
                                    type = float, default = 0.5)


	#Flags to turn on and off simulation functions
        parser.add_argument('-c','--count_subs', action='store_true', default=False, help = 'boolean flag to turn on substitution counting. ' +
                                    'This will slow down simulations')
        parser.add_argument('-o','--output_gens', action='store_true', default=False, help = 'boolean flag to output every 100th generation. ' +
                                    'This can be helpful in tracking simulation progression')
        parser.add_argument('-B','--backup', action='store_true', default=False, help = 'boolean flag to turn on backups of the simulations, ' +
                                    'allowing a restart of simulations if required. This will increase space and time complexity')
        parser.add_argument('-P', '--polymorphisms', action='store_true', default=False, help = 'boolean flag to turn on the creation of file ' +
                                    'specifying all polymorphic and fixed states at the end of a branch')
        parser.add_argument('-S', '--calculate_selection', action='store_true', default=False, help = 'boolean flag that turns on calculations ' +
                                    'of selection by counting synonymous and non-synonymous fixed substitutions')


        #Get arguments from user
        arguments = parser.parse_args()
        
        return (arguments)
        
    
    
    #Command to make string version of mutation matrix from csv file
    def make_mutation_matrix(self, mutation_matrix):
        mut_mat = pandas.read_csv(mutation_matrix, names = ["A","C","G","T"])

        nrow = mut_mat.shape[0]

        #Check that mutational matrices are either 4 by 4 or 4 by 64
        if ((nrow != 4 and nrow !=64) or (mut_mat.shape[1] != 4)):
            print("Mutational matrices must be either 4 by 4 or 4 by 64. Representing mutations from " +
                "nucleotide to nucleotide or tri-nucleotide to nucleotide, respectfully.")
            sys.exit(0)

        mut_mat = mut_mat.to_numpy()

        #Check to make sure that mutations from nucleotide to itself are 0
        if(nrow == 4):
            diag_sum = sum(np.diag(mut_mat))
        else:
            col_1 = mut_mat[:,0]
            col_2 = mut_mat[:,1]
            col_3 = mut_mat[:,2]
            col_4 = mut_mat[:,3]
            diag_sum = sum(col_1[0:4] + col_1[16:20] + col_1[32:36] + col_1[48:52] +
                            col_2[4:8] + col_2[20:24] + col_2[36:40] + col_2[52:56] +
                            col_3[8:12] + col_3[24:28] + col_3[40:44] + col_3[56:60] +
                            col_4[12:16] + col_4[28:32] + col_4[44:48] + col_4[60:64])

        if(diag_sum != 0):
            print("All mutations from a nucleotide to itself must be 0. ie. in 4 by 4 " +
                "mutation matrices, all diagonals must be 0 and in 4 by 64 mutation matrices, " +
                "the first 4 rows in column 1 must be 0, the second 4 rows in column 2 must be 0, etc.")
            sys.exit(0)


        #Make and return string of the mutational matrix
        mut_mat_str = "matrix(c(" + str(list(mut_mat.flatten()))[1:-1] + "), ncol = 4, byrow = T)"
        return mut_mat_str    
        
        
        
   #Go through arguments and make sure that the user hasn't provided any arguments that cannot be processed     
    def check_arguments (self, arguments):
        
        arguments_pass = True
    
        #Check to see if user has selected high performance computing. If using ensure that time and partition are given
        hpc = arguments.high_performance_computing
        if (hpc and (arguments.partition == None or arguments.time == None)):
            print("When using high performance computing, partition and time data must be provided. Closing program.")
            arguments_pass = False

        #Check to make sure gene count and coding ratio are valid
        if (arguments.gene_count < 0 or arguments.gene_count > arguments.genome_length):
            print("Number of genes must be greater than 0 and less than the length of the genome. Closing program.")
            arguments_pass = False;

        if (arguments.coding_ratio < 0 or arguments.coding_ratio > 1.0):
            print("Coding ratio must be greater than 0 and less than or equal to 1. Please re-enter as a ratio (0, 1]. Closing program.")
            arguments_pass = False;


        #Ensure that if users the user has only one gene present if using a fasta file
        if (arguments.fasta_file != None and (arguments.coding_ratio != 1.0 or arguments.gene_count != 1)):
            print("When specifying an ancestral sequence with a fasta file, the sequence of only one fully coding gene should be provided. " +
                    "Closing program.")
            arguments_pass = False;            
            
        return(arguments_pass)
        
    
    #Read input tree to make output file names
    def process_filenames(self, tree_name, backup):
    
        # Find where data needs to be output to, set up documents and folders accordingly accordingly
        output_file_start = os.getcwd() + '/' + tree_name.split('.')[0]
        split_starting_output = output_file_start.split('/')
        
        #Find names of directories
        output_files_directory = "/".join(split_starting_output[0:(len(split_starting_output)-1)]) + "/slimScripts"
        
        #Make output file directories
        try:
            os.mkdir(output_files_directory)
        # Catch recreation of folder, allow user to specify whether files should be overwritten
        except OSError:
            cont = input ("The required directories already exits, program files will be overwritten, continue (y/n)?")
            while(cont != "y"):
                if(cont == "n"):
            	    sys.exit(0)
                cont = input("Continue? Please enter y or n")
                
        
        #Set up where files will be output to 
        output_filename = output_files_directory + "/" + split_starting_output[-1]
        
        
        #Make backup folder
        if(backup):
            backup_files_directory = "/".join(split_starting_output[0:(len(split_starting_output)-1)]) + "/backupFiles"
            try:
                os.mkdir(backup_files_directory)
            except OSError:
                print("using same backup folder")
            
            return((output_filename, output_file_start, backup_files_directory))
        else:
            return((output_filename, output_file_start, None))
        
        
    
    
    
    #Read through user input and process into a dictionary of starting parameters
    def make_param_dict(self, arguments):
    
        # Set up the starting parameters
        param_dict = vars(arguments)
        
        #Figure out burn in time
        param_dict["burn_in"] = arguments.burn_in_multiplier * arguments.population_size
        del param_dict["burn_in_multiplier"]
      
        #Figure out if jukes_cantor
        param_dict["jukes_cantor"] = arguments.mutation_matrix == None

        
        return(param_dict)
        



    #Save input given by user to a file so that it can be viewed/remembered after running the simulation
    def save_input(self, param_dict):
        
        file_start = param_dict["filenames"][1]
        
        # If this is a jukes cantor model write add theta to dictionary
        if (param_dict["jukes_cantor"]):
            param_dict["theta"] = 4*param_dict["mutation_rate"]*param_dict["population_size"]
        
        # Save starting parameters for later reference
        parameter_file = open(file_start + "_parameters.yaml", 'w')
        yaml.dump(param_dict,parameter_file)
        parameter_file.close()
        
       
        
    
    #Return the dictionary of starting parameters
    def get_params(self):
        return(self.param_dict)
        
        



if __name__ == '__main__':
    input_reader = readInput()
    input_reader.process_input()