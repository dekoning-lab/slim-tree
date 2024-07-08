#Program to take in a tree in Newick file and output SLiM code to run a simulation from that tree

import utils
import os

class SLiMTree:

    #Main script to run other commands
    def __init__(self):

        start_params = self.read_input()
        start_params = self.process_fitness(start_params)
        
        
        #Get the data for each clade in the tree
        clades = utils.cladeReader.cladeReader(start_params)
        clade_dict_list = clades.get_clade_dict_list()

        #Write and run the slim code
        self.write_slim_code(clade_dict_list, start_params)
        
        
        # Start the SLiM code to run
        if(start_params["high_performance_computing"]):
            os.system("sbatch \"" + start_params["filenames"][0] + "_p1.sh\"")
        else:
           os.system("slim \"" + start_params["filenames"][0] + "_p1.slim\"")


    #Read the input from the user
    def read_input(self):
    
        #Read and process input from user commands
        input_reader = utils.readInput.readInput()
        input_reader.process_input()
        start_params = input_reader.get_params()
        
        
        #Save input to file
        input_reader.save_input(start_params)
        return(start_params)
        
        
        
        
    #Process fitnesses from the starting parameters
    def process_fitness(self, start_params):
        #Read in the stationary distributions, processes fitness profiles and find ancestral sequences
        fitness_finder = utils.findFitness.findFitness(start_params["codon_stationary_distributions"])
        if(start_params["aa_fitness_distributions"] != None):
            fitness_finder.process_existing_fitness_file(start_params["aa_fitness_distributions"])
        elif (start_params["jukes_cantor"]):
            fitness_finder.find_optimal_fitnesses(start_params["mutation_rate"], 
                        start_params["population_size"], start_params["high_performance_computing"], 
                        start_params["partition"], start_params["time"])
        else: 
            fitness_finder.find_optimal_fitnesses_mu_mat(start_params["mutation_matrix"][0], 
                        start_params["population_size"], start_params["high_performance_computing"], 
                        start_params["partition"], start_params["time"])
                        
        start_params["fitness_profiles"], start_params["min_fitness"]  = fitness_finder.process_fitness_dists()
        
        
        #Get coding regions and find the ancestral sequence and assign fitness profiles to coding regions
        if(start_params["fasta_file"] == None):
            start_params["coding_seqs"] = utils.findCoding.findCoding(start_params["genome_length"], start_params["coding_ratio"],
                        start_params["gene_count"]).get_coding_regions()
            start_params["fitness_profile_nums"] = fitness_finder.define_fitness_profiles(True,
                    start_params["coding_seqs"], start_params["genome_length"])
            start_params["ancestral_sequence"] = fitness_finder.find_ancestral(start_params["coding_seqs"], 
                        start_params["fitness_profile_nums"])
        else:
            start_params["ancestral_sequence"], start_params["genome_length"] = fitness_finder.find_ancestral_fasta(start_params["fasta_file"])
            start_params["coding_seqs"] = utils.findCoding.findCoding(start_params["genome_length"]).get_coding_regions()
            start_params["fitness_profile_nums"] = fitness_finder.define_fitness_profiles(False,
                    start_params["coding_seqs"], start_params["genome_length"])
            
            
            
        # Find the scaling factor for mostly neutral fitnesses
        start_params["scaling_value"] = fitness_finder.find_fitness_scaling(start_params["fitness_profile_nums"], 
                    start_params["coding_ratio"] != 1)
                        
                        
        #If dN/dS is being calculated find the denominators
        if(start_params["calculate_selection"]):
            sel_denom = utils.calculateSelectionDenominators.calculateSelectionDenominators(fitness_finder.get_stationary_mat(),
                        start_params["fitness_profile_nums"], start_params["mutation_rate"], start_params["mutation_matrix"])
            start_params["dn_denom"] = sel_denom.get_dn()
            start_params["ds_denom"] = sel_denom.get_ds()
        
        return(start_params)


  
    #Using the dictionary of the data from the clades, write slim code to represent the phylogeny
    def write_slim_code (self, clade_dict_list, start_params):

        #Open SLiM writer based on tool type and write the initialize statement
        if(start_params["high_performance_computing"]):
            slim_writer = utils.writeSLiMHPC.writeSLiMHPC(start_params)
        else:
            slim_writer = utils.writeSLiM.writeSLiM(start_params)

        #Write a script for each clade which will be run sequentially
        if (start_params["nonWF"]):
            #If this is a non-Wright-Fisher model, use a different write_subpop function
            for clade in clade_dict_list:
                slim_writer.write_subpop_nonwf(clade);
        else:            
            for clade in clade_dict_list:
                slim_writer.write_subpop(clade)

        if(not start_params["high_performance_computing"]):
            slim_writer.close_file()




if __name__=='__main__':
    SLiMTree()
