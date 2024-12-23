#Script which takes in stationary distributions and processes the file. Also finds or processes the fitness file
import sys, subprocess, os, random, re, time
import pandas as pd
import numpy as np
from Bio.Seq import Seq
from Bio.Seq import CodonTable
from Bio import SeqIO
from copy import deepcopy

class findFitness:

    def __init__(self, stationary_dist_file, neutral):
    
        #All possible codons
        self.codons = ['TTT', 'TTC','TTA', 'TTG', 'CTT', 
                        'CTC', 'CTA', 'CTG', 'ATT', 'ATC', 
                        'ATA', 'ATG', 'GTT', 'GTC', 'GTA', 
                        'GTG', 'TCT', 'TCC', 'TCA', 'TCG',
                        'CCT', 'CCC', 'CCA', 'CCG', 'ACT', 
                        'ACC', 'ACA', 'ACG', 'GCT', 'GCC', 
                        'GCA', 'GCG', 'TAT', 'TAC', 'CAT', 
                        'CAC', 'CAA', 'CAG', 'AAT', 'AAC', 
                        'AAA', 'AAG', 'GAT', 'GAC', 'GAA', 
                        'GAG', 'TGT', 'TGC', 'TGG','CGT', 
                        'CGC', 'CGA', 'CGG', 'AGT', 'AGC', 
                        'AGA', 'AGG', 'GGT', 'GGC', 'GGA', 'GGG']
        
        
        if (neutral):
            self.stationary_mat = pd.DataFrame(index=np.arange(61), columns=np.arange(0))
            self.stationary_mat.insert(0, 0, 1/61)
            #Every amino acid is represented
            self.stationary_mat.index = self.codons
            self.validify_stationary_distribution(self.stationary_mat)
        else:
            self.stationary_dist_file = stationary_dist_file
            self.stationary_mat = pd.read_csv(stationary_dist_file, header = None, index_col = 0)
            self.ndists = self.stationary_mat.shape[1]
            self.validify_stationary_distribution(self.stationary_mat)
            
            
            
            #Read csv of codon numbers
            slim_codon_nums = os.path.join( os.path.dirname( __file__ ), '..' ) + '/fitnessDataFiles/slim_codon_nums.csv'
            self.slim_codons = pd.read_csv(slim_codon_nums, header = 0, index_col = 0).transpose().to_dict('list')
        
        
    
    
    #Function to make sure that the stationary distributions are provided in the correct format
    def validify_stationary_distribution(self, stationary_mat):
        #Translate to amino acids
        
        try:
            codons = Seq("".join(list(stationary_mat.index)))
            self.AAs = [*str(codons.translate())]
        except TypeError: #Make sure that codons are provided in the first row
            print("Please ensure the first row of your stationary distributions is the codon names. Exiting.")
            sys.exit(0)
        except CodonTable.TranslationError: #Check to make sure only valid codons are provided
            print("Please ensure that your stationary distribution only has valid codons. Exiting.")
            sys.exit(0)
        
        #Check to make sure stop codons are not provided
        if("*" in self.AAs):
            print("Do not include stop codons in your stationary distribution. Exiting.")
            sys.exit(0)
            
        #Check to make sure that all codons are represented in the stationary distribution
        all_codons_set = set(self.codons)
        given_codons_set = set (list(stationary_mat.index))
        
        if(all_codons_set != given_codons_set):
            print("Please ensure that every codon is represented in your stationary distributions. Exiting.")
            sys.exit(0)
            
    
    
    
    #If an existing fitness file is given, read in the fitness file, verify the file and put amino acids in the right order
    def process_existing_fitness_file(self, fitness_file):
        #Read fitness matrix
        fitness_mat = pd.read_csv(fitness_file, header = None, index_col = 0)
        
        #Verify that the fitness matrix is set up properly
        self.validify_fitness_file(fitness_mat)
        
        #Reorder columns in order required for the software
        fitness_mat = fitness_mat.sort_index(axis=0)
        
        self.fitness_mat = fitness_mat
        
        
        
    #Function verifies that there are the same number of stationary distributions as fitnesses and that fitness files are in amino acids,
    #while stationary distributions are in codons
    def validify_fitness_file(self, fitness_mat):
    
        #Verify that there are the same number of stationary distributions and fitness values
        nfitnesses = fitness_mat.shape[1]
        
        if(self.ndists!=nfitnesses):
            print("The same number of fitness profiles and stationary distributions must be provided. Exiting.")
            sys.exit(0)
            
        
        #Verify that the fitness profiles are in terms of amino acids and that all are provided
        npossibilities = fitness_mat.shape[0]
        if(npossibilities != 21 or set(fitness_mat.index) != set(self.AAs + ["X"]) ):
            print("Fitness data files must be in terms of amino acids. There should be the 20 amino acids and stops (ie. 21 rows). Exiting.")
            sys.exit(0)
            
            
            
    #Function which finds fitnesses from stationary distributions using R script if fitness file not given.
    #Writes fitnesses to new file so they may be reused in the future
    def find_optimal_fitnesses(self, mutation_rate, population_size, hpc, partition, time_p):
        print("Finding fitnesses for stationary distributions")
        fitness_mat =  os.getcwd() + "/table_fitness_dists.csv"
        
        #If existing fitness mat exists remove it, ensures HPC waits if file is already there
        if(os.path.isfile(fitness_mat)): os.remove(fitness_mat) 
        
        
        #If hpc make file to run R script to find fitness profiles and run batch file
        if (hpc):
            batch_file = open("find_fitness.sh", "w")
            batch_file.write(("#!/bin/sh\n\n#SBATCH -J find_fitness \n#SBATCH -t " + str(time_p) +
                "\n#SBATCH -p "  + str(partition) + 
                "\n#SBATCH -o fitness.out\n#SBATCH -e fitness.err" +
                "\n#SBATCH -n 10" + 
                "\nRscript " + os.path.dirname(os.path.realpath(__file__)) + "/fitness_profile_finder.R" +
                " -f " + self.stationary_dist_file + " -N " + str(population_size) + " -v " + str(mutation_rate)+ " -o " + fitness_mat))
            batch_file.close()
            
            subprocess.run(["sbatch",  "find_fitness.sh"])
            
            #Waits until batch file finished before continuing - need to throw message if R script fails
            while not (os.path.isfile(fitness_mat) and os.path.getsize("fitness.err") == 0): 
                time.sleep(1)
            
        else:      
            #Run R script to find fitness profiles
            subprocess.call(["Rscript", os.path.dirname(os.path.realpath(__file__)) + "/fitness_profile_finder.R", 
                    "-f", self.stationary_dist_file, "-N", str(population_size), "-v", str(mutation_rate), "-o", fitness_mat])
        
        #Check if fitness matrix exists, if it does not exist error is in R script
        if(os.path.isfile(fitness_mat)):
            self.fitness_mat = pd.read_csv(fitness_mat, header = None, index_col = 0)
        else:
            print("There seems to be a bug in running the R script to get fitnesses. Check your packages. Closing program.")
            sys.exit(0)
        
    
    #Function to find fitnesses if a non-jukes-cantor matrix is supplied - uses the average mutation rate of the matrix
    def find_optimal_fitnesses_mu_mat(self, mutation_matrix, population_size, hpc, partition, time_p, test = False):
        #Find mean of all non-diagonal elements
        sum_all = mutation_matrix.values.sum()
        nrow = mutation_matrix.shape[0]
        
        mu =  (sum_all - np.diag(mutation_matrix).sum())/(nrow**2 - nrow)

        #Find the fitnesses if not a test, if a test just return the mean to test
        if(test):
            return(mu)
        else:
            self.find_optimal_fitnesses(mu, population_size, hpc, partition, time_p)




    #Process the fitness distributions so that they may be used by slim-tree
    def process_fitness_dists(self):
 
        # Add an extra distribution for non-coding neutral alleles - neutral, 20 AA + stop codon
        self.stationary_mat["neutral"] = 1
        self.fitness_mat["neutral"] = 1


        #Convert fitness profiles to a dictionary for easier integration into eidos
        fitness_profiles = self.fitness_mat.transpose().to_dict('list')
        
        
        #Find the smallest fitness values
        fitness_distributions = self.fitness_mat.values.tolist()
        min_fitness = min(np.array(fitness_distributions).flatten())
        
        return(fitness_profiles, min_fitness)
        
        
        
    #Set up the fitness profile for each position in the genome
    def define_fitness_profiles(self, randomize_fitness_profiles, coding_poses, genome_length):


        # Set up fitness profiles
        if(randomize_fitness_profiles): #User did not provide a sequence
            fitness_profile_nums = []
            neut_fitness_pos = self.fitness_mat.shape[1] - 1
            for coding_pos in range(len(coding_poses)):
                fitness_profile_nums = (fitness_profile_nums + [neut_fitness_pos] + #starting codon is correct
                            random.choices(range(neut_fitness_pos),k=coding_poses[coding_pos,1] - coding_poses[coding_pos,0] - 1))

                if (coding_pos != len(coding_poses) - 1):
                    fitness_profile_nums = fitness_profile_nums + [neut_fitness_pos]* (coding_poses[coding_pos+1,0] - 
                                    coding_poses[coding_pos,1])
                else:
                    fitness_profile_nums = fitness_profile_nums + [neut_fitness_pos]* (genome_length - 
                                    coding_poses[coding_pos,1])
        
        else: #Use user provided fitness profiles either provided or taken from the stationary dists
            
            if(self.ndists != genome_length): 
                print("Please ensure that when using a fasta file, the same number of fitness profiles are provided as the length " +
                "of the genome in the fasta file. Exiting.")
                sys.exit(0)

            # A fitness profile is given for every position in genome
            fitness_profile_nums = list(range(genome_length))
            
        return(fitness_profile_nums)
        
        
        
    
    
    # Randomly select an ancestral sequence from the codon seqence files if a fasta file is not provided
    def find_ancestral(self, coding_regions, profile_nums):
        
        #Set up ORFs
        start_codon_nums = coding_regions[:,0]
        stop_codon_nums = coding_regions[:,1]

        #Methionine - start codon
        start_codon = "14"

        #Stop codons
        stop_codons = ["48", "50", "56"]

        #Middle codons - chosen according to distribution of alleles
        codons = []

        for dist_num in profile_nums:
                weights = list(self.stationary_mat.iloc[:,dist_num])
                codon_name = random.choices(self.stationary_mat.index, weights = weights, k = 1)[0]
                codons.append(str(self.slim_codons[codon_name][0]))

        #Replace start and stop codons with start and stop codons
        for codon_num in start_codon_nums: codons[codon_num] = start_codon
        for codon_num in stop_codon_nums: codons[codon_num] = random.choice(stop_codons)

        return (codons)



    #Randomly select the ancestral sequence from a neutral models
    def find_ancestral_neutral(self, genome_length):
        return list(str(element) for element in np.random.randint(low = 0,high=63,size=genome_length))

        
    #Check to make sure file is in terms of nucleotides
    def dna_check(self,sequence):
        dna = set('ACTG')
        return all(base.upper() in dna for base in sequence)
        
        
        
    # Find the ancestral sequence when a fasta file is given by the user 
    def find_ancestral_fasta(self, fasta_file):
        
        #Read in fasta file
        try:
            num_record = 0
            for record in SeqIO.parse(fasta_file, "fasta"):
                ans_seq = str(record.seq).upper()
                num_record += 1

        except:
            print("Please provide ancestral sequence file in fasta format. Exiting.")
            sys.exit(0)
        
        
        #Make sure only one sequence uploaded
        if(num_record == 0):
            #Impropper formatting 
            print("Please provide ancestral sequence file in fasta format. Exiting.")
            sys.exit(0)
        elif(num_record > 1):
            #Too many sequences
            print("You have provided more than one sequence in your fasta file. " +
                    "Please provide only one sequence. Exiting.")
            sys.exit(0)
            
        
        #Check to make sure file is in terms of nucleotides
        if(not self.dna_check(ans_seq)):
            print("Please ensure that your fasta file is in terms of nucleotides, not amino acids. Exiting.")
            sys.exit(0)
        
        #Find the length of the genome in terms of amino acids 
        genome_length = len(ans_seq)/3
        
        #Convert the fasta file to slim codon numbers by splitting into codons
        codon_names = re.findall('...', ans_seq)
        
        codons = []
        for codon_name in codon_names:
            codons.append(str(self.slim_codons[codon_name][0]))
            
        return (codons, int(genome_length))
        
    
    
    #Convert stationary distribution from codons to amino acids
    def convert_stat_dist(self):
        combined_stat_dists = deepcopy(self.stationary_mat)
        combined_stat_dists.index = self.AAs
        combined_stat_dists = combined_stat_dists.groupby(combined_stat_dists.index).sum()
        
        return(combined_stat_dists)
    
    
    
    
    # Calculates expected fitness from fitness profiles for scaling in non-wright-fisher models
    def find_fitness_scaling(self, fitness_profile_nums, multiple_genes):
        #Set up lists and data frames
        stationary_dists = self.convert_stat_dist()
        expected_fitness_profiles = []
        

        # Find the expected value for each fitness profile
        for num_dist in range(len(stationary_dists.columns)-1):
            mean = 0
            stationary = stationary_dists.iloc[:,num_dist]
            fitness = self.fitness_mat.iloc[:,num_dist]
            for num_profile in range(len(fitness)-1):
                mean += stationary.iloc[num_profile] * fitness.iloc[num_profile]

            expected_fitness_profiles.append(mean)

        # Add fitness profile with mean of 1 to account for the neutral areas
        expected_fitness_profiles.append(1)


        # Find the expected value for each site in the genome based on it's respective fitness profile
        expected_fitnesses = []
        for fitness_profile in fitness_profile_nums:
            expected_fitnesses.append(expected_fitness_profiles[fitness_profile])
            
    
        # Find the expected value of all sites
        scaling_value = np.sum(expected_fitnesses)
        
        return(scaling_value)
    
    
    #Return the stationary matrix for use in other scripts
    def get_stationary_mat(self):
        return(self.stationary_mat)
        
        
    


if __name__ == '__main__':
    fit = findFitness(sys.argv[1])
    
    if(sys.argv[2] != None):
        fit.process_existing_fitness_file(sys.argv[2])
