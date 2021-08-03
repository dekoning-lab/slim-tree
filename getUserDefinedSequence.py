#Program to take in gb file and fasta file and create ancestral sequence, coding ranges and fitness profiles
#Required packages:
#Bio, random, pandas, numpy, sys

import random, pandas, sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

class getUserDefinedSequence:
    
    #Initialize a gb file and fasta file given by the user to find the ancestral sequence and coding regions from. Also initialize stationary distributions
    def __init__(self, fasta_file, stationary_dists = None, fitness_dists = None):
        self.gb_file = gb_file
        self.fasta_file = fasta_file
        self.stationary_dists = stationary_dists
        self.fitness_dists = fitness_dists
        
        if (self.stationary_dists != None):
            self.stationary_dists = self.stationary_dists.transpose()
            self.stationary_dists = self.stationary_dists.drop("neutral")
    

    #Find the ancestral sequence from the fasta file given by the user
    def get_ancestral_sequence(self):
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            sequence = str(record.seq)
        
        try:    
            return sequence
        except:
            print("Please provide fasta file in fasta format. Program closing.")
            sys.exit(0)
            
        
            
    
    #Find the fitness profiles for a given coding sequence
    def find_fitness_profiles(self, sequence, coding_features):
        fitness_profiles = []
        expected_fitnesses = []
        
        #Filter through each amino acid and choose fitness profile appropriately
        for feature in coding_features:
            coding_seq = list(Seq(sequence[feature[0]: feature[1]]).translate())
            
            for aa in coding_seq:
                #print(aa)
                if (aa == "*"):
                    profile = [len(self.stationary_dists.index)]
                else:
                    scaled_data = self.stationary_dists[aa]/sum(self.stationary_dists[aa])
                    profile = random.choices(range(len(self.stationary_dists.index)), weights = scaled_data, k = 1)
                    expected_profile_mean = sum(scaled_data * self.fitness_dists[aa][:-1])
                    expected_fitnesses.append(expected_profile_mean)
                    expected_fitnesses.append(expected_profile_mean)
                fitness_profiles = fitness_profiles + profile
                
            
            #sys.exit(0)
        self.expected_fitness = np.prod(expected_fitnesses)
        
        return fitness_profiles
        
        
    #Find the expected fitness of each amino acid for scaling in non-WF models
    def get_fitness_scaling(self):
        return self.expected_fitness
    



#Testing code for whent the class is not imported
if __name__ == "__main__":
    stationary_distributions = pandas.io.parsers.read_csv("fitnessDataFiles/table_stationary_distributions.csv") 
    get_seq = getUserDefinedSequence("../carsonella.gb", "../carsonella.fasta", stationary_distributions)
    coding_feats = get_seq.get_coding_features()
    ans_seq = get_seq.get_ancestral_sequence()
    get_seq.find_fitness_profiles(ans_seq, coding_feats)
    
