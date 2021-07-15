#Program to take in gb file and fasta file and create ancestral sequence, coding ranges and fitness profiles
#Required packages:
#Bio, random, pandas, numpy, sys

import random, pandas, sys
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq

class getUserDefinedSequence:
    
    #Initialize a gb file and fasta file given by the user to find the ancestral sequence and coding regions from. Also initialize stationary distributions
    def __init__(self, gb_file, fasta_file, stationary_dists = None, fitness_dists = None):
        self.gb_file = gb_file
        self.fasta_file = fasta_file
        self.stationary_dists = stationary_dists
        self.fitness_dists = fitness_dists
        
        if (self.stationary_dists != None):
            self.stationary_dists = self.stationary_dists.transpose()
            self.stationary_dists = self.stationary_dists.drop("neutral")

    #Get the coding ranges of sequences from the gb file
    def get_coding_features(self):
        #Find the sequence to be taken
        for record in SeqIO.parse(self.gb_file, "gb"):
            features = record.features
            

    
        #Get the coding sequence numbers from the sequences
        cds = []

        try:    
            for feature in features:
                if(feature.type=="CDS"):
                    cds.append(str(feature.location))
        except:
            print("Please provide gb file in genbank format. Program closing.")
            sys.exit(0)
            
        #Find the coding numbers as ints and format according to SLiM-Tree requirements
        orfs = []
        for seq in cds:
            if(seq.find("+") != -1):
                 
                
                if(seq.find("join") != -1):
                    for seq_part in seq.split(","):
                        digits = []
                        for num in seq_part.split(":"):
                            digits.append(int(''.join(s for s in list(num) if s.isdigit())))
                        if(((digits[1] - digits[0]) ) % 3 != 0): #If connected sequences - may not be exactly the correct codon count - add sequences to the end
                            digits[1] = digits[1] + 3 - (((digits[1] - digits[0]) ) % 3)
                        orfs.append(digits)
                
                else:
                    digits = []
                    for num in seq.split(":"):
                        digits.append(int(''.join(s for s in list(num) if s.isdigit())))
                
                    orfs.append(digits)

        orfs = np.stack(np.sort(np.array(orfs), axis = 0))

        return orfs
    

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
    
