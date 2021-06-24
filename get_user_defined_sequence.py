#Program to take in gb file and fasta file and create ancestral sequence, coding ranges and fitness profiles
#Required packages:
#Bio, random

import random
from Bio import SeqIO
from Bio.Seq import Seq

class get_user_defined_sequence:
    
    #Initialize a gb file and fasta file given by the user to find the ancestral sequence and coding regions from
    def __init__(self, gb_file, fasta_file):
        self.gb_file = gb_file
        self.fasta_file = fasta_file

    #Get the coding ranges of sequences from the gb file
    def get_coding_features(self):
        #Find the sequence to be taken
        for record in SeqIO.parse(self.gb_file, "gb"):
            features = record.features
    
        #Get the coding sequence numbers from the sequences
        cds = []                                                                                                          
        for feature in features:
            if(feature.type=="CDS"):
                cds.append(str(feature.location))
        
        #Find the coding numbers as ints and format according to SLiM-Tree requirements
        orfs = []
        for seq in cds:
            if(seq.find("+") == -1):
                pass
            digits = []                                                                                                   
            for num in seq.split(":"):
                digits.append(''.join(s for s in list(num) if s.isdigit()))
    
            print(digits)
            if(len(orfs) == 0 or orfs[-1][1] < digits[0] ): #If two orfs overlap - combine into a single orf (this has to happen for SLiM)
                orfs.append(digits)
            else:
                orfs[-1][1] = digits [1]
                
        return orfs
    

    #Find the ancestral sequence from the fasta file given by the user
    def get_ancestral_sequence(self):
        for record in SeqIO.parse(self.fasta_file, "fasta"):
            sequence = str(record.seq)
            self.seq = record.seq
        
        if (len(sequence) % 3 == 0): #If the sequence is divisible by 3, return the sequence
            return sequence
        else: #If the sequence is not divisible by 3, return the sequence plus additional nucleotides (SLiM requires the sequence to be divisible by 3)
            num_missing_seqs = 3 - (len(sequence) % 3)
            additional_seqs = "".join(random.choices( ["A", "C", "G", "T"], k=num_missing_seqs))
            return sequence + additional_seqs
            
    
    #Find the fitness profiles for a given coding sequence
    def find_fitness_profiles(self, sequence):
        pass
        #print(Seq(sequence).translate()[0:444])
    



#Testing code for whent the class is not imported
if __name__ == "__main__":
    get_seq = get_user_defined_sequence("../carsonella.gb", "../carsonella.fasta")
    get_seq.get_coding_features()
    ans_seq = get_seq.get_ancestral_sequence()
    get_seq.find_fitness_profiles(ans_seq)
    