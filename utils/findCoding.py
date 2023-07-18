#Script which finds coding sequences from user input

import math
import numpy as np

class findCoding:
    def __init__(self, genome_length, coding_ratio = None, gene_count = None):
        self.genome_length = genome_length
        self.coding_ratio = coding_ratio
        self.gene_count = gene_count
        
        if(coding_ratio == None):
            self.get_coding_seqs_fasta()
        else:
            self.get_coding_seqs()
        
    
    
    
    
    #Return the range of coding sequences for the set number of genes.
    def get_coding_seqs(self):
        
        if(self.gene_count == 0 or self.coding_ratio == 0):
            return None
        
        percent_coding = math.ceil(int(self.genome_length) * self.coding_ratio) #Gives approximate number of coding amino acids
        avg_coding_length = math.ceil(percent_coding / self.gene_count) #Gives avg length of coding region
        avg_noncoding_length = 0
        if (self.gene_count != 1):
            avg_noncoding_length = math.floor((self.genome_length - percent_coding) / (self.gene_count - 1)) #Average length of non-coding regions by subtracting number of coding aa from total aa

        coding_regions = []
        current_aa = 0

        for i in range (self.gene_count):
            coding_regions.append(current_aa)
            coding_regions.append(min(current_aa + avg_coding_length, self.genome_length - 1)) #Ensures that you are not surpasssing the length of the genome
            current_aa = current_aa + avg_noncoding_length + avg_coding_length  #Accounts for the non-coding region + coding region added previously
            
            #Make sure that the genome will not be longer than the genome
            if (current_aa + avg_coding_length > self.genome_length - 1):
                current_aa = self.genome_length - avg_coding_length - 1

        self.coding_regions = np.stack(np.array_split(coding_regions, self.gene_count))
        
        
        
        
    #If user defined sequence the coding seqs will be decided just by the genome length
    def get_coding_seqs_fasta(self):
        self.coding_regions = np.stack(np.array_split([int(0), int(self.genome_length)], 1))

    
    #Returns coding regions
    def get_coding_regions(self):
        return(self.coding_regions)
   

