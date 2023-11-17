

import statistics, sys
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Seq import CodonTable
from collections import Counter

#Class which calculates the selection denominators for dn ds analysis
class calculateSelectionDenominators:
    
    def __init__(self, stationary_distributions, fitness_profile_nums, mu, short_mu_mat):
        
        #Set up variables
        self.stationary_distributions = stationary_distributions
        
        #Get the list of codons 
        self.codons = list(stationary_distributions.index)
        self.ncodons = len(self.stationary_distributions.index)
        
        #Find synonymous and non-synonymous positions in stationary dists
        self.find_syn_codons(stationary_distributions)
        
        
        #Find rates of mutation
        if(short_mu_mat == None): #If there is a not a mutation matrix, we need to make 4 by 4 mutation matrix from rate of mutation
            short_mu_mat = np.full((4, 4), mu/3)
        else: #Grab the correct matrix - ie. the matrix that is in array form
            short_mu_mat = short_mu_mat[0]
            
        self.mu_mat = self.find_rate_mut(short_mu_mat)
        
        self.dN_denom, self.dS_denom = self.calculate_selection_denominators(fitness_profile_nums)
    

    # Grab the possible codons and split into synonymous and non-synonymous from the stationary distributions
    def find_syn_codons(self, stationary_distributions):
    
        #Find amino acids associated with the codons
        codons = Seq("".join(list(stationary_distributions.index)))
        AAs = np.array([*str(codons.translate())])
        
        self.syn_subs = []
        
        #Recurse through codons and find the positions of synonymous codons
        for codon_num in range(self.ncodons):
            #Create boolean lists of matching and not matching AAs
            amino_acid = AAs[codon_num]
            AAs_to_match = AAs.copy()
            
            #Ensure that same amino acid isn't included
            AAs_to_match[codon_num] = None
            matching_AAs = AAs_to_match == amino_acid
            
            self.syn_subs.append(matching_AAs.tolist())
            
            
    
    
    
    # Find all possible codons with only 1 mutation and find mutation rate for that codon
    def find_rate_mut(self, short_mu_mat):
        #Set up mutation matrix
        mu_mat = np.full((self.ncodons, self.ncodons), 0.0)
        
        #Set up conversion from nucleotide character to position in short mutation matrix
        convert_nucleotide = {"A": 0, "C": 1, "G" : 2, "T" : 3}
        
        #Recurse through codons and find mutation rate
        for cod1 in range(self.ncodons):
            cod1_split = set(enumerate([*self.codons[cod1]]))
            
            for cod2 in range(self.ncodons):
                cod2_split = set(enumerate([*self.codons[cod2]]))
                differences = list(cod1_split.symmetric_difference(cod2_split))
                num_dif = len(differences)/2
                
                if(num_dif != 1):
                    continue
               
                old_codon_num = convert_nucleotide[list(differences[1])[1]]
                new_codon_num = convert_nucleotide[list(differences[0])[1]]
                
                mu_mat[cod1, cod2] = short_mu_mat[old_codon_num, new_codon_num]
                
        return (mu_mat)
            
            
            
    
    # Loop through each amino acid in the stationary distributions and calculates piQ
    # ratio is the ratio of amino acids in the genome that are from the distribution
    def get_dist_ds_dn(self, profile_num, profile_quantity):
        
        stationary_dist = self.stationary_distributions.iloc[:,profile_num]
        denom_ds = 0
        denom_dn = 0
    
    
        #Recurse through codons and find pi_Qij
        for i in range(self.ncodons):       
            for j in range(self.ncodons):
            
                Qij = self.mu_mat[i][j]
                pi_Qij = stationary_dist[i]*Qij*profile_quantity
                
                if(self.syn_subs[i][j]):
                    denom_ds += pi_Qij
                else:
                    denom_dn += pi_Qij
        
        # print(ratio)
        return([denom_ds, denom_dn])
    
    
    
    # Calculates the denominator of the dN and dS true equations based on the fitness profiles numbers
    def calculate_selection_denominators(self, fitness_profile_nums):    
        
        ndists = range(self.stationary_distributions.shape[1]-1)
        
        # Find the ratio of each profile in this genome
        num_each_profile = Counter(fitness_profile_nums)
        num_each_profile = pd.DataFrame(sorted([(i, num_each_profile[i] ) for i in num_each_profile]))
          
        # Find the value of ds for each of the stationary distributions
        ds_dn = np.array(list(map(lambda x: self.get_dist_ds_dn(x, num_each_profile.loc[x,1]), ndists)))

        
        #Sum and multiply by 3 to get the denominators of dn and ds
        return (sum(ds_dn)*3)
        
        
    #Return dn to public
    def get_dn(self):
        return(self.dN_denom)
    
    
    #Return ds to public
    def get_ds(self):
        return(self.dS_denom)
    
    
    
    
    
    
    
    
    