

import statistics, sys
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Seq import CodonTable
from collections import Counter

#Class which calculates the selection denominators for dn ds analysis
class calculateSelectionDenominators:
    
    def __init__(self, stationary_distributions, fitness_profile_nums, mu, mu_mat):
        
        #Set up variables
        self.stationary_distributions = stationary_distributions
        self.ncodons = len(self.stationary_distributions.index)
        
        #Find synonymous and non-synonymous positions in stationary dists
        self.find_syn_codons(stationary_distributions)
        self.find_num_muts(stationary_distributions)
        
        #Find rates of mutation
        jc = mu_mat == None #If there is a mutation matrix then it's not a jukes-cantor model
        self.mu_mat = self.find_mu_mat(jc, mu_mat, mu)
        
        
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
            
            
    
    # Find all possible codons with more than 1 mutation
    def find_num_muts(self, stationary_distributions):
        #Get the list of codons 
        codons = list(stationary_distributions.index)
        
        #Recurse through codons and find which ones have more than 1 mutation
        self.num_muts = []
        for cod1 in codons:
            muts_list = []
            for cod2 in codons:
                cod1_split = [*cod1]
                cod2_split = [*cod2]
                num_same = len(set(enumerate(cod1_split)).intersection(set(enumerate(cod2_split))))
                muts_list.append(3-num_same)
            
            self.num_muts.append(muts_list)
            
            
         
    
    # Converts a mutation matrix or single value (for jukes-cantor), to a list of average forward and reverse mutation rates
    def find_mu_mat (self, jc, mu_mat, mu):
        
    
        # Make mutation matrix of the mu value if Jukes-Cantor
        if (jc):
            mu_mat = np.full((61, 61), mu/3)
        
        
        #Need to convert other mutation matrices to larger matrix
        

        
        
        
        return(mu_mat)
    

        
    
    
    # Loop through each amino acid in the stationary distributions and calculates piQ
    # ratio is the ratio of amino acids in the genome that are from the distribution
    def get_dist_ds_dn(self, profile_num, profile_quantity):
        
        stationary_dist = self.stationary_distributions.iloc[:,profile_num]
        denom_ds = 0
        denom_dn = 0
    
    
        #Recurse through codons and find pi_Qij
        for i in range(self.ncodons):       
            for j in range(self.ncodons):
                #Ensure that only a single substitution is being added
                if (self.num_muts[i][j] != 1):
                    continue
           
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
        print((sum(ds_dn)*3))
        sys.exit(0)
        
        #Sum and multiply by 3 to get the denominators of dn and ds
        return (sum(ds_dn))
        
        
    #Return dn to public
    def get_dn(self):
        return(self.dN_denom)
    
    
    #Return ds to public
    def get_ds(self):
        return(self.dS_denom)
    
    
    
    
    
    
    
    
    