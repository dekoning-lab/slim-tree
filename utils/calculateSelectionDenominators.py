

import statistics, sys
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio.Seq import CodonTable
from collections import Counter

#Class which calculates the selection denominators for dn ds analysis
class calculateSelectionDenominators:
    
    #fitness_profile_nums holds the profile number of every codon in the genome. sites lists the codon positions whose
    #mutational opportunity is counted (see selection_sites); None counts every codon
    def __init__(self, stationary_distributions, fitness_profile_nums, mu, short_mu_mat, sites = None):
        
        #Set up variables
        self.stationary_distributions = stationary_distributions
        
        #Get the list of codons 
        self.codons = list(self.stationary_distributions.index)
        self.ncodons = len(self.stationary_distributions.index)
        
        #Find synonymous and non-synonymous positions in stationary dists
        syn_subs = self.find_syn_codons(stationary_distributions)
        
        
        #Find rates of mutation
        if(short_mu_mat is None): #If there is a not a mutation matrix, we need to make 4 by 4 mutation matrix from rate of mutation
            short_mu_mat = np.full((4, 4), mu/3)
        else: #Grab the correct matrix - ie. the matrix that is in array form
            short_mu_mat = short_mu_mat[0]
            
        self.mu_mat = self.find_rate_mut(short_mu_mat)

        #Only the counted sites contribute to the denominators
        if(sites is not None):
            fitness_profile_nums = [fitness_profile_nums[site] for site in sites]

        self.dS_denom, self.dN_denom = self.calculate_selection_denominators(fitness_profile_nums, syn_subs)
    

    # Grab the possible codons and split into synonymous and non-synonymous from the stationary distributions
    def find_syn_codons(self, stationary_distributions):
    
        #Find amino acids associated with the codons
        codons = Seq("".join(list(self.codons)))
        AAs = np.array([*str(codons.translate())])
        
        syn_subs = []
        
        #Recurse through codons and find the positions of synonymous codons
        for codon_num in range(self.ncodons):
            #Create boolean lists of matching and not matching AAs
            amino_acid = AAs[codon_num]
            AAs_to_match = AAs.copy()
            
            matching_AAs = AAs_to_match == amino_acid
            syn_subs.append(matching_AAs.tolist())
            
        return(syn_subs)   
    
    
    
    # Find all possible codons with only 1 mutation and find mutation rate for that codon
    def find_rate_mut(self, short_mu_mat):
        #Set up mutation matrix
        mu_mat = np.full((self.ncodons, self.ncodons), 0.0)
        
        #Set up conversion from nucleotide character to position in short mutation matrix
        convert_nucleotide = {"A": 0, "C": 1, "G" : 2, "T" : 3}
        
        #Recurse through codons and find mutation rate
        for cod1num in range(self.ncodons):
            cod1 = self.codons[cod1num]
            cod1_split = set(enumerate([*cod1]))
            
            
            for cod2num in range(self.ncodons):
                cod2 = self.codons[cod2num]
                cod2_split = set(enumerate([*cod2]))
                differences = list(cod1_split.symmetric_difference(cod2_split))
                num_dif = len(differences)/2
                
                if(num_dif != 1):
                    continue
               
                #Get the position of the mutation in the codon
                pos_mut = differences[0][0]
                old_codon_num = convert_nucleotide[cod1[pos_mut]]
                new_codon_num = convert_nucleotide[cod2[pos_mut]]
                
                mu_mat[cod1num, cod2num] = short_mu_mat[old_codon_num, new_codon_num]
                
        return (mu_mat)
            
            
            
    
    # Loop through each amino acid in the stationary distributions and calculates piQ
    # The ratio of amino acids in the genome that are from the distribution
    def get_dist_ds_dn(self, profile_num, profile_quantity, syn_subs):
        
        stationary_dist = self.stationary_distributions.iloc[:,profile_num]
        denom_ds = 0
        denom_dn = 0
        
        #Recurse through codons and find pi_Qij
        for i in range(self.ncodons):       
            for j in range(self.ncodons):
            
                Qij = self.mu_mat[i][j]
                pi_Qij = stationary_dist.iloc[i]*Qij*profile_quantity
                if(syn_subs[i][j]):
                    denom_ds += pi_Qij
                else:
                    denom_dn += pi_Qij

        return([denom_ds, denom_dn])
    
    
    
    # Calculates the denominators of dS and dN: the numbers of synonymous and nonsynonymous sites in the counted codons
    def calculate_selection_denominators(self, fitness_profile_nums, syn_subs):    
        
        # Find how many sites use each profile, keyed by profile number. A profile that no site uses is simply absent
        num_each_profile = sorted(Counter(fitness_profile_nums).items())

        # No counted sites: no mutational opportunity
        if(len(num_each_profile) == 0):
            return (np.zeros(2))
        
        # Find the value of ds and dn for each of the profiles that is used
        ds_dn = np.array([self.get_dist_ds_dn(profile_num, quantity, syn_subs) for profile_num, quantity in num_each_profile])

        
        #Sum the mutational opportunities over sites: sum of pi_i * mu_ij for synonymous and for nonsynonymous changes
        opportunity = sum(ds_dn)
        if(sum(opportunity) == 0):
            return (np.zeros(2))

        #Goldman & Yang (1994): the proportions of the opportunity that are synonymous and nonsynonymous (their rho_s
        #and rho_a at V = infinity, which sum to 1), times 3, are the numbers of synonymous and nonsynonymous nucleotide
        #sites per codon. Multiplied by the number of counted codons they are the numbers of sites in the sequence, so
        #that count / denominator is Ks or Ka: substitutions per synonymous or per nonsynonymous site
        return (3 * len(fitness_profile_nums) * opportunity / sum(opportunity))
        
        
    #Return dn to public
    def get_dn(self):
        return(self.dN_denom)
    
    
    #Return ds to public
    def get_ds(self):
        return(self.dS_denom)



#Codon positions that count towards dN/dS: the codons strictly inside a gene. Each gene's start codon and terminal stop
#codon, and all non-coding codons, are left out: they have no synonymous or non-synonymous changes to score under the
#profiles (they use SLiM-Tree's "neutral" profile, a column of ones that is not a distribution over codons).
#coding_seqs holds one [start codon, terminal stop codon] pair per gene. In neutral runs (-N) there are no
#genes, so every codon counts
def selection_sites(coding_seqs, genome_length, neutral_evolution = False):
    if(neutral_evolution or coding_seqs is None):
        return(list(range(int(genome_length))))

    sites = []
    for start, end in np.asarray(coding_seqs).reshape(-1, 2):
        sites += list(range(int(start) + 1, int(end)))
    return(sites)
