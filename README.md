# SLiM-Tree

• Here, we present SLiM-Tree a simulation tool that applies population genetics techniques to phylogenetic timescales without relying on conventional models (under realistic models of sequence-fitness relationships). It is a flexible software package to automate seting up pure population genetics simulations along a phylogeny with a realistic model of molecular evolution.

• How it works: It Employs SLiM to create a platform which can evolve populations with or without using Wright-Fisher model, allowing users to explore by relaxing simplified assumptions and models.

For a full description of SLiM-Tree usage please refer to the user manual. 

• Requirements: SLiMTree requires installation of Python3 and slim (https://messerlab.org/slim/). Required python packages are sys, argparse, BioPython, matplotlib, random, pandas, numpy, os, json, string and math. If using protein based fitness effects - java and c are also required.

• Usage: slim-tree [-h] [-fd AA_FITNESS_DISTRIBUTIONS] [-hpc] [-p PARTITION] [-t TIME] [-w] [-n POPULATION_SIZE] [-b BURN_IN_MULTIPLIER]
                 [-r RECOMBINATION_RATE] [-v MUTATION_RATE] [-m MUTATION_MATRIX] [-d TREE_DATA_FILE] [-g GENOME_LENGTH] [-G GENE_COUNT]
                 [-C CODING_RATIO] [-f FASTA_FILE] [-k SAMPLE_SIZE] [-sr SPLIT_RATIO] [-c] [-o] [-B] [-P] [-S]
                 input_tree codon_stationary_distributions


• How to run SLiM_Tree: run the command python3 ../slim-tree/ <input_tree> <codon_stationary_distributions> 


Additional arguments include:
	
  	-h: --help  
   		show this help message and exit
   	
    	-hpc: --high_performance_computing 
     		boolean flag to turn on slim-tree high performance computing. Slurm is required

	-fd AA_FITNESS_DISTRIBUTIONS: --aa_fitness_distributions AA_FITNESS_DISTRIBUTIONS 
 					file containing a amino acid fitnesses
	  
	-p PARTITION: --partition PARTITION 
 			partition to run Slurm on - required if using high performance computing
	  
	-t TIME: --time TIME  
 		maximum time to run each simulation for - suggested time is the maximum time available for a partition -required if using high
   		performance computing.

   	-w: --nonWF           
    	    boolean flag to specify that a non-wright-fisher model should be used in lieu of a wright-fisher model.
  
 	-n POPULATION_SIZE: --population_size POPULATION_SIZE 
  			    starting population size for the simulation, default = 100

  	-b BURN_IN_MULTIPLIER: --burn_in_multiplier BURN_IN_MULTIPLIER 
   				value to multiply population size by for burn in, default = 10

	-r RECOMBINATION_RATE: --recombination_rate RECOMBINATION_RATE
 				recombination rate, default = 2.5e-8

 	-v MUTATION_RATE: --mutation_rate MUTATION_RATE
                          starting mutation rate for the simulation, default = 2.5e-6

	-m MUTATION_MATRIX: --mutation_matrix MUTATION_MATRIX
                            CSV file specifying a mutation rate matrix, matrix should be either 4 by 4 or 4 by 64 specifying rates from nucleotide
                            to nucleotide and tri-nucleotide to nucleotide respectfully. Nucleotides and tri-nucleotides should be in alphabetical
                            order with no headers. If mutation rate matrix is supplied, mutation rate will be ignored
  
   	-d TREE_DATA_FILE: --tree_data_file TREE_DATA_FILE
                           file to change the population size for specific branches using YAML formatting. When using HPC, other parameters may
                           also be changed.
			
	-g GENOME_LENGTH: --genome_length GENOME_LENGTH
                          length of the genome - amino acids, default = 300
			
  	-G GENE_COUNT: --gene_count GENE_COUNT
                	number of genes to be simulated by the model, default = 1
			
  	-C CODING_RATIO: --coding_ratio CODING_RATIO
                        ratio of the genome which is coding, default = 1.0
			
  	-f FASTA_FILE: --fasta_file FASTA_FILE
                        fasta file containing ancestral sequence (amino acids), replaces random creation of ancestral sequence. Fitness
                        profiles for each amino acid are required
			
  	-k SAMPLE_SIZE: --sample_size SAMPLE_SIZE
                        size of sample obtained from each population at a tree tip at the end of the simulations.Input 'all' for the every
                        member of the tree tip samples and consensus for the consensus sequence of the population at each tip. default = all
			
 	-sr SPLIT_RATIO: --split_ratio SPLIT_RATIO
                         proportion of a population that goes into the first daughter branch at a tree branching point in non-wright fisher
                         models. must be ratio between 0 and 1.0. default = 0.5
			
  	-c: --count_subs      
   	    boolean flag to turn on substitution counting. This will slow down simulations
   
  	-o: --output_gens
   	    boolean flag to output every 100th generation. This can be helpful in tracking simulation progression
   
  	-B: --backup
   	    boolean flag to turn on backups of the simulations, allowing a restart of simulations if required. This will increase space and 
	    time complexity
			
  	-P: --polymorphisms
   	    boolean flag to turn on the creation of file specifying all polymorphic and fixed states at the end of a branch
   
  	-S: --calculate_selection
            boolean flag that turns on calculations of selection by counting synonymous and non-synonymous fixed substitutions


The folder DataPostProcessing contains scripts that can be used for post processing of the output data and the folder. 
