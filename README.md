# SLiM-Tree

• Here, we presents SLiM-Tree, a flexible software package developed to automate the use of SLiM (https://messerlab.org/slim/) to perform pure population genetics simulations along a phylogeny, using a realistic model of molecular evolution, without the standard suite of simplifying assumptions used in phylogenetics (e.g., mutation limited evolution/ weak mutation or infinite sites, limited mutation among competing allelic types between fixations, no polymorphism, etc.).


• SLiM-Tree requires installation of SLiM (https://github.com/MesserLab/SLiM/releases/download/v3.7.1/SLiM_Manual.pdf, last version 4.0.1), python3 (version 3.10 or later) and R (version 4.2.2 or later). It is also important to avoid using "conda" for these installations, as it can cause significant issues with version compatibility.

• For Python3 and R, their dependencies listed below should also be installed:

	•Required python packages: sys, argparse, BioPython, matplotlib, random, pandas, numpy, os, json, string, Pyyaml, math, statistics, collections, copy, subprocess, re, time, shutil, builtins, filecmp, contextlib, unittest, io and csv.
 	•Required R packages:  dplyr, BB, data.table, optparse, seqinr, doParallel and Rfast.

• The SLiM-Tree package can be installed by cloning the repository:

	$ git clone git@github.com:dekoning-lab/slim-tree.git
	$ git checkout redev (will be dev by the time we want to publish)

• How to run SLiM_Tree: run the command $python3 ../slim-tree/ <input_tree> <codon_stationary_distributions.csv> 
• Additional parameters and flags are available and detailed in the Parameters section. For a full description of SLiM-Tree usage please refer to the user manual. 

• Additional parameters include:
	
  	-h: 
   	Show this help message and exit
   	
	-hpc: 
 	Boolean flag to turn on slim-tree high performance computing. Slurm is required

	-fd aa_fitness_distributions.csv:
 	The SLiM-Tree is designed to generate a table of amino acid fitness distributions using the provided stationary distributions. However, you also have the option to input your table of fitness distribution in .csv format using `-fd' flag.
	  
	-p partition_name:
 	Partition to run Slurm on - required if using high performance computing (-hpc)
	  
	-t time(00:00:00):
 	Maximum time to run each simulation for, required if using high performance computing (-hpc)

   	-w --nonWF:
    	Boolean flag to specify that a non-wright-fisher model should be used in lieu of the default wright-fisher model. 
  
 	-n population_size (default = 100):
  	Starting population size for the simulation. 

  	-b burn_in_multiplier (default = 10):
   	Value to multiply population size by for burn in.

	-r recombination_rate (default = 2.5e-8):
 	Recombination rate. 

 	-v mutation_rate (default =2.5e-6):
        Starting mutation rate for the simulation.

	-m mutation_matrix.csv
	A file in .csv format, specifying a mutation rate matrix. Matrix should be either 4*4 or 4*64, specifying rates from nucleotide to nucleotide and tri-nucleotide to nucleotide respectfully. Nucleotides and tri-nucleotides should be in alphabetical order with no headers. If mutation rate matrix is supplied, mutation rate will be ignored
  
   	-d tree_data_file.yaml:
        File to change the population size (-n) or mutation rate (-v) for specific branches using .yaml formatting. 
	When using -hpc, other parameters may also be changed.
			
	-g genome_length (default = 300):
        Length of the genome in amino acids.
			
  	-G gene_count (default = 1):
        Number of genes to be simulated by the model
			
  	-C coding_ratio (default = 1.0):
        Ratio of the genome which is coding.
			
  	-f fasta_file.fasta:
        Fasta file containing ancestral sequence (amino acids), replaces random creation of ancestral sequence. Fitness profiles for each amino acid are required
			
  	-k sampling_method (default = all):
        Size of sample obtained from each population at a tree tip at the end of the simulations.Input "consensus" for the consensus sequence of the population at each tip. 
			
 	-sr split_ratio (default = 0.5):
        Proportion of a population that goes into the first daughter branch at a tree branching point in non-wright fisher models. must be ratio between 0 and 1.0. 
			
  	-c:      
   	Boolean flag to turn on substitution counting. This will slow down simulations
   
  	-o:
   	Boolean flag to output every 100th generation. This can be helpful in tracking simulation progression
   
  	-B: 
   	Boolean flag to turn on backups of the simulations, allowing a restart of simulations if required. This will increase space and time complexity. 
			
  	-P:
   	Boolean flag to turn on the creation of file specifying all polymorphic and fixed states at the end of a branch
   
  	-S:
        Boolean flag that turns on calculations of selection by counting synonymous and non-synonymous fixed substitutions


• Usage: slim-tree <input_tree> <codon_stationary_distributions.csv> [-h] [-fd aa_fitness_distributions.csv] [-hpc] [-p partition] [-t time] [-w --nonWF] [-n population_size] [-b burn_in_multiplier] [-r recombination_rate] [-v mutation_rate] [-m mutation_matrix] [-d tree_data_file.yaml] [-g genome_length] [-G gene_count] [-C coding_ratio] [-f fasta_file.fasta] [-k sampling_method] [-sr split_ratio] [-c count_substitutions] [-o output_gene] [-B backup] [-P polymorphisms] [-S calculate_selection_dn/ds]


The folder DataPostProcessing contains scripts that can be used for post processing of the output data and the folder. 
