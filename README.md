# SLiM-Tree

SLiM-Tree is a simulation tool which automates pure population genetics simulations over phylogenetic timescales under realistic models of sequence-fitness relationships. For a full description of SLiM-Tree usage please refer to the user manual. 

SLiMTree requires installation of Python3 and slim (https://messerlab.org/slim/). Required python packages are sys, argparse, BioPython, matplotlib, random, pandas, numpy, os, json, string and math

To run SLiMTree run the command python3 SLiMTree -i <your input tree - in newick format>.


Additional arguments include:
	
  	-h: help - displays a help message showing options 

	-T: tool - specify whether you are using SLiM-Tree or SLiM-Tree-HPC (for parallelization)
	  
	-p: partition - the partition which scripts should be written to for SLiM-Tree-HPC
	  
	-t: time - the maximum amount of time that scripts should be run for before timing out on SLiM-Tree-HPC
	  
	-n: population size - the size of each population
	  
	-v: the mutation rate
	  
	-g: the length of the genome in codons
	  
	-r: recombination rate
	  
	-b: burn in multiplier - a value which will be multiplied by the population size to determine the length of the burn-in period to establish mutations

        -k: sample size - the size of the sample of genomes to be output at the end of the simulation, users can specify number, 'all' for all genomes or 'consensus' for the consensus sequence
	
	-sr: the ratio of the population that is split into each subpopulation in non-WF models of evolution
	  
	-G: the number of genes to be used in the genome - must be 1 if using protein based fitness effects
	  
	-C: the percentage of the genome which is coding - must be 1 if using protein based fitness effects
	
	-hp: boolean specifying whether to model haploidy instead of diploidy
	  
	-c: boolean specifying whether supstitutions should be counted
	  
	-o: boolean specifying whether the generation of the simulation should be output (helps to keep track of longer simulations)
	  
	-B: boolean specifying whether simulations should be backed up
	  
	-w: boolean specifying whether a Wright-Fisher model should be used
	
	-P: boolean specifying whether polymorphic states and percentages should be output
	
	-s: boolean specifying whether the user provides a sequence
	
	-f: link to fasta file for user provided sequences
	
	-gb: link to genbank file for user profided fitnesses
	
	-R: boolean specifying whether fitness profiles should be randomized, there must be an equal number of fitness profiles to the genome length if false
	
	-fc: boolean specifying whether fitness effects should be calculated using site-heterogeneous fitness profiles. If false, fitness effects will be calculated using protein based structure and a pdb file must be provided
	
	-pdb: link to pdb file containing structure to model
	
	-pdbs: link to folder containing structurally diverse pdbs of approximately the same size as the supplied pdb to be used as a distribution of epistatic interactions
	
	-cid: the chain id of the chain to be modelled in the main pdb
	
	-cids: file containing the chain ids to be used in the distribution of epistatic interactions
	
	-ct: the maximum distance that amino acids may be apart to be in contact with eachother
	
	-jc: boolean specifying whether a Jukes-Cantor mutational matrix should be used. If false, users must specify a csv file containing mutational matrix.
	
	-m: csv file containing mutational matrix defining the mutation rate from one nucleotide to another in alphabetical order (ie. A, C, G, T)
	


Users may also specify another file under -d with changes to parameters for different branches. 

The folder DataPostProcessing contains scripts that were can be used for post processing of the output data and the folder. 
