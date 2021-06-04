# SLiM-Tree

SLiM-Tree is a simulation tool which automates pure population genetics simulations over phylogenetic timescales under realistic models of sequence-fitness relationships. For a full description of SLiM-Tree usage please refer to the user manual. 

SLiMTree requires installation of Python3 and slim (https://messerlab.org/slim/). Required python packages are sys, argparse, BioPython, matplotlib, random, pandas, numpy, os, json, string and math

To run SLiMTree run the command python3 SLiMTree -i <your input tree - in newick format>.


Additional arguments include:

      -T: tool - specify whether you are using SLiM-Tree or SLiM-Tree-HPC (for parallelization)
	  
	  -p: partition - the partition which scripts should be written to for SLiM-Tree-HPC
	  
	  -t: time - the maximum amount of time that scripts should be run for before timing out on SLiM-Tree-HPC
	  
	  -n: population size - the size of each population
	  
	  -v: the mutation rate
	  
	  -g: the length of the genome in codons
	  
	  -r: recombination rate
	  
	  -b: burn in multiplier - a value which will be multiplied by the population size to determine the length of the burn-in period to establish mutations

      -k: sample size - the size of the sample of genomes to be output at the end of the simulation
	  
	  -G: number of genes
	  
	  -C: the percentage of the genome which is coding
	  
	  -c: boolean specifying whether supstitutions should be counted
	  
	  -o: boolean specifying whether the generation of the simulation should be output (helps to keep track of longer simulations)
	  
	  -B: boolean specifying whether simulations should be backed up
	  
	  -w: boolean specifying whether a Wright-Fisher model should be used


Users may also specify another file under -d with changes to parameters for different branches. 

The folder DataPostProcessing contains scripts that were can be used for post processing of the output data and the folder. ExampleBatchScripts contains scripts used to run SLiM-Tree for some of the simulations we ran.
