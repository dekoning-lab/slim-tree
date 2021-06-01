# SLiMTree

SLiMTree is a Wright-Fisher simulation tool that simulates evolution along a phylogeny using SLiM written by the Messer Lab. 

SLiMTree requires installation of Python3 and slim (https://messerlab.org/slim/). Required python packages are sys, argparse, BioPython, matplotlib, random, csv, numpy, os and json

To run SLiMTree run the command python3 SLiMTree -i <your input tree - in newick format>.


Additional arguments include:

      -b: burn in multiplier - a value which will be multiplied by the population size to determine the brun in before the tree splits

      -g: the length of the genome in amino acids

      -k: sample size - the size of the sample of populations output at the end of the simulation

      -n: population size - the size of each population

      -r: recombination rate

      -v: the mutation rate



Users may also specify another file under -d which gives this data on specific populations.

Note: for the time being SLiMTree is only able to run on high performance computing clusters as different branches are run through different nodes. 


The folder DataPostProcessing contains scripts that were can be used for post processing of the output data and the folder ExampleBatchScripts contains scripts used to run SLiM for my thesis.

Finally, the folder CompiledData contains compiled raw substitution and branch length counts from my thesis project. 
