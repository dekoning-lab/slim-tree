#!/bin/sh

#SBATCH -J SLiM_Simulation_p8
#SBATCH -t 04:00:00
#SBATCH -p apophis
#SBATCH -o p8.out
#SBATCH -e p8.err
#SBATCH -n 1

slim '/work/dk_lab/work/afarinesh/SLiM_Tree_paper_EMB/EX_2/slimScripts/tree_p8.slim'