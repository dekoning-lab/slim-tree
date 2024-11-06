#!/bin/sh

#SBATCH -J SLiM_Simulation_p7
#SBATCH -t 04:00:00
#SBATCH -p apophis
#SBATCH -o p7.out
#SBATCH -e p7.err
#SBATCH -n 1

slim '/work/dk_lab/work/afarinesh/SLiM_Tree_paper_EMB/EX_1/slimScripts/tree_p7.slim'