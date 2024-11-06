#!/bin/sh

#SBATCH -J SLiM_Simulation_p6
#SBATCH -t 04:00:00
#SBATCH -p apophis
#SBATCH -o p6.out
#SBATCH -e p6.err
#SBATCH -n 1

slim '/work/dk_lab/work/afarinesh/SLiM_Tree_paper_EMB/EX_4/slimScripts/tree_p6.slim'