#!/bin/sh

#SBATCH -J SLiM_Simulation_p2
#SBATCH -t 04:00:00
#SBATCH -p apophis
#SBATCH -o p2.out
#SBATCH -e p2.err
#SBATCH -n 1

slim '/work/dk_lab/work/afarinesh/SLiM_Tree_paper_MBE/EX_10/slimScripts/tree_p2.slim'