#!/bin/sh

#SBATCH -J SLiM_Simulation_p1
#SBATCH -t 04:00:00
#SBATCH -p apophis
#SBATCH -o p1.out
#SBATCH -e p1.err
#SBATCH -n 1

slim '/work/dk_lab/work/afarinesh/SLiM_Tree_paper_MBE/EX_11/slimScripts/tree_p1.slim'