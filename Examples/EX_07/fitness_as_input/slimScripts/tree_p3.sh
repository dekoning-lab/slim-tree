#!/bin/sh

#SBATCH -J SLiM_Simulation_p3
#SBATCH -t 04:00:00
#SBATCH -p apophis
#SBATCH -o p3.out
#SBATCH -e p3.err
#SBATCH -n 1

slim '/work/dk_lab/work/afarinesh/SLiM_Tree_paper_MBE/EX_07/slimScripts/tree_p3.slim'