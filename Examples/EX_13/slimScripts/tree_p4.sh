#!/bin/sh

#SBATCH -J SLiM_Simulation_p4
#SBATCH -t 04:00:00
#SBATCH -p apophis
#SBATCH -o p4.out
#SBATCH -e p4.err
#SBATCH -n 1

slim '/work/dk_lab/work/afarinesh/SLiM_Tree_paper_MBE/EX_14/slimScripts/tree_p4.slim'