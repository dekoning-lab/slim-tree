#!/bin/sh

#SBATCH -J SLiM_Simulation_p9
#SBATCH -t 04:00:00
#SBATCH -p apophis
#SBATCH -o p9.out
#SBATCH -e p9.err
#SBATCH -n 1

slim '/work/dk_lab/work/afarinesh/SLiM_Tree_paper_EMB/EX_4/slimScripts/tree_p9.slim'