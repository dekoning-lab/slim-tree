#!/bin/sh

#SBATCH -J SLiM_Simulation_p5
#SBATCH -t 04:00:00
#SBATCH -p apophis
#SBATCH -o p5.out
#SBATCH -e p5.err
#SBATCH -n 1

slim '/work/dk_lab/work/afarinesh/SLiM_Tree_paper_EMB/EX_3/slimScripts/tree_p5.slim'