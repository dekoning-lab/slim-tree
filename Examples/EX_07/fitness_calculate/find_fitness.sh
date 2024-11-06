#!/bin/sh

#SBATCH -J find_fitness 
#SBATCH -t 04:00:00
#SBATCH -p apophis
#SBATCH -o fitness.out
#SBATCH -e fitness.err
#SBATCH -n 10
Rscript /home/afarinesh.panahy/software/slim-tree/utils/fitness_profile_finder.R -f table_stationary_distributions.csv -N 100 -v 0.0015 -o /work/dk_lab/work/afarinesh/SLiM_Tree_paper_MBE/EX_07/fitness_calculate/table_fitness_dists.csv