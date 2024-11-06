#!/bin/bash

# Set the job name and specify the cluster's job scheduler directives
# Replace "your_job_name" with your desired job name
# Replace the directives below with the appropriate directives for your cluster
# These directives are for Slurm as an example

#SBATCH --job-name=ex6
#SBATCH --partition=apophis
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=05:00:00

# Navigate to the directory where your R script and Perl script are located
cd /work/dk_lab/work/afarinesh/SLiM_Tree_paper_EMB/EX_6


python3 /home/afarinesh.panahy/software/slim-tree tree table_stationary_distributions.csv -hpc -t 04:00:00 -p apophis -v 0.0015 -fd table_fitness_dists.csv -n 250 


# Exit
exit 0 
