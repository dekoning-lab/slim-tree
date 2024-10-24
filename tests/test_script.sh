#!/bin/sh

#SBATCH -J run_script
#SBATCH -t 12:00:00
#SBATCH -p cpu2019
#SBATCH -o script.out
#SBATCH -e script.err
#SBATCH -n 1

module load R/4.2.0
python3 ~/slim-tree/ ~/slim-tree/tests/test_tree.txt ~/slim-tree/tests/table_stationary_distributions.csv -hpc -t 12:00:00 -p cpu2019 -o -S <<< 'y'